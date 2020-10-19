#include "bin-packing.h"
#include "population.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define CLOCK2SEC(START, END) \
        (((double)END - START) * ((double)1.0 / CLOCKS_PER_SEC))

#if defined (DEBUG_BIN)
static void print_bin(const bin_t *bin) {
        printf("fill: %Lf\n"
               "count: %zu\n"
               "item_indices:\n",
               bin->fill, bin->count);
        for (size_t i=0; i<bin->count; i++) {
                printf("%zu ", bin->item_indices[i]);
        }
        putchar('\n');
}
static void print_chrom(const chrom_t *chrom) {
        printf("fitness: %lf\n"
               "bin_cap: %Lf\n"
               "num_bins: %zu\n"
               "bins:\n",
               chrom->fitness, chrom->bin_cap, chrom->num_bins);
        for (size_t i=0; i<chrom->num_bins; i++) {
                printf("bin %zu:\n", i);
                print_bin(chrom->bins[i]);
        }
}
#endif

static struct llarray *result_bin_alloc(const bin_t *bin,
                                        const long double *item_sizes,
                                        size_t num_items) {
        struct llarray *arr = malloc(offsetof(struct llarray, elems)
                                     + (bin->count * sizeof(*arr->elems)));
        arr->num_elems = bin->count;
        for (size_t i=0; i<arr->num_elems; i++) {
                arr->elems[i] = item_sizes[bin->item_indices[i]];
        }
        return arr;
}
static result_t *result_alloc(const chrom_t *best_chrom,
                              const long double *item_sizes, size_t num_items) {
        assert(best_chrom != NULL);
        result_t *res = malloc(offsetof(result_t, bins)
                               + (best_chrom->num_bins * sizeof(*res->bins)));
        *res = (result_t){.fitness = best_chrom->fitness,
                          .num_bins = best_chrom->num_bins};
        for (size_t i=0; i<res->num_bins; i++) {
                res->bins[i] = result_bin_alloc(best_chrom->bins[i],
                                                item_sizes, num_items);
        }
        return res;
}
void result_free(result_t *res) {
        for (size_t i=0; i<res->num_bins; i++) {
                free(res->bins[i]);
        }
        free(res);
}

typedef pop_t tourn_t;
static tourn_t *tournament_select(const pop_t *pop, size_t mating_pool_size,
                                  double tournament_p,
                                  unsigned tournament_size) {
        /* initialize container for tournament */
        tourn_t *mp = pop_alloc(mating_pool_size);
        /* fill mating pool through tournament selection */
        for (size_t i=0; i<mating_pool_size; i++) {
                mp->chroms[i] = pop->chroms[rand() % pop->num_chroms];
                /* apply selection based on chosen tournament size */
                for (unsigned j=1; j<tournament_size; j++) {
                        size_t k = rand() % pop->num_chroms;
                        if (mp->chroms[i]->fitness
                            < pop->chroms[k]->fitness) {
                                mp->chroms[i] = pop->chroms[k];
                        }
                }
        }
        return mp;
}
static inline void tourn_free(tourn_t *mating_pool) {
        free(mating_pool);
}
static const chrom_t *find_elite(const pop_t *pop) {
        const chrom_t *elite = pop->chroms[0];
        for (size_t i=1; i<pop->num_chroms; i++) {
                if (pop->chroms[i]->fitness > elite->fitness) {
                        elite = pop->chroms[i];
                }
        }
        return elite;
}
static pop_t *child_pop(const tourn_t *mating_pool, size_t population_size,
                        const chrom_t *elite_chrom,
                        const long double *item_sizes, size_t num_items) {
        pop_t *child = pop_alloc(population_size);
        child->chroms[0] = chrom_copy(elite_chrom);
        for (size_t i=1; i<child->num_chroms; i++) {
                size_t i1 = rand() % mating_pool->num_chroms;
                size_t i2 = rand() % mating_pool->num_chroms;
                child->chroms[i] = chrom_cx(mating_pool->chroms[i1],
                                            mating_pool->chroms[i2],
                                            item_sizes, num_items);
        }
        return child;
}
static void mutate_pop(pop_t *pop, double mutation_rate,
                       const long double *item_sizes, size_t num_items) {
        /* starts at 1 because 0 contains the elite chromosome from
         * previous generations */
        for (size_t i=1; i<pop->num_chroms; i++) {
                chrom_mutate(pop->chroms[i], mutation_rate,
                             item_sizes, num_items);
        }
}
static inline void print_stats(size_t gen_num, const chrom_t *best_chrom,
                               double secs) {
        printf("%zu\t %zu\t %lf\t %lf\n",
               gen_num, best_chrom->num_bins, best_chrom->fitness, secs);
}
static int bin_cmp(const void *a, const void *b) {
        const bin_t *av = *(bin_t **)a;
        const bin_t *bv = *(bin_t **)b;
        if (av->fill < bv->fill) {
                return -1;
        } else if (av->fill > bv->fill) {
                return 1;
        } else {
                return 0;
        }
}
static void inversion(pop_t *pop) {
        /* start at 1 because 0 is the elite chromosome */
        for (size_t i=1; i<pop->num_chroms; i++) {
                qsort(pop->chroms[i]->bins, pop->chroms[i]->num_bins,
                      sizeof(*pop->chroms[i]->bins), bin_cmp);
        }
}

result_t *bin_packing(const prob_set_t *ps) {
        /* verify values before use */
        assert(ps->item_sizes != NULL);
        assert(ps->num_items > 0);
        assert(ps->bin_capacity > 0);
        for (size_t i=0; i<ps->num_items; i++) {
                assert(ps->item_sizes[i] <= ps->bin_capacity);
        }
        assert(ps->max_generations > 0);
        assert(ps->population_size > 0);
        assert(ps->mating_pool_size > 0);
        assert((ps->max_mutation_rate >= 0.0)
               && (ps->max_mutation_rate <= 1.0));
        assert((ps->tournament_p >= 0.0) && (ps->tournament_p <= 1.0));
        assert(ps->tournament_size > 0);

        clock_t start = clock();
        printf("gen #\t # bins\t fitness\t cum. sec\n");
        pop_t *pop = pop_rand_init(ps->bin_capacity, ps->population_size,
                                   ps->item_sizes, ps->num_items);
        const chrom_t *best = find_elite(pop);
        clock_t end = clock();
        print_stats(1, best, CLOCK2SEC(start, end));
        for (size_t gen=1; gen<ps->max_generations; gen++) {
                if ((best->fitness >= nextafter(1.0, 0.0))
                    || (best->num_bins <= ps->terminal_num_bins)
                    || (CLOCK2SEC(start, end) >= ps->max_secs)) {
                        break;
                }
                tourn_t *t = tournament_select(pop, ps->mating_pool_size,
                                               ps->tournament_p,
                                               ps->tournament_size);
                pop_t *child = child_pop(t, ps->population_size, best,
                                         ps->item_sizes, ps->num_items);
                tourn_free(t);
                if (ps->use_inversion_operator) {
                        inversion(child);
                }
                mutate_pop(child, ps->max_mutation_rate,
                           ps->item_sizes, ps->num_items);
                const chrom_t *new_best = find_elite(child);
                end = clock();
                print_stats(gen + 1, new_best, CLOCK2SEC(start, end));
#ifdef DEBUG_BIN
                print_chrom(new_best);
#endif
                pop_free(pop);
                pop = child;
                best = new_best;
        }
        result_t *res = result_alloc(best, ps->item_sizes, ps->num_items);
        pop_free(pop);
        return res;
}

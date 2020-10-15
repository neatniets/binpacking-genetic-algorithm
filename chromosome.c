#include "chromosome.h"
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#ifdef DEBUG
#include <stdio.h>
#endif

#define FITNESS_K       2

static bin_t *bin_alloc(void) {
        bin_t *bin = malloc(sizeof(*bin));
        *bin = (bin_t){.fill = 0,
                       .count = 0,
                       .item_indices = NULL};
        return bin;
}
static void bin_free(bin_t *bin) {
        if (bin == NULL) {
                return;
        }
        free(bin->item_indices);
        free(bin);
}
static bin_t *bin_copy(const bin_t *bin) {
        bin_t *copy = bin_alloc();
        *copy = (bin_t){.fill = bin->fill,
                        .count = bin->count,
                        .item_indices = malloc(bin->count
                                               * sizeof(*bin->item_indices))
                       };
        memcpy(copy->item_indices,
               bin->item_indices,
               bin->count * sizeof(*copy->item_indices));
        return copy;
}
static void bin_add(bin_t *bin, size_t index, long long value) {
        bin->fill += value;
        bin->item_indices = realloc(bin->item_indices,
                                    (bin->count + 1)
                                    * sizeof(*bin->item_indices));
        bin->item_indices[bin->count] = index;
        bin->count++;
}

static chrom_t *chrom_alloc(size_t bin_cap) {
        chrom_t *chrom = malloc(sizeof(*chrom));
        *chrom = (chrom_t){.fitness = 0,
                           .bin_cap = bin_cap,
                           .num_bins = 0,
                           .bins = NULL};
        return chrom;
}
static void chrom_new_bin(chrom_t *chrom) {
        chrom->bins = realloc(chrom->bins,
                              (chrom->num_bins + 1)
                              * sizeof(*chrom->bins));
        chrom->bins[chrom->num_bins] = bin_alloc();
        chrom->num_bins++;
}
static void chrom_add_bin(chrom_t *chrom, bin_t *bin) {
        chrom_new_bin(chrom);
        bin_free(chrom->bins[chrom->num_bins - 1]);
        chrom->bins[chrom->num_bins - 1] = bin;
}
static void chrom_del_bin(chrom_t *chrom, size_t bin_index) {
        bin_free(chrom->bins[bin_index]);
        memcpy(chrom->bins + bin_index,
               chrom->bins + bin_index + 1,
               sizeof(*chrom->bins) * (chrom->num_bins - (bin_index + 1)));
        /* chrom->bins = realloc(chrom->bins,
         *                       sizeof(*chrom->bins)
         *                       * (chrom->num_bins - 1)); */
        chrom->num_bins--;
}
static void eval_fitness(chrom_t *chrom, int fitness_k) {
        chrom->fitness = 0;
        for (size_t i=0; i<chrom->num_bins; i++) {
                chrom->fitness += pow((double)chrom->bins[i]->fill
                                      / chrom->bin_cap,
                                      fitness_k)
                                  / chrom->num_bins;
        }
}

static void fit(chrom_t *chrom, size_t index, long long value) {
        for (size_t i=0; i<chrom->num_bins; i++) {
                if (chrom->bins[i]->fill + value <= chrom->bin_cap) {
#ifdef DEBUG
                        printf("adding to bin %zu:\n"
                               "index: %zu\tsize: %lld\n",
                               i, index, value);
#endif
                        bin_add(chrom->bins[i], index, value);
                        return;
                }
        }
#ifdef DEBUG
        printf("adding to new bin:\n"
               "index: %zu\tsize: %lld\n",
               index, value);
#endif
        chrom_new_bin(chrom);
        bin_add(chrom->bins[chrom->num_bins - 1], index, value);
}
static void first_fit(chrom_t *chrom, const long long *item_sizes,
                      bool *is_item_used, size_t num_items,
                      size_t start_pos) {
        for (size_t loop = 0, i = start_pos, end = num_items;
             loop < 2;
             loop++, i = 0, end = start_pos) {
                for (; i < end; i++) {
                        if (!is_item_used[i]) {
                                fit(chrom, i, item_sizes[i]);
                                is_item_used[i] = true;
                        }
                }
        }
}
chrom_t *rand_first_fit(const long long *item_sizes, size_t num_items,
                        size_t bin_cap) {
        chrom_t *chrom = chrom_alloc(bin_cap);
        bool *is_item_used = calloc(num_items, sizeof(*is_item_used));
        first_fit(chrom, item_sizes, is_item_used, num_items,
                  rand() % num_items);
        free(is_item_used);
        eval_fitness(chrom, FITNESS_K);
        return chrom;
}
void chrom_free(chrom_t *chrom) {
        if (chrom == NULL) {
                return;
        }
        for (size_t i=0; i<chrom->num_bins; i++) {
                bin_free(chrom->bins[i]);
        }
        free(chrom->bins);
        free(chrom);
}

chrom_t *chrom_copy(const chrom_t *chrom) {
        chrom_t *copy = chrom_alloc(chrom->bin_cap);
        *copy = (chrom_t){.fitness = chrom->fitness,
                          .num_bins = chrom->num_bins,
                          .bins = malloc(chrom->num_bins
                                         * sizeof(*copy->bins))};
        for (size_t i=0; i<copy->num_bins; i++) {
                copy->bins[i] = bin_copy(chrom->bins[i]);
        }
        return copy;
}

/** Returns true if used items conflict with items in bin */
static bool check4conflict(const bool *is_item_used, const bin_t *bin) {
        for (size_t i=0; i<bin->count; i++) {
                if (is_item_used[bin->item_indices[i]]) {
                        return true;
                }
        }
        return false;
}
static void mark_used(bool *is_item_used, const bin_t *bin) {
        for (size_t i=0; i<bin->count; i++) {
                is_item_used[bin->item_indices[i]] = true;
        }
}
static void mark_unused(bool *is_item_used, const bin_t *bin) {
        for (size_t i=0; i<bin->count; i++) {
                is_item_used[bin->item_indices[i]] = false;
        }
}
chrom_t *chrom_cx(const chrom_t *parent1, const chrom_t *parent2,
                  const long long *item_sizes, size_t num_items) {
        size_t p2_start = rand() % (parent2->num_bins);
        size_t p2_count = rand() % (parent2->num_bins - p2_start) + 1;
        size_t p1_pos = rand() % (parent1->num_bins + 1);
#ifdef DEBUG
        printf("\nparent2 start: %zu\n"
               "parent2 count: %zu\n"
               "parent1 pos: %zu\n\n",
               p2_start, p2_count, p1_pos);
#endif
        chrom_t *child = chrom_alloc(parent1->bin_cap);
        bool *is_item_used = calloc(num_items, sizeof(*is_item_used));
        /* marking items from the chosen bins from parent2 as used */
        for (size_t i=p2_start; i<p2_start+p2_count; i++) {
                mark_used(is_item_used, parent2->bins[i]);
        }
        /* 1.) adding bins w/o conflicts from parent1 prior to p1_pos
         * 2.) adding bins from parent2
         * 3.) adding bins w/o conflicts from parent1 after p1_pos */
        for (size_t p1i = 0, p2i = p2_start, loop = 0,
             p1end = p1_pos, p2end = p2_start+p2_count;
             loop < 2;
             loop++, p1i = p1_pos, p1end = parent1->num_bins) {
                /* 1.) & 3.) */
                for (; p1i < p1end; p1i++) {
                        if (!check4conflict(is_item_used,
                                            parent1->bins[p1i])) {
#ifdef DEBUG
                                printf("bin %zu in parent1 has no"
                                       " conflicts; adding bin\n",
                                       p1i);
#endif
                                chrom_add_bin(child,
                                              bin_copy(parent1->bins[p1i]));
                                mark_used(is_item_used, parent1->bins[p1i]);
                        }
                }
                /* 2.) */
                for (; p2i < p2end; p2i++) {
#ifdef DEBUG
                        printf("adding bin %zu from parent2\n", p2i);
#endif
                        chrom_add_bin(child, bin_copy(parent2->bins[p2i]));
                }
        }
        /* first-fit the items not in any bins in child */
        first_fit(child, item_sizes, is_item_used, num_items, 0);
        free(is_item_used);
        eval_fitness(child, FITNESS_K);
        return child;
}
void chrom_mutate(chrom_t *chrom, double mutation_rate,
                  const long long *item_sizes, size_t num_items) {
        assert((mutation_rate >= 0.0) && (mutation_rate <= 1.0));
#ifdef DEBUG
        printf("mutation_rate: %lf\n", mutation_rate);
#endif
        if (mutation_rate == 0.0) {
                return;
        }
        bool *is_item_used = malloc(sizeof(*is_item_used) * num_items);
        for (size_t i=0; i<num_items; i++) {
                is_item_used[i] = true;
        }
        /* 'mutate' (delete) each bin with a probability specified by the
         * mutation rate */
        for (size_t i=0; i<chrom->num_bins; i++) {
#ifndef DEBUG
                if (rand() * (1 / (double)RAND_MAX) <= mutation_rate) {
#else
                double rand_mut = rand() * (1 / (double)RAND_MAX);
                printf("randomly generated value: %lf\n", rand_mut);
                if (rand_mut <= mutation_rate) {
                        printf("deleting bin %zu\n", i);
#endif
                        mark_unused(is_item_used, chrom->bins[i]);
                        chrom_del_bin(chrom, i);
                        i--; // have to account for bins shifting back
                }
        }
        /* first fit items from deleted bins */
        first_fit(chrom, item_sizes, is_item_used, num_items, 0);
        free(is_item_used);
        eval_fitness(chrom, FITNESS_K);
}

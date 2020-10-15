#include "population.h"
#include <stdio.h>
#include <stdlib.h>

#define POP_SZ          10
#define ARR_SZ          20
#define TEST_CAP        1000

static void print_bin(const bin_t *bin,
                      const long double *item_sizes, size_t num_items);
static void print_chrom(const chrom_t *chrom,
                        const long double *item_sizes, size_t num_items);
static void print_pop(const pop_t *pop,
                      const long double *item_sizes, size_t num_items);

int main(void) {
        srand(3);
        long double *arr = malloc(ARR_SZ * sizeof(*arr));
        printf("test set:\n");
        for (size_t i=0; i<ARR_SZ; i++) {
                arr[i] = rand() % (TEST_CAP / 2) + 1;
                printf("%Lf ", arr[i]);
        }
        putchar('\n');
        printf("rand pop\n");
        pop_t *pop = pop_rand_init(TEST_CAP, POP_SZ, arr, ARR_SZ);
        printf("pop:\n");
        print_pop(pop, arr, ARR_SZ);
        printf("free\n");
        pop_free(pop);
        free(arr);
        return 0;
}

static void print_bin(const bin_t *bin,
                      const long double *item_sizes, size_t num_items) {
        printf("fill: %Lf\n"
               "count: %zu\n"
               "items:\n",
               bin->fill, bin->count);
        for (size_t i=0; i<bin->count; i++) {
                printf("index: %zu\tsize: %Lf\n",
                       bin->item_indices[i],
                       item_sizes[bin->item_indices[i]]);
        }
}
static void print_chrom(const chrom_t *chrom,
                        const long double *item_sizes, size_t num_items) {
        printf("fitness: %lf\n"
               "bin_cap: %Lf\n"
               "num_bins: %zu\n"
               "bins:\n",
               chrom->fitness, chrom->bin_cap, chrom->num_bins);
        for (size_t i=0; i<chrom->num_bins; i++) {
                printf("bin %zu:\n", i);
                print_bin(chrom->bins[i], item_sizes, num_items);
        }
}
static void print_pop(const pop_t *pop,
                      const long double *item_sizes, size_t num_items) {
        printf("num_chroms: %zu\n"
               "chroms:\n",
               pop->num_chroms);
        for (size_t i=0; i<pop->num_chroms; i++) {
                printf("chrom %zu:\n", i);
                print_chrom(pop->chroms[i], item_sizes, num_items);
        }
        putchar('\n');
}

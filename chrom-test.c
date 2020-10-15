#include "chromosome.h"
#include <stdio.h>
#include <stdlib.h>

#define ARR_SZ          20
#define TEST_CAP        1000
#define MUT_RATE        (0.75)

static void print_bin(const bin_t *bin,
                      const long long *item_sizes, size_t num_items);
static void print_chrom(const chrom_t *chrom,
                        const long long *item_sizes, size_t num_items);

int main(void) {
        srand(3);
        long long *arr = malloc(ARR_SZ * sizeof(*arr));
        printf("test set:\n");
        for (size_t i=0; i<ARR_SZ; i++) {
                arr[i] = rand() % (TEST_CAP / 2) + 1;
                printf("%lld ", arr[i]);
        }
        putchar('\n');
        printf("rand\n");
        chrom_t *chrom = rand_first_fit(arr, ARR_SZ, TEST_CAP);
        printf("chrom 1:\n");
        print_chrom(chrom, arr, ARR_SZ);
        putchar('\n');
        printf("rand\n");
        chrom_t *chrom2 = rand_first_fit(arr, ARR_SZ, TEST_CAP);
        printf("chrom 2:\n");
        print_chrom(chrom2, arr, ARR_SZ);
        putchar('\n');
        printf("crossover\n");
        chrom_t *child = chrom_cx(chrom, chrom2, arr, ARR_SZ);
        printf("child chrom:\n");
        print_chrom(child, arr, ARR_SZ);
        putchar('\n');
        printf("mutate child\n");
        chrom_mutate(child, MUT_RATE, arr, ARR_SZ);
        printf("mutated chrom:\n");
        print_chrom(child, arr, ARR_SZ);
        putchar('\n');
        printf("copy chrom\n");
        chrom_t *copy = chrom_copy(child);
        printf("free chrom 1\n");
        chrom_free(chrom);
        printf("free chrom 2\n");
        chrom_free(chrom2);
        printf("free child chrom\n");
        chrom_free(child);
        printf("copied chrom:\n");
        print_chrom(copy, arr, ARR_SZ);
        putchar('\n');
        printf("free copy\n");
        chrom_free(copy);
        free(arr);
        return 0;
}

static void print_bin(const bin_t *bin,
                      const long long *item_sizes, size_t num_items) {
        printf("fill: %lld\n"
               "count: %zu\n"
               "items:\n",
               bin->fill, bin->count);
        for (size_t i=0; i<bin->count; i++) {
                printf("index: %zu\tsize: %lld\n",
                       bin->item_indices[i],
                       item_sizes[bin->item_indices[i]]);
        }
}
static void print_chrom(const chrom_t *chrom,
                        const long long *item_sizes, size_t num_items) {
        printf("fitness: %lf\n"
               "bin_cap: %zu\n"
               "num_bins: %zu\n"
               "bins:\n",
               chrom->fitness, chrom->bin_cap, chrom->num_bins);
        for (size_t i=0; i<chrom->num_bins; i++) {
                printf("bin %zu:\n", i);
                print_bin(chrom->bins[i], item_sizes, num_items);
        }
}

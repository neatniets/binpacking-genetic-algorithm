#include "bin-packing.h"
#include <stdio.h>
#include <stdlib.h>

#define NUM_PASSES      25

static void falk_main(void);

int main(void) {
        falk_main();
        return 0;
}

static void falk_main_solve(void) {
        /* skip problem identifier */
        scanf(" %*s");
        size_t bin_capacity, num_items, optimal_num_bins;
        scanf(" %zu %zu %zu",
              &bin_capacity, &num_items, &optimal_num_bins);
        long long *item_sizes = malloc(num_items
                                       * sizeof(*item_sizes));
        for (size_t i=0; i<num_items; i++) {
                scanf(" %lld", item_sizes+i);
        }
        prob_set_t ps = {.item_sizes = item_sizes,
                         .num_items = num_items,
                         .bin_capacity = bin_capacity,
                         .max_generations = 10,
                         .population_size = 100,
                         .mating_pool_size = 100,
                         .max_mutation_rate = 0.5,
                         .tournament_p = 1.0,
                         .tournament_size = 2,
                         .fitness_k = 2,
                         .use_adaptive_mutation = false};
        printf("OPTIMAL NUMBER OF BINS: %zu\n", optimal_num_bins);
        for (size_t i=0; i<NUM_PASSES; i++) {
                printf("PASS #%zu:\n", i);
                result_free(bin_packing(&ps));
        }
        free(item_sizes);
}
static void falk_main(void) {
        size_t num_problems;
        scanf(" %zu", &num_problems);
        for (size_t i=0; i<num_problems; i++) {
                printf("PROBLEM #%zu:\n", i);
                falk_main_solve();
        }
}

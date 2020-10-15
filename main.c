#include "bin-packing.h"
#include <stdio.h>
#include <stdlib.h>

#define ARR_SZ          1000
#define CAP             1000
#define MAX_GEN         1000
#define POP_SZ          100
#define MP_SZ           POP_SZ
#define MUT_RT          0.0
#define TOURN_P         1.0
#define TOURN_SZ        2
#define FIT_K           2

int main(void) {
        long long *arr = malloc(ARR_SZ * sizeof(*arr));
        for (size_t i=0; i<ARR_SZ; i++) {
                arr[i] = rand() % CAP + 1;
        }
        prob_set_t ps = {.item_sizes = arr,
                         .num_items = ARR_SZ,
                         .bin_capacity = CAP,
                         .max_generations = MAX_GEN,
                         .population_size = POP_SZ,
                         .mating_pool_size = MP_SZ,
                         .max_mutation_rate = MUT_RT,
                         .tournament_p = TOURN_P,
                         .tournament_size = TOURN_SZ,
                         .fitness_k = FIT_K,
                         .use_adaptive_mutation = false};
        result_free(bin_packing(&ps));
        free(arr);
        return 0;
}

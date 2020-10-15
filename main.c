#include "bin-packing.h"
#include <stdio.h>

const long long arr[] = {9, 8, 7, 6, 4, 3, 2, 1,
                         1, 2, 3, 4, 6, 7, 8, 9,
                         9, 8, 7, 6, 4, 3, 2, 1};
#define ARR_SZ          (sizeof(arr)/sizeof(*arr))
#define CAP             20
#define MAX_GEN         100
#define POP_SZ          ARR_SZ
#define MP_SZ           POP_SZ
#define MUT_RT          0.5
#define TOURN_P         1.0
#define TOURN_SZ        2
#define FIT_K           2

int main(void) {
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
        return 0;
}

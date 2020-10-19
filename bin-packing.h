#ifndef BIN_PACKING_H
#define BIN_PACKING_H

#include <stddef.h>
#include <stdbool.h>

typedef struct problem_set prob_set_t;
struct problem_set {
        const long double *item_sizes;
        size_t num_items;
        size_t bin_capacity;
        size_t max_generations;
        size_t terminal_num_bins;
        double max_secs;
        size_t population_size;
        size_t mating_pool_size;
        double max_mutation_rate;
        double tournament_p;
        unsigned tournament_size;
        bool use_inversion_operator;
};

struct llarray {
        size_t num_elems;
        long double elems[];
};
typedef struct result result_t;
struct result {
        double fitness;
        size_t num_bins;
        struct llarray *bins[];
};

void result_free(result_t *res);

result_t *bin_packing(const prob_set_t *ps);

#endif /* !BIN_PACKING_H */

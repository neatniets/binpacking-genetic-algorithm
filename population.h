#ifndef POPULATION_H
#define POPULATION_H

#include "chromosome.h"

typedef struct population pop_t;
struct population {
        const size_t num_chroms;
        chrom_t *chroms[];
};

pop_t *pop_alloc(size_t pop_size);
pop_t *pop_rand_init(long double bin_capacity, size_t pop_size,
                     const long double *item_sizes, size_t num_items);
void pop_free(pop_t *pop);

#endif /* !POPULATION_H */

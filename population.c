#include "population.h"
#include <stdlib.h>
#include <string.h>

pop_t *pop_alloc(size_t pop_size) {
        pop_t *pop = malloc(offsetof(pop_t, chroms)
                            + (pop_size * sizeof(*pop->chroms)));
        memcpy((size_t *)&pop->num_chroms, &pop_size,
               sizeof(pop->num_chroms));
        memset(&pop->chroms, 0, pop->num_chroms * sizeof(*pop->chroms));
        return pop;
}
pop_t *pop_rand_init(long double bin_capacity, size_t pop_size,
                     const long double *item_sizes, size_t num_items) {
        pop_t *pop = pop_alloc(pop_size);
        for (size_t i=0; i<pop->num_chroms; i++) {
                pop->chroms[i] = rand_first_fit(item_sizes, num_items,
                                                bin_capacity);
        }
        return pop;
}
void pop_free(pop_t *pop) {
        if (pop == NULL) {
                return;
        }
        for (size_t i=0; i<pop->num_chroms; i++) {
                chrom_free(pop->chroms[i]);
        }
        free(pop);
}

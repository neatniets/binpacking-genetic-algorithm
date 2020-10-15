#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include <stddef.h>

typedef struct bin bin_t;
struct bin {
        long double fill;
        size_t count;
        size_t *item_indices;
};

typedef struct chromosome chrom_t;
struct chromosome {
        double fitness;
        long double bin_cap;
        size_t num_bins;
        bin_t **bins;
};

chrom_t *rand_first_fit(const long double *item_sizes, size_t num_items,
                        long double bin_cap);
void chrom_free(chrom_t *chrom);

chrom_t *chrom_copy(const chrom_t *chrom);

chrom_t *chrom_cx(const chrom_t *parent1, const chrom_t *parent2,
                  const long double *item_sizes, size_t num_items);
void chrom_mutate(chrom_t *chrom, double mutation_rate,
                  const long double *item_sizes, size_t num_items);

#endif /* !CHROMOSOME_H */

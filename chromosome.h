#ifndef CHROMOSOME_H
#define CHROMOSOME_H

typedef struct bin bin_t;
struct bin {
        long long fill;
        size_t count;
        size_t *item_indices;
};

typedef struct chromosome chrom_t;
struct chromosome {
        double fitness;
        long long bin_cap;
        size_t num_bins;
        bin_t **bins;
};

chrom_t *rand_first_fit(const long long *item_sizes, size_t num_items);
void chrom_free(chrom_t *chrom);

chrom_t *chrom_cx(const chrom_t *parent1, const chrom_t *parent2,
                  const long long *item_sizes, size_t num_items);
void chrom_mutate(chrom_t *chrom);

#endif /* !CHROMOSOME_H */

#include "chromosome.h"
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define FITNESS_K       2

static bin_t *bin_alloc(void) {
        bin_t *bin = malloc(sizeof(*bin));
        *bin = (bin_t){.fill = 0,
                       .count = 0,
                       .item_indices = NULL};
        return bin;
}
static void bin_free(bin_t *bin) {
        free(bin->item_indices);
        free(bin);
}
static bin_t *bin_copy(const bin_t *bin) {
        bin_t *copy = bin_alloc();
        *copy = (bin_t){.fill = bin->fill,
                        .count = bin->count,
                        .item_indices = malloc(bin->count
                                               * sizeof(*bin->item_indices))
                       };
        memcpy(copy->item_indices,
               bin->item_indices,
               bin->count * sizeof(*copy->item_indices));
        return copy;
}
static void bin_add(bin_t *bin, size_t index, long long value) {
        bin->fill += value;
        bin->item_indices = realloc(bin->item_indices,
                                    (bin->count + 1)
                                    * sizeof(*bin->item_indices));
        bin->item_indices[bin->count] = index;
        bin->count++;
}

static chrom_t *chrom_alloc(size_t bin_cap) {
        chrom_t *chrom = malloc(sizeof(*chrom));
        *chrom = (chrom_t){.fitness = 0,
                           .bin_cap = bin_cap,
                           .num_bins = 0,
                           .bins = NULL};
        return chrom;
}
static void chrom_new_bin(chrom_t *chrom) {
        chrom->bins = realloc(chrom->bins,
                              (chrom->num_bins + 1)
                              * sizeof(*chrom->bins));
        chrom->bins[chrom->num_bins] = bin_alloc();
        chrom->num_bins++;
}
static void chrom_add_bin(chrom_t *chrom, bin_t *bin) {
        chrom_new_bin(chrom);
        bin_free(chrom->bins[chrom->num_bins - 1]);
        chrom->bins[chrom->num_bins - 1] = bin;
}
static void eval_fitness(chrom_t *chrom, int fitness_k) {
        chrom->fitness = 0;
        for (size_t i=0; i<chrom->num_bins; i++) {
                chrom->fitness += pow((double)chrom->bins[i]->fill
                                      / chrom->bin_cap,
                                      fitness_k)
                                  / chrom->num_bins;
        }
}

static void fit(chrom_t *chrom, size_t index, long long value) {
        for (size_t i=0; i<chrom->num_bins; i++) {
                if (chrom->bins[i]->fill + value <= chrom->bin_cap) {
                        bin_add(chrom->bins[i], index, value);
                        return;
                }
        }
        chrom_new_bin(chrom);
        bin_add(chrom->bins[chrom->num_bins - 1], index, value);
}
static void first_fit(chrom_t *chrom, const long long *item_sizes,
                      bool *is_item_used, size_t num_items,
                      size_t start_pos) {
        for (size_t loop = 0, i = start_pos, end = num_items;
             loop < 2;
             loop++, i = 0, end = start_pos) {
                if (!is_item_used[i]) {
                        fit(chrom, i, item_sizes[i]);
                        is_item_used[i] = true;
                }
        }
}
chrom_t *rand_first_fit(const long long *item_sizes, size_t num_items,
                        size_t bin_cap) {
        chrom_t *chrom = chrom_alloc(bin_cap);
        bool *is_item_used = calloc(num_items, sizeof(*is_item_used));
        first_fit(chrom, item_sizes, is_item_used, num_items,
                  rand() % num_items);
        free(is_item_used);
        eval_fitness(chrom, FITNESS_K);
        return chrom;
}
void chrom_free(chrom_t *chrom) {
        for (size_t i=0; i<chrom->num_bins; i++) {
                bin_free(chrom->bins[i]);
        }
        free(chrom->bins);
        free(chrom);
}

chrom_t *chrom_cx(const chrom_t *parent1, const chrom_t *parent2,
                  const long long *item_sizes, size_t num_items);
void chrom_mutate(chrom_t *chrom);

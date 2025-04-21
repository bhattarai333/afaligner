#ifndef HELPER_H
#define HELPER_H

#include <stddef.h>
#include <sys/types.h>

typedef struct D_matrix_element {
    double distance;
    ssize_t prev_i;
    ssize_t prev_j;
} D_matrix_element;

// Helper functions
double* get_coarsed_sequence(double* s, size_t n, size_t l);
size_t* get_window(size_t n, size_t m, size_t* path_buffer, size_t path_len, int radius);
void update_window(size_t* window, size_t n, size_t m, ssize_t i, ssize_t j);
double euclid_distance(double* x, double* y, size_t l);
D_matrix_element get_best_candidate(D_matrix_element* candidates, size_t n);
void reverse_path(size_t* path, ssize_t path_len);

#endif // HELPER_H
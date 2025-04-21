#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include "helper.h"
#include "logger.h"

sparse_matrix* create_sparse_matrix(size_t rows, size_t cols) {
    log_function_entry("create_sparse_matrix");
    sparse_matrix* matrix = malloc(sizeof(sparse_matrix));
    if (!matrix) {
        log_error("Failed to allocate sparse matrix");
        return NULL;
    }
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->elements = NULL;
    log_debug("Created sparse matrix with dimensions %zux%zu", rows, cols);
    log_function_exit("create_sparse_matrix", 0);
    return matrix;
}

void set_sparse_value(sparse_matrix* matrix, size_t i, size_t j, D_matrix_element value) {
    log_function_entry("set_sparse_value");
    sparse_element* elem = malloc(sizeof(sparse_element));
    if (!elem) {
        log_error("Failed to allocate sparse matrix element at (%zu, %zu)", i, j);
        return;
    }
    elem->i = i;
    elem->j = j;
    elem->value = value;
    elem->next = matrix->elements;
    matrix->elements = elem;
    log_debug("Set value at (%zu, %zu) with distance %f", i, j, value.distance);
    log_function_exit("set_sparse_value", 0);
}

D_matrix_element get_sparse_value(sparse_matrix* matrix, size_t i, size_t j) {
    log_function_entry("get_sparse_value");
    sparse_element* current = matrix->elements;
    while (current) {
        if (current->i == i && current->j == j) {
            log_debug("Found value at (%zu, %zu) with distance %f", i, j, current->value.distance);
            log_function_exit("get_sparse_value", current->value.distance);
            return current->value;
        }
        current = current->next;
    }
    D_matrix_element default_elem = {DBL_MAX, -1, -1};
    log_debug("No value found at (%zu, %zu), returning default", i, j);
    log_function_exit("get_sparse_value", DBL_MAX);
    return default_elem;
}

void free_sparse_matrix(sparse_matrix* matrix) {
    log_function_entry("free_sparse_matrix");
    sparse_element* current = matrix->elements;
    while (current) {
        sparse_element* next = current->next;
        free(current);
        current = next;
    }
    free(matrix);
    log_debug("Freed sparse matrix and all elements");
    log_function_exit("free_sparse_matrix", 0);
}

double *get_coarsed_sequence(double *s, size_t n, size_t l) {
    log_function_entry("get_coarsed_sequence");
    log_info("Creating coarsed sequence with n=%zu, l=%zu", n, l);

    size_t coarsed_sequence_len = n / 2;
    double *coarsed_sequence = malloc(coarsed_sequence_len * l * sizeof(double));

    if (coarsed_sequence == NULL) {
        log_error("Failed to allocate memory for coarsed sequence");
        return NULL;
    }

    log_debug("Computing average values for coarsed sequence");
    for (size_t i = 0; 2 * i + 1 < n ; i++) {
        for (size_t j = 0; j < l; j++) {
            coarsed_sequence[l*i+j] = (s[l*(2*i)+j] + s[l*(2*i+1)+j]) / 2;
        }
    }

    log_debug("Coarsed sequence created successfully");
    log_function_exit("get_coarsed_sequence", 0);
    return coarsed_sequence;
}

size_t *get_window(size_t n, size_t m, size_t *path_buffer, size_t path_len, int radius) {
    log_function_entry("get_window");
    log_info("Computing window with n=%zu, m=%zu, path_len=%zu, radius=%d", n, m, path_len, radius);

    size_t *window = malloc(2*n*sizeof(size_t));
    if (window == NULL) {
        log_error("Failed to allocate memory for window");
        return NULL;
    }

    log_debug("Initializing window boundaries");
    for (size_t i = 0; i < n; i++) {
        window[2*i] = m;
        window[2*i+1] = 0;
    }

    log_debug("Computing window constraints from path");
    for (size_t k = 0; k < path_len; k++) {
        size_t i = path_buffer[2*k];
        size_t j = path_buffer[2*k+1];

        for (ssize_t x = -radius; x < radius + 1; x++) {
            update_window(window, n, m, 2*(i + x), 2*(j - radius));
            update_window(window, n, m, 2*(i + x) + 1, 2*(j - radius));
            update_window(window, n, m, 2*(i + x), 2*(j + radius + 1) + 1);
            update_window(window, n, m, 2*(i + x) + 1, 2*(j + radius + 1) + 1);
        }
    }

    log_debug("Window computation completed");
    log_function_exit("get_window", 0);
    return window;
}

void update_window(size_t *window, size_t n, size_t m, ssize_t i, ssize_t j) {
    log_function_entry("update_window");
    log_debug("Updating window at i=%zd, j=%zd", i, j);

    if (i < 0 || i >= n) {
        log_debug("Index i out of bounds, skipping update");
        log_function_exit("update_window", 0);
        return;
    }

    if (j < 0) {
        log_debug("Clamping j to 0 (was %zd)", j);
        j = 0;
    }
    if (j > m - 1) {
        log_debug("Clamping j to %zu (was %zd)", m - 1, j);
        j = m - 1;
    }

    if (j < window[2*i]) {
        log_debug("Updating lower bound at i=%zd from %zu to %zd", i, window[2*i], j);
        window[2*i] = j;
    }
    if (j >= window[2*i+1]) {
        log_debug("Updating upper bound at i=%zd from %zu to %zd", i, window[2*i+1], j + 1);
        window[2*i+1] = j + 1;
    }

    log_function_exit("update_window", 0);
}

double euclid_distance(double *x, double *y, size_t l) {
    log_function_entry("euclid_distance");

    double sum = 0;
    for (size_t i = 0; i < l; i++) {
        double v = x[i] - y[i];
        sum += v * v;
    }
    double result = sqrt(sum);

    log_debug("Computed Euclidean distance: %f", result);
    log_function_exit("euclid_distance", result);
    return result;
}

double get_distance_array(D_matrix_element *D_matrix, size_t n, size_t m, size_t *window, size_t i, size_t j) {
    if (window != NULL) {
        if (i >= n || j >= m || j < window[2*i] || j >= window[2*i+1]) {
            return DBL_MAX;
        }
    } else if (i >= n || j >= m) {
        return DBL_MAX;
    }
    return D_matrix[i*m + j].distance;
}

// In dtwbd.c, replace the get_distance function with:
static double get_distance_sparse(sparse_matrix* matrix, size_t n, size_t m, size_t* window, size_t i, size_t j) {
    if (window != NULL) {
        if (i >= n || j >= m || j < window[2*i] || j >= window[2*i+1]) {
            log_function_exit("get_distance_sparse", DBL_MAX);
            return DBL_MAX;
        }
    } else if (i >= n || j >= m) {
        log_function_exit("get_distance_sparse", DBL_MAX);
        return DBL_MAX;
    }

    D_matrix_element elem = get_value(matrix, i, j);
    return elem.distance;
}



D_matrix_element get_best_candidate(D_matrix_element *candidates, size_t n) {
    log_function_entry("get_best_candidate");
    log_debug("Evaluating %zu candidates", n);

    double min_distance = DBL_MAX;
    D_matrix_element best_candidate;

    for (size_t i = 0; i < n; i++) {
        if (candidates[i].distance < min_distance) {
            min_distance = candidates[i].distance;
            best_candidate = candidates[i];
            log_debug("New best candidate found at index %zu with distance %f", i, min_distance);
        }
    }

    log_debug("Best candidate selected with distance %f", best_candidate.distance);
    log_function_exit("get_best_candidate", best_candidate.distance);
    return best_candidate;
}

void reverse_path(size_t *path, ssize_t path_len) {
    log_function_entry("reverse_path");
    log_info("Reversing path of length %zd", path_len);

    for (size_t i = 0, j = path_len - 1; i < j; i++, j--) {
        log_debug("Swapping elements at positions %zu and %zu", i, j);
        size_t tmp_s = path[2*i];
        size_t tmp_t = path[2*i+1];
        path[2*i] = path[2*j];
        path[2*i+1] = path[2*j+1];
        path[2*j] = tmp_s;
        path[2*j+1] = tmp_t;
    }

    log_debug("Path reversal completed");
    log_function_exit("reverse_path", path_len);
}


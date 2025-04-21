#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include "helper.h"

double *get_coarsed_sequence(double *s, size_t n, size_t l) {
    size_t coarsed_sequence_len = n / 2;
    double *coarsed_sequence = malloc(coarsed_sequence_len * l * sizeof(double));

    for (size_t i = 0; 2 * i + 1 < n ; i++) {
        for (size_t j = 0; j < l; j++) {
            coarsed_sequence[l*i+j] = (s[l*(2*i)+j] + s[l*(2*i+1)+j]) / 2;
        }
    }

    return coarsed_sequence;
}

size_t *get_window(size_t n, size_t m, size_t *path_buffer, size_t path_len, int radius) {
    size_t *window = malloc(2*n*sizeof(size_t));

    for (size_t i = 0; i < n; i++) {
        window[2*i] = m;
        window[2*i+1] = 0;
    }

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

    return window;
}

void update_window(size_t *window, size_t n, size_t m, ssize_t i, ssize_t j) {
    if (i < 0 ||

i >= n) return;

    if (j < 0) {
        j = 0;
    }
    if (j > m - 1) {
        j = m - 1;
    }
    if (j < window[2*i]) {
        window[2*i] = j;
    }
    if (j >= window[2*i+1]) {
        window[2*i+1] = j + 1;
    }
}

double euclid_distance(double *x, double *y, size_t l) {
    double sum = 0;
    for (size_t i = 0; i < l; i++) {
        double v = x[i] - y[i];
        sum += v * v;
    }
    return sqrt(sum);
}

double get_distance(D_matrix_element *D_matrix, size_t n, size_t m, size_t *window, size_t i, size_t j) {
    if (i < 0 || i >= n || j < 0 ||

j >= m) {
        return DBL_MAX;
    }

    if (window == NULL ||

(j >= window[2*i] && j < window[2*i+1])) {
        return D_matrix[i*m+j].distance;
    }

    return DBL_MAX;
}

D_matrix_element get_best_candidate(D_matrix_element *candidates, size_t n) {
    double min_distance = DBL_MAX;
    D_matrix_element best_candidate;

    for (size_t i = 0; i < n; i++) {
        if (candidates[i].distance < min_distance) {
            min_distance = candidates[i].distance;
            best_candidate = candidates[i];
        }
    }

    return best_candidate;
}

void reverse_path(size_t *path, ssize_t path_len) {
    for (size_t i = 0, j = path_len - 1; i < j; i++, j--) {
        size_t tmp_s = path[2*i];
        size_t tmp_t = path[2*i+1];
        path[2*i] = path[2*j];
        path[2*i+1] = path[2*j+1];
        path[2*j] = tmp_s;
        path[2*j+1] = tmp_t;
    }
}
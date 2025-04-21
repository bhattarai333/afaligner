#define _POSIX_C_SOURCE 200809L
#include <stdlib.h>
#include <stddef.h>
#include <sys/types.h>
#include <float.h>
#include <stdbool.h>
#include "helper.h"
#include "logger.h"

ssize_t DTWBD(double *s, double *t, size_t n, size_t m, size_t l,
              double skip_penalty, size_t *window, double *path_distance,
              size_t *path_buffer) {
    log_function_entry("DTWBD");
    log_info("Starting DTWBD with n=%zu, m=%zu, l=%zu, skip_penalty=%f", n, m, l, skip_penalty);

    // Allocate D_matrix
    D_matrix_element *D_matrix = malloc(n * m * sizeof(D_matrix_element));
    if (D_matrix == NULL) {
        log_error("Failed to allocate D_matrix");
        return -1;
    }

    // Initialize variables
    size_t end_i = 0;
    size_t end_j = 0;
    double min_path_distance = DBL_MAX;
    ssize_t path_len = 0;
    bool match = false;

    // Initialize D_matrix
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            if (window == NULL ||

(j >= window[2*i] && j < window[2*i+1])) {
                D_matrix[i*m + j].distance = euclid_distance(&s[i*l], &t[j*l], l);
                D_matrix[i*m + j].prev_i = -1;
                D_matrix[i*m + j].prev_j = -1;
            } else {
                D_matrix[i*m + j].distance = DBL_MAX;
                D_matrix[i*m + j].prev_i = -1;
                D_matrix[i*m + j].prev_j = -1;
            }
        }
    }

    // Dynamic programming
    for (size_t i = 0; i < n; i++) {
        for (size_t j = window ? window[2*i] : 0;
             j < (window ? window[2*i+1] : m); j++) {

            D_matrix_element candidates[3];
            size_t num_candidates = 0;

            // Regular step
            if (i > 0 && j > 0) {
                candidates[num_candidates].distance =
                    get_distance(D_matrix, n, m, window, i-1, j-1) +
                    D_matrix[i*m + j].distance;
                candidates[num_candidates].prev_i = i-1;
                candidates[num_candidates].prev_j = j-1;
                num_candidates++;
            }

            // Skip in s
            if (i > 0) {
                candidates[num_candidates].distance =
                    get_distance(D_matrix, n, m, window, i-1, j) +
                    skip_penalty;
                candidates[num_candidates].prev_i = i-1;
                candidates[num_candidates].prev_j = j;
                num_candidates++;
            }

            // Skip in t
            if (j > 0) {
                candidates[num_candidates].distance =
                    get_distance(D_matrix, n, m, window, i, j-1) +
                    skip_penalty;
                candidates[num_candidates].prev_i = i;
                candidates[num_candidates].prev_j = j-1;
                num_candidates++;
            }

            if (num_candidates > 0) {
                D_matrix_element best = get_best_candidate(candidates, num_candidates);
                D_matrix[i*m + j].distance += best.distance;
                D_matrix[i*m + j].prev_i = best.prev_i;
                D_matrix[i*m + j].prev_j = best.prev_j;

                if (i == n-1 && j == m-1) {
                    match = true;
                    end_i = i;
                    end_j = j;
                    min_path_distance = D_matrix[i*m + j].distance;
                }
            }
        }
    }

    if (match) {
        log_info("Path found with distance %f at (%zu, %zu)",
                 min_path_distance, end_i, end_j);

        // Reconstruct path
        size_t curr_i = end_i;
        size_t curr_j = end_j;
        path_len = 0;

        while (curr_i >= 0 && curr_j >= 0) {
            path_buffer[2*path_len] = curr_i;
            path_buffer[2*path_len+1] = curr_j;
            path_len++;

            ssize_t next_i = D_matrix[curr_i*m + curr_j].prev_i;
            ssize_t next_j = D_matrix[curr_i*m + curr_j].prev_j;

            if (next_i < 0 ||

next_j < 0) break;

            curr_i = next_i;
            curr_j = next_j;
        }

        reverse_path(path_buffer, path_len);
        *path_distance = min_path_distance;
    }

    free(D_matrix);
    log_function_exit("DTWBD", path_len);
    return path_len;
}
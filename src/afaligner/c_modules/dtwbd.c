#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <stdio.h>

#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;

__declspec(dllimport) size_t FastDTWBD();
__declspec(dllimport) size_t DTWBD();
#endif

typedef struct {
    double distance;
    ssize_t prev_i;
    ssize_t prev_j;
} D_matrix_element;

// Helper function to convert 2D coordinates to 1D index in banded matrix
static inline size_t get_band_index(size_t i, size_t j, size_t band_width, size_t j_offset) {
    return i * band_width + (j - j_offset);
}

// Helper function to check if a point is within the band
static inline bool is_in_band(size_t i, size_t j, size_t j_start, size_t band_width) {
    return (j >= j_start) && (j < j_start + band_width);
}

double euclid_distance(double *x, double *y, size_t l) {
    double sum = 0;
    for (size_t i = 0; i < l; i++) {
        double v = x[i] - y[i];
        sum += v * v;
    }
    return sqrt(sum);
}

void get_coarsed_sequence(double *sequence, size_t n, size_t l, size_t radius, double *coarsed) {
    for (size_t i = 0; i < n/2; i++) {
        for (size_t j = 0; j < l; j++) {
            coarsed[i*l + j] = (sequence[(2*i)*l + j] + sequence[(2*i+1)*l + j]) / 2;
        }
    }
}

void get_window(size_t *path, ssize_t path_len, size_t radius, size_t n, size_t m, size_t **window) {
    *window = malloc(sizeof(size_t) * n * 2);
    if (!*window) {
        return;
    }

    for (size_t i = 0; i < n; i++) {
        (*window)[2*i] = m;
        (*window)[2*i+1] = 0;
    }

    for (ssize_t i = 0; i < path_len; i++) {
        size_t x = path[2*i];
        size_t y = path[2*i+1];

        for (size_t j = (x > radius ? x - radius : 0);
             j <= (x + radius < n ? x + radius : n - 1); j++) {
            size_t min_y = y > radius ? y - radius : 0;
            size_t max_y = y + radius < m ? y + radius : m - 1;

            if (min_y < (*window)[2*j])
                (*window)[2*j] = min_y;
            if (max_y + 1 > (*window)[2*j+1])
                (*window)[2*j+1] = max_y + 1;
        }
    }
}

void update_window(size_t *window, size_t n, size_t m, size_t radius) {
    for (size_t i = 0; i < n; i++) {
        if (window[2*i+1] - window[2*i] > 2 * radius + 1) {
            size_t center = (window[2*i] + window[2*i+1]) / 2;
            window[2*i] = center > radius ? center - radius : 0;
            window[2*i+1] = center + radius + 1 < m ? center + radius + 1 : m;
        }
    }
}

D_matrix_element get_best_candidate(D_matrix_element *candidates, size_t n) {
    double min_distance = DBL_MAX;
    D_matrix_element best_candidate = {DBL_MAX, -1, -1};

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

double get_distance(D_matrix_element *D_matrix, size_t n, size_t band_width, size_t j_offset,
                   size_t i, size_t j) {
    if (i < 0 || i >= n ||

j < 0) {
        return DBL_MAX;
    }

    if (!is_in_band(i, j, j_offset, band_width)) {
        return DBL_MAX;
    }

    size_t band_idx = get_band_index(i, j, band_width, j_offset);
    return D_matrix[band_idx].distance;
}

ssize_t DTWBD(
    double *s, double *t,
    size_t n, size_t m, size_t l,
    double skip_penalty,
    size_t *window,
    double *path_distance,
    size_t *path_buffer
) {
    // Calculate band width based on window size
    size_t band_width;
    if (window) {
        size_t max_width = 0;
        for (size_t i = 0; i < n; i++) {
            size_t width = window[2*i+1] - window[2*i];
            if (width > max_width) {
                max_width = width;
            }
        }
        band_width = max_width;
    } else {
        band_width = m;
    }

    // Allocate banded matrix
    D_matrix_element *D_matrix = malloc(sizeof(D_matrix_element) * n * band_width);
    if (D_matrix == NULL) {
        fprintf(stderr, "ERROR: malloc() failed when allocating D_matrix\n");
        return -1;
    }

    // Initialize variables for minimum path tracking
    double min_last_row_distance = DBL_MAX;
    size_t min_last_row_j = 0;

    // Fill the matrix
    for (size_t i = 0; i < n; i++) {
        size_t j_start = window ? window[2*i] : 0;
        size_t j_end = window ? window[2*i+1] : m;

        for (size_t j = j_start; j < j_end; j++) {
            size_t band_idx = get_band_index(i, j, band_width, j_start);

            // Calculate distance between current points
            double point_distance = euclid_distance(&s[i*l], &t[j*l], l);

            // Initialize candidates array
            D_matrix_element candidates[3] = {
                {DBL_MAX, -1, -1},
                {DBL_MAX, -1, -1},
                {DBL_MAX, -1, -1}
            };

            // Get distances from previous cells
            if (i > 0 && j > 0) {
                double prev_distance = get_distance(D_matrix, n, band_width,
                                                 window ? window[2*(i-1)] : 0,
                                                 i-1, j-1);
                if (prev_distance != DBL_MAX) {
                    candidates[0].distance = prev_distance + point_distance;
                    candidates[0].prev_i = i-1;
                    candidates[0].prev_j = j-1;
                }
            }

            if (i > 0) {
                double prev_distance = get_distance(D_matrix, n, band_width,
                                                 window ? window[2*(i-1)] : 0,
                                                 i-1, j);
                if (prev_distance != DBL_MAX) {
                    candidates[1].distance = prev_distance + point_distance + skip_penalty;
                    candidates[1].prev_i = i-1;
                    candidates[1].prev_j = j;
                }
            }

            if (j > 0) {
                double prev_distance = get_distance(D_matrix, n, band_width,
                                                 j_start,
                                                 i, j-1);
                if (prev_distance != DBL_MAX) {
                    candidates[2].distance = prev_distance + point_distance + skip_penalty;
                    candidates[2].prev_i = i;
                    candidates[2].prev_j = j-1;
                }
            }

            // Special case for first cell
            if (i == 0 && j == 0) {
                candidates[0].distance = point_distance;
                candidates[0].prev_i = -1;
                candidates[0].prev_j = -1;
            }

            // Get best candidate
            D_matrix[band_idx] = get_best_candidate(candidates, 3);

            // Update minimum distance in last row
            if (i == n-1 && D_matrix[band_idx].distance < min_last_row_distance) {
                min_last_row_distance = D_matrix[band_idx].distance;
                min_last_row_j = j;
            }
        }
    }

    // Reconstruct path
    ssize_t path_len = 0;
    size_t curr_i = n-1;
    size_t curr_j = min_last_row_j;

    while (curr_i != (size_t)-1 && curr_j != (size_t)-1) {
        path_buffer[2*path_len] = curr_i;
        path_buffer[2*path_len+1] = curr_j;
        path_len++;

        size_t j_start = window ? window[2*curr_i] : 0;
        size_t band_idx = get_band_index(curr_i, curr_j, band_width, j_start);

        size_t next_i = D_matrix[band_idx].prev_i;
        size_t next_j = D_matrix[band_idx].prev_j;

        curr_i = next_i;
        curr_j = next_j;
    }

    // Reverse path since we reconstructed it backwards
    reverse_path(path_buffer, path_len);

    // Set return values
    *path_distance = min_last_row_distance;

    // Clean up
    free(D_matrix);

    return path_len;
}
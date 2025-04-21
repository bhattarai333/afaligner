#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "dtwbd.h"
#include "sparse_matrix.h"
#include "uthash.h"
#include <stdbool.h>
#include "fastdtwbd.h"
#include "logger.h"
#include <stdlib.h>


#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif


ssize_t DTWBD(
    double *s, size_t n,
    double *t, size_t m,
    size_t dim,
    double skip_penalty,
    size_t *window,
    size_t *path_buffer,
    double *path_distance
) {
    SparseMatrix *D_sparse = create_sparse_matrix();
    double d;

    // Logging start of DTWBD function
    log_info("Starting DTWBD function");

    // Fill the distance matrix D_sparse
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            if (window && (j < window[2 * i] || j >= window[2 * i + 1])) {
                continue;
            }

            d = euclid_distance(&s[i * dim], &t[j * dim], dim);

            double min_prev_distance = DBL_MAX;
            ssize_t prev_i = -1, prev_j = -1;

            // Check previous elements in the matrix
            D_matrix_element *e;
            e = get_element(D_sparse, i - 1, j);
            if (e && e->distance < min_prev_distance) {
                min_prev_distance = e->distance;
                prev_i = i - 1;
                prev_j = j;
            }

            e = get_element(D_sparse, i, j - 1);
            if (e && e->distance < min_prev_distance) {
                min_prev_distance = e->distance;
                prev_i = i;
                prev_j = j - 1;
            }

            e = get_element(D_sparse, i - 1, j - 1);
            if (e && e->distance < min_prev_distance) {
                min_prev_distance = e->distance;
                prev_i = i - 1;
                prev_j = j - 1;
            }

            // Update the current matrix element
            D_matrix_element cur;
            cur.distance = d + (min_prev_distance == DBL_MAX ? 0 : min_prev_distance);
            cur.prev_i = prev_i;
            cur.prev_j = prev_j;

            set_element(D_sparse, i, j, cur);

            // Log the current values for debugging
            log_debug("i: %zu, j: %zu, distance: %.4f, min_prev_distance: %.4f, prev_i: %zd, prev_j: %zd",
                      i, j, cur.distance, min_prev_distance, cur.prev_i, cur.prev_j);
        }
    }

    double min_path_distance = DBL_MAX;
    size_t end_i = 0, end_j = 0;
    bool match = false;

    // Find the minimum path distance
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            D_matrix_element *e = get_element(D_sparse, i, j);
            if (!e) continue;

            double cur_path_distance = e->distance + skip_penalty * (n - i + m - j - 2);

            // Log each comparison for path distance
            log_debug("Comparing path distance at i: %zu, j: %zu, cur_path_distance: %.4f",
                      i, j, cur_path_distance);

            if (cur_path_distance < min_path_distance) {
                min_path_distance = cur_path_distance;
                end_i = i;
                end_j = j;
                match = true;
            }
        }
    }

    ssize_t path_len = 0;

    if (match) {
        D_matrix_element *e;
        *path_distance = min_path_distance;

        // Log the minimum path distance and end points
        log_info("Found match. Min path distance: %.4f, end_i: %zu, end_j: %zu",
                 min_path_distance, end_i, end_j);

        for (ssize_t i = end_i, j = end_j; i != -1;) {
            e = get_element(D_sparse, i, j);
            if (!e) break;
            path_buffer[2 * path_len] = i;
            path_buffer[2 * path_len + 1] = j;
            path_len++;

            ssize_t next_i = e->prev_i;
            ssize_t next_j = e->prev_j;
            i = next_i;
            j = next_j;

            // Log each step of the path reconstruction
            log_debug("Reconstructing path: i: %zd, j: %zd", i, j);
        }

        reverse_path(path_buffer, path_len);

        // Log path reversal
        log_info("Path reconstruction complete, length: %zd", path_len);
    } else {
        log_info("No matching path found");
    }

    // Free the sparse matrix
    free_sparse_matrix(D_sparse);

    return path_len;
}

double euclid_distance(double *x, double *y, size_t l) {
    double sum = 0;

    // Log the input vectors and their length
    log_debug("Calculating Euclidean distance. Length: %zu", l);

    for (size_t i = 0; i < l; i++) {
        double v = x[i] - y[i];
        double squared_diff = v * v;
        sum += squared_diff;

        // Log the difference and squared difference for each dimension
        log_debug("i: %zu, x[i]: %.4f, y[i]: %.4f, v: %.4f, squared_diff: %.4f",
                  i, x[i], y[i], v, squared_diff);
    }

    double distance = sqrt(sum);

    // Log the final computed distance
    log_info("Euclidean distance calculated: %.4f", distance);

    return distance;
}

void reverse_path(size_t *path, ssize_t path_len) {
    // Log the original path
    log_debug("Reversing path. Path length: %zu", path_len);
    log_debug("Original path: ");
    for (ssize_t i = 0; i < path_len; i++) {
        log_debug("Path[%zu]: (%zu, %zu)", i, path[2*i], path[2*i + 1]);
    }

    for (size_t i = 0, j = path_len - 1; i < j; i++, j--) {
        // Log the swap details
        log_debug("Swapping indices i: %zu, j: %zu", i, j);
        log_debug("Before swap: path[%zu]: (%zu, %zu) <-> path[%zu]: (%zu, %zu)",
                  i, path[2*i], path[2*i+1], j, path[2*j], path[2*j+1]);

        // Perform the swap
        size_t tmp_s = path[2 * i];
        size_t tmp_t = path[2 * i + 1];
        path[2 * i] = path[2 * j];
        path[2 * i + 1] = path[2 * j + 1];
        path[2 * j] = tmp_s;
        path[2 * j + 1] = tmp_t;

        // Log the swap result
        log_debug("After swap: path[%zu]: (%zu, %zu) <-> path[%zu]: (%zu, %zu)",
                  i, path[2*i], path[2*i+1], j, path[2*j], path[2*j+1]);
    }

    // Log the final reversed path
    log_debug("Reversed path: ");
    for (ssize_t i = 0; i < path_len; i++) {
        log_debug("Path[%zu]: (%zu, %zu)", i, path[2*i], path[2*i + 1]);
    }
}

ssize_t FastDTWBD(
    double *s, double *t,
    size_t n, size_t m,
    size_t l,
    double skip_penalty,
    int radius,
    double *path_distance,
    size_t *path_buffer
) {
    ssize_t path_len;
    size_t min_sequence_len = 2 * (radius + 1) + 1;

    log_debug("Starting FastDTWBD with parameters: n=%zu, m=%zu, l=%zu, skip_penalty=%lf, radius=%d", n, m, l, skip_penalty, radius);

    // Base case
    if (n < min_sequence_len || m < min_sequence_len) {
        log_debug("Base case reached, calling DTWBD.");
        return DTWBD(s, n, t, m, l, skip_penalty, NULL, path_buffer, path_distance);
    }

    // Create coarsed sequences
    log_debug("Creating coarsed sequences for s and t.");
    double *coarsed_s = get_coarsed_sequence(s, n, l);
    if (!coarsed_s) {
        log_error("Failed to allocate coarsed_s.");
        return -1;  // Check allocation
    }

    double *coarsed_t = get_coarsed_sequence(t, m, l);
    if (!coarsed_t) {
        log_error("Failed to allocate coarsed_t.");
        free(coarsed_s);
        return -1;  // Check allocation
    }

    // Recursive call
    log_debug("Calling FastDTWBD recursively with coarsed sequences.");
    path_len = FastDTWBD(coarsed_s, coarsed_t, n / 2, m / 2, l, skip_penalty, radius, path_distance, path_buffer);

    if (path_len > 0) {
        log_debug("Path length from recursive call: %zd", path_len);

        // Create window and call DTWBD
        size_t *window = get_window(n, m, path_buffer, path_len, radius);
        if (window) {
            log_debug("Window created, calling DTWBD with the window.");
            path_len = DTWBD(s, n, t, m, l, skip_penalty, window, path_buffer, path_distance);
            free(window);
        } else {
            log_warn("Window creation failed.");
        }
    } else {
        log_warn("Recursive call returned an invalid path length.");
    }

    // Cleanup
    log_debug("Cleaning up, freeing coarsed sequences.");
    free(coarsed_s);
    free(coarsed_t);

    return path_len;
}


double *get_coarsed_sequence(double *s, size_t n, size_t l) {
    size_t coarsed_sequence_len = n / 2;

    log_debug("Allocating memory for coarsed sequence of length: %zu", coarsed_sequence_len);
    double *coarsed_sequence = malloc(coarsed_sequence_len * l * sizeof(double));

    if (!coarsed_sequence) {
        log_error("Memory allocation for coarsed sequence failed.");
        return NULL;
    }

    log_debug("Memory allocation successful for coarsed sequence.");

    // Create the coarsed sequence
    for (size_t i = 0; 2 * i + 1 < n; i++) {
        #log_debug("Creating coarsed sequence element %zu", i);

        for (size_t j = 0; j < l; j++) {
            coarsed_sequence[l * i + j] = (s[l * (2 * i) + j] + s[l * (2 * i + 1) + j]) / 2;
        }
    }

    log_debug("Coarsed sequence creation complete.");
    return coarsed_sequence;
}


size_t *get_window(size_t n, size_t m, size_t *path_buffer, size_t path_len, int radius) {
    log_debug("Allocating memory for window of size: %zu", 2 * n);
    size_t *window = malloc(2 * n * sizeof(size_t));

    if (!window) {
        log_error("Memory allocation for window failed.");
        return NULL;
    }

    log_debug("Memory allocation successful for window.");

    // Initialize window with max and min values
    for (size_t i = 0; i < n; i++) {
        window[2 * i] = m;    // maximum value for lower limit
        window[2 * i + 1] = 0;  // minimum value for upper limit
    }

    log_debug("Initialized window limits.");

    // Update window limits based on the path buffer
    for (size_t k = 0; k < path_len; k++) {
        size_t i = path_buffer[2 * k];
        size_t j = path_buffer[2 * k + 1];

        log_debug("Processing path element (%zu, %zu) at index %zu", i, j, k);

        for (ssize_t x = -radius; x < radius + 1; x++) {
            // update lower window limit
            update_window(window, n, m, 2 * (i + x), 2 * (j - radius));
            update_window(window, n, m, 2 * (i + x) + 1, 2 * (j - radius));

            // update upper window limit
            update_window(window, n, m, 2 * (i + x), 2 * (j + radius + 1) + 1);
            update_window(window, n, m, 2 * (i + x) + 1, 2 * (j + radius + 1) + 1);
        }
    }

    log_debug("Window computation complete, returning window.");

    return window;
}


void update_window(size_t *window, size_t n, size_t m, ssize_t i, ssize_t j) {
    log_debug("Updating window at index %zd with value %zd", i, j);

    if (i < 0 || i >= n) {
        log_debug("Index %zd is out of bounds (0 to %zu). No update performed.", i, n - 1);
        return;
    }

    if (j < 0) {
        log_debug("Index %zd is less than 0, adjusting to 0.", j);
        j = 0;
    }
    if (j > m - 1) {
        log_debug("Index %zd is greater than or equal to m (%zu), adjusting to %zu.", j, m, m - 1);
        j = m - 1;
    }

    if (j < window[2*i]) {
        log_debug("Updating lower limit of window at index %zd: %zu -> %zu", i, window[2*i], j);
        window[2*i] = j;
    }
    if (j >= window[2*i+1]) {
        log_debug("Updating upper limit of window at index %zd: %zu -> %zu", i, window[2*i+1], j + 1);
        window[2*i+1] = j + 1;
    }
}


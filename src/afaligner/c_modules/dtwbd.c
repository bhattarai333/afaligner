#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <stdio.h>
#include <time.h>

// Add this to export the functions
#define EXPORT __attribute__((visibility("default")))

typedef struct {
    double distance;
    ssize_t prev_i;
    ssize_t prev_j;
} D_matrix_element;

// Export the main functions
EXPORT ssize_t FastDTWBD(
    double *s, double *t,
    size_t n, size_t m, size_t l,
    double skip_penalty,
    size_t radius,
    double *path_distance,
    size_t *path_buffer
);

EXPORT ssize_t DTWBD(
    double *s, double *t,
    size_t n, size_t m, size_t l,
    double skip_penalty,
    size_t *window,
    double *path_distance,
    size_t *path_buffer
);


ssize_t FastDTWBD(
    double *s, double *t,
    size_t n, size_t m, size_t l,
    double skip_penalty,
    size_t radius,
    double *path_distance,
    size_t *path_buffer
);

// Global file pointer for logging
static FILE* log_fp = NULL;

// Initialize logging
static void init_logging() {
    if (log_fp == NULL) {
        log_fp = fopen("./output/afaligner.log", "a");
        if (log_fp == NULL) {
            fprintf(stderr, "Failed to open log file ./output/afaligner.log\n");
            return;
        }

        // Write header with timestamp
        time_t now = time(NULL);
        char timestamp[64];
        strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", localtime(&now));
        fprintf(log_fp, "\n=== New session started at %s ===\n", timestamp);
        fflush(log_fp);
    }
}

// Close logging
static void close_logging() {
    if (log_fp != NULL) {
        fclose(log_fp);
        log_fp = NULL;
    }
}

// Modify the logging macros
#define DEBUG_LOG(fmt, ...) do { \
    init_logging(); \
    if (log_fp) { \
        fprintf(log_fp, "[DEBUG] %s:%d: " fmt "\n", __func__, __LINE__, ##__VA_ARGS__); \
        fflush(log_fp); \
    } \
} while(0)

#define ERROR_LOG(fmt, ...) do { \
    init_logging(); \
    if (log_fp) { \
        fprintf(log_fp, "[ERROR] %s:%d: " fmt "\n", __func__, __LINE__, ##__VA_ARGS__); \
        fflush(log_fp); \
    } \
} while(0)

// Add cleanup function that should be called when the library is unloaded
__attribute__((destructor)) static void cleanup_logging() {
    close_logging();
}


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
    DEBUG_LOG("Window allocation size: %zu bytes", sizeof(size_t) * n * 2);
    if (!*window) {
        return;
    }

    // Add bounds checking for radius
    if (radius > n || radius > m) {
        ERROR_LOG("Radius %zu exceeds sequence dimensions n=%zu, m=%zu", radius, n, m);
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
    if (i < 0 || i >= n || j < 0) {
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
    DEBUG_LOG("Started DTWBD with n=%zu, m=%zu, l=%zu, skip_penalty=%f", n, m, l, skip_penalty);

    // Input validation
    if (!s || !t || !path_distance || !path_buffer) {
        ERROR_LOG("Null pointer passed as argument");
        return -1;
    }

    // Calculate band width based on window size
    size_t band_width;
    if (window) {
        size_t max_width = 0;
        for (size_t i = 0; i < n; i++) {
            size_t width = window[2*i+1] - window[2*i];
            DEBUG_LOG("Window[%zu]: start=%zu, end=%zu, width=%zu",
                     i, window[2*i], window[2*i+1], width);
            if (window[2*i+1] > m) {
                ERROR_LOG("Window end (%zu) exceeds sequence length (%zu)", window[2*i+1], m);
                return -1;
            }
            if (width > max_width) {
                max_width = width;
            }
        }
        band_width = max_width;
        DEBUG_LOG("Using windowed band_width=%zu", band_width);
    } else {
        band_width = m;
        DEBUG_LOG("Using full band_width=%zu", band_width);
    }

    // Allocate banded matrix
    DEBUG_LOG("Allocating D_matrix: size=%zu", sizeof(D_matrix_element) * n * band_width);
    D_matrix_element *D_matrix = malloc(sizeof(D_matrix_element) * n * band_width);
    if (D_matrix == NULL) {
        ERROR_LOG("malloc() failed when allocating D_matrix");
        return -1;
    }

    // Initialize variables for minimum path tracking
    double min_last_row_distance = DBL_MAX;
    size_t min_last_row_j = 0;

    // Fill the matrix
    DEBUG_LOG("Starting matrix fill");
    for (size_t i = 0; i < n; i++) {
        size_t j_start = window ? window[2*i] : 0;
        size_t j_end = window ? window[2*i+1] : m;

        DEBUG_LOG("Processing row %zu: j_start=%zu, j_end=%zu", i, j_start, j_end);

        if (j_end > m) {
            ERROR_LOG("j_end (%zu) exceeds sequence length (%zu)", j_end, m);
            free(D_matrix);
            return -1;
        }

        for (size_t j = j_start; j < j_end; j++) {
            size_t band_idx = get_band_index(i, j, band_width, j_start);

            if (band_idx >= n * band_width) {
                ERROR_LOG("band_idx (%zu) exceeds matrix size (%zu)", band_idx, n * band_width);
                free(D_matrix);
                return -1;
            }

            // Calculate distance between current points
            double point_distance = euclid_distance(&s[i*l], &t[j*l], l);
            DEBUG_LOG("Point distance at (%zu,%zu) = %f", i, j, point_distance);

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
            DEBUG_LOG("Best distance at (%zu,%zu) = %f", i, j, D_matrix[band_idx].distance);

            // Update minimum distance in last row
            if (i == n-1 && D_matrix[band_idx].distance < min_last_row_distance) {
                min_last_row_distance = D_matrix[band_idx].distance;
                min_last_row_j = j;
                DEBUG_LOG("New minimum in last row: j=%zu, distance=%f", j, min_last_row_distance);
            }
        }
    }

    // Reconstruct path
    DEBUG_LOG("Starting path reconstruction from i=%zu, j=%zu", n-1, min_last_row_j);
    ssize_t path_len = 0;
    size_t curr_i = n-1;
    size_t curr_j = min_last_row_j;

    while (curr_i != (size_t)-1 && curr_j != (size_t)-1) {
        if (path_len >= n + m) {
            ERROR_LOG("Path length exceeded maximum possible length");
            free(D_matrix);
            return -1;
        }

        path_buffer[2*path_len] = curr_i;
        path_buffer[2*path_len+1] = curr_j;
        path_len++;

        size_t j_start = window ? window[2*curr_i] : 0;
        size_t band_idx = get_band_index(curr_i, curr_j, band_width, j_start);

        if (band_idx >= n * band_width) {
            ERROR_LOG("Invalid band_idx during path reconstruction");
            free(D_matrix);
            return -1;
        }

        size_t next_i = D_matrix[band_idx].prev_i;
        size_t next_j = D_matrix[band_idx].prev_j;

        DEBUG_LOG("Path step %zd: (%zu,%zu) -> (%zu,%zu)",
                 path_len, curr_i, curr_j, next_i, next_j);

        curr_i = next_i;
        curr_j = next_j;
    }

    // Reverse path since we reconstructed it backwards
    DEBUG_LOG("Reversing path of length %zd", path_len);
    reverse_path(path_buffer, path_len);

    // Set return values
    *path_distance = min_last_row_distance;
    DEBUG_LOG("Final path distance: %f", min_last_row_distance);

    // Clean up
    free(D_matrix);
    DEBUG_LOG("DTWBD completed successfully");

    return path_len;
}

ssize_t FastDTWBD(
    double *s, double *t,
    size_t n, size_t m, size_t l,
    double skip_penalty,
    size_t radius,
    double *path_distance,
    size_t *path_buffer
) {
    DEBUG_LOG("=== FastDTWBD Start ===");
    DEBUG_LOG("Parameters: n=%zu, m=%zu, l=%zu, radius=%zu, skip_penalty=%f",
              n, m, l, radius, skip_penalty);

        // Add at the start of the function
    if (radius >= n || radius >= m) {
        ERROR_LOG("Radius %zu too large for sequence dimensions n=%zu, m=%zu", radius, n, m);
        radius = (n < m ? n : m) / 2;  // Set to half of smaller dimension
        DEBUG_LOG("Adjusted radius to %zu", radius);
    }

    // Input validation
    if (!s || !t) {
        ERROR_LOG("Null input sequences");
        return -1;
    }
    if (!path_distance || !path_buffer) {
        ERROR_LOG("Null output parameters");
        return -1;
    }
    if (n == 0 || m == 0 || l == 0) {
        ERROR_LOG("Zero-length input parameters");
        return -1;
    }

    // Base case for recursion
    if (n < 2 || m < 2) {
        DEBUG_LOG("Base case reached, calling DTWBD directly");
        return DTWBD(s, t, n, m, l, skip_penalty, NULL, path_distance, path_buffer);
    }

    // Calculate size of coarsed sequences
    size_t coarsed_n = n / 2;
    size_t coarsed_m = m / 2;
    DEBUG_LOG("Coarsed sequence sizes: coarsed_n=%zu, coarsed_m=%zu", coarsed_n, coarsed_m);

    // Memory allocation size checks
    size_t coarsed_s_size = sizeof(double) * coarsed_n * l;
    size_t coarsed_t_size = sizeof(double) * coarsed_m * l;
    size_t coarsed_path_size = sizeof(size_t) * (coarsed_n + coarsed_m) * 2;

    if (coarsed_s_size == 0 || coarsed_t_size == 0 || coarsed_path_size == 0) {
        ERROR_LOG("Invalid allocation sizes calculated");
        return -1;
    }

    DEBUG_LOG("Allocating coarsed sequences: s_size=%zu bytes, t_size=%zu bytes",
              coarsed_s_size, coarsed_t_size);

    // Allocate memory for coarsed sequences
    double *coarsed_s = malloc(coarsed_s_size);
    if (!coarsed_s) {
        ERROR_LOG("Failed to allocate coarsed_s");
        return -1;
    }

    double *coarsed_t = malloc(coarsed_t_size);
    if (!coarsed_t) {
        ERROR_LOG("Failed to allocate coarsed_t");
        free(coarsed_s);
        return -1;
    }

    DEBUG_LOG("Creating coarsed sequences");
    // Create coarsed sequences
    get_coarsed_sequence(s, n, l, radius, coarsed_s);
    get_coarsed_sequence(t, m, l, radius, coarsed_t);

    DEBUG_LOG("Allocating coarsed path of size %zu bytes", coarsed_path_size);
    // Allocate memory for coarsed path
    size_t *coarsed_path = malloc(coarsed_path_size);
    if (!coarsed_path) {
        ERROR_LOG("Failed to allocate coarsed_path");
        free(coarsed_s);
        free(coarsed_t);
        return -1;
    }

    // Recursive call
    DEBUG_LOG("Making recursive call with coarsed sequences");
    double coarsed_distance;
    ssize_t coarsed_path_len = FastDTWBD(
        coarsed_s, coarsed_t,
        coarsed_n, coarsed_m, l,
        skip_penalty, radius,
        &coarsed_distance,
        coarsed_path
    );

    // Clean up coarsed sequences early
    DEBUG_LOG("Freeing coarsed sequences");
    free(coarsed_s);
    free(coarsed_t);

    if (coarsed_path_len < 0) {
        ERROR_LOG("Recursive call failed with length %zd", coarsed_path_len);
        free(coarsed_path);
        return -1;
    }

    DEBUG_LOG("Recursive call succeeded, path_len=%zd", coarsed_path_len);

    // Calculate window from coarsed path
    DEBUG_LOG("Calculating window from coarsed path");
    size_t *window = NULL;
    get_window(coarsed_path, coarsed_path_len, radius, n, m, &window);

    // Clean up coarsed path
    free(coarsed_path);

    if (!window) {
        ERROR_LOG("Failed to create window");
        return -1;
    }

    // Update window to ensure minimum width
    DEBUG_LOG("Updating window bounds");
    update_window(window, n, m, radius);

    // Final DTW with window
    DEBUG_LOG("Performing final DTWBD with window");
    ssize_t path_len = DTWBD(s, t, n, m, l, skip_penalty, window, path_distance, path_buffer);

    // Clean up window
    DEBUG_LOG("Freeing window");
    free(window);

    if (path_len < 0) {
        ERROR_LOG("Final DTWBD call failed with length %zd", path_len);
        return -1;
    }

    DEBUG_LOG("FastDTWBD completed successfully, path_len=%zd", path_len);
    return path_len;
}
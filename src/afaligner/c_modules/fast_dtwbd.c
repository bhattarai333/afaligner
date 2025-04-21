#define _POSIX_C_SOURCE 200809L  // For ssize_t on some systems
#include <stdlib.h>
#include <stddef.h>
#include <sys/types.h>
#include "helper.h"
#include "logger.h"
#include <stdbool.h>

// Add export declaration
#ifdef _WIN32
    #define EXPORT __declspec(dllexport)
#else
    #define EXPORT __attribute__((visibility("default")))
#endif

// Add the EXPORT macro to the function declaration
EXPORT ssize_t FastDTWBD(double *s, double *t, size_t n, size_t m, size_t l,
                        double skip_penalty, int radius, double *path_distance,
                        size_t *path_buffer) {
    log_function_entry("FastDTWBD");
    log_info("Parameters: n=%zu, m=%zu, l=%zu, skip_penalty=%f, radius=%d",
             n, m, l, skip_penalty, radius);

    ssize_t path_len;
    size_t min_sequence_len = 2 * (radius + 1) + 1;

    if (n < min_sequence_len || m < min_sequence_len) {
        log_debug("Sequences too short, falling back to regular DTWBD");
        path_len = DTWBD(s, t, n, m, l, skip_penalty, NULL, path_distance, path_buffer);
        log_function_exit("FastDTWBD", path_len);
        return path_len;
    }

    double *coarsed_s = get_coarsed_sequence(s, n, l);
    if (!coarsed_s) {
        log_error("Failed to allocate memory for coarsed_s");
        return -1;
    }

    double *coarsed_t = get_coarsed_sequence(t, m, l);
    if (!coarsed_t) {
        log_error("Failed to allocate memory for coarsed_t");
        free(coarsed_s);
        return -1;
    }

    log_debug("Recursive FastDTWBD call with coarsed sequences");
    path_len = FastDTWBD(coarsed_s, coarsed_t, n/2, m/2, l, skip_penalty,
                         radius, path_distance, path_buffer);

    if (path_len < 0) {
        log_error("Recursive FastDTWBD call failed");
        free(coarsed_s);
        free(coarsed_t);
        return -1;
    }

    size_t *window = get_window(n, m, path_buffer, path_len, radius);
    if (!window) {
        log_error("Failed to allocate memory for window");
        free(coarsed_s);
        free(coarsed_t);
        return -1;
    }

    log_debug("Computing final DTWBD with window constraints");
    path_len = DTWBD(s, t, n, m, l, skip_penalty, window, path_distance, path_buffer);

    free(coarsed_s);
    free(coarsed_t);
    free(window);

    log_info("FastDTWBD completed with path_length=%zd, path_distance=%f",
             path_len, *path_distance);
    log_function_exit("FastDTWBD", path_len);
    return path_len;
}
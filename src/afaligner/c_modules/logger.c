#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <sys/stat.h>
#include "logger.h"
#include <stdbool.h>

static FILE* log_file = NULL;
static const char* level_strings[] = {
    "DEBUG",
    "INFO",
    "WARN",
    "ERROR"
};

static void ensure_logger_initialized() {
    if (!logger_initialized) {
        init_logger(DEFAULT_LOG_FILE);
        logger_initialized = true;
    }
}

static void create_directory(const char* path) {
    char tmp[256];
    char *p = NULL;
    size_t len;

    snprintf(tmp, sizeof(tmp), "%s", path);
    len = strlen(tmp);
    if (tmp[len - 1] == '/')
        tmp[len - 1] = 0;
    for (p = tmp + 1; *p; p++) {
        if (*p == '/') {
            *p = 0;
#ifdef _WIN32
            _mkdir(tmp);
#else
            mkdir(tmp, S_IRWXU);
#endif
            *p = '/';
        }
    }
#ifdef _WIN32
    _mkdir(tmp);
#else
    mkdir(tmp, S_IRWXU);
#endif
}

int init_logger(const char* log_path) {
    // Create output directory if it doesn't exist
    create_directory("./output");

    log_file = fopen(log_path, "a");
    if (!log_file) {
        fprintf(stderr, "Failed to open log file: %s\n", log_path);
        return -1;
    }

    // Write header to log file
    time_t now = time(NULL);
    fprintf(log_file, "\n\n=== New Session Started at %s===\n", ctime(&now));
    fflush(log_file);
    return 0;
}

void close_logger(void) {
    if (log_file) {
        time_t now = time(NULL);
        fprintf(log_file, "=== Session Ended at %s===\n", ctime(&now));
        fclose(log_file);
        log_file = NULL;
    }
}

static void log_message(LogLevel level, const char* format, va_list args) {
    ensure_logger_initialized();
    if (!log_file) return;

    time_t now = time(NULL);
    struct tm* timeinfo = localtime(&now);
    char timestamp[20];
    strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", timeinfo);

    fprintf(log_file, "[%s] [%s] ", timestamp, level_strings[level]);
    vfprintf(log_file, format, args);
    fprintf(log_file, "\n");
    fflush(log_file);
}

void log_debug(const char* format, ...) {
    ensure_logger_initialized();
    va_list args;
    va_start(args, format);
    log_message(LOG_DEBUG, format, args);
    va_end(args);
}

void log_info(const char* format, ...) {
    ensure_logger_initialized();
    va_list args;
    va_start(args, format);
    log_message(LOG_INFO, format, args);
    va_end(args);
}

void log_warn(const char* format, ...) {
    ensure_logger_initialized();
    va_list args;
    va_start(args, format);
    log_message(LOG_WARN, format, args);
    va_end(args);
}

void log_error(const char* format, ...) {
    ensure_logger_initialized();
    va_list args;
    va_start(args, format);
    log_message(LOG_ERROR, format, args);
    va_end(args);
}

void log_function_entry(const char* function_name) {
    ensure_logger_initialized();
    log_debug("Entering function: %s", function_name);
}

void log_function_exit(const char* function_name, ssize_t result) {
    ensure_logger_initialized();
    log_debug("Exiting function: %s (result: %zd)", function_name, result);
}
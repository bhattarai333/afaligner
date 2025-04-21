#include "logger.h"
#include <stdio.h>
#include <stdarg.h>
#include <time.h>

#define LOG_FILE_PATH "./output/afaligner.log"

static void log_with_prefix(const char *prefix, const char *format, va_list args) {
    FILE *log_file = fopen(LOG_FILE_PATH, "a");
    if (!log_file) return;

    time_t now = time(NULL);
    struct tm *t = localtime(&now);

    fprintf(log_file, "[%02d:%02d:%02d] [%s] ",
            t->tm_hour, t->tm_min, t->tm_sec, prefix);
    vfprintf(log_file, format, args);
    fprintf(log_file, "\n");

    fclose(log_file);
}

void log_info(const char *format, ...) {
#if LOG_INFO_ENABLED
    va_list args;
    va_start(args, format);
    log_with_prefix("INFO", format, args);
    va_end(args);
#endif
}

void log_debug(const char *format, ...) {
#if LOG_DEBUG_ENABLED
    va_list args;
    va_start(args, format);
    log_with_prefix("DEBUG", format, args);
    va_end(args);
#endif
}

void log_error(const char *format, ...) {
#if LOG_DEBUG_ENABLED
    va_list args;
    va_start(args, format);
    fprintf(stderr, "ERROR: ");
    vfprintf(stderr, format, args);
    va_end(args);
    fprintf(stderr, "\n");
#endif
}

void log_memory_alloc(const void *ptr, size_t size, const char *func, int line) {
#if LOG_MEM_ENABLED
    FILE *log_file = fopen(LOG_FILE_PATH, "a");
    if (!log_file) return;
    fprintf(log_file, "[MEMORY] Allocated %zu bytes at %p (%s:%d)\n", size, ptr, func, line);
    fclose(log_file);
#endif
}

void log_memory_free(const void *ptr, const char *func, int line) {
#if LOG_MEM_ENABLED
    FILE *log_file = fopen(LOG_FILE_PATH, "a");
    if (!log_file) return;
    fprintf(log_file, "[MEMORY] Freed memory at %p (%s:%d)\n", ptr, func, line);
    fclose(log_file);
#endif
}

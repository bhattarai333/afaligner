#ifndef LOGGER_H
#define LOGGER_H

#include <stddef.h>

// Enable or disable logging categories
#define LOG_INFO_ENABLED 1
#define LOG_MEM_ENABLED  1
#define LOG_DEBUG_ENABLED 1

#if defined(_WIN32) || defined(__WIN32__)
    #ifdef BUILDING_FASTDTWBD
        #define EXPORT __declspec(dllexport)
    #else
        #define EXPORT __declspec(dllimport)
    #endif
#else
    #ifdef BUILDING_FASTDTWBD
        #define EXPORT __attribute__((visibility("default")))
    #else
        #define EXPORT
    #endif
#endif

// Logging functions
EXPORT void log_info(const char *format, ...);
EXPORT void log_debug(const char *format, ...);
EXPORT void log_memory_alloc(const void *ptr, size_t size, const char *func, int line);
EXPORT void log_memory_free(const void *ptr, const char *func, int line);

// Helper macros for memory logging
#define LOG_ALLOC(ptr, size) log_memory_alloc(ptr, size, __func__, __LINE__)
#define LOG_FREE(ptr)        log_memory_free(ptr, __func__, __LINE__)

#endif // LOGGER_H

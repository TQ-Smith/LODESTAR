
#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
    char date_and_time[51];
    va_list args;
    time_t current_time;
    FILE* file;
} Logger;

static Logger* log = NULL;

void _init_log(const char* file);
void _log(const char* prefix, const char* fmt, ...);
void _close_log();

#define INIT_LOG(file) _init_log(file)
#define LOG_INFO(...) _log("[INFO]", __VA_ARGS__)
#define LOG_WARNING(...) _log("[WARNING]", __VA_ARGS__)
#define LOG_ERROR(...) _log("[ERROR]", __VA_ARGS__)
#define CLOSE_LOG() _close_log()

#endif

// Adapted from: https://tuttlem.github.io/2012/12/08/simple-logging-in-c.html

#include "Logger.h"

#include <pthread.h>

pthread_mutex_t logLock = PTHREAD_MUTEX_INITIALIZER;

void _init_log(const char* file_name) {
    pthread_mutex_lock(&logLock);
    log = (Logger*) malloc(sizeof(Logger));
    if (file_name == NULL)
        log -> file = stderr;
    else {
        log -> file = fopen(file_name, "w");
        if (log -> file == NULL) {
            log -> file = stderr;
            LOG_ERROR("Unable to open %s for logging! Using stderr ...\n", file_name);
        }
    }
    log -> current_time = time(NULL);
    tzset();
    pthread_mutex_unlock(&logLock);
}

void _log(const char* prefix, const char* fmt, ...) {
    pthread_mutex_lock(&logLock);
    if (log == NULL) {
        pthread_mutex_unlock(&logLock);
        return;
    }
    strftime(log -> date_and_time, sizeof(log -> date_and_time) - 1, "%a %b %d %T %Z %Y", localtime(&log -> current_time));
    fprintf(log -> file, "%s: %s ", log -> date_and_time, prefix);
    va_start(log -> args, fmt);
    vfprintf(log -> file, fmt, log -> args);
    va_end(log -> args);
    pthread_mutex_unlock(&logLock);
}

void _close_log() {
    pthread_mutex_lock(&logLock);
    if (log -> file == NULL) {
        pthread_mutex_unlock(&logLock);
        return;
    }
    if (log -> file != stderr)
        fclose(log -> file);
    free(log);
    log = NULL;
    pthread_mutex_unlock(&logLock);
}

/*
int main() {
    INIT_LOG("test.log");
    LOG_INFO("Test %d: Info message.\n", 1);
    LOG_WARNING("Test %d: Warning message.\n", 2);
    LOG_ERROR("Test %d: Error message.\n", 3);
    CLOSE_LOG();
}
*/

#ifndef _THREAD_POOL_
#define _THREAD_POOL_

#include <pthread.h>

#include <stdbool.h>

typedef struct ThreadPoolWork {

    void (*routine) ();
    void* arg;
    struct ThreadPoolWork* next;

} ThreadPoolWork_t;

typedef struct ThreadPool {

    int numThreads;
    int maxQueueSize;

    int doNotBlockWhenFull;

    pthread_t *threads;
    int curQueueSize;
    ThreadPoolWork_t* queueHead;
    ThreadPoolWork_t* queueTail;
    pthread_mutex_t queueLock;
    pthread_cond_t queueNotEmpty;
    pthread_cond_t queueNotFull;
    pthread_cond_t queueEmpty;

    bool queueClosed;
    bool shutdown;

} ThreadPool_t;

ThreadPool_t* init_thread_pool(int numWorkerThreads, int maxQueueSize, int doNotBlockWhenFull);

int thread_pool_add_work(ThreadPool_t* pool, void* routine, void* arg);

int thread_pool_destroy(ThreadPool_t* pool, bool finish);

#endif
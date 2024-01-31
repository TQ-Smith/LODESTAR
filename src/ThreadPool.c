
#include "ThreadPool.h"

#include <stdlib.h>

#include <stdio.h>

void* thread_loop(ThreadPool_t* pool) {

    ThreadPoolWork_t* work;

    while (true) {
        
        pthread_mutex_lock(&(pool -> queueLock));
        while((pool -> curQueueSize == 0) && !(pool -> shutdown))
            pthread_cond_wait(&(pool -> queueNotEmpty), &(pool -> queueLock));

        if (pool -> shutdown) {
            pthread_mutex_unlock(&(pool -> queueLock));
            pthread_exit(NULL);
        }

        work = pool -> queueHead;
        pool -> curQueueSize--;

        if (pool -> curQueueSize == 0)
            pool -> queueHead = pool -> queueTail = NULL;
        else
            pool -> queueHead = work -> next;

        if (!(pool -> doNotBlockWhenFull) && (pool -> curQueueSize == (pool -> maxQueueSize - 1)))
            pthread_cond_broadcast(&(pool -> queueNotFull));

        if (pool -> curQueueSize == 0)
            pthread_cond_signal(&(pool -> queueEmpty));

        pthread_mutex_unlock(&(pool -> queueLock));

        (*(work -> routine))(work -> arg);

        free(work);

    }

}

ThreadPool_t* init_thread_pool(int numWorkerThreads, int maxQueueSize, bool doNotBlockWhenFull) {

    ThreadPool_t* pool = (ThreadPool_t*) calloc(1, sizeof(ThreadPool_t));

    pool -> numThreads = numWorkerThreads;
    pool -> maxQueueSize = maxQueueSize;
    pool -> doNotBlockWhenFull = doNotBlockWhenFull;
    pool -> threads = (pthread_t*) calloc(numWorkerThreads, sizeof(pthread_t));
    pool -> curQueueSize = 0;
    pool -> queueHead = NULL;
    pool -> queueTail = NULL;
    pool -> queueClosed = false;
    pool -> shutdown = false;

    pthread_mutex_init(&(pool -> queueLock), NULL);
    pthread_cond_init(&(pool -> queueNotEmpty), NULL);
    pthread_cond_init(&(pool -> queueNotFull), NULL);
    pthread_cond_init(&(pool -> queueEmpty), NULL);

    for (int i = 0; i < numWorkerThreads; i++)
        pthread_create(&(pool -> threads[i]), NULL, (void * (*)(void *)) thread_loop, (void*) pool);

    return pool;

}

int thread_pool_add_work(ThreadPool_t* pool, void* routine, void* arg) {

    ThreadPoolWork_t* work;

    pthread_mutex_lock(&pool -> queueLock);

    if ((pool -> curQueueSize == pool -> maxQueueSize) && pool -> doNotBlockWhenFull) {
        pthread_mutex_unlock(&pool -> queueLock);
        return -1;
    }

    while ((pool -> curQueueSize == pool -> maxQueueSize) && (!(pool -> shutdown || pool -> queueClosed)))
        pthread_cond_wait(&pool -> queueNotFull, &pool -> queueLock);

    if (pool -> shutdown || pool -> queueClosed) {
        pthread_mutex_unlock(&pool -> queueLock);
        return -1;
    }

    work = (ThreadPoolWork_t*) calloc(1, sizeof(ThreadPool_t));
    work -> routine = routine;
    work -> arg = arg;
    work -> next = NULL;
    if (pool -> curQueueSize == 0) {
        pool -> queueTail = pool -> queueHead = work;
        pthread_cond_broadcast(&pool -> queueNotEmpty);
    } else {
        (pool -> queueTail) -> next = work;
        pool -> queueTail = work;
    }
    pool -> curQueueSize++;
    pthread_mutex_unlock(&pool -> queueLock);

    return 1;

}

void thread_pool_wait(ThreadPool_t* pool) {

    pthread_mutex_lock(&(pool -> queueLock));

    pthread_mutex_unlock(&(pool -> queueLock));

}

int thread_pool_destroy(ThreadPool_t* pool, bool finish) {

    ThreadPoolWork_t* temp;

    pthread_mutex_lock(&(pool -> queueLock));

    if (pool -> queueClosed || pool -> shutdown) {
        pthread_mutex_unlock(&(pool -> queueLock));
        return 0;
    }

    pool -> queueClosed = true;

    if (finish)
        while (pool -> curQueueSize != 0)
            pthread_cond_wait(&(pool -> queueEmpty), &(pool -> queueLock));
    
    pool -> shutdown = true;

    pthread_mutex_unlock(&(pool -> queueLock));
    pthread_cond_broadcast(&(pool -> queueNotEmpty));
    pthread_cond_broadcast(&(pool -> queueNotFull));

    for (int i = 0; i < pool -> numThreads; i++)
        pthread_join(pool -> threads[i], NULL);

    free(pool -> threads);
    while (pool -> queueHead != NULL) {
        temp = pool -> queueHead -> next;
        pool -> queueHead = pool -> queueHead -> next;
        free(temp);
    }
    free(pool);

    return 0;

}

void task(void *arg){
	printf("Thread #%u working on %d\n", (int) pthread_self(), (int) (long) arg);
}

int main() {

    ThreadPool_t* pool = init_thread_pool(4, 100, false);

    printf("\nAdding 50 tasks to threadpool:\n");
	for (int i = 0; i < 50; i++)
		thread_pool_add_work(pool, task, (void*) (long) i);

    thread_pool_destroy(pool, false);

    printf("Destroyed Thread Pool!\n");

}
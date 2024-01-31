
// File: ThreadPool.c
// Date: 31 January 2024
// Author: TQ Smith
// Purpose: A minimalistic thread pool adapted from Pthreads Programming by Nichols, Buttlar, Farrell.
//          Other resource used: https://nachtimwald.com/2019/04/12/thread-pool-in-c/

#include "ThreadPool.h"

#include <stdlib.h>

#include <stdio.h>

void* thread_loop(ThreadPool_t* pool) {

    ThreadPoolWork_t* work;

    while (true) {
        
        pthread_mutex_lock(&(pool -> queueLock));

        // Wait for something to be added to the queue.
        while((pool -> curQueueSize == 0) && !(pool -> shutdown))
            pthread_cond_wait(&(pool -> queueNotEmpty), &(pool -> queueLock));
        
        // If shutdown was indicated, release access and exit thread.
        if (pool -> shutdown) {
            pthread_mutex_unlock(&(pool -> queueLock));
            pthread_exit(NULL);
        }

        // Get the request at the head of the queue, decrement size of queue, and adjust head.
        work = pool -> queueHead;
        pool -> curQueueSize--;

        if (pool -> curQueueSize == 0)
            pool -> queueHead = pool -> queueTail = NULL;
        else
            pool -> queueHead = work -> next;
        
        // Set conditional if queue is full and the adding thread should block when full.
        if (!(pool -> doNotBlockWhenFull) && (pool -> curQueueSize == (pool -> maxQueueSize - 1)))
            pthread_cond_broadcast(&(pool -> queueNotFull));
        
        // Indicate the queue is empty.
        if (pool -> curQueueSize == 0)
            pthread_cond_signal(&(pool -> queueEmpty));
        
        // Release access to the queue.
        pthread_mutex_unlock(&(pool -> queueLock));

        // Execute function.
        (*(work -> routine))(work -> arg);

        // Deallocate the processed request.
        free(work);

    }

}

ThreadPool_t* init_thread_pool(int numWorkerThreads, int maxQueueSize, bool doNotBlockWhenFull) {

    // Allocate our thread pool structure.
    ThreadPool_t* pool = (ThreadPool_t*) calloc(1, sizeof(ThreadPool_t));

    // Set our default values.
    pool -> numThreads = numWorkerThreads;
    pool -> maxQueueSize = maxQueueSize;
    pool -> doNotBlockWhenFull = doNotBlockWhenFull;
    pool -> threads = (pthread_t*) calloc(numWorkerThreads, sizeof(pthread_t));
    pool -> curQueueSize = 0;
    pool -> queueHead = NULL;
    pool -> queueTail = NULL;
    pool -> queueClosed = false;
    pool -> shutdown = false;

    // Initalize our mutex and condition variables.
    pthread_mutex_init(&(pool -> queueLock), NULL);
    pthread_cond_init(&(pool -> queueNotEmpty), NULL);
    pthread_cond_init(&(pool -> queueNotFull), NULL);
    pthread_cond_init(&(pool -> queueEmpty), NULL);

    // Create our threads.
    for (int i = 0; i < numWorkerThreads; i++)
        pthread_create(&(pool -> threads[i]), NULL, (void * (*)(void *)) thread_loop, (void*) pool);

    // Return our created pool.
    return pool;

}

int thread_pool_add_work(ThreadPool_t* pool, void* routine, void* arg) {

    ThreadPoolWork_t* work;

    // Gain access to the pull.
    pthread_mutex_lock(&pool -> queueLock);

    // If queue is full and we aren't blocking on adding, then release access to
    //  queue and return -1.
    if ((pool -> curQueueSize == pool -> maxQueueSize) && pool -> doNotBlockWhenFull) {
        pthread_mutex_unlock(&pool -> queueLock);
        return -1;
    }

    // Block until there is room in the pool to add the request.
    while ((pool -> curQueueSize == pool -> maxQueueSize) && (!(pool -> shutdown || pool -> queueClosed)))
        pthread_cond_wait(&pool -> queueNotFull, &pool -> queueLock);
    
    // If the pool is being destroyed, release access and return -1.
    if (pool -> shutdown || pool -> queueClosed) {
        pthread_mutex_unlock(&pool -> queueLock);
        return -1;
    }

    // Create our work request and add our work to the queue.
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

    // Release access to the pool.
    pthread_mutex_unlock(&pool -> queueLock);

    // Successfully added request. Return success.
    return 0;

}

void thread_pool_wait(ThreadPool_t* pool) {

    // Gain access to the queue.
    pthread_mutex_lock(&(pool -> queueLock));

    // Wait for all the threads and request to finish processing.
    while (pool -> curQueueSize != 0)
        pthread_cond_wait(&(pool -> queueEmpty), &(pool -> queueLock));
    
    // Release access to the pool.
    pthread_mutex_unlock(&(pool -> queueLock));

}

void destroy_thread_pool(ThreadPool_t* pool, bool finish) {

    ThreadPoolWork_t* temp;

    // Gain access to the queue.
    pthread_mutex_lock(&(pool -> queueLock));

    // If pool is already shutdown, release access.
    if (pool -> queueClosed || pool -> shutdown)
        pthread_mutex_unlock(&(pool -> queueLock));
    
    // Restrict adding to the queue.
    pool -> queueClosed = true;

    // If flag is set, wait for queue to empty and threads to finish processing.
    if (finish)
        while (pool -> curQueueSize != 0)
            pthread_cond_wait(&(pool -> queueEmpty), &(pool -> queueLock));
    
    // Set flag that pool shutdown is occuring.
    pool -> shutdown = true;

    // Release access to the queue.
    pthread_mutex_unlock(&(pool -> queueLock));

    // Signal other thread to stop processing.
    pthread_cond_broadcast(&(pool -> queueNotEmpty));
    pthread_cond_broadcast(&(pool -> queueNotFull));

    // Join all threads together.
    for (int i = 0; i < pool -> numThreads; i++)
        pthread_join(pool -> threads[i], NULL);
    
    // Free threads.
    free(pool -> threads);
    
    // If there are still request in the queue (finish is false)
    //  then deallocate the queue.
    while (pool -> queueHead != NULL) {
        temp = pool -> queueHead -> next;
        pool -> queueHead = pool -> queueHead -> next;
        free(temp);
    }

    // Deallocae main structure.
    free(pool);

}

// Used to test the thread pool.

/*
#include <unistd.h>

void task(void *arg){
	printf("Thread #%u working on %d\n", (int) pthread_self(), (int) (long) arg);
    sleep(1);
}

int main() {

    ThreadPool_t* pool = init_thread_pool(4, 10, false);

    printf("\nAdding 50 tasks to threadpool:\n");
	for (int i = 0; i < 50; i++)
		thread_pool_add_work(pool, task, (void*) (long) i);

    thread_pool_wait(pool);

    destroy_thread_pool(pool, false);

    printf("Destroyed Thread Pool!\n");

}
*/
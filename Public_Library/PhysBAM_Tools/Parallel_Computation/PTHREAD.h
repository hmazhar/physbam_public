//#####################################################################
// Copyright 2011, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef _PTHREAD_H_
#define _PTHREAD_H_

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>

#include <pthread.h>
#ifdef __APPLE__

struct pthread_barrier_t
{
    pthread_mutex_t lock;
    pthread_cond_t cond;
    int threads;
    int blocked_threads;
};

inline int pthread_barrier_init(pthread_barrier_t* barrier,void* foo,int number_of_threads)
{
    int ret1=pthread_mutex_init(&barrier->lock,0);
    int ret2=pthread_cond_init(&barrier->cond,0);
    if(ret1 !=0 || ret2 !=0){PHYSBAM_FATAL_ERROR("failed to allocate mutex or condition in pthread_barrier_init");}
    barrier->threads=number_of_threads;
    barrier->blocked_threads=0;
    return 0;
}

inline int pthread_barrier_destroy(pthread_barrier_t* barrier)
{
    pthread_mutex_destroy(&barrier->lock);
    pthread_cond_destroy(&barrier->cond);
    return 0;
}

inline int pthread_barrier_wait (pthread_barrier_t *barrier)
{
    pthread_mutex_lock(&barrier->lock);
    barrier->blocked_threads++;
    if(barrier->blocked_threads==barrier->threads){
        barrier->blocked_threads=0;pthread_cond_broadcast(&barrier->cond);
    }
    else pthread_cond_wait(&barrier->cond,&barrier->lock);
    pthread_mutex_unlock(&barrier->lock);
    return 0;
}

#endif
#endif

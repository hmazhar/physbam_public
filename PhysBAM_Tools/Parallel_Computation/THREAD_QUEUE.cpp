//#####################################################################
// Copyright 2007, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
using namespace PhysBAM;

#ifdef USE_PTHREADS
//#####################################################################
// Constructor
//#####################################################################
THREAD_QUEUE::THREAD_QUEUE(const int thread_count,const bool set_affinity)
    :threads(thread_count),active_threads(thread_count),inactive_threads(0),queue(65535)
{
    pthread_attr_init(&attr);
    pthread_cond_init(&done_condition,0);
    pthread_cond_init(&todo_condition,0);
    pthread_mutex_init(&queue_lock,0);
    for(int i=1;i<=threads.m;i++){
        pthread_create(&threads(i),0,Thread_Routine,this);
        if(set_affinity){
            cpu_set_t cpuset;CPU_ZERO(&cpuset);CPU_SET(i-1,&cpuset);
            pthread_setaffinity_np(threads(i),sizeof(cpu_set_t),&cpuset);}}
}
//#####################################################################
// Destructor
//#####################################################################
THREAD_QUEUE::~THREAD_QUEUE()
{
    EXITER* exiter=new EXITER[threads.m];
    for(int i=1;i<=threads.m;i++) Queue(&exiter[i-1]);
    for(int i=1;i<=threads.m;i++) pthread_join(threads(i),NULL);
    
    pthread_cond_destroy(&done_condition);
    pthread_cond_destroy(&todo_condition);
    pthread_mutex_destroy(&queue_lock);
    pthread_attr_destroy(&attr);
    delete[] exiter;
}
//#####################################################################
// Queue
//#####################################################################
void THREAD_QUEUE::Queue(TASK* task)
{
    pthread_mutex_lock(&queue_lock);
    if(inactive_threads) pthread_cond_signal(&todo_condition);
    queue.Enqueue(task);
    pthread_mutex_unlock(&queue_lock);
}
//#####################################################################
// Wait
//#####################################################################
void THREAD_QUEUE::Wait()
{
    pthread_mutex_lock(&queue_lock);
    while(!queue.Empty() || active_threads!=0) pthread_cond_wait(&done_condition,&queue_lock);
    pthread_mutex_unlock(&queue_lock);
}
//#####################################################################
// Number_Of_Threads
//#####################################################################
int THREAD_QUEUE::Number_Of_Threads()
{
    return threads.Size();
}
#else
THREAD_QUEUE::THREAD_QUEUE(const int thread_count,const bool set_affinity) {PHYSBAM_FATAL_ERROR();}
THREAD_QUEUE::~THREAD_QUEUE() {PHYSBAM_FATAL_ERROR();}
void THREAD_QUEUE::Queue(TASK* task) {PHYSBAM_FATAL_ERROR();}
void THREAD_QUEUE::Wait() {PHYSBAM_FATAL_ERROR();}
int THREAD_QUEUE::Number_Of_Threads() {return 1;}
#endif

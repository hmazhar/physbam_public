//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INT_ITERATOR_THREADED
//#####################################################################
#ifndef __INT_ITERATOR_THREADED__
#define __INT_ITERATOR_THREADED__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/INTERVAL.h>
#include <PhysBAM_Tools/Parallel_Computation/ITERATOR_TASK.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
namespace PhysBAM{

/*template<class TYPE>
class INT_ITERATOR_THREADED
{
public:
    int row_jump;
    int start_index,end_index;
    THREAD_QUEUE* thread_queue;
    TYPE& threaded_class;

    INT_ITERATOR_THREADED(const int start_index_input,const int end_index_input,THREAD_QUEUE* thread_queue_input,TYPE& threaded_class_input)
        :row_jump(1),start_index(start_index_input),end_index(end_index_input),thread_queue(thread_queue_input),threaded_class(threaded_class_input)
    {
        if(thread_queue) row_jump=(end_index-start_index+1)/thread_queue->Number_Of_Threads()+1;
    }

    void Run()
    {
        if(!thread_queue){threaded_class.Run(start_index,end_index);return;}
#ifdef USE_PTHREADS
        int min_value=start_index,max_value=end_index;
        int local_start_index,local_end_index;
        for(int i=min_value;i<=max_value;i+=row_jump){
            local_start_index=min(i,max_value);local_end_index=min(i+row_jump-1,max_value);
            INT_ITERATOR_TASK<TYPE>* task=new INT_ITERATOR_TASK<TYPE>(threaded_class,local_start_index,local_end_index);
            thread_queue->Queue(task);}
        thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }
};*/

template<class TYPE>
class INT_ITERATOR_THREADED_ALPHA
{
public:
    int row_jump;
    int start_index,end_index;
    ARRAY<INTERVAL<int> > intervals;
    THREAD_QUEUE* thread_queue;

    INT_ITERATOR_THREADED_ALPHA(THREAD_QUEUE* thread_queue_input)
        :end_index(thread_queue_input->Number_Of_Threads()),thread_queue(thread_queue_input)
    {}

    INT_ITERATOR_THREADED_ALPHA(const int start_index_input,const int end_index_input,THREAD_QUEUE* thread_queue_input)
        :row_jump(1),start_index(start_index_input),end_index(end_index_input),thread_queue(thread_queue_input)
    {
        if(thread_queue){
            int number_of_intervals=thread_queue->Number_Of_Threads();
            row_jump=(end_index-start_index+1)/number_of_intervals+1;
            int min_value=start_index,max_value=end_index;
            int local_start_index,local_end_index;
            for(int i=min_value;i<=max_value;i+=row_jump){
                local_start_index=min(i,max_value);local_end_index=min(i+row_jump-1,max_value);
                intervals.Append(INTERVAL<int>(local_start_index,local_end_index));}}
        else intervals.Append(INTERVAL<int>(start_index,end_index));
    }

    void Run(TYPE& my_class,void (TYPE::*func)(int))
    {
        if(!thread_queue){(my_class.*func)(1);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=end_index;i++){
            INT_ITERATOR_TASK_0<TYPE>* task=new INT_ITERATOR_TASK_0<TYPE>(my_class,func,i);
            thread_queue->Queue(task);}
        thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    void Run(TYPE& my_class,void (TYPE::*func)(int,int))
    {
        if(!thread_queue){(my_class.*func)(start_index,end_index);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=intervals.m;i++){
            INT_ITERATOR_TASK_0<TYPE>* task=new INT_ITERATOR_TASK_0<TYPE>(my_class,func,intervals(i).min_corner,intervals(i).max_corner);
            thread_queue->Queue(task);}
        thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    void Run(TYPE& my_class,void (TYPE::*func)(int,int,int))
    {
        if(!thread_queue){(my_class.*func)(start_index,end_index,1);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=intervals.m;i++){
            INT_ITERATOR_TASK_0<TYPE>* task=new INT_ITERATOR_TASK_0<TYPE>(my_class,func,intervals(i).min_corner,intervals(i).max_corner,i);
            thread_queue->Queue(task);}
        thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1> void Run(TYPE& my_class,void (TYPE::*func)(T1,int),T1 arg1)
    {
        if(!thread_queue){(my_class.*func)(arg1,1);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=end_index;i++){
            INT_ITERATOR_TASK_1<TYPE,T1>* task=new INT_ITERATOR_TASK_1<TYPE,T1>(my_class,func,i,arg1);
            thread_queue->Queue(task);}
        thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1> void Run(TYPE& my_class,void (TYPE::*func)(T1,int,int),T1 arg1)
    {
        if(!thread_queue){(my_class.*func)(arg1,start_index,end_index);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=intervals.m;i++){
            INT_ITERATOR_TASK_1<TYPE,T1>* task=new INT_ITERATOR_TASK_1<TYPE,T1>(my_class,func,intervals(i).min_corner,intervals(i).max_corner,arg1);
            thread_queue->Queue(task);}
        thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1> void Run(TYPE& my_class,void (TYPE::*func)(T1,int,int,int),T1 arg1)
    {
        if(!thread_queue){(my_class.*func)(arg1,start_index,end_index,1);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=intervals.m;i++){
            INT_ITERATOR_TASK_1<TYPE,T1>* task=new INT_ITERATOR_TASK_1<TYPE,T1>(my_class,func,intervals(i).min_corner,intervals(i).max_corner,i,arg1);
            thread_queue->Queue(task);}
        thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1,class T2> void Run(TYPE& my_class,void (TYPE::*func)(T1,T2,int,int),T1 arg1,T2 arg2)
    {
        if(!thread_queue){(my_class.*func)(arg1,arg2,start_index,end_index);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=intervals.m;i++){
            INT_ITERATOR_TASK_2<TYPE,T1,T2>* task=new INT_ITERATOR_TASK_2<TYPE,T1,T2>(my_class,func,intervals(i).min_corner,intervals(i).max_corner,arg1,arg2);
            thread_queue->Queue(task);}
        thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1,class T2> void Run(TYPE& my_class,void (TYPE::*func)(T1,T2,int,int,int),T1 arg1,T2 arg2)
    {
        if(!thread_queue){(my_class.*func)(arg1,arg2,start_index,end_index,1);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=intervals.m;i++){
            INT_ITERATOR_TASK_2<TYPE,T1,T2>* task=new INT_ITERATOR_TASK_2<TYPE,T1,T2>(my_class,func,intervals(i).min_corner,intervals(i).max_corner,i,arg1,arg2);
            thread_queue->Queue(task);}
        thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1,class T2,class T3> void Run(TYPE& my_class,void (TYPE::*func)(T1,T2,T3,int,int),T1 arg1,T2 arg2,T3 arg3)
    {
        if(!thread_queue){(my_class.*func)(arg1,arg2,arg3,start_index,end_index);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=intervals.m;i++){
            INT_ITERATOR_TASK_3<TYPE,T1,T2,T3>* task=new INT_ITERATOR_TASK_3<TYPE,T1,T2,T3>(my_class,func,intervals(i).min_corner,intervals(i).max_corner,arg1,arg2,arg3);
            thread_queue->Queue(task);}
        thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1,class T2,class T3> void Run(TYPE& my_class,void (TYPE::*func)(T1,T2,T3,int,int,int),T1 arg1,T2 arg2,T3 arg3)
    {
        if(!thread_queue){(my_class.*func)(arg1,arg2,arg3,start_index,end_index,1);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=intervals.m;i++){
            INT_ITERATOR_TASK_3<TYPE,T1,T2,T3>* task=new INT_ITERATOR_TASK_3<TYPE,T1,T2,T3>(my_class,func,intervals(i).min_corner,intervals(i).max_corner,i,arg1,arg2,arg3);
            thread_queue->Queue(task);}
        thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1,class T2,class T3,class T4> void Run(TYPE& my_class,void (TYPE::*func)(T1,T2,T3,T4,int,int),T1 arg1,T2 arg2,T3 arg3,T4 arg4)
    {
        if(!thread_queue){(my_class.*func)(arg1,arg2,arg3,arg4,start_index,end_index);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=intervals.m;i++){
            INT_ITERATOR_TASK_4<TYPE,T1,T2,T3,T4>* task=new INT_ITERATOR_TASK_4<TYPE,T1,T2,T3,T4>(my_class,func,intervals(i).min_corner,intervals(i).max_corner,arg1,arg2,arg3,arg4);
            thread_queue->Queue(task);}
        thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1,class T2,class T3,class T4,class T5> void Run(TYPE& my_class,void (TYPE::*func)(T1,T2,T3,T4,T5,int,int),T1 arg1,T2 arg2,T3 arg3,T4 arg4,T5 arg5)
    {
        if(!thread_queue){(my_class.*func)(arg1,arg2,arg3,arg4,arg5,start_index,end_index);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=intervals.m;i++){
            INT_ITERATOR_TASK_5<TYPE,T1,T2,T3,T4,T5>* task=new INT_ITERATOR_TASK_5<TYPE,T1,T2,T3,T4,T5>(my_class,func,intervals(i).min_corner,intervals(i).max_corner,arg1,arg2,arg3,arg4,arg5);
            thread_queue->Queue(task);}
        thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }
};
}
#endif

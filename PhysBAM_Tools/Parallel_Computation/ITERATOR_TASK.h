//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ITERATOR_TASK
//#####################################################################
#ifdef USE_PTHREADS
#ifndef __ITERATOR_TASK__
#define __ITERATOR_TASK__

#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
namespace PhysBAM{

template<class TYPE,class T_ITERATOR,class T>
class ITERATOR_TASK:public THREAD_QUEUE::TASK
{   
public:
    TYPE& threaded_task;
    T_ITERATOR iterator;
    T dt,time;

    ITERATOR_TASK(TYPE& threaded_task_input,const T_ITERATOR& iterator_input,const T dt_input,const T time_input)
        :threaded_task(threaded_task_input),iterator(iterator_input),dt(dt_input),time(time_input)
    {}
    
    void Run()
    {
        threaded_task.Run(iterator,dt,time);
    }
};

template<class TYPE>
class INT_ITERATOR_TASK:public THREAD_QUEUE::TASK
{   
public:
    TYPE& threaded_task;
    int start_index,end_index,tid;

    INT_ITERATOR_TASK(TYPE& threaded_task_input,const int start_index_input,const int end_index_input,const int tid_input)
        :threaded_task(threaded_task_input),start_index(start_index_input),end_index(end_index_input),tid(tid_input)
    {}
    
    void Run()
    {
        threaded_task.Run(start_index,end_index,tid);
    }
};

template<class TYPE>
class INT_ITERATOR_TASK_0:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*tidonlyfunc)(int);
    void (TYPE::*func)(int,int);
    void (TYPE::*tidfunc)(int,int,int);
    int start_index,end_index,tid;

    INT_ITERATOR_TASK_0(TYPE& class_input,void (TYPE::*func_input)(int),const int tid_input)
        :my_class(class_input),tidonlyfunc(func_input),func(0),tidfunc(0),tid(tid_input)
    {}

    INT_ITERATOR_TASK_0(TYPE& class_input,void (TYPE::*func_input)(int,int),const int start_index_input,const int end_index_input)
        :my_class(class_input),tidonlyfunc(0),func(func_input),tidfunc(0),start_index(start_index_input),end_index(end_index_input)
    {}

    INT_ITERATOR_TASK_0(TYPE& class_input,void (TYPE::*func_input)(int,int,int),const int start_index_input,const int end_index_input,const int tid_input)
        :my_class(class_input),tidonlyfunc(0),func(0),tidfunc(func_input),start_index(start_index_input),end_index(end_index_input),tid(tid_input)
    {}
 
    void Run()
    {
        if(tidonlyfunc) (my_class.*tidonlyfunc)(tid);
        else if(func) (my_class.*func)(start_index,end_index);
        else (my_class.*tidfunc)(start_index,end_index,tid);
    }
};

template<class TYPE,class T1>
class INT_ITERATOR_TASK_1:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*tidonlyfunc)(T1,int);
    void (TYPE::*func)(T1,int,int);
    void (TYPE::*tidfunc)(T1,int,int,int);
    int start_index,end_index,tid;
    T1 arg1;

    INT_ITERATOR_TASK_1(TYPE& class_input,void (TYPE::*func_input)(T1,int),const int tid_input,T1 arg1_input)
        :my_class(class_input),tidonlyfunc(func_input),func(0),tidfunc(0),tid(tid_input),arg1(arg1_input)
    {}

    INT_ITERATOR_TASK_1(TYPE& class_input,void (TYPE::*func_input)(T1,int,int),const int start_index_input,const int end_index_input,T1 arg1_input)
        :my_class(class_input),tidonlyfunc(0),func(func_input),tidfunc(0),start_index(start_index_input),end_index(end_index_input),arg1(arg1_input)
    {}

    INT_ITERATOR_TASK_1(TYPE& class_input,void (TYPE::*func_input)(T1,int,int,int),const int start_index_input,const int end_index_input,const int tid_input,T1 arg1_input)
        :my_class(class_input),tidonlyfunc(0),func(0),tidfunc(func_input),start_index(start_index_input),end_index(end_index_input),tid(tid_input),arg1(arg1_input)
    {}
 
    void Run()
    {
        if(tidonlyfunc) (my_class.*tidonlyfunc)(arg1,tid);
        else if(func) (my_class.*func)(arg1,start_index,end_index);
        else (my_class.*tidfunc)(arg1,start_index,end_index,tid);
    }
};

template<class TYPE,class T1,class T2>
class INT_ITERATOR_TASK_2:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T1,T2,int,int);
    void (TYPE::*tidfunc)(T1,T2,int,int,int);
    int start_index,end_index,tid;
    T1 arg1;T2 arg2;

    INT_ITERATOR_TASK_2(TYPE& class_input,void (TYPE::*func_input)(T1,T2,int,int),const int start_index_input,const int end_index_input,T1 arg1_input,T2 arg2_input)
        :my_class(class_input),func(func_input),tidfunc(0),start_index(start_index_input),end_index(end_index_input),arg1(arg1_input),arg2(arg2_input)
    {}

    INT_ITERATOR_TASK_2(TYPE& class_input,void (TYPE::*func_input)(T1,T2,int,int,int),const int start_index_input,const int end_index_input,const int tid_input,T1 arg1_input,T2 arg2_input)
        :my_class(class_input),func(0),tidfunc(func_input),start_index(start_index_input),end_index(end_index_input),tid(tid_input),arg1(arg1_input),arg2(arg2_input)
    {}
    
    void Run()
    {
        if(func) (my_class.*func)(arg1,arg2,start_index,end_index);
        else (my_class.*tidfunc)(arg1,arg2,start_index,end_index,tid);
    }
};

template<class TYPE,class T1,class T2,class T3>
class INT_ITERATOR_TASK_3:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T1,T2,T3,int,int);
    void (TYPE::*tidfunc)(T1,T2,T3,int,int,int);
    int start_index,end_index,tid;
    T1 arg1;T2 arg2;T3 arg3;

    INT_ITERATOR_TASK_3(TYPE& class_input,void (TYPE::*func_input)(T1,T2,T3,int,int),const int start_index_input,const int end_index_input,T1 arg1_input,T2 arg2_input,T3 arg3_input)
        :my_class(class_input),func(func_input),tidfunc(0),start_index(start_index_input),end_index(end_index_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input)
    {}

    INT_ITERATOR_TASK_3(TYPE& class_input,void (TYPE::*func_input)(T1,T2,T3,int,int,int),const int start_index_input,const int end_index_input,const int tid_input,T1 arg1_input,T2 arg2_input,T3 arg3_input)
        :my_class(class_input),func(0),tidfunc(func_input),start_index(start_index_input),end_index(end_index_input),tid(tid_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input)
    {}
    
    void Run()
    {
        if(func) (my_class.*func)(arg1,arg2,arg3,start_index,end_index);
        else (my_class.*tidfunc)(arg1,arg2,arg3,start_index,end_index,tid);
    }
};

template<class TYPE,class T1,class T2,class T3,class T4>
class INT_ITERATOR_TASK_4:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T1,T2,T3,T4,int,int);
    int start_index,end_index;
    T1 arg1;T2 arg2;T3 arg3;T4 arg4;

    INT_ITERATOR_TASK_4(TYPE& class_input,void (TYPE::*func_input)(T1,T2,T3,T4,int,int),const int start_index_input,const int end_index_input,T1 arg1_input,T2 arg2_input,T3 arg3_input,T4 arg4_input)
        :my_class(class_input),func(func_input),start_index(start_index_input),end_index(end_index_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input),arg4(arg4_input)
    {}
    
    void Run()
    {
        (my_class.*func)(arg1,arg2,arg3,arg4,start_index,end_index);
    }
};

template<class TYPE,class T1,class T2,class T3,class T4,class T5>
class INT_ITERATOR_TASK_5:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T1,T2,T3,T4,T5,int,int);
    int start_index,end_index;
    T1 arg1;T2 arg2;T3 arg3;T4 arg4;T5 arg5;

    INT_ITERATOR_TASK_5(TYPE& class_input,void (TYPE::*func_input)(T1,T2,T3,T4,T5,int,int),const int start_index_input,const int end_index_input,T1 arg1_input,T2 arg2_input,T3 arg3_input,T4 arg4_input,T5 arg5_input)
        :my_class(class_input),func(func_input),start_index(start_index_input),end_index(end_index_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input),arg4(arg4_input),arg5(arg5_input)
    {}
    
    void Run()
    {
        (my_class.*func)(arg1,arg2,arg3,arg4,arg5,start_index,end_index);
    }
};

template<class TYPE,class T1,class T2,class T3,class T4,class T5,class T6>
class INT_ITERATOR_TASK_6:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T1,T2,T3,T4,T5,T6,int,int);
    int start_index,end_index;
    T1 arg1;T2 arg2;T3 arg3;T4 arg4;T5 arg5;T6 arg6;;

    INT_ITERATOR_TASK_6(TYPE& class_input,void (TYPE::*func_input)(T1,T2,T3,T4,T5,T6,int,int),const int start_index_input,const int end_index_input,T1 arg1_input,T2 arg2_input,T3 arg3_input,T4 arg4_input,T5 arg5_input,T6 arg6_input)
        :my_class(class_input),func(func_input),start_index(start_index_input),end_index(end_index_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input),arg4(arg4_input),arg5(arg5_input),arg6(arg6_input)
    {}
    
    void Run()
    {
        (my_class.*func)(arg1,arg2,arg3,arg4,arg5,arg6,start_index,end_index);
    }
};

template<class TYPE,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
class INT_ITERATOR_TASK_7:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T1,T2,T3,T4,T5,T6,T7,int,int);
    int start_index,end_index;
    T1 arg1;T2 arg2;T3 arg3;T4 arg4;T5 arg5;T6 arg6;T7 arg7;

    INT_ITERATOR_TASK_7(TYPE& class_input,void (TYPE::*func_input)(T1,T2,T3,T4,T5,T6,T7,int,int),const int start_index_input,const int end_index_input,T1 arg1_input,T2 arg2_input,T3 arg3_input,T4 arg4_input,T5 arg5_input,T6 arg6_input,T7 arg7_input)
        :my_class(class_input),func(func_input),start_index(start_index_input),end_index(end_index_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input),arg4(arg4_input),arg5(arg5_input),arg6(arg6_input),arg7(arg7_input)
    {}
    
    void Run()
    {
        (my_class.*func)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,start_index,end_index);
    }
};

template<class TYPE,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
class INT_ITERATOR_TASK_8:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T1,T2,T3,T4,T5,T6,T7,T8,int,int);
    int start_index,end_index;
    T1 arg1;T2 arg2;T3 arg3;T4 arg4;T5 arg5;T6 arg6;T7 arg7;T8 arg8;

    INT_ITERATOR_TASK_8(TYPE& class_input,void (TYPE::*func_input)(T1,T2,T3,T4,T5,T6,T7,T8,int,int),const int start_index_input,const int end_index_input,T1 arg1_input,T2 arg2_input,T3 arg3_input,T4 arg4_input,T5 arg5_input,T6 arg6_input,T7 arg7_input,T8 arg8_input)
        :my_class(class_input),func(func_input),start_index(start_index_input),end_index(end_index_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input),arg4(arg4_input),arg5(arg5_input),arg6(arg6_input),arg7(arg7_input),arg8(arg8_input)
    {}
    
    void Run()
    {
        (my_class.*func)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,start_index,end_index);
    }
};

template<class TYPE,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8, class T9, class T10, class T11, class T12, class T13>
class INT_ITERATOR_TASK_13:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,int,int);
    int start_index,end_index;
    T1 arg1;T2 arg2;T3 arg3;T4 arg4;T5 arg5;T6 arg6;T7 arg7;T8 arg8;T9 arg9;T10 arg10;T11 arg11;T12 arg12;T13 arg13;

    INT_ITERATOR_TASK_13(TYPE& class_input,void (TYPE::*func_input)(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,int,int),const int start_index_input,const int end_index_input,T1 arg1_input,T2 arg2_input,T3 arg3_input,T4 arg4_input,T5 arg5_input,T6 arg6_input,T7 arg7_input,T8 arg8_input,T9 arg9_input,T10 arg10_input,T11 arg11_input,T12 arg12_input,T13 arg13_input)
        :my_class(class_input),func(func_input),start_index(start_index_input),end_index(end_index_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input),arg4(arg4_input),arg5(arg5_input),arg6(arg6_input),arg7(arg7_input),arg8(arg8_input),arg9(arg9_input),arg10(arg10_input),arg11(arg11_input),arg12(arg12_input),arg13(arg13_input)
    {}
    
    void Run()
    {
        (my_class.*func)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,start_index,end_index);
    }
};

template<class TYPE,class T_ITERATOR>
class ITERATOR_TASK_0:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T_ITERATOR&);
    T_ITERATOR iterator;

    ITERATOR_TASK_0(TYPE& class_input,void (TYPE::*func_input)(T_ITERATOR&),const T_ITERATOR& iterator_input)
        :my_class(class_input),func(func_input),iterator(iterator_input)
    {}
    
    void Run()
    {
        (my_class.*func)(iterator);
    }
};

template<class TYPE,class T_ITERATOR,class T1>
class ITERATOR_TASK_1:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T_ITERATOR&,T1);
    T_ITERATOR iterator;
    T1 arg1;

    ITERATOR_TASK_1(TYPE& class_input,void (TYPE::*func_input)(T_ITERATOR&,T1),const T_ITERATOR& iterator_input,T1 arg1_input)
        :my_class(class_input),func(func_input),iterator(iterator_input),arg1(arg1_input)
    {}
    
    void Run()
    {
        (my_class.*func)(iterator,arg1);
    }
};

template<class TYPE,class T_ITERATOR,class T1,class T2>
class ITERATOR_TASK_2:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T_ITERATOR&,T1,T2);
    T_ITERATOR iterator;
    T1 arg1;T2 arg2;

    ITERATOR_TASK_2(TYPE& class_input,void (TYPE::*func_input)(T_ITERATOR&,T1,T2),const T_ITERATOR& iterator_input,T1 arg1_input,T2 arg2_input)
        :my_class(class_input),func(func_input),iterator(iterator_input),arg1(arg1_input),arg2(arg2_input)
    {}
    
    void Run()
    {
        (my_class.*func)(iterator,arg1,arg2);
    }
};

template<class TYPE,class T_ITERATOR,class T1,class T2,class T3>
class ITERATOR_TASK_3:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T_ITERATOR&,T1,T2,T3);
    T_ITERATOR iterator;
    T1 arg1;T2 arg2;T3 arg3;

    ITERATOR_TASK_3(TYPE& class_input,void (TYPE::*func_input)(T_ITERATOR&,T1,T2,T3),const T_ITERATOR& iterator_input,T1 arg1_input,T2 arg2_input,T3 arg3_input)
        :my_class(class_input),func(func_input),iterator(iterator_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input)
    {}
    
    void Run()
    {
        (my_class.*func)(iterator,arg1,arg2,arg3);
    }
};

template<class TYPE,class T_ITERATOR,class T1,class T2,class T3,class T4>
class ITERATOR_TASK_4:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T_ITERATOR&,T1,T2,T3,T4);
    T_ITERATOR iterator;
    T1 arg1;T2 arg2;T3 arg3;T4 arg4;

    ITERATOR_TASK_4(TYPE& class_input,void (TYPE::*func_input)(T_ITERATOR&,T1,T2,T3,T4),const T_ITERATOR& iterator_input,T1 arg1_input,T2 arg2_input,T3 arg3_input,T4 arg4_input)
        :my_class(class_input),func(func_input),iterator(iterator_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input),arg4(arg4_input)
    {}
    
    void Run()
    {
        (my_class.*func)(iterator,arg1,arg2,arg3,arg4);
    }
};

template<class TYPE,class T_ITERATOR,class T1,class T2,class T3,class T4,class T5>
class ITERATOR_TASK_5:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T_ITERATOR&,T1,T2,T3,T4,T5);
    T_ITERATOR iterator;
    T1 arg1;T2 arg2;T3 arg3;T4 arg4;T5 arg5;

    ITERATOR_TASK_5(TYPE& class_input,void (TYPE::*func_input)(T_ITERATOR&,T1,T2,T3,T4,T5),const T_ITERATOR& iterator_input,T1 arg1_input,T2 arg2_input,T3 arg3_input,T4 arg4_input,T5 arg5_input)
        :my_class(class_input),func(func_input),iterator(iterator_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input),arg4(arg4_input),arg5(arg5_input)
    {}
    
    void Run()
    {
        (my_class.*func)(iterator,arg1,arg2,arg3,arg4,arg5);
    }
};

template<class TYPE,class T_ITERATOR,class T1,class T2,class T3,class T4,class T5,class T6>
class ITERATOR_TASK_6:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T_ITERATOR&,T1,T2,T3,T4,T5,T6);
    T_ITERATOR iterator;
    T1 arg1;T2 arg2;T3 arg3;T4 arg4;T5 arg5;T6 arg6;

    ITERATOR_TASK_6(TYPE& class_input,void (TYPE::*func_input)(T_ITERATOR&,T1,T2,T3,T4,T5,T6),const T_ITERATOR& iterator_input,T1 arg1_input,T2 arg2_input,T3 arg3_input,T4 arg4_input,T5 arg5_input,T6 arg6_input)
        :my_class(class_input),func(func_input),iterator(iterator_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input),arg4(arg4_input),arg5(arg5_input),arg6(arg6_input)
    {}
    
    void Run()
    {
        (my_class.*func)(iterator,arg1,arg2,arg3,arg4,arg5,arg6);
    }
};

template<class TYPE,class T_ITERATOR,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
class ITERATOR_TASK_7:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T_ITERATOR&,T1,T2,T3,T4,T5,T6,T7);
    T_ITERATOR iterator;
    T1 arg1;T2 arg2;T3 arg3;T4 arg4;T5 arg5;T6 arg6;T7 arg7;

    ITERATOR_TASK_7(TYPE& class_input,void (TYPE::*func_input)(T_ITERATOR&,T1,T2,T3,T4,T5,T6,T7),const T_ITERATOR& iterator_input,T1 arg1_input,T2 arg2_input,T3 arg3_input,T4 arg4_input,T5 arg5_input,T6 arg6_input,T7 arg7_input)
        :my_class(class_input),func(func_input),iterator(iterator_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input),arg4(arg4_input),arg5(arg5_input),arg6(arg6_input),arg7(arg7_input)
    {}
    
    void Run()
    {
        (my_class.*func)(iterator,arg1,arg2,arg3,arg4,arg5,arg6,arg7);
    }
};

template<class TYPE,class T_ITERATOR,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
class ITERATOR_TASK_8:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T_ITERATOR&,T1,T2,T3,T4,T5,T6,T7,T8);
    T_ITERATOR iterator;
    T1 arg1;T2 arg2;T3 arg3;T4 arg4;T5 arg5;T6 arg6;T7 arg7;T8 arg8;

    ITERATOR_TASK_8(TYPE& class_input,void (TYPE::*func_input)(T_ITERATOR&,T1,T2,T3,T4,T5,T6,T7,T8),const T_ITERATOR& iterator_input,T1 arg1_input,T2 arg2_input,T3 arg3_input,T4 arg4_input,T5 arg5_input,T6 arg6_input,T7 arg7_input,T8 arg8_input)
        :my_class(class_input),func(func_input),iterator(iterator_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input),arg4(arg4_input),arg5(arg5_input),arg6(arg6_input),arg7(arg7_input),arg8(arg8_input)
    {}
    
    void Run()
    {
        (my_class.*func)(iterator,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8);
    }
};

template<class TYPE,class T_ITERATOR,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9>
class ITERATOR_TASK_9:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T_ITERATOR&,T1,T2,T3,T4,T5,T6,T7,T8,T9);
    T_ITERATOR iterator;
    T1 arg1;T2 arg2;T3 arg3;T4 arg4;T5 arg5;T6 arg6;T7 arg7;T8 arg8;T9 arg9;

    ITERATOR_TASK_9(TYPE& class_input,void (TYPE::*func_input)(T_ITERATOR&,T1,T2,T3,T4,T5,T6,T7,T8,T9),const T_ITERATOR& iterator_input,T1 arg1_input,T2 arg2_input,T3 arg3_input,T4 arg4_input,T5 arg5_input,T6 arg6_input,T7 arg7_input,T8 arg8_input,T9 arg9_input)
        :my_class(class_input),func(func_input),iterator(iterator_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input),arg4(arg4_input),arg5(arg5_input),arg6(arg6_input),arg7(arg7_input),arg8(arg8_input),arg9(arg9_input)
    {}
    
    void Run()
    {
        (my_class.*func)(iterator,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9);
    }
};

template<class TYPE,class T_ITERATOR,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9,class T10,class T11>
class ITERATOR_TASK_11:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T_ITERATOR&,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11);
    T_ITERATOR iterator;
    T1 arg1;T2 arg2;T3 arg3;T4 arg4;T5 arg5;T6 arg6;T7 arg7;T8 arg8;T9 arg9;T10 arg10;T11 arg11;

    ITERATOR_TASK_11(TYPE& class_input,void (TYPE::*func_input)(T_ITERATOR&,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11),const T_ITERATOR& iterator_input,T1 arg1_input,T2 arg2_input,T3 arg3_input,T4 arg4_input,T5 arg5_input,T6 arg6_input,T7 arg7_input,T8 arg8_input,T9 arg9_input,T10 arg10_input,T11 arg11_input)
        :my_class(class_input),func(func_input),iterator(iterator_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input),arg4(arg4_input),arg5(arg5_input),arg6(arg6_input),arg7(arg7_input),arg8(arg8_input),arg9(arg9_input),arg10(arg10_input),arg11(arg11_input)
    {}
    
    void Run()
    {
        (my_class.*func)(iterator,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11);
    }
};

template<class TYPE,class T_ITERATOR,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9,class T10,class T11,class T12>
class ITERATOR_TASK_12:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T_ITERATOR&,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12);
    T_ITERATOR iterator;
    T1 arg1;T2 arg2;T3 arg3;T4 arg4;T5 arg5;T6 arg6;T7 arg7;T8 arg8;T9 arg9;T10 arg10;T11 arg11;T12 arg12;

    ITERATOR_TASK_12(TYPE& class_input,void (TYPE::*func_input)(T_ITERATOR&,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12),const T_ITERATOR& iterator_input,T1 arg1_input,T2 arg2_input,T3 arg3_input,T4 arg4_input,T5 arg5_input,T6 arg6_input,T7 arg7_input,T8 arg8_input,T9 arg9_input,T10 arg10_input,T11 arg11_input,T12 arg12_input)
        :my_class(class_input),func(func_input),iterator(iterator_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input),arg4(arg4_input),arg5(arg5_input),arg6(arg6_input),arg7(arg7_input),arg8(arg8_input),arg9(arg9_input),arg10(arg10_input),arg11(arg11_input),arg12(arg12_input)
    {}
    
    void Run()
    {
        (my_class.*func)(iterator,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12);
    }
};

template<class TYPE,class T_ITERATOR,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9,class T10,class T11,class T12,class T13>
class ITERATOR_TASK_13:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T_ITERATOR&,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13);
    T_ITERATOR iterator;
    T1 arg1;T2 arg2;T3 arg3;T4 arg4;T5 arg5;T6 arg6;T7 arg7;T8 arg8;T9 arg9;T10 arg10;T11 arg11;T12 arg12;T13 arg13;

    ITERATOR_TASK_13(TYPE& class_input,void (TYPE::*func_input)(T_ITERATOR&,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13),const T_ITERATOR& iterator_input,T1 arg1_input,T2 arg2_input,T3 arg3_input,T4 arg4_input,T5 arg5_input,T6 arg6_input,T7 arg7_input,T8 arg8_input,T9 arg9_input,T10 arg10_input,T11 arg11_input,T12 arg12_input,T13 arg13_input)
        :my_class(class_input),func(func_input),iterator(iterator_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input),arg4(arg4_input),arg5(arg5_input),arg6(arg6_input),arg7(arg7_input),arg8(arg8_input),arg9(arg9_input),arg10(arg10_input),arg11(arg11_input),arg12(arg12_input),arg13(arg13_input)
    {}
    
    void Run()
    {
        (my_class.*func)(iterator,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13);
    }
};

template<class TYPE,class T_ITERATOR,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9,class T10,class T11,class T12,class T13,class T14>
class ITERATOR_TASK_14:public THREAD_QUEUE::TASK
{   
public:
    TYPE& my_class;
    void (TYPE::*func)(T_ITERATOR&,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14);
    T_ITERATOR iterator;
    T1 arg1;T2 arg2;T3 arg3;T4 arg4;T5 arg5;T6 arg6;T7 arg7;T8 arg8;T9 arg9;T10 arg10;T11 arg11;T12 arg12;T13 arg13;T14 arg14;

    ITERATOR_TASK_14(TYPE& class_input,void (TYPE::*func_input)(T_ITERATOR&,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14),const T_ITERATOR& iterator_input,T1 arg1_input,T2 arg2_input,T3 arg3_input,T4 arg4_input,T5 arg5_input,T6 arg6_input,T7 arg7_input,T8 arg8_input,T9 arg9_input,T10 arg10_input,T11 arg11_input,T12 arg12_input,T13 arg13_input,T14 arg14_input)
        :my_class(class_input),func(func_input),iterator(iterator_input),arg1(arg1_input),arg2(arg2_input),arg3(arg3_input),arg4(arg4_input),arg5(arg5_input),arg6(arg6_input),arg7(arg7_input),arg8(arg8_input),arg9(arg9_input),arg10(arg10_input),arg11(arg11_input),arg12(arg12_input),arg13(arg13_input),arg14(arg14_input)
    {}
    
    void Run()
    {
        (my_class.*func)(iterator,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14);
    }
};
}
#endif
#endif

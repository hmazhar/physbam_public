//#####################################################################
// Copyright 2010, Michael Lentine, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DOMAIN_ITERATOR_THREADED
//#####################################################################
#ifndef __DOMAIN_ITERATOR_THREADED__
#define __DOMAIN_ITERATOR_THREADED__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/ITERATOR_TASK.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <climits>
namespace PhysBAM{

template<class TYPE,class TV>
class DOMAIN_ITERATOR_THREADED_ALPHA
{
    typedef VECTOR<int,TV::dimension> TV_INT;
    typedef typename TV::SCALAR T;

public:
    THREAD_QUEUE* thread_queue;
    ARRAY<RANGE<TV_INT> > domains;
    int domain_type;
    int number_of_domains;
    bool do_wait;
    TV_INT split_per_dimension;

    DOMAIN_ITERATOR_THREADED_ALPHA(const RANGE<TV_INT>& domain,THREAD_QUEUE* thread_queue_input,const int axis=1,const int overlap_rows=0,const int domain_type_input=1,const int number_of_domains_multiplier=4,const int row_jump_override=0,const bool do_wait_input=true)
        :thread_queue(thread_queue_input),domain_type(domain_type_input),do_wait(do_wait_input)
    {
        int total_rows=domain.max_corner(axis)-domain.min_corner(axis)+1;
        if(thread_queue && total_rows>0){
            number_of_domains=number_of_domains_multiplier*thread_queue->Number_Of_Threads();
            if(domain_type==1){
                if(row_jump_override<=0){
                    number_of_domains=min(number_of_domains,total_rows); // have to have at least 1 row per domain
                    int row_jump=total_rows/number_of_domains;
                    int row_remainder=total_rows%number_of_domains;
                    domains.Resize(number_of_domains);
                    int previous_domain_max=domain.min_corner(axis)-1;
                    for(int i=1;i<=number_of_domains;i++){
                        domains(i)=domain;
                        domains(i).min_corner(axis)=previous_domain_max+1;
                        previous_domain_max+=row_jump;
                        if(i<=row_remainder)previous_domain_max+=1; // if this is one of the early domains and we have more remainders, add it to this one
                        domains(i).max_corner(axis)=previous_domain_max+overlap_rows;
                        domains(i).max_corner(axis)=min(domains(i).max_corner(axis),domain.max_corner(axis));}
                    assert(previous_domain_max==domain.max_corner(axis));} // if this isn't true, something is wrong with the above logic
                else PHYSBAM_FATAL_ERROR();} // not implemented yet
            else if(domain_type==2){
                ARRAY<ARRAY<int> > boundaries;
                split_per_dimension=Split_Range_Wrapper(split_per_dimension,domain,boundaries);
                //TODO: sync split_per_dimension and update if diff
                domains.Resize(split_per_dimension.Product());
                GRID<TV> process_grid(split_per_dimension,RANGE<TV>::Centered_Box());
                int count=0;
                for(typename GRID<TV>::NODE_ITERATOR iterator(process_grid);iterator.Valid();iterator.Next()){
                    TV_INT coordinates=iterator.Node_Index();
                    TV_INT start,end;
                    for(int axis=1;axis<=TV::dimension;axis++){
                        start[axis]=boundaries(axis)(coordinates[axis])+1;
                        end[axis]=boundaries(axis)(coordinates[axis]+1);}
                    domains(++count).min_corner=start+domain.min_corner-TV_INT::Constant_Vector(overlap_rows);
                    domains(count).max_corner=end+domain.min_corner+TV_INT::Constant_Vector(overlap_rows);
                    for(int axis=1;axis<=TV::dimension;axis++){
                        domains(count).min_corner(axis)=max(domains(count).min_corner(axis),domain.min_corner(axis));
                        domains(count).max_corner(axis)=min(domains(count).max_corner(axis),domain.max_corner(axis));}}}
            else PHYSBAM_FATAL_ERROR();}
        else{domains.Resize(1);domains(1)=domain;number_of_domains=1;}
    }

    TV_INT Split_Range_Wrapper(TV_INT& processes_per_dimension,const RANGE<TV_INT>& global_range,ARRAY<ARRAY<int> >& boundaries)
    {
        TV_INT range=global_range.max_corner-global_range.min_corner;
        return Split_Range(range,processes_per_dimension,boundaries);
    }

    static bool Minimize_2D_Surface_Area(const int number_of_domains,const int x,const int y,int& count_x)
    {
        count_x=max(1,(int)sqrt((T)number_of_domains*x/y));
        if(number_of_domains%count_x==0){
            if(number_of_domains%(count_x+1)==0 && y*count_x+x*number_of_domains/count_x >= y*(count_x+1)+x*number_of_domains/(count_x+1)) count_x++;}
        else if(number_of_domains%(count_x+1)==0) count_x++;
        else return false;
        return true;
    }

    VECTOR<int,1> Split_Range(const VECTOR<int,1>& global_range,const VECTOR<int,1>& processes_per_dimension,ARRAY<ARRAY<int> >& boundaries)
    {
        PHYSBAM_ASSERT(!processes_per_dimension.x || processes_per_dimension.x==number_of_domains);
        boundaries.Resize(1);
        Split_Dimension(global_range.x,number_of_domains,boundaries(1));
        return VECTOR<int,1>(number_of_domains);
    }

    VECTOR<int,2> Split_Range(const VECTOR<int,2>& global_range,const VECTOR<int,2>& processes_per_dimension,ARRAY<ARRAY<int> >& boundaries)
    {
        int x=global_range.x,y=global_range.y;VECTOR<int,2> count;
        if(processes_per_dimension!=VECTOR<int,2>()) count=processes_per_dimension;
        else{ // try to figure out counts by minimizing surface area between processes
            if(!Minimize_2D_Surface_Area(number_of_domains,x,y,count.x)){
                LOG::cerr<<"Don't know how to divide domain in both directions."<<std::endl;PHYSBAM_NOT_IMPLEMENTED();}
            count.y=number_of_domains/count.x;}
        PHYSBAM_ASSERT(count.x*count.y==number_of_domains);
        boundaries.Resize(2);
        Split_Dimension(x,count.x,boundaries(1));
        Split_Dimension(y,count.y,boundaries(2));
        return count;
    }

    VECTOR<int,3> Split_Range(const VECTOR<int,3>& global_range,const VECTOR<int,3>& processes_per_dimension,ARRAY<ARRAY<int> >& boundaries)
    {
        int x=global_range.x,y=global_range.y,z=global_range.z;VECTOR<int,3> count;
        if(processes_per_dimension!=VECTOR<int,3>()) count=processes_per_dimension;
        else{ // try to figure out counts by minimizing surface area between processes
            T minimum_surface_area=FLT_MAX;VECTOR<int,3> test_count;
            for(test_count.z=1;test_count.z<=number_of_domains;test_count.z++) if(number_of_domains%test_count.z==0){
                if(Minimize_2D_Surface_Area(number_of_domains/test_count.z,x,y,test_count.x)){
                    test_count.y=number_of_domains/(test_count.x*test_count.z);
                    T surface_area=test_count.x*(y*z)+test_count.y*(x*z)+test_count.z*(x*y);
                    if(surface_area<minimum_surface_area){count=test_count;minimum_surface_area=surface_area;}}
                if(Minimize_2D_Surface_Area(number_of_domains/test_count.z,x,y,test_count.y)){
                    test_count.x=number_of_domains/(test_count.y*test_count.z);
                    T surface_area=test_count.x*(y*z)+test_count.y*(x*z)+test_count.z*(x*y);
                    if(surface_area<minimum_surface_area){count=test_count;minimum_surface_area=surface_area;}}}
            if(minimum_surface_area==INT_MAX){LOG::cerr<<"Don't know how to divide domain in all directions."<<std::endl;PHYSBAM_NOT_IMPLEMENTED();}}
        PHYSBAM_ASSERT(count.x*count.y*count.z==number_of_domains);
        boundaries.Resize(3);
        Split_Dimension(x,count.x,boundaries(1));
        Split_Dimension(y,count.y,boundaries(2));
        Split_Dimension(z,count.z,boundaries(3));
        return count;
    }

    void Split_Dimension(const int x,const int processes,ARRAY<int>& boundaries)
    {
        int range_over_processes=(x+1)/processes;
        int remainder=(x+1)%processes;
        boundaries.Resize(processes+1);boundaries(1)=-1;
        for(int p=1;p<=processes;p++){
            boundaries(p+1)=boundaries(p)+range_over_processes;
            if(p<=remainder)boundaries(p+1)+=1;}
    }

    void Run(TYPE& my_class,void (TYPE::*func)(RANGE<TV_INT>&))
    {
        if(!thread_queue){(my_class.*func)(domains(1));return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=domains.m;i++){
            ITERATOR_TASK_0<TYPE,RANGE<TV_INT> >* task=new ITERATOR_TASK_0<TYPE,RANGE<TV_INT> >(my_class,func,domains(i));
            thread_queue->Queue(task);}
        if(do_wait)thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1> void Run(TYPE& my_class,void (TYPE::*func)(RANGE<TV_INT>&,T1),T1 arg1)
    {
        if(!thread_queue){(my_class.*func)(domains(1),arg1);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=domains.m;i++){
            ITERATOR_TASK_1<TYPE,RANGE<TV_INT>,T1>* task=new ITERATOR_TASK_1<TYPE,RANGE<TV_INT>,T1>(my_class,func,domains(i),arg1);
            thread_queue->Queue(task);}
        if(do_wait)thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1,class T2> void Run(TYPE& my_class,void (TYPE::*func)(RANGE<TV_INT>&,T1,T2),T1 arg1,T2 arg2)
    {
        if(!thread_queue){(my_class.*func)(domains(1),arg1,arg2);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=domains.m;i++){
            ITERATOR_TASK_2<TYPE,RANGE<TV_INT>,T1,T2>* task=new ITERATOR_TASK_2<TYPE,RANGE<TV_INT>,T1,T2>(my_class,func,domains(i),arg1,arg2);
            thread_queue->Queue(task);}
        if(do_wait)thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1,class T2,class T3> void Run(TYPE& my_class,void (TYPE::*func)(RANGE<TV_INT>&,T1,T2,T3),T1 arg1,T2 arg2,T3 arg3)
    {
        if(!thread_queue){(my_class.*func)(domains(1),arg1,arg2,arg3);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=domains.m;i++){
            ITERATOR_TASK_3<TYPE,RANGE<TV_INT>,T1,T2,T3>* task=new ITERATOR_TASK_3<TYPE,RANGE<TV_INT>,T1,T2,T3>(my_class,func,domains(i),arg1,arg2,arg3);
            thread_queue->Queue(task);}
        if(do_wait)thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }
 
    template<class T1,class T2,class T3,class T4> void Run(TYPE& my_class,void (TYPE::*func)(RANGE<TV_INT>&,T1,T2,T3,T4),T1 arg1,T2 arg2,T3 arg3,T4 arg4)
    {
        if(!thread_queue){(my_class.*func)(domains(1),arg1,arg2,arg3,arg4);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=domains.m;i++){
            ITERATOR_TASK_4<TYPE,RANGE<TV_INT>,T1,T2,T3,T4>* task=new ITERATOR_TASK_4<TYPE,RANGE<TV_INT>,T1,T2,T3,T4>(my_class,func,domains(i),arg1,arg2,arg3,arg4);
            thread_queue->Queue(task);}
        if(do_wait)thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1,class T2,class T3,class T4,class T5> void Run(TYPE& my_class,void (TYPE::*func)(RANGE<TV_INT>&,T1,T2,T3,T4,T5),T1 arg1,T2 arg2,T3 arg3,T4 arg4,T5 arg5)
    {
        if(!thread_queue){(my_class.*func)(domains(1),arg1,arg2,arg3,arg4,arg5);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=domains.m;i++){
            ITERATOR_TASK_5<TYPE,RANGE<TV_INT>,T1,T2,T3,T4,T5>* task=new ITERATOR_TASK_5<TYPE,RANGE<TV_INT>,T1,T2,T3,T4,T5>(my_class,func,domains(i),arg1,arg2,arg3,arg4,arg5);
            thread_queue->Queue(task);}
        if(do_wait)thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1,class T2,class T3,class T4,class T5,class T6> void Run(TYPE& my_class,void (TYPE::*func)(RANGE<TV_INT>&,T1,T2,T3,T4,T5,T6),T1 arg1,T2 arg2,T3 arg3,T4 arg4,T5 arg5,T6 arg6)
    {
        if(!thread_queue){(my_class.*func)(domains(1),arg1,arg2,arg3,arg4,arg5,arg6);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=domains.m;i++){
            ITERATOR_TASK_6<TYPE,RANGE<TV_INT>,T1,T2,T3,T4,T5,T6>* task=new ITERATOR_TASK_6<TYPE,RANGE<TV_INT>,T1,T2,T3,T4,T5,T6>(my_class,func,domains(i),arg1,arg2,arg3,arg4,arg5,arg6);
            thread_queue->Queue(task);}
        if(do_wait)thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1,class T2,class T3,class T4,class T5,class T6,class T7> void Run(TYPE& my_class,void (TYPE::*func)(RANGE<TV_INT>&,T1,T2,T3,T4,T5,T6,T7),T1 arg1,T2 arg2,T3 arg3,T4 arg4,T5 arg5,T6 arg6,T7 arg7)
    {
        if(!thread_queue){(my_class.*func)(domains(1),arg1,arg2,arg3,arg4,arg5,arg6,arg7);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=domains.m;i++){
            ITERATOR_TASK_7<TYPE,RANGE<TV_INT>,T1,T2,T3,T4,T5,T6,T7>* task=new ITERATOR_TASK_7<TYPE,RANGE<TV_INT>,T1,T2,T3,T4,T5,T6,T7>(my_class,func,domains(i),arg1,arg2,arg3,arg4,arg5,arg6,arg7);
            thread_queue->Queue(task);}
        if(do_wait)thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8> void Run(TYPE& my_class,void (TYPE::*func)(RANGE<TV_INT>&,T1,T2,T3,T4,T5,T6,T7,T8),T1 arg1,T2 arg2,T3 arg3,T4 arg4,T5 arg5,T6 arg6,T7 arg7,T8 arg8)
    {
        if(!thread_queue){(my_class.*func)(domains(1),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=domains.m;i++){
            ITERATOR_TASK_8<TYPE,RANGE<TV_INT>,T1,T2,T3,T4,T5,T6,T7,T8>* task=new ITERATOR_TASK_8<TYPE,RANGE<TV_INT>,T1,T2,T3,T4,T5,T6,T7,T8>(my_class,func,domains(i),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8);
            thread_queue->Queue(task);}
        thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9> void Run(TYPE& my_class,void (TYPE::*func)(RANGE<TV_INT>&,T1,T2,T3,T4,T5,T6,T7,T8,T9),T1 arg1,T2 arg2,T3 arg3,T4 arg4,T5 arg5,T6 arg6,T7 arg7,T8 arg8,T9 arg9)
    {
        if(!thread_queue){(my_class.*func)(domains(1),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=domains.m;i++){
            ITERATOR_TASK_9<TYPE,RANGE<TV_INT>,T1,T2,T3,T4,T5,T6,T7,T8,T9>* task=new ITERATOR_TASK_9<TYPE,RANGE<TV_INT>,T1,T2,T3,T4,T5,T6,T7,T8,T9>(my_class,func,domains(i),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9);
            thread_queue->Queue(task);}
        if(do_wait)thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9,class T10,class T11> void Run(TYPE& my_class,void (TYPE::*func)(RANGE<TV_INT>&,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11),T1 arg1,T2 arg2,T3 arg3,T4 arg4,T5 arg5,T6 arg6,T7 arg7,T8 arg8,T9 arg9,T10 arg10,T11 arg11)
    {
        if(!thread_queue){(my_class.*func)(domains(1),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=domains.m;i++){
            ITERATOR_TASK_11<TYPE,RANGE<TV_INT>,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>* task=new ITERATOR_TASK_11<TYPE,RANGE<TV_INT>,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11>(my_class,func,domains(i),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11);
            thread_queue->Queue(task);}
        if(do_wait)thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9,class T10,class T11,class T12> void Run(TYPE& my_class,void (TYPE::*func)(RANGE<TV_INT>&,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12),T1 arg1,T2 arg2,T3 arg3,T4 arg4,T5 arg5,T6 arg6,T7 arg7,T8 arg8,T9 arg9,T10 arg10,T11 arg11,T12 arg12)
    {
        if(!thread_queue){(my_class.*func)(domains(1),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=domains.m;i++){
            ITERATOR_TASK_12<TYPE,RANGE<TV_INT>,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12>* task=new ITERATOR_TASK_12<TYPE,RANGE<TV_INT>,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12>(my_class,func,domains(i),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12);
            thread_queue->Queue(task);}
        if(do_wait)thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9,class T10,class T11,class T12,class T13> void Run(TYPE& my_class,void (TYPE::*func)(RANGE<TV_INT>&,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13),T1 arg1,T2 arg2,T3 arg3,T4 arg4,T5 arg5,T6 arg6,T7 arg7,T8 arg8,T9 arg9,T10 arg10,T11 arg11,T12 arg12,T13 arg13)
    {
        if(!thread_queue){(my_class.*func)(domains(1),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=domains.m;i++){
            ITERATOR_TASK_13<TYPE,RANGE<TV_INT>,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13>* task=new ITERATOR_TASK_13<TYPE,RANGE<TV_INT>,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13>(my_class,func,domains(i),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13);
            thread_queue->Queue(task);}
        if(do_wait)thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }

    template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9,class T10,class T11,class T12,class T13,class T14> void Run(TYPE& my_class,void (TYPE::*func)(RANGE<TV_INT>&,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14),T1 arg1,T2 arg2,T3 arg3,T4 arg4,T5 arg5,T6 arg6,T7 arg7,T8 arg8,T9 arg9,T10 arg10,T11 arg11,T12 arg12,T13 arg13,T14 arg14)
    {
        if(!thread_queue){(my_class.*func)(domains(1),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14);return;}
#ifdef USE_PTHREADS
        for(int i=1;i<=domains.m;i++){
            ITERATOR_TASK_14<TYPE,RANGE<TV_INT>,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14>* task=new ITERATOR_TASK_14<TYPE,RANGE<TV_INT>,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14>(my_class,func,domains(i),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14);
            thread_queue->Queue(task);}
        if(do_wait)thread_queue->Wait();
#else
        PHYSBAM_FATAL_ERROR();
#endif
    }
};
}
#endif

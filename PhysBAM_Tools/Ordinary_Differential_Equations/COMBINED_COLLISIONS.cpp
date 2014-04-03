//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/COMBINED_COLLISIONS.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/COMBINED_COLLISIONS_NONLINEAR_FUNCTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS<TV>::
COMBINED_COLLISIONS(IMPULSE* accumulator_input)
    :accumulator(accumulator_input),test_system(false)
{
    accumulator->Clear();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS<TV>::
~COMBINED_COLLISIONS()
{
}
//#####################################################################
// Function Add_Collider
//#####################################################################
template<class TV> void COMBINED_COLLISIONS<TV>::
Add_Collider(COLLIDER* collider)
{
    colliders.Append(collider);
}
//#####################################################################
// Function Apply_Combined_Impulse
//#####################################################################
template<class TV> void COMBINED_COLLISIONS<TV>::
Apply_Combined_Impulse(T dt,T time)
{
    accumulator->Resize();
    collider_info.Resize(colliders.m);
    accumulator->Setup_Discover_State(dt,time);
    for(int i=1;i<=colliders.m;i++){
        colliders(i)->Discover(dt,time);
        int n=colliders(i)->Count();
        collider_info(i).scale.Resize(n);
        for(int j=1;j<=n;j++)
            colliders(i)->Precompute(j);

        ARRAYS_COMPUTATIONS::Fill(collider_info(i).scale,(T)1);}

    Refresh_Referenced_List();

    accumulator->Setup_Impulse_State(referenced_list,dt,time);

    if(test_system) Test_System();
    Compute_Impulse_Using_Scale(dt,time);
    for(int i=1;i<=5;i++){
        Update_Scales_With_Diagnosis(dt,time);
        Compute_Impulse_Using_Scale(dt,time);}
    Apply_Computed_Impulses(dt,time);
    Notify_Moved();
    accumulator->Finish_Step(referenced_list,dt,time);
}
//#####################################################################
// Function Apply_Combined_Impulse
//#####################################################################
template<class TV> void COMBINED_COLLISIONS<TV>::
Apply_Combined_Projections(T dt,T time)
{
    // affect velocities only, not positions
    // alternately, affect both, but save and restore
    accumulator->Resize();
    collider_info.Resize(colliders.m);
    accumulator->Setup_Discover_State(dt,time);
    for(int i=1;i<=colliders.m;i++){
        colliders(i)->Discover(dt,time);
        int n=colliders(i)->Count();
        collider_info(i).scale.Resize(n);
        for(int j=1;j<=n;j++)
            colliders(i)->Precompute(j);

        ARRAYS_COMPUTATIONS::Fill(collider_info(i).scale,(T)1);}

    Refresh_Referenced_List();

    accumulator->Setup_Impulse_State(referenced_list,dt,time);

    accumulator->Clear(referenced_list);
    for(int iterations=1;iterations<=10;iterations++){
        for(int i=1;i<=colliders.m;i++)
            for(int j=1;j<=collider_info(i).scale.m;j++)
                colliders(i)->Accumulate_Impulse(j,*accumulator,collider_info(i).scale(j),dt,time);}
    Recompute_Active_Set(dt,time);
    Apply_Computed_Impulses(dt,time);
    Notify_Moved();
}
//#####################################################################
// Function Refresh_Referenced_List
//#####################################################################
template<class TV> void COMBINED_COLLISIONS<TV>::
Refresh_Referenced_List()
{
    referenced_particles.Remove_All();
    referenced_list.Remove_All();

    for(int i=1;i<=colliders.m;i++)
        for(int j=1;j<=collider_info(i).scale.m;j++)
            colliders(i)->Affected_Bodies(j,referenced_particles);

    referenced_particles.Get_Keys(referenced_list);
}
//#####################################################################
// Function Update_Scales_With_Diagnosis
//#####################################################################
template<class TV> void COMBINED_COLLISIONS<TV>::
Update_Scales_With_Diagnosis(T dt,T time)
{
    for(int i=1;i<=colliders.m;i++){
        for(int j=1;j<=collider_info(i).scale.m;j++){
            T progress=colliders(i)->Diagnose_Impulse(j,*accumulator,dt,time);
            if(progress<(T).01) progress=0;
            if(progress)
                collider_info(i).scale(j)/=progress;
            else
                collider_info(i).scale(j)=0;}}
}
//#####################################################################
// Function Compute_Impulse_Using_Scale
//#####################################################################
template<class TV> void COMBINED_COLLISIONS<TV>::
Compute_Impulse_Using_Scale(T dt,T time)
{
    accumulator->Clear(referenced_list);
    for(int i=1;i<=colliders.m;i++)
        for(int j=1;j<=collider_info(i).scale.m;j++)
            colliders(i)->Accumulate_Impulse(j,*accumulator,collider_info(i).scale(j),dt,time);
}
//#####################################################################
// Function Apply_Computed_Impulses
//#####################################################################
template<class TV> void COMBINED_COLLISIONS<TV>::
Apply_Computed_Impulses(T dt,T time)
{
    accumulator->Apply(referenced_list,dt,time);
}
//#####################################################################
// Function Solve
//#####################################################################
template<class TV> void COMBINED_COLLISIONS<TV>::
Solve(T dt,T time,int iterations)
{
    for(int i=1;i<=iterations;i++)
        Apply_Combined_Impulse(dt,time);
}
//#####################################################################
// Function Solve
//#####################################################################
template<class TV> void COMBINED_COLLISIONS<TV>::
Iterate_Projections(T dt,T time,int iterations)
{
    // run as projections, see which ones actually got used, run again with only those projections, rinse and repeat
    for(int i=1;i<=iterations;i++){
        Apply_Combined_Projections(dt,time);
        PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Pull in step %d",i),2,2);
    }
}
//#####################################################################
// Function Recompute_Active_Set
//#####################################################################
template<class TV> void COMBINED_COLLISIONS<TV>::
Recompute_Active_Set(T dt,T time)
{
    for(int i=1;i<=colliders.m;i++){
        int n=colliders(i)->Count();
        for(int j=1;j<=n;j++)
            colliders(i)->Diagnose_Impulse(j,*accumulator,dt,time);}
}
//#####################################################################
// Function Notify_Moved
//#####################################################################
template<class TV> void COMBINED_COLLISIONS<TV>::
Notify_Moved()
{
    for(int i=1;i<=colliders.m;i++)
        for(int j=1;j<=collider_info(i).scale.m;j++)
            colliders(i)->Moved(j);
}
//#####################################################################
// Function Test_System
//#####################################################################
template<class TV> void COMBINED_COLLISIONS<TV>::
Test_System()
{
    COMBINED_COLLISIONS_NONLINEAR_FUNCTION<TV> f(*this,1,0);
    f.Test_System();
}
//#####################################################################
// Function Find_Scales_With_CG
//#####################################################################
template<class TV> void COMBINED_COLLISIONS<TV>::
Find_Scales_With_CG(T dt,T time)
{
    ITERATIVE_SOLVER<T> is;
    is.tolerance=(T)1e-5;

    COMBINED_COLLISIONS_NONLINEAR_FUNCTION<TV> f(*this,dt,time);
    COMBINED_COLLISIONS_PARAMETER_SPACE<TV> x(*this);
    int n=0;
    for(int i=1;i<=x.parameters.m;i++) n+=x.parameters(i).m;
    is.Conjugate_Gradient(f,x,2,min(20,n));
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS<TV>::IMPULSE::
IMPULSE()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS<TV>::IMPULSE::
~IMPULSE()
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS<TV>::COLLIDER::
COLLIDER()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS<TV>::COLLIDER::
~COLLIDER()
{
}
template class COMBINED_COLLISIONS<VECTOR<float,1> >;
template class COMBINED_COLLISIONS<VECTOR<float,2> >;
template class COMBINED_COLLISIONS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMBINED_COLLISIONS<VECTOR<double,1> >;
template class COMBINED_COLLISIONS<VECTOR<double,2> >;
template class COMBINED_COLLISIONS<VECTOR<double,3> >;
#endif

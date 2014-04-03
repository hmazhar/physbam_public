//#####################################################################
// Copyright 2002, Ron Fedkiw, Eran Guendelman, Robert Bridson
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TURBULENCE  
//##################################################################### 
//
// Generates a periodic turbulent field from a random complex sequence.
//
//#####################################################################
#ifndef __TURBULENCE__
#define __TURBULENCE__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UNIFORM_FORWARD.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class T>
class TURBULENCE
{
private:
    RANDOM_NUMBERS<T> *random,random_default;
    int incompressible; // (1) incompressible, (0) compressible
    T k_inertial; // sets the scale of the turbulence
    T epsilon;    // dissipation
    T constant;  // Kolmogorov constant
    T rescaled_average_velocity; // rescales the computed velocity to have this as the average velocity
public:
    T time_start,time_end;

public:
    TURBULENCE() 
    {
        random=&random_default;
        Enforce_Incompressibility();
        Set_Lowest_Angular_Frequency();
        Set_Dissipation();
        Set_Constant();
        Set_Rescaled_Average_Velocity();
    }
    
    void Set_Custom_Random(RANDOM_NUMBERS<T>& random_input)
    {random=&random_input;}

    void Enforce_Incompressibility()
    {incompressible=1;}

    void Do_Not_Enforce_Incompressibility()
    {incompressible=0;}

    void Set_Lowest_Angular_Frequency(const T k_inertial_input=4) // k_inertial=2*pi*|frequency|
    {k_inertial=k_inertial_input;}

    void Set_Dissipation(const T epsilon_input=1)
    {epsilon=epsilon_input;}
    
    void Set_Constant(const T constant_input=1.5)
    {constant=constant_input;}

    void Set_Rescaled_Average_Velocity(const T velocity_magnitude=1)
    {rescaled_average_velocity=velocity_magnitude;}

//#####################################################################
    void Generate_Random_Turbulence(const GRID<VECTOR<T,2> >& grid,ARRAY<T,VECTOR<int,2> >& u,ARRAY<T,VECTOR<int,2> >& v) const;
    void Generate_Random_Turbulence(const GRID<VECTOR<T,3> >& grid,ARRAY<T,VECTOR<int,3> >& u,ARRAY<T,VECTOR<int,3> >& v,ARRAY<T,VECTOR<int,3> >& w) const;
//#####################################################################
};
}    
#endif

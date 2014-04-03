//#####################################################################
// Copyright 2009-2010, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __WATER_DRIVER__
#define __WATER_DRIVER__
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
namespace PhysBAM{


template<class TV> class WATER_EXAMPLE;

template<class TV>
class WATER_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef EXTRAPOLATION_UNIFORM<GRID<TV>,T> T_EXTRAPOLATION_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_BASE T_ARRAYS_BASE;
    
protected:
    int current_frame;
    T time;
    int output_number;

    WATER_EXAMPLE<TV>& example;
public:
    WATER_DRIVER(WATER_EXAMPLE<TV>& example);
    virtual ~WATER_DRIVER();

    void Scalar_Advance(const T dt,const T time);
    void Convect(const T dt,const T time);
    void Add_Forces(const T dt,const T time);
    void Project(const T dt,const T time);
    void Execute_Main_Program();
    void Initialize();
    void Advance_To_Target_Time(const T target_time);
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files(const int frame);
    void Write_Substep(const std::string& title,const int substep,const int level=0);

//#####################################################################
};
}
#endif

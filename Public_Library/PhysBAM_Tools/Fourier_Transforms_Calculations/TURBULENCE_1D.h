//#####################################################################
// Copyright 2007. Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TURBULENCE_1D  
//##################################################################### 
#ifndef __TURBULENCE_1D__
#define __TURBULENCE_1D__

#include <PhysBAM_Tools/Fourier_Transforms_Calculations/TURBULENCE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
namespace PhysBAM{

template<class T>
class TURBULENCE_1D:public TURBULENCE<T>
{
    typedef VECTOR<T,1> TV;
public:
    using TURBULENCE<T>::time_start;using TURBULENCE<T>::time_end;

public:
    TURBULENCE_1D()
    {}
    
    TURBULENCE_1D(const GRID<TV>& grid_input)
    {}

    void Initialize_Grid(const GRID<TV>& grid_input)
    {PHYSBAM_NOT_IMPLEMENTED();}

    void Generate_Initial_Turbulence(const T time_start_input=0,const T time_end_input=1)
    {PHYSBAM_NOT_IMPLEMENTED();}
    
    void Advance_Turbulence()
    {PHYSBAM_NOT_IMPLEMENTED();}

    VECTOR<T,1> Turbulent_Velocity(const VECTOR<T,1>& X,const T fraction) const
    {PHYSBAM_NOT_IMPLEMENTED();}

    T Turbulent_Face_Velocity(const int axis,const VECTOR<T,1>& X,const T fraction) const
    {PHYSBAM_NOT_IMPLEMENTED();}
    
    T Turbulent_U_Velocity(const VECTOR<T,1>& X,const T fraction) const
    {PHYSBAM_NOT_IMPLEMENTED();}

//#####################################################################
};
}    
#endif


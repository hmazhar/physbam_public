//#####################################################################
// Copyright 2002, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BLACKBODY
//##################################################################### 
#ifndef __BLACKBODY__
#define __BLACKBODY__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/CIE_XYZ.h>
namespace PhysBAM{

template<class T>
class BLACKBODY
{
public:
    CIE_XYZ<T> cie; 

    BLACKBODY() 
    {}

//#####################################################################
    void Calculate_Radiance_Spectrum(const T temperature,const GRID<VECTOR<T,1> >& grid,ARRAY<T,VECTOR<int,1> >& radiance_spectrum) const;
    VECTOR<T,3> Calculate_XYZ(const T temperature) const; 
//#####################################################################
};   
}
#endif

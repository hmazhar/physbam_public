//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_FACE_VELOCITY_DONOR_CELL_RLE
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __ADVECTION_FACE_VELOCITY_DONOR_CELL_RLE__
#define __ADVECTION_FACE_VELOCITY_DONOR_CELL_RLE__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
namespace PhysBAM{

template<class T> class RLE_GRID_2D;
template<class T> class RLE_GRID_3D;

template<class T>
class ADVECTION_FACE_VELOCITY_DONOR_CELL_RLE
{
public:
    bool clamp_divergence_fix;

    ADVECTION_FACE_VELOCITY_DONOR_CELL_RLE()
        :clamp_divergence_fix(false)
    {}

//#####################################################################
    void Euler_Step(const RLE_GRID_2D<T>& grid,ARRAY<T>& Z,const ARRAY<T>& Z_ghost,const ARRAY<T>& V,const T dt);
    void Euler_Step(const RLE_GRID_3D<T>& grid,ARRAY<T>& Z,const ARRAY<T>& Z_ghost,const ARRAY<T>& V,const T dt);
//#####################################################################
};
}
#endif
#endif

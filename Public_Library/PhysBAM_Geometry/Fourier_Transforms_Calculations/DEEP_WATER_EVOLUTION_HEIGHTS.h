//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __DEEP_WATER_EVOLUTION_HEIGHTS__
#define __DEEP_WATER_EVOLUTION_HEIGHTS__

#include <PhysBAM_Tools/Fourier_Transforms_Calculations/DEEP_WATER_EVOLUTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>

namespace PhysBAM{

template<class TV>
class DEEP_WATER_EVOLUTION_HEIGHTS:public DEEP_WATER_EVOLUTION<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<T,TV::dimension+1> TV_FULL;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_T;
    typedef typename GRID<TV>::NODE_ITERATOR NODE_ITERATOR;
    typedef DEEP_WATER_EVOLUTION<TV> BASE;
public:
    using BASE::grid;using BASE::h;using BASE::Xh;using BASE::lambda;using BASE::Set_H_Hats_From_Height;

    void Intersect_With_Geometry(const RIGID_GEOMETRY<TV_FULL>& rigid_geometry,const T collision_thickness);
//#####################################################################
};

}
#endif

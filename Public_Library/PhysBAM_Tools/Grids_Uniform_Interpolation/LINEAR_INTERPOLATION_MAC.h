//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __LINEAR_INTERPOLATION_MAC__
#define __LINEAR_INTERPOLATION_MAC__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
namespace PhysBAM{

template<class TV,class T2>
class LINEAR_INTERPOLATION_MAC
{
public:
    VECTOR<GRID<TV>,TV::m> grids;
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T2> interp;

    LINEAR_INTERPOLATION_MAC(GRID<TV>& grid)
    {
        for(int d=1;d<=TV::m;d++) grids(d)=grid.Get_Axis_X_Face_Grid(d);
    }

    VECTOR<T2,TV::m> Clamped_To_Array(const ARRAY<T2,FACE_INDEX<TV::m> >& V,const TV& X) const
    {VECTOR<T2,TV::m> Y;for(int d=1;d<=TV::m;d++) Y(d)=interp.Clamped_To_Array(grids(d),V.Component(d),X);return Y;}

//#####################################################################
};
}
#endif

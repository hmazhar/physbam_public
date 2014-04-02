//#####################################################################
// Copyright 2003-2005, Ronald Fedkiw, Eran Guendelman, Duc Nguyen.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_DYADIC
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_DYADIC__
#define __LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_DYADIC__

#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_DYADIC<T_GRID>
class LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_DYADIC:public LINEAR_INTERPOLATION_DYADIC<T_GRID,T2,T_FACE_LOOKUP>
{
    typedef typename T_GRID::SCALAR T;
public:
    template<class T3> struct REBIND{typedef LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_DYADIC<T_GRID,T3,T_FACE_LOOKUP> TYPE;};

    LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_DYADIC()
    {}
       
    T2 Clamped_To_Array(const OCTREE_GRID<T>& grid,const ARRAY<T2>& u,const VECTOR<T,3>& X) const
    {VECTOR<T,3> X_new(clamp(X.x,grid.uniform_grid.xmin,grid.uniform_grid.xmax),clamp(X.y,grid.uniform_grid.ymin,grid.uniform_grid.ymax),clamp(X.z,grid.uniform_grid.zmin,grid.uniform_grid.zmax));
    return LINEAR_INTERPOLATION_DYADIC<T,T2,OCTREE_GRID<T> >::Clamped_To_Array(grid,u,X_new);}

//#####################################################################
};
}
#endif
#endif

//#####################################################################
// Copyright 2003-2006, Ronald Fedkiw, Eran Guendelman, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM
//#####################################################################
#ifndef __LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM__
#define __LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<T_GRID>
class LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM:public LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef VECTOR<int,TV::dimension> TV_INT;
    typedef typename T_GRID::UNIFORM_GRID T_UNIFORM_GRID;typedef typename GRID_ARRAYS_POLICY<T_UNIFORM_GRID>::ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;
public:
    template<class T3> struct REBIND{typedef LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM<T_GRID,T3,T_FACE_LOOKUP> TYPE;};
    using LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::Clamped_To_Array;

    LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM()
    {}

    T2 Clamped_To_Array(const T_UNIFORM_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const
    {RANGE<TV_INT> domain=u.Domain_Indices();TV X_new(clamp(X,grid.X(domain.min_corner),grid.X(domain.max_corner)));
    return LINEAR_INTERPOLATION_UNIFORM<T_UNIFORM_GRID,T2>::Clamped_To_Array(grid,u,X_new);}

//#####################################################################
};
}
#endif

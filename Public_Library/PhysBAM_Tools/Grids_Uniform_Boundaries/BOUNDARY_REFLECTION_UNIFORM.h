//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_REFLECTION_UNIFORM
//#####################################################################
#ifndef __BOUNDARY_REFLECTION_UNIFORM__
#define __BOUNDARY_REFLECTION_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T2>
class BOUNDARY_REFLECTION_UNIFORM:public BOUNDARY_UNIFORM<T_GRID,T2>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef VECTOR<bool,T_GRID::dimension> TV_BOOL;
    typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,T_GRID::dimension> TV_SIDES;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;
    typedef typename T_ARRAYS_BASE::template REBIND<T2>::TYPE T_ARRAYS_T2;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
public:
    typedef BOUNDARY_UNIFORM<T_GRID,T2> BASE;
    using BASE::Set_Constant_Extrapolation;using BASE::Constant_Extrapolation;using BASE::Fill_Single_Ghost_Region;

    BOUNDARY_REFLECTION_UNIFORM(const TV_SIDES& constant_extrapolation=TV_SIDES())
    {
        Set_Constant_Extrapolation(constant_extrapolation);
    }

//#####################################################################
    virtual void Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) PHYSBAM_OVERRIDE;
    virtual void Fill_Single_Ghost_Region(const T_GRID& grid,T_ARRAYS_T2& u_ghost,const RANGE<TV_INT>& region,const int side,const T dt,const T time,const int number_of_ghost_cells=3) const PHYSBAM_OVERRIDE ;
    virtual void Apply_Boundary_Condition(const T_GRID& grid,T_ARRAYS_T2& u,const T time) PHYSBAM_OVERRIDE {} // do nothing
//#####################################################################
};
}
#endif

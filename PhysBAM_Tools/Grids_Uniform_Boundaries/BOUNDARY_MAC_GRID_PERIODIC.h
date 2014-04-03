//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_MAC_GRID_PERIODIC
//#####################################################################
#ifndef __BOUNDARY_MAC_GRID_PERIODIC__
#define __BOUNDARY_MAC_GRID_PERIODIC__

#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T2>
class BOUNDARY_MAC_GRID_PERIODIC:public BOUNDARY_UNIFORM<T_GRID,T2>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;typedef typename T_ARRAYS_BASE::template REBIND<T2>::TYPE T_ARRAYS_T2;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename REBIND<T_FACE_ARRAYS_SCALAR,T2>::TYPE T_FACE_ARRAYS_T2;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
public:

    BOUNDARY_MAC_GRID_PERIODIC() 
    {}

//#####################################################################
    void Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) PHYSBAM_OVERRIDE;
    void Fill_Ghost_Cells_Face(const T_GRID& grid,const T_FACE_ARRAYS_T2& u,T_FACE_ARRAYS_T2& u_ghost,const T time,const int number_of_ghost_cells=3) PHYSBAM_OVERRIDE;
    void Apply_Boundary_Condition(const T_GRID& grid,T_ARRAYS_T2& u,const T time) PHYSBAM_OVERRIDE {} // do nothing
//#####################################################################
};
}
#endif

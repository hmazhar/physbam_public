//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ARRAYS_UTILITIES__
#define __ARRAYS_UTILITIES__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T_GRID,class T2>
class ARRAYS_UTILITIES
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::VECTOR_T TV;typedef typename TV::template REBIND<T2>::TYPE TV_T2;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_DIMENSION_T2;
    typedef typename T_ARRAYS_DIMENSION_T2::template REBIND<TV_T2>::TYPE T_ARRAYS_DIMENSION_VECTOR_T2;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename REBIND<T_FACE_ARRAYS_SCALAR,T2>::TYPE T_FACE_ARRAYS_T2;
    typedef typename REBIND<typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR,int>::TYPE T_ARRAYS_INT;typedef typename REBIND<typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR,bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
public:

//#####################################################################
    static void Make_Ghost_Mask_From_Active_Mask(const T_GRID& grid,const T_ARRAYS_BOOL& input_mask,T_ARRAYS_BOOL& output_mask,const int stencil_width,const int number_of_ghost_cells=0);
    static void Compute_Face_Data_From_Cell_Data(const T_GRID& face_grid,T_FACE_ARRAYS_T2& face_array,const T_ARRAYS_DIMENSION_T2& cell_array,const int number_of_ghost_cells=0);
    static void Compute_Gradient_At_Faces_From_Cell_Data(const T_GRID& face_grid,T_FACE_ARRAYS_T2& grad_face_array,const T_ARRAYS_DIMENSION_T2& cell_array,const int number_of_ghost_cells=0);
    static void Compute_Gradient_At_Cells_From_Face_Data(const T_GRID& face_grid,T_ARRAYS_DIMENSION_VECTOR_T2& grad_cell_array,const T_FACE_ARRAYS_T2& face_array,const int number_of_ghost_cells=0);
    static void Compute_Divergence_At_Cells_From_Face_Data(const T_GRID& face_grid,T_ARRAYS_DIMENSION_T2& div_cell_array,const T_FACE_ARRAYS_T2& face_array,const int number_of_ghost_cells=0);
//#####################################################################
};
}
#endif

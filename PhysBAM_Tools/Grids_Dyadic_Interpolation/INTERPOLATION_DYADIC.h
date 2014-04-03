//#####################################################################
// Copyright 2003-2005, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERPOLATION_DYADIC
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __INTERPOLATION_DYADIC__
#define __INTERPOLATION_DYADIC__

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Dyadic/BLOCK_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/FACE_LOOKUP_DYADIC.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_FORWARD.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_DYADIC<T_GRID>
class INTERPOLATION_DYADIC:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::CELL T_CELL;
    typedef typename T_GRID::INDEX INDEX;
public:
    template<class T3> struct REBIND{typedef INTERPOLATION_DYADIC<T_GRID,T3,T_FACE_LOOKUP> TYPE;};

    virtual ~INTERPOLATION_DYADIC()
    {}

    T2 Clamped_To_Array_Node(const T_GRID& grid,const ARRAY<T2>& u,const TV& X) const
    {return From_Leaf_Cell_Node(grid,grid.Leaf_Cell(X),u,X);}

    T2 Clamped_To_Array_Cell(const T_GRID& grid,const ARRAY<T2>& u,const ARRAY<T2>* u_node,const TV& X) const
    {const T_CELL* base_cell=grid.Base_Cell(X);if(base_cell) return From_Block_Cell(grid,BLOCK_DYADIC<T_GRID>(grid,base_cell),u,X);
    const T_CELL* leaf_cell=grid.Leaf_Cell(X);assert(leaf_cell&&u_node);return From_Leaf_Cell_Cell(grid,leaf_cell,u,*u_node,X);}

    TV Clamped_To_Array_Face(const T_GRID& grid,const ARRAY<T>& u,const ARRAY<TV>* u_node,const TV& X) const
    {const T_CELL* base_cell=grid.Base_Cell(X);if(base_cell) return From_Block_Face(grid,BLOCK_DYADIC<T_GRID>(grid,base_cell),T_FACE_LOOKUP(u),X);
    const T_CELL* leaf_cell=grid.Leaf_Cell(X);assert(leaf_cell&&u_node);return From_Leaf_Cell_Face(grid,leaf_cell,u,*u_node,X);}

    T Clamped_To_Array_Face_Component(const int axis,const T_GRID& grid,const ARRAY<TV>& u,const ARRAY<TV>* u_node,const TV& X) const
    {const T_CELL* base_cell=grid.Base_Cell(X);if(base_cell) return From_Block_Face_Component(axis,grid,BLOCK_DYADIC<T_GRID>(grid,base_cell),T_FACE_LOOKUP(u),X);
    const T_CELL* leaf_cell=grid.Leaf_Cell(X);assert(leaf_cell&&u_node);return From_Leaf_Cell_Face_Component(axis,grid,leaf_cell,u,*u_node,X);}

    T2 From_Close_Cell_Node(const T_GRID& grid,const T_CELL* cell,const ARRAY<T2>& u,const TV& X) const
    {return From_Leaf_Cell_Node(grid,grid.Leaf_Cell_By_Neighbor_Path(cell,X),u,X);}

    T2 From_Close_Cell_Cell(const T_GRID& grid,const T_CELL* cell,const ARRAY<T2>& u,const ARRAY<T2>* u_node,const TV& X) const
    {const T_CELL* base_cell=grid.Base_Cell_By_Neighbor_Path(cell,X);if(base_cell) return From_Block_Cell(grid,BLOCK_DYADIC<T_GRID>(grid,base_cell),u,X);
    const T_CELL* leaf_cell=grid.Leaf_Cell_By_Neighbor_Path(cell,X);assert(leaf_cell&&u_node);return From_Leaf_Cell_Cell(grid,leaf_cell,u,*u_node,X);}

    TV From_Close_Cell_Face(const T_GRID& grid,const T_CELL* cell,const ARRAY<T>& u,const ARRAY<TV>* u_node,const TV& X) const
    {const T_CELL* base_cell=grid.Base_Cell_By_Neighbor_Path(cell,X);if(base_cell) return From_Block_Face(grid,BLOCK_DYADIC<T_GRID>(grid,base_cell),u,X);
    const T_CELL* leaf_cell=grid.Leaf_Cell_By_Neighbor_Path(cell,X);assert(leaf_cell&&u_node);return From_Leaf_Cell_Face(grid,leaf_cell,u,*u_node,X);}

//#####################################################################
    virtual ARRAY<PAIR<INDEX,T> > Clamped_To_Array_Cell_Weights(const T_GRID& grid,const ARRAY<T2>& u,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T2 From_Leaf_Cell_Node(const T_GRID& grid,const T_CELL* leaf_cell,const ARRAY<T2>& u,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual TV From_Leaf_Cell_Face(const T_GRID& grid,const T_CELL* leaf_cell,const ARRAY<T>& u,const ARRAY<TV>& u_node,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T From_Leaf_Cell_Face_Component(const int axis,const T_GRID& grid,const T_CELL* leaf_cell,const ARRAY<T>& u,const ARRAY<TV>& u_node,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T2 From_Leaf_Cell_Cell(const T_GRID& grid,const T_CELL* leaf_cell,const ARRAY<T2>& u,const ARRAY<T2>& u_node,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual TV From_Block_Face(const T_GRID& grid,const BLOCK_DYADIC<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T From_Block_Face_Component(const int axis,const T_GRID& grid,const BLOCK_DYADIC<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
        {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T2 From_Block_Cell(const T_GRID& grid,const BLOCK_DYADIC<T_GRID>& block,const ARRAY<T2>& u,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};
}
#endif
#endif

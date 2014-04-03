//#####################################################################
// Copyright 2003-2005, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INTERPOLATION_DYADIC
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __LINEAR_INTERPOLATION_DYADIC__
#define __LINEAR_INTERPOLATION_DYADIC__

#include <PhysBAM_Tools/Grids_Dyadic/BLOCK_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/FACE_LOOKUP_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_BINTREE_HELPER.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_OCTREE_HELPER.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_QUADTREE_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_1D_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_2D_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_3D_HELPER.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_DYADIC<T_GRID>
class LINEAR_INTERPOLATION_DYADIC:public INTERPOLATION_DYADIC<T_GRID,T2,T_FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::CELL T_CELL;typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_HELPER T_LINEAR_INTERPOLATION_HELPER;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_MAC_HELPER T_LINEAR_INTERPOLATION_MAC_HELPER;
public:
    template<class T3> struct REBIND{typedef LINEAR_INTERPOLATION_DYADIC<T_GRID,T3,T_FACE_LOOKUP> TYPE;};

    T2 From_Leaf_Cell_Node(const T_GRID& grid,const T_CELL* cell,const ARRAY<T2>& u,const TV& X) const PHYSBAM_OVERRIDE
    {return T_LINEAR_INTERPOLATION_HELPER::Interpolate_Nodes(grid,cell,u,X);}

    TV From_Leaf_Cell_Face(const T_GRID& grid,const T_CELL* cell,const ARRAY<T>& u,const ARRAY<TV>& u_node,const TV& X) const PHYSBAM_OVERRIDE
    {return T_LINEAR_INTERPOLATION_HELPER::Interpolate_Faces(grid,cell,u,u_node,X);}

    T From_Leaf_Cell_Face_Component(const int axis,const T_GRID& grid,const T_CELL* cell,const ARRAY<T>& u,const ARRAY<TV>& u_node,const TV& X) const PHYSBAM_OVERRIDE
    {return T_LINEAR_INTERPOLATION_HELPER::Interpolate_Face(grid,axis,cell,u,u_node,X);}

    T2 From_Leaf_Cell_Cell(const T_GRID& grid,const T_CELL* cell,const ARRAY<T2>& u,const ARRAY<T2>& u_node,const TV& X) const PHYSBAM_OVERRIDE
    {return T_LINEAR_INTERPOLATION_HELPER::Interpolate_Cells(grid,cell,u,u_node,X);}

    TV From_Block_Face(const T_GRID& grid,const BLOCK_DYADIC<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
    {return T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Face(block,u,X);}

    T From_Block_Face_Component(const int axis,const T_GRID& grid,const BLOCK_DYADIC<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const PHYSBAM_OVERRIDE
    {return T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Face_Component(axis+1,block,u,X);}

    T2 From_Block_Cell(const T_GRID& grid,const BLOCK_DYADIC<T_GRID>& block,const ARRAY<T2>& u,const TV& X) const PHYSBAM_OVERRIDE
    {return T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Cell(block,u,X);}

    ARRAY<PAIR<typename T_GRID::INDEX,T> > Clamped_To_Array_Cell_Weights(const T_GRID& grid,const ARRAY<T2>& u,const TV& X) const PHYSBAM_OVERRIDE
    {return T_LINEAR_INTERPOLATION_HELPER::Interpolate_Cell_Weights(grid,u,X);}
//#####################################################################
};
}
#endif
#endif

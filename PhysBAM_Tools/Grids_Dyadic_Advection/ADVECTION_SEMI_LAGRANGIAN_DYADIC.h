//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_DYADIC
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __ADVECTION_SEMI_LAGRANGIAN_DYADIC__
#define __ADVECTION_SEMI_LAGRANGIAN_DYADIC__


#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/AVERAGING_DYADIC.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> // T_AVERAGING=AVERAGING_DYADIC<T,T,T_GRID>, T_INTERPOLATION=LINEAR_INTERPOLATION_DYADIC<T_GRID,T2>
class ADVECTION_SEMI_LAGRANGIAN_DYADIC:public ADVECTION<T_GRID,T2,typename T_AVERAGING::FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::CELL T_CELL;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_HELPER T_LINEAR_INTERPOLATION_DYADIC_HELPER;typedef typename T_AVERAGING::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;typedef typename REBIND<T_BOUNDARY,T2>::TYPE T_BOUNDARY_T2;
public:
    template<class T3> struct REBIND{typedef ADVECTION_SEMI_LAGRANGIAN_DYADIC<T_GRID,T3,T_AVERAGING,typename T_INTERPOLATION::template REBIND<T3>::TYPE> TYPE;};
    template<class T_INTERPOLATION_2> struct REBIND_INTERPOLATION{typedef ADVECTION_SEMI_LAGRANGIAN_DYADIC<T_GRID,T2,T_AVERAGING,T_INTERPOLATION_2> TYPE;};

    ADVECTION_SEMI_LAGRANGIAN_DYADIC()
    {}

    void Update_Advection_Equation_Node(const T_GRID& grid,ARRAY<T2>& Z,const ARRAY<T2>& Z_ghost,
        const ARRAY<TV>& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const ARRAY<T2>* Z_min_ghost=0,const ARRAY<T2>* Z_max_ghost=0,ARRAY<T2>* Z_min=0,ARRAY<T2>* Z_max=0)
    {assert(!Z_min && !Z_max);
    T_INTERPOLATION interpolation;
    for(DYADIC_GRID_ITERATOR_NODE<T_GRID> iterator(grid,grid.Map_Regular_Nodes());iterator.Valid();iterator.Next()){
        TV sample_location=iterator.Location()-dt*V(iterator.Node_Index());
        Z(iterator.Node_Index())=interpolation.From_Close_Cell_Node(grid,iterator.Deepest_Cell(),Z_ghost,sample_location);}}

    void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,ARRAY<T2>& Z,const ARRAY<T2>& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const ARRAY<T2>* Z_min_ghost,const ARRAY<T2>* Z_max_ghost,ARRAY<T2>* Z_min,ARRAY<T2>* Z_max)
    {assert(!Z_min && !Z_max);
    T_INTERPOLATION interpolation;T_AVERAGING averaging;
    ARRAY<T2> Z_node_ghost(grid.number_of_nodes,false);T_LINEAR_INTERPOLATION_DYADIC_HELPER::Interpolate_From_Cells_To_Nodes(grid,Z_ghost,Z_node_ghost);
    for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid,0);iterator.Valid();iterator.Next()){
        TV sample_location=iterator.Location()-dt*averaging.Face_To_Cell_Vector(grid,iterator.Cell_Pointer(),face_velocities);
        Z(iterator.Cell_Index())=interpolation.From_Close_Cell_Cell(grid,iterator.Cell_Pointer(),Z_ghost,&Z_node_ghost,sample_location);}}

    void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,ARRAY<T>& Z,const T_FACE_LOOKUP& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,ARRAY<T>* Z_min,ARRAY<T>* Z_max)
    {assert(!Z_min && !Z_max);
    // TODO: make this use the interpolation class. Not trivial to do so though because of the T2 type and the fact that only one component should be interpolated
    T_AVERAGING averaging;T_INTERPOLATION interpolation;
    // average face values to nodes (these are only used away from the interface/objects)
    ARRAY<TV> Z_node_ghost(grid.number_of_nodes,false);T_LINEAR_INTERPOLATION_DYADIC_HELPER::Interpolate_From_Faces_To_Nodes(grid,Z_ghost.V_face,Z_node_ghost);
    ARRAY<TV> V_node(grid.number_of_nodes,false);T_LINEAR_INTERPOLATION_DYADIC_HELPER::Interpolate_From_Faces_To_Nodes(grid,face_velocities.V_face,V_node);
    for(DYADIC_GRID_ITERATOR_FACE<T_GRID> iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){int face=iterator.Face_Index();
        TV sample_location=iterator.Location()-dt*averaging.Face_To_Face_Vector(grid,face,face_velocities,V_node);
        T_CELL* base_cell=grid.Base_Cell_By_Neighbor_Path(iterator.Deepest_Cell(),sample_location);
        if(base_cell) Z(face)=interpolation.From_Block_Face_Component(iterator.Axis(),grid,BLOCK_DYADIC<T_GRID>(grid,base_cell),Z_ghost.Starting_Point_Face(face),sample_location);
        else{
            T_CELL* leaf_cell=grid.Leaf_Cell_By_Neighbor_Path(iterator.Deepest_Cell(),sample_location);
            Z(face)=interpolation.From_Leaf_Cell_Face_Component(iterator.Axis(),grid,leaf_cell,Z_ghost.V_face,Z_node_ghost,sample_location);}}}

//#####################################################################
};
}
#endif
#endif

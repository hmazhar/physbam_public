//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_DYADIC  
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_DYADIC.h>
//#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_BINTREE_HELPER.h>
//#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
//#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_OCTREE_HELPER.h>
//#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_QUADTREE_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T2> BOUNDARY_DYADIC<T_GRID,T2>::
BOUNDARY_DYADIC()
{}
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T2> BOUNDARY_DYADIC<T_GRID,T2>::
BOUNDARY_DYADIC(const TV_SIDES& constant_extrapolation)
{
    Set_Constant_Extrapolation(constant_extrapolation);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID,class T2> BOUNDARY_DYADIC<T_GRID,T2>::
~BOUNDARY_DYADIC()
{}
//#####################################################################
// Hack
//#####################################################################
// type safe hack for getting around excessive instantiation
template<class T,class T2> struct Hack{static T& Convert(T2& data){PHYSBAM_FATAL_ERROR();}};
template<class T> struct Hack<T,T>{static T& Convert(T& data){return data;}};
//#####################################################################
// Function Fill_Ghost_Cells_Node
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_DYADIC<T_GRID,T2>::
Fill_Ghost_Cells_Node(const T_GRID& grid,const ARRAY<T2>& u,ARRAY<T2>& u_ghost,const T time)
{
    ARRAY<T2>::Copy(u,u_ghost);
    for(int side=1;side<=T_GRID::number_of_faces_per_cell;side++)Fill_Individual_Side_Ghost_Cells_Node(grid,side,u_ghost,time);    
}
//#####################################################################
// Function Fill_Ghost_Cells_Cell
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_DYADIC<T_GRID,T2>::
Fill_Ghost_Cells_Cell(const T_GRID& grid,const ARRAY<T2>& u,ARRAY<T2>& u_ghost,const T time)
{
    ARRAY<T2>::Copy(u,u_ghost);
    for(int side=1;side<=T_GRID::number_of_faces_per_cell;side++)Fill_Individual_Side_Ghost_Cells_Cell(grid,side,u_ghost,time);
}
//#####################################################################
// Function Fill_Ghost_Cells_Face
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_DYADIC<T_GRID,T2>::
Fill_Ghost_Cells_Face(const T_GRID& grid,const ARRAY<T>& u,ARRAY<T>& u_ghost,const T time)
{
    ARRAY<T>::Copy(u,u_ghost);
    for(int side=1;side<=T_GRID::number_of_faces_per_cell;side++)Fill_Individual_Side_Ghost_Cells_Face(grid,side,u_ghost,time);   
}
//#####################################################################
// Function Fill_Individual_Side_Ghost_Cells_Node
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_DYADIC<T_GRID,T2>::
Fill_Individual_Side_Ghost_Cells_Node(const T_GRID& grid,const int side,ARRAY<T2>& u_ghost,const T time)
{
       PHYSBAM_FATAL_ERROR("Dyadic Fill_Individual_Side_Ghost_Cells_Node not yet implemented");
//     int component=(side-1)/2+1;
//     RANGE<TV> ghost_domain(grid.uniform_grid.Ghost_Domain(grid.number_of_ghost_cells)),domain(grid.uniform_grid.domain);ghost_domain.Change_Size((T)-1e-3*grid.Minimum_Edge_Length());
//     TV clamp_box_min=ghost_domain.Minimum_Corner(),clamp_box_max=ghost_domain.Maximum_Corner();
//     clamp_box_max[component]=domain.Maximum_Corner()[component];clamp_box_min[component]=domain.Minimum_Corner()[component];
//     RANGE<TV> clamp_box(clamp_box_min,clamp_box_max);
// 
//     if(use_fixed_boundary) for(DYADIC_GRID_ITERATOR_NODE<T_GRID> iterator(grid,grid.Map_Individual_Side_Ghost_Nodes(side));iterator.Valid();iterator.Next()) 
//         u_ghost(iterator.Node_Index())=fixed_boundary_value;
//     else if(clamp_below&&clamp_above) for(DYADIC_GRID_ITERATOR_NODE<T_GRID> iterator(grid,grid.Map_Individual_Side_Ghost_Nodes(side));iterator.Valid();iterator.Next()){
//         T2 boundary_value=T_LINEAR_INTERPOLATION_HELPER::Interpolate_Nodes(grid,u_ghost,clamp_box.Clamp(iterator.Location()));
//         u_ghost(iterator.Node_Index())=clamp(boundary_value,lower_threshold,upper_threshold);}
//     else if(clamp_below) for(DYADIC_GRID_ITERATOR_NODE<T_GRID> iterator(grid,grid.Map_Individual_Side_Ghost_Nodes(side));iterator.Valid();iterator.Next()){
//         T2 boundary_value=T_LINEAR_INTERPOLATION_HELPER::Interpolate_Nodes(grid,u_ghost,clamp_box.Clamp(iterator.Location()));
//         u_ghost(iterator.Node_Index())=clamp_min(boundary_value,lower_threshold);}
//     else if(clamp_above) for(DYADIC_GRID_ITERATOR_NODE<T_GRID> iterator(grid,grid.Map_Individual_Side_Ghost_Nodes(side));iterator.Valid();iterator.Next()){
//         T2 boundary_value=T_LINEAR_INTERPOLATION_HELPER::Interpolate_Nodes(grid,u_ghost,clamp_box.Clamp(iterator.Location()));
//         u_ghost(iterator.Node_Index())=clamp_max(boundary_value,upper_threshold);}
//     else for(DYADIC_GRID_ITERATOR_NODE<T_GRID> iterator(grid,grid.Map_Individual_Side_Ghost_Nodes(side));iterator.Valid();iterator.Next()){
//         T2 boundary_value=T_LINEAR_INTERPOLATION_HELPER::Interpolate_Nodes(grid,u_ghost,clamp_box.Clamp(iterator.Location()));
//         u_ghost(iterator.Node_Index())=boundary_value;}
}
//#####################################################################
// Function Fill_Individual_Side_Ghost_Cells_Cell
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_DYADIC<T_GRID,T2>::
Fill_Individual_Side_Ghost_Cells_Cell(const T_GRID& grid,const int side,ARRAY<T2>& u_ghost,const T time)
{
    if(use_fixed_boundary){for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid,grid.number_of_ghost_cells,T_UNIFORM_GRID::GHOST_REGION,side);iterator.Valid();iterator.Next()) 
        u_ghost(iterator.Cell_Index())=fixed_boundary_value;}
    else{
        int inside_cell=side&1;ARRAY<VECTOR<T_CELL*,T_GRID::number_of_neighbors_per_cell> >& neighbors=grid.Neighbors();
        for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid,grid.number_of_ghost_cells,T_UNIFORM_GRID::GHOST_REGION,side);iterator.Valid();iterator.Next())u_ghost(iterator.Cell_Index())=T2();
        for(DYADIC_GRID_ITERATOR_FACE<T_GRID> iterator(grid,grid.Map_Individual_Domain_Faces(side));iterator.Valid();iterator.Next()){
            T_CELL* cell=iterator.Cell(inside_cell);T size=cell->Face_Size();T2 value=u_ghost(cell->Cell());
            while(neighbors(cell->Cell())(side)){cell=neighbors(cell->Cell())(side);assert(!cell->Has_Children());u_ghost(cell->Cell())+=(size/cell->Face_Size())*value;}}}
    if(clamp_below&&clamp_above)for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid,grid.number_of_ghost_cells,T_UNIFORM_GRID::GHOST_REGION,side);iterator.Valid();iterator.Next()){
        u_ghost(iterator.Cell_Index())=clamp(u_ghost(iterator.Cell_Index()),lower_threshold,upper_threshold);}
    else if(clamp_below) for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid,grid.number_of_ghost_cells,T_UNIFORM_GRID::GHOST_REGION,side);iterator.Valid();iterator.Next()){
        u_ghost(iterator.Cell_Index())=clamp_min(u_ghost(iterator.Cell_Index()),lower_threshold);}
    else if(clamp_above) for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid,grid.number_of_ghost_cells,T_UNIFORM_GRID::GHOST_REGION,side);iterator.Valid();iterator.Next()){
        u_ghost(iterator.Cell_Index())=clamp_max(u_ghost(iterator.Cell_Index()),upper_threshold);}
}
//#####################################################################
// Function Fill_Individual_Side_Ghost_Cells_Face
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_DYADIC<T_GRID,T2>::
Fill_Individual_Side_Ghost_Cells_Face(const T_GRID& grid,const int side,ARRAY<T>& u_ghost,const T time)
{
    TV lower,upper,fixed;
    if(use_fixed_boundary) fixed=Hack<TV,T2>::Convert(fixed_boundary_value);
    if(clamp_below) lower=Hack<TV,T2>::Convert(lower_threshold);if(clamp_above) upper=Hack<TV,T2>::Convert(upper_threshold);
    if(use_fixed_boundary){for(DYADIC_GRID_ITERATOR_FACE<T_GRID> iterator(grid,grid.Map_Individual_Side_Ghost_Faces(side));iterator.Valid();iterator.Next()) 
        u_ghost(iterator.Face_Index())=fixed[iterator.Axis()+1];}
    else{
        int axis=(side+1)/2,inside_cell=side&1;T tolerance=(T).25*grid.Minimum_Edge_Length();ARRAY<VECTOR<T_CELL*,T_GRID::number_of_neighbors_per_cell> >& neighbors=grid.Neighbors();
        for(DYADIC_GRID_ITERATOR_FACE<T_GRID> iterator(grid,grid.Map_Individual_Side_Ghost_Faces(side));iterator.Valid();iterator.Next())u_ghost(iterator.Face_Index())=0;
        for(DYADIC_GRID_ITERATOR_FACE<T_GRID> iterator(grid,grid.Map_Individual_Domain_Faces(side));iterator.Valid();iterator.Next()){
            T_CELL* cell=iterator.Cell(inside_cell);T size=cell->Face_Size();T dx=cell->DX().x;
            T parallel_value=u_ghost(cell->Face(side-1));
            T locations[T_GRID::number_of_faces_per_cell-2],values[T_GRID::number_of_faces_per_cell-2];
            static const int axis_to_other_axis[3][2]={{2,3},{1,3},{1,2}};
            for(int i=0;i<T_GRID::dimension-1;i++){int other_axis=axis_to_other_axis[axis-1][i];
                int index1=2*i,index2=2*i+1,side1=2*other_axis-1,side2=2*other_axis;
                locations[index1]=grid.Face_Location(side1-1,cell)[other_axis];
                locations[index2]=grid.Face_Location(side2-1,cell)[other_axis];
                T_CELL *neighbor1=neighbors(cell->Cell())(side1),*neighbor2=neighbors(cell->Cell())(side2);
                values[index1]=!neighbor1||!neighbor1->Has_Children()?u_ghost(cell->Face(side1-1)):0;
                values[index2]=!neighbor2||neighbor2->Depth_Of_This_Cell()<cell->Depth_Of_This_Cell()?u_ghost(cell->Face(side2-1)):0;}
            while(neighbors(cell->Cell())(side)){
                cell=neighbors(cell->Cell())(side);assert(!cell->Has_Children());
                u_ghost(cell->Face(side-1))+=(size/cell->Face_Size())*parallel_value;
                T weight=T_GRID::dimension==3?dx/cell->DX().x:1;
                for(int i=0;i<T_GRID::dimension-1;i++){int other_axis=axis_to_other_axis[axis-1][i];
                    for(int j=0;j<=1;j++){int index=2*i+j,face_index=2*other_axis-2+j;
                        if(abs(grid.Face_Location(face_index,cell)[other_axis]-locations[index])<tolerance) u_ghost(cell->Face(face_index))+=weight*values[index];}}}}}
    if(clamp_below&&clamp_above) for(DYADIC_GRID_ITERATOR_FACE<T_GRID> iterator(grid,grid.Map_Individual_Side_Ghost_Faces(side));iterator.Valid();iterator.Next()){
        u_ghost(iterator.Face_Index())=clamp(u_ghost(iterator.Face_Index()),lower[iterator.Axis()+1],upper[iterator.Axis()+1]);}
    else if(clamp_below) for(DYADIC_GRID_ITERATOR_FACE<T_GRID> iterator(grid,grid.Map_Individual_Side_Ghost_Faces(side));iterator.Valid();iterator.Next()){
        u_ghost(iterator.Face_Index())=clamp_min(u_ghost(iterator.Face_Index()),lower[iterator.Axis()+1]);}
    else if(clamp_above) for(DYADIC_GRID_ITERATOR_FACE<T_GRID> iterator(grid,grid.Map_Individual_Side_Ghost_Faces(side));iterator.Valid();iterator.Next()){
        u_ghost(iterator.Face_Index())=clamp_max(u_ghost(iterator.Face_Index()),upper[iterator.Axis()+1]);}
}
//#####################################################################
#define INSTANTIATION_HELPER(T_GRID) \
    template class BOUNDARY_DYADIC<T_GRID,T_GRID::SCALAR>; \
    template class BOUNDARY_DYADIC<T_GRID,T_GRID::VECTOR_T>; \
    template class BOUNDARY_DYADIC<T_GRID,MATRIX_POLICY<T_GRID::VECTOR_T>::SYMMETRIC_MATRIX>;

#if COMPILE_WITH_BINTREE_SUPPORT
INSTANTIATION_HELPER(BINTREE_GRID<float>)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(BINTREE_GRID<double>)
#endif
#endif

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
INSTANTIATION_HELPER(OCTREE_GRID<float>)
INSTANTIATION_HELPER(QUADTREE_GRID<float>)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(OCTREE_GRID<double>)
INSTANTIATION_HELPER(QUADTREE_GRID<double>)
#endif
#endif

#endif

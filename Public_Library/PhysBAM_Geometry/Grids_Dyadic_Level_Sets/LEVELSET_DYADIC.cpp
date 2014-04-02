//#####################################################################
// Copyright 2003-2005, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
// #if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_CHILDREN.h>
#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_BINTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_CHILDREN.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_CHILDREN.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/FACE_LOOKUP_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Collisions/GRID_BASED_COLLISION_GEOMETRY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/FAST_MARCHING_METHOD_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_DYADIC.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
template<class T,class T_GRID> typename LEVELSET<T,T_GRID>::T_LINEAR_INTERPOLATION_SCALAR LEVELSET<T,T_GRID>::interpolation_default;
template<class T,class T_GRID> typename REBIND<typename LEVELSET<T,T_GRID>::T_LINEAR_INTERPOLATION_SCALAR,typename T_GRID::VECTOR_T>::TYPE LEVELSET<T,T_GRID>::normal_interpolation_default;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> LEVELSET_DYADIC<T_GRID>::
LEVELSET_DYADIC(T_GRID& grid_input,ARRAY<T>& phi_input)
    :grid(grid_input),phi(phi_input),normals(0),curvature(0)
{
    Set_Band_Width();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> LEVELSET_DYADIC<T_GRID>::
~LEVELSET_DYADIC()
{
    delete normals;delete curvature;
}
//#####################################################################
// Function Collision_Aware_Phi
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_DYADIC<T_GRID>::
Collision_Aware_Phi(const TV& location) const
{
    assert(collision_body_list);
    return LINEAR_INTERPOLATION_COLLIDABLE_CELL_DYADIC<T_GRID,T>(*collision_body_list,&valid_mask_current,collidable_phi_replacement_value).Clamped_To_Array_Cell(grid,phi,0,location);
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_DYADIC<T_GRID>::
CFL(const ARRAY<T>& face_velocities) const
{
    T max_V_norm=0;
    for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid);iterator.Valid();iterator.Next()){int cell=iterator.Cell_Index();
        T local_V_norm=0;
        for(int axis=1;axis<=T_GRID::dimension;axis++)
            local_V_norm+=maxabs(face_velocities(grid.First_Face_Index_In_Cell(axis-1,cell)),face_velocities(grid.Second_Face_Index_In_Cell(axis-1,cell)));
        max_V_norm=max(max_V_norm,local_V_norm);}
    T dt_convection=max_V_norm/grid.Minimum_Edge_Length();
    T dt_curvature=curvature_motion?sigma*T_GRID::number_of_faces_per_cell/sqr(grid.Minimum_Edge_Length()):0;
    T dt_overall=dt_convection+dt_curvature;
    return 1/max(dt_overall,1/max_time_step);
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
template<class T_GRID> void LEVELSET_DYADIC<T_GRID>::
Compute_Normals(const T time)
{
    grid.Fully_Refined_Block();
    ARRAY<T> phi_ghost(grid.number_of_cells);boundary->Fill_Ghost_Cells_Cell(grid,phi,phi_ghost,time);
    ARRAY<T> phi_nodes_ghost(grid.number_of_nodes);LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID>::Interpolate_From_Cells_To_Nodes(grid,phi_ghost,phi_nodes_ghost);
    if(!normals) normals=new ARRAY<TV>(grid.number_of_cells,false);else normals->Resize(grid.number_of_cells,false,false);
    for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid,2);iterator.Valid();iterator.Next()){
        TV location=iterator.Location();
        for(int axis=1;axis<=T_GRID::dimension;axis++){
            TV offset;offset[axis]=iterator.DX()[axis];
            (*normals)(iterator.Cell_Index())[axis]=(interpolation->From_Close_Cell_Cell(grid,iterator.Cell_Pointer(),phi_ghost,&phi_nodes_ghost,location+offset)-
                                                      interpolation->From_Close_Cell_Cell(grid,iterator.Cell_Pointer(),phi_ghost,&phi_nodes_ghost,location-offset))/(2*offset[axis]);}
        (*normals)(iterator.Cell_Index()).Normalized();}
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T_GRID> typename T_GRID::VECTOR_T LEVELSET_DYADIC<T_GRID>::
Normal(const TV& location,const bool use_precomputed_normals_if_possible,const ARRAY<T>* phi_nodes) const
{
    const CELL* base_cell=grid.Base_Cell(location);if(base_cell) return Normal(BLOCK_DYADIC<T_GRID>(grid,base_cell),location,use_precomputed_normals_if_possible,phi_nodes);
    const CELL* leaf_cell=grid.Leaf_Cell(location);if(leaf_cell) return Normal_From_Leaf_Cell(leaf_cell,location,phi_nodes);
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T_GRID> typename T_GRID::VECTOR_T LEVELSET_DYADIC<T_GRID>::
Normal_From_Leaf_Cell(const CELL* leaf_cell,const TV& location,const ARRAY<T>* phi_nodes) const
{
    TV normal;
    for(int axis=1;axis<=T_GRID::dimension;axis++){
        TV offset;offset[axis]=leaf_cell->DX()[axis];
        normal[axis]=(Phi_From_Close_Cell(leaf_cell,location+offset,phi_nodes)-Phi_From_Close_Cell(leaf_cell,location-offset,phi_nodes))/(2*offset[axis]);}
    return normal.Normalized();
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T_GRID> typename T_GRID::VECTOR_T LEVELSET_DYADIC<T_GRID>::
Normal(const BLOCK_DYADIC<T_GRID>& block,const TV& location,const bool use_precomputed_normals_if_possible,const ARRAY<T>* phi_nodes) const
{
    if(normals && use_precomputed_normals_if_possible){
        LINEAR_INTERPOLATION_DYADIC<T_GRID,TV> interpolation;
        return interpolation.From_Block_Cell(grid,block,*normals,location).Normalized();}
    else{
        TV normal;
        for(int axis=1;axis<=T_GRID::dimension;axis++){
            TV offset;offset[axis]=grid.Minimum_Edge_Length();
            normal[axis]=(Phi_From_Close_Cell(block.Base_Cell(),location+offset,phi_nodes)-Phi_From_Close_Cell(block.Base_Cell(),location-offset,phi_nodes))/(2*offset[axis]);}
        return normal.Normalized();}
}
//#####################################################################
// Function Extended_Normal
//#####################################################################
template<class T_GRID> typename T_GRID::VECTOR_T LEVELSET_DYADIC<T_GRID>::
Extended_Normal(const TV& location,const bool use_precomputed_normals_if_possible,const ARRAY<T>* phi_nodes) const
{
    TV clamped_location=grid.Clamp(location);
    const CELL* base_cell=grid.Base_Cell(clamped_location);if(base_cell)return Normal(BLOCK_DYADIC<T_GRID>(grid,base_cell),location,use_precomputed_normals_if_possible,phi_nodes);
    const CELL* leaf_cell=grid.Leaf_Cell(clamped_location);if(leaf_cell)return Normal_From_Leaf_Cell(leaf_cell,location,phi_nodes);
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Extended_Normal
//#####################################################################
template<class T_GRID> typename T_GRID::VECTOR_T LEVELSET_DYADIC<T_GRID>::
Extended_Normal_From_Leaf_Cell(const CELL* leaf_cell,const TV& location,const ARRAY<T>* phi_nodes) const
{
    TV normal;
    for(int axis=1;axis<=T_GRID::dimension;axis++){
        TV offset;offset[axis]=leaf_cell->DX()[axis];
        normal[axis]=(Extended_Phi_From_Close_Cell(leaf_cell,location+offset,phi_nodes)-Extended_Phi_From_Close_Cell(leaf_cell,location-offset,phi_nodes))/(2*offset[axis]);}
    return normal.Normalized();
}
//#####################################################################
// Function Create_Octree_Levelset
//#####################################################################
namespace PhysBAM{template<class T_GRID> struct FIND_CELLS_TO_REFINE_HELPER{T_GRID* grid;ARRAY<typename T_GRID::CELL*>* cells_to_refine;typename T_GRID::SCALAR (*phi_function)(const typename T_GRID::VECTOR_T& location);};}
template<class T_GRID> static void Find_Cells_To_Refine(void* data,const typename T_GRID::CELL* cell)
{
    FIND_CELLS_TO_REFINE_HELPER<T_GRID>* helper=(FIND_CELLS_TO_REFINE_HELPER<T_GRID>*)data;
    typedef typename T_GRID::SCALAR T;
    if(cell->Depth_Of_This_Cell()>=helper->grid->maximum_depth)return;
    T dx=(T).5*cell->Diagonal_Length();
    if(abs(helper->phi_function(helper->grid->Node_Location(cell->Cell())))<dx){
        helper->cells_to_refine->Append((typename T_GRID::CELL*)cell);return;}
}
template<class T_GRID> void LEVELSET_DYADIC<T_GRID>::
Create_Dyadic_Levelset(const int max_depth,const int min_depth,T (*phi_function)(const TV& location),const bool verbose)
{
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(verbose) LOG::cout<<"Creating Dyadic Structure"<<std::endl;
#endif
    while(true){
        ARRAY<CELL*> cells_to_refine;
        FIND_CELLS_TO_REFINE_HELPER<T_GRID> helper;helper.grid=&grid;helper.cells_to_refine=&cells_to_refine;helper.phi_function=phi_function;
        MAP_MESH::Map_Cells(grid.uniform_grid,grid.cells,grid.number_of_ghost_cells,&helper,Find_Cells_To_Refine<T_GRID>);
        if(cells_to_refine.m==0)break;
        for(int i=1;i<=cells_to_refine.m;i++){
            if(cells_to_refine(i)->Has_Children()) PHYSBAM_FATAL_ERROR();
            cells_to_refine(i)->Create_Children(grid.number_of_cells,0,grid.number_of_nodes,0,grid.number_of_faces,0,&grid);}}
    phi.Resize(grid.number_of_cells);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(verbose) LOG::cout<<"Number of nodes="<<phi.m<<std::endl;
#endif
    for(int i=1;i<=grid.number_of_cells;i++)
        phi(i)=phi_function(grid.Cell_Location(i));
}
//#####################################################################
// Function Whitney_Criteria
//#####################################################################
template<class T_GRID> bool LEVELSET_DYADIC<T_GRID>::
Whitney_Criteria(const CELL* cell,const T half_band_width) const
{
    return abs(phi(cell->Cell()))<=(half_band_width+(T).5)*cell->DX().Magnitude();
}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class T_GRID> void LEVELSET_DYADIC<T_GRID>::
Fast_Marching_Method(const T time,const T stopping_distance,const ARRAY<int>* seed_indices)
{       
    Get_Signed_Distance_Using_FMM(phi,time,stopping_distance,seed_indices);
}
//#####################################################################
// Function Get_Signed_Distance_Using_FMM
//#####################################################################
template<class T_GRID> void LEVELSET_DYADIC<T_GRID>::
Get_Signed_Distance_Using_FMM(ARRAY<T>& signed_distance,const T time,const T stopping_distance,const ARRAY<int>* seed_indices)
{       
    ARRAY<T> phi_ghost(grid.number_of_cells);boundary->Fill_Ghost_Cells_Cell(grid,phi,phi_ghost,time); 
    if(refine_fmm_initialization_with_iterative_solver) ARRAY<T>::Copy(phi_ghost,phi);
    FAST_MARCHING_METHOD_DYADIC<T_GRID> fmm(*this);
    fmm.Fast_Marching_Method_Cells(phi_ghost,stopping_distance,seed_indices);
    ARRAY<T>::Copy(phi_ghost,signed_distance);
    boundary->Apply_Boundary_Condition(grid,signed_distance,time);
}
//#####################################################################
// Function Fast_Marching_Method_Outside_Band
//#####################################################################
template<class T_GRID> void LEVELSET_DYADIC<T_GRID>::
Fast_Marching_Method_Outside_Band(const T half_band_width,const T time,const T stopping_distance)
{              
    ARRAY<T> phi_ghost(grid.number_of_cells);boundary->Fill_Ghost_Cells_Cell(grid,phi,phi_ghost,time); 
    if(refine_fmm_initialization_with_iterative_solver) ARRAY<T>::Copy(phi_ghost,phi);
    FAST_MARCHING_METHOD_DYADIC<T_GRID> fmm(*this);
    fmm.Fast_Marching_Method_Cells(phi_ghost,stopping_distance);
    for(int i=1;i<=grid.number_of_nodes;i++) if(abs(phi_ghost(i)) > half_band_width) phi(i)=phi_ghost(i);
    boundary->Apply_Boundary_Condition(grid,phi,time);
}
//#####################################################################
// Function Approximate_Negative_Material
//#####################################################################
// calculates the approximate volume using Heaviside functions
namespace PhysBAM{template<class T_GRID> struct APPROXIMATE_MATERIAL_HELPER{const LEVELSET_DYADIC<T_GRID>* levelset;typename T_GRID::SCALAR *total_size,*interface_half_width;};}
template<class T_GRID> static void Approximate_Negative_Material_Helper(void* data,const typename T_GRID::CELL* cell)
{
    APPROXIMATE_MATERIAL_HELPER<T_GRID>* helper=(APPROXIMATE_MATERIAL_HELPER<T_GRID>*)data;
    typedef typename T_GRID::SCALAR T;
    ARRAY<T>& phi=helper->levelset->phi;T& interface_half_width=*helper->interface_half_width;
    *helper->total_size+=LEVELSET_UTILITIES<T>::Heaviside(-phi(cell->Cell())*T_GRID::one_over_number_of_nodes_per_cell,interface_half_width)*cell->Cell_Size();
}
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_DYADIC<T_GRID>::
Approximate_Negative_Material(const T interface_thickness,const T time) const
{
    T total_size=0;T interface_half_width=interface_thickness*grid.Minimum_Edge_Length()/2;
    APPROXIMATE_MATERIAL_HELPER<T_GRID> helper;helper.levelset=this;helper.total_size=&total_size;helper.interface_half_width=&interface_half_width;
    MAP_MESH::Map_Cells(grid.uniform_grid,grid.cells,0,&helper,Approximate_Negative_Material_Helper<T_GRID>);
    return total_size;
}
//#####################################################################
// Function Approximate_Positive_Material
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LEVELSET_DYADIC<T_GRID>::
Approximate_Positive_Material(const T interface_thickness,const T time) const
{
    return grid.Domain().Size()-Approximate_Negative_Material(interface_thickness,time);
}
//#####################################################################
#if COMPILE_WITH_BINTREE_SUPPORT
template class LEVELSET_DYADIC<BINTREE_GRID<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_DYADIC<BINTREE_GRID<double> >;
#endif
#endif

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template class LEVELSET_DYADIC<OCTREE_GRID<float> >;
template class LEVELSET_DYADIC<QUADTREE_GRID<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_DYADIC<OCTREE_GRID<double> >;
template class LEVELSET_DYADIC<QUADTREE_GRID<double> >;
#endif
#endif

#endif

//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SCATTERED_INTERPOLATION
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Particles_Interpolation/SCATTERED_INTERPOLATION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> SCATTERED_INTERPOLATION<T_GRID>::
SCATTERED_INTERPOLATION()
{
    Set_Radius_Of_Influence(0);
    Use_Tent_Weights();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> SCATTERED_INTERPOLATION<T_GRID>::
~SCATTERED_INTERPOLATION()
{}
//#####################################################################
// Function Transfer_To_Grid
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_T2> void SCATTERED_INTERPOLATION<T_GRID>::
Transfer_To_Grid_Helper(ARRAY_VIEW<const TV> domain,ARRAY_VIEW<const typename T_ARRAYS_T2::ELEMENT> range,const T_GRID& grid,T_ARRAYS_T2& grid_data) const
{
    typedef typename T_ARRAYS_T2::ELEMENT T2;
    T_ARRAYS_ARRAY_INT points_in_cell;Bin_Domain_Values(domain,grid,points_in_cell);
    if(use_tent_weights) Transfer_With_Tent_Weights(domain,range,grid,grid_data,points_in_cell);
    else if(use_distance_averaged_weights) Transfer_With_Distance_Averaged_Weights(domain,range,grid,grid_data,points_in_cell);
}
//#####################################################################
// Function Bin_Domain_Values
//#####################################################################
template<class T_GRID> void SCATTERED_INTERPOLATION<T_GRID>::
Bin_Domain_Values(ARRAY_VIEW<const TV> domain_values,const T_GRID& grid,T_ARRAYS_ARRAY_INT& points_in_cell)
{
    points_in_cell.Resize(grid.Get_MAC_Grid().Domain_Indices());
    for(int k=1;k<=domain_values.Size();k++)points_in_cell(grid.Clamp_To_Cell(domain_values(k))).Append(k);
}
//#####################################################################
// Function Transfer_With_Distance_Averaged_Weights
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_T2> void SCATTERED_INTERPOLATION<T_GRID>::
Transfer_With_Distance_Averaged_Weights(ARRAY_VIEW<const TV> domain,ARRAY_VIEW<const typename T_ARRAYS_T2::ELEMENT> range,const T_GRID& grid,T_ARRAYS_T2& grid_data,
    const T_ARRAYS_ARRAY_INT& points_in_cell) const
{
    typedef typename T_ARRAYS_T2::ELEMENT T2;
    assert(!grid.Is_MAC_Grid());
    for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();TV X=iterator.Location();
        grid_data(node)=T2();
        RANGE<TV_INT> index_box=Grid_Influence_Bounds(grid,X);
        int count=0;for(CELL_ITERATOR cell_iterator(grid,index_box);cell_iterator.Valid();cell_iterator.Next()) count+=points_in_cell(cell_iterator.Cell_Index()).m;
        bool found_data=false;ARRAY<T> weights(count);int weight_index=0;T product=(T)1;
        for(CELL_ITERATOR cell_iterator(grid,index_box);cell_iterator.Valid();cell_iterator.Next()){TV_INT cell=cell_iterator.Cell_Index();
            const ARRAY<int>& point_indices=points_in_cell(cell);
            for(int indirect_index=1;indirect_index<=point_indices.m;indirect_index++){
                weight_index++;int k=point_indices(indirect_index);
                T distance_squared=(domain(k)-X).Magnitude_Squared();
                if(distance_squared<=radius_of_influence_squared){
                    T distance=sqrt(distance_squared);found_data=true;weights(weight_index)=product;product*=distance;for(int l=1;l<weight_index;l++) if(weights(l)) weights(l)*=distance;}}}
        if(found_data){
            weight_index=0;T normalization=(T)0;
            for(CELL_ITERATOR cell_iterator(grid,index_box);cell_iterator.Valid();cell_iterator.Next()){TV_INT cell=cell_iterator.Cell_Index();
                const ARRAY<int>& point_indices=points_in_cell(cell);
                for(int indirect_index=1;indirect_index<=point_indices.m;indirect_index++) if(weights(++weight_index) > 0){
                    int k=point_indices(indirect_index);grid_data(node)+=weights(weight_index)*range(k);normalization+=weights(weight_index);}}
            grid_data(node)/=normalization;}}
}
//#####################################################################
// Function Transfer_With_Distance_Averaged_Weights
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_T2> void SCATTERED_INTERPOLATION<T_GRID>::
Transfer_With_Tent_Weights(ARRAY_VIEW<const TV> domain,ARRAY_VIEW<const typename T_ARRAYS_T2::ELEMENT> range,const T_GRID& grid,T_ARRAYS_T2& grid_data,
    const T_ARRAYS_ARRAY_INT& points_in_cell) const
{

    typedef typename T_ARRAYS_T2::ELEMENT T2;
    assert(!grid.Is_MAC_Grid());
    for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();TV X=iterator.Location();
        grid_data(node)=T2();int used_particles=0;
        RANGE<TV_INT> index_box=Grid_Influence_Bounds(grid,X);
        for(CELL_ITERATOR cell_iterator(grid,index_box);cell_iterator.Valid();cell_iterator.Next()){TV_INT cell=cell_iterator.Cell_Index();
            const ARRAY<int>& point_indices=points_in_cell(cell);
            for(int indirect_index=1;indirect_index<=point_indices.m;indirect_index++){
                int k=point_indices(indirect_index);
                T distance_squared=(TV(domain(k))-X).Magnitude_Squared();
                if(distance_squared<=radius_of_influence_squared){
                    T distance=sqrt(distance_squared);T weight=1-distance/radius_of_influence;grid_data(node)+=weight*range(k);used_particles++;}}}
        if(used_particles>0)grid_data(node)/=(T)used_particles;}
}
//#####################################################################
// Function Transfer_With_Distance_Averaged_Weights
//#####################################################################
template<class T_GRID> void SCATTERED_INTERPOLATION<T_GRID>::
Instantiate()
{
    // This is a work around for Visual C++ compiler.
    typename REBIND<T_ARRAYS_SCALAR,VECTOR<T,3> >::TYPE grid_data;
    Transfer_To_Grid_Helper(ARRAY_VIEW<const TV>(0,0),ARRAY_VIEW<const VECTOR<T,3> >(0,0),T_GRID(),grid_data);
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template class SCATTERED_INTERPOLATION<GRID<VECTOR<T,1> > >;   \
    template class SCATTERED_INTERPOLATION<GRID<VECTOR<T,3> > >;
INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif

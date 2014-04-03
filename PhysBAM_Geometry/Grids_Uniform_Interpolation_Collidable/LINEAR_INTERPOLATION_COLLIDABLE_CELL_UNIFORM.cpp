//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Eran Guendelman, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM.h>
using namespace PhysBAM;
template<class T_GRID,class T2> LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<T_GRID,T2>::
LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM(const T_GRID_BASED_COLLISION_GEOMETRY& body_list_input,const T_ARRAYS_BOOL* cell_valid_value_mask_input,
    const T2& default_cell_replacement_value_input,const bool extrapolate_to_invalid_cell_values_input)
    :body_list(body_list_input),cell_valid_value_mask(cell_valid_value_mask_input),default_cell_replacement_value(default_cell_replacement_value_input),
    extrapolate_to_invalid_cell_values(extrapolate_to_invalid_cell_values_input)
{
}
template<class T_GRID,class T2> LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<T_GRID,T2>::
~LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM()
{
}
template<class T_GRID,class T2> T2 LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<T_GRID,T2>::
Trace_And_Get_Value(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X,const TV_INT& cell,TV_INT* cells,int& valid_mask,const int mask_value) const
{
    TV intersection_point;
    COLLISION_GEOMETRY_ID body_id;
    int aggregate_id;
    assert(cell_valid_value_mask);
    if(body_list.Cell_Center_Intersection(cell,cells,T_GRID::number_of_cells_per_block,X,body_id,aggregate_id,intersection_point) || !(*cell_valid_value_mask)(cell)) return default_cell_replacement_value;
    else{valid_mask|=mask_value;return u(cell);}
}
template<class T_GRID,class T2> T2 LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<T_GRID,T2>::
Invalid_Value_Replacement(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X,const TV_INT& cell,const TV& intersection_point,
    const COLLISION_GEOMETRY_ID body_id,const int aggregate_id,bool& valid,const T ray_t_max) const
{
    valid=false;
    return default_cell_replacement_value;
}
template<class T_GRID,class T2> void LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<T_GRID,T2>::
Extrapolate_To_Invalid_Values(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X,const TV_INT& cell,T2* values,int& valid_mask) const
{
    static const int all_valid_mask=(1<<T_GRID::number_of_cells_per_block)-1,neighbor_mask=T_GRID::number_of_cells_per_block-1;
    if(valid_mask==0) for(int t=0;t<T_GRID::number_of_cells_per_block;t++)values[t]=default_cell_replacement_value;
    else if(valid_mask!=all_valid_mask){ // not everything is valid
        for(int t=0;t<T_GRID::number_of_cells_per_block;t++)if(!(valid_mask&(1<<t))){
            int ring;
            for(ring=1;ring<T_GRID::dimension;ring++){ // one ring in 2d, one and then two ring in 3d
                int ring_mask=ring==2?neighbor_mask:0;
                T2 sum=T2();int count=0;
                for(int n=0;n<T_GRID::dimension;n++){
                    int neighbor=t^(1<<n)^ring_mask;
                    if(valid_mask&(1<<neighbor)){count++;sum+=values[neighbor];}}
                if(count>0){values[t]=sum/(T)count;break;}} // average valid ring neighbors
            if(ring==T_GRID::dimension) values[t]=values[t^neighbor_mask];} // final ring (must be valid because valid_mask!=0)
        valid_mask=all_valid_mask;}
}
template<class T_GRID,class T2> T2 LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<T_GRID,T2>::
From_Base_Node(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X,const TV_INT& cell,bool& interpolation_valid,VECTOR<T2,2>* extrema) const
{
    BLOCK_UNIFORM<T_GRID> block(grid,cell+TV_INT::All_Ones_Vector());
    assert(block.Inside(X));
    if(!body_list.Occupied_Block(block)){interpolation_valid=true;LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2> linear;
        if(extrema) *extrema=linear.Extrema_From_Base_Node(grid,u,u,X,cell);return linear.From_Base_Node(grid,u,X,cell);}
    COLLISION_GEOMETRY_ID body_id(0);int aggregate_id=0;
    if(body_list.Inside_Any_Simplex_Of_Any_Body(X,body_id,aggregate_id)) return Invalid_Value_Replacement(grid,u,X,cell,X,body_id,aggregate_id,interpolation_valid);
    int valid_mask=0;
    T2 values[T_GRID::number_of_nodes_per_cell];TV_INT nodes_in_cell[T_GRID::number_of_nodes_per_cell];grid.Nodes_In_Cell_From_Minimum_Corner_Node(cell,nodes_in_cell);
    TV_INT cells_in_block[T_GRID::number_of_cells_per_block];block.All_Cell_Indices(cells_in_block);
    for(int k=0;k<T_GRID::number_of_nodes_per_cell;k++)values[k]=Trace_And_Get_Value(grid,u,X,nodes_in_cell[k],cells_in_block,valid_mask,1<<k);
    if(extrapolate_to_invalid_cell_values) Extrapolate_To_Invalid_Values(grid,u,X,cell,values,valid_mask);
    interpolation_valid=(valid_mask!=0);
    if(extrema){for(int k=0;k<T_GRID::number_of_nodes_per_cell;k++){extrema->x=Componentwise_Min(values[k],extrema->x);extrema->y=Componentwise_Max(values[k],extrema->y);}return T2();}
    else return LINEAR_INTERPOLATION<T,T2>::Linear(values,(X-grid.X(cell))*grid.one_over_dX);
}
template<class T_GRID,class T2> VECTOR<T2,2> LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<T_GRID,T2>::
Extrema_Clamped_To_Array(const T_GRID& grid,const T_ARRAYS_T2& u_min,const T_ARRAYS_T2& u_max,const TV& X) const
{
    return Extrema_From_Base_Node(grid,u_min,u_max,X,INTERPOLATION_UNIFORM<T_GRID,T2>::Clamped_Index_End_Minus_One(grid,u_min,X));
}
template<class T_GRID,class T2> VECTOR<T2,2> LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<T_GRID,T2>::
Extrema_From_Base_Node(const T_GRID& grid,const T_ARRAYS_T2& u_min,const T_ARRAYS_T2& u_max,const TV& X,const TV_INT& cell) const
{
    VECTOR<T2,2> extrema_min,extrema_max;
    bool unused;
    From_Base_Node(grid,u_min,X,cell,unused,&extrema_min);
    From_Base_Node(grid,u_min,X,cell,unused,&extrema_max);
    return VECTOR<T2,2>(extrema_min.x,extrema_max.y);
}
template<class T_GRID,class T2> T2 LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<T_GRID,T2>::
Clamped_To_Array(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const
{
    return From_Base_Node(grid,u,X,INTERPOLATION_UNIFORM<T_GRID,T2>::Clamped_Index_End_Minus_One(grid,u,X));
}
// TODO: this only works for cells, which would ideally be enforced at compile time
template<class T_GRID,class T2> T2 LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<T_GRID,T2>::
From_Base_Node(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X,const TV_INT& cell) const
{
    bool interpolation_valid;
    return From_Base_Node(grid,u,X,cell,interpolation_valid);
}
template class LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<VECTOR<float,1> >,float>;
template class LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<VECTOR<float,2> >,float>;
template class LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<VECTOR<float,3> >,float>;
template class LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<VECTOR<float,3> >,VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<VECTOR<double,1> >,double>;
template class LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<VECTOR<double,2> >,double>;
template class LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<VECTOR<double,3> >,double>;
template class LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<VECTOR<double,3> >,VECTOR<double,3> >;
#endif

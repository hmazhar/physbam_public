//#####################################################################
// Copyright 2005, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Collisions/GRID_BASED_COLLISION_GEOMETRY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/FAST_MARCHING_METHOD_DYADIC.h>
#include <PhysBAM_Geometry/Level_Sets/FAST_MARCHING.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> FAST_MARCHING_METHOD_DYADIC<T_GRID>::
FAST_MARCHING_METHOD_DYADIC(const LEVELSET_DYADIC<T_GRID>& levelset)
    :grid(levelset.grid),levelset(levelset),currently_using_cells(false),cell_pointer_from_index(levelset.grid.Cell_Pointer_From_Index()),
     node_neighbors(levelset.grid.Node_Neighbors()),cell_neighbors(levelset.grid.Neighbors()),phi_nodes_for_initialize(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> FAST_MARCHING_METHOD_DYADIC<T_GRID>::
~FAST_MARCHING_METHOD_DYADIC()
{}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class T_GRID> void FAST_MARCHING_METHOD_DYADIC<T_GRID>::
Fast_Marching_Method_Cells(ARRAY<T>& phi_ghost,const T stopping_distance,const ARRAY<int>* seed_indices,const bool add_seed_indices_for_ghost_cells)
{
    assert(!seed_indices); // TODO: make seed indices work with cell based fast marching
    // march the fine cells by the interface
    currently_using_cells=true;levelset.grid.Fully_Refined_Block();
    ARRAY<T> phi_nodes(grid.number_of_nodes);LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID>::Interpolate_From_Cells_To_Nodes(grid,phi_ghost,phi_nodes);phi_nodes_for_initialize=&phi_nodes;
    LEVELSET_DYADIC<T_GRID> levelset_temp(grid,phi_ghost);levelset_for_initialize=&levelset_temp;
    Fast_Marching_Method(phi_ghost,stopping_distance,0);
    // interpolate phi values to the nodes
    phi_nodes.Fill(0);
    ARRAY<T> weight(grid.number_of_nodes);
    ARRAY<int> number_of_valid_cells_touching(grid.number_of_nodes);
    for(CELL_ITERATOR iterator(grid,grid.number_of_ghost_cells);iterator.Valid();iterator.Next()){
        T cell_weight=1/iterator.Cell_Pointer()->DX().x;
        T weight_times_phi=cell_weight*phi_ghost(iterator.Cell_Index());
        for(int i=0;i<T_GRID::number_of_nodes_per_cell;i++){
            phi_nodes(iterator.Cell_Pointer()->Node(i))+=weight_times_phi;
            weight(iterator.Cell_Pointer()->Node(i))+=cell_weight;}
        if(iterator.Cell_Pointer()->Depth_Of_This_Cell()==grid.maximum_depth)
            for(int i=0;i<T_GRID::number_of_nodes_per_cell;i++)number_of_valid_cells_touching(iterator.Cell_Pointer()->Node(i))++;}
    // seed the node part of the march
    ARRAY<int> seed_nodes;seed_nodes.Preallocate(number_of_valid_cells_touching.Count_Matches(T_GRID::number_of_cells_per_node));
    for(NODE_ITERATOR iterator(grid,grid.Map_All_Nodes());iterator.Valid();iterator.Next()){
        phi_nodes(iterator.Node_Index())/=weight(iterator.Node_Index());
        if(number_of_valid_cells_touching(iterator.Node_Index())==T_GRID::number_of_cells_per_node)seed_nodes.Append(iterator.Node_Index());}
    // fast march the nodes
    currently_using_cells=false;
    Fast_Marching_Method(phi_nodes,stopping_distance,&seed_nodes,true);
    // intepolate the nodal values back onto the cells
    for(CELL_ITERATOR iterator(grid,grid.number_of_ghost_cells);iterator.Valid();iterator.Next()){
        if(iterator.Cell_Pointer()->Depth_Of_This_Cell()<grid.maximum_depth){
            phi_ghost(iterator.Cell_Index())=0;
            for(int i=0;i<T_GRID::number_of_nodes_per_cell;i++)phi_ghost(iterator.Cell_Index())+=phi_nodes(iterator.Cell_Pointer()->Node(i));
            phi_ghost(iterator.Cell_Index())*=T_GRID::one_over_number_of_nodes_per_cell;}}
}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class T_GRID> void FAST_MARCHING_METHOD_DYADIC<T_GRID>::
Fast_Marching_Method_Nodes(ARRAY<T>& phi_ghost_nodes,const T stopping_distance,const ARRAY<int>* seed_indices,const bool add_seed_indices_for_ghost_cells,ARRAY<T>* phi_ghost_cells)
{
    currently_using_cells=false;levelset.grid.Fully_Refined_Block();
    ARRAY<T>* phi_cells=phi_ghost_cells;
    if(!phi_cells){
        phi_cells=new ARRAY<T>(grid.number_of_cells);
        LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID>::Interpolate_From_Nodes_To_Cells(grid,phi_ghost_nodes,*phi_cells);phi_nodes_for_initialize=&phi_ghost_nodes;}
    LEVELSET_DYADIC<T_GRID> levelset_temp(grid,*phi_cells);levelset_for_initialize=&levelset_temp;phi_nodes_for_initialize=&phi_ghost_nodes;
    Fast_Marching_Method(phi_ghost_nodes,stopping_distance,seed_indices,add_seed_indices_for_ghost_cells);
    if(!phi_ghost_cells) delete phi_cells;
}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class T_GRID> void FAST_MARCHING_METHOD_DYADIC<T_GRID>::
Fast_Marching_Method(ARRAY<T>& phi_ghost,const T stopping_distance,const ARRAY<int>* seed_indices,const bool add_seed_indices_for_ghost_cells)
{
    int heap_length=0;
    ARRAY<bool,VECTOR<int,1> > done(0,Number_Of_Computational_Nodes()); // starts at 0 so that the 0th node is never done - Update_Close_Point() relies on this
    ARRAY<int> close_k(Number_Of_Computational_Nodes());
    ARRAY<int> heap(Number_Of_Computational_Nodes(),false);

    Initialize_Interface(phi_ghost,done,close_k,heap,heap_length,seed_indices,add_seed_indices_for_ghost_cells);

    assert(T_GRID::number_of_faces_per_cell==T_GRID::number_of_neighbors_per_node); // this has to be true for node and cell fast marching to work
    while(heap_length != 0){
        int index=heap(1); // smallest point is on top of heap
        if(stopping_distance && abs(phi_ghost(index)) > stopping_distance){ // exit early
            for(int index=1;index<=Number_Of_Computational_Nodes();index++)if(Use_Computational_Node(index))if(!done(index)){
                phi_ghost(index)=LEVELSET_UTILITIES<T>::Sign(phi_ghost(index))*stopping_distance;}
            return;}
        done(index)=true;close_k(index)=0; // add to done, remove from close
        FAST_MARCHING<T>::Down_Heap(phi_ghost,close_k,heap,heap_length);heap_length--; // remove point from heap

        if(levelset.collision_body_list){
            for(int i=1;i<=T_GRID::number_of_faces_per_cell;i++){
                int neighbor_index=Neighbor(i,index);
                if(neighbor_index!=0 && !done(neighbor_index) && Neighbor_Visible(i,index)) Update_Or_Add_Neighbor(phi_ghost,done,close_k,heap,heap_length,neighbor_index);}}
        else
            for(int i=1;i<=T_GRID::number_of_faces_per_cell;i++){
                int neighbor_index=Neighbor(i,index);
                if(neighbor_index!=0 && !done(neighbor_index)) Update_Or_Add_Neighbor(phi_ghost,done,close_k,heap,heap_length,neighbor_index);}}
}
//#####################################################################
// Function Update_Or_Add_Neighbor
//#####################################################################
template<class T_GRID> inline void FAST_MARCHING_METHOD_DYADIC<T_GRID>::
Update_Or_Add_Neighbor(ARRAY<T>& phi_ghost,ARRAY<bool,VECTOR<int,1> >& done,ARRAY<int>& close_k,ARRAY<int>& heap,int& heap_length,int neighbor)
{
    if(close_k(neighbor)){
        Update_Close_Point(phi_ghost,done,neighbor);
        FAST_MARCHING<T>::Up_Heap(phi_ghost,close_k,heap,close_k(neighbor));}
    else{close_k(neighbor)=1; // add to close 
        Update_Close_Point(phi_ghost,done,neighbor);
        heap_length++;heap(heap_length)=neighbor;
        FAST_MARCHING<T>::Up_Heap(phi_ghost,close_k,heap,heap_length);}
}
//#####################################################################
// Function Initialize_Interface
//#####################################################################
// pass heap_length by reference
template<class T_GRID> void FAST_MARCHING_METHOD_DYADIC<T_GRID>::
Initialize_Interface(ARRAY<T>& phi,ARRAY<bool,VECTOR<int,1> >& done,ARRAY<int>& close_k,ARRAY<int>& heap,int& heap_length,const ARRAY<int>* seed_indices,const bool add_seed_indices_for_ghost_cells)
{
    if(seed_indices){
        for(int i=1;i<=seed_indices->m;i++)Add_To_Initial(done,close_k,(*seed_indices)(i));
        if(add_seed_indices_for_ghost_cells){
            assert(!currently_using_cells); // this only works for nodes right now
            for(int side=1;side<=T_GRID::number_of_faces_per_cell;side++){
                for(NODE_ITERATOR iterator(grid,grid.Map_Individual_Side_Ghost_And_Boundary_Nodes(side));iterator.Valid();iterator.Next()){
                    for(int axis=1;axis<=T_GRID::dimension;axis++){
                        int node2=Neighbor(2*axis,iterator.Node_Index());
                        if(node2&&LEVELSET_UTILITIES<T>::Interface(phi(iterator.Node_Index()),phi(node2))){
                            if(!done(iterator.Node_Index()))Add_To_Initial(done,close_k,iterator.Node_Index());
                            if(!done(node2))Add_To_Initial(done,close_k,node2);}}}}}}
    else{
        T sqr_epsilon_dx_dy=sqr(levelset.small_number*sqr(grid.Minimum_Edge_Length())),sqr_epsilon_dx_dy_dz=sqr_epsilon_dx_dy*sqr(grid.Minimum_Edge_Length());
        ARRAY<T> phi_new(Number_Of_Computational_Nodes(),false);
        phi_new.Fill(2*grid.Minimum_Edge_Length()); // ok positive, since minmag is used below

        if(levelset.collision_body_list){
            for(int i=1;i<=Number_Of_Computational_Nodes();i++)if(Use_Computational_Node(i))for(int axis=1;axis<=T_GRID::dimension;axis++){
                int node2=Neighbor(2*axis,i);if(node2){
                    if(!Neighbor_Visible(2*axis,i)){
                        if(phi(i)<=0) Add_To_Initial(done,close_k,i);
                        if(phi(node2)<=0) Add_To_Initial(done,close_k,node2);}
                    else if(LEVELSET_UTILITIES<T>::Interface(phi(i),phi(node2))){
                        Add_To_Initial(done,close_k,i);Add_To_Initial(done,close_k,node2);}}}}
        else{
            for(int i=1;i<=Number_Of_Computational_Nodes();i++)if(Use_Computational_Node(i))for(int axis=1;axis<=T_GRID::dimension;axis++){
                int node2=Neighbor(2*axis,i);if(node2){
                    if(LEVELSET_UTILITIES<T>::Interface(phi(i),phi(node2))){
                        Add_To_Initial(done,close_k,i);Add_To_Initial(done,close_k,node2);}}}}
        
        T fmm_initialization_iterative_tolerance_absolute=levelset.fmm_initialization_iterative_tolerance*grid.Minimum_Edge_Length();

        if(levelset.collision_body_list){
            COLLISION_GEOMETRY_ID body_id;
            for(int index=1;index<=Number_Of_Computational_Nodes();index++)if(Use_Computational_Node(index))if(done(index)){
                T value[T_GRID::dimension]; // the phi value to use in the given direction
                int number_of_axis=0; // the number of axis that we want to use later
                bool really_clamp_phi_with_collision_bodies=levelset.clamp_phi_with_collision_bodies&&phi(index)<=0;
                T abs_phi=abs(phi(index));
                TV location=Computational_Node_Location(index);
                for(int axis=1;axis<=T_GRID::dimension;axis++){
                    int low=Neighbor(2*axis-1,index),high=Neighbor(2*axis,index);
                    bool use_low=false,use_high=false;T value_low=0,value_high=0;
                    T dx_low=low?location[axis]-Computational_Node_Location(low)[axis]:0,dx_high=high?Computational_Node_Location(high)[axis]-location[axis]:0;
                    if(low){
                        if(!Neighbor_Visible(2*axis-1,index)){
                            RAY<TV> ray(location,-TV::Axis_Vector(axis),true);ray.t_max=dx_low;ray.semi_infinite=false;
                            levelset.collision_body_list->Closest_Non_Intersecting_Point_Of_Any_Body(ray,body_id);
                            T phi_object=max(ray.t_max,levelset.small_number);
                            use_low=true;if(really_clamp_phi_with_collision_bodies) value_low=min(abs_phi,phi_object);else value_low=phi_object;}
                        else if(done(low)&&LEVELSET_UTILITIES<T>::Interface(phi(index),phi(low))){
                            use_low=true;value_low=LEVELSET_UTILITIES<T>::Theta(phi(index),phi(low))*dx_low;}}
                    if(high){
                        if(!Neighbor_Visible(2*axis,index)){
                            RAY<TV> ray(location,TV::Axis_Vector(axis),true);ray.t_max=dx_high;ray.semi_infinite=false;
                            levelset.collision_body_list->Closest_Non_Intersecting_Point_Of_Any_Body(ray,body_id);
                            T phi_object=max(ray.t_max,levelset.small_number);
                            use_high=true;if(really_clamp_phi_with_collision_bodies) value_high=min(abs_phi,phi_object);else value_high=phi_object;}
                        else if(done(high)&&LEVELSET_UTILITIES<T>::Interface(phi(index),phi(high))){
                            use_high=true;value_high=LEVELSET_UTILITIES<T>::Theta(phi(index),phi(high))*dx_high;}}
                    if(!use_low){if(use_high) value[number_of_axis]=value_high;else number_of_axis--;}
                    else if(!use_high) value[number_of_axis]=value_low;
                    else value[number_of_axis]=min(value_low,value_high);
                    number_of_axis++;}
                assert(number_of_axis);
                if(number_of_axis==1) phi_new(index)=value[0];
                else if(number_of_axis==2){T d2=sqr(value[0])+sqr(value[1]);if(d2 > sqr_epsilon_dx_dy) phi_new(index)=value[0]*value[1]/sqrt(d2);else phi_new(index)=min(value[0],value[1]);}
                else{assert(T_GRID::dimension==3); // 2d should never get to this point
                    T value_xy=value[0]*value[1],value_xz=value[0]*value[2],value_yz=value[1]*value[2],d2=sqr(value_xy)+sqr(value_xz)+sqr(value_yz);
                    if(d2 > sqr_epsilon_dx_dy_dz) phi_new(index)=value_xy*value[2]/sqrt(d2);else phi_new(index)=min(value[0],value[1],value[2]);}
                if(levelset.refine_fmm_initialization_with_iterative_solver){
                    TV vertex=location;T phi_value=phi(index);
                    CELL* leaf_cell=grid.Leaf_Cell(vertex);
                    for(int iterations=1;iterations<=levelset.fmm_initialization_iterations;iterations++){
                        if(abs(phi_value)>10*grid.Minimum_Edge_Length())break; // stop if it looks like it's blowing up
                        vertex-=phi_value*Normal(leaf_cell,vertex);
                        leaf_cell=grid.Leaf_Cell_By_Neighbor_Path(leaf_cell,vertex);if(!leaf_cell)break; // we are outside the grid
                        phi_value=Phi(leaf_cell,vertex);
                        if(abs(phi_value)<=fmm_initialization_iterative_tolerance_absolute){
                            T new_phi_value=(vertex-location).Magnitude();
                            if((new_phi_value-phi_new(index)) /max(fmm_initialization_iterative_tolerance_absolute,phi_new(index))<levelset.fmm_initialization_iterative_drift_fraction && 
                               new_phi_value>0)
                                phi_new(index)=new_phi_value;
                            break;}}}
                phi_new(index)*=LEVELSET_UTILITIES<T>::Sign(phi(index));}}
        else{
            for(int index=1;index<=Number_Of_Computational_Nodes();index++)if(Use_Computational_Node(index))if(done(index)){
                T value[T_GRID::dimension]; // the phi value to use in the given direction
                int number_of_axis=0; // the number of axis that we want to use later
                TV location=Computational_Node_Location(index);
                for(int axis=1;axis<=T_GRID::dimension;axis++){
                    int low=Neighbor(2*axis-1,index),high=Neighbor(2*axis,index);
                    bool use_low=(done(low)&&LEVELSET_UTILITIES<T>::Interface(phi(index),phi(low))),use_high=(done(high)&&LEVELSET_UTILITIES<T>::Interface(phi(index),phi(high)));
                    if(!use_low){if(use_high)value[number_of_axis]=LEVELSET_UTILITIES<T>::Theta(phi(index),phi(high))*(Computational_Node_Location(high)[axis]-location[axis]);else number_of_axis--;}
                    else if(!use_high) value[number_of_axis]=LEVELSET_UTILITIES<T>::Theta(phi(index),phi(low))*(location[axis]-Computational_Node_Location(low)[axis]);
                    else value[number_of_axis]=min(LEVELSET_UTILITIES<T>::Theta(phi(index),phi(high))*(Computational_Node_Location(high)[axis]-location[axis]),
                                                   LEVELSET_UTILITIES<T>::Theta(phi(index),phi(low))*(location[axis]-Computational_Node_Location(low)[axis]));
                    number_of_axis++;}
                if(number_of_axis==1) phi_new(index)=value[0];
                else if(number_of_axis==2){T d2=sqr(value[0])+sqr(value[1]);if(d2 > sqr_epsilon_dx_dy) phi_new(index)=value[0]*value[1]/sqrt(d2);else phi_new(index)=min(value[0],value[1]);}
                else{assert(T_GRID::dimension==3); // 2d should never get to this point
                    T value_xy=value[0]*value[1],value_xz=value[0]*value[2],value_yz=value[1]*value[2],d2=sqr(value_xy)+sqr(value_xz)+sqr(value_yz);
                    if(d2 > sqr_epsilon_dx_dy_dz) phi_new(index)=value_xy*value[2]/sqrt(d2);else phi_new(index)=min(value[0],value[1],value[2]);}
                if(levelset.refine_fmm_initialization_with_iterative_solver){
                    TV vertex=location;T phi_value=phi(index);
                    CELL* leaf_cell=grid.Leaf_Cell(vertex);
                    for(int iterations=1;iterations<=levelset.fmm_initialization_iterations;iterations++){
                        if(abs(phi_value)>10*grid.Minimum_Edge_Length())break; // stop if it looks like it's blowing up
                        vertex-=phi_value*Normal(leaf_cell,vertex);
                        leaf_cell=grid.Leaf_Cell_By_Neighbor_Path(leaf_cell,vertex);if(!leaf_cell)break; // we are outside the grid
                        phi_value=Phi(leaf_cell,vertex);
                        if(abs(phi_value)<=fmm_initialization_iterative_tolerance_absolute){
                            T new_phi_value=(vertex-location).Magnitude();
                            if((new_phi_value-phi_new(index))/max(fmm_initialization_iterative_tolerance_absolute,phi_new(index))<levelset.fmm_initialization_iterative_drift_fraction && 
                               new_phi_value>0)
                                phi_new(index)=new_phi_value;
                            break;}}}
                phi_new(index)*=LEVELSET_UTILITIES<T>::Sign(phi(index));}}
        for(int index=1;index<=Number_Of_Computational_Nodes();index++)if(Use_Computational_Node(index))if(done(index)) phi(index)=phi_new(index);}   // initialize done points

    // initialize close points
    for(int i=1;i<=Number_Of_Computational_Nodes();i++)if(Use_Computational_Node(i))if(close_k(i)){
        Update_Close_Point(phi,done,i);
        heap_length++;heap(heap_length)=i;
        FAST_MARCHING<T>::Up_Heap(phi,close_k,heap,heap_length);}
}
//#####################################################################
// Function Update_Close_Point
//##################################################################### 
// needs done=0 around the outside of the domain
template<class T_GRID> void FAST_MARCHING_METHOD_DYADIC<T_GRID>::
Update_Close_Point(ARRAY<T>& phi,const ARRAY<bool,VECTOR<int,1> >& done,const int index)
{
    T value[T_GRID::dimension]; // the phi value to use in the given direction
    T dx[T_GRID::dimension]; // the edge length in the given direction
    int number_of_axis=0; // the number of axis that we want to use later

    TV node_location_i=Computational_Node_Location(index);
    // check each principal axis
    for(int axis=1;axis<=T_GRID::dimension;axis++){
        int low=Neighbor(2*axis-1,index),high=Neighbor(2*axis,index);
        bool check_low=done(low),check_high=done(high);
        if(levelset.collision_body_list){
            if(check_low && !Neighbor_Visible(2*axis-1,index)) check_low=false;
            if(check_high && !Neighbor_Visible(2*axis,index)) check_high=false;}
        if(!check_low){
            if(check_high){value[number_of_axis]=phi(high);dx[number_of_axis]=Computational_Node_Location(high)[axis]-node_location_i[axis];}
            else number_of_axis--;}
        else if(!check_high){value[number_of_axis]=phi(low);dx[number_of_axis]=node_location_i[axis]-Computational_Node_Location(low)[axis];}
        else{
            if(abs(phi(low))<=abs(phi(high))){value[number_of_axis]=phi(low);dx[number_of_axis]=node_location_i[axis]-Computational_Node_Location(low)[axis];}
            else{value[number_of_axis]=phi(high);dx[number_of_axis]=Computational_Node_Location(high)[axis]-node_location_i[axis];}}
        number_of_axis++;}

    phi(index)=FAST_MARCHING<T>::template Solve_Close_Point<T_GRID::dimension>(phi(index),number_of_axis,value,dx);
}
//#####################################################################
// Function Add_To_Initial
//##################################################################### 
template<class T_GRID> void FAST_MARCHING_METHOD_DYADIC<T_GRID>::
Add_To_Initial(ARRAY<bool,VECTOR<int,1> >& done,ARRAY<int>& close_k,const int index)
{
    done(index)=true;close_k(index)=0; // add to done, remove from close 
    // add neighbors to close if not done 
    if(levelset.collision_body_list){
        for(int i=1;i<=T_GRID::number_of_faces_per_cell;i++){int neighbor_index=Neighbor(i,index);
            if(neighbor_index!=0&&!done(neighbor_index)&&Neighbor_Visible(i,index)) close_k(neighbor_index)=1;}}
    else
        for(int i=1;i<=T_GRID::number_of_faces_per_cell;i++){int neighbor_index=Neighbor(i,index);
            if(neighbor_index!=0&&!done(neighbor_index)) close_k(neighbor_index)=1;}
}
//#####################################################################
template class FAST_MARCHING_METHOD_DYADIC<OCTREE_GRID<float> >;
template class FAST_MARCHING_METHOD_DYADIC<QUADTREE_GRID<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FAST_MARCHING_METHOD_DYADIC<OCTREE_GRID<double> >;
template class FAST_MARCHING_METHOD_DYADIC<QUADTREE_GRID<double> >;
#endif
#endif

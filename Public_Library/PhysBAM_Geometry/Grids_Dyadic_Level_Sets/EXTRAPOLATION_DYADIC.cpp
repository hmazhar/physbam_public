//#####################################################################
// Copyright 2003-2005, Ron Fedkiw, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/EXTRAPOLATION_DYADIC.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T2> EXTRAPOLATION_DYADIC<T_GRID,T2>::
EXTRAPOLATION_DYADIC(T_GRID& grid_input)
    :grid(grid_input),currently_using_cells(false),currently_using_faces(false),
    cell_pointer_from_index(grid_input.Cell_Pointer_From_Index()),node_neighbors(grid_input.Node_Neighbors()),cell_neighbors(grid_input.Neighbors()),collision_aware_extrapolation(false),
    neighbors_visible(0),face_iterator_for_neighbors(grid_input,grid_input.Map_All_Faces()),
    boundary_scalar_default(*new T_BOUNDARY_SCALAR),boundary_extrapolate_default(*new T_BOUNDARY_EXTRAPOLATE)
{
    Set_Band_Width();
    Set_Small_Number();
    Set_Isobaric_Fix_Width();
    Set_Custom_Seed_Indices();
    Set_Custom_Seed_Done();
    tolerance=small_number*grid.Minimum_Edge_Length();
    boundary_scalar=&boundary_scalar_default;
    boundary_extrapolate=&boundary_extrapolate_default;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID,class T2> EXTRAPOLATION_DYADIC<T_GRID,T2>::
~EXTRAPOLATION_DYADIC()
{
    delete &boundary_scalar_default;delete &boundary_extrapolate_default;
}
//#####################################################################
// Function Extrapolate_Nodes
//#####################################################################
template<class T_GRID,class T2> void EXTRAPOLATION_DYADIC<T_GRID,T2>::
Extrapolate_Nodes(ARRAY<T> phi_nodes,ARRAY<T2>& u,const bool zero_outside_band,const T time)
{
    ARRAY<T> phi_nodes_ghost(grid.number_of_nodes,false);
    ARRAY<T2> u_ghost(grid.number_of_nodes,false); 
    boundary_scalar->Fill_Ghost_Cells_Node(grid,phi_nodes,phi_nodes_ghost,time);
    boundary_extrapolate->Fill_Ghost_Cells_Node(grid,u,u_ghost,time);
    if(zero_outside_band)for(NODE_ITERATOR iterator(grid,grid.Map_All_Nodes());iterator.Valid();iterator.Next())
        if(phi_nodes_ghost(iterator.Node_Index())>=band_width)u(iterator.Node_Index())=T2();
    if(isobaric_fix_width) phi_nodes_ghost+=isobaric_fix_width;
  
    Extrapolate(phi_nodes_ghost,u_ghost,time);
    ARRAY<T2>::Copy(u_ghost,u);
    boundary_extrapolate->Apply_Boundary_Condition(grid,u,time);
}
//#####################################################################
// Function Extrapolate_Faces
//#####################################################################
template<class T_GRID,class T2> void EXTRAPOLATION_DYADIC<T_GRID,T2>::
Extrapolate_Faces(ARRAY<T> phi_faces,ARRAY<T>& u,const bool zero_outside_band,const T time,const bool extrapolate_nodes)
{
    ARRAY<T> phi_faces_ghost(grid.number_of_faces,false);
    ARRAY<T> u_ghost(grid.number_of_faces,false); 
    boundary_scalar->Fill_Ghost_Cells_Face(grid,phi_faces,phi_faces_ghost,time);
    boundary_scalar->Fill_Ghost_Cells_Face(grid,u,u_ghost,time);
    if(zero_outside_band)for(FACE_ITERATOR iterator(grid,grid.Map_All_Faces());iterator.Valid();iterator.Next())
        if(phi_faces_ghost(iterator.Face_Index())>=band_width)u(iterator.Face_Index())=0;
    if(isobaric_fix_width) phi_faces_ghost+=isobaric_fix_width;

    currently_using_faces=true;
    Extrapolate(phi_faces_ghost,u_ghost,time);

    // interpolate phi values to the nodes
    if(extrapolate_nodes){
        ARRAY<TV> u_nodes(grid.number_of_nodes);
        ARRAY<TV> weight(grid.number_of_nodes);
        ARRAY<int> number_of_valid_cells_touching(grid.number_of_nodes);
        for(FACE_ITERATOR iterator(grid,grid.Map_All_Faces());iterator.Valid();iterator.Next()){
            T face_weight=1/iterator.Face_DX();
            T weight_times_u=face_weight*u_ghost(iterator.Face_Index());
            for(int i=0;i<T_GRID::number_of_nodes_per_face;i++){
                u_nodes(iterator.Face_Node(i))[iterator.Axis()+1]+=weight_times_u;
                weight(iterator.Face_Node(i))[iterator.Axis()+1]+=face_weight;}}
        for(CELL_ITERATOR iterator(grid,grid.number_of_ghost_cells);iterator.Valid();iterator.Next()){
            if(iterator.Cell_Pointer()->Depth_Of_This_Cell()==grid.maximum_depth)
                for(int i=0;i<T_GRID::number_of_nodes_per_cell;i++)number_of_valid_cells_touching(iterator.Cell_Pointer()->Node(i))++;}
        // seed the node part of the extrapolation
        ARRAY<int> seed_nodes;seed_nodes.Preallocate(number_of_valid_cells_touching.Count_Matches(T_GRID::number_of_cells_per_node));
        for(NODE_ITERATOR iterator(grid,grid.Map_All_Nodes());iterator.Valid();iterator.Next()){
            u_nodes(iterator.Node_Index())/=weight(iterator.Node_Index());
            if(number_of_valid_cells_touching(iterator.Node_Index())==T_GRID::number_of_cells_per_node)seed_nodes.Append(iterator.Node_Index());}
        // extrapolate the nodes
        currently_using_faces=false;
        Set_Custom_Seed_Indices(seed_indices);
        ARRAY<T> phi_nodes;LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID>::Interpolate_From_Faces_To_Nodes(grid,phi_faces,phi_nodes);
        Extrapolate(phi_nodes,u_nodes,time);
        // intepolate the nodal values back onto the faces
        for(FACE_ITERATOR iterator(grid,grid.Map_All_Faces());iterator.Valid();iterator.Next()){
            if(iterator.Other_Cell()->Depth_Of_This_Cell()<grid.maximum_depth){
                u_ghost(iterator.Face_Index())=0;
                for(int i=0;i<T_GRID::number_of_nodes_per_face;i++)u_ghost(iterator.Face_Index())+=u_nodes(iterator.Face_Node(i))[iterator.Axis()+1];
                u_ghost(iterator.Face_Index())*=T_GRID::one_over_number_of_nodes_per_face;}}}

    ARRAY<T>::Copy(u_ghost,u);
    boundary_scalar->Apply_Boundary_Condition(grid,u,time);
}
//#####################################################################
// Function Extrapolate_Cells
//#####################################################################
template<class T_GRID,class T2> void EXTRAPOLATION_DYADIC<T_GRID,T2>::
Extrapolate_Cells(ARRAY<T> phi,ARRAY<T2>& u,const bool zero_outside_band,const T time,const bool extrapolate_nodes)
{
    ARRAY<T> phi_ghost(grid.number_of_cells,false);
    ARRAY<T2> u_ghost(grid.number_of_cells,false); 
    boundary_scalar->Fill_Ghost_Cells_Cell(grid,phi,phi_ghost,time);
    boundary_extrapolate->Fill_Ghost_Cells_Cell(grid,u,u_ghost,time);
    if(zero_outside_band)for(CELL_ITERATOR iterator(grid,grid.number_of_ghost_cells);iterator.Valid();iterator.Next())
        if(phi_ghost(iterator.Cell_Index())>=band_width)u(iterator.Cell_Index())=T2();
    if(isobaric_fix_width) phi_ghost+=isobaric_fix_width;
  
    currently_using_cells=true;
    Extrapolate(phi_ghost,u_ghost,time);

    // interpolate phi values to the nodes
    if(extrapolate_nodes){
        ARRAY<T2> u_nodes(grid.number_of_nodes);
        ARRAY<T> weight(grid.number_of_nodes);
        ARRAY<int> number_of_valid_cells_touching(grid.number_of_nodes);
        for(CELL_ITERATOR iterator(grid,grid.number_of_ghost_cells);iterator.Valid();iterator.Next()){
            T cell_weight=1/iterator.Cell_Pointer()->DX().x;
            T2 weight_times_u=cell_weight*u_ghost(iterator.Cell_Index());
            for(int i=0;i<T_GRID::number_of_nodes_per_cell;i++){
                u_nodes(iterator.Cell_Pointer()->Node(i))+=weight_times_u;
                weight(iterator.Cell_Pointer()->Node(i))+=cell_weight;}
            if(iterator.Cell_Pointer()->Depth_Of_This_Cell()==grid.maximum_depth)
                for(int i=0;i<T_GRID::number_of_nodes_per_cell;i++)number_of_valid_cells_touching(iterator.Cell_Pointer()->Node(i))++;}
        // seed the node part of the extrapolation
        ARRAY<int> seed_nodes;seed_nodes.Preallocate(number_of_valid_cells_touching.Count_Matches(T_GRID::number_of_cells_per_node));
        for(NODE_ITERATOR iterator(grid,grid.Map_All_Nodes());iterator.Valid();iterator.Next()){
            u_nodes(iterator.Node_Index())/=weight(iterator.Node_Index());
            if(number_of_valid_cells_touching(iterator.Node_Index())==T_GRID::number_of_cells_per_node)seed_nodes.Append(iterator.Node_Index());}
        // extrapolate the nodes
        currently_using_cells=false;
        Set_Custom_Seed_Indices(seed_indices);
        ARRAY<T> phi_nodes;LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID>::Interpolate_From_Cells_To_Nodes(grid,phi,phi_nodes);
        Extrapolate(phi_nodes,u_nodes,time);
        // intepolate the nodal values back onto the cells
        for(CELL_ITERATOR iterator(grid,grid.number_of_ghost_cells);iterator.Valid();iterator.Next()){
            if(iterator.Cell_Pointer()->Depth_Of_This_Cell()<grid.maximum_depth){
                u_ghost(iterator.Cell_Index())=T2();
                for(int i=0;i<T_GRID::number_of_nodes_per_cell;i++)u_ghost(iterator.Cell_Index())+=u_nodes(iterator.Cell_Pointer()->Node(i));
                u_ghost(iterator.Cell_Index())*=T_GRID::one_over_number_of_nodes_per_cell;}}}

    ARRAY<T2>::Copy(u_ghost,u);
    boundary_extrapolate->Apply_Boundary_Condition(grid,u,time);
}
//#####################################################################
// Function Extrapolate
//#####################################################################
template<class T_GRID,class T2> template<class T3> void EXTRAPOLATION_DYADIC<T_GRID,T2>::
Extrapolate(ARRAY<T>& phi_ghost,ARRAY<T3>& u_ghost,const T time)
{
    int heap_length=0;
    ARRAY<bool,VECTOR<int,1> > done(0,Number_Of_Computational_Nodes()); // starts at 0 so that the 0th node is never done - Update_Close_Point() relies on this
    ARRAY<bool> close(Number_Of_Computational_Nodes());
    ARRAY<int> heap(Number_Of_Computational_Nodes());
    Initialize(phi_ghost,done,close,heap,heap_length);

    while(heap_length && phi_ghost(heap(1)) <= band_width+isobaric_fix_width){
        int index=Remove_Root_From_Heap(phi_ghost,heap,heap_length,close);
        done(index)=true;
        Update_Close_Point(u_ghost,phi_ghost,done,index);

        if(collision_aware_extrapolation){
            for(int i=1;i<=T_GRID::number_of_faces_per_cell;i++){int neighbor_index=Neighbor(i,index);
                if(neighbor_index!=0 && !done(neighbor_index) && !close(neighbor_index) && Neighbor_Visible(i,index)) 
                    Add_To_Heap(phi_ghost,heap,heap_length,close,neighbor_index);}}
        else
            for(int i=1;i<=T_GRID::number_of_faces_per_cell;i++){int neighbor_index=Neighbor(i,index);
                if(neighbor_index!=0 && !done(neighbor_index) && !close(neighbor_index)) 
                    Add_To_Heap(phi_ghost,heap,heap_length,close,neighbor_index);}}
}
//#####################################################################
// Function Initialize
//#####################################################################
// pass heap_length by reference. Phi is cell centered, face, or nodal depending on which type of extrapolation is being done
template<class T_GRID,class T2> void EXTRAPOLATION_DYADIC<T_GRID,T2>::
Initialize(const ARRAY<T>& phi,ARRAY<bool,VECTOR<int,1> >& done,ARRAY<bool>& close,ARRAY<int>& heap,int& heap_length)
{
    assert(!seed_indices||!seed_done);
    if(seed_indices){for(int i=1;i<=seed_indices->m;i++) done((*seed_indices)(i))=true;}
    else if(seed_done){for(int i=1;i<=seed_done->m;i++) if((*seed_done)(i)) done(i)=true;}
    else if(currently_using_cells){for(CELL_ITERATOR iterator(grid,grid.number_of_ghost_cells);iterator.Valid();iterator.Next())
        if(iterator.Cell_Pointer()->Depth_Of_This_Cell()==grid.maximum_depth&&phi(iterator.Cell_Index())<=0)done(iterator.Cell_Index())=true;}
    else if(currently_using_faces){for(FACE_ITERATOR iterator(grid,grid.Map_All_Faces());iterator.Valid();iterator.Next())
        if(iterator.Deepest_Cell()->Depth_Of_This_Cell()==grid.maximum_depth&&phi(iterator.Face_Index())<=0)done(iterator.Face_Index())=true;}
    else{for(NODE_ITERATOR iterator(grid,grid.Map_All_Nodes());iterator.Valid();iterator.Next()){
        if(phi(iterator.Node_Index())<=0)done(iterator.Node_Index())=true;}}

    // find neighbors of done nodes which have positive phi
    if(collision_aware_extrapolation){
        if(currently_using_cells){
            for(CELL_ITERATOR iterator(grid,grid.number_of_ghost_cells);iterator.Valid();iterator.Next())if(done(iterator.Cell_Index())){
                for(int i=1;i<=T_GRID::number_of_neighbors_per_cell;i++){
                    int neighbor_index=cell_neighbors(iterator.Cell_Index())(i)->Cell();
                    if(neighbor_index!=0 && cell_pointer_from_index(neighbor_index)->Depth_Of_This_Cell()==grid.maximum_depth) // only use cells at the finest level
                        if(!done(neighbor_index) && !close(neighbor_index) && phi(neighbor_index)>0 && Neighbor_Visible(i,iterator.Cell_Index()))
                            Add_To_Heap(phi,heap,heap_length,close,neighbor_index);}}}
        else if(currently_using_faces){
            for(FACE_ITERATOR iterator(grid,grid.Map_All_Faces());iterator.Valid();iterator.Next())if(done(iterator.Face_Index())){
                for(int i=1;i<=T_GRID::number_of_neighbors_per_face;i++){
                    int neighbor_index=iterator.Minimal_Neighbor(i-1);
                    if(neighbor_index!=0 && !done(neighbor_index) && !close(neighbor_index) && phi(neighbor_index)>0 && Neighbor_Visible(i,iterator.Face_Index()))
                        Add_To_Heap(phi,heap,heap_length,close,neighbor_index);}}}
        else{
            for(NODE_ITERATOR iterator(grid,grid.Map_All_Nodes());iterator.Valid();iterator.Next())if(done(iterator.Node_Index())){
                for(int i=1;i<=T_GRID::number_of_neighbors_per_node;i++){
                    int neighbor_index=node_neighbors(iterator.Node_Index())(i);
                    if(neighbor_index!=0 && !done(neighbor_index) && !close(neighbor_index) && phi(neighbor_index)>0 && Neighbor_Visible(i,iterator.Node_Index())) 
                        Add_To_Heap(phi,heap,heap_length,close,neighbor_index);}}}}
    else{
        if(currently_using_cells){
            for(CELL_ITERATOR iterator(grid,grid.number_of_ghost_cells);iterator.Valid();iterator.Next())if(done(iterator.Cell_Index())){
                for(int i=1;i<=T_GRID::number_of_neighbors_per_cell;i++){
                    int neighbor_index=cell_neighbors(iterator.Cell_Index())(i)->Cell();
                    if(neighbor_index!=0 && cell_pointer_from_index(neighbor_index)->Depth_Of_This_Cell()==grid.maximum_depth) // only use cells at the finest level
                        if(!done(neighbor_index) && !close(neighbor_index) && phi(neighbor_index)>0) Add_To_Heap(phi,heap,heap_length,close,neighbor_index);}}}
        else if(currently_using_faces){
            for(FACE_ITERATOR iterator(grid,grid.Map_All_Faces());iterator.Valid();iterator.Next())if(done(iterator.Face_Index())){
                for(int i=1;i<=T_GRID::number_of_neighbors_per_face;i++){
                    int neighbor_index=iterator.Minimal_Neighbor(i-1);
                    if(neighbor_index!=0 && !done(neighbor_index) && !close(neighbor_index) && phi(neighbor_index)>0)
                        Add_To_Heap(phi,heap,heap_length,close,neighbor_index);}}}
        else{
            for(NODE_ITERATOR iterator(grid,grid.Map_All_Nodes());iterator.Valid();iterator.Next())if(done(iterator.Node_Index())){
                for(int i=1;i<=T_GRID::number_of_neighbors_per_node;i++){
                    int neighbor_index=node_neighbors(iterator.Node_Index())(i);
                    if(neighbor_index!=0 && !done(neighbor_index) && !close(neighbor_index) && phi(neighbor_index)>0) Add_To_Heap(phi,heap,heap_length,close,neighbor_index);}}}}
}
//#####################################################################
// Function Update_Close_Point
//##################################################################### 
// needs done=0 around the outside of the domain
template<class T_GRID,class T2> template<class T3> void EXTRAPOLATION_DYADIC<T_GRID,T2>::
Update_Close_Point(ARRAY<T3>& u,const ARRAY<T>& phi,const ARRAY<bool,VECTOR<int,1> >& done,const int i)
{      
    T3 value[T_GRID::dimension]; // the velocity value to use in the given direction
    T dx[T_GRID::dimension]; // the edge length in the given direction
    T phix_dx[T_GRID::dimension]; // the difference in phi value for the direction
    int number_of_axis=0; // the number of axis that we want to use later

    TV node_location_i=Computational_Node_Location(i);
    // check each principal axis
    for(int axis=1;axis<=T_GRID::dimension;axis++){
        int low=Neighbor(2*axis-1,i),high=Neighbor(2*axis,i);
        bool check_low=done(low),check_high=done(high);
        if(collision_aware_extrapolation){
            if(check_low && !Neighbor_Visible(2*axis-1,i)) check_low=false;
            if(check_high && !Neighbor_Visible(2*axis,i)) check_high=false;}
        if(!check_low){
            if(check_high){value[number_of_axis]=u(high);phix_dx[number_of_axis]=phi(i)-phi(high);dx[number_of_axis]=Computational_Node_Location(high)[axis]-node_location_i[axis];}
            else number_of_axis--;}
        else if(!check_high){value[number_of_axis]=u(low);phix_dx[number_of_axis]=phi(i)-phi(low);dx[number_of_axis]=node_location_i[axis]-Computational_Node_Location(low)[axis];}
        else{
            if(phi(low)<=phi(high)){value[number_of_axis]=u(low);phix_dx[number_of_axis]=phi(i)-phi(low);dx[number_of_axis]=node_location_i[axis]-Computational_Node_Location(low)[axis];}
            else{value[number_of_axis]=u(high);phix_dx[number_of_axis]=phi(i)-phi(high);dx[number_of_axis]=Computational_Node_Location(high)[axis]-node_location_i[axis];}}
        number_of_axis++;}
       
    assert(number_of_axis);
    if(number_of_axis==1) u(i)=value[0];
    else if(T_GRID::dimension==2 || number_of_axis==2){
        T a=phix_dx[0],b=phix_dx[1],denominator=a+b,fraction=(T).5; // a should be multiplied by sqr(dy/dx) (or whatever axis there are in the calculation), but we only support square grids, so this is ok.
        if(denominator > small_number) fraction=max((T)0,min((T)1,a/denominator));
        u(i)=fraction*value[0]+(1-fraction)*value[1];}
    else{assert(T_GRID::dimension==3); // should only get here in 3D
        T a=phix_dx[0],b=phix_dx[1],c=phix_dx[2],denominator=a+b+c; // a,b should be multiplied by sqr(dy/dx), but we only support square grids, so this is ok.
        T fraction_1=(T)one_third,fraction_2=(T)one_third;
        if(denominator > small_number){fraction_1=max((T)0,min((T)1,a/denominator));fraction_2=max((T)0,min((T)1-fraction_1,b/denominator));}
        u(i)=fraction_1*value[0]+fraction_2*value[1]+((T)1-fraction_1-fraction_2)*value[2];}
}
//#####################################################################
template class EXTRAPOLATION_DYADIC<OCTREE_GRID<float>,float>;
template class EXTRAPOLATION_DYADIC<QUADTREE_GRID<float>,float>;
template class EXTRAPOLATION_DYADIC<OCTREE_GRID<float>,VECTOR<float,3> >;
template class EXTRAPOLATION_DYADIC<QUADTREE_GRID<float>,VECTOR<float,2> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EXTRAPOLATION_DYADIC<OCTREE_GRID<double>,double>;
template class EXTRAPOLATION_DYADIC<QUADTREE_GRID<double>,double>;
template class EXTRAPOLATION_DYADIC<OCTREE_GRID<double>,VECTOR<double,3> >;
template class EXTRAPOLATION_DYADIC<QUADTREE_GRID<double>,VECTOR<double,2> >;
#endif
#endif

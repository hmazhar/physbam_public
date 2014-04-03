//#####################################################################
// Copyright 2005, Geoffrey Irving
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Collisions/GRID_BASED_COLLISION_GEOMETRY_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/FAST_MARCHING_METHOD_RLE.h>
#include <PhysBAM_Geometry/Level_Sets/FAST_MARCHING.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> FAST_MARCHING_METHOD_RLE<T_GRID>::
FAST_MARCHING_METHOD_RLE(const LEVELSET_RLE<T_GRID>& levelset)
    :levelset(levelset),grid(levelset.grid),neighbors(grid.Short_Cell_Neighbors()),
    neighbors_visible(levelset.collision_body_list?&levelset.collision_body_list->cell_neighbors_visible:0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> FAST_MARCHING_METHOD_RLE<T_GRID>::
~FAST_MARCHING_METHOD_RLE()
{}
//#####################################################################
// Function Fast_Marching
//#####################################################################
template<class T_GRID> void FAST_MARCHING_METHOD_RLE<T_GRID>::
Fast_Marching(ARRAY<T>& phi,const T stopping_distance,const bool outside_only) const
{
    if(&phi!=&levelset.phi) PHYSBAM_FATAL_ERROR("phi should be levelset.phi so that we can use levelset.Phi and levelset.Normal");
    ARRAY<bool,VECTOR<int,1> > done(0,grid.number_of_cells); // starts at 0 so that the 0th cell is never done - Update_Close_Point() relies on this
    Initialize_Interface(phi,done,outside_only);

    int heap_length=0;
    ARRAY<int> close_k(grid.number_of_cells);
    ARRAY<int> heap(grid.number_of_cells,false);
    T far_negative=-min(grid.negative_bandwidth,stopping_distance),far_positive=min(grid.positive_bandwidth,stopping_distance);
    for(int c=1;c<=grid.number_of_cells;c++){
        if(outside_only && phi(c)<0) continue;
        if(done(c)){
            for(int n=1;n<=T_GRID::number_of_neighbors_per_cell;n++){int neighbor=Neighbor_If_Visible(n,c);
                if(neighbor&&!done(neighbor)&&!close_k(neighbor)) Update_Or_Add_Neighbor(phi,done,close_k,heap,heap_length,neighbor);}}
        else if(!close_k(c)) phi(c)=phi(c)<=0?far_negative:far_positive;}

    while(heap_length){
        int cell=heap(1); // smallest point is on top of heap
        if(abs(phi(cell)) > stopping_distance) break; // exit early
        done(cell)=true;close_k(cell)=0; // add to done, remove from close
        FAST_MARCHING<T>::Down_Heap(phi,close_k,heap,heap_length);heap_length--; // remove point from heap

        for(int n=1;n<=T_GRID::number_of_neighbors_per_cell;n++){int neighbor=Neighbor_If_Visible(n,cell);
            if(neighbor&&!done(neighbor)) Update_Or_Add_Neighbor(phi,done,close_k,heap,heap_length,neighbor);}}
}
//#####################################################################
// Function Slow_Marching
//#####################################################################
template<class T_GRID> void FAST_MARCHING_METHOD_RLE<T_GRID>::
Slow_Marching(ARRAY<T>& phi,const T stopping_distance,const T tolerance) const
{
    if(&phi!=&levelset.phi) PHYSBAM_FATAL_ERROR("phi should be levelset.phi so that we can use levelset.Phi and levelset.Normal");
    ARRAY<bool,VECTOR<int,1> > done(0,grid.number_of_cells); // starts at 0 so that the 0th cell is never done - Update_Close_Point() relies on this
    Initialize_Interface(phi,done);

    T absolute_tolerance=tolerance*grid.Minimum_Edge_Length();
    ARRAY<unsigned char> needs_update_in_iteration(grid.number_of_cells);
    ARRAY<bool,VECTOR<int,1> > touched(done);

    // mark cells near interface
    T far_negative=-min(grid.negative_bandwidth,stopping_distance),far_positive=min(grid.positive_bandwidth,stopping_distance);
    for(int cell=1;cell<=grid.number_of_cells;cell++){
        if(!done(cell)) phi(cell)=phi(cell)<=0?far_negative:far_positive;
        else for(int n=1;n<=T_GRID::number_of_neighbors_per_cell;n++){int neighbor=Neighbor_If_Visible(n,cell);
            if(neighbor&&!done(neighbor)) needs_update_in_iteration(neighbor)=1;}}

    // loop until nothing changes
    for(unsigned char iteration=1;;iteration++){
        bool cell_marked=false;
        for(int cell=1;cell<=grid.number_of_cells;cell++)if(needs_update_in_iteration(cell)==iteration){
            T old_phi=phi(cell);touched(cell)=true;
            Update_Close_Point(phi,touched,cell);
            if(abs(old_phi-phi(cell))>absolute_tolerance && abs(phi(cell))<stopping_distance)
                for(int n=1;n<=T_GRID::number_of_neighbors_per_cell;n++){int neighbor=Neighbor_If_Visible(n,cell);
                    if(neighbor && !done(neighbor) && abs(phi(cell))<abs(phi(neighbor))){needs_update_in_iteration(neighbor)=iteration+1;cell_marked=true;}}}
        if(!cell_marked) break;}
}
//#####################################################################
// Function Update_Or_Add_Neighbor
//#####################################################################
template<class T_GRID> inline void FAST_MARCHING_METHOD_RLE<T_GRID>::
Update_Or_Add_Neighbor(ARRAY<T>& phi,ARRAY<bool,VECTOR<int,1> >& done,ARRAY<int>& close_k,ARRAY<int>& heap,int& heap_length,const int neighbor) const
{
    if(close_k(neighbor)){
        Update_Close_Point(phi,done,neighbor);
        FAST_MARCHING<T>::Up_Heap(phi,close_k,heap,close_k(neighbor));}
    else{close_k(neighbor)=1; // add to close
        Update_Close_Point(phi,done,neighbor);
        heap_length++;heap(heap_length)=neighbor;
        FAST_MARCHING<T>::Up_Heap(phi,close_k,heap,heap_length);}
}
//#####################################################################
// Function Initialize_Interface
//#####################################################################
template<class T_GRID> void FAST_MARCHING_METHOD_RLE<T_GRID>::
Initialize_Interface(ARRAY<T>& phi,ARRAY<bool,VECTOR<int,1> >& done,const bool outside_only) const
{
    T sqr_epsilon_cell_size=sqr(levelset.small_number*grid.uniform_grid.Cell_Size()),sqr_epsilon_face_size[3];
    if(T_GRID::dimension==2){sqr_epsilon_face_size[2]=sqr_epsilon_cell_size;}
    else if(T_GRID::dimension==3) for(int axis=1;axis<=3;axis++) sqr_epsilon_face_size[axis-1]=sqr_epsilon_cell_size/sqr(grid.uniform_grid.dX[axis]);
    T fmm_initialization_iterative_tolerance_absolute=levelset.fmm_initialization_iterative_tolerance*grid.uniform_grid.Minimum_Edge_Length();
    // find interface cells
    T_GRID::template Face_Loop<Find_Interface_Cells>(*this,phi,done);
    // set interface cell phi values
    ARRAY<T> phi_new(grid.number_of_cells,false);
    phi_new.Fill(2*grid.uniform_grid.dX.Max()); // ok positive, since minmag is used below
    for(CELL_ITERATOR cell(grid,grid.number_of_ghost_cells);cell;cell++){int c=cell.Cell();if(done(c)){
        if(outside_only && phi(c)<0) continue;
        T value[T_GRID::dimension]; // the phi value to use in the given direction
        int number_of_axis=0; // the number of axis that we want to use later
        int missing_axis=3; // used in number_of_axis==2 case only, so it gives you which axis is missing (==3 for 2d)
        bool really_clamp_phi_with_collision_bodies=levelset.clamp_phi_with_collision_bodies&&phi(c)<=0;
        T abs_phi=abs(phi(c));
        TV X=cell.X();
        for(int axis=1;axis<=T_GRID::dimension;axis++){T dx=grid.uniform_grid.dX[axis];
            int low=neighbors(c)(2*axis-1),high=neighbors(c)(2*axis);
            bool use_low=false,use_high=false;T value_low=0,value_high=0;
            if(low){
                if(levelset.collision_body_list && !Neighbor_Visible(axis,low,c)){
                    RAY<TV> ray(X,-TV::Axis_Vector(axis),true);ray.t_max=dx;ray.semi_infinite=false;
                    COLLISION_GEOMETRY_ID body_id;levelset.collision_body_list->Closest_Non_Intersecting_Point_Of_Any_Body(ray,body_id);
                    T phi_object=max(ray.t_max,levelset.small_number);
                    use_low=true;if(really_clamp_phi_with_collision_bodies) value_low=min(abs_phi,phi_object);else value_low=phi_object;}
                else if(done(low)&&LEVELSET_UTILITIES<T>::Interface(phi(c),phi(low))){
                    use_low=true;value_low=LEVELSET_UTILITIES<T>::Theta(phi(c),phi(low))*dx;}}
            if(high){
                if(levelset.collision_body_list && !Neighbor_Visible(axis,c,high)){
                    RAY<TV> ray(X,TV::Axis_Vector(axis),true);ray.t_max=dx;ray.semi_infinite=false;
                    COLLISION_GEOMETRY_ID body_id;levelset.collision_body_list->Closest_Non_Intersecting_Point_Of_Any_Body(ray,body_id);
                    T phi_object=max(ray.t_max,levelset.small_number);
                    use_high=true;if(really_clamp_phi_with_collision_bodies) value_high=min(abs_phi,phi_object);else value_high=phi_object;}
                else if(done(high)&&LEVELSET_UTILITIES<T>::Interface(phi(c),phi(high))){
                    use_high=true;value_high=LEVELSET_UTILITIES<T>::Theta(phi(c),phi(high))*dx;}}
            if(!use_low){if(use_high) value[number_of_axis]=value_high;else{missing_axis=axis;number_of_axis--;}}
            else if(!use_high) value[number_of_axis]=value_low;
            else value[number_of_axis]=min(value_low,value_high);
            number_of_axis++;}
        assert(number_of_axis);
        if(number_of_axis==1) phi_new(c)=value[0];
        else if(number_of_axis==2){T d2=sqr(value[0])+sqr(value[1]);
            if(d2 > sqr_epsilon_face_size[missing_axis-1]) phi_new(c)=value[0]*value[1]/sqrt(d2);else phi_new(c)=min(value[0],value[1]);}
        else{assert(T_GRID::dimension==3); // 2d should never get to this point
            T value_xy=value[0]*value[1],value_xz=value[0]*value[2],value_yz=value[1]*value[2],d2=sqr(value_xy)+sqr(value_xz)+sqr(value_yz);
            if(d2 > sqr_epsilon_cell_size) phi_new(c)=value_xy*value[2]/sqrt(d2);else phi_new(c)=min(value[0],value[1],value[2]);}
        if(levelset.refine_fmm_initialization_with_iterative_solver){
            assert(levelset.normals);
            TV vertex_minus_X=-phi(c)*(*levelset.normals)(c);
            for(int iterations=1;iterations<=levelset.fmm_initialization_iterations;iterations++){
                TV vertex=X+vertex_minus_X;
                T_BLOCK block(grid,vertex);if(!block) break;
                T phi_value=levelset.Phi(block,vertex);
                if(abs(phi_value)<=fmm_initialization_iterative_tolerance_absolute){
                    T new_phi_value=(vertex_minus_X).Magnitude();
                    if((new_phi_value-phi_new(c))/max(fmm_initialization_iterative_tolerance_absolute,phi_new(c))<levelset.fmm_initialization_iterative_drift_fraction && new_phi_value>0)
                        phi_new(c)=new_phi_value;
                    break;}
                if(abs(phi_value)>10*grid.uniform_grid.Minimum_Edge_Length()) break; // stop if it looks like it's blowing up
                TV normal=levelset.Normal(block,vertex);
                vertex_minus_X=(TV::Dot_Product(vertex_minus_X,normal)-phi_value)*normal;}}
        phi_new(c)*=LEVELSET_UTILITIES<T>::Sign(phi(c));}}
    // initialize done phi
    for(int c=1;c<=grid.number_of_cells;c++)if(done(c)) phi(c)=phi_new(c);
}
//#####################################################################
// Function Find_Interface_Cells
//#####################################################################
template<class T_GRID> template<class T_FACE> void FAST_MARCHING_METHOD_RLE<T_GRID>::
Find_Interface_Cells::Apply(const FAST_MARCHING_METHOD_RLE<T_GRID>& fast_marching,ARRAY<T>& phi,ARRAY<bool,VECTOR<int,1> >& done)
{
    const T_GRID& grid=fast_marching.grid;
    for(T_FACE face(grid,grid.number_of_ghost_cells,false);face;face++)if(face.Both_Cells_Short()){int c1=face.cell1.Cell(),c2=face.cell2.Cell();
        if(fast_marching.levelset.collision_body_list && !fast_marching.Neighbor_Visible(face.Axis(),c1,c2)){
            if(fast_marching.levelset.clamp_phi_with_collision_bodies && fast_marching.levelset.collision_body_list->Refine_Due_To_Objects(c1,c2)){
                if(phi(c1)<=0) done(c1)=true;
                if(phi(c2)<=0) done(c2)=true;}}
        else if(LEVELSET_UTILITIES<T>::Interface(phi(c1),phi(c2))) done(c1)=done(c2)=true;}
}
//#####################################################################
// Function Update_Close_Point
//##################################################################### 
template<class T_GRID> void FAST_MARCHING_METHOD_RLE<T_GRID>::
Update_Close_Point(ARRAY<T>& phi,const ARRAY<bool,VECTOR<int,1> >& done,const int cell) const
{
    T value[T_GRID::dimension]; // the phi value to use in the given direction
    T dx[T_GRID::dimension]; // the edge length in the given direction
    int number_of_axis=0; // the number of axis that we want to use later

    // check each principal axis
    for(int axis=1;axis<=T_GRID::dimension;axis++){
        int low=neighbors(cell)(2*axis-1),high=neighbors(cell)(2*axis);
        bool check_low=done(low),check_high=done(high);
        if(levelset.collision_body_list){
            if(check_low && !Neighbor_Visible(axis,low,cell)) check_low=false;
            if(check_high && !Neighbor_Visible(axis,cell,high)) check_high=false;}
        dx[number_of_axis]=grid.uniform_grid.dX[axis];
        if(!check_low){
            if(check_high) value[number_of_axis]=phi(high);
            else number_of_axis--;}
        else if(!check_high) value[number_of_axis]=phi(low);
        else value[number_of_axis]=minmag(phi(low),phi(high));
        number_of_axis++;}

    phi(cell)=FAST_MARCHING<T>::template Solve_Close_Point<T_GRID::dimension>(phi(cell),number_of_axis,value,dx);
}
//#####################################################################
template class FAST_MARCHING_METHOD_RLE<RLE_GRID_2D<float> >;
template class FAST_MARCHING_METHOD_RLE<RLE_GRID_3D<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FAST_MARCHING_METHOD_RLE<RLE_GRID_2D<double> >;
template class FAST_MARCHING_METHOD_RLE<RLE_GRID_3D<double> >;
#endif
#endif

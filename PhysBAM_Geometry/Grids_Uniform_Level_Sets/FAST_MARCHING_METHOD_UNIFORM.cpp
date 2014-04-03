//#####################################################################
// Copyright 2005-2007, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#endif
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/FAST_MARCHING.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> FAST_MARCHING_METHOD_UNIFORM<T_GRID>::
FAST_MARCHING_METHOD_UNIFORM(const T_LEVELSET& levelset,const int ghost_cells_input,THREAD_QUEUE* thread_queue_input)
    :levelset(levelset),ghost_cells(ghost_cells_input),thread_queue(thread_queue_input)
{
    cell_grid=levelset.grid.Is_MAC_Grid()?levelset.grid:levelset.grid.Get_MAC_Grid_At_Regular_Positions();
    RANGE<TV_INT> domain_indices=cell_grid.Domain_Indices().Thickened(ghost_cells);dimension_start=domain_indices.Minimum_Corner();dimension_end=domain_indices.Maximum_Corner();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> FAST_MARCHING_METHOD_UNIFORM<T_GRID>::
~FAST_MARCHING_METHOD_UNIFORM()
{}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class T_GRID> void FAST_MARCHING_METHOD_UNIFORM<T_GRID>::
Fast_Marching_Method(T_ARRAYS_SCALAR& phi_ghost,const T stopping_distance,const ARRAY<TV_INT>* seed_indices,const bool add_seed_indices_for_ghost_cells)
{
    int heap_length=0;
    T_ARRAYS_BOOL done(cell_grid.Domain_Indices(ghost_cells+1));
    T_ARRAYS_INT close_k(cell_grid.Domain_Indices(ghost_cells+1)); // extra cell so that it is the same size as the done array for optimizations
    ARRAY<TV_INT> heap(cell_grid.Domain_Indices().Thickened(ghost_cells).Size(),false); // a generic version of heap((m+6)*(n+6)*(mn+6),false);
    Initialize_Interface(phi_ghost,done,close_k,heap,heap_length,seed_indices,add_seed_indices_for_ghost_cells);

    while(heap_length != 0){
        TV_INT index=heap(1); // smallest point is on top of heap
        if(stopping_distance && abs(phi_ghost(index)) > stopping_distance){ // exit early
            for(CELL_ITERATOR iterator(cell_grid);iterator.Valid();iterator.Next()) if(!done(iterator.Cell_Index())){
                phi_ghost(iterator.Cell_Index())=LEVELSET_UTILITIES<T>::Sign(phi_ghost(iterator.Cell_Index()))*stopping_distance;}
            return;}
        done(index)=true;close_k(index)=0; // add to done, remove from close
        FAST_MARCHING<T>::Down_Heap(phi_ghost,close_k,heap,heap_length);heap_length--; // remove point from heap

        if(levelset.collision_body_list)
            for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                if(index[axis] != dimension_start[axis] && !done(index-axis_vector) && Neighbor_Visible(axis,index-axis_vector))
                    Update_Or_Add_Neighbor(phi_ghost,done,close_k,heap,heap_length,index-axis_vector);
                if(index[axis] != dimension_end[axis] && !done(index+axis_vector) && Neighbor_Visible(axis,index))
                    Update_Or_Add_Neighbor(phi_ghost,done,close_k,heap,heap_length,index+axis_vector);}
        else
            for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                if(index[axis] != dimension_start[axis] && !done(index-axis_vector))
                    Update_Or_Add_Neighbor(phi_ghost,done,close_k,heap,heap_length,index-axis_vector);
                if(index[axis] != dimension_end[axis] && !done(index+axis_vector))
                    Update_Or_Add_Neighbor(phi_ghost,done,close_k,heap,heap_length,index+axis_vector);}}
    //DOMAIN_ITERATOR_THREADED_ALPHA<FAST_MARCHING_METHOD_UNIFORM<T_GRID>,TV>(cell_grid.Domain_Indices(ghost_cells),thread_queue,1,ghost_cells,2,1).template Run<T_ARRAYS_SCALAR&,T,const ARRAY<TV_INT>*,bool>(*this,&FAST_MARCHING_METHOD_UNIFORM<T_GRID>::Fast_Marching_Method_Threaded,phi_ghost,stopping_distance,seed_indices,add_seed_indices_for_ghost_cells);
}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class T_GRID> void FAST_MARCHING_METHOD_UNIFORM<T_GRID>::
Fast_Marching_Method_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_SCALAR& phi_ghost,const T stopping_distance,const ARRAY<TV_INT>* seed_indices,const bool add_seed_indices_for_ghost_cells)
{
    int heap_length=0;
    T_ARRAYS_SCALAR phi_new(domain);
    for(CELL_ITERATOR iterator(cell_grid,domain);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();phi_new(cell)=phi_ghost(cell);}
    T_ARRAYS_BOOL done(domain.Thickened(1));
    T_ARRAYS_INT close_k(domain.Thickened(1)); // extra cell so that it is the same size as the done array for optimizations
    ARRAY<TV_INT> heap(domain.Size(),false); // a generic version of heap((m+6)*(n+6)*(mn+6),false);
    Initialize_Interface(domain,phi_new,done,close_k,heap,heap_length,seed_indices,add_seed_indices_for_ghost_cells);

    while(heap_length != 0){
        TV_INT index=heap(1); // smallest point is on top of heap
        if(stopping_distance && abs(phi_new(index)) > stopping_distance){ // exit early
            for(CELL_ITERATOR iterator(cell_grid,domain);iterator.Valid();iterator.Next()) if(!done(iterator.Cell_Index())){
                phi_new(iterator.Cell_Index())=LEVELSET_UTILITIES<T>::Sign(phi_new(iterator.Cell_Index()))*stopping_distance;}
            return;}
        done(index)=true;close_k(index)=0; // add to done, remove from close
        FAST_MARCHING<T>::Down_Heap(phi_new,close_k,heap,heap_length);heap_length--; // remove point from heap

        if(levelset.collision_body_list)
            for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                if(index[axis] != domain.min_corner[axis] && !done(index-axis_vector) && Neighbor_Visible(axis,index-axis_vector))
                    Update_Or_Add_Neighbor(phi_new,done,close_k,heap,heap_length,index-axis_vector);
                if(index[axis] != domain.max_corner[axis] && !done(index+axis_vector) && Neighbor_Visible(axis,index))
                    Update_Or_Add_Neighbor(phi_new,done,close_k,heap,heap_length,index+axis_vector);}
        else
            for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                if(index[axis] != domain.min_corner[axis] && !done(index-axis_vector))
                    Update_Or_Add_Neighbor(phi_new,done,close_k,heap,heap_length,index-axis_vector);
                if(index[axis] != domain.max_corner[axis] && !done(index+axis_vector))
                    Update_Or_Add_Neighbor(phi_new,done,close_k,heap,heap_length,index+axis_vector);}}

    RANGE<TV_INT> interior_domain=domain.Thickened(-ghost_cells);
    for(CELL_ITERATOR iterator(cell_grid,interior_domain);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();phi_ghost(cell)=phi_new(cell);}
}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class T_GRID> void FAST_MARCHING_METHOD_UNIFORM<T_GRID>::
Fast_Marching_Method(T_ARRAYS_SCALAR& phi_ghost,T_ARRAYS_BOOL& done,const T stopping_distance,const bool add_seed_indices_for_ghost_cells)
{
    int heap_length=0;
    T_ARRAYS_INT close_k(cell_grid.Domain_Indices(ghost_cells+1)); // extra cell so that it is the same size as the done array for optimizations
    ARRAY<TV_INT> heap(cell_grid.Domain_Indices().Thickened(ghost_cells).Size(),false); // a generic version of heap((m+6)*(n+6)*(mn+6),false);
    Initialize_Interface(phi_ghost,done,close_k,heap,heap_length,add_seed_indices_for_ghost_cells);

    while(heap_length != 0){
        TV_INT index=heap(1); // smallest point is on top of heap
        if(stopping_distance && abs(phi_ghost(index)) > stopping_distance){ // exit early
            for(CELL_ITERATOR iterator(cell_grid);iterator.Valid();iterator.Next()) if(!done(iterator.Cell_Index())){
                phi_ghost(iterator.Cell_Index())=LEVELSET_UTILITIES<T>::Sign(phi_ghost(iterator.Cell_Index()))*stopping_distance;}
            return;}
        done(index)=true;close_k(index)=0; // add to done, remove from close
        FAST_MARCHING<T>::Down_Heap(phi_ghost,close_k,heap,heap_length);heap_length--; // remove point from heap

        if(levelset.collision_body_list)
            for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                if(index[axis] != dimension_start[axis] && !done(index-axis_vector) && Neighbor_Visible(axis,index-axis_vector))
                    Update_Or_Add_Neighbor(phi_ghost,done,close_k,heap,heap_length,index-axis_vector);
                if(index[axis] != dimension_end[axis] && !done(index+axis_vector) && Neighbor_Visible(axis,index))
                    Update_Or_Add_Neighbor(phi_ghost,done,close_k,heap,heap_length,index+axis_vector);}
        else
            for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                if(index[axis] != dimension_start[axis] && !done(index-axis_vector))
                    Update_Or_Add_Neighbor(phi_ghost,done,close_k,heap,heap_length,index-axis_vector);
                if(index[axis] != dimension_end[axis] && !done(index+axis_vector))
                    Update_Or_Add_Neighbor(phi_ghost,done,close_k,heap,heap_length,index+axis_vector);}}
}
//#####################################################################
// Function Update_Or_Add_Neighbor
//#####################################################################
template<class T_GRID> inline void FAST_MARCHING_METHOD_UNIFORM<T_GRID>::
Update_Or_Add_Neighbor(T_ARRAYS_SCALAR& phi_ghost,T_ARRAYS_BOOL& done,T_ARRAYS_INT& close_k,ARRAY<TV_INT>& heap,int& heap_length,const TV_INT& neighbor)
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
template<class T_GRID> void FAST_MARCHING_METHOD_UNIFORM<T_GRID>::
Initialize_Interface(RANGE<TV_INT>& domain,T_ARRAYS_SCALAR& phi_ghost,T_ARRAYS_BOOL& done,T_ARRAYS_INT& close_k,ARRAY<TV_INT>& heap,int& heap_length,const ARRAY<TV_INT>* seed_indices,const bool add_seed_indices_for_ghost_cells)
{ 
    RANGE<TV_INT> interior_domain=domain.Thickened(-ghost_cells);
    if(seed_indices){
        for(int i=1;i<=seed_indices->m;i++)Add_To_Initial(done,close_k,(*seed_indices)(i));
        if(add_seed_indices_for_ghost_cells){
            for(int axis=1;axis<=TV::dimension;axis++) for(int side=1;side<=2;side++){
                RANGE<TV_INT> ghost_domain=domain;if(side==1) ghost_domain.max_corner(axis)=interior_domain.min_corner(axis)-1;else ghost_domain.min_corner(axis)=interior_domain.max_corner(axis)+1;
                for(CELL_ITERATOR iterator(cell_grid,ghost_domain);iterator.Valid();iterator.Next()){TV_INT index=iterator.Cell_Index();
                    for(int i=1;i<=T_GRID::number_of_neighbors_per_cell;i++){TV_INT neighbor_index(iterator.Cell_Neighbor(i));
                        if(domain.Lazy_Inside(neighbor_index) && LEVELSET_UTILITIES<T>::Interface(phi_ghost(index),phi_ghost(neighbor_index))){
                            if(!done(index))Add_To_Initial(done,close_k,index);
                            if(!done(neighbor_index))Add_To_Initial(done,close_k,neighbor_index);}}}}}}
    else{
        T_ARRAYS_SCALAR phi_new(domain.Thickened(1),false); // same size as done and close_k for array accelerations
        phi_new.Fill(2*cell_grid.dX.Max()); // ok positive, since minmag is used below
        if(levelset.collision_body_list){
            for(int axis=1;axis<=TV::dimension;axis++){
                domain.min_corner(axis)++;
                for(FACE_ITERATOR iterator(cell_grid,domain,axis);iterator.Valid();iterator.Next()){TV_INT index1=iterator.First_Cell_Index(),index2=iterator.Second_Cell_Index();
                    if(!Neighbor_Visible(iterator.Axis(),index1)){
                        if(phi_ghost(index1)<=0) Add_To_Initial(done,close_k,index1);
                        if(phi_ghost(index2)<=0) Add_To_Initial(done,close_k,index2);}
                    else if(LEVELSET_UTILITIES<T>::Interface(phi_ghost(index1),phi_ghost(index2))){
                        Add_To_Initial(done,close_k,index1);Add_To_Initial(done,close_k,index2);}}
                domain.min_corner(axis)--;}}
        else{
            for(int axis=1;axis<=TV::dimension;axis++){
                domain.min_corner(axis)++;
                for(FACE_ITERATOR iterator(cell_grid,domain,axis);iterator.Valid();iterator.Next()){TV_INT index1=iterator.First_Cell_Index(),index2=iterator.Second_Cell_Index();
                    if(LEVELSET_UTILITIES<T>::Interface(phi_ghost(index1),phi_ghost(index2))){
                        Add_To_Initial(done,close_k,index1);Add_To_Initial(done,close_k,index2);}}
                domain.min_corner(axis)--;}}

        Initialize_Interface_Threaded(domain,phi_ghost,phi_new,done);
        for(CELL_ITERATOR iterator(cell_grid,domain);iterator.Valid();iterator.Next()) if(done(iterator.Cell_Index())) phi_ghost(iterator.Cell_Index())=phi_new(iterator.Cell_Index());} // initialize done points

    // initialize close points
    for(CELL_ITERATOR iterator(cell_grid,domain);iterator.Valid();iterator.Next()) if(close_k(iterator.Cell_Index())){
        Update_Close_Point(phi_ghost,done,iterator.Cell_Index());
        heap_length++;heap(heap_length)=iterator.Cell_Index();
        FAST_MARCHING<T>::Up_Heap(phi_ghost,close_k,heap,heap_length);}
}
//#####################################################################
// Function Initialize_Interface
//#####################################################################
// pass heap_length by reference
template<class T_GRID> void FAST_MARCHING_METHOD_UNIFORM<T_GRID>::
Initialize_Interface(T_ARRAYS_SCALAR& phi_ghost,T_ARRAYS_BOOL& done,T_ARRAYS_INT& close_k,ARRAY<TV_INT>& heap,int& heap_length,const ARRAY<TV_INT>* seed_indices,const bool add_seed_indices_for_ghost_cells)
{ 
    if(seed_indices){
        for(int i=1;i<=seed_indices->m;i++)Add_To_Initial(done,close_k,(*seed_indices)(i));
        if(add_seed_indices_for_ghost_cells){RANGE<TV_INT> ghost_domain=cell_grid.Domain_Indices().Thickened(ghost_cells);
            for(CELL_ITERATOR iterator(cell_grid,ghost_cells,T_GRID::GHOST_REGION);iterator.Valid();iterator.Next()){TV_INT index=iterator.Cell_Index();
                for(int i=1;i<=T_GRID::number_of_neighbors_per_cell;i++){TV_INT neighbor_index(iterator.Cell_Neighbor(i));
                    if(ghost_domain.Lazy_Inside(neighbor_index) && LEVELSET_UTILITIES<T>::Interface(phi_ghost(index),phi_ghost(neighbor_index))){
                        if(!done(index))Add_To_Initial(done,close_k,index);
                        if(!done(neighbor_index))Add_To_Initial(done,close_k,neighbor_index);}}}}}
    else{
        T_ARRAYS_SCALAR phi_new(cell_grid.Domain_Indices(ghost_cells+1),false); // same size as done and close_k for array accelerations
        phi_new.Fill(2*cell_grid.dX.Max()); // ok positive, since minmag is used below
        if(levelset.collision_body_list){
            for(FACE_ITERATOR iterator(cell_grid,ghost_cells,T_GRID::INTERIOR_REGION);iterator.Valid();iterator.Next()){TV_INT index1=iterator.First_Cell_Index(),index2=iterator.Second_Cell_Index();
                if(!Neighbor_Visible(iterator.Axis(),index1)){
                    if(phi_ghost(index1)<=0) Add_To_Initial(done,close_k,index1);
                    if(phi_ghost(index2)<=0) Add_To_Initial(done,close_k,index2);}
                else if(LEVELSET_UTILITIES<T>::Interface(phi_ghost(index1),phi_ghost(index2))){
                    Add_To_Initial(done,close_k,index1);Add_To_Initial(done,close_k,index2);}}}
        else{
            for(FACE_ITERATOR iterator(cell_grid,ghost_cells,T_GRID::INTERIOR_REGION);iterator.Valid();iterator.Next()){TV_INT index1=iterator.First_Cell_Index(),index2=iterator.Second_Cell_Index();
                if(LEVELSET_UTILITIES<T>::Interface(phi_ghost(index1),phi_ghost(index2))){
                    Add_To_Initial(done,close_k,index1);Add_To_Initial(done,close_k,index2);}}}

        DOMAIN_ITERATOR_THREADED_ALPHA<FAST_MARCHING_METHOD_UNIFORM<T_GRID>,TV>(cell_grid.Domain_Indices(ghost_cells),thread_queue).template Run<T_ARRAYS_SCALAR&,T_ARRAYS_SCALAR&,T_ARRAYS_BOOL&>(*this,&FAST_MARCHING_METHOD_UNIFORM<T_GRID>::Initialize_Interface_Threaded,phi_ghost,phi_new,done);

        for(CELL_ITERATOR iterator(cell_grid,ghost_cells);iterator.Valid();iterator.Next()) if(done(iterator.Cell_Index())) phi_ghost(iterator.Cell_Index())=phi_new(iterator.Cell_Index());} // initialize done points

    // initialize close points
    for(CELL_ITERATOR iterator(cell_grid,ghost_cells);iterator.Valid();iterator.Next()) if(close_k(iterator.Cell_Index())){
        Update_Close_Point(phi_ghost,done,iterator.Cell_Index());
        heap_length++;heap(heap_length)=iterator.Cell_Index();
        FAST_MARCHING<T>::Up_Heap(phi_ghost,close_k,heap,heap_length);}
}
//#####################################################################
// Function Initialize_Interface
//#####################################################################
// pass heap_length by reference
template<class T_GRID> void FAST_MARCHING_METHOD_UNIFORM<T_GRID>::
Initialize_Interface_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_SCALAR& phi_ghost,T_ARRAYS_SCALAR& phi_new,T_ARRAYS_BOOL& done)
{ 
    T_LEVELSET levelset_ghost(cell_grid,phi_ghost);

    T sqr_epsilon_cell_size=sqr(levelset.small_number*cell_grid.Cell_Size()),sqr_epsilon_face_size[3];
    if(T_GRID::dimension==2){sqr_epsilon_face_size[2]=sqr_epsilon_cell_size;}
    else if(T_GRID::dimension==3) for(int axis=1;axis<=3;axis++) sqr_epsilon_face_size[axis-1]=sqr_epsilon_cell_size/sqr(cell_grid.dX[axis]);
    
    T fmm_initialization_iterative_tolerance_absolute=levelset.fmm_initialization_iterative_tolerance*cell_grid.Minimum_Edge_Length();    
 
    if(levelset.collision_body_list){
        COLLISION_GEOMETRY_ID body_id;
        for(CELL_ITERATOR iterator(cell_grid,domain);iterator.Valid();iterator.Next()) if(done(iterator.Cell_Index())){TV_INT index=iterator.Cell_Index();
            T value[T_GRID::dimension]={0}; // the phi value to use in the given direction
            int number_of_axis=0; // the number of axis that we want to use later
            int missing_axis=3; // used in number_of_axis==2 case only, so it gives you which axis is missing (==3 for 2d)
            bool really_clamp_phi_with_collision_bodies=levelset.clamp_phi_with_collision_bodies&&phi_ghost(index)<=0;
            T abs_phi=abs(phi_ghost(index));
            TV location=iterator.Location();
            for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis),low=index-axis_vector,high=index+axis_vector;T dx=cell_grid.dX[axis];
                bool use_low=false,use_high=false;T value_low=0,value_high=0;
                if(index[axis]>dimension_start[axis]){
                    if(!Neighbor_Visible(axis,low)){
                        RAY<TV> ray(location,-TV::Axis_Vector(axis),true);ray.t_max=dx;ray.semi_infinite=false;
                        levelset.collision_body_list->Closest_Non_Intersecting_Point_Of_Any_Body(ray,body_id);
                        T phi_object=max(ray.t_max,levelset.small_number);
                        use_low=true;if(really_clamp_phi_with_collision_bodies) value_low=min(abs_phi,phi_object);else value_low=phi_object;}
                    else if(done(low)&&LEVELSET_UTILITIES<T>::Interface(phi_ghost(index),phi_ghost(low))){
                        use_low=true;value_low=LEVELSET_UTILITIES<T>::Theta(phi_ghost(index),phi_ghost(low))*dx;}}
                if(index[axis]<dimension_end[axis]){
                    if(!Neighbor_Visible(axis,index)){
                        RAY<TV> ray(location,TV::Axis_Vector(axis),true);ray.t_max=dx;ray.semi_infinite=false;
                        levelset.collision_body_list->Closest_Non_Intersecting_Point_Of_Any_Body(ray,body_id);
                        T phi_object=max(ray.t_max,levelset.small_number);
                        use_high=true;if(really_clamp_phi_with_collision_bodies) value_high=min(abs_phi,phi_object);else value_high=phi_object;}
                    else if(done(high)&&LEVELSET_UTILITIES<T>::Interface(phi_ghost(index),phi_ghost(high))){
                        use_high=true;value_high=LEVELSET_UTILITIES<T>::Theta(phi_ghost(index),phi_ghost(high))*dx;}}
                if(!use_low){if(use_high) value[number_of_axis]=value_high;else{missing_axis=axis;number_of_axis--;}}
                else if(!use_high) value[number_of_axis]=value_low;
                else value[number_of_axis]=min(value_low,value_high);
                number_of_axis++;}
            assert(number_of_axis);
            if(number_of_axis==1) phi_new(index)=value[0];
            else if(number_of_axis==2){T d2=sqr(value[0])+sqr(value[1]);
                if(d2 > sqr_epsilon_face_size[missing_axis-1]) phi_new(index)=value[0]*value[1]/sqrt(d2);else phi_new(index)=min(value[0],value[1]);}
            else{PHYSBAM_ASSERT(T_GRID::dimension==3); // 2d should never get to this point
                T value_xy=value[0]*value[1],value_xz=value[0]*value[2],value_yz=value[1]*value[2],d2=sqr(value_xy)+sqr(value_xz)+sqr(value_yz);
                if(d2 > sqr_epsilon_cell_size) phi_new(index)=value_xy*value[2]/sqrt(d2);else phi_new(index)=min(value[0],value[1],value[2]);}
            if(levelset.refine_fmm_initialization_with_iterative_solver){
                TV vertex=location;T phi_value=phi_ghost(index);
                for(int iterations=1;iterations<=levelset.fmm_initialization_iterations;iterations++){
                    if(abs(phi_value)>10*cell_grid.Minimum_Edge_Length())break; // stop if it looks like it's blowing up
                    vertex-=phi_value*levelset_ghost.Normal(vertex);phi_value=levelset_ghost.Phi(vertex);
                    if(abs(phi_value)<=fmm_initialization_iterative_tolerance_absolute){
                        T new_phi_value=(vertex-location).Magnitude();
                        if((new_phi_value-phi_new(index))/max(fmm_initialization_iterative_tolerance_absolute,phi_new(index))<levelset.fmm_initialization_iterative_drift_fraction && new_phi_value>0)
                            phi_new(index)=new_phi_value;
                        break;}}}
            phi_new(index)*=LEVELSET_UTILITIES<T>::Sign(phi_ghost(index));}}
    else{
        for(CELL_ITERATOR iterator(cell_grid,domain);iterator.Valid();iterator.Next()) if(done(iterator.Cell_Index())){TV_INT index=iterator.Cell_Index();
            T value[T_GRID::dimension]={0}; // the phi value to use in the given direction
            int number_of_axis=0; // the number of axis that we want to use later
            int missing_axis=3; // used in number_of_axis==2 case only, so it gives you which axis is missing (==3 for 2d)
            TV location=iterator.Location();
            for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis),low=index-axis_vector,high=index+axis_vector;T dx=cell_grid.dX[axis];
                bool use_low=(done(low)&&LEVELSET_UTILITIES<T>::Interface(phi_ghost(index),phi_ghost(low))),use_high=(done(high)&&LEVELSET_UTILITIES<T>::Interface(phi_ghost(index),phi_ghost(high)));
                if(!use_low){if(use_high)value[number_of_axis]=LEVELSET_UTILITIES<T>::Theta(phi_ghost(index),phi_ghost(high))*dx;else{missing_axis=axis;number_of_axis--;}}
                else if(!use_high) value[number_of_axis]=LEVELSET_UTILITIES<T>::Theta(phi_ghost(index),phi_ghost(low))*dx;
                else value[number_of_axis]=min(LEVELSET_UTILITIES<T>::Theta(phi_ghost(index),phi_ghost(high)),LEVELSET_UTILITIES<T>::Theta(phi_ghost(index),phi_ghost(low)))*dx;
                number_of_axis++;}
            if(number_of_axis==1) phi_new(index)=value[0];
            else if(number_of_axis==2){T d2=sqr(value[0])+sqr(value[1]);if(d2 > sqr_epsilon_face_size[missing_axis-1]) phi_new(index)=value[0]*value[1]/sqrt(d2);else phi_new(index)=min(value[0],value[1]);}
            else{PHYSBAM_ASSERT(T_GRID::dimension==3); // 2d should never get to this point
                T value_xy=value[0]*value[1],value_xz=value[0]*value[2],value_yz=value[1]*value[2],d2=sqr(value_xy)+sqr(value_xz)+sqr(value_yz);
                if(d2 > sqr_epsilon_cell_size) phi_new(index)=value_xy*value[2]/sqrt(d2);else phi_new(index)=min(value[0],value[1],value[2]);}
            if(levelset.refine_fmm_initialization_with_iterative_solver){TV location=iterator.Location();
                TV vertex=location;T phi_value=phi_ghost(index);
                for(int iterations=1;iterations<=levelset.fmm_initialization_iterations;iterations++){
                    if(abs(phi_value)>10*cell_grid.Minimum_Edge_Length())break; // stop if it looks like it's blowing up
                    vertex-=phi_value*levelset_ghost.Normal(vertex);phi_value=levelset_ghost.Phi(vertex);
                    if(abs(phi_value)<=fmm_initialization_iterative_tolerance_absolute){
                        T new_phi_value=(vertex-location).Magnitude();
                        if((new_phi_value-phi_new(index))/max(fmm_initialization_iterative_tolerance_absolute,phi_new(index))<levelset.fmm_initialization_iterative_drift_fraction && new_phi_value>0)
                            phi_new(index)=new_phi_value;
                        break;}}}
            phi_new(index)*=LEVELSET_UTILITIES<T>::Sign(phi_ghost(index));}}
}
//#####################################################################
// Function Initialize_Interface
//#####################################################################
// pass heap_length by reference
template<class T_GRID> void FAST_MARCHING_METHOD_UNIFORM<T_GRID>::
Initialize_Interface(T_ARRAYS_SCALAR& phi_ghost,T_ARRAYS_BOOL& done,T_ARRAYS_INT& close_k,ARRAY<TV_INT>& heap,int& heap_length,const bool add_seed_indices_for_ghost_cells)
{
    T_LEVELSET levelset_ghost(cell_grid,phi_ghost);

    for(CELL_ITERATOR iterator(cell_grid,ghost_cells);iterator.Valid();iterator.Next()) if(done(iterator.Cell_Index())) Add_To_Initial(done,close_k,iterator.Cell_Index());
    if(add_seed_indices_for_ghost_cells){RANGE<TV_INT> ghost_domain=cell_grid.Domain_Indices().Thickened(ghost_cells);
        for(CELL_ITERATOR iterator(cell_grid,ghost_cells,T_GRID::GHOST_REGION);iterator.Valid();iterator.Next()){TV_INT index=iterator.Cell_Index();
            for(int i=1;i<=T_GRID::number_of_neighbors_per_cell;i++){TV_INT neighbor_index(iterator.Cell_Neighbor(i));
                if(ghost_domain.Lazy_Inside(neighbor_index) && LEVELSET_UTILITIES<T>::Interface(phi_ghost(index),phi_ghost(neighbor_index))){
                    if(!done(index))Add_To_Initial(done,close_k,index);
                    if(!done(neighbor_index))Add_To_Initial(done,close_k,neighbor_index);}}}}

   // initialize close points
    for(CELL_ITERATOR iterator(cell_grid,ghost_cells);iterator.Valid();iterator.Next()) if(close_k(iterator.Cell_Index())){
        Update_Close_Point(phi_ghost,done,iterator.Cell_Index());
        heap_length++;heap(heap_length)=iterator.Cell_Index();
        FAST_MARCHING<T>::Up_Heap(phi_ghost,close_k,heap,heap_length);}
}
//#####################################################################
// Function Update_Close_Point
//##################################################################### 
// needs done=0 around the outside of the domain
template<class T_GRID> void FAST_MARCHING_METHOD_UNIFORM<T_GRID>::
Update_Close_Point(T_ARRAYS_SCALAR& phi_ghost,const T_ARRAYS_BOOL& done,const TV_INT& index)
{
    T value[T_GRID::dimension]={}; // the phi value to use in the given direction
    T dx[T_GRID::dimension]; // the edge length in the given direction
    int number_of_axis=0; // the number of axis that we want to use later

    // check each principal axis
    for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis),low=index-axis_vector,high=index+axis_vector;
        bool check_low=done(low),check_high=done(high);
        if(levelset.collision_body_list){
            if(check_low && !Neighbor_Visible(axis,low)) check_low=false;
            if(check_high && !Neighbor_Visible(axis,index)) check_high=false;}
        dx[number_of_axis]=cell_grid.dX[axis];
        if(!check_low){
            if(check_high)value[number_of_axis]=phi_ghost(high);
            else number_of_axis--;}
        else if(!check_high) value[number_of_axis]=phi_ghost(low);
        else value[number_of_axis]=minmag(phi_ghost(low),phi_ghost(high));
        number_of_axis++;}

    phi_ghost(index)=FAST_MARCHING<T>::template Solve_Close_Point<T_GRID::dimension>(phi_ghost(index),number_of_axis,value,dx);
}
//#####################################################################
// Function Add_To_Initial
//##################################################################### 
template<class T_GRID> void FAST_MARCHING_METHOD_UNIFORM<T_GRID>::
Add_To_Initial(T_ARRAYS_BOOL& done,T_ARRAYS_INT& close_k,const TV_INT& index)
{
    done(index)=true;close_k(index)=0; // add to done, remove from close 
    // add neighbors to close if not done 
    if(levelset.collision_body_list){
        for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
            if(!done(index-axis_vector) && Neighbor_Visible(axis,index-axis_vector)) close_k(index-axis_vector)=1;
            if(!done(index+axis_vector) && Neighbor_Visible(axis,index)) close_k(index+axis_vector)=1;}}
    else{
        for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
            if(!done(index-axis_vector)) close_k(index-axis_vector)=1;
            if(!done(index+axis_vector)) close_k(index+axis_vector)=1;}}
}
//#####################################################################
template class FAST_MARCHING_METHOD_UNIFORM<GRID<VECTOR<float,1> > >;
template class FAST_MARCHING_METHOD_UNIFORM<GRID<VECTOR<float,2> > >;
template class FAST_MARCHING_METHOD_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FAST_MARCHING_METHOD_UNIFORM<GRID<VECTOR<double,1> > >;
template class FAST_MARCHING_METHOD_UNIFORM<GRID<VECTOR<double,2> > >;
template class FAST_MARCHING_METHOD_UNIFORM<GRID<VECTOR<double,3> > >;
#endif

//#####################################################################
// Copyright 2011, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace CUT_CELL_COMPUTATIONS
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/OBJECTS_IN_CELL.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/CUT_CELL.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
namespace PhysBAM{
namespace CUT_CELL_COMPUTATIONS{
//#####################################################################
// Helper Functions
//#####################################################################
namespace {
template<class T,int d>
bool Is_Occluded_Cell(VECTOR<T,d> centroid,VECTOR<T,d> cell_center,OBJECTS_IN_CELL<GRID<VECTOR<T,d> >,COLLISION_GEOMETRY_ID>& objects_in_cell,const COLLISION_GEOMETRY_COLLECTION<VECTOR<T,d> >& collision_geometry_collection,VECTOR<int,d>& cell_index,VECTOR<int,d>& neighbor_index)
{ARRAY<COLLISION_GEOMETRY_ID> collision_objects;objects_in_cell.Get_Objects_For_Cells(cell_index,neighbor_index,collision_geometry_collection.bodies.m,collision_objects);return collision_geometry_collection.Intersection_Between_Points(centroid,cell_center,&collision_objects);}
template<class T,int d>
bool Is_Occluded_Cell_Center(VECTOR<T,d> centroid,VECTOR<T,d> cell_center,OBJECTS_IN_CELL<GRID<VECTOR<T,d> >,COLLISION_GEOMETRY_ID>& objects_in_cell,const COLLISION_GEOMETRY_COLLECTION<VECTOR<T,d> >& collision_geometry_collection,VECTOR<int,d>& cell_index)
{ARRAY<COLLISION_GEOMETRY_ID> collision_objects;objects_in_cell.Get_Objects_For_Cell(cell_index,collision_objects);return collision_geometry_collection.Intersection_Between_Points(centroid,cell_center,&collision_objects);}
}
//#####################################################################
// Function Compute_Cut_Geometries 1-D
//#####################################################################
template<class T>
void Compute_Cut_Geometries(const GRID<VECTOR<T,1> >& grid,const int num_ghost_cells,typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<T,1> > >::GRID_BASED_COLLISION_GEOMETRY& collision_bodies_affecting_fluid,ARRAY<CUT_CELLS<T,1>*,VECTOR<int,1> >& cut_cells){
    typedef VECTOR<T,1> TV;
    typedef VECTOR<int,1> TV_INT;
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid,num_ghost_cells);iterator.Valid();iterator.Next()){
        TV_INT index=iterator.Cell_Index();COLLISION_GEOMETRY_ID body_id;
        ARRAY<COLLISION_GEOMETRY_ID> collision_objects;collision_bodies_affecting_fluid.objects_in_cell.Get_Objects_For_Cell(index,collision_objects);
        if(!collision_objects.Size()) cut_cells(index)=0;
        else{
            cut_cells(index) = new CUT_CELLS<T,1>();cut_cells(index)->dominant_element=0;
            RANGE<TV> cell_volume=iterator.Bounding_Box();
            RAY<TV> l_to_r_ray(cell_volume.min_corner,cell_volume.max_corner-cell_volume.min_corner,true);l_to_r_ray.t_max=l_to_r_ray.direction.Normalize();l_to_r_ray.semi_infinite=false;
            RAY<TV> r_to_l_ray(cell_volume.max_corner,cell_volume.min_corner-cell_volume.max_corner,true);r_to_l_ray.t_max=r_to_l_ray.direction.Normalize();r_to_l_ray.semi_infinite=false;
            if(collision_bodies_affecting_fluid.Intersection_With_Any_Simplicial_Object(l_to_r_ray,body_id,&collision_objects)){
                int poly_index=cut_cells(index)->geometry.Append(POLYGON<TV>(RANGE<TV>(cell_volume.min_corner, l_to_r_ray.Point(l_to_r_ray.t_max))));
                TV centroid=ARRAYS_COMPUTATIONS::Average(cut_cells(index)->geometry(poly_index).X);

                cut_cells(index)->visibility.Resize(poly_index);
                if(!Is_Occluded_Cell_Center<T,1>(centroid,grid.Center(index),collision_bodies_affecting_fluid.objects_in_cell,collision_bodies_affecting_fluid.collision_geometry_collection,index)){
                    cut_cells(index)->dominant_element=poly_index;
                    cut_cells(index)->visibility(poly_index).Append(index);}
                if(!cut_cells(index)->dominant_element)
                    for(int node=1;!cut_cells(index)->dominant_element && node<=cut_cells(index)->geometry(poly_index).X.Size();++node)
                        if(!Is_Occluded_Cell_Center<T,1>(cut_cells(index)->geometry(poly_index).X(node),grid.Center(index),collision_bodies_affecting_fluid.objects_in_cell,collision_bodies_affecting_fluid.collision_geometry_collection,index)){
                            cut_cells(index)->dominant_element=poly_index;
                            cut_cells(index)->visibility(poly_index).Append(index);}

                TV_INT neighbor_cell_index=index-TV_INT::Axis_Vector(1);
                if(!Is_Occluded_Cell<T,1>(centroid,grid.Center(neighbor_cell_index),collision_bodies_affecting_fluid.objects_in_cell,collision_bodies_affecting_fluid.collision_geometry_collection,index,neighbor_cell_index))
                    cut_cells(index)->visibility(poly_index).Append(neighbor_cell_index);}
            if(collision_bodies_affecting_fluid.Intersection_With_Any_Simplicial_Object(r_to_l_ray,body_id,&collision_objects)){
                int poly_index=cut_cells(index)->geometry.Append(POLYGON<TV>(RANGE<TV>(r_to_l_ray.Point(r_to_l_ray.t_max),cell_volume.max_corner)));
                TV centroid=ARRAYS_COMPUTATIONS::Average(cut_cells(index)->geometry(poly_index).X);

                cut_cells(index)->visibility.Resize(poly_index);
                if(!Is_Occluded_Cell_Center<T,1>(centroid,grid.Center(index),collision_bodies_affecting_fluid.objects_in_cell,collision_bodies_affecting_fluid.collision_geometry_collection,index)){
                    cut_cells(index)->dominant_element=poly_index;
                    cut_cells(index)->visibility(poly_index).Append(index);}
                if(!cut_cells(index)->dominant_element)
                    for(int node=1;!cut_cells(index)->dominant_element && node<=cut_cells(index)->geometry(poly_index).X.Size();++node)
                        if(!Is_Occluded_Cell_Center<T,1>(cut_cells(index)->geometry(poly_index).X(node),grid.Center(index),collision_bodies_affecting_fluid.objects_in_cell,collision_bodies_affecting_fluid.collision_geometry_collection,index)){
                            cut_cells(index)->dominant_element=poly_index;
                            cut_cells(index)->visibility(poly_index).Append(index);}

                TV_INT neighbor_cell_index=index+TV_INT::Axis_Vector(1);
                if(!Is_Occluded_Cell<T,1>(centroid,grid.Center(neighbor_cell_index),collision_bodies_affecting_fluid.objects_in_cell,collision_bodies_affecting_fluid.collision_geometry_collection,index,neighbor_cell_index))
                    cut_cells(index)->visibility(poly_index).Append(neighbor_cell_index);}
            if(!cut_cells(index)->geometry.Size()){delete cut_cells(index);cut_cells(index)=0;}}}
}
//#####################################################################
// Function Compute_Cut_Geometries 2-D
//#####################################################################
template<class T>
void Compute_Cut_Geometries(const GRID<VECTOR<T,2> >& grid,const int num_ghost_cells,typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<T,2> > >::GRID_BASED_COLLISION_GEOMETRY& collision_bodies_affecting_fluid,ARRAY<CUT_CELLS<T,2>*,VECTOR<int,2> >& cut_cells)
{
    PHYSBAM_NOT_IMPLEMENTED();
    // ARRAY<TV_INT> neighbor_offsets(8);
    // neighbor_offsets(1)=TV_INT( 1,-1); neighbor_offsets(2)=TV_INT( 1, 0); neighbor_offsets(3)=TV_INT( 1, 1);
    // neighbor_offsets(4)=TV_INT( 0,-1);                                    neighbor_offsets(5)=TV_INT( 0, 1);
    // neighbor_offsets(6)=TV_INT(-1,-1); neighbor_offsets(7)=TV_INT(-1, 0); neighbor_offsets(8)=TV_INT(-1, 1);

    // typedef VECTOR<T,1> TV;
    // typedef VECTOR<int,1> TV_INT;
    // for(typename GRID<TV>::CELL_ITERATOR iterator(grid,num_ghost_cells);iterator.Valid();iterator.Next()){
    //     TV_INT index=iterator.Cell_Index();COLLISION_GEOMETRY_ID body_id;
    //     ARRAY<COLLISION_GEOMETRY_ID> collision_objects;collision_bodies_affecting_fluid.objects_in_cell.Get_Objects_For_Cell(index,collision_objects);
    //     if(!collision_objects.Size()) cut_cells(index)=0;
    //     else{
    //         LOG::cout<<"There are "<<collision_objects.Size()<<" objects which cut the cell of interest"<<std::endl;
    //         cut_cells(index) = new CUT_CELLS<T,2>();cut_cells(index)->dominant_element=0;
    //         POLYGON<TV> full_cell_volume(iterator.Bounding_Box());

    //         if(!cut_cells(index)->geometry.Size()){delete cut_cells(index);cut_cells(index)=0;}}}
}
//#####################################################################
// Function Compute_Cut_Geometries 3-D
//#####################################################################
template<class T>
void Compute_Cut_Geometries(const GRID<VECTOR<T,3> >& grid,const int num_ghost_cells,typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<T,3> > >::GRID_BASED_COLLISION_GEOMETRY& collision_bodies_affecting_fluid,ARRAY<CUT_CELLS<T,3>*,VECTOR<int,3> >& cut_cells)
{PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
template void Compute_Cut_Geometries(const GRID<VECTOR<float,1> >&,const int,typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<float,1> > >::GRID_BASED_COLLISION_GEOMETRY&,ARRAY<CUT_CELLS<float,1>*,VECTOR<int,1> >&);
template void Compute_Cut_Geometries(const GRID<VECTOR<float,2> >&,const int,typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<float,2> > >::GRID_BASED_COLLISION_GEOMETRY&,ARRAY<CUT_CELLS<float,2>*,VECTOR<int,2> >&);
template void Compute_Cut_Geometries(const GRID<VECTOR<float,3> >&,const int,typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<float,3> > >::GRID_BASED_COLLISION_GEOMETRY&,ARRAY<CUT_CELLS<float,3>*,VECTOR<int,3> >&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void Compute_Cut_Geometries(const GRID<VECTOR<double,1> >&,const int,typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<double,1> > >::GRID_BASED_COLLISION_GEOMETRY&,ARRAY<CUT_CELLS<double,1>*,VECTOR<int,1> >&);
template void Compute_Cut_Geometries(const GRID<VECTOR<double,2> >&,const int,typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<double,2> > >::GRID_BASED_COLLISION_GEOMETRY&,ARRAY<CUT_CELLS<double,2>*,VECTOR<int,2> >&);
template void Compute_Cut_Geometries(const GRID<VECTOR<double,3> >&,const int,typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<double,3> > >::GRID_BASED_COLLISION_GEOMETRY&,ARRAY<CUT_CELLS<double,3>*,VECTOR<int,3> >&);
#endif
}
}

//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __POTENTIAL_COLLISIONS__
#define __POTENTIAL_COLLISIONS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>

#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>

namespace PhysBAM
{

template<class TV_INT> class RANGE;

//TODO: This should be in GET_POTENTIAL_COLLISIONS.cpp but doing so causes gcc to crash
namespace POTENTIAL_COLLISIONS
{
    extern OPERATION_HASH<COLLISION_GEOMETRY_ID> already_added;
    template<class T_COLLISION_GEOMETRY,class T_ARRAY,class ID> 
    void Get_Potential_Collisions(const typename T_COLLISION_GEOMETRY::VECTOR_T::template REBIND<int>::TYPE& voxel,const HASHTABLE<typename T_COLLISION_GEOMETRY::VECTOR_T::template REBIND<int>::TYPE,ARRAY<ID>*>& hashtable,const ARRAY<ID>& bodies_not_in_partition,T_ARRAY& collision_bodies,ARRAY<ID>& object_indices)
    {
        object_indices.Remove_All();
        ARRAY<ID>* occupancy_list=0;
        if(hashtable.Get(voxel,occupancy_list)) for(int t=1;t<=occupancy_list->m;t++) object_indices.Append((*occupancy_list)(t));
        for(int i=1;i<=bodies_not_in_partition.m;i++) object_indices.Append(bodies_not_in_partition(i));
    }
    template<class T_COLLISION_GEOMETRY,class T_ARRAY,class ID>
    void Get_Potential_Collisions(const typename T_COLLISION_GEOMETRY::VECTOR_T::template REBIND<int>::TYPE& voxel,const HASHTABLE<typename T_COLLISION_GEOMETRY::VECTOR_T::template REBIND<int>::TYPE,ARRAY<ID>*>& hashtable,const ARRAY<ID>& bodies_not_in_partition,T_ARRAY& collision_bodies,ARRAY<T_COLLISION_GEOMETRY*>& objects)
    {
        objects.Remove_All();
        ARRAY<ID>* occupancy_list=0;
        if(hashtable.Get(voxel,occupancy_list)) for(int t=1;t<=occupancy_list->m;t++) objects.Append(collision_bodies((*occupancy_list)(t)));
        for(int i=1;i<=bodies_not_in_partition.m;i++) objects.Append(collision_bodies(bodies_not_in_partition(i)));
    }
    template<class T_COLLISION_GEOMETRY,class T_ARRAY,class ID>
    void Get_Potential_Collisions(const ID index,const RANGE<typename T_COLLISION_GEOMETRY::VECTOR_T::template REBIND<int>::TYPE>& range,const HASHTABLE<typename T_COLLISION_GEOMETRY::VECTOR_T::template REBIND<int>::TYPE,ARRAY<ID>*>& hashtable,const ARRAY<ID>& bodies_not_in_partition,T_ARRAY& collision_bodies,ARRAY<ID>& object_indices,bool only_higher_index)
    {
        typedef typename T_COLLISION_GEOMETRY::VECTOR_T TV;
        typedef typename TV::template REBIND<int>::TYPE TV_INT;

        object_indices.Remove_All();
        already_added.Initialize(collision_bodies.Size());
        GRID<TV> unused;
        for(typename GRID<TV>::CELL_ITERATOR iterator(unused,range);iterator.Valid();iterator.Next()){TV_INT voxel=iterator.Cell_Index();
            ARRAY<ID>* occupancy_list=0;
            if(hashtable.Get(voxel,occupancy_list)) for(int t=1;t<=occupancy_list->m;t++){
                ID k=(*occupancy_list)(t);
                if(k>index || (!only_higher_index && k!=index))
                    if(!already_added.Is_Marked_Current(k)){already_added.Mark(k);object_indices.Append(k);}}}
        for(int i=1;i<=bodies_not_in_partition.m;i++){ // add bodies that aren't in spatial partition - don't need append unique
            ID k=bodies_not_in_partition(i);if(k>index || (!only_higher_index && k!=index)) object_indices.Append(k);}
    }
    template<class T_COLLISION_GEOMETRY,class T_ARRAY,class ID>
    void Get_Potential_Collisions(const ID index,const RANGE<typename T_COLLISION_GEOMETRY::VECTOR_T::template REBIND<int>::TYPE>& range,const HASHTABLE<typename T_COLLISION_GEOMETRY::VECTOR_T::template REBIND<int>::TYPE,ARRAY<ID>*>& hashtable,const ARRAY<ID>& bodies_not_in_partition,T_ARRAY& collision_bodies,ARRAY<T_COLLISION_GEOMETRY*>& objects,bool only_higher_index)
    {
        typedef typename T_COLLISION_GEOMETRY::VECTOR_T TV;
        typedef typename TV::template REBIND<int>::TYPE TV_INT;
        
        objects.Remove_All();
        already_added.Initialize(collision_bodies.Size());
        GRID<TV> unused;
        for(typename GRID<TV>::CELL_ITERATOR iterator(unused,range);iterator.Valid();iterator.Next()){TV_INT voxel=iterator.Cell_Index();
            ARRAY<ID>* occupancy_list=0;
            if(hashtable.Get(voxel,occupancy_list)) for(int t=1;t<=occupancy_list->m;t++){
                ID k=(*occupancy_list)(t);
                if(k>index || (!only_higher_index && k!=index))
                    if(!already_added.Is_Marked_Current(k)){already_added.Mark(k);objects.Append(collision_bodies(k));}}}
        for(int i=1;i<=bodies_not_in_partition.m;i++){ // add bodies that aren't in spatial partition - don't need append unique
            ID k=bodies_not_in_partition(i);if(k>index || (!only_higher_index && k!=index)) objects.Append(collision_bodies(k));}
    }
}
}
#endif

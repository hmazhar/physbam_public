//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COLLISION_GEOMETRY_SPATIAL_PARTITION
//#####################################################################
#ifndef __COLLISION_GEOMETRY_SPATIAL_PARTITION__
#define __COLLISION_GEOMETRY_SPATIAL_PARTITION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Collision_Detection_Computations/GET_POTENTIAL_COLLISIONS.h>
#include <PhysBAM_Geometry/Collisions/COLLISIONS_GEOMETRY_FORWARD.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <climits> // for INT_MAX (for old gcc)

namespace PhysBAM{

template<class T_COLLISION_GEOMETRY,class T_ARRAY,class ID>
class COLLISION_GEOMETRY_SPATIAL_PARTITION:public NONCOPYABLE
{
private:
    typedef typename T_COLLISION_GEOMETRY::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;

    T_ARRAY& collision_bodies;
    ARRAY<RANGE<TV_INT>,ID> voxel_range;
    HASHTABLE<TV_INT,ARRAY<ID>*> hashtable;
    int reinitialize_counter;
public:
    T voxel_size,one_over_voxel_size;
    ARRAY<ID> bodies_not_in_partition;
    T collision_body_thickness;

    COLLISION_GEOMETRY_SPATIAL_PARTITION(T_ARRAY& collision_bodies_input,T collision_body_thickness_input=0);
    ~COLLISION_GEOMETRY_SPATIAL_PARTITION();

    void Get_Potential_Collisions(const ID index,ARRAY<ID>& object_indices,bool only_higher_index=false) const
    {POTENTIAL_COLLISIONS::Get_Potential_Collisions<T_COLLISION_GEOMETRY,T_ARRAY,ID>(index,voxel_range(index),hashtable,bodies_not_in_partition,collision_bodies,object_indices,only_higher_index);}

    void Get_Potential_Collisions_Using_Current_Position(const ID index,ARRAY<ID>& object_indices,bool only_higher_index=false) const
    {POTENTIAL_COLLISIONS::Get_Potential_Collisions<T_COLLISION_GEOMETRY,T_ARRAY,ID>(index,Voxel_Range(index),hashtable,bodies_not_in_partition,collision_bodies,object_indices,only_higher_index);}

    template<class T_RETURN>
    void Get_Potential_Collisions(const TV& location,ARRAY<T_RETURN>& objects) const
    {POTENTIAL_COLLISIONS::Get_Potential_Collisions<T_COLLISION_GEOMETRY,T_ARRAY,ID>(Voxel(location),hashtable,bodies_not_in_partition,collision_bodies,objects);}

    template<class T_RETURN>
    void Get_Potential_Collisions(const RANGE<TV>& box,ARRAY<T_RETURN>& objects) const
    {POTENTIAL_COLLISIONS::Get_Potential_Collisions<T_COLLISION_GEOMETRY,T_ARRAY,ID>(ID(-1),Voxel(box),hashtable,bodies_not_in_partition,collision_bodies,objects,false);}

    void Set_Collision_Body_Thickness(const T collision_body_thickness_input)
    {collision_body_thickness=collision_body_thickness_input;Reinitialize();}

    ID Number_Of_Bodies()
    {return collision_bodies.Size();}

private:
    bool Out_Of_Range(const TV& location) const
    {for(int i=1;i<=TV::m;i++) if(location(i)*one_over_voxel_size<=INT_MIN/2||location(i)*one_over_voxel_size>=INT_MAX/2) return true; return false;}

    RANGE<TV_INT> Voxel_Range(const ID index) const
    {return Voxel(collision_bodies(index)->Axis_Aligned_Bounding_Box().Thickened(collision_body_thickness));}

    TV_INT Voxel(const TV& location) const
    {if(Out_Of_Range(location)) throw FLOATING_POINT_ERROR(std::string("Floating Point Exception: Location is too far away or voxel size is too small"));
    return TV_INT(floor(location*one_over_voxel_size));}

    RANGE<TV_INT> Voxel(const RANGE<TV>& box) const
    {if(Out_Of_Range(box.min_corner) || Out_Of_Range(box.max_corner)) throw FLOATING_POINT_ERROR(std::string("Floating Point Exception: Box is too big or voxel size is too small"));
    return RANGE<TV_INT>(Voxel(box.min_corner),Voxel(box.max_corner));}

    int Number_Of_Voxels_Occupied(const ID index) const
    {return (voxel_range(index).Edge_Lengths()+TV_INT::All_Ones_Vector()).Product();}

public:
//#####################################################################
    RANGE<TV> Scene_Bounding_Box();
    T Scene_Bounding_Box_Size();
    T Average_Bounding_Box_Size();
    T Maximum_Bounding_Box_Size();
    void Compute_Voxel_Size(const SPATIAL_PARTITION_VOXEL_SIZE_HEURISTIC heuristic,const int number_of_boxes,const T voxel_size_scale_factor);
    void Reinitialize();
    void Print_Initial_Statistics() const;
    void Update_Body(const ID index);
    void Insert_Into_Hashtable(const ID index);
    bool Remove_From_Cell(const TV_INT& voxel,const ID index);
    void Add_To_Cell(const TV_INT& voxel,const ID index);
    void Remove_If_Not_Still_Present(const RANGE<TV_INT>& old_locations,const RANGE<TV_INT>& new_locations,const ID index);
    void Add_If_Newly_Present(const RANGE<TV_INT>& old_locations,const RANGE<TV_INT>& new_locations,const ID index);
//#####################################################################
};
}
#endif

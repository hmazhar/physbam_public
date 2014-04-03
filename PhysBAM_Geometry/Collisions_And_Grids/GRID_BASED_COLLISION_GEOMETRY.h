//#####################################################################
// Copyright 2005-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_BASED_COLLISION_GEOMETRY
//#####################################################################
#ifndef __GRID_BASED_COLLISION_GEOMETRY__
#define __GRID_BASED_COLLISION_GEOMETRY__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/OBJECTS_IN_CELL.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic_Arrays/GRID_ARRAYS_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE_Arrays/GRID_ARRAYS_POLICY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/LINEAR_INTERPOLATION_RLE.h>
#endif
namespace PhysBAM{

template <class T_GRID>
class GRID_BASED_COLLISION_GEOMETRY:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::INDEX T_INDEX;typedef typename T_GRID::BLOCK T_BLOCK;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_ARRAYS_BOOL::template REBIND<VECTOR<bool,T_GRID::dimension> >::TYPE T_ARRAYS_BOOL_DIMENSION;
    typedef typename T_FACE_ARRAYS_BOOL::template REBIND<VECTOR<bool,T_GRID::dimension> >::TYPE T_FACE_ARRAYS_BOOL_DIMENSION;
    typedef typename IF<TV::dimension==2,T,typename IF<TV::dimension==1,ONE,TV>::TYPE>::TYPE T_WEIGHTS;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX T_SIMPLEX;
public:
    COLLISION_GEOMETRY_COLLECTION<TV> collision_geometry_collection;
    // TODO: Add flag to disable collision bodies for fluid.

    T_GRID& grid;
    T collision_thickness;
    int number_of_ghost_cells;

    OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID> objects_in_cell;
    T_ARRAYS_BOOL occupied_blocks;
    T_ARRAYS_BOOL swept_occupied_blocks;
    T_ARRAYS_BOOL_DIMENSION cell_neighbors_visible; // length T_GRID::dimension, order: right top back (for dyadic, tree is fully refined where this would have any effect)
    T_FACE_ARRAYS_BOOL_DIMENSION face_neighbors_visible; // length is T_GRID::dimension, order: right top back (for dyadic, tree is fully refined where this would have any effect)

    GRID_BASED_COLLISION_GEOMETRY(T_GRID& grid_input);
    virtual ~GRID_BASED_COLLISION_GEOMETRY();

    void Save_State(const int state_index,const T time=0)
    {for(COLLISION_GEOMETRY_ID i(1);i<=collision_geometry_collection.bodies.m;i++) if(Is_Active(i)) collision_geometry_collection.bodies(i)->Save_State(state_index,time);}

    void Restore_State(const int state_index)
    {for(COLLISION_GEOMETRY_ID i(1);i<=collision_geometry_collection.bodies.m;i++) if(Is_Active(i)) collision_geometry_collection.bodies(i)->Restore_State(state_index);}

    void Average_States(const int state1, const int state2,const int result_state,const T interpolation_distance)
    {for(COLLISION_GEOMETRY_ID i(1);i<=collision_geometry_collection.bodies.m;i++) if(Is_Active(i)) collision_geometry_collection.bodies(i)->Average_States(state1,state2,result_state,interpolation_distance);}

    void Delete_State(const int state_index)
    {for(COLLISION_GEOMETRY_ID i(1);i<=collision_geometry_collection.bodies.m;i++) if(Is_Active(i)) collision_geometry_collection.bodies(i)->Delete_State(state_index);}
    
    bool Intersection_With_Any_Simplicial_Object(RAY<TV>& ray,COLLISION_GEOMETRY_ID& body_id,const ARRAY<COLLISION_GEOMETRY_ID>* objects=0) const
    {return collision_geometry_collection.Intersection_With_Any_Simplicial_Object(ray,body_id,objects);}
    
    void Add_Bodies(COLLISION_GEOMETRY_COLLECTION<TV>& collision_geometry_list)
    {collision_geometry_collection.Add_Bodies(collision_geometry_list);}

    void Remove_Body(COLLISION_GEOMETRY_ID id)
    {collision_geometry_collection.Remove_Body(id);}

    bool Is_Active(COLLISION_GEOMETRY_ID id) const
    {return collision_geometry_collection.Is_Active(id);}

//#####################################################################
    void Add_Bodies(RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection);
    virtual void Rasterize_Objects();
    bool Earliest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,T_WEIGHTS& weights,COLLISION_GEOMETRY_ID& body_id,int& simplex_id,
        const ARRAY<COLLISION_GEOMETRY_ID>* objects=0) const;
    bool Latest_Crossover(const TV& start_X,const TV& end_X,const T dt,COLLISION_GEOMETRY_ID& body_id,int& simplex_id,TV& initial_hit_point,const ARRAY<COLLISION_GEOMETRY_ID>* objects=0) const;
    bool Latest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,T_WEIGHTS& weights,COLLISION_GEOMETRY_ID& body_id,int& simplex_id,
        POINT_SIMPLEX_COLLISION_TYPE& returned_collision_type,const ARRAY<COLLISION_GEOMETRY_ID>* objects=0) const;
    bool Any_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,const ARRAY<COLLISION_GEOMETRY_ID>* objects=0) const;
    void Update_Intersection_Acceleration_Structures(const bool use_swept_simplex_hierarchy,const int state1=0,const int state2=0);
    bool Get_Body_Penetration(const TV& start_X,const TV& end_X,const T contour_value,const T dt,COLLISION_GEOMETRY_ID& body_id,int& simplex_id,
        T& start_phi,T& end_phi,TV& end_body_normal,TV& body_velocity,const ARRAY<COLLISION_GEOMETRY_ID>* objects=0) const;
    bool Push_Out_Point(TV& X,const T collision_distance,const bool check_particle_crossover,bool& particle_crossover,const ARRAY<COLLISION_GEOMETRY_ID>* objects=0) const;
    bool Occupied_Block(const T_BLOCK& block) const;
    bool Swept_Occupied_Block(const T_BLOCK& block) const;
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT    
    void Read_State(TYPED_ISTREAM& input,const int state_index);
    void Write_State(TYPED_OSTREAM& output,const int state_index) const;
#endif    
//#####################################################################
};
}
#endif

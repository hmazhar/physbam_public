//#####################################################################
// Copyright 2002-2008, Zhaosheng Bao, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_COLLISION_GEOMETRY_3D
//##################################################################### 
#ifndef __RIGID_COLLISION_GEOMETRY_3D__
#define __RIGID_COLLISION_GEOMETRY_3D__

#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY.h>
namespace PhysBAM{

template<class T> class TRIANGLE_HIERARCHY;
template<class TV> class RIGID_COLLISION_GEOMETRY;

template<class T_input>
class RIGID_COLLISION_GEOMETRY<VECTOR<T_input,3> >:public RIGID_COLLISION_GEOMETRY_BASE<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
    typedef RIGID_COLLISION_GEOMETRY_BASE<TV> BASE;
    typedef typename TV::SPIN T_SPIN;
public:
    using BASE::saved_states;using BASE::rigid_geometry;using BASE::collision_thickness;

    RIGID_COLLISION_GEOMETRY(RIGID_GEOMETRY<TV>& rigid_geometry_input);
public:
    virtual ~RIGID_COLLISION_GEOMETRY();

//#####################################################################
    bool Earliest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,TV& weights,int& triangle_id) const PHYSBAM_OVERRIDE;
    bool Latest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,TV& weights,int& triangle_id,POINT_SIMPLEX_COLLISION_TYPE& returned_collision_type) const PHYSBAM_OVERRIDE;
    bool Any_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt) const PHYSBAM_OVERRIDE;
    void Get_Simplex_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,const bool with_body_motion,const T extra_thickness,const T body_thickness_factor) const PHYSBAM_OVERRIDE;
    void Update_Intersection_Acceleration_Structures(const bool use_swept_triangle_hierarchy,const int state1=0,const int state2=0) PHYSBAM_OVERRIDE;
    TRIANGLE_3D<T> World_Space_Simplex(const int triangle_id,const bool use_saved_state=false) const PHYSBAM_OVERRIDE;
    TRIANGLE_3D<T> World_Space_Simplex(const int triangle_id,const FRAME<TV>& frame) const;
//#####################################################################
};
}
#endif

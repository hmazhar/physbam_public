//#####################################################################
// Copyright 2004, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_RAY
//#####################################################################
#ifndef __RENDERING_RAY__
#define __RENDERING_RAY__

#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
namespace PhysBAM{

template<class T> class RENDERING_RAY_DEBUG;
template<class T> class RENDERING_OBJECT;

template<class T>
class RENDERING_RAY
{
public:
    RAY<VECTOR<T,3> > ray;
    int recursion_depth;
    T ray_contribution;
    RENDERING_RAY_DEBUG<T>* debug_ray;
    const RENDERING_OBJECT<T>* current_object;
    enum {DUMMY_RAY,COLOR_RAY,SHADOW_RAY,PHOTON_RAY,PHOTON_GATHER_RAY,UNKNOWN_RAY} ray_type;

    RENDERING_RAY(const RAY<VECTOR<T,3> >& ray,const T ray_contribution_input,const RENDERING_OBJECT<T>* current_object_input)
        :ray(ray),recursion_depth(0),ray_contribution(ray_contribution_input),debug_ray(0),current_object(current_object_input),ray_type(UNKNOWN_RAY)
    {}

    RENDERING_RAY()
        :recursion_depth(0),ray_contribution(1),debug_ray(0),current_object(0),ray_type(UNKNOWN_RAY)
    {}

//#####################################################################
};
}
#endif

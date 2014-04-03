//#####################################################################
// Copyright 2004, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_RAY_DEBUG
//#####################################################################
#ifndef __RENDERING_RAY_DEBUG__
#define __RENDERING_RAY_DEBUG__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY.h>
namespace PhysBAM{

template<class T> class PHOTON;
template<class T> class IRRADIANCE_SAMPLE;

template<class T>
class RENDERING_RAY_DEBUG:public NONCOPYABLE
{
public:
    RENDERING_RAY<T> ray;
    RENDERING_RAY_DEBUG* parent;
    ARRAY<RENDERING_RAY_DEBUG*> children;
    ARRAY<std::string> comments;
    ARRAY<PHOTON<T>*> photons_used;
    ARRAY<int> irradiance_cache_samples_used;
    bool hit_object;
    VECTOR<T,3> same_side_normal;

    RENDERING_RAY_DEBUG()
        :parent(0),hit_object(false)
    {}

    RENDERING_RAY_DEBUG(const RENDERING_RAY<T>& ray_input)
        :ray(ray_input),parent(0),hit_object(false)
    {}

    ~RENDERING_RAY_DEBUG()
    {for(int i=1;i<=children.m;i++)delete children(i);}

    void Add_Child(RENDERING_RAY<T>& ray_to_add)
    {RENDERING_RAY_DEBUG* child=new RENDERING_RAY_DEBUG(ray_to_add);
    child->parent=this;
    ray_to_add.debug_ray=child->ray.debug_ray=child;
    children.Append(child);}

    void Add_Comment(const std::string& comment_string)
    {comments.Append(comment_string);}

    void Print(std::ostream& output,const int spaces=0)
    {for(int space=1;space<=spaces;space++)output<<" ";
    std::string ray_type_str;
    switch(ray.ray_type){
        case RENDERING_RAY<T>::DUMMY_RAY:ray_type_str="DUMMY";break;
        case RENDERING_RAY<T>::COLOR_RAY:ray_type_str="COLOR";break;
        case RENDERING_RAY<T>::SHADOW_RAY:ray_type_str="SHADOW";break;
        case RENDERING_RAY<T>::UNKNOWN_RAY:ray_type_str="UNKNOWN";break;}
    output<<"RAY "<<ray_type_str<<" endpoint: "<<ray.ray.endpoint<<" direction: "<<ray.ray.direction<<std::endl;
    for(int i=1;i<=children.m;i++)children(i)->Print(spaces+5);}

//#####################################################################
};
}
#endif

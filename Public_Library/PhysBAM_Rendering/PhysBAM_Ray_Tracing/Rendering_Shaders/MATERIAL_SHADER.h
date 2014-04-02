//#####################################################################
// Copyright 2003-2004, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATERIAL_SHADER
//#####################################################################
#ifndef __MATERIAL_SHADER__
#define __MATERIAL_SHADER__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class T,int d> class VECTOR;
template<class T> class RENDER_WORLD;
template<class T> class RENDERING_RAY;
template<class T> class RENDERING_OBJECT;
template<class T> class RENDERING_LIGHT;
template<class T> class PHOTON_MAP;

template<class T>
class MATERIAL_SHADER:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;
protected:
    RENDER_WORLD<T>& world;
public:

    MATERIAL_SHADER(RENDER_WORLD<T>& world_input)
        :world(world_input)
    {}

    virtual ~MATERIAL_SHADER()
    {}

//#####################################################################
    virtual TV Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal) const=0;
    virtual TV Shade_Surface_Using_Indirect_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal) const {return TV();}
    virtual TV Shade_Light_Ray(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal,const RENDERING_LIGHT<T>& light,
        const RENDERING_RAY<T>& full_ray) const {return TV();}
    virtual TV Shade_Surface_Using_Approximate_Full_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal) const {return TV();}
    virtual void Receive_Photon(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& intersection_point,const TV& same_side_normal,const TV& power,
        const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,const int diffuse_bounces,const int specular_bounces) const {}
//#####################################################################
};
}
#endif

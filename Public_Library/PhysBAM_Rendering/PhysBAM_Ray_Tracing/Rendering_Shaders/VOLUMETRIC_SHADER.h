//#####################################################################
// Copyright 2003-2005, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VOLUMETRIC_SHADER
//#####################################################################
#ifndef __VOLUMETRIC_SHADER__
#define __VOLUMETRIC_SHADER__

#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/PHOTON_MAP.h>
namespace PhysBAM{

template<class T> class RENDER_WORLD;
template<class T> class RENDERING_LIGHT;

template<class T>
class VOLUMETRIC_SHADER:public NONCOPYABLE
{
public:
    bool supports_photon_mapping;
    RENDER_WORLD<T>& world;

    VOLUMETRIC_SHADER(RENDER_WORLD<T>& world_input)
        :supports_photon_mapping(false),world(world_input)
    {}

    virtual ~VOLUMETRIC_SHADER()
    {}

    static T Isotropic_Phase_Function(const VECTOR<T,3>& incoming,const VECTOR<T,3>& outgoing)
    {return T(1)/T(4*pi);}

    static T Henyey_Greenstein_Phase_Function(const VECTOR<T,3>& incoming,const VECTOR<T,3>& outgoing,const T g)
    {T cos_theta=VECTOR<T,3>::Dot_Product(incoming,outgoing);
    return ((T)1-sqr(g))/(T(4*pi)*pow((T)1+sqr(g)-(T)2*g*cos_theta,(T)1.5));}

    VECTOR<T,3> Inverse_Isotropic_Phase_Direction(const VECTOR<T,3>& incoming)
    {return world.random.template Get_Direction<VECTOR<T,3> >();}

    VECTOR<T,3> Inverse_Henyey_Greenstein_Phase_Direction(const VECTOR<T,3>& incoming,const T g)
    {T u=world.random.Get_Uniform_Number(T(0),T(1));
        T cos_theta=(1/(2*g))*(1+sqr(g)-sqr((1-sqr(g))/(1-g+2*g*u)));
    VECTOR<T,3> random_direction=world.random.template Get_Direction<VECTOR<T,3> >();
    VECTOR<T,3> orthogonal_direction=(random_direction-random_direction).Normalized();
    T orthogonal_length=tan(acos(cos_theta));
    return (incoming+orthogonal_direction*orthogonal_length).Normalized();}

//#####################################################################
    virtual VECTOR<T,3> Attenuate_Color(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& object,const VECTOR<T,3>& color)=0;
    virtual VECTOR<T,3> Attenuate_Light(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& object,const RENDERING_LIGHT<T>& light,const VECTOR<T,3>& light_color)=0;
    virtual VECTOR<T,3> Attenuate_Photon(const RENDERING_RAY<T>& ray, const RENDERING_OBJECT<T>& object, const VECTOR<T,3>& photon_power, bool& should_throw){return VECTOR<T,3>(0,0,0);}
    virtual bool Scatter_Photon_Ray(const RENDERING_OBJECT<T>& object,RENDERING_RAY<T>& modified_ray,VECTOR<T,3>& photon_power,const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,const T fixed_step_size){return false;}
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_POINT_LIGHT__
#define __RENDERING_POINT_LIGHT__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Lights/RENDERING_LIGHT.h>
namespace PhysBAM{

template<class T>
class RENDERING_POINT_LIGHT:public RENDERING_LIGHT<T>
{
public:
    using RENDERING_LIGHT<T>::world;using RENDERING_LIGHT<T>::position;using RENDERING_LIGHT<T>::color;using RENDERING_LIGHT<T>::brightness;using RENDERING_LIGHT<T>::global_photon_random;
    using RENDERING_LIGHT<T>::caustic_photon_random;using RENDERING_LIGHT<T>::volume_photon_random;using RENDERING_LIGHT<T>::sample_points_random;
    using RENDERING_LIGHT<T>::supports_global_photon_mapping;using RENDERING_LIGHT<T>::supports_caustic_photon_mapping;using RENDERING_LIGHT<T>::supports_volume_photon_mapping;

    bool distance_attenuation;

    RENDERING_POINT_LIGHT(const VECTOR<T,3>& position_input,const VECTOR<T,3>& color_input,const T brightness_input,bool distance_attenuation_input,RENDER_WORLD<T>& world_input,
        const bool supports_global_photons,const bool supports_caustic_photons,const bool supports_volume_photons)
        :RENDERING_LIGHT<T>(position_input,color_input,brightness_input,world_input,supports_global_photons,supports_caustic_photons,supports_volume_photons),
        distance_attenuation(distance_attenuation_input)
    {}

    void Sample_Points(const VECTOR<T,3>& surface_position,const VECTOR<T,3>& surface_normal,ARRAY<RAY<VECTOR<T,3> > >& sample_array)const PHYSBAM_OVERRIDE
    {sample_array.Resize(1);sample_array(1)=RAY<VECTOR<T,3> >(SEGMENT_3D<T>(surface_position,position));}

    VECTOR<T,3> Emitted_Light(const RENDERING_RAY<T>& ray) const PHYSBAM_OVERRIDE
    {if(distance_attenuation)return brightness*color*T(one_over_four_pi)/sqr(ray.ray.t_max);
    else return brightness*color;}

    int Emit_Photons(RENDERING_RAY<T>& parent_ray,PHOTON_MAP<T>& photon_map,const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type) const  PHYSBAM_OVERRIDE
    {int photons_emitted=0;
    while(photon_map.Light_Emission_Quota_Remains()){
        RANDOM_NUMBERS<T>* random=0;
        if(type==PHOTON_MAP<T>::GLOBAL_PHOTON_MAP) random=&global_photon_random;
        else if(type==PHOTON_MAP<T>::CAUSTIC_PHOTON_MAP) random=&caustic_photon_random;
        else if(type==PHOTON_MAP<T>::VOLUME_PHOTON_MAP) random=&volume_photon_random;
        world.random.Set_Seed((int)abs(random->Get_Uniform_Number((T)0,(T)1000000)));
        photons_emitted++;
        VECTOR<T,3> direction(world.random.template Get_Direction<VECTOR<T,3> >());
        RENDERING_RAY<T> photon_ray(RAY<VECTOR<T,3> >(position,direction),1,parent_ray.current_object);
        world.Cast_Photon(photon_ray,parent_ray,color*brightness,type,0,0);}
    return photons_emitted;}

//#####################################################################
};
}
#endif

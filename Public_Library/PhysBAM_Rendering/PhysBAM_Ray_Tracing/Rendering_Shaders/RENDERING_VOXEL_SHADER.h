//#####################################################################
// Copyright 2003-2005, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_VOXEL_SHADER__
#define __RENDERING_VOXEL_SHADER__

#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/BLACKBODY.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_LEVELSET_MULTIPLE_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/VOLUMETRIC_SHADER.h>
namespace PhysBAM{

template<class TV> class IMPLICIT_OBJECT;

template<class T>
class RENDERING_VOXEL_SHADER:public VOLUMETRIC_SHADER<T>
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    using VOLUMETRIC_SHADER<T>::world;using VOLUMETRIC_SHADER<T>::supports_photon_mapping;

    TV absorption;
    TV absorption_shadow;
    TV scattering;
    T inscattering_amplification;
    bool isotropic_scattering;
    T phase_function_g;
    BLACKBODY<T> blackbody;
    MATRIX<T,3> world_xyz_to_display_xyz;
    T emission_amplification;
    bool use_lms_scaling;
    const IMPLICIT_OBJECT<TV>* empty_implicit_surface;
    const RENDERING_LEVELSET_MULTIPLE_OBJECT<LEVELSET_MULTIPLE<GRID<TV> > >* levelset_multiple_object;
    int non_empty_region;
    bool use_empty_implicit_surface_for_light_attenuation;
    
    bool use_temperature_remap;
    bool use_blackbody_ramp;
    bool use_constant_emission_color;
    TV constant_emission_color;
    ARRAY<T,VECTOR<int,1> > temperature_remap;
    GRID<VECTOR<T,1> > temperature_remap_grid;
    LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,T> temperature_remap_interpolation;

    GRID<VECTOR<T,1> > blackbody_ramp_grid;
    ARRAY<TV,VECTOR<int,1> > blackbody_ramp;
    LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,TV> blackbody_interpolation;

    RENDERING_VOXEL_SHADER(const TV& absorption_input,const TV& scattering_input,const T inscattering_amplification_input,
        const T emission_amplification_input,const T white_point_temperature,const bool use_lms_scaling_input,
        const bool use_blackbody_ramp_input,const bool use_constant_emission_color_input,const TV constant_emission_color_input,
        RENDER_WORLD<T>& world_input);
    virtual ~RENDERING_VOXEL_SHADER();

    T Phase(const TV& incoming, const TV& outgoing)
    {if(isotropic_scattering)return VOLUMETRIC_SHADER<T>::Isotropic_Phase_Function(incoming,outgoing);
    else return VOLUMETRIC_SHADER<T>::Henyey_Greenstein_Phase_Function(incoming,outgoing,phase_function_g);}

    TV Inverse_Phase_Direction(const TV& incoming)
    {if(isotropic_scattering) return Inverse_Isotropic_Phase_Direction(incoming);
    else return Inverse_Henyey_Greenstein_Phase_Direction(incoming,phase_function_g);}
    
    void Use_Empty_Implicit_Surface(const IMPLICIT_OBJECT<TV>& empty_implicit_surface_input,const bool use_empty_implicit_surface_for_light_attenuation_input=true)
    {empty_implicit_surface=&empty_implicit_surface_input;supports_photon_mapping=false;use_empty_implicit_surface_for_light_attenuation=use_empty_implicit_surface_for_light_attenuation_input;}

    void Use_Levelset_Multiple_Region(const RENDERING_LEVELSET_MULTIPLE_OBJECT<LEVELSET_MULTIPLE<GRID<TV> > >* levelset_multiple_object_input,int region_in){
      levelset_multiple_object=levelset_multiple_object_input;
      non_empty_region=region_in;}

    void Use_Temperature_Remap(ARRAY<T,VECTOR<int,1> >& temperature_remap_input)
    {temperature_remap=temperature_remap_input;use_temperature_remap=true;
    temperature_remap_grid.Initialize(temperature_remap.counts.x,0,3000);}

    void Set_Absorption_Shadow(const TV& absorption_shadow_input)
    {absorption_shadow=absorption_shadow_input;}

//#####################################################################
    TV Attenuate_Color(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& object,const TV& color) PHYSBAM_OVERRIDE;
    TV Attenuate_Light(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& object,const RENDERING_LIGHT<T>& light,const TV& light_color) PHYSBAM_OVERRIDE;
    bool Scatter_Photon_Ray(const RENDERING_OBJECT<T>& object,RENDERING_RAY<T>& ray,TV& photon_power,const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type,const T fixed_step_size) PHYSBAM_OVERRIDE;
    TV Attenuate_Photon(const RENDERING_RAY<T>& ray, const RENDERING_OBJECT<T>& object, const TV& photon_power, bool& should_throw) PHYSBAM_OVERRIDE {return TV(photon_power);}
//#####################################################################
};
}
#endif

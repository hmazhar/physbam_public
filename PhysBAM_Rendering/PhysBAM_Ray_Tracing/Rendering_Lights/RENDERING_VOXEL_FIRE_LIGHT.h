//#####################################################################
// Copyright 2004-2007, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_VOXEL_FIRE_LIGHT
//#####################################################################
#ifndef __RENDERING_VOXEL_FIRE_LIGHT__
#define __RENDERING_VOXEL_FIRE_LIGHT__

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Log/PROGRESS_INDICATOR.h>
#include <PhysBAM_Tools/Random_Numbers/PIECEWISE_CONSTANT_PDF.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Lights/RENDERING_LIGHT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_UNIFORM_VOXELS.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/RENDERING_VOXEL_SHADER.h>
namespace PhysBAM{

template<class T>
class RENDERING_VOXEL_FIRE_LIGHT:public RENDERING_LIGHT<T>
{
    typedef VECTOR<T,3> TV;
private:
    using RENDERING_LIGHT<T>::world;using RENDERING_LIGHT<T>::supports_global_photon_mapping;using RENDERING_LIGHT<T>::supports_caustic_photon_mapping;
    using RENDERING_LIGHT<T>::supports_volume_photon_mapping;using RENDERING_LIGHT<T>::global_photon_random;using RENDERING_LIGHT<T>::caustic_photon_random;
    using RENDERING_LIGHT<T>::volume_photon_random;using RENDERING_LIGHT<T>::sample_points_random;

    RENDERING_UNIFORM_VOXELS<T>& fire_voxels;
    RENDERING_VOXEL_SHADER<T>& fire_shader;
    T total_average_power;
    VECTOR<T,3> total_power;
    ARRAY<T,VECTOR<int,3> > pdf;
    ARRAY<PIECEWISE_CONSTANT_PDF<T> ,VECTOR<int,2> > z_cdf;
    ARRAY<PIECEWISE_CONSTANT_PDF<T> ,VECTOR<int,1> > y_cdf;
    PIECEWISE_CONSTANT_PDF<T> x_cdf;
public:
    RENDERING_VOXEL_FIRE_LIGHT(RENDERING_UNIFORM_VOXELS<T>& fire_voxels_input,RENDERING_VOXEL_SHADER<T>& fire_shader_input,RENDER_WORLD<T>& world_input,const bool supports_global_photons,
        const bool supports_caustic_photons,const bool supports_volume_photons)
        :RENDERING_LIGHT<T>(VECTOR<T,3>(0,0,0),VECTOR<T,3>(0,0,0),0,world_input,supports_global_photons,supports_caustic_photons,supports_volume_photons),
        fire_voxels(fire_voxels_input),fire_shader(fire_shader_input)
    {
        // build pdf
        LOG::SCOPE scope("building fire light pdf","Building Fire_Light PDF");
        GRID<TV>& grid=fire_voxels.grid;
        pdf.Resize(1,fire_voxels.grid.counts.x,1,fire_voxels.grid.counts.y,1,fire_voxels.grid.counts.z);
        T cell_volume=fire_voxels.grid.dX.x*fire_voxels.grid.dX.y*fire_voxels.grid.dX.z;
        total_average_power=0;total_power=VECTOR<T,3>(0,0,0);
        for(int i=1;i<=fire_voxels.grid.counts.x;i++)for(int j=1;j<=fire_voxels.grid.counts.y;j++)for(int ij=1;ij<=fire_voxels.grid.counts.z;ij++){
            T cell_power;VECTOR<T,3> cell_power_spectrum;
            if(fire_shader.empty_implicit_surface && fire_shader.empty_implicit_surface->Lazy_Inside(fire_voxels.grid.X(i,j,ij))){cell_power=0;cell_power_spectrum=VECTOR<T,3>(0,0,0);}
            else{
                VECTOR<T,3> xyz=fire_shader.world_xyz_to_display_xyz*fire_shader.blackbody.Calculate_XYZ(fire_voxels.Source_Term(2,fire_voxels.grid.X(i,j,ij)));
                cell_power_spectrum=fire_shader.blackbody.cie.XYZ_To_RGB(xyz)*fire_shader.emission_amplification;
                T luminance=xyz.y*fire_shader.emission_amplification;
                cell_power=T(4*T(pi)*cell_volume)*luminance;
                cell_power_spectrum*=T(4*T(pi)*cell_volume);}
            total_power+=cell_power_spectrum;
            total_average_power+=cell_power;pdf(i,j,ij)=cell_power;}
        // generate cdf for Pr(z|x,y)
        PROGRESS_INDICATOR progress(grid.counts.x);
        LOG::cout<<"    Generating CDF for Pr(z|x,y)";
        z_cdf.Resize(1,grid.counts.x,1,grid.counts.y);
        for(int i=1;i<=grid.counts.x;i++){
            progress.Progress();
            for(int j=1;j<=grid.counts.y;j++){
                z_cdf(i,j).Initialize(grid.counts.z);
                for(int ij=1;ij<=grid.counts.z;ij++)z_cdf(i,j).pdf(ij)=pdf(i,j,ij);
                z_cdf(i,j).Compute_Cumulative_Distribution_Function();}}
        // generate cdf for Pr(y|x)
        progress.Initialize(grid.counts.x);
        LOG::cout<<"    Generating CDF for Pr(y|x)";
        y_cdf.Resize(1,grid.counts.x);
        for(int i=1;i<=grid.counts.x;i++){
            progress.Progress();
            y_cdf(i).Initialize(grid.counts.y);
            for(int j=1;j<=grid.counts.y;j++)y_cdf(i).pdf(j)=z_cdf(i,j).normalization_constant;
            y_cdf(i).Compute_Cumulative_Distribution_Function();}
        // generate cdf for Pr(x)
        LOG::cout<<"    Generating CDF for Pr(x)..."<<std::endl;
        x_cdf.Initialize(grid.counts.x);
        for(int i=1;i<=grid.counts.x;i++)x_cdf.pdf(i)=y_cdf(i).normalization_constant;
        x_cdf.Compute_Cumulative_Distribution_Function();
        LOG::cout<<"x_cdf normalization "<<x_cdf.normalization_constant<<" total_power "<<total_average_power<<std::endl;
    }

    virtual ~RENDERING_VOXEL_FIRE_LIGHT()
    {}

    void Sample_Points(const VECTOR<T,3>& surface_position,const VECTOR<T,3>& surface_normal,ARRAY<RAY<VECTOR<T,3> > >& sample_array)const PHYSBAM_OVERRIDE
    {sample_array.Resize(1);sample_array(1)=RAY<VECTOR<T,3> >(surface_position,VECTOR<T,3>(1,0,0),true);} // just give dummy return to avoid division by zero

    VECTOR<T,3> Sample_Point_In_Volume(const T xi_x,const T xi_y,const T xi_z,T& probability_of_location)const
    {PAIR<int,T> sample_x=x_cdf.Sample(xi_x);int i=sample_x.x;
    PAIR<int,T> sample_y=y_cdf(i).Sample(xi_y);int j=sample_y.x;
    PAIR<int,T> sample_z=z_cdf(i,j).Sample(xi_z);int ij=sample_z.x;
    probability_of_location=pdf(i,j,ij)/total_average_power;
    return VECTOR<T,3>(sample_x.y,sample_y.y,sample_z.y)-VECTOR<T,3>((T).5,(T).5,(T).5)+fire_voxels.grid.X(i,j,ij);}

    VECTOR<T,3> Emitted_Light(const RENDERING_RAY<T>& ray) const PHYSBAM_OVERRIDE
    {return VECTOR<T,3>(0,0,0);}

    int Emit_Photons(RENDERING_RAY<T>& parent_ray,PHOTON_MAP<T>& photon_map,const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type) const  PHYSBAM_OVERRIDE
    {int number_emitted=0;
    while(photon_map.Light_Emission_Quota_Remains()){
        RANDOM_NUMBERS<T>* random=0;
        if(type==PHOTON_MAP<T>::GLOBAL_PHOTON_MAP) random=&global_photon_random;
        else if(type==PHOTON_MAP<T>::CAUSTIC_PHOTON_MAP) random=&caustic_photon_random;
        else if(type==PHOTON_MAP<T>::VOLUME_PHOTON_MAP) random=&volume_photon_random;
        world.random.Set_Seed(abs((int)random->Get_Uniform_Number((T)0,(T)1)));

        T xi_x=world.random.Get_Uniform_Number((T)0,(T)1),xi_y=world.random.Get_Uniform_Number((T)0,(T)1),xi_z=world.random.Get_Uniform_Number((T)0,(T)1);
        T probability;
        VECTOR<T,3> sample_location=Sample_Point_In_Volume(xi_x,xi_y,xi_z,probability);VECTOR<T,3> sample_direction(world.random.template Get_Direction<VECTOR<T,3> >());
        RENDERING_RAY<T> photon_emit_ray(RAY<VECTOR<T,3> >(sample_location,sample_direction,true),1,parent_ray.current_object);
        world.Cast_Photon(photon_emit_ray,parent_ray,total_power,type,0,0);
        number_emitted++;}
    return number_emitted;}

//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2002-2009, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_UNIFORM_VOXELS
//#####################################################################
#ifndef __RENDERING_UNIFORM_VOXELS__
#define __RENDERING_UNIFORM_VOXELS__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_VOXELS.h>
namespace PhysBAM{

template<class T>
class RENDERING_UNIFORM_VOXELS:public RENDERING_VOXELS<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_VOXELS<T>::box;using RENDERING_VOXELS<T>::small_number;using RENDERING_VOXELS<T>::precompute_single_scattering;

    GRID<TV>& grid;
    GRID<TV>& coarse_grid;
    ARRAY<ARRAY<T,VECTOR<int,3> >*> data; // defined at each grid point
    ARRAY<T> data_scale;
    ARRAY<T> data_offset;
    ARRAY<bool> data_clamp_low_value;
    ARRAY<T> data_lowest_value;
    ARRAY<bool> data_clamp_high_value;
    ARRAY<T> data_highest_value;
    ARRAY<ARRAY<VECTOR<T,3> ,VECTOR<int,3> >*> precomputed_light;
    ARRAY<ARRAY<bool,VECTOR<int,3> >*> precomputed_light_valid;
    ARRAY<VECTOR<int,3> > map_from_accessor_index_to_my_index;
    T volumetric_step;
    INTERPOLATION_UNIFORM<GRID<TV>,T>* voxel_source_interpolation;
    INTERPOLATION_UNIFORM<GRID<TV>,VECTOR<T,3> >* voxel_light_interpolation;
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> default_voxel_source_interpolation;
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,VECTOR<T,3> > default_voxel_light_interpolation;
    int number_of_smoothing_steps;

    RENDERING_UNIFORM_VOXELS(GRID<TV>& grid_input,ARRAY<T,VECTOR<int,3> >& data_input,const T volumetric_step)
        :grid(grid_input),coarse_grid(grid_input),volumetric_step(volumetric_step),number_of_smoothing_steps(0)
    {
        box=grid_input.domain;
        voxel_light_interpolation=&default_voxel_light_interpolation;voxel_source_interpolation=&default_voxel_source_interpolation;
        data.Append(&data_input);
    }

    RENDERING_UNIFORM_VOXELS(GRID<TV>& grid_input,GRID<TV>& coarse_grid_input,ARRAY<T,VECTOR<int,3> >& data_input,const T volumetric_step)
        :grid(grid_input),coarse_grid(coarse_grid_input),volumetric_step(volumetric_step),number_of_smoothing_steps(0)
    {
        box=grid_input.domain;
        voxel_light_interpolation=&default_voxel_light_interpolation;voxel_source_interpolation=&default_voxel_source_interpolation;
        data.Append(&data_input);
    }

    void Print_Parameters()
    {for(int i=1;i<=data.Size();i++){
        LOG::cout<<"data_scale("<<i<<")="<<data_scale(i)<<" data_offset("<<i<<")="<<data_offset(i)<<std::endl;
        LOG::cout<<"data_clamp_low_value("<<i<<")="<<data_clamp_low_value(i)<<" data_lowest_value("<<i<<")="<<data_lowest_value(i)<<std::endl;
        LOG::cout<<"data_clamp_high_value("<<i<<")="<<data_clamp_high_value(i)<<" data_highest_value("<<i<<")="<<data_highest_value(i)<<std::endl;
        LOG::cout<<std::endl;}}

    T Volumetric_Integration_Step(const RAY<VECTOR<T,3> > &ray,const T xi) const PHYSBAM_OVERRIDE
    {return xi*volumetric_step;}

    T Source_Term(const int source_term_index,const VECTOR<T,3>& location) const PHYSBAM_OVERRIDE
    {T value=voxel_source_interpolation->Clamped_To_Array(grid,*data(source_term_index),location);
    value=data_scale(source_term_index)*(value+data_offset(source_term_index));
    if(data_clamp_low_value(source_term_index)) value=max(value,data_lowest_value(source_term_index));
    if(data_clamp_high_value(source_term_index)) value=min(value,data_highest_value(source_term_index));
    return value;}

    void Get_Node_Locations(ARRAY<VECTOR<T,3> >& locations) PHYSBAM_OVERRIDE
    {locations.Resize(coarse_grid.counts.Product());map_from_accessor_index_to_my_index.Resize(locations.m);int index=1;
    for(int i=1;i<=coarse_grid.counts.x;i++)for(int j=1;j<=coarse_grid.counts.y;j++)for(int ij=1;ij<=coarse_grid.counts.z;ij++){
        map_from_accessor_index_to_my_index(index)=VECTOR<int,3>(i,j,ij);
        locations(index)=coarse_grid.X(i,j,ij);index++;}}

    bool Use_Precomputed_Light_Data(const VECTOR<T,3>& location,const int light_index) const PHYSBAM_OVERRIDE
    {if(!precompute_single_scattering)return false;
    VECTOR<int,3> index=INTERPOLATION_UNIFORM<GRID<TV>,VECTOR<T,3> >::Clamped_Index_End_Minus_One(coarse_grid,*precomputed_light(light_index),location);
    int i=index.x,j=index.y,ij=index.z;
    if(coarse_grid.Outside(location))return false;
    bool i_j_ij=(*precomputed_light_valid(light_index))(i,j,ij),i_j_ij1=(*precomputed_light_valid(light_index))(i,j,ij+1),i_j1_ij=(*precomputed_light_valid(light_index))(i,j+1,ij),
        i_j1_ij1=(*precomputed_light_valid(light_index))(i,j+1,ij+1),i1_j_ij=(*precomputed_light_valid(light_index))(i+1,j,ij),i1_j_ij1=(*precomputed_light_valid(light_index))(i+1,j,ij+1),
        i1_j1_ij=(*precomputed_light_valid(light_index))(i+1,j+1,ij),i1_j1_ij1=(*precomputed_light_valid(light_index))(i+1,j+1,ij+1);
        return i_j_ij&&i_j_ij1&&i_j1_ij&&i_j1_ij1&&i1_j_ij&&i1_j_ij1&&i1_j1_ij&&i1_j1_ij1;}
    
    void Set_Precomputed_Light_Data(const int location_index,const int light_index,const VECTOR<T,3>& light_value) PHYSBAM_OVERRIDE
    {(*precomputed_light(light_index))(map_from_accessor_index_to_my_index(location_index))=light_value;}

    void Set_Precomputed_Light_Valid(const int location_index,const int light_index,const bool value) PHYSBAM_OVERRIDE
    {VECTOR<int,3> index=map_from_accessor_index_to_my_index(location_index);
    (*precomputed_light_valid(light_index))(index)=value;}

    VECTOR<T,3> Precomputed_Light_Data(const VECTOR<T,3>& location,const int light) const PHYSBAM_OVERRIDE
    {return voxel_light_interpolation->Clamped_To_Array(coarse_grid,*precomputed_light(light),location);}

    void Set_Custom_Source_Interpolation(INTERPOLATION_UNIFORM<GRID<TV>,T>* interpolation)
    {voxel_source_interpolation=interpolation;}

    void Set_Custom_Light_Interpolation(INTERPOLATION_UNIFORM<GRID<TV>,VECTOR<T,3> >* interpolation)
    {voxel_light_interpolation=interpolation;}

protected:
    void Prepare_For_Precomputation(RENDER_WORLD<T>& world) PHYSBAM_OVERRIDE
    {precomputed_light.Resize(world.Lights().m);precomputed_light_valid.Resize(world.Lights().m);
    for(int i=1;i<=precomputed_light.m;i++)precomputed_light(i)=new ARRAY<VECTOR<T,3> ,VECTOR<int,3> >(1,coarse_grid.counts.x,1,coarse_grid.counts.y,1,coarse_grid.counts.z);
    for(int i=1;i<=precomputed_light.m;i++){precomputed_light_valid(i)=new ARRAY<bool,VECTOR<int,3> >(1,coarse_grid.counts.x,1,coarse_grid.counts.y,1,coarse_grid.counts.z);precomputed_light_valid(i)->Fill(false);}}

//#####################################################################
    void Postprocess_Light_Field() PHYSBAM_OVERRIDE;
//#####################################################################
};   
}
#endif


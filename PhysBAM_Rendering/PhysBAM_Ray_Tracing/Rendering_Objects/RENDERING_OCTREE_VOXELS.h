#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Frank Losasso, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_OCTREE_VOXELS
//#####################################################################
#ifndef __RENDERING_OCTREE_VOXELS__
#define __RENDERING_OCTREE_VOXELS__

#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_PLANE_INTERSECTION.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_VOXELS.h>
namespace PhysBAM{

template<class T> class RENDERING_LIGHT;

template<class T>
class RENDERING_OCTREE_VOXELS:public RENDERING_VOXELS<T>
{
public:
    using RENDERING_VOXELS<T>::box;using RENDERING_VOXELS<T>::small_number;using RENDERING_VOXELS<T>::precompute_single_scattering;

    OCTREE_GRID<T>& octree_grid;
    ARRAY<ARRAY<T>*> data; // defined at each grid point
    mutable OCTREE_CELL<T>* current_cell;
    mutable OCTREE_CELL<T>* old_cell;
    ARRAY<ARRAY<bool>* > precomputed_light_valid;
    ARRAY<ARRAY<bool>* > precomputed_light_ever_valid;
    ARRAY<ARRAY<VECTOR<T,3> >* > precomputed_light;
    INTERPOLATION_DYADIC<OCTREE_GRID<T>,T>* voxel_source_interpolation;
    INTERPOLATION_DYADIC<OCTREE_GRID<T>,VECTOR<T,3> >* voxel_light_interpolation;
    LINEAR_INTERPOLATION_DYADIC<OCTREE_GRID<T>,T> default_voxel_source_interpolation;
    LINEAR_INTERPOLATION_DYADIC<OCTREE_GRID<T>,VECTOR<T,3> > default_voxel_light_interpolation;
    ARRAY<VECTOR<T,3> > node_locations;
    T volumetric_step;

    RENDERING_OCTREE_VOXELS(OCTREE_GRID<T>& octree_grid_input,ARRAY<T>& data_input,T step_input)
        :octree_grid(octree_grid_input),current_cell(0),old_cell(0)
    {
        box=octree_grid_input.Domain();
        voxel_light_interpolation=&default_voxel_light_interpolation;voxel_source_interpolation=&default_voxel_source_interpolation;
        data.Append(&data_input);
        volumetric_step=step_input;
    }

    ~RENDERING_OCTREE_VOXELS()
    {}

    void Get_Hits_With_Plane(const PLANE<T>& plane,RAY<VECTOR<T,3> > ray,T &min_hit,T &second_min_hit) const
    {if(INTERSECTION::Intersects(ray,plane) && ray.t_max<min_hit){second_min_hit=min_hit;min_hit=ray.t_max;}}

     T Volumetric_Integration_Step(const RAY<VECTOR<T,3> >& ray, const T xi) const PHYSBAM_OVERRIDE
    {VECTOR<T,3> original_location=ray.Point(0); //or some small number?
    if(!current_cell||octree_grid.Inside_Cell(current_cell,original_location)||current_cell->Has_Children())
        current_cell=octree_grid.Leaf_Cell(original_location);
    old_cell=current_cell;
    VECTOR<T,3> cell_size=current_cell->DX();
    T step_range=min(cell_size.x,cell_size.y,cell_size.z);
    VECTOR<T,3> center=current_cell->Center();VECTOR<T,3> DX_OVER_TWO=cell_size/2;
    T min_hit=FLT_MAX,second_min_hit=FLT_MAX;
    Get_Hits_With_Plane(PLANE<T>(VECTOR<T,3>(-1,0,0),center+VECTOR<T,3>(-DX_OVER_TWO.x,0,0)), ray, min_hit, second_min_hit);
    Get_Hits_With_Plane(PLANE<T>(VECTOR<T,3>(1,0,0),center+VECTOR<T,3>( DX_OVER_TWO.x,0,0)), ray, min_hit, second_min_hit);
    Get_Hits_With_Plane(PLANE<T>(VECTOR<T,3>(0,-1,0),center+VECTOR<T,3>(0,-DX_OVER_TWO.y,0)), ray, min_hit, second_min_hit);
    Get_Hits_With_Plane(PLANE<T>(VECTOR<T,3>(0,1,0),center+VECTOR<T,3>(0,DX_OVER_TWO.y,0)), ray, min_hit, second_min_hit);
    Get_Hits_With_Plane(PLANE<T>(VECTOR<T,3>(0,0,-1),center+VECTOR<T,3>(0,0,-DX_OVER_TWO.z)), ray, min_hit, second_min_hit);
    Get_Hits_With_Plane(PLANE<T>(VECTOR<T,3>(0,0,1),center+VECTOR<T,3>(0,0,DX_OVER_TWO.z)), ray, min_hit, second_min_hit);
    step_range=min(step_range,second_min_hit+small_number);
    //if(step_range>0.00048828126 && Octree_Cell_In_Thick_Smoke(*current_octree_cell)) step_range=(T)0.00048828125;
    T volumetric_integration_step=((T).2+xi*(T).6)*step_range;
    while(1){// get the new location and update the current_cell
        VECTOR<T,3> new_location=ray.Point(volumetric_integration_step);
        current_cell=octree_grid.Leaf_Cell_By_Neighbor_Path(current_cell,new_location);
        if(volumetric_integration_step<=current_cell->DX().x) break;
        volumetric_integration_step*=0.5;}
    return volumetric_integration_step*volumetric_step;}

    T Source_Term(const int source_term_index,const VECTOR<T,3>& location) const PHYSBAM_OVERRIDE
    {return voxel_source_interpolation->Clamped_To_Array_Cell(octree_grid,*data(source_term_index),&(*data(source_term_index)),location);}

    bool Use_Precomputed_Light_Data(const VECTOR<T,3>& location,const int light) const PHYSBAM_OVERRIDE
    {if(!precompute_single_scattering)return false;
    OCTREE_CELL<T>* cell=octree_grid.Leaf_Cell(location);
    for(int i=0;i<8;i++)if(!(*precomputed_light_valid(light))(cell->Node(i)))return false;
    return true;}

    VECTOR<T,3> Precomputed_Light_Data(const VECTOR<T,3>& location,const int light) const PHYSBAM_OVERRIDE
    {return voxel_light_interpolation->Clamped_To_Array_Cell(octree_grid,*precomputed_light(light),&(*precomputed_light(light)),location);} // TODO:  Someone should sanity check this.

    void Set_Custom_Source_Interpolation(INTERPOLATION_DYADIC<OCTREE_GRID<T>,T>* interpolation)
    {voxel_source_interpolation=interpolation;}

    void Set_Custom_Light_Interpolation(INTERPOLATION_DYADIC<OCTREE_GRID<T>,VECTOR<T,3> >* interpolation)
    {voxel_light_interpolation=interpolation;}

    void Get_Node_Locations(ARRAY<VECTOR<T,3> >& locations) PHYSBAM_OVERRIDE
    {node_locations.Resize(octree_grid.number_of_nodes,false,false);
    for(int i=1;i<=octree_grid.number_of_nodes;i++) node_locations(i)=octree_grid.Node_Location(i);
    locations=node_locations;}

    void Set_Precomputed_Light_Data(const int location_index,const int light_index,const VECTOR<T,3>& light_value) PHYSBAM_OVERRIDE
    {(*precomputed_light(light_index))(location_index)=light_value;}

    void Set_Precomputed_Light_Valid(const int location_index,const int light_index,const bool value) PHYSBAM_OVERRIDE
    {if((*precomputed_light_ever_valid(light_index))(location_index)) (*precomputed_light_valid(light_index))(location_index)=true;}

protected:
    void Prepare_For_Precomputation(RENDER_WORLD<T>& world) PHYSBAM_OVERRIDE
    {LOG::cout<<"-->RENDERING_OCTREE_VOXELS::Prepare_For_Precomputation"<<std::endl;
        ARRAY<RENDERING_LIGHT<T>*> lights=world.Lights();
    precomputed_light_valid.Resize(lights.m);precomputed_light.Resize(lights.m);precomputed_light_ever_valid.Resize(lights.m);
    LOG::cout<<"Number of lights: "<<lights.m<<std::endl;
    for(int i=1;i<=lights.m;i++){
        precomputed_light_valid(i)=new ARRAY<bool>(octree_grid.number_of_nodes);
        precomputed_light_ever_valid(i)=new ARRAY<bool>(octree_grid.number_of_nodes,false);
        precomputed_light_ever_valid(i)->Fill(true);
        precomputed_light(i)=new ARRAY<VECTOR<T,3> >(octree_grid.number_of_nodes);}
    LOG::cout<<"<--RENDERING_OCTREE_VOXELS::Prepare_For_Precomputation"<<std::endl;}

//#####################################################################
};
}
#endif
#endif

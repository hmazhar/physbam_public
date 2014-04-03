#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2006, Jiayi Chong, Igor Neverov, Andrew Selle, Mike Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE, using the Hierarchical Rendering Technique for Translucent Materials
//#####################################################################
#ifndef __SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE__
#define __SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE__

#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Collisions/POINT_REPULSION.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE.h>
namespace PhysBAM{

template<class T>
class SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;
public:
    ARRAY<SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE<T> > samples;
    OCTREE_GRID<T> octree_grid;
    OCTREE_CELL<T> *root_cell;
    ARRAY<ARRAY<int> > octree_cell_samples;
    ARRAY<SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE<T> > aggregate_samples;
    RANGE<TV> bounding_box;
    int max_irradiance_leaf_samples;
    T error_criterion;
    const T large_number;
    
    SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE(int input_max_irradiance_leaf_samples=8,T input_error_criterion=0.1)
        :root_cell(0),large_number(1000000)
    {max_irradiance_leaf_samples=input_max_irradiance_leaf_samples;error_criterion=input_error_criterion;}

    virtual ~SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE()
    {}
    
    T Approximate_Maximum_Solid_Angle(T area,TV position1,TV position2) 
    {if(position1!=position2)return (T)(area/(position1-position2).Magnitude_Squared());
    LOG::cerr<<"Returning large number!"<<std::endl;
    return large_number;}

//#####################################################################
    template<class RW> void Initialize_Sample_Locations_From_File(const std::string& file,TRIANGULATED_SURFACE<T>& triangulated_surface);
    void Initialize_Sample_Locations_From_Vertices(const TRIANGULATED_SURFACE<T>& triangulated_surface);
    void Build_Octree_Cells(OCTREE_CELL<T>* cell);
    void Build_Aggregate_Samples(OCTREE_CELL<T>* dummy_cell);
    void Build_Octree();
    void Process_Irradiance_Candidates(ARRAY<SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE<T> >& candidate_list,const TV& point,OCTREE_CELL<T>& cell,ARRAY<int>& sample_numbers);
    void Find_Irradiance_Candidates(ARRAY<SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE<T> >& candidate_list,const TV& point);
//#####################################################################
};
}
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Read_Write/Rendering_Shaders/READ_WRITE_SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE.h>
#endif
#endif

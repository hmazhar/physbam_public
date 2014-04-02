#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __IRRADIANCE_CACHE__
#define __IRRADIANCE_CACHE__

#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Random_Numbers/PIECEWISE_CONSTANT_PDF.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/IRRADIANCE_SAMPLE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
namespace PhysBAM{

template<class T> class RENDER_WORLD;
template<class T> class RENDERING_OBJECT;
template<class T> class RENDERING_RAY;

template<class T>
class IRRADIANCE_CACHE:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;
public:
    OCTREE_GRID<T> octree_grid;
    OCTREE_CELL<T> *root_cell;
    ARRAY<ARRAY<IRRADIANCE_SAMPLE<T> >,VECTOR<int,1> > samples;
    T error_threshold;
    bool cache_enabled;
    int samples_stored;
    int sqrt_of_final_gather_samples;
    int final_gather_samples;
    T one_over_final_gather_samples;
    T min_radius_of_influence,max_radius_of_influence;
    GRID<VECTOR<T,2> > importance_sample_grid;
    PIECEWISE_CONSTANT_PDF<T> pdf;
    bool use_photon_map_multiple_importance_sampling;
    RANDOM_NUMBERS<T> random;

    IRRADIANCE_CACHE()
        :cache_enabled(false),use_photon_map_multiple_importance_sampling(false)
    {
        const int photon_map_pdf_phi_samples=16,photon_map_pdf_theta_samples=4;
        importance_sample_grid.Initialize(photon_map_pdf_phi_samples,photon_map_pdf_theta_samples,0,1,0,1,true);
        pdf.Initialize(photon_map_pdf_phi_samples*photon_map_pdf_theta_samples);
    }

    void Set_Final_Gather_Samples(const int final_gather_samples_input)
    {sqrt_of_final_gather_samples=(int)sqrt((T)final_gather_samples_input);final_gather_samples=sqr(sqrt_of_final_gather_samples);one_over_final_gather_samples=(T)1/(T)final_gather_samples;}
                                  
    void Enable_Cache(const RANGE<TV>& bounding_box,const T error_threshold_input)
    {samples.Resize(1,8);
    RANGE<TV> pseudo_bounding_box(bounding_box.Minimum_Corner(),bounding_box.Maximum_Corner()+bounding_box.Edge_Lengths());
    octree_grid.Initialize(GRID<TV>(3,3,3,pseudo_bounding_box),10,0,true,false);root_cell=octree_grid.cells(1,1,1);
    cache_enabled=true;
    error_threshold=error_threshold_input;
    min_radius_of_influence=(T).001*pow(octree_grid.Domain().Size(),(T)one_third);max_radius_of_influence=(T).125*pow(octree_grid.Domain().Size(),(T)one_third);}

//#####################################################################
    bool Irradiance_Estimate(const TV& location,const TV& direction,TV& interpolated_irradiance,const RENDERING_RAY<T>& ray);
    void Irradiance_Estimate_Octree(OCTREE_CELL<T>* cell,const TV& location,const TV& direction,const RENDERING_RAY<T>& ray,int& found_samples,
        TV& weighted_irradiance_accumulator,T& weight_accumulator,const PLANE<T>& plane);
    void Store_Irradiance_Estimate(const TV& location,const TV& direction,const TV& irradiance,const T harmonic_mean_distance);
    void Store_Irradiance_Estimate_In_Octree(OCTREE_CELL<T>* cell,const RANGE<TV>& box,const T box_diagonal_squared,const IRRADIANCE_SAMPLE<T>& sample);
    TV Compute_Indirect_Light(RENDER_WORLD<T>& world,const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const TV& same_side_position,const TV& same_side_normal);
//#####################################################################
};
}
#endif
#endif

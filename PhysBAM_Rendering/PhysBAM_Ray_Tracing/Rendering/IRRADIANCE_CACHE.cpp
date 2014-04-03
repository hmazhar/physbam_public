#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IRRADIANCE_CACHE
//##################################################################### 
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_CELL.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/IRRADIANCE_CACHE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/PHOTON.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY_DEBUG.h>
using namespace PhysBAM;
//#####################################################################
// Class Irradiance_Estimate
//##################################################################### 
template<class T> bool IRRADIANCE_CACHE<T>::
Irradiance_Estimate(const TV& location,const TV& direction,TV& interpolated_irradiance,const RENDERING_RAY<T>& ray)
{
    if(cache_enabled){
        int found_samples=0;TV weighted_irradiance_accumulator(0,0,0);T weight_accumulator=0;PLANE<T> tangent_plane(direction,location);
        Irradiance_Estimate_Octree(root_cell,location,direction,ray,found_samples,weighted_irradiance_accumulator,weight_accumulator,tangent_plane);  
        if(found_samples>0){interpolated_irradiance=weighted_irradiance_accumulator/weight_accumulator;return true;}
        return false;}
    return false;
    interpolated_irradiance=TV(0,0,0);return true; // irradiance cache disabled return 0 for indirect light and treat it as valid so no indirect montecarlo...
}

template<class T> void IRRADIANCE_CACHE<T>::
Irradiance_Estimate_Octree(OCTREE_CELL<T>* cell,const TV& location,const TV& direction,const RENDERING_RAY<T>& ray,int& found_samples,TV& weighted_irradiance_accumulator,
                                T& weight_accumulator,const PLANE<T>& tangent_plane)
{
    // check each irradiance sample to see how much it contributes
    ARRAY<IRRADIANCE_SAMPLE<T> >& sample_list=samples(cell->Cell());
    for(int i=1;i<=sample_list.m;i++){
        IRRADIANCE_SAMPLE<T>& sample=sample_list(i);
        // Culling tests very similar to PBRT's  (Greg Humphrey and Matt Pharr)
        T distance_squared=(location-sample.location).Magnitude_Squared();
        if(distance_squared>sqr(sample.harmonic_mean_distance))continue; // reject too far
        if(TV::Dot_Product(location-sample.location,(sample.direction+direction))<T(-.01))continue; // reject big difference in orientation
        T error_denominator=sample.harmonic_mean_distance*TV::Dot_Product(direction,sample.direction);
        if(error_denominator!=0){
            T error=distance_squared/(error_denominator);
            if(error<1){
                T weight=sqr(1-error);weighted_irradiance_accumulator+=sample.irradiance*weight;weight_accumulator+=weight;found_samples++;}}}
    // traverse down to cells that the point is in.
    if(cell->Has_Children())for(int k=0;k<8;k++){
        OCTREE_CELL<T>* child=cell->Child(k);
        if(octree_grid.Inside_Cell(child,location)){
            Irradiance_Estimate_Octree(child,location,direction,ray,found_samples,weighted_irradiance_accumulator,weight_accumulator,tangent_plane);return;}}
}
//#####################################################################
// Class Store_Irradiance_Estimate
//##################################################################### 
template<class T> void IRRADIANCE_CACHE<T>::
Store_Irradiance_Estimate(const TV& location,const TV& direction,const TV& irradiance,const T harmonic_mean_distance)
{
    if(cache_enabled){
        T influence_distance=clamp(harmonic_mean_distance,min_radius_of_influence,max_radius_of_influence);
        influence_distance*=error_threshold;
        IRRADIANCE_SAMPLE<T> sample(location,direction,irradiance,influence_distance);
        TV radius_vector(influence_distance,influence_distance,influence_distance);
        RANGE<TV> sample_box(location-radius_vector,location+radius_vector);T box_diagonal_squared=sqr(2*influence_distance);
        Store_Irradiance_Estimate_In_Octree(root_cell,sample_box,box_diagonal_squared,sample);}
}
//#####################################################################
// Class Store_Irradiance_Estimate_In_Octree
//##################################################################### 
template<class T> void IRRADIANCE_CACHE<T>::
Store_Irradiance_Estimate_In_Octree(OCTREE_CELL<T>* cell,const RANGE<TV>& box,const T box_diagonal_squared,const IRRADIANCE_SAMPLE<T>& sample)
{
    // Make sure the box overlaps this
    TV center=cell->Center(),dx_over_two=cell->DX()*T(0.5),lower=center-dx_over_two,upper=center+dx_over_two;
    if(!RANGE<TV>(lower,upper).Lazy_Intersection(box))return;
    if(cell->DX().Magnitude_Squared()<box_diagonal_squared||cell->Depth_Of_This_Cell()>octree_grid.maximum_depth){
        if(samples.counts.x<cell->Cell())samples.Resize(1,octree_grid.number_of_cells*2);samples(cell->Cell()).Append(sample);return;}
    if(!cell->Has_Children())cell->Create_Children(octree_grid.number_of_cells,0,octree_grid.number_of_nodes,0,octree_grid.number_of_faces,0,&octree_grid);
    for(int k=0;k<8;k++)Store_Irradiance_Estimate_In_Octree(cell->Child(k),box,box_diagonal_squared,sample);
}
//#####################################################################
// Class Compute_Indirect_Light
//##################################################################### 
template<class T> VECTOR<T,3> IRRADIANCE_CACHE<T>::
Compute_Indirect_Light(RENDER_WORLD<T>& world,const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
                                    const RENDERING_OBJECT<T>& intersection_object,const TV& same_side_position,const TV& same_side_normal)
{
    TV indirect_irradiance(0,0,0);
    if(!Irradiance_Estimate(same_side_position,same_side_normal,indirect_irradiance,ray)){
        // setup matrix to transform the canonical hemisphere directions
        TV temporary_vector=same_side_normal.Orthogonal_Vector();
        TV u=TV::Cross_Product(temporary_vector,same_side_normal);u.Normalize();
        TV v=TV::Cross_Product(same_side_normal,u);
        MATRIX<T,3> normal_transformer(u,v,same_side_normal);
        if(use_photon_map_multiple_importance_sampling){
            T sum_of_reciprical_distances=0;int hit_rays=0;
            MATRIX<T,3> normal_transformer_inverse=normal_transformer.Inverse();
            int photon_count;T max_distance_of_photons_found;
            ARRAY<PHOTON<T>*> photons(world.number_of_photons_for_estimate);ARRAY<T> distances(world.number_of_photons_for_estimate);
            world.global_photon_map.Locate_Photons(same_side_position,world.max_photon_distance,photons,distances,photon_count,max_distance_of_photons_found);
            T total_power=0;
            for(int photon_index=1;photon_index<=photon_count;photon_index++){
                PHOTON<T>& photon=*photons(photon_index);TV direction=normal_transformer_inverse*photon.direction;T phi=atan(direction.y/direction.x);
                VECTOR<T,2> sample_coordinate(1-direction.z*direction.z,phi/(2*T(pi))); // 1-z*z corresponds to 1-cos(theta)^2
                T photon_power=photon.power.Max(); // need to convert to scalar
                VECTOR<int,2> index=importance_sample_grid.Clamp_To_Cell(sample_coordinate);pdf.pdf((index.y-1)*16+index.x)+=photon_power;total_power+=photon_power;}
            for(int i=1;i<=importance_sample_grid.counts.x*importance_sample_grid.counts.y;i++)if(pdf.pdf(i)==0)pdf.pdf(i)+=T(1e-3)*total_power;
            pdf.Compute_Cumulative_Distribution_Function();
            for(int sample_number=1;sample_number<=final_gather_samples;sample_number++){
                PAIR<int,T> sample=pdf.Sample(random.Get_Uniform_Number((T)0,(T)1));
                int i=sample.x%importance_sample_grid.counts.x;int j=sample.x/importance_sample_grid.counts.x;
                VECTOR<T,2> sample_position=importance_sample_grid.X(i,j)+importance_sample_grid.dX*VECTOR<T,2>(random.Get_Uniform_Number((T)0,(T)1),random.Get_Uniform_Number((T)0,(T)1))-VECTOR<T,2>(0.5,0.5);
                T theta=asin(sqrt(sample_position.x));T phi=2*T(pi)*sample_position.y;
                TV pseudo_direction(cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta));
                TV reflected_direction=normal_transformer*pseudo_direction;
                RENDERING_RAY<T> ray_to_light(RAY<VECTOR<T,3> >(same_side_position,reflected_direction),1,&exiting_object);
                TV biased_radiance=clamp(TV::Dot_Product(ray_to_light.ray.direction,same_side_normal),T(0),T(1))*world.Cast_Ray_For_Photon_Gather(ray_to_light,ray);
                TV unbiased_radiance=biased_radiance*T(total_power/(pdf.pdf(sample.x)*importance_sample_grid.counts.x*importance_sample_grid.counts.y));
                indirect_irradiance+=unbiased_radiance;
                if(ray_to_light.ray.semi_infinite==false){hit_rays++;sum_of_reciprical_distances+=T(1)/T(ray_to_light.ray.t_max);}}
            indirect_irradiance*=one_over_final_gather_samples*T(pi);
            // store
            T harmonic_mean_distance;if(hit_rays==0)harmonic_mean_distance=0;else harmonic_mean_distance=T(hit_rays)/T(sum_of_reciprical_distances);
            Store_Irradiance_Estimate(same_side_position,same_side_normal,indirect_irradiance,harmonic_mean_distance);}
        else{
            T sum_of_reciprical_distances=0;int hit_rays=0;
            const int N=sqrt_of_final_gather_samples;const int M=sqrt_of_final_gather_samples;const T one_over_N=T(1)/T(N);const T one_over_M=T(1)/T(M);
            // compute array of radiance and distances (and accumualted return value)
            ARRAY<TV,VECTOR<int,2> > radiance(1,M,1,N);ARRAY<T,VECTOR<int,2> > distances(1,M,1,N);
            // compute irradiance and reciprocal mean distance
            for(int i=1;i<=N;i++)for(int j=1;j<=M;j++){
                T theta=asin(sqrt((j-random.Get_Uniform_Number(T(0),T(1)))*one_over_M));
                T phi=T(2)*T(pi)*(i-random.Get_Uniform_Number(T(0),T(1)))*one_over_N;
                TV pseudo_direction(cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta));
                TV reflected_direction=normal_transformer*pseudo_direction;
                RENDERING_RAY<T> ray_to_light(RAY<VECTOR<T,3> >(same_side_position,reflected_direction),1,&exiting_object);
                radiance(j,i)=clamp(TV::Dot_Product(ray_to_light.ray.direction,same_side_normal),T(0),T(1))*world.Cast_Ray_For_Photon_Gather(ray_to_light,ray);
                indirect_irradiance+=radiance(j,i);
                if(ray_to_light.ray.semi_infinite==false){hit_rays++;sum_of_reciprical_distances+=T(1)/T(ray_to_light.ray.t_max);distances(j,i)=ray_to_light.ray.t_max;}
                else distances(j,i)=FLT_MAX;}
            indirect_irradiance*=one_over_M*one_over_N*T(pi);
            // store
            T harmonic_mean_distance;if(hit_rays==0)harmonic_mean_distance=0;else harmonic_mean_distance=T(hit_rays)/T(sum_of_reciprical_distances);
            Store_Irradiance_Estimate(same_side_position,same_side_normal,indirect_irradiance,harmonic_mean_distance);}
    }
    return indirect_irradiance;
}
//##################################################################### 

template class IRRADIANCE_CACHE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class IRRADIANCE_CACHE<double>;
#endif
#endif

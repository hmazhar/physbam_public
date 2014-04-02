//#####################################################################
// Copyright 2004, Ronald Fedkiw, Frank Losasso, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Structure of this class is inspired by Henrik Wann Jackassen's "Realistic Synthesis Using Photon Mapping" Appendix B
#include <PhysBAM_Tools/Arrays_Computations/HEAPIFY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/integer_log.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_3D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/PHOTON_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY_DEBUG.h>
using namespace PhysBAM;
//#####################################################################
// Function Store_Photon
//#####################################################################
template<class T> bool PHOTON_MAP<T>::
Store_Photon(const TV& location,const TV& direction,const TV& power,const int depth)
{
    if(!Light_Emission_Quota_Remains())return false;
    photons_stored++;
    photons(photons_stored).location=location;
    photons(photons_stored).direction=direction;
    photons(photons_stored).power=power;
    photons(photons_stored).kdtree_split_axis=PHOTON<T>::KDTREE_LEAF;
    bounding_box.Enlarge_To_Include_Point(location);
    return true;
}
//#####################################################################
// Function Locate_Photons
//#####################################################################
template<class T> void PHOTON_MAP<T>::
Locate_Photons(const TV& location,const T max_distance_squared,ARRAY<PHOTON<T>*>& photons_found,ARRAY<T>& distance_squared_of_photons_found,
               int& number_photons_found,T& max_distance_squared_of_found_photons)
{
    number_photons_found=0;
    if(photons.m==0)return;
    T temp_max_distance=max_distance_squared;
    Locate_Photons_Helper(1,location,temp_max_distance,number_photons_found,photons_found,distance_squared_of_photons_found);
    if(number_photons_found<=photons_found.m){
        max_distance_squared_of_found_photons=0;
        for(int i=1;i<=photons_found.m;i++) if(distance_squared_of_photons_found(i)>max_distance_squared_of_found_photons)
            max_distance_squared_of_found_photons=distance_squared_of_photons_found(i);}
    else max_distance_squared_of_found_photons=temp_max_distance;
    number_photons_found=min(photons_found.m,number_photons_found);
}
//#####################################################################
// Function Locate_Photons
//#####################################################################
template<class T> void PHOTON_MAP<T>::
Locate_Photons_Helper(const int photon_index,const TV& location,T& max_distance_squared,int& number_of_photons_found,ARRAY<PHOTON<T>*>& photons_found,
                      ARRAY<T>& distance_squared_of_photons_found)
{
    PHOTON<T>& photon=photons(photon_index);
    if(photon_index<photons_stored_over_two){ // make sure node has children
        T axis_distance=location[photon.kdtree_split_axis]-photon.location[photon.kdtree_split_axis];
        if(axis_distance>0){ // point belongs on right subtree is dominant
            Locate_Photons_Helper(2*photon_index+1,location,max_distance_squared,number_of_photons_found,photons_found,distance_squared_of_photons_found);
            if(sqr(axis_distance)<max_distance_squared)
                Locate_Photons_Helper(2*photon_index,location,max_distance_squared,number_of_photons_found,photons_found,distance_squared_of_photons_found);}
        else{ // point belongs on left subtree
            Locate_Photons_Helper(2*photon_index,location,max_distance_squared,number_of_photons_found,photons_found,distance_squared_of_photons_found);
            if(sqr(axis_distance)<max_distance_squared)
                Locate_Photons_Helper(2*photon_index+1,location,max_distance_squared,number_of_photons_found,photons_found,distance_squared_of_photons_found);}}
    T distance_squared_to_query_location=(photon.location-location).Magnitude_Squared();
    if(distance_squared_to_query_location<=max_distance_squared){ //&&photon.depth>1){
        if(number_of_photons_found<photons_found.m){
            number_of_photons_found++;
            photons_found(number_of_photons_found)=&photon;distance_squared_of_photons_found(number_of_photons_found)=distance_squared_to_query_location;}  
        else{
            if(number_of_photons_found==photons_found.m){ARRAYS_COMPUTATIONS::Heapify(distance_squared_of_photons_found,photons_found);number_of_photons_found++;}
            //photons_found(1)=&photon;distance_squared_of_photons_found(1)=distance_squared_to_query_location;
            // can't use heapify here, because we need to break at the position we place the new photon
            int current_index=1;
            int left,right,index_of_largest;
            for(;;){
                left=2*current_index;right=2*current_index+1;index_of_largest=current_index;
                if(left>photons_found.m)break;
                else if(right>photons_found.m||distance_squared_of_photons_found(left)>distance_squared_of_photons_found(right)) index_of_largest=left; 
                else index_of_largest=right;
                if(distance_squared_to_query_location>distance_squared_of_photons_found(index_of_largest))break; // we found a place to put our new guy
                exchange(distance_squared_of_photons_found(current_index),distance_squared_of_photons_found(index_of_largest));
                exchange(photons_found(current_index),photons_found(index_of_largest));
                current_index=index_of_largest;}
            distance_squared_of_photons_found(current_index)=distance_squared_to_query_location;
            photons_found(current_index)=&photon;
            max_distance_squared=distance_squared_of_photons_found(1);}} // copy the new bounding distance
}
//#####################################################################
// Function Irradiance_Estimate
//#####################################################################
template<class T> VECTOR<T,3> PHOTON_MAP<T>::
Irradiance_Estimate(const TV& location,const TV& normal,const T max_distance_squared,const int number_of_photons,const RENDERING_RAY<T>& ray,const PHOTON_MAP_TYPE type)
{
    int found_samples;
    TV irradiance;
    ARRAY<PHOTON<T>*> nearby_photons(number_of_photons);
    ARRAY<T> photon_squared_distances(number_of_photons);
    T one_over_max_distance,actual_max_distance_squared,max_distance_cubed;
    Locate_Photons(location,max_distance_squared,nearby_photons,photon_squared_distances,found_samples,actual_max_distance_squared);
    if(ray.debug_ray) ray.debug_ray->Add_Comment(STRING_UTILITIES::string_sprintf("Got %d photons",found_samples));
    if(found_samples<4) return TV();
    one_over_max_distance=T(1)/sqrt(actual_max_distance_squared);
    max_distance_cubed=actual_max_distance_squared*sqrt(actual_max_distance_squared);
    if(ray.debug_ray) ray.debug_ray->Add_Comment(STRING_UTILITIES::string_sprintf("Actual Max Distance Squared %f",actual_max_distance_squared));
    //const T alpha=0.918,beta=1.953;
    T one_over_k=T(1)/T(1.1);
    for(int i=1;i<=found_samples;i++){
        PHOTON<T>& photon=*nearby_photons(i);
        //T filter_weight=alpha*(1-(1-exp(-beta*photon_squared_distances(i)/(2*actual_max_distance_squared)))/(1-exp(-beta)));
        T filter_weight=1-sqrt(photon_squared_distances(i))*one_over_max_distance*one_over_k;
        if(ray.debug_ray)ray.debug_ray->photons_used.Append(nearby_photons(i));
        if(TV::Dot_Product(photon.direction,normal)<0)irradiance+=photon.power*filter_weight;}
    TV total_irradiance;
    if(type==PHOTON_MAP<T>::VOLUME_PHOTON_MAP) total_irradiance=irradiance/((T(4)/T(3))*T(pi)*max_distance_cubed)/(1-T(two_thirds)*one_over_k);
    else total_irradiance=irradiance/(T(pi)*actual_max_distance_squared)/(1-T(two_thirds)*one_over_k);
    if(ray.debug_ray) ray.debug_ray->Add_Comment(STRING_UTILITIES::string_sprintf("Irradiance Estimate from Photon Map %f %f %f",total_irradiance.x,total_irradiance.y,total_irradiance.z));
    return total_irradiance;
}
//#####################################################################
// Function Max_Number_Of_Photons
//#####################################################################
template<class T> int PHOTON_MAP<T>::
Max_Number_Of_Photons() const
{
    return photons.m;
}
//#####################################################################
// Function Prepare_Photon_Map_For_Rendering
//#####################################################################
template<class T> void PHOTON_MAP<T>::
Prepare_Photon_Map_For_Rendering()
{
    Construct_Balanced_KD_Tree();
}
//#####################################################################
// Function Construct_Balanced_KD_Tree 
//#####################################################################
template<class T> void PHOTON_MAP<T>::
Construct_Balanced_KD_Tree()
{
    if(photons_stored<=0)return;
    ARRAY<int> permutation_array(photons.m),kdtree_array(photons.m);
    for(int i=1;i<=photons.m;i++)permutation_array(i)=kdtree_array(i)=i;
    Balance_Sub_KD_Tree(1,1,photons.m,permutation_array,kdtree_array,bounding_box);
    // reorder array (current_index position is always been saved so you can write new value there)
    int current_index=1;int index_of_photon_to_relocate=1;
    PHOTON<T> photon_to_relocate=photons(index_of_photon_to_relocate);
    for(int dummy_counter=1;dummy_counter<=photons.m;dummy_counter++){
        int swap_to=kdtree_array(current_index);kdtree_array(current_index)=-1;
        if(swap_to==index_of_photon_to_relocate){
            photons(current_index)=photon_to_relocate;
            if(dummy_counter<photons.m){
                while(index_of_photon_to_relocate<=photons.m&&kdtree_array(index_of_photon_to_relocate)==-1)index_of_photon_to_relocate++;
                photon_to_relocate=photons(index_of_photon_to_relocate);current_index=index_of_photon_to_relocate;}}
        else{
            photons(current_index)=photons(swap_to);
            current_index=swap_to;}}
    photons_stored_over_two=photons_stored/2-1;
}
//#####################################################################
// Function Balance_KD_Tree
//#####################################################################
template<class T> void PHOTON_MAP<T>::
Balance_Sub_KD_Tree(const int index,const int first_index,const int last_index,ARRAY<int>& permutation_array,ARRAY<int>& kdtree_array,RANGE<TV> box)
{
    if(last_index==first_index){kdtree_array(index)=permutation_array(first_index);return;}
    // Find the partition that will generate a left balanced kd-tree
    int partition_index;
    int elements=last_index-first_index+1;int filled_subtree_height=integer_log(elements+1)-1;int filled_subtree_elements=(1<<filled_subtree_height)-1;
    int extra_elements=elements-2*filled_subtree_elements-1;
    if(extra_elements>=filled_subtree_elements+1)partition_index=first_index+filled_subtree_elements*2+1;
    else partition_index=first_index+filled_subtree_elements+extra_elements;
    // Choose axis to partition on
    T dx,dy,dz;box.Edge_Lengths().Get(dx,dy,dz);
    typename PHOTON<T>::PHOTON_KDTREE_SPLIT_AXIS axis;
    if(dx>dy&&dy>dz)axis=PHOTON<T>::KDTREE_SPLIT_X;
    else if(dy>dx)axis=PHOTON<T>::KDTREE_SPLIT_Y;
    else axis=PHOTON<T>::KDTREE_SPLIT_Z;
    // Now partition about the axis
    Median_Split(partition_index,first_index,last_index,permutation_array,axis);
    photons(permutation_array(partition_index)).kdtree_split_axis=axis;
    kdtree_array(index)=permutation_array(partition_index);
    // Build left and right subtrees
    TV lower_bound=box.Minimum_Corner(),upper_bound=box.Maximum_Corner(),split_photon_position=photons(kdtree_array(index)).location;
    if(partition_index>first_index){
        TV modified_upper_bound=upper_bound;modified_upper_bound[axis]=split_photon_position[axis];
        Balance_Sub_KD_Tree(2*index,first_index,partition_index-1,permutation_array,kdtree_array,RANGE<TV>(lower_bound,modified_upper_bound));}
    if(partition_index<last_index){
        TV modified_lower_bound=lower_bound;modified_lower_bound[axis]=split_photon_position[axis];
        Balance_Sub_KD_Tree(2*index+1,partition_index+1,last_index,permutation_array,kdtree_array,RANGE<TV>(modified_lower_bound,upper_bound));}
}
//#####################################################################
// Function Median_Split
//#####################################################################
template<class T> void PHOTON_MAP<T>::
Median_Split(const int required_partition_index,const int first_index,const int last_index,ARRAY<int>& permutation_array,const typename PHOTON<T>::PHOTON_KDTREE_SPLIT_AXIS axis)
{
    int left=first_index,right=last_index;
    while(right>left){
        int obtained_partition_index=Partition_Sub_Array(left,right,permutation_array,axis);
        if(obtained_partition_index>required_partition_index)right=obtained_partition_index-1;
        else if(obtained_partition_index<required_partition_index)left=obtained_partition_index+1;
        else break;} // i=required_partition_index so we're done
}
//#####################################################################
// Function Partition_Sub_Array (routine from Sedgewick's C++ Algo pg 319)
//#####################################################################
template<class T> int PHOTON_MAP<T>::
Partition_Sub_Array(const int left_index,const int right_index,ARRAY<int>& permutation_array,const typename PHOTON<T>::PHOTON_KDTREE_SPLIT_AXIS axis)
{
    int i=left_index-1,j=right_index;T partition_value=photons(permutation_array(right_index)).location[axis];int temp;
    while(1){
        while(photons(permutation_array(++i)).location[axis]<partition_value); // find left side swap candidate
        while(partition_value<photons(permutation_array(--j)).location[axis]&&j!=left_index); // find right side swap candidate
        if(i>=j)break;
        temp=permutation_array(i);permutation_array(i)=permutation_array(j);permutation_array(j)=temp;}
    temp=permutation_array(i);permutation_array(i)=permutation_array(right_index);permutation_array(right_index)=temp;
    return i; // return the split index
}
//#####################################################################
// Function Print_Photon_List
//#####################################################################
template<class T> void PHOTON_MAP<T>::
Print_Photon_List()
{
    LOG::cout<<"PHOTON LIST"<<std::endl;
    for(int i=1;i<=photons.m;i++)LOG::cout<<"i="<<i<<" Photon "<<photons(i).location<<" Power "<<photons(i).power<<" Split on "<<photons(i).kdtree_split_axis<<std::endl;
}
//#####################################################################
// Function Print_Photon_Tree
//#####################################################################
template<class T> void PHOTON_MAP<T>::
Print_Photon_Tree(const int depth,const int index)
{
    if(index>photons.m)return;
    if(index==1)LOG::cout<<"PHOTON TREE"<<std::endl;
    for(int i=1;i<=depth;i++)LOG::cout<<" ";
    LOG::cout<<"Photon "<<photons(index).location<<" Split on "<<photons(index).kdtree_split_axis<<std::endl;
    Print_Photon_Tree(depth+5,index*2);
    Print_Photon_Tree(depth+5,index*2+1);
}
//#####################################################################
// Function Begin_Light_Emission
//#####################################################################
template<class T> void PHOTON_MAP<T>::
Begin_Light_Emission(const int number_of_photons_for_light)
{
    if(light_emission_start_index!=-1)return;
    light_emission_photon_quota=number_of_photons_for_light;
    light_emission_start_index=photons_stored+1;
}
//#####################################################################
// Function End_Light_Emission
//#####################################################################
template<class T> void PHOTON_MAP<T>::
End_Light_Emission(const int number_of_photons_emitted)
{
    if(light_emission_start_index==-1||number_of_photons_emitted==0)return;
    T power_scale=T(1)/T(number_of_photons_emitted);
    for(int i=light_emission_start_index;i<=photons_stored;i++)photons(i).power*=power_scale;
    light_emission_start_index=-1;
}
//#####################################################################
// Function Light_Emission_Quota_Remains
//#####################################################################
template<class T> bool PHOTON_MAP<T>::
Light_Emission_Quota_Remains()
{
    if(light_emission_start_index==-1)return false;
    int number=photons_stored-light_emission_start_index+1;
    if(number>=light_emission_photon_quota)return false;else return true;
}
//#####################################################################
template class PHOTON_MAP<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PHOTON_MAP<double>;
#endif

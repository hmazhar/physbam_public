//#####################################################################
// Copyright 2004, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PHOTON_MAP__
#define __PHOTON_MAP__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/PHOTON.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY.h>
namespace PhysBAM{

template<class T>
class PHOTON_MAP:public NONCOPYABLE
{
  typedef VECTOR<T,3> TV;
public:
    ARRAY<PHOTON<T> > photons;
    int photons_stored,photons_stored_over_two;
    int light_emission_start_index;
    int light_emission_photon_quota;
public:
    RANGE<TV> bounding_box;
    enum PHOTON_MAP_TYPE{CAUSTIC_PHOTON_MAP,VOLUME_PHOTON_MAP,GLOBAL_PHOTON_MAP};

    PHOTON_MAP(const int maximum_photons)
    {
        Resize_Photon_Map(maximum_photons);
    }

    void Resize_Photon_Map(const int new_maximum_photons)
    {photons.Resize(new_maximum_photons);photons_stored=0;light_emission_start_index=-1;}
    
    TV Photon_Position(const int index) const
    {return photons(index).location;}
    
    TV Photon_Power(const int index) const
    {return photons(index).power;}
    
    TV Photon_Direction(const int index) const
    {return photons(index).direction;}

//#####################################################################
    bool Store_Photon(const TV& location,const TV& direction,const TV& power,const int depth);
    void Locate_Photons(const TV& location,const T max_distance,ARRAY<PHOTON<T>*>& photons_found,ARRAY<T>& distance_squared_of_photons_found,
        int& number_photons_found,T& max_distance_squared_of_found_photons);
    void Locate_Photons_Helper(const int photon_index,const TV& location,T& max_distance_squared,int& number_of_photons_found,ARRAY<PHOTON<T>*>& photons_found,
        ARRAY<T>& distance_squared_of_photons_found);
    TV Irradiance_Estimate(const TV& location,const TV& normal,const T max_distance,const int number_of_photons,const RENDERING_RAY<T>& ray,const PHOTON_MAP_TYPE type);
    void Prepare_Photon_Map_For_Rendering();
    int Max_Number_Of_Photons() const; 
private:
    void Construct_Balanced_KD_Tree();
    void Balance_Sub_KD_Tree(const int index,const int first_index,const int last_index,ARRAY<int>& permutation_array,ARRAY<int>& kdtree_array,RANGE<TV> box);
    void Median_Split(const int required_partition_index,const int first_index,const int last_index,ARRAY<int>& permutation_array,const typename PHOTON<T>::PHOTON_KDTREE_SPLIT_AXIS axis);
    int Partition_Sub_Array(const int left_index,const int right_index,ARRAY<int>& permutation_array,const typename PHOTON<T>::PHOTON_KDTREE_SPLIT_AXIS axis);
    void Print_Photon_List();
    void Print_Photon_Tree(const int depth=0,const int index=1);
public:
    void Begin_Light_Emission(const int number_of_photons_for_light);
    void End_Light_Emission(const int number_of_photons_emitted);
    bool Light_Emission_Quota_Remains();
//#####################################################################
};
}
#endif

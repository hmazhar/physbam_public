//#####################################################################
// Copyright 2003-2004, Ron Fedkiw, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_3D_SPATIAL_PARTITION
//##################################################################### 
//
// A class for quickly getting pairs of particles within an interaction distance. Call Reinitialize to set it up initially and also whenever the particles move around.
// There are two ways to get nearby pairs of particles:
// (1) Given a particle index, you can get a list of other particles within the interaction radius using Get_Particles_Within_Interaction_Radius().
// (2) If you want to loop through all interacting pairs and you don't care in what order, a much more efficient way is:
//      Reset_Pair_Finder();while(Get_Next_Particles_Potentially_Within_Interaction_Radius(index1,nearby_particle_indices)){}
//      The particles in "nearby_particle_indices" only *might* be within radius of index1, so you still need to double check if they are.
//      Doing things this way, you will iterate through all distinct pairs.
//
//##################################################################### 
#ifndef __PARTICLE_3D_SPATIAL_PARTITION__
#define __PARTICLE_3D_SPATIAL_PARTITION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
namespace PhysBAM{

template<class T>
class PARTICLE_3D_SPATIAL_PARTITION:public NONCOPYABLE
{  
private:
    typedef VECTOR<T,3> TV;

    GEOMETRY_PARTICLES<TV>& particles;
    T interaction_radius,interaction_radius_squared,one_over_interaction_radius;
    HASHTABLE<VECTOR<int,3>,ARRAY<int>*> hashtable;
    HASHTABLE_ITERATOR<VECTOR<int,3>,ARRAY<int>*> hashtable_iterator;
    int reinitialize_counter;
    int number_neighbor_lists;
    ARRAY<int> *main_list,*neighbor_lists[28];
    int main_list_index,neighbor_lists_index[28];

public:
    PARTICLE_3D_SPATIAL_PARTITION(GEOMETRY_PARTICLES<TV>& particles_input,T interaction_radius_input)
        :particles(particles_input),interaction_radius(interaction_radius_input),interaction_radius_squared(interaction_radius*interaction_radius),
         hashtable_iterator(hashtable),reinitialize_counter(0)
    {
        assert(interaction_radius>0);
        one_over_interaction_radius=1/interaction_radius;
    }

    ~PARTICLE_3D_SPATIAL_PARTITION()
    {hashtable.Delete_Pointers_Stored_In_Table();}

    void Reinitialize() // ensures that each hashtable entry is a list of indices in increasing order
    {if(reinitialize_counter%10 == 0){hashtable.Delete_Pointers_Stored_In_Table();hashtable.Remove_All();}else hashtable.Reset_List_Arrays_Stored_In_Table();
    reinitialize_counter++;
    for(int k=1;k<=particles.array_collection->Size();k++){
        VECTOR<int,3> voxel=Voxel(k);
        ARRAY<int>* occupancy_list=0;
        if(!hashtable.Get(voxel,occupancy_list)){occupancy_list=new ARRAY<int>();occupancy_list->Preallocate(5);hashtable.Insert(voxel,occupancy_list);}
        occupancy_list->Append(k);}}

    void Get_Particles_Within_Interaction_Radius(const int particle_index,ARRAY<int>& interacting_particle_indices,bool only_higher_indices=false)
    {VECTOR<int,3> voxel=Voxel(particle_index);
    interacting_particle_indices.Remove_All();
    for(int i=voxel.x-1;i<=voxel.x+1;i++) for(int j=voxel.y-1;j<=voxel.y+1;j++) for(int k=voxel.z-1;k<=voxel.z+1;k++){
        ARRAY<int>* occupancy_list=0;
        if(hashtable.Get(VECTOR<int,3>(i,j,k),occupancy_list)) for(int t=1;t<=occupancy_list->m;t++){
            int index=(*occupancy_list)(t);
            if((index > particle_index || (!only_higher_indices  && index != particle_index))&&(particles.X(particle_index)-particles.X(index)).Magnitude_Squared() <= interaction_radius_squared) 
                interacting_particle_indices.Append(index);}}}

    void Reset_Pair_Finder()
    {hashtable_iterator.Reset();Initialize_Next_Voxels_To_Process();}

    bool Get_Next_Particles_Potentially_Within_Interaction_Radius(int &particle_index,ARRAY<int>& interacting_particle_indices)
    {interacting_particle_indices.Remove_All();
    if(main_list_index > main_list->m && !Initialize_Next_Voxels_To_Process()) return false;
    particle_index=(*main_list)(main_list_index);
    for(int t=main_list_index+1;t<=main_list->m;t++) interacting_particle_indices.Append((*main_list)(t));
    for(int i=0;i<number_neighbor_lists;i++){
        assert(neighbor_lists[i]);
        while(neighbor_lists_index[i] <= neighbor_lists[i]->m && particle_index > (*neighbor_lists[i])(neighbor_lists_index[i])) neighbor_lists_index[i]++;
        // append everything from neighbor_lists_index[i] on (as long as it's in the radius)
        for(int ii=neighbor_lists_index[i];ii<=neighbor_lists[i]->m;ii++) interacting_particle_indices.Append((*neighbor_lists[i])(ii));}
    main_list_index++;
    return true;}

private:
    VECTOR<int,3> Voxel(const int particle_index)
    {const TV& position=particles.X(particle_index);
    return VECTOR<int,3>((int)floor(position.x*one_over_interaction_radius),(int)floor(position.y*one_over_interaction_radius),(int)floor(position.z*one_over_interaction_radius));}

    bool Initialize_Next_Voxels_To_Process()
    {VECTOR<int,3> voxel;if(!Next_Voxel_To_Process(voxel,main_list)) return false;
    assert(main_list->m > 0);   // should at least have one particle in it
    main_list_index=1;int index=0;
    for(int i=voxel.x-1;i<=voxel.x+1;i++) for(int j=voxel.y-1;j<=voxel.y+1;j++) for(int k=voxel.z-1;k<=voxel.z+1;k++){
        if(i==voxel.x&&j==voxel.y&&k==voxel.z) continue;
        ARRAY<int>* occupancy_list=0;
        if(hashtable.Get(VECTOR<int,3>(i,j,k),occupancy_list)){neighbor_lists[index]=occupancy_list;neighbor_lists_index[index]=1;index++;}}
    number_neighbor_lists=index;
    return true;}

    bool Next_Voxel_To_Process(VECTOR<int,3>& voxel,ARRAY<int>*& occupancy_list)
    {for(;;){
        if(!hashtable_iterator.Valid()) return false;
        occupancy_list=hashtable_iterator.Data();if(occupancy_list->m > 0) break;
        hashtable_iterator.Next();}
    voxel=hashtable_iterator.Key();
    hashtable_iterator.Next();
    return true;}

//#####################################################################
};
}
#endif

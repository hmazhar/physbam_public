//#####################################################################
// Copyright 2005, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class AXIAL_COLLISION_HASH
//#####################################################################
#ifndef __AXIAL_COLLISION_HASH__
#define __AXIAL_COLLISION_HASH__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

template<class T>
class AXIAL_COLLISION_HASH:public NONCOPYABLE
{
public:
    int primitives;
    int axis_clusters;
    VECTOR<ARRAY<T>,3> axial_primitive_centroids;
    VECTOR<ARRAY<RANGE<VECTOR<T,1> > >,3> axial_primitive_bounding_boxes;
    VECTOR<ARRAY<RANGE<VECTOR<T,1> > >,3> axial_intervals;
    ARRAY<ARRAY<int>*,VECTOR<int,3> > cluster_grid_hash;
    ARRAY<VECTOR<int,2> > cluster_ranges;
    ARRAY<VECTOR<int,3> > primitive_to_axial_clusters;
    
    AXIAL_COLLISION_HASH()
    {}

    virtual ~AXIAL_COLLISION_HASH()
    {}

    void Initialize(const int primitives_input,const int axis_clusters_input)
    {primitives=primitives_input;axis_clusters=axis_clusters_input;
    for(int axis=1;axis<=3;axis++){
        axial_primitive_centroids[axis].Resize(primitives,false);
        axial_primitive_bounding_boxes[axis].Resize(primitives,false);
        axial_intervals[axis].Resize(axis_clusters,false);}
    cluster_grid_hash.Resize(1,axis_clusters_input,1,axis_clusters_input,1,axis_clusters_input,false);
    primitive_to_axial_clusters.Resize(primitives);
    // generate ranges
    cluster_ranges.Resize(axis_clusters,false,false);
    int quotient=primitives/axis_clusters,remainder=primitives%axis_clusters;
    for(int i=1;i<=axis_clusters;i++) cluster_ranges(i)=VECTOR<int,2>((i-1)*quotient+min(i-1,remainder)+1,i*quotient+min(i,remainder));}

    void Update_Grid_Hash()
    {ARRAY<int> primitive_indices(primitives,false);
    for(int axis=1;axis<=3;axis++){
        for(int i=1;i<=primitives;i++) primitive_indices(i)=i;
        Sort(primitive_indices,Indirect_Comparison(axial_primitive_centroids[axis]));
        ARRAY<RANGE<VECTOR<T,1> > >::Copy(RANGE<VECTOR<T,1> >(FLT_MAX,-FLT_MAX),axial_intervals[axis]);
        for(int cluster=1;cluster<=axis_clusters;cluster++) for(int i=cluster_ranges(cluster).x;i<=cluster_ranges(cluster).y;i++){
            int primitive_index=primitive_indices(i);
            axial_intervals[axis](cluster).Enlarge_To_Include_Box(axial_primitive_bounding_boxes[axis](primitive_index));
            primitive_to_axial_clusters(primitive_index)[axis]=cluster;}}
    // Build cartesian product
    for(int i=1;i<=axis_clusters;i++) for(int j=1;j<=axis_clusters;j++) for(int ij=1;ij<=axis_clusters;ij++) if(cluster_grid_hash(i,j,ij)){delete cluster_grid_hash(i,j,ij);cluster_grid_hash(i,j,ij)=0;}
    for(int i=1;i<=primitives;i++){
        VECTOR<int,3>& index=primitive_to_axial_clusters(i);
        if(!cluster_grid_hash(index)) cluster_grid_hash(index)=new ARRAY<int>;
        cluster_grid_hash(index)->Append(i);}}

//#####################################################################
};   
}
#endif

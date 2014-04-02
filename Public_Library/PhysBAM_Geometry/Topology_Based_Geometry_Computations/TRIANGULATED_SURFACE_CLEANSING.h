//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGULATED_SURFACE_CLEANSING__
#define __TRIANGULATED_SURFACE_CLEANSING__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_3D_SPATIAL_PARTITION.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS{

template<class T> void
Close_Surface(TRIANGULATED_SURFACE<T>& ts,const bool merge_coincident_vertices,const T merge_coincident_vertices_threshold,const bool fill_holes,const bool verbose)
{
    typedef VECTOR<T,3> TV;
    PHYSBAM_ASSERT(merge_coincident_vertices_threshold<1);
    bool incident_elements_defined=(ts.mesh.incident_elements!=0);if(!incident_elements_defined) ts.mesh.Initialize_Incident_Elements();
    bool boundary_mesh_defined=(ts.mesh.boundary_mesh!=0);if(!boundary_mesh_defined) ts.mesh.Initialize_Boundary_Mesh();
    bool node_on_boundary_defined=(ts.mesh.node_on_boundary!=0);if(!node_on_boundary_defined) ts.mesh.Initialize_Node_On_Boundary();

    // for each node on boundary, get minimum length of boundary segments adjacent to it - also keep track of longest boundary segment
    ARRAY<T> minimum_incident_boundary_segment_length(ts.mesh.boundary_mesh->number_nodes,false);
    ARRAYS_COMPUTATIONS::Fill(minimum_incident_boundary_segment_length,(T)FLT_MAX);
    T maximum_boundary_segment_length=0;
    for(int i=1;i<=ts.mesh.boundary_mesh->elements.m;i++){
        int node1,node2;ts.mesh.boundary_mesh->elements(i).Get(node1,node2);
        T length=(ts.particles.X(node1)-ts.particles.X(node2)).Magnitude();
        minimum_incident_boundary_segment_length(node1)=min(minimum_incident_boundary_segment_length(node1),length);
        minimum_incident_boundary_segment_length(node2)=min(minimum_incident_boundary_segment_length(node2),length);
        maximum_boundary_segment_length=max(maximum_boundary_segment_length,length);}

    if(merge_coincident_vertices && maximum_boundary_segment_length>0){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(verbose) LOG::cout<<"MERGING VERTICES (threshold="<<merge_coincident_vertices_threshold<<")"<<std::endl;
#endif
        int number_merged=0;
        PARTICLE_3D_SPATIAL_PARTITION<T> particle_spatial_partition(ts.particles,2*maximum_boundary_segment_length);
        particle_spatial_partition.Reinitialize();particle_spatial_partition.Reset_Pair_Finder();
        int index1;ARRAY<int> nearby_particle_indices;nearby_particle_indices.Preallocate(100);
        while(particle_spatial_partition.Get_Next_Particles_Potentially_Within_Interaction_Radius(index1,nearby_particle_indices)) if((*ts.mesh.node_on_boundary)(index1)){
            TV position1=ts.particles.X(index1);
            int closest_index2=0;T closest_distance_squared=FLT_MAX;
            for(int k=1;k<=nearby_particle_indices.m;k++){
                int index2=nearby_particle_indices(k);if(!(*ts.mesh.node_on_boundary)(index2)) continue;
                T distance_squared=(ts.particles.X(index2)-position1).Magnitude_Squared();
                if(distance_squared < closest_distance_squared){closest_index2=index2;closest_distance_squared=distance_squared;}}
            if(closest_index2){ // merge index1 to closest_index2
                T real_threshold=merge_coincident_vertices_threshold*min(minimum_incident_boundary_segment_length(index1),minimum_incident_boundary_segment_length(closest_index2));
                if(closest_distance_squared<=sqr(real_threshold)){
                    number_merged++;
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
                    if(verbose) LOG::cout<<index1<<"->"<<closest_index2<<" "<<std::flush;
#endif
                    for(int j=1;j<=(*ts.mesh.incident_elements)(index1).m;j++){
                        int t=(*ts.mesh.incident_elements)(index1)(j);
                        if(ts.mesh.elements(t)(1)==index1)ts.mesh.elements(t)(1)=closest_index2;
                        else if(ts.mesh.elements(t)(2)==index1)ts.mesh.elements(t)(2)=closest_index2;
                        else if(ts.mesh.elements(t)(3)==index1)ts.mesh.elements(t)(3)=closest_index2;}
                    (*ts.mesh.incident_elements)(closest_index2).Append_Elements((*ts.mesh.incident_elements)(index1));}}} // dynamically update the incident triangles list
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(verbose) LOG::cout<<std::endl<<number_merged<<" vertices merged"<<std::endl<<std::endl;
#endif
        ts.Refresh_Auxiliary_Structures();}

    if(fill_holes){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(verbose) LOG::cout<<"HOLE FILLING"<<std::endl;
#endif
        bool connected_segments_defined=(ts.mesh.boundary_mesh->connected_segments!=0);if(!connected_segments_defined) ts.mesh.boundary_mesh->Initialize_Connected_Segments();
        ARRAY<ARRAY<VECTOR<int,2> > >& connected_segments=*ts.mesh.boundary_mesh->connected_segments;
        for(int i=1;i<=connected_segments.m;i++){
            TV centroid;
            for(int j=1;j<=connected_segments(i).m;j++){int node1,node2;connected_segments(i)(j).Get(node1,node2);centroid+=ts.particles.X(node1)+ts.particles.X(node2);}
            centroid/=(T)(2*connected_segments(i).m); // assuming we count each node exactly twice, this gives us the average
            int new_particle_index=ts.particles.array_collection->Add_Element();ts.mesh.number_nodes++;ts.particles.X(new_particle_index)=centroid;
            if(ts.particles.store_velocity){
                TV velocity;
                for(int j=1;j<=connected_segments(i).m;j++){int node1,node2;connected_segments(i)(j).Get(node1,node2);velocity+=ts.particles.V(node1)+ts.particles.V(node2);}
                ts.particles.V(new_particle_index)=velocity/(T)(2*connected_segments(i).m);}
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(verbose){LOG::cout<<"Adding particle "<<new_particle_index<<": "<<ts.particles.X(new_particle_index)<<std::endl;LOG::cout<<"Adding triangles: "<<std::flush;}
#endif
            for(int j=1;j<=connected_segments(i).m;j++){ // assumes segment orientation is consistent with triangle orientation!
                int node1,node2;connected_segments(i)(j).Get(node1,node2);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
                if(verbose) LOG::cout<<"("<<node2<<","<<node1<<","<<new_particle_index<<") "<<std::flush;
#endif
                ts.mesh.elements.Append(VECTOR<int,3>(node2,node1,new_particle_index));}
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(verbose) LOG::cout<<std::endl<<std::endl;
#endif
        }
        if(!connected_segments_defined){delete ts.mesh.boundary_mesh->connected_segments;ts.mesh.boundary_mesh->connected_segments=0;}}

    if(!incident_elements_defined){delete ts.mesh.incident_elements;ts.mesh.incident_elements=0;}
    if(!boundary_mesh_defined){delete ts.mesh.boundary_mesh;ts.mesh.boundary_mesh=0;}
    if(!node_on_boundary_defined){delete ts.mesh.node_on_boundary;ts.mesh.node_on_boundary=0;}

    ts.Discard_Valence_Zero_Particles_And_Renumber();  // refreshes auxiliary structures too 
}
}
}
#endif

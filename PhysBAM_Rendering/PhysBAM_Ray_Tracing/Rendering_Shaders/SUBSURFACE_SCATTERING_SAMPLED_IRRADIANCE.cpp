#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2006, Jiayi Chong, Igor Neverov, Andrew Selle, Mike Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE, using the Hierarchical Rendering Technique for Translucent Materials
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Math_Tools/READ_WRITE_RANGE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE.h>
using namespace PhysBAM;

//#####################################################################
// Function Initialize_Sample_Locations_From_File
//#####################################################################
template<class T> template<class RW> void SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<T>::
Initialize_Sample_Locations_From_File(const std::string& file,TRIANGULATED_SURFACE<T>& triangulated_surface)
{
    ARRAY<VECTOR<T,3> > point_weights;ARRAY<T> original_triangle_areas;ARRAY<int> sample_triangles;
    {std::istream* input=FILE_UTILITIES::Safe_Open_Input(file);
    POINT_REPULSION<T>::template Read<RW>(*input,point_weights,sample_triangles,original_triangle_areas);delete input;}
    PHYSBAM_ASSERT(point_weights.m==sample_triangles.m);PHYSBAM_ASSERT(triangulated_surface.mesh.elements.m==original_triangle_areas.m);
    samples.Resize(point_weights.m);
    T surface_area_per_sample=triangulated_surface.Total_Area()/samples.m; // assume points evenly distributed on surface
    ARRAY<T> new_triangle_areas(triangulated_surface.mesh.elements.m);
    for(int i=1;i<=triangulated_surface.mesh.elements.m;++i){
        int node1,node2,node3;triangulated_surface.mesh.elements(i).Get(node1,node2,node3);
        new_triangle_areas(i)=TRIANGLE_3D<T>::Area(triangulated_surface.particles.X(node1),triangulated_surface.particles.X(node2),triangulated_surface.particles.X(node3));}
    T min_scale,max_scale;min_scale=max_scale=new_triangle_areas(sample_triangles(1))/original_triangle_areas(sample_triangles(1));
    for(int i=1;i<=samples.m;++i){
        int node1,node2,node3;triangulated_surface.mesh.elements(sample_triangles(i)).Get(node1,node2,node3);
        samples(i).position=TRIANGLE_3D<T>::Point_From_Barycentric_Coordinates(point_weights(i),triangulated_surface.particles.X(node1),triangulated_surface.particles.X(node2),triangulated_surface.particles.X(node3));
        // scale surface_area_by_sample by the change in area from the original area (from stretching)
        T scale_factor=new_triangle_areas(sample_triangles(i))/original_triangle_areas(sample_triangles(i));
        if(scale_factor<min_scale)min_scale=scale_factor;if(scale_factor>max_scale)max_scale=scale_factor;
        samples(i).area=surface_area_per_sample*scale_factor;
        samples(i).normal=triangulated_surface.Normal(samples(i).position,sample_triangles(i));}
    LOG::cerr<<"Samples loaded from file: min_scale="<<min_scale<<", max_scale="<<max_scale<<std::endl;
}
//#####################################################################
// Function Initialize_Sample_Locations_From_Vertices
//#####################################################################
template<class T> void SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<T>::
Initialize_Sample_Locations_From_Vertices(const TRIANGULATED_SURFACE<T>& triangulated_surface)
{
    const bool half_samples=false;
    if(half_samples){
        ARRAY<int> points;ARRAY<bool> done(triangulated_surface.mesh.elements.m,false);ARRAYS_COMPUTATIONS::Fill(done,false);
        TRIANGLE_MESH& mesh=triangulated_surface.mesh;GEOMETRY_PARTICLES<VECTOR<T,3> >& particles=triangulated_surface.particles;
        mesh.Initialize_Incident_Elements();
        for(int i=1;i<=triangulated_surface.particles.array_collection->Size();i++){
            if(done(i)) continue;else done(i)=true;
            T area=0;
            VECTOR<T,3> normal(0,0,0);
            ARRAY<int>& incident=(*triangulated_surface.mesh.incident_elements)(i);
            for(int j=1;j<=incident.m;j++){
                int node1,node2,node3;mesh.elements(incident(j)).Get(node1,node2,node3);
                TRIANGLE_3D<T> triangle(particles.X(node1),particles.X(node2),particles.X(node3));
                area+=triangle.Area();
                done(node1)=done(node2)=done(node3)=true;
                normal+=triangle.normal;}
            normal/=(T)incident.m;normal.Normalize();
            samples.Append(SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE<T>(particles.X(i),normal,area));}}
    else{
        //PHYSBAM_ASSERT(triangulated_surface.mesh.number_nodes==triangulated_surface.vertex_normals->m);
        triangulated_surface.mesh.Initialize_Incident_Elements();
        ARRAY<int> referenced;
        triangulated_surface.mesh.elements.Flattened().Get_Unique(referenced);
        samples.Resize(referenced.m);
        for(int ii=1;ii<=samples.m;++ii){
            int p=referenced(ii);
            T vertex_area=0;
            for(int j=1;j<=(*triangulated_surface.mesh.incident_elements)(p).m;++j){
                int node1,node2,node3;triangulated_surface.mesh.elements((*triangulated_surface.mesh.incident_elements)(p)(j)).Get(node1,node2,node3);
                vertex_area+=TRIANGLE_3D<T>::Area(triangulated_surface.particles.X(node1),triangulated_surface.particles.X(node2),triangulated_surface.particles.X(node3));}
            if((*triangulated_surface.mesh.incident_elements)(p).m<1) LOG::cout<<"sorry vertex "<<p<<" is broken"<<std::endl;
            vertex_area/=3;
            samples(ii).position=triangulated_surface.particles.X(p);samples(ii).area=vertex_area;
            samples(ii).normal=(*triangulated_surface.vertex_normals)(p);
            assert((*triangulated_surface.mesh.incident_elements)(p).m>0);}}
}
//#####################################################################
// Function Build_Octree_Cells
//#####################################################################
template<class T> void SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<T>::
Build_Octree_Cells(OCTREE_CELL<T>* cell)
{
    if(octree_cell_samples(cell->Cell()).m<=max_irradiance_leaf_samples) return;
    cell->Create_Children(octree_grid.number_of_cells,0,octree_grid.number_of_nodes,0,octree_grid.number_of_faces,0,&octree_grid);
    octree_cell_samples.Resize(octree_cell_samples.m+8);
    ARRAY<int>& samples_inside=octree_cell_samples(cell->Cell());
    for(int i=1;i<=samples_inside.m;i++){
        VECTOR<T,3> location=samples(samples_inside(i)).position;
        octree_cell_samples(octree_grid.Inside_Offspring(cell,location)->Cell()).Append(samples_inside(i));}
    for(int i=0;i<8;i++) Build_Octree_Cells(cell->Child(i));
}
//#####################################################################
// Function Build_Aggregate_Samples
//#####################################################################
template<class T> void SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<T>::
Build_Aggregate_Samples(OCTREE_CELL<T>* dummy_cell)
{
    for(int i=1;i<=aggregate_samples.m;++i){
         if(octree_cell_samples(i).m==0)continue;
         SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE<T> accumulated_sample;
         accumulated_sample.position=accumulated_sample.transmitted_irradiance=accumulated_sample.normal=accumulated_sample.transmitted_irradiance_product=VECTOR<T,3>(0,0,0);accumulated_sample.area=0;aggregate_samples(i)=accumulated_sample;
         T irradiance_magnitude_sum=0;
         for(int j=1;j<=octree_cell_samples(i).m;++j){
             int child_index=octree_cell_samples(i)(j);
             T irradiance_magnitude=samples(child_index).transmitted_irradiance.Magnitude();
             accumulated_sample.area+=samples(child_index).area;
             accumulated_sample.transmitted_irradiance+=samples(child_index).transmitted_irradiance;
             accumulated_sample.transmitted_irradiance_product+=samples(child_index).transmitted_irradiance_product;
             accumulated_sample.position+=samples(child_index).position*irradiance_magnitude;
             irradiance_magnitude_sum+=irradiance_magnitude;}
             if(irradiance_magnitude_sum>0){
                 accumulated_sample.position=accumulated_sample.position/irradiance_magnitude_sum;}
             else{
                 accumulated_sample.position=VECTOR<T,3>(0,0,0); 
                 for(int k=1;k<=octree_cell_samples(i).m;++k){accumulated_sample.position+=samples(octree_cell_samples(i)(k)).position;}
                 accumulated_sample.position/=(T)octree_cell_samples(i).m;}
             aggregate_samples(i)=accumulated_sample;}
}
//#####################################################################
// Function Build_Octree
//#####################################################################
template<class T> void SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<T>::
Build_Octree()
{
    LOG::cerr<<"Creating Octree, error_criterion="<<error_criterion<<" bounding_box="<<bounding_box<<std::endl;
    //bounding_box.Change_Size(bounding_box.Edge_Lengths().Max());
    RANGE<TV> pseudo_bounding_box(bounding_box.Minimum_Corner(),bounding_box.Maximum_Corner()+bounding_box.Edge_Lengths());
    octree_grid.Initialize(GRID<TV>(3,3,3,pseudo_bounding_box),0,0,false,false);root_cell=octree_grid.cells(1,1,1);
    octree_cell_samples.Resize(8);octree_cell_samples(1).Resize(samples.m);
    for(int i=1;i<=samples.m;i++) octree_cell_samples(1)(i)=i;
    LOG::cout<<"Num Irradiance Sample points: "<<samples.m<<std::endl;
    LOG::Time("Refine Octree");
    Build_Octree_Cells(root_cell);
    LOG::Time("Build aggregate irradiance samples"); 
    aggregate_samples.Resize(octree_grid.number_of_cells);
    Build_Aggregate_Samples(root_cell);
}
//#####################################################################
// Function Process_Irradiance_Candidates
//#####################################################################
template<class T> void SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<T>::
Process_Irradiance_Candidates(ARRAY<SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE<T> >& candidate_list,const VECTOR<T,3>& point,OCTREE_CELL<T>& cell,ARRAY<int>& sample_numbers)
{
    const SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE<T>& cell_sample=aggregate_samples(cell.Cell());
    if(octree_grid.Inside_Cell(&cell,point) || Approximate_Maximum_Solid_Angle(cell_sample.area,point,cell_sample.position)>error_criterion){
        if(!cell.Has_Children()){
            ARRAY<int>& candidates_in_cell=octree_cell_samples(cell.Cell());
            for(int i=1;i<=candidates_in_cell.m;i++){sample_numbers.Append(candidates_in_cell(i));candidate_list.Append(samples(candidates_in_cell(i)));}
            return;}
        else for(int i=0;i<8;++i){Process_Irradiance_Candidates(candidate_list,point,cell.children->children[i],sample_numbers);}}
    else{sample_numbers.Append(-cell.Cell());candidate_list.Append(cell_sample);}
}
//#####################################################################
// Function Find_Irradiance_Candidates
//#####################################################################
template<class T> void SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<T>::
Find_Irradiance_Candidates(ARRAY<SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE<T> >& candidate_list,const VECTOR<T,3>& point)
{
    ARRAY<int> sample_numbers;Process_Irradiance_Candidates(candidate_list,point,*root_cell,sample_numbers);
}
//#####################################################################
template class SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<float>;
template void SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<float>::Initialize_Sample_Locations_From_File<float>(const std::string& file,TRIANGULATED_SURFACE<float>&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<double>::Initialize_Sample_Locations_From_File<float>(const std::string& file,TRIANGULATED_SURFACE<double>&);
template void SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<double>::Initialize_Sample_Locations_From_File<double>(const std::string& file,TRIANGULATED_SURFACE<double>&);
template class SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<double>;
#endif
#endif

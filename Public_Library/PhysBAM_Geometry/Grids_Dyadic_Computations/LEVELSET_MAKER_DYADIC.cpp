//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/FLOOD_FILL_GRAPH.h>
#include <PhysBAM_Tools/Data_Structures/GRAPH.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/SEGMENT_3D_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Computations/LEVELSET_MAKER_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/EXTRAPOLATION_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/FAST_MARCHING_METHOD_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Function Coarsen_Tree
//#####################################################################
template<class T> static bool
Coarsen_Tree(OCTREE_CELL<T>* cell,ARRAY<T>& phi,const T refinement_distance)
{
    if(cell->Has_Children()){
        bool can_coarsen=true;
        for(int i=0;i<8;i++) can_coarsen&=Coarsen_Tree(cell->Child(i),phi,refinement_distance);
        if(can_coarsen){
            T cell_phi=0;for(int i=0;i<8;i++) cell_phi+=phi(cell->Child(i)->Cell());cell_phi*=(T).125;
            phi(cell->Cell())=cell_phi;cell->Delete_Children();}
        else return false;}
    T dx=(T)1.01*refinement_distance+(T).5*(T)root_three*cell->DX().x;
    return abs(phi(cell->Cell()))>dx;
}
//#####################################################################
// Function Compute_Level_Set
//#####################################################################
// Assumes triangulated surface is closed
template<class T> bool LEVELSET_MAKER_DYADIC<T>::
Compute_Level_Set(TRIANGULATED_SURFACE<T>& triangulated_surface,OCTREE_GRID<T>& grid,ARRAY<T>& phi,ARRAY<T>* phi_nodes,ARRAY<VECTOR<T,3> >* velocity,const T coarsen_bandwidth,int (*depth_function)(int,const RANGE<VECTOR<T,3> >&,void*),void* depth_data)
{
    if(remove_degenerate_triangles_area_threshold){
        if(verbose) LOG::Time("Remove Degenerate Triangles");
        triangulated_surface.Remove_Degenerate_Triangles(remove_degenerate_triangles_area_threshold);}

    // initialize acceleration structures
    if(verbose) LOG::Time("Initialize Acceleration Structures");
    bool triangle_list_defined=triangulated_surface.triangle_list!=0;if(!triangulated_surface.triangle_list) triangulated_surface.Update_Triangle_List();
    bool hierarchy_defined=triangulated_surface.hierarchy!=0;if(!triangulated_surface.hierarchy) triangulated_surface.Initialize_Hierarchy();
    bool bounding_box_defined=triangulated_surface.bounding_box!=0;if(!triangulated_surface.bounding_box) triangulated_surface.Update_Bounding_Box();
    bool adjacent_triangles_defined=triangulated_surface.mesh.adjacent_elements!=0;if(!adjacent_triangles_defined) triangulated_surface.mesh.Initialize_Adjacent_Elements();
    bool incident_triangles_defined=triangulated_surface.mesh.incident_elements!=0;if(!incident_triangles_defined) triangulated_surface.mesh.Initialize_Incident_Elements();
 
    bool compute_velocity=velocity && triangulated_surface.particles.store_velocity;

    if(use_fmm && compute_velocity && extrapolate_velocity && fmm_one_sided_band_width && fmm_one_sided_band_width<velocity_extrapolation_one_sided_band_width+1){
        LOG::cerr<<"Extending FMM band width to one more than velocity extrapolation band width"<<std::endl;
        fmm_one_sided_band_width=velocity_extrapolation_one_sided_band_width+1;}
    T fmm_stopping_distance=fmm_one_sided_band_width*grid.Minimum_Edge_Length();

    const T surface_thickness_over_two=Surface_Thickness_Over_Two(grid),surface_padding_for_flood_fill=Surface_Padding_For_Flood_Fill(grid);

    // build the octree topology
    if(verbose) LOG::Time("Refining Octree");
    ARRAY<OCTREE_CELL<T>*> refined_cells;
    for(int t=1;t<=triangulated_surface.mesh.elements.m;t++){
        const TRIANGLE_3D<T>& triangle=(*triangulated_surface.triangle_list)(t);
        TRIANGLE_3D<T> enlarged_triangle=triangle;enlarged_triangle.Change_Size(surface_padding_for_flood_fill);
        RANGE<VECTOR<T,3> > triangle_bounding_box=enlarged_triangle.Bounding_Box();
        triangle_bounding_box.Change_Size(surface_thickness_over_two);
        refined_cells.Remove_All();
        triangle_bounding_box.Change_Size((T).5*grid.Minimum_Edge_Length());
        if(depth_function) grid.Refine_Cells_Intersecting_Box(triangle_bounding_box,refined_cells,depth_function(t,triangle_bounding_box,depth_data));
        else grid.Refine_Cells_Intersecting_Box(triangle_bounding_box,refined_cells);}
    refined_cells.Clean_Memory();grid.Tree_Topology_Changed();grid.Node_Iterator_Data();

    if(compute_velocity) velocity->Resize(grid.number_of_cells);

    phi.Resize(grid.number_of_cells,false);phi_nodes->Resize(grid.number_of_nodes,false);ARRAYS_COMPUTATIONS::Fill(phi,(T)FLT_MAX);
 
    bool need_flood_fill=compute_signed_distance_function || compute_heaviside_function;
    if(need_flood_fill && verbose) LOG::Time("Building graph for flood fill");
    GRAPH graph(grid.number_of_cells);ARRAY<int> triangle_intersection_list;
    if(need_flood_fill) for(DYADIC_GRID_ITERATOR_FACE<OCTREE_GRID<T> > iterator(grid,grid.Map_All_Faces());iterator.Valid();iterator.Next()){
        OCTREE_CELL<T> *cell=iterator.Deepest_Cell(),*other_cell=iterator.Other_Cell();
        if(!other_cell) continue; // outer boundary face
        VECTOR<T,3> x1=cell->Center(),x2=other_cell->Center();
        SEGMENT_3D<T> segment=SEGMENT_3D<T>(x1,x2);
        triangulated_surface.hierarchy->Intersection_List(RANGE<VECTOR<T,3> >::Bounding_Box(x1,x2),triangle_intersection_list,surface_padding_for_flood_fill*2);
        bool add=true;
        for(int item=1;item<=triangle_intersection_list.m;item++){
            int t=triangle_intersection_list(item);
            const TRIANGLE_3D<T>& triangle=(*triangulated_surface.triangle_list)(t);
            TRIANGLE_3D<T> enlarged_triangle=triangle;enlarged_triangle.Change_Size(surface_padding_for_flood_fill);
        if(INTERSECTION::Intersects(segment,enlarged_triangle,surface_thickness_over_two)){add=false;break;}}
        if(add) graph.Add_Undirected_Edge(cell->Cell(),other_cell->Cell());
        triangle_intersection_list.Remove_All();}
    RANGE<VECTOR<T,3> > domain=grid.uniform_grid.domain;domain.Change_Size((T)1e-5);
    
    ARRAYS_COMPUTATIONS::Fill(graph.valid_nodes,false);
    for(DYADIC_GRID_ITERATOR_CELL<OCTREE_GRID<T> > iterator(grid,3);iterator.Valid();iterator.Next()) if(domain.Lazy_Inside(iterator.Location())) graph.valid_nodes(iterator.Cell_Index())=true;

    bool store_closest_triangle_index=need_flood_fill && !only_boundary_region_is_outside;
    ARRAY<int> closest_triangle_index;
    if(store_closest_triangle_index) closest_triangle_index.Resize(grid.number_of_cells);

    bool store_initialized_indices=use_fmm;
    if(store_initialized_indices){initialized_indices_octree.Exact_Resize(0);initialized_indices_octree.Preallocate(20);}

    if(verbose) LOG::Time("Rasterizing Triangles");
    ARRAY<OCTREE_CELL<T>*> intersecting_cells;
    BOX<VECTOR<T,3> > grid_domain(grid.Domain());
    for(int t=1;t<=triangulated_surface.mesh.elements.m;t++){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(t%1000==0 && verbose) LOG::cout<<"Doing triangle "<<t<<" of "<<triangulated_surface.mesh.elements.m<<std::endl;
#endif
        const TRIANGLE_3D<T>& triangle=(*triangulated_surface.triangle_list)(t);
        TRIANGLE_3D<T> enlarged_triangle=triangle;enlarged_triangle.Change_Size(surface_padding_for_flood_fill);
        RANGE<VECTOR<T,3> > triangle_bounding_box=enlarged_triangle.Bounding_Box();
        if(!grid_domain.Lazy_Intersection(triangle_bounding_box)) continue;
        triangle_bounding_box.Change_Size(surface_thickness_over_two+(T).5*grid.Minimum_Edge_Length());
        intersecting_cells.Remove_All();
        grid.Get_Cells_Intersecting_Box(triangle_bounding_box,intersecting_cells);

        for(int k=1;k<=intersecting_cells.m;k++){
            int i=intersecting_cells(k)->Cell();
            VECTOR<T,3> grid_position=intersecting_cells(k)->Center(),weights,closest_point=triangle.Closest_Point(grid_position,weights);
            T distance_squared=(grid_position-closest_point).Magnitude_Squared();
            if(phi(i)==FLT_MAX || distance_squared<sqr(phi(i))){
                if(store_initialized_indices && phi(i)==FLT_MAX)initialized_indices_octree.Append(i);
                phi(i)=sqrt(distance_squared);
                if(store_closest_triangle_index) closest_triangle_index(i)=t;
                if(compute_velocity){
                    int node1,node2,node3;triangulated_surface.mesh.elements(t).Get(node1,node2,node3);
                    (*velocity)(i)=weights.x*triangulated_surface.particles.V(node1)+weights.y*triangulated_surface.particles.V(node2)+weights.z*triangulated_surface.particles.V(node3);}}}}
     intersecting_cells.Clean_Memory();

    if((compute_signed_distance_function || compute_unsigned_distance_function) && use_fmm && fmm_stopping_distance)
        for(int i=1;i<=grid.number_of_cells;i++) phi(i)=min(phi(i),fmm_stopping_distance);
    else if(compute_heaviside_function) ARRAYS_COMPUTATIONS::Fill(phi,grid.Minimum_Edge_Length());

    if(use_orthogonal_vote) PHYSBAM_NOT_IMPLEMENTED(); // not supported for now
    else if(need_flood_fill){ // Need flood fill to determine sign (inside/outside)
        if(verbose) LOG::Time("Flood Fill");
        ARRAY<int> colors(grid.number_of_cells);FLOOD_FILL_GRAPH flood_fill;//flood_fill.Optimize_Fill_For_Single_Cell_Regions(true);
        int number_of_colors=flood_fill.Flood_Fill(graph,colors);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(verbose) LOG::cout<<"(got "<<number_of_colors<<" colors)... "<<std::endl;
#endif
        if(number_of_colors==1){ // there is only one color. check if the whole domain is inside or outside then return
            if(triangulated_surface.Inside(grid.uniform_grid.domain.min_corner)) 
                ARRAYS_COMPUTATIONS::Fill(phi,(T)-FLT_MAX);else ARRAYS_COMPUTATIONS::Fill(phi,(T)FLT_MAX);
            return false;}
        
        if(verbose) LOG::Time("Classifying colors");
        ARRAY<bool,VECTOR<int,1> > color_is_inside(1,number_of_colors);
        ARRAY<T> color_maximum_distance(number_of_colors,false);ARRAYS_COMPUTATIONS::Fill(color_maximum_distance,(T)-1);
        ARRAY<OCTREE_CELL<T>*> color_representatives(number_of_colors);
        for(DYADIC_GRID_ITERATOR_CELL<OCTREE_GRID<T> > iterator(grid,3);iterator.Valid();iterator.Next()){int cell_index=iterator.Cell_Index();
            if(graph.valid_nodes(cell_index) && closest_triangle_index(cell_index) && color_maximum_distance(colors(cell_index))<phi(cell_index)){
                color_maximum_distance(colors(cell_index))=phi(cell_index);color_representatives(colors(cell_index))=iterator.Cell_Pointer();}}
        for(int color=1;color<=number_of_colors;color++){
            if(color_maximum_distance(color)<0){LOG::cerr<<"Error: could not determine inside/outside for color "<<color<<std::endl;}
            else color_is_inside(color)=triangulated_surface.Inside_Relative_To_Triangle(color_representatives(color)->Center(),
                closest_triangle_index(color_representatives(color)->Cell()),surface_thickness_over_two);}
        for(DYADIC_GRID_ITERATOR_CELL<OCTREE_GRID<T> > iterator(grid,3);iterator.Valid();iterator.Next()){
            int cell_index=iterator.Cell_Index();
            if(graph.valid_nodes(cell_index) && color_is_inside(colors(cell_index))) phi(cell_index)*=-1;}}

    if(positive_boundary_band){
        LOG::Time(STRING_UTILITIES::string_sprintf("Enforcing positive boundary of %d",positive_boundary_band));
        T clip=((T)positive_boundary_band+(T).5)*grid.Minimum_Edge_Length();
        for(DYADIC_GRID_ITERATOR_CELL<OCTREE_GRID<T> > iterator(grid,3);iterator.Valid();iterator.Next()){
            if(grid_domain.Signed_Distance(iterator.Location()) >= -clip) phi(iterator.Cell_Index())=max(phi(iterator.Cell_Index()),grid.Minimum_Edge_Length());}}

    for(DYADIC_GRID_ITERATOR_CELL<OCTREE_GRID<T> > iterator(grid,3);iterator.Valid();iterator.Next())
        phi(iterator.Cell_Index())=clamp(phi(iterator.Cell_Index()),-10*grid.Minimum_Edge_Length(),10*grid.Minimum_Edge_Length()); // clamp away from FLT_MAX to avoid
                                                                                                                                           // floating point exceptions
    BOUNDARY_DYADIC<OCTREE_GRID<T>,T>().Fill_Ghost_Cells_Cell(grid,phi,phi,0);
    if(use_fmm && (compute_unsigned_distance_function || compute_signed_distance_function)){
        LEVELSET_OCTREE<T> levelset(grid,phi);
        if(phi_nodes){
            if(verbose) LOG::Time(STRING_UTILITIES::string_sprintf("Fast Marching Nodes (one sided band width=%f)",fmm_one_sided_band_width));
            LINEAR_INTERPOLATION_DYADIC_HELPER<OCTREE_GRID<T> >::Interpolate_From_Cells_To_Nodes(grid,phi,*phi_nodes);
            FAST_MARCHING_METHOD_DYADIC<OCTREE_GRID<T> > fmm(levelset);
            if(compute_unsigned_distance_function) fmm.Fast_Marching_Method_Nodes(*phi_nodes,fmm_one_sided_band_width,&initialized_indices_octree,false,&phi);
            else fmm.Fast_Marching_Method_Nodes(*phi_nodes,fmm_one_sided_band_width,0,false,&phi);
            LINEAR_INTERPOLATION_DYADIC_HELPER<OCTREE_GRID<T> >::Interpolate_From_Nodes_To_Cells(grid,*phi_nodes,phi);}
        else{
            if(verbose) LOG::Time(STRING_UTILITIES::string_sprintf("Fast Marching Cells (one sided band width=%f)",fmm_one_sided_band_width));
            if(compute_unsigned_distance_function) levelset.Fast_Marching_Method(0,fmm_stopping_distance,&initialized_indices_octree);
            else if(compute_signed_distance_function) levelset.Fast_Marching_Method(0,fmm_stopping_distance);}}

    if(coarsen_bandwidth){
        if(verbose) LOG::Time("Coarsening");
        T refinement_distance=(T).5*coarsen_bandwidth*grid.Minimum_Edge_Length();
        int old_number_of_cells=grid.number_of_cells,old_number_of_nodes=grid.number_of_nodes,old_number_of_faces=grid.number_of_faces;
        for(int i=1;i<=grid.uniform_grid.counts.x;i++) for(int j=1;j<=grid.uniform_grid.counts.y;j++) for(int ij=1;ij<=grid.uniform_grid.counts.z;ij++) 
            Coarsen_Tree(grid.cells(i,j,ij),phi,refinement_distance);

        ARRAY<int> cell_mapping_array,node_mapping_array,face_mapping_array;
        grid.Compact_Array_Indices(&cell_mapping_array,&node_mapping_array,&face_mapping_array);
        ARRAY<T> temp(max(old_number_of_cells,old_number_of_nodes,old_number_of_faces),false);  
        ARRAY<T>::Compact_Array_Using_Compaction_Array(phi,cell_mapping_array,&temp);
        phi.Resize(grid.number_of_cells);
        grid.Tree_Topology_Changed();}

    if(compute_velocity && extrapolate_velocity){
        if(!compute_signed_distance_function || !use_fmm){LOG::cerr<<"Can only extrapolate velocity if computing signed distance function"<<std::endl;}
        else{
            if(verbose) LOG::Time(STRING_UTILITIES::string_sprintf("Extrapolating velocity (one sided band width=%f)",velocity_extrapolation_one_sided_band_width));
            EXTRAPOLATION_DYADIC<OCTREE_GRID<T>,VECTOR<T,3> > extrapolation(grid);
            extrapolation.Set_Custom_Seed_Indices(&initialized_indices_octree);
            extrapolation.Set_Band_Width(velocity_extrapolation_one_sided_band_width);
            extrapolation.Extrapolate_Cells(phi,*velocity,false);
            phi*=-1;
            extrapolation.Extrapolate_Cells(phi,*velocity,false);
            phi*=-1;}}

    if(verbose) LOG::Stop_Time();

    // delete acceleration structures if defined in this function
    if(!incident_triangles_defined){delete triangulated_surface.mesh.incident_elements;triangulated_surface.mesh.incident_elements=0;}
    if(!adjacent_triangles_defined){delete triangulated_surface.mesh.adjacent_elements;triangulated_surface.mesh.adjacent_elements=0;}
    if(!bounding_box_defined){delete triangulated_surface.bounding_box;triangulated_surface.bounding_box=0;}
    if(!hierarchy_defined){delete triangulated_surface.hierarchy;triangulated_surface.hierarchy=0;}
    if(!triangle_list_defined){delete triangulated_surface.triangle_list;triangulated_surface.triangle_list=0;}
    return true;
}
//#####################################################################
template class LEVELSET_MAKER_DYADIC<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_MAKER_DYADIC<double>;
#endif
#endif

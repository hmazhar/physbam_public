//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Data_Structures/FLOOD_FILL_GRAPH.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FLOOD_FILL_3D.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/SEGMENT_3D_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Function Compute_Level_Set
//#####################################################################
// Assumes triangulated surface is closed
template<class T> bool LEVELSET_MAKER_UNIFORM<T>::
Compute_Level_Set(TRIANGULATED_SURFACE<T>& triangulated_surface,GRID<TV>& grid,ARRAY<T,TV_INT>& phi,ARRAY<TV,TV_INT>* velocity)
{
    PHYSBAM_ASSERT(grid.counts.x && grid.counts.y && grid.counts.z,"Attempting to rasterize surface onto grid with no indices");
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
    if(compute_velocity) velocity->Resize(grid.Domain_Indices());

    phi.Resize(grid.Domain_Indices(),false);phi.Fill(FLT_MAX);
   
    if(use_fmm && compute_velocity && extrapolate_velocity && fmm_one_sided_band_width && fmm_one_sided_band_width<velocity_extrapolation_one_sided_band_width+1){
        LOG::cerr<<"Extending FMM band width to one more than velocity extrapolation band width"<<std::endl;
        fmm_one_sided_band_width=velocity_extrapolation_one_sided_band_width+1;}
    T fmm_stopping_distance=fmm_one_sided_band_width*grid.dX.Max();

    bool need_flood_fill=compute_signed_distance_function || compute_heaviside_function;
    ARRAY<bool,TV_INT> edge_is_blocked_x,edge_is_blocked_y,edge_is_blocked_z;
    if(need_flood_fill){
        edge_is_blocked_x.Resize(2,grid.counts.x,1,grid.counts.y,1,grid.counts.z);
        edge_is_blocked_y.Resize(1,grid.counts.x,2,grid.counts.y,1,grid.counts.z);
        edge_is_blocked_z.Resize(1,grid.counts.x,1,grid.counts.y,2,grid.counts.z);}

    bool store_closest_triangle_index=need_flood_fill && !only_boundary_region_is_outside;
    ARRAY<int,TV_INT> closest_triangle_index;
    if(store_closest_triangle_index) closest_triangle_index.Resize(grid.Domain_Indices());

    bool store_initialized_indices=use_fmm;
    if(store_initialized_indices){initialized_indices.Exact_Resize(0);initialized_indices.Preallocate(20);}

    const T surface_thickness_over_two=Surface_Thickness_Over_Two(grid),surface_padding_for_flood_fill=Surface_Padding_For_Flood_Fill(grid);

    const BOX<TV>& grid_domain=grid.domain;
    if(verbose) LOG::Time("Rasterizing Triangles");
    for(int t=1;t<=triangulated_surface.mesh.elements.m;t++){
        const TRIANGLE_3D<T>& triangle=(*triangulated_surface.triangle_list)(t);
        TRIANGLE_3D<T> enlarged_triangle=triangle;if(surface_padding_for_flood_fill) enlarged_triangle.Change_Size(surface_padding_for_flood_fill);
        RANGE<TV> triangle_bounding_box=enlarged_triangle.Bounding_Box();
        triangle_bounding_box.Change_Size(surface_thickness_over_two);
        if(!grid_domain.Lazy_Intersection(triangle_bounding_box)) continue;
        TV_INT min_index=grid.Clamped_Index(triangle_bounding_box.Minimum_Corner()),max_index=grid.Clamped_Index_End_Minus_One(triangle_bounding_box.Maximum_Corner())+TV_INT(1,1,1);
        for(int i=min_index.x;i<=max_index.x;i++) for(int j=min_index.y;j<=max_index.y;j++) for(int k=min_index.z;k<=max_index.z;k++){
            TV grid_position=grid.X(i,j,k),weights,closest_point=triangle.Closest_Point(grid_position,weights);
            T distance_squared=(grid_position-closest_point).Magnitude_Squared();
            if(phi(i,j,k)==FLT_MAX || distance_squared<sqr(phi(i,j,k))){
                if(store_initialized_indices && phi(i,j,k)==FLT_MAX)initialized_indices.Append(TV_INT(i,j,k));
                phi(i,j,k)=sqrt(distance_squared);
                if(store_closest_triangle_index) closest_triangle_index(i,j,k)=t;
                if(compute_velocity){
                    int node1,node2,node3;triangulated_surface.mesh.elements(t).Get(node1,node2,node3);
                    (*velocity)(i,j,k)=weights.x*triangulated_surface.particles.V(node1)+weights.y*triangulated_surface.particles.V(node2)+weights.z*triangulated_surface.particles.V(node3);}}}
        if(need_flood_fill){
            for(int i=min_index.x+1;i<=max_index.x;i++) for(int j=min_index.y;j<=max_index.y;j++) for(int k=min_index.z;k<=max_index.z;k++)
                if(!edge_is_blocked_x(i,j,k)) edge_is_blocked_x(i,j,k)=INTERSECTION::Intersects(SEGMENT_3D<T>(grid.X(i,j,k),grid.X(i-1,j,k)),enlarged_triangle,surface_thickness_over_two);
            for(int i=min_index.x;i<=max_index.x;i++) for(int j=min_index.y+1;j<=max_index.y;j++) for(int k=min_index.z;k<=max_index.z;k++)
                if(!edge_is_blocked_y(i,j,k)) edge_is_blocked_y(i,j,k)=INTERSECTION::Intersects(SEGMENT_3D<T>(grid.X(i,j,k),grid.X(i,j-1,k)),enlarged_triangle,surface_thickness_over_two);
            for(int i=min_index.x;i<=max_index.x;i++) for(int j=min_index.y;j<=max_index.y;j++) for(int k=min_index.z+1;k<=max_index.z;k++)
                if(!edge_is_blocked_z(i,j,k)) edge_is_blocked_z(i,j,k)=INTERSECTION::Intersects(SEGMENT_3D<T>(grid.X(i,j,k),grid.X(i,j,k-1)),enlarged_triangle,surface_thickness_over_two);}}

    if((compute_signed_distance_function || compute_unsigned_distance_function) && use_fmm && fmm_stopping_distance)
        for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) for(int k=1;k<=grid.counts.z;k++) phi(i,j,k)=min(phi(i,j,k),fmm_stopping_distance);
    else if(compute_heaviside_function) phi.Fill(grid.dX.Max());

    if(write_debug_data){
        GRID<TV> output_grid=grid;
        if(!grid.Is_MAC_Grid())
            output_grid=GRID<TV>(grid.counts,RANGE<TV>(grid.Xmin()-grid.dX/2,grid.Xmax()+grid.dX/2),true);
        // TODO: put this back if you need it
        /*FILE_UTILITIES::Write_To_File<T>("grid.debug",output_grid);
          FILE_UTILITIES::Write_To_File<T>("triangulated_surface.debug",triangulated_surface);*/}

    if(use_orthogonal_vote){
        // TODO: put this back if you need it
        /*if(write_debug_data){FILE_UTILITIES::Write_To_File<T>("edge_is_blocked.debug",edge_is_blocked_x,edge_is_blocked_y,edge_is_blocked_z);}*/
        if(verbose) LOG::Time("Computing sign using orthogonal vote");
        ARRAY<bool,TV_INT> is_inside(grid.Domain_Indices());
        for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) for(int k=1;k<=grid.counts.z;k++){
            if(closest_triangle_index(i,j,k)) is_inside(i,j,k)=triangulated_surface.Inside_Relative_To_Triangle(grid.X(i,j,k),closest_triangle_index(i,j,k),surface_thickness_over_two);}
        ARRAY<char,TV_INT> vote(grid.Domain_Indices());
        for(int j=1;j<=grid.counts.y;j++) for(int k=1;k<=grid.counts.z;k++) Process_Segment(grid.counts.x,edge_is_blocked_x,is_inside,TV_INT(1,j,k),1,vote);
        for(int i=1;i<=grid.counts.x;i++) for(int k=1;k<=grid.counts.z;k++) Process_Segment(grid.counts.y,edge_is_blocked_y,is_inside,TV_INT(i,1,k),2,vote);
        for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) Process_Segment(grid.counts.z,edge_is_blocked_z,is_inside,TV_INT(i,j,1),3,vote);
        is_inside.Clean_Memory();
        for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) for(int k=1;k<=grid.counts.z;k++) if(vote(i,j,k)>=3) phi(i,j,k)*=-1;

        if(keep_only_largest_inside_region){
            ARRAY<int,TV_INT> colors(grid.Domain_Indices());ARRAY<bool,TV_INT> null_edge_is_blocked(grid.Domain_Indices(1));
            for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) for(int k=1;k<=grid.counts.z;k++) if(phi(i,j,k)>0) colors(i,j,k)=-1; // make outside regions uncolorable
            FLOOD_FILL_3D flood_fill;flood_fill.Optimize_Fill_For_Single_Cell_Regions(true);
            int number_of_colors=flood_fill.Flood_Fill(colors,null_edge_is_blocked,null_edge_is_blocked,null_edge_is_blocked);
            // TODO: put this back if you need it
            //if(write_debug_data){FILE_UTILITIES::Write_To_File<T>("colors.debug",colors);}
            ARRAY<int> region_size(number_of_colors);
            for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) for(int k=1;k<=grid.counts.z;k++) if(colors(i,j,k)>0) region_size(colors(i,j,k))++;
            int max_region_size=ARRAYS_COMPUTATIONS::Max(region_size);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(verbose) LOG::cout<<"Keeping only largest inside region (max region size = "<<max_region_size<<")... "<<std::flush;
#endif
            // flip smaller regions back to positive sign
            for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) for(int k=1;k<=grid.counts.z;k++) if(colors(i,j,k)>0 && region_size(colors(i,j,k))<max_region_size) phi(i,j,k)*=-1;}
    }
    else if(need_flood_fill){ // Need flood fill to determine sign (inside/outside)
        if(verbose) LOG::Time("Flood Fill");
        ARRAY<int,TV_INT> colors(grid.Domain_Indices());FLOOD_FILL_3D flood_fill;flood_fill.Optimize_Fill_For_Single_Cell_Regions(true);
        int number_of_colors=flood_fill.Flood_Fill(colors,edge_is_blocked_x,edge_is_blocked_y,edge_is_blocked_z);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(verbose) LOG::cout<<"(got "<<number_of_colors<<" colors)... "<<std::endl;
#endif
        if(number_of_colors==1 && !phi_offset){ // there is only one color. check if the whole domain is inside or outside then return
            if(triangulated_surface.Inside(grid.X(1,1,1)))phi.Fill(-FLT_MAX);else phi.Fill(FLT_MAX);return false;}
        // TODO: put this back if you need it
        /*if(write_debug_data){
            FILE_UTILITIES::Write_To_File<T>("colors.debug",colors);
            FILE_UTILITIES::Write_To_File<T>("edge_is_blocked.debug",edge_is_blocked_x,edge_is_blocked_y,edge_is_blocked_z);}
        if(write_debug_path){
            ARRAY<TV_INT> path_nodes;
            bool path_exists=flood_fill.Path_Between_Nodes(RANGE<TV_INT>(1,grid.counts.x,1,grid.counts.y,1,grid.counts.z),path_start_node,path_end_node,
                                                           edge_is_blocked_x,edge_is_blocked_y,edge_is_blocked_z,&path_nodes);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(verbose){LOG::cout<<"Path between "<<path_start_node<<" and "<<path_end_node<<" "<<(path_exists?"exists":"doesn't exist")<<std::endl;}
#endif
            ARRAY<bool,TV_INT> path(grid.Domain_Indices());
            for(int i=1;i<=path_nodes.m;i++){path(path_nodes(i))=true;}
            FILE_UTILITIES::Write_To_File<T>("path.debug",path);}*/
        ARRAY<bool,VECTOR<int,1> > color_is_inside(1,number_of_colors);
        if(only_boundary_region_is_outside){
            ARRAY<bool> color_touches_boundary(number_of_colors);
            if(verbose) LOG::Time("Marking boundary region as outside");
            flood_fill.Identify_Colors_Touching_Boundary(number_of_colors,colors,edge_is_blocked_x,edge_is_blocked_y,edge_is_blocked_z,color_touches_boundary);
            if(verbose && color_touches_boundary.Number_True()>1) LOG::cerr<<"Warning: Got "<<color_touches_boundary.Number_True()<<" colors touching boundary"<<std::endl;
            for(int i=1;i<=number_of_colors;i++) color_is_inside(i)=!color_touches_boundary(i);}
        else{
            ARRAY<T> color_maximum_distance(number_of_colors,false);ARRAYS_COMPUTATIONS::Fill(color_maximum_distance,(T)-1);
            ARRAY<TV_INT> color_representatives(number_of_colors);
            for(int i=1;i<=grid.counts.x;i++)for(int j=1;j<=grid.counts.y;j++)for(int k=1;k<=grid.counts.z;k++)if(closest_triangle_index(i,j,k) && color_maximum_distance(colors(i,j,k))<phi(i,j,k)){
                color_maximum_distance(colors(i,j,k))=phi(i,j,k);color_representatives(colors(i,j,k))=TV_INT(i,j,k);}
            for(int color=1;color<=number_of_colors;color++){
                if(color_maximum_distance(color)<0){LOG::cerr<<"Error: could not determine inside/outside for color "<<color<<std::endl;PHYSBAM_FATAL_ERROR();}
                else color_is_inside(color)=triangulated_surface.Inside_Relative_To_Triangle(grid.X(color_representatives(color)),
                    closest_triangle_index(color_representatives(color)),surface_thickness_over_two);}}
        if(keep_only_largest_inside_region){
            ARRAY<int> region_size(number_of_colors);
            for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) for(int k=1;k<=grid.counts.z;k++) if(colors(i,j,k)>0 && color_is_inside(colors(i,j,k))) region_size(colors(i,j,k))++;
            int max_region_size=ARRAYS_COMPUTATIONS::Max(region_size);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(verbose) LOG::cout<<"Keeping only largest inside region (max region size = "<<max_region_size<<")... "<<std::flush;
#endif
            for(int i=1;i<=number_of_colors;i++) if(color_is_inside(i) && region_size(i)<max_region_size) color_is_inside(i)=false;}
        if(flip_sign_if_corners_are_inside){ // If the majority of corners are labelled as inside then we flip signs
            int num_corners_inside=(int)color_is_inside(colors(1,1,1))+(int)color_is_inside(colors(1,1,grid.counts.z))+
                                   (int)color_is_inside(colors(1,grid.counts.y,1))+(int)color_is_inside(colors(1,grid.counts.y,grid.counts.z))+
                                   (int)color_is_inside(colors(grid.counts.x,1,1))+(int)color_is_inside(colors(grid.counts.x,1,grid.counts.z))+
                                   (int)color_is_inside(colors(grid.counts.x,grid.counts.y,1))+(int)color_is_inside(colors(grid.counts.x,grid.counts.y,grid.counts.z));
            if(num_corners_inside>4){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
                if(verbose) LOG::cout<<"Majority of corners are inside -- flipping sign!"<<std::endl;
#endif
                for(int i=1;i<=number_of_colors;i++) color_is_inside(i)=!color_is_inside(i);}}
        for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) for(int k=1;k<=grid.counts.z;k++) if(color_is_inside(colors(i,j,k))) phi(i,j,k)*=-1;}

    if(positive_boundary_band){
        RANGE<TV> clip(grid.X_plus_half(TV_INT()+positive_boundary_band),grid.X_plus_half(grid.counts-positive_boundary_band));
        for(int j=1;j<=grid.counts.y;j++) for(int k=1;k<=grid.counts.z;k++){
            for(int i=1;i<=positive_boundary_band+1;i++) phi(i,j,k)=max(phi(i,j,k),clip.min_corner.x-grid.Axis_X(i,1));
            for(int i=grid.counts.x-positive_boundary_band;i<=grid.counts.x;i++) phi(i,j,k)=max(phi(i,j,k),grid.Axis_X(i,1)-clip.max_corner.x);}

        for(int i=1;i<=grid.counts.x;i++) for(int k=1;k<=grid.counts.z;k++){
            for(int j=1;j<=positive_boundary_band+1;j++) phi(i,j,k)=max(phi(i,j,k),clip.min_corner.y-grid.Axis_X(j,2));
            for(int j=grid.counts.y-positive_boundary_band;j<=grid.counts.y;j++) phi(i,j,k)=max(phi(i,j,k),grid.Axis_X(j,2)-clip.max_corner.y);}

        for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++){
            for(int k=1;k<=positive_boundary_band+1;k++) phi(i,j,k)=max(phi(i,j,k),clip.min_corner.z-grid.Axis_X(k,3));
            for(int k=grid.counts.z-positive_boundary_band;k<=grid.counts.z;k++) phi(i,j,k)=max(phi(i,j,k),grid.Axis_X(k,3)-clip.max_corner.z);}}

    if(use_fmm && (compute_unsigned_distance_function || compute_signed_distance_function)){
        for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) for(int k=1;k<=grid.counts.z;k++) 
            phi(i,j,k)=clamp(phi(i,j,k),-10*grid.min_dX,10*grid.min_dX); // clamp away from FLT_MAX to avoid floating point exceptions
        if(verbose) LOG::Time(STRING_UTILITIES::string_sprintf("Fast Marching (one sided band width=%f)",fmm_one_sided_band_width));
        GRID<TV> grid_copy=grid;LEVELSET_3D<GRID<TV> > levelset(grid_copy,phi);
        if(compute_unsigned_distance_function) levelset.Fast_Marching_Method(0,fmm_stopping_distance,&initialized_indices);
        else if(compute_signed_distance_function) levelset.Fast_Marching_Method(0,fmm_stopping_distance,phi_offset?&initialized_indices:0);}

    if(compute_velocity && extrapolate_velocity){
        if(!compute_signed_distance_function || !use_fmm){LOG::cerr<<"Can only extrapolate velocity if computing signed distance function"<<std::endl;}
        else{
            if(verbose) LOG::Time(STRING_UTILITIES::string_sprintf("Extrapolating velocity (one sided band width=%f)",velocity_extrapolation_one_sided_band_width));
            GRID<TV> grid_copy=grid;
            EXTRAPOLATION_UNIFORM<GRID<TV>,TV> extrapolation(grid_copy,phi,*velocity,3);
            extrapolation.Set_Custom_Seed_Indices(&initialized_indices);
            extrapolation.Set_Band_Width(velocity_extrapolation_one_sided_band_width);
            extrapolation.Extrapolate(0);
            phi*=-1;
            extrapolation.Extrapolate(0);
            phi*=-1;}}

    if(phi_offset){
        phi-=phi_offset;
        if(use_fmm && compute_signed_distance_function)
            LEVELSET_3D<GRID<TV> >(grid,phi).Fast_Marching_Method(0,fmm_stopping_distance);}

    // TODO: put this back if you need it
    /*if(write_debug_data){
        GRID<TV> grid_copy=grid;LEVELSET_3D<GRID<TV> > levelset(grid_copy,phi);
        FILE_UTILITIES::Write_To_File<T>("levelset.debug",levelset);
        if(compute_velocity)FILE_UTILITIES::Write_To_File<T>("velocity.debug",velocity);}*/

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(verbose) LOG::cout<<"Done"<<std::endl;
#endif

    // delete acceleration structures if defined in this function
    if(!incident_triangles_defined){delete triangulated_surface.mesh.incident_elements;triangulated_surface.mesh.incident_elements=0;}
    if(!adjacent_triangles_defined){delete triangulated_surface.mesh.adjacent_elements;triangulated_surface.mesh.adjacent_elements=0;}
    if(!bounding_box_defined){delete triangulated_surface.bounding_box;triangulated_surface.bounding_box=0;}
    if(!hierarchy_defined){delete triangulated_surface.hierarchy;triangulated_surface.hierarchy=0;}
    if(!triangle_list_defined){delete triangulated_surface.triangle_list;triangulated_surface.triangle_list=0;}
    return true;
}
//#####################################################################
// Function Process_Segment
//#####################################################################
static void Process_Segment(const int m,const ARRAY<bool,VECTOR<int,3> >& edge_is_blocked,const ARRAY<bool,VECTOR<int,3> >& is_inside,const VECTOR<int,3>& start_index,const int axis,ARRAY<char,VECTOR<int,3> >& vote)
{
    typedef VECTOR<int,3> TV_INT;
    TV_INT index=start_index,increment;increment[axis]=1;
    int segment_start=1;bool segment_starts_inside=false;
    for(index[axis]=1;index[axis]<=m;index[axis]++){
        if(index[axis]>1 && edge_is_blocked(index)){segment_start=index[axis];segment_starts_inside=is_inside(index);}
        if(index[axis]==m || (index[axis]<m && edge_is_blocked(index+increment))){
            int segment_end=index[axis];bool segment_ends_inside=index[axis]<m?is_inside(index):false;
            int vote_increment=(int)segment_starts_inside+(int)segment_ends_inside;
            TV_INT t=start_index;
            for(t[axis]=segment_start;t[axis]<=segment_end;t[axis]++) vote(t)+=vote_increment;}}
}
//#####################################################################
template class LEVELSET_MAKER_UNIFORM<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_MAKER_UNIFORM<double>;
#endif

//#####################################################################
// Copyright 2004-2005, Jiayi Chong, Geoffrey Irving, Igor Neverov, Andy Selle, Michael Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POINT_REPULSION
//##################################################################### 
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Collisions/POINT_REPULSION.h>
using namespace PhysBAM; 
//##########################################################################################
// Function Initialize
//##########################################################################################
template<class T> void POINT_REPULSION<T>::
Initialize(int point_number)
{
    surface.Use_Face_Normals();
    surface.Update_Vertex_Normals();
    surface.Update_Triangle_List();
    surface.mesh.Initialize_Adjacent_Elements();
    Initialize_Radii_Of_Curvature();
    surface_area=surface.Total_Area();
    if(surface.bounding_box==0)surface.Update_Bounding_Box();
    surface_linear_size=surface.bounding_box->Edge_Lengths().Magnitude();
    overpopulation_fraction=(T)1.65;
    Set_Standard_Number_Points(point_number);
    Set_Number_Points(point_number);
    Initialize_Sampling();
    LOG::cout<<"Initializing Vicinity Triangles"<<std::endl;Initialize_Vicinity_Triangles();
    LOG::cout<<"Initializing Segment Mesh"<<std::endl;surface.mesh.Initialize_Segment_Mesh();
}
//##########################################################################################
// Function Initialize_Vicinity_Triangles
//##########################################################################################
template<class T> void POINT_REPULSION<T>::
Initialize_Vicinity_Triangles()
{
    if(surface.mesh.adjacent_elements==0) PHYSBAM_FATAL_ERROR("Need adjacent elements to initialize vinicity_triangles");
    const ARRAY<ARRAY<int> >& adjacent_triangles=*surface.mesh.adjacent_elements;
    vicinity_triangles.Resize(surface.mesh.elements.m); 
    ARRAY<int> candidates,next_candidates;
    average_number_vicinity_triangles=0;int progress=0,new_progress=0;
    LOG::cout<<"Vicinity triangle list size: "<<vicinity_triangles.m<<std::endl;
    for(int t=1;t<=vicinity_triangles.m;++t){
        ARRAY<int>& vicinity_triangles_t=vicinity_triangles(t);vicinity_triangles_t.Resize(0);
        vicinity_triangles_t.Append(t);
        candidates=adjacent_triangles(t);T r=Repulsion_Radius(t);int cycle=1;
        while(candidates.m){
            for(int i=1;i<=candidates.m;++i){
                int candidate=candidates(i);
                if(cycle==1 || Acceptable_Vicinity_Triangle(t,r,candidate)){
                    vicinity_triangles_t.Append_Unique(candidate);
                    const ARRAY<int>& adjacent_to_candidate=adjacent_triangles(candidate);
                    for(int j=1;j<=adjacent_to_candidate.m;++j)     
                        if(vicinity_triangles_t.Find(adjacent_to_candidate(j))==0)next_candidates.Append_Unique(adjacent_to_candidate(j));}}
            ++cycle;candidates=next_candidates;next_candidates.Resize(0);}
        average_number_vicinity_triangles+=vicinity_triangles_t.m;
        new_progress=(int)100.0*t/vicinity_triangles.m;if(new_progress>progress){progress=new_progress;LOG::cout<<progress<<"%"<<" ";}}
    LOG::cout<<std::endl;average_number_vicinity_triangles/=vicinity_triangles.m;
}
//##########################################################################################
// Function Update_Neighbor_Points
//##########################################################################################
template<class T> void POINT_REPULSION<T>::
Update_Neighbor_Points()
{
    if(vicinity_triangles.m!=surface.mesh.elements.m) PHYSBAM_FATAL_ERROR("Need vicinity_triangles to update neighbor points");
    for(int i=1;i<=points.m;++i){
        POINT_REPULSION_DATA<T>& data=points(i);VECTOR<T,3> position=data.position,normal=surface.Normal(position,data.triangle); 
        T radius=sqr(Repulsion_Radius(data.triangle)); 
        data.neighbors.Resize(0);const ARRAY<int>& candidate_triangles=vicinity_triangles(data.triangle);
        for(int t=1;t<=candidate_triangles.m;++t){
            const ARRAY<int>& candidate_points=points_in_triangle(candidate_triangles(t));
            for(int j=1;j<=candidate_points.m;++j){
                int candidate=candidate_points(j); 
                if(i!=candidate && Acceptable_Neighbor_Point(position,radius,normal,candidate))data.neighbors.Append(candidate);}}}
}
//##########################################################################################
// Function Move_Points
//##########################################################################################
template<class T> void POINT_REPULSION<T>::
Move_Points()
{
    VECTOR<T,3> position,normal,accumulated_force;int new_triangle;
    ARRAY<VECTOR<T,3> ,VECTOR<int,1> > new_points(1,points.m);
    ARRAY<int,VECTOR<int,1> > new_triangles(1,points.m);
    for(int i=1;i<=points.m;++i){new_points(i)=points(i).position;new_triangles(i)=points(i).triangle;}
    for(int i=1;i<=points.m;++i){
        POINT_REPULSION_DATA<T>& data=points(i);position=data.position;
        normal=surface.Normal(position,data.triangle); 
        accumulated_force=VECTOR<T,3>(0,0,0);
        T radius=Repulsion_Radius(data.triangle),radii=radii_of_curvature(data.triangle);
        for(int j=1;j<=data.neighbors.m;++j){
            VECTOR<T,3> reference_point=points(data.neighbors(j)).position;reference_point-=VECTOR<T,3>::Dot_Product(reference_point-position,normal)*normal;
            T distance=(reference_point-position).Magnitude();accumulated_force+=(position-reference_point)*(max<T>(0,1-distance/radius)/distance);}
        if(accumulated_force.Max_Abs()==0) continue;
        T force=accumulated_force.Magnitude(),f1=min(force*(T).9*radius,(T).1*radius,(T).5*radii);
        VECTOR<T,3> in_plane_displacement=(attenuation_factor*f1/force)*accumulated_force;
        RAY<VECTOR<T,3> > ray(position+in_plane_displacement+normal,-normal); 
        if(Intersection_With_Triangles(ray,vicinity_triangles(data.triangle),new_triangle)){
            new_points(i)=ray.Point(ray.t_max);new_triangles(i)=new_triangle;}}
    for(int i=1;i<=points.m;++i){points(i).position=new_points(i);points(i).triangle=new_triangles(i);}
    attenuation_factor*=(T).95;Update_Points_In_Triangle();Update_Neighbor_Points();Update_Stats();
}
//##########################################################################################
// Function Update_Stats
//##########################################################################################
template<class T> void POINT_REPULSION<T>::
Update_Stats()
{
    T ratio_sum=0,number_sum=0,number_radius=0;int number_nodes_with_neighbors=0;
    for(int i=1;i<=points.m;++i){
        POINT_REPULSION_DATA<T>& data=points(i);number_sum+=data.neighbors.m; number_radius+=Repulsion_Radius(data.triangle);
        if(data.neighbors.m){
            ++number_nodes_with_neighbors;T local_sum=0;
            for(int j=1;j<=data.neighbors.m;++j)local_sum+=(data.position-points(data.neighbors(j)).position).Magnitude();
            ratio_sum+=(local_sum/data.neighbors.m)/Repulsion_Radius(data.triangle);}}
    average_distance_to_neighbor_over_repulsion_radius=ratio_sum/number_nodes_with_neighbors;
    average_number_neighbors=number_sum/points.m;average_repulsion_radius=number_radius/points.m;
}
//#####################################################################
template class POINT_REPULSION<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class POINT_REPULSION<double>;
#endif

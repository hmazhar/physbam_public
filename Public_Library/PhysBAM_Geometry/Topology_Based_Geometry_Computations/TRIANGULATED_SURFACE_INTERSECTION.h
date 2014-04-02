//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGULATED_SURFACE_INTERSECTION__
#define __TRIANGULATED_SURFACE_INTERSECTION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_TRIANGULATED_SURFACE_INTERSECTION.h>

namespace PhysBAM{
template<class T> class TRIANGULATED_SURFACE;
template<class TV> class RAY;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
template<class T>
bool Check_For_Self_Intersection(TRIANGULATED_SURFACE<T>& ts,const T thickness_over_2,const bool update_bounding_boxes,ARRAY<VECTOR<int,2> >* intersecting_segment_triangle_pairs)
{  
    bool segment_mesh_defined=ts.mesh.segment_mesh!=0;if(!segment_mesh_defined) ts.mesh.Initialize_Segment_Mesh();
    bool intersection=ts.Segment_Triangle_Intersection(*ts.mesh.segment_mesh,ts.particles.X,thickness_over_2,update_bounding_boxes,intersecting_segment_triangle_pairs);
    if(!segment_mesh_defined){delete ts.mesh.segment_mesh;ts.mesh.segment_mesh=0;}
    return intersection;
}
template<class T>
bool Find_First_Segment_Triangle_Intersection(TRIANGULATED_SURFACE<T>& ts,const SEGMENT_MESH& test_segment_mesh,ARRAY_VIEW<const VECTOR<T,3> > X,const T thickness_over_2,const int max_coarsening_attempts,
    const bool update_bounding_boxes)
{
    if(ts.Segment_Triangle_Intersection(test_segment_mesh,X,thickness_over_2,update_bounding_boxes)){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        LOG::cout<<"SELF INTERSECTIONS !"<<std::endl;
#endif
        return true;}
    else{
        for(int loops=1;loops<=max_coarsening_attempts;loops++){
            T distance=(1<<loops)*thickness_over_2;
            if(ts.Segment_Triangle_Intersection(test_segment_mesh,X,distance,false)){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
                LOG::cout<<"collision at a proximity < "<<distance<<std::endl;
#endif
                return true;}
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            else LOG::cout<<"ok at a proximity = "<<distance<<std::endl;
#endif
    }}
    return false;
}
template<class T>
bool Segment_Triangle_Intersection(TRIANGULATED_SURFACE<T>& ts,const SEGMENT_MESH& test_segment_mesh,ARRAY_VIEW<const VECTOR<T,3> > X,const T thickness_over_2,const bool update_bounding_boxes,
    ARRAY<VECTOR<int,2> >* intersecting_segment_triangle_pairs)
{  
    bool intersection=false;
    ARRAY<ARRAY<int> > triangles_near_edges(test_segment_mesh.elements.m);ts.Get_Triangles_Near_Edges(triangles_near_edges,test_segment_mesh,X,thickness_over_2,update_bounding_boxes);
    for(int e=1;e<=test_segment_mesh.elements.m;e++){
        SEGMENT_3D<T> segment(X(test_segment_mesh.elements(e)(1)),X(test_segment_mesh.elements(e)(2)));
        for(int k=1;k<=triangles_near_edges(e).m;k++){int t=triangles_near_edges(e)(k);
            TRIANGLE_3D<T> triangle(X(ts.mesh.elements(t)(1)),X(ts.mesh.elements(t)(2)),X(ts.mesh.elements(t)(3)));
            if(INTERSECTION::Intersects(segment,triangle,thickness_over_2)){intersection=true;
                if(intersecting_segment_triangle_pairs) intersecting_segment_triangle_pairs->Append(VECTOR<int,2>(e,t));else return true;}}}
    return intersection;
}
template<class T>
void Get_Triangles_Near_Edges(TRIANGULATED_SURFACE<T>& ts,ARRAY<ARRAY<int> >& triangles_near_edges,const SEGMENT_MESH& test_segment_mesh,ARRAY_VIEW<const VECTOR<T,3> > X,const T thickness_over_2,
    const bool update_bounding_boxes)
{  
    bool hierarchy_defined=ts.hierarchy!=0;if(!hierarchy_defined) ts.Initialize_Hierarchy(false);
    if(!hierarchy_defined || update_bounding_boxes) ts.hierarchy->Update_Boxes(X);

    for(int k=1;k<=test_segment_mesh.elements.m;k++){
        int node1,node2;test_segment_mesh.elements(k).Get(node1,node2);
        RANGE<VECTOR<T,3> > box(X(node1));box.Enlarge_To_Include_Point(X(node2));
        ts.hierarchy->Intersection_List(box,triangles_near_edges(k),thickness_over_2);
        for(int kk=1;kk<=triangles_near_edges(k).m;kk++){int t=triangles_near_edges(k)(kk); // remove elements that contain a node of the edge being tested
            for(int i=1;i<=3;i++) if(ts.mesh.elements(t)(i) == node1 || ts.mesh.elements(t)(i) == node2){triangles_near_edges(k).Remove_Index_Lazy(kk);kk--;break;}}}
    
    if(!hierarchy_defined){delete ts.hierarchy;ts.hierarchy=0;}
}


}
}
#endif

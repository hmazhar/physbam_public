//#####################################################################
// Copyright 2003-2005, Duc Nguyen, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENTED_CURVE_2D  
//##################################################################### 
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/SEGMENT_2D_SEGMENT_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
using namespace PhysBAM;
//####################################################################
// Register this class as read-write
//####################################################################
namespace PhysBAM{
bool Register_Segmented_Curve_2d(){
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<SEGMENTED_CURVE_2D<float> >();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<SEGMENTED_CURVE_2D<double> >();
#endif
    return true;
}
static bool registered=Register_Segmented_Curve_2d();
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> SEGMENTED_CURVE_2D<T>::
SEGMENTED_CURVE_2D(SEGMENT_MESH& segment_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input)
    :SEGMENTED_CURVE<TV>(segment_mesh_input,particles_input)
{
    PHYSBAM_ASSERT(registered);
}
//#####################################################################
// Function Calculate_Signed_Distance
//#####################################################################
template<class T> T SEGMENTED_CURVE_2D<T>:: 
Calculate_Signed_Distance(const TV& location,T thickness_over_two) const
{
    T distance;
    Closest_Point_On_Curve(location,thickness_over_two,0,&distance);
    if(Inside(location,thickness_over_two)) distance*=-1;
    return distance;
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool SEGMENTED_CURVE_2D<T>::  
Inside(const TV& location,T thickness_over_two) const
{
    // early pruning
    if(bounding_box){if(bounding_box->Outside(location,thickness_over_two)) return false;}

    // Not totally robust yet
    T total_signed_angle=0;
    for(int i=1;i<=mesh.elements.m;i++){
        SEGMENT_2D<T> segment=(segment_list)?(*segment_list)(i):Get_Element(i);
        TV v1=segment.x1-location;
        TV v2=segment.x2-location;
        VECTOR<T,1> cross=TV::Cross_Product(v1,v2);
        T signed_angle=TV::Angle_Between(v1,v2);
        if(cross.x<0) signed_angle*=-1;
        total_signed_angle+=signed_angle;}
    T total_angle=abs(total_signed_angle);

    return (abs(total_angle-2*pi)<abs(total_signed_angle-0));
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool SEGMENTED_CURVE_2D<T>::  
Outside(const TV& location,T thickness_over_two) const
{
    return !Inside(location,thickness_over_two);
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool SEGMENTED_CURVE_2D<T>::  
Boundary(const TV& location,T thickness_over_two) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Inside_Any_Segment
//#####################################################################
template<class T> bool SEGMENTED_CURVE_2D<T>::  
Inside_Any_Segment(const TV& point,int& segment_id,const T thickness_over_two) const
{
    for(int i=1;i<=mesh.elements.m;i++){
        SEGMENT_2D<T> segment=(segment_list)?(*segment_list)(i):Get_Element(i);
        if(segment.Inside(point,thickness_over_two)){segment_id=i;return true;}}
    return false;
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,2> SEGMENTED_CURVE_2D<T>::
Normal(const int aggregate) const
{
    SEGMENT_2D<T> segment=(segment_list)?(*segment_list)(aggregate):Get_Element(aggregate);
    return segment.Normal();
}
//#####################################################################
// Function Find_First_Segment_Segment_Intersection
//#####################################################################
template<class T> bool SEGMENTED_CURVE_2D<T>::
Find_First_Segment_Segment_Intersection(const SEGMENT_MESH& test_segment_mesh,ARRAY_VIEW<const TV> X,const T thickness_over_2,const int max_coarsening_attempts,
    const bool update_bounding_boxes)
{
    if(Segment_Segment_Intersection(test_segment_mesh,X,thickness_over_2,update_bounding_boxes)){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        LOG::cout<<"SELF INTERSECTIONS !"<<std::endl;
#endif
        return true;}
    else{
        for(int loops=1;loops<=max_coarsening_attempts;loops++){
            T distance=pow((T)2,loops)*thickness_over_2;
            if(Segment_Segment_Intersection(test_segment_mesh,X,distance,false)){
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
//#####################################################################
// Function Segment_Segment_Intersection
//#####################################################################
template<class T> bool SEGMENTED_CURVE_2D<T>::
Segment_Segment_Intersection(const SEGMENT_MESH& test_segment_mesh,ARRAY_VIEW<const TV> X,const T thickness_over_2,const bool update_bounding_boxes,
    ARRAY<VECTOR<int,2> >* intersecting_segment_segment_pairs)
{  
    bool intersection=false;
    ARRAY<ARRAY<int> > segments_near_segments(test_segment_mesh.elements.m);Get_Segments_Near_Segments(segments_near_segments,test_segment_mesh,X,thickness_over_2,update_bounding_boxes);
    for(int s1=1;s1<=test_segment_mesh.elements.m;s1++){
        SEGMENT_2D<T> segment1(X(test_segment_mesh.elements(s1)(1)),X(test_segment_mesh.elements(s1)(2)));
        for(int k=1;k<=segments_near_segments(s1).m;k++){int s2=segments_near_segments(s1)(k);
            SEGMENT_2D<T> segment2(X(mesh.elements(s2)(1)),X(mesh.elements(s2)(2)));
            if(INTERSECTION::Intersects(segment1,segment2,thickness_over_2)){intersection=true;
                if(intersecting_segment_segment_pairs) intersecting_segment_segment_pairs->Append(VECTOR<int,2>(s1,s2));else return true;}}}
    return intersection;
}
//#####################################################################
// Function Get_Triangles_Near_Edges
//#####################################################################
template<class T> void SEGMENTED_CURVE_2D<T>::
Get_Segments_Near_Segments(ARRAY<ARRAY<int> >& segments_near_segments,const SEGMENT_MESH& test_segment_mesh,ARRAY_VIEW<const TV> X,const T thickness_over_2,
    const bool update_bounding_boxes)
{  
    bool hierarchy_defined=hierarchy!=0;if(!hierarchy_defined) Initialize_Hierarchy(false);
    if(!hierarchy_defined || update_bounding_boxes) hierarchy->Update_Boxes(X);

    for(int k=1;k<=test_segment_mesh.elements.m;k++){
        int node1,node2;test_segment_mesh.elements(k).Get(node1,node2);
        RANGE<TV> box(X(node1));box.Enlarge_To_Include_Point(X(node2));
        hierarchy->Intersection_List(box,segments_near_segments(k),thickness_over_2);
        for(int kk=1;kk<=segments_near_segments(k).m;kk++){int s=segments_near_segments(k)(kk); // remove elements that contain a node of the segment being tested
            for(int i=1;i<=2;i++) if(mesh.elements(s)(i) == node1 || mesh.elements(s)(i) == node2){segments_near_segments(k).Remove_Index_Lazy(kk);kk--;break;}}}
    
    if(!hierarchy_defined){delete hierarchy;hierarchy=0;}
}
//####################################################################
template class SEGMENTED_CURVE_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SEGMENTED_CURVE_2D<double>;
#endif

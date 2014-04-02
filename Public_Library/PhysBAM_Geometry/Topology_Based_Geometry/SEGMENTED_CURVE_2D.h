//#####################################################################
// Copyright 2003-2007, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Duc Nguyen, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENTED_CURVE_2D
//##################################################################### 
#ifndef __SEGMENTED_CURVE_2D__
#define __SEGMENTED_CURVE_2D__

#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class T_input>
class SEGMENTED_CURVE_2D:public SEGMENTED_CURVE<VECTOR<T_input,2> >
{
    typedef T_input T;typedef VECTOR<T,2> TV;
public:
    typedef SEGMENTED_CURVE<TV> BASE;
    using BASE::mesh;using BASE::particles;using BASE::segment_list;using BASE::bounding_box;using BASE::Update_Bounding_Box;using BASE::hierarchy;
    using BASE::Get_Element;using BASE::Closest_Point_On_Curve;using BASE::Total_Length;using BASE::Initialize_Hierarchy;

    SEGMENTED_CURVE_2D(SEGMENT_MESH& segment_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input);

    bool Inside_Any_Simplex(const TV& point,int& segment_id,const T thickness_over_two=0) const
    {return Inside_Any_Segment(point,segment_id,thickness_over_two);}

    TV Normal(const TV& location,const int aggregate) const
    {return Normal(aggregate);}

//#####################################################################
    T Calculate_Signed_Distance(const TV& location,T thickness_over_two=0) const;
    bool Boundary(const TV& location,T thickness_over_two=0) const;
    bool Inside(const TV& location,T thickness_over_two=0) const;
    bool Outside(const TV& location,T thickness_over_two=0) const;
    bool Inside_Any_Segment(const TV& point,int& segment_id,const T thickness_over_two=0) const;
    bool Segment_Segment_Intersection(const SEGMENT_MESH& test_segment_mesh,ARRAY_VIEW<const TV> X,const T thickness_over_2=0,const bool update_bounding_boxes=true,
        ARRAY<VECTOR<int,2> >* intersecting_segment_segment_pairs=0);
    bool Find_First_Segment_Segment_Intersection(const SEGMENT_MESH& test_segment_mesh,ARRAY_VIEW<const TV> X,const T thickness_over_2,const int max_coarsening_attempts=5,
        const bool update_bounding_boxes=true);
    void Get_Segments_Near_Segments(ARRAY<ARRAY<int> >& segments_near_segments,const SEGMENT_MESH& test_segment_mesh,ARRAY_VIEW<const TV> X,const T thickness_over_2,
        const bool update_bounding_boxes);
    TV Normal(const int aggregate) const;
//#####################################################################
};
}
#endif

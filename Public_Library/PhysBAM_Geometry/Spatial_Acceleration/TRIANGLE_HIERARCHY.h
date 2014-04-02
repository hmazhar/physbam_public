//#####################################################################
// Copyright 2002-2005, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Andrew Selle, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_HIERARCHY
//##################################################################### 
#ifndef __TRIANGLE_HIERARCHY__
#define __TRIANGLE_HIERARCHY__

#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
namespace PhysBAM{

template<class T_input>
class TRIANGLE_HIERARCHY:public BOX_HIERARCHY<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef BOX_HIERARCHY<VECTOR<T,3> > BASE;
    using BASE::leaves;using BASE::root;using BASE::parents;using BASE::children;using BASE::box_hierarchy;using BASE::box_radius;using BASE::Leaf;using BASE::Update_Nonleaf_Boxes;
    using BASE::Intersection_List;

    TRIANGLE_MESH& triangle_mesh;
    GEOMETRY_PARTICLES<VECTOR<T,3> >& particles;
    ARRAY<TRIANGLE_3D<T> >* triangle_list;
    ARRAY<ARRAY<int> > triangles_in_group;
    int triangles_per_group;
    
    TRIANGLE_HIERARCHY(TRIANGLE_MESH& triangle_mesh_input,GEOMETRY_PARTICLES<VECTOR<T,3> >& particles_input,const bool update_boxes=true,const int triangles_per_group=0);
    TRIANGLE_HIERARCHY(TRIANGLE_MESH& triangle_mesh_input,GEOMETRY_PARTICLES<VECTOR<T,3> >& particles_input,ARRAY<TRIANGLE_3D<T> >& triangle_list_input,const bool update_boxes=true,const int triangles_per_group=0);
    virtual ~TRIANGLE_HIERARCHY();

    void Update_Boxes(const T extra_thickness=0)
    {Update_Leaf_Boxes(extra_thickness);Update_Nonleaf_Boxes();}

    template<class T_ARRAY_TV>
    void Update_Boxes(const T_ARRAY_TV& X,const T extra_thickness=0) // use X instead of the current particle positions
    {Update_Leaf_Boxes(X,extra_thickness);Update_Nonleaf_Boxes();}

    template<class T_ARRAY_TV>
    void Update_Boxes(const T_ARRAY_TV& start_X,const T_ARRAY_TV& end_X,const T extra_thickness=0) // bound triangles moving from start_X to end_X
    {Update_Leaf_Boxes(start_X,end_X,extra_thickness);Update_Nonleaf_Boxes();}

    void Update_Boxes(const FRAME<TV>& start_frame,const FRAME<TV>& end_frame,const T extra_thickness=0) // for moving triangles
    {Update_Leaf_Boxes(start_frame,end_frame,extra_thickness);Update_Nonleaf_Boxes();}

    void Update_Leaf_Boxes(const T extra_thickness=0)
    {Calculate_Bounding_Boxes(box_hierarchy);if(extra_thickness) Thicken_Leaf_Boxes(extra_thickness);}

    template<class T_ARRAY_TV>
    void Update_Leaf_Boxes(const T_ARRAY_TV& X,const T extra_thickness=0) // use X instead of the current particle positions
    {Calculate_Bounding_Boxes(box_hierarchy,X);if(extra_thickness) Thicken_Leaf_Boxes(extra_thickness);}
    
    template<class T_ARRAY_TV>
    void Update_Leaf_Boxes(const T_ARRAY_TV& start_X,const T_ARRAY_TV& end_X,const T extra_thickness=0) // bound triangles moving from start_X to end_X
    {Calculate_Bounding_Boxes(box_hierarchy,start_X,end_X);if(extra_thickness) Thicken_Leaf_Boxes(extra_thickness);}

    void Update_Leaf_Boxes(const FRAME<TV>& start_frame,const FRAME<TV>& end_frame,const T extra_thickness=0) // for moving triangles
    {Calculate_Bounding_Boxes(box_hierarchy,start_frame,end_frame);if(extra_thickness) Thicken_Leaf_Boxes(extra_thickness);}

    void Intersection_List(const VECTOR<T,3>& point,ARRAY<int>& intersection_list,const T thickness_over_two=0) const PHYSBAM_OVERRIDE
    {if(triangles_per_group){
        ARRAY<int> group_list;group_list.Preallocate(10);
        Intersection_List(root,point,group_list,thickness_over_two);
        for(int i=1;i<=group_list.m;i++) intersection_list.Append_Elements(triangles_in_group(group_list(i)));}
    else Intersection_List(root,point,intersection_list,thickness_over_two);}
    
    void Intersection_List(const RANGE<TV>& test_box,ARRAY<int>& intersection_list,const T thickness_over_two=0) const PHYSBAM_OVERRIDE
    {if(triangles_per_group){
        ARRAY<int> group_list;group_list.Preallocate(10);
        Intersection_List(root,test_box,group_list,thickness_over_two);
        for(int i=1;i<=group_list.m;i++) intersection_list.Append_Elements(triangles_in_group(group_list(i)));}
    else Intersection_List(root,test_box,intersection_list,thickness_over_two);}
        
    void Intersection_List(const ORIENTED_BOX<TV>& test_box,ARRAY<int>& intersection_list) const PHYSBAM_OVERRIDE
    {if(triangles_per_group){
        ARRAY<int> group_list;group_list.Preallocate(10);
        Intersection_List(root,test_box,group_list);
        for(int i=1;i<=group_list.m;i++) intersection_list.Append_Elements(triangles_in_group(group_list(i)));}
    else Intersection_List(root,test_box,intersection_list);}
    
    void Intersection_List(const PLANE<T>& test_plane,ARRAY<int>& intersection_list,const T thickness_over_two=0) const PHYSBAM_OVERRIDE
    {if(triangles_per_group){
        ARRAY<int> group_list;group_list.Preallocate(10);
        Intersection_List(root,test_plane,group_list,thickness_over_two);
        for(int i=1;i<=group_list.m;i++) intersection_list.Append_Elements(triangles_in_group(group_list(i)));}
    else Intersection_List(root,test_plane,intersection_list,thickness_over_two);}

    void Intersection_List(const IMPLICIT_OBJECT<VECTOR<T,3> >& implicit_surface,const MATRIX<T,3>& rotation,const VECTOR<T,3>& translation,ARRAY<int>& intersection_list,
        const T contour_value=0) const
    {if(triangles_per_group){
        ARRAY<int> group_list;group_list.Preallocate(10);
        Intersection_List(root,implicit_surface,rotation,translation,group_list,contour_value);
        for(int i=1;i<=group_list.m;i++) intersection_list.Append_Elements(triangles_in_group(group_list(i)));}
    else Intersection_List(root,implicit_surface,rotation,translation,intersection_list,contour_value);}

    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,ARRAY_VIEW<const TV> X)
    {Calculate_Bounding_Boxes_Helper(bounding_boxes,X);}

    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,INDIRECT_ARRAY<ARRAY_VIEW<const TV> > X)
    {Calculate_Bounding_Boxes_Helper(bounding_boxes,X);}

    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,ARRAY_VIEW<const TV> start_X,ARRAY_VIEW<const TV> end_X)
    {Calculate_Bounding_Boxes_Helper(bounding_boxes,start_X,end_X);}

    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,INDIRECT_ARRAY<ARRAY_VIEW<const TV> > start_X,INDIRECT_ARRAY<ARRAY_VIEW<const TV> > end_X)
    {Calculate_Bounding_Boxes_Helper(bounding_boxes,start_X,end_X);}

//#####################################################################
    void Initialize_Hierarchy_Using_KD_Tree() PHYSBAM_OVERRIDE;
    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes);
    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,const FRAME<TV>& start_frame,const FRAME<TV>& end_frame);
    void Calculate_Bounding_Box_Radii(const ARRAY<RANGE<TV> >& bounding_boxes,ARRAY<T>& radius) PHYSBAM_OVERRIDE;
    bool Intersection(RAY<VECTOR<T,3> >& ray,const T thickness_over_two=0,const bool use_ray_bounding_box=false) const;
    bool Intersection(RAY<VECTOR<T,3> >& ray,VECTOR<T,3>& weights,const T thickness_over_two=0,const bool use_ray_bounding_box=false) const;
    // for internal use - but octrees use them as well so they're not private
    bool Intersection(const int box,RAY<VECTOR<T,3> >& ray,const T thickness,const T thickness_over_two,const bool use_ray_bounding_box=false) const;
private:
    template<class T_ARRAY_TV> void Calculate_Bounding_Boxes_Helper(ARRAY<RANGE<TV> >& bounding_boxes,T_ARRAY_TV X);
    template<class T_ARRAY_TV> void Calculate_Bounding_Boxes_Helper(ARRAY<RANGE<TV> >& bounding_boxes,T_ARRAY_TV start_X,T_ARRAY_TV end_X);
//#####################################################################
};   
}
#endif



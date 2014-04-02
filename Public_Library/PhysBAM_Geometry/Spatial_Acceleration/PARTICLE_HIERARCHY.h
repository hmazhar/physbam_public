//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_HIERARCHY
//#####################################################################
#ifndef __PARTICLE_HIERARCHY__
#define __PARTICLE_HIERARCHY__

#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SPATIAL_ACCELERATION_FORWARD.h>
namespace PhysBAM{

template<class TV,class T_ARRAY>
class PARTICLE_HIERARCHY:public BOX_HIERARCHY<TV>
{
private:
    typedef typename TV::SCALAR T;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::HYPERPLANE T_HYPERPLANE;
public:
    typedef BOX_HIERARCHY<TV> BASE;
    using BASE::leaves;using BASE::root;using BASE::parents;using BASE::children;using BASE::box_hierarchy;using BASE::box_radius;using BASE::Leaf;using BASE::Update_Nonleaf_Boxes;
    using BASE::Intersection_List;

    const T_ARRAY& X;
    ARRAY<ARRAY<int> > particles_in_group;
    int particles_per_group;

    PARTICLE_HIERARCHY(const T_ARRAY& X_input,const bool update_boxes=true,const int particles_per_group_input=10);
    virtual ~PARTICLE_HIERARCHY();

    void Update_Boxes(const T extra_thickness=0)
    {Update_Leaf_Boxes(extra_thickness);Update_Nonleaf_Boxes();}

    void Update_Boxes(const ARRAY<TV>& X,const T extra_thickness=0)
    {Update_Leaf_Boxes(X,extra_thickness);Update_Nonleaf_Boxes();}

    void Update_Leaf_Boxes(const T extra_thickness=0)
    {Calculate_Bounding_Boxes(box_hierarchy);if(extra_thickness) Thicken_Leaf_Boxes(extra_thickness);}

    template<class T_ARRAY_TV> void Update_Leaf_Boxes(const T_ARRAY_TV& X,const T extra_thickness=0)
    {Calculate_Bounding_Boxes(box_hierarchy,X);if(extra_thickness) Thicken_Leaf_Boxes(extra_thickness);}

    void Intersection_List(const TV& point,ARRAY<int>& intersection_list,const T thickness_over_two=0) const PHYSBAM_OVERRIDE
    {if(particles_per_group){
        ARRAY<int> group_list;group_list.Preallocate(10);
        BASE::Intersection_List(root,point,group_list,thickness_over_two);
        for(int i=1;i<=group_list.m;i++) intersection_list.Append_Elements(particles_in_group(group_list(i)));}
    else BASE::Intersection_List(root,point,intersection_list,thickness_over_two);}

    void Intersection_List(const RANGE<TV>& test_box,ARRAY<int>& intersection_list,const T thickness_over_two=0) const PHYSBAM_OVERRIDE
    {if(particles_per_group){
        ARRAY<int> group_list;group_list.Preallocate(10);
        BASE::Intersection_List(root,test_box,group_list,thickness_over_two);
        for(int i=1;i<=group_list.m;i++) intersection_list.Append_Elements(particles_in_group(group_list(i)));}
    else BASE::Intersection_List(root,test_box,intersection_list,thickness_over_two);}

    void Intersection_List(const ORIENTED_BOX<TV>& test_box,ARRAY<int>& intersection_list) const PHYSBAM_OVERRIDE
    {if(particles_per_group){
        ARRAY<int> group_list;group_list.Preallocate(10);
        BASE::Intersection_List(root,test_box,group_list);
        for(int i=1;i<=group_list.m;i++) intersection_list.Append_Elements(particles_in_group(group_list(i)));}
    else BASE::Intersection_List(root,test_box,intersection_list);}

    void Intersection_List(const T_HYPERPLANE& test_plane,ARRAY<int>& intersection_list,const T thickness_over_two=0) const PHYSBAM_OVERRIDE
    {if(particles_per_group){
        ARRAY<int> group_list;group_list.Preallocate(10);
        BASE::Intersection_List(root,test_plane,group_list,thickness_over_two);
        for(int i=1;i<=group_list.m;i++) intersection_list.Append_Elements(particles_in_group(group_list(i)));}
    else BASE::Intersection_List(root,test_plane,intersection_list,thickness_over_two);}

    void Intersection_List(const IMPLICIT_OBJECT<TV>& implicit_object,const MATRIX<T,TV::dimension>& rotation,const TV& translation,ARRAY<int>& intersection_list,const T contour_value=0) const PHYSBAM_OVERRIDE
    {if(particles_per_group){
        ARRAY<int> group_list;group_list.Preallocate(10);
        BASE::Intersection_List(root,implicit_object,rotation,translation,group_list,contour_value);
        for(int i=1;i<=group_list.m;i++) intersection_list.Append_Elements(particles_in_group(group_list(i)));}
    else BASE::Intersection_List(root,implicit_object,rotation,translation,intersection_list,contour_value);}

//#####################################################################
    void Initialize_Hierarchy_Using_KD_Tree() PHYSBAM_OVERRIDE;
    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes) const;
    template<class T_ARRAY_TV> void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,const T_ARRAY_TV& X) const;
    void Calculate_Bounding_Box_Radii(const ARRAY<RANGE<TV> >& bounding_boxes,ARRAY<T>& radius) PHYSBAM_OVERRIDE;
//#####################################################################
};
//#####################################################################
// Function Calculate_Bounding_Boxes
//#####################################################################
template<class TV,class T_ARRAY> template<class T_ARRAY_TV> void PARTICLE_HIERARCHY<TV,T_ARRAY>::
Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,const T_ARRAY_TV& X) const
{
    STATIC_ASSERT((IS_SAME<TV,typename T_ARRAY_TV::ELEMENT>::value));
    if(particles_per_group) for(int k=1;k<=leaves;k++){
        if(particles_in_group(k).m) bounding_boxes(k)=RANGE<TV>::Bounding_Box(X.Subset(particles_in_group(k)));}
    else for(int k=1;k<=leaves;k++) bounding_boxes(k).Reset_Bounds(X(k));
}
//#####################################################################
}
#endif

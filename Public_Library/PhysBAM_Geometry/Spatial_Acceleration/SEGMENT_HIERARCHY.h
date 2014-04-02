//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Andrew Selle, Joseph Teran, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_HIERARCHY
//##################################################################### 
#ifndef __SEGMENT_HIERARCHY__
#define __SEGMENT_HIERARCHY__

#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
namespace PhysBAM{
class SEGMENT_MESH;

template<class TV>
class SEGMENT_HIERARCHY:public BOX_HIERARCHY<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::SEGMENT T_SEGMENT;
public:
    typedef BOX_HIERARCHY<TV> BASE;
    using BASE::leaves;using BASE::root;using BASE::parents;using BASE::children;using BASE::box_hierarchy;using BASE::box_radius;using BASE::Update_Nonleaf_Boxes;using BASE::Intersection_List;

    SEGMENT_MESH& segment_mesh;
    GEOMETRY_PARTICLES<TV>& particles;
    ARRAY<T_SEGMENT>* segment_list;
    
    SEGMENT_HIERARCHY(SEGMENT_MESH& segment_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input,const bool update_boxes=true);
    SEGMENT_HIERARCHY(SEGMENT_MESH& segment_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input,ARRAY<T_SEGMENT>& segment_list_input,const bool update_boxes=true);
    virtual ~SEGMENT_HIERARCHY();

    void Update_Boxes(const T extra_thickness=0)
    {Update_Leaf_Boxes(extra_thickness);Update_Nonleaf_Boxes();}

    template<class T_ARRAY_TV>
    void Update_Boxes(const T_ARRAY_TV& X,const T extra_thickness=0) // use X instead of the current particle positions
    {Update_Leaf_Boxes(X,extra_thickness);Update_Nonleaf_Boxes();}

    template<class T_ARRAY_TV>
    void Update_Boxes(const T_ARRAY_TV& start_X,const T_ARRAY_TV& end_X,const T extra_thickness=0) // bound segments moving from start_X to end_X
    {Update_Leaf_Boxes(start_X,end_X,extra_thickness);Update_Nonleaf_Boxes();}

    void Update_Boxes(const FRAME<TV>& start_frame,const FRAME<TV>& end_frame,const T extra_thickness=0) // for a moving rigid body
    {Update_Leaf_Boxes(start_frame,end_frame,extra_thickness);Update_Nonleaf_Boxes();}

    void Update_Leaf_Boxes(const T extra_thickness=0)
    {Calculate_Bounding_Boxes(box_hierarchy);if(extra_thickness) Thicken_Leaf_Boxes(extra_thickness);}

    template<class T_ARRAY_TV>
    void Update_Leaf_Boxes(const T_ARRAY_TV& X,const T extra_thickness=0) // use X instead of the current particle positions
    {Calculate_Bounding_Boxes(box_hierarchy,X);if(extra_thickness) Thicken_Leaf_Boxes(extra_thickness);}
    
    template<class T_ARRAY_TV>
    void Update_Leaf_Boxes(const T_ARRAY_TV& start_X,const T_ARRAY_TV& end_X,const T extra_thickness=0) // bound segments moving from start_X to end_X
    {Calculate_Bounding_Boxes(box_hierarchy,start_X,end_X);if(extra_thickness) Thicken_Leaf_Boxes(extra_thickness);}

    void Update_Leaf_Boxes(const FRAME<TV>& start_frame,const FRAME<TV>& end_frame,const T extra_thickness=0) // for a moving rigid body
    {Calculate_Bounding_Boxes(box_hierarchy,start_frame,end_frame);if(extra_thickness) Thicken_Leaf_Boxes(extra_thickness);}

    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes)
    {Calculate_Bounding_Boxes(bounding_boxes,particles.X);}

    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,ARRAY_VIEW<const TV> X)
    {Calculate_Bounding_Boxes_Helper(bounding_boxes,X);}

    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,ARRAY_VIEW<const TV> start_X,ARRAY_VIEW<const TV> end_X)
    {Calculate_Bounding_Boxes_Helper(bounding_boxes,start_X,end_X);}

    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,INDIRECT_ARRAY<ARRAY_VIEW<const TV> > X)
    {Calculate_Bounding_Boxes_Helper(bounding_boxes,X);}

//#####################################################################
    void Initialize_Hierarchy_Using_KD_Tree() PHYSBAM_OVERRIDE;
    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,const FRAME<TV>& start_frame,const FRAME<TV>& end_frame);
    void Calculate_Bounding_Box_Radii(const ARRAY<RANGE<TV> >& bounding_boxes,ARRAY<T>& radius) PHYSBAM_OVERRIDE;
private:
    template<class T_ARRAY_TV> void Calculate_Bounding_Boxes_Helper(ARRAY<RANGE<TV> >& bounding_boxes,T_ARRAY_TV X);
    template<class T_ARRAY_TV> void Calculate_Bounding_Boxes_Helper(ARRAY<RANGE<TV> >& bounding_boxes,T_ARRAY_TV start_X,T_ARRAY_TV end_X);
//#####################################################################
};   
}
#endif

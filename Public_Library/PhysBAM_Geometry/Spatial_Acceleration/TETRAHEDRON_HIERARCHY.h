//#####################################################################
// Copyright 2004, Zhaosheng Bao, Ron Fedkiw, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TETRAHEDRON_HIERARCHY
//##################################################################### 
#ifndef __TETRAHEDRON_HIERARCHY__
#define __TETRAHEDRON_HIERARCHY__

#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology/TETRAHEDRON_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class T_input>
class TETRAHEDRON_HIERARCHY:public BOX_HIERARCHY<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef BOX_HIERARCHY<TV> BASE;
    using BASE::leaves;using BASE::root;using BASE::parents;using BASE::children;using BASE::box_hierarchy;using BASE::box_radius;using BASE::Update_Nonleaf_Boxes;

    TETRAHEDRON_MESH& tetrahedron_mesh;
    GEOMETRY_PARTICLES<TV>& particles;
    ARRAY<TETRAHEDRON<T> >* tetrahedron_list;
    
    TETRAHEDRON_HIERARCHY(TETRAHEDRON_MESH& tetrahedron_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input,const bool update_boxes=true);
    TETRAHEDRON_HIERARCHY(TETRAHEDRON_MESH& tetrahedron_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input,ARRAY<TETRAHEDRON<T> >& tetrahedron_list_input,
        const bool update_boxes=true);
    virtual ~TETRAHEDRON_HIERARCHY();

    void Update_Boxes(const T extra_thickness=0)
    {Update_Leaf_Boxes(extra_thickness);Update_Nonleaf_Boxes();}

    void Update_Boxes(const ARRAY<TV>& X,const T extra_thickness=0) // use X instead of the current particle positions
    {Update_Leaf_Boxes(X,extra_thickness);Update_Nonleaf_Boxes();}

    void Update_Boxes(const ARRAY<TV>& start_X,const ARRAY<TV>& end_X,const T extra_thickness=0) // bound tetrahedrons moving from start_X to end_X
    {Update_Leaf_Boxes(start_X,end_X,extra_thickness);Update_Nonleaf_Boxes();}

    void Update_Boxes(const FRAME<TV>& start_frame,const FRAME<TV>& end_frame,const T extra_thickness=0) // for a moving rigid body
    {Update_Leaf_Boxes(start_frame,end_frame,extra_thickness);Update_Nonleaf_Boxes();}

    void Update_Leaf_Boxes(const T extra_thickness=0)
    {Calculate_Bounding_Boxes(box_hierarchy);if(extra_thickness) Thicken_Leaf_Boxes(extra_thickness);}

    void Update_Leaf_Boxes(const ARRAY<TV>& X,const T extra_thickness=0) // use X instead of the current particle positions
    {Calculate_Bounding_Boxes(box_hierarchy,X);if(extra_thickness) Thicken_Leaf_Boxes(extra_thickness);}
    
    void Update_Leaf_Boxes(const ARRAY<TV>& start_X,const ARRAY<TV>& end_X,const T extra_thickness=0) // bound tetrahedrons moving from start_X to end_X
    {Calculate_Bounding_Boxes(box_hierarchy,start_X,end_X);if(extra_thickness) Thicken_Leaf_Boxes(extra_thickness);}

    void Update_Leaf_Boxes(const FRAME<TV>& start_frame,const FRAME<TV>& end_frame,const T extra_thickness=0) // for a moving rigid body
    {Calculate_Bounding_Boxes(box_hierarchy,start_frame,end_frame);if(extra_thickness) Thicken_Leaf_Boxes(extra_thickness);}

    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes)
    {Calculate_Bounding_Boxes(bounding_boxes,particles.X);}
    
//#####################################################################
    void Initialize_Hierarchy_Using_KD_Tree() PHYSBAM_OVERRIDE;
    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,ARRAY_VIEW<const TV> X);
    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,ARRAY_VIEW<const TV> start_X,ARRAY_VIEW<const TV> end_X);
    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,const FRAME<TV>& start_frame,const FRAME<TV>& end_frame);
    void Calculate_Bounding_Box_Radii(const ARRAY<RANGE<TV> >& bounding_boxes,ARRAY<T>& radius) PHYSBAM_OVERRIDE;
//#####################################################################
};   
}
#endif


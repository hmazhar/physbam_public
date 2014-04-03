//#####################################################################
// Copyright 2002-2004, Zhaosheng Bao.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_HIERARCHY_2D
//##################################################################### 
#ifndef __TRIANGLE_HIERARCHY_2D__
#define __TRIANGLE_HIERARCHY_2D__

#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
namespace PhysBAM{

template<class T_input>
class TRIANGLE_HIERARCHY_2D:public BOX_HIERARCHY<VECTOR<T_input,2> >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;
public:
    typedef BOX_HIERARCHY<TV> BASE;
    using BASE::leaves;using BASE::root;using BASE::parents;using BASE::children;using BASE::box_hierarchy;using BASE::box_radius;using BASE::Update_Nonleaf_Boxes;

    TRIANGLE_MESH& triangle_mesh;
    GEOMETRY_PARTICLES<TV>& particles;
    ARRAY<TRIANGLE_2D<T> >* triangle_list;
    
    TRIANGLE_HIERARCHY_2D(TRIANGLE_MESH& triangle_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input,const bool update_boxes=true)
        :triangle_mesh(triangle_mesh_input),particles(particles_input),triangle_list(0)
    {  
        Initialize_Hierarchy_Using_KD_Tree();
        if(update_boxes) Update_Boxes();
    }
    
    TRIANGLE_HIERARCHY_2D(TRIANGLE_MESH& triangle_mesh_input,GEOMETRY_PARTICLES<TV>& particles_input,ARRAY<TRIANGLE_2D<T> >& triangle_list_input,
        const bool update_boxes=true)
        :triangle_mesh(triangle_mesh_input),particles(particles_input),triangle_list(&triangle_list_input)
    {  
        Initialize_Hierarchy_Using_KD_Tree();
        if(update_boxes) Update_Boxes();
    }

    virtual ~TRIANGLE_HIERARCHY_2D()
    {}

    void Update_Boxes()
    {Update_Leaf_Boxes();Update_Nonleaf_Boxes();}

    template<class T_ARRAY_TV>
    void Update_Boxes(const T_ARRAY_TV& X) // use X instead of the current particle positions
    {Update_Leaf_Boxes(X);Update_Nonleaf_Boxes();}

    template<class T_ARRAY_TV>
    void Update_Boxes(const T_ARRAY_TV& start_X,const T_ARRAY_TV& end_X) // bound triangles moving from start_X to end_X
    {Update_Leaf_Boxes(start_X,end_X);Update_Nonleaf_Boxes();}

    void Update_Leaf_Boxes()
    {Calculate_Bounding_Boxes(box_hierarchy);}

    template<class T_ARRAY_TV>
    void Update_Leaf_Boxes(const T_ARRAY_TV& X) // use X instead of the current particle positions
    {Calculate_Bounding_Boxes(box_hierarchy,X);}
    
    template<class T_ARRAY_TV>
    void Update_Leaf_Boxes(const T_ARRAY_TV& start_X,const T_ARRAY_TV& end_X) // bound triangles moving from start_X to end_X
    {Calculate_Bounding_Boxes(box_hierarchy,start_X,end_X);}

    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes)
    {Calculate_Bounding_Boxes(bounding_boxes,particles.X);}
       
//#####################################################################
    void Initialize_Hierarchy_Using_KD_Tree() PHYSBAM_OVERRIDE;
    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,ARRAY_VIEW<const TV> X);
    void Calculate_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,ARRAY_VIEW<const TV> start_X,ARRAY_VIEW<const TV> end_X);
    void Calculate_Bounding_Box_Radii(const ARRAY<RANGE<TV> >& bounding_boxes,ARRAY<T>& radius) PHYSBAM_OVERRIDE;
//#####################################################################
};   
}
#endif


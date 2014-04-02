//#####################################################################
// Copyright 2003-2008, Eran Guendelman, Sergey Koltakov, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_SELECTION
//##################################################################### 
#ifndef __OPENGL_SELECTION__
#define __OPENGL_SELECTION__

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Topology/SIMPLEX_MESH.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
namespace PhysBAM{

class OPENGL_SELECTION:public NONCOPYABLE
{
public:
    enum TYPE { GRID_CELL_1D,GRID_CELL_2D, GRID_NODE_2D, GRID_CELL_3D, GRID_CELL_LIST_3D, GRID_NODE_3D, GRID_NODE_LIST_3D, RLE_CELL_2D, RLE_CELL_3D, BINTREE_CELL, BINTREE_NODE, QUADTREE_CELL, QUADTREE_NODE, OCTREE_NODE, OCTREE_CELL, 
                OCTREE_FACE, POINTS_2D, POINTS_3D, COMPONENT_PARTICLES_2D, COMPONENT_PARTICLES_3D, 
                POINT_SIMPLICES_1D,
                SEGMENTED_CURVE_VERTEX_2D, SEGMENTED_CURVE_SEGMENT_2D, SEGMENTED_CURVE_VERTEX_3D, SEGMENTED_CURVE_SEGMENT_3D, SEGMENTED_CURVE_3D,
                TRIANGULATED_SURFACE_VERTEX, TRIANGULATED_SURFACE_SEGMENT, TRIANGULATED_SURFACE_TRIANGLE, 
                TRIANGULATED_AREA_VERTEX, TRIANGULATED_AREA_SEGMENT, TRIANGULATED_AREA_TRIANGLE,
                TETRAHEDRALIZED_VOLUME_VERTEX, TETRAHEDRALIZED_VOLUME_TETRAHEDRON, COMPONENT_RIGID_BODIES_2D, COMPONENT_RIGID_BODIES_3D, 
                ARTICULATED_RIGID_BODIES_JOINT_2D, ARTICULATED_RIGID_BODIES_JOINT_3D, ARTICULATED_RIGID_BODIES_MUSCLE_2D, MUSCLE_3D, MUSCLE_SURFACE_3D,
                COMPONENT_DEFORMABLE_COLLECTION_1D,
                COMPONENT_DEFORMABLE_OBJECT_2D, COMPONENT_DEFORMABLE_OBJECT_2D_INTERACTIVE, 
                COMPONENT_DEFORMABLE_COLLECTION_3D, COMPONENT_DEFORMABLE_COLLECTION_3D_INTERACTIVE,
                COMPONENT_HEIGHTFIELD_1D, COMPONENT_HEIGHTFIELD_2D, COMPONENT_CURVE_VERTEX_2D, COMPONENT_CURVE_SEGMENT_2D, COMPONENT_MUSCLES_2D,DEBUG_PARTICLES_2D };

    TYPE type;
    OPENGL_OBJECT *object;
    float min_depth, max_depth;
    bool hide;

    OPENGL_SELECTION(TYPE type,OPENGL_OBJECT* object);
    virtual ~OPENGL_SELECTION();
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const=0;
    virtual TYPE Actual_Type() const;

//#####################################################################
// Useful functions to draw highlighted primitives
//#####################################################################
    template<class TV> static void Draw_Highlighted_Vertex(const TV& position,int id=0,const OPENGL_COLOR& color=OPENGL_PREFERENCES::selection_highlight_color);
    template<class TV> static void Draw_Highlighted_Segment(const TV& x1,const TV& x2,int id=0,const OPENGL_COLOR& color=OPENGL_PREFERENCES::selection_highlight_color);
    template<class TV> static void Draw_Highlighted_Curve(const ARRAY<VECTOR<TV,2> >& X,int id=0,const OPENGL_COLOR& color=OPENGL_PREFERENCES::selection_highlight_color);
    template<class TV> static void Draw_Highlighted_Triangle_Boundary(const TV& x1,const TV& x2,const TV& x3,int id=0,const OPENGL_COLOR& color=OPENGL_PREFERENCES::selection_highlight_color);
    template<class T> static void Draw_Highlighted_Tetrahedron_Boundary(const VECTOR<T,3>& x1,const VECTOR<T,3>& x2,const VECTOR<T,3>& x3,const VECTOR<T,3>& x4,int id=0,
        const OPENGL_COLOR& color=OPENGL_PREFERENCES::selection_highlight_color);
    template<class T> static void Draw_Highlighted_Quad(const VECTOR<T,2>& x00,const VECTOR<T,2>& x11,const OPENGL_COLOR& color=OPENGL_PREFERENCES::selection_highlight_color);
    template<class T> static void Draw_Highlighted_Quad(const VECTOR<T,3>& node1,const VECTOR<T,3>& node2,const VECTOR<T,3>& node3,const VECTOR<T,3>& node4,
        const OPENGL_COLOR& color=OPENGL_PREFERENCES::selection_highlight_color);
    template<class T> static void Draw_Highlighted_Box(const VECTOR<T,3>& x000,const VECTOR<T,3>& x111,const OPENGL_COLOR& color=OPENGL_PREFERENCES::selection_highlight_color);
    template<class TV,int d> static void Draw_Vertices_For_Selection(const SIMPLEX_MESH<d>& mesh,const GEOMETRY_PARTICLES<TV>& particles);
};

}

#endif

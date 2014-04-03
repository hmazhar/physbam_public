//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_TETRAHEDRALIZED_VOLUME
//#####################################################################
#ifndef __OPENGL_TETRAHEDRALIZED_VOLUME__
#define __OPENGL_TETRAHEDRALIZED_VOLUME__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Topology/TETRAHEDRON_MESH.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T>
class OPENGL_TETRAHEDRALIZED_VOLUME:public OPENGL_OBJECT
{
    typedef VECTOR<T,3> TV;
public:
    OPENGL_MATERIAL material;
    TETRAHEDRON_MESH* mesh;
    GEOMETRY_PARTICLES<VECTOR<T,3> >* particles;
    int current_tetrahedron;
    int current_node;
    int current_boundary_triangle;
    int minimum_valence;
    ARRAY<int> subset;
    ARRAY<int> subset_triangles;
    ARRAY<int> subset_particles;
    bool boundary_only;
    ARRAY<T> spectrum;
    bool draw_subsets;
    int cutaway_mode;
    T cutaway_fraction;
    TETRAHEDRON_MESH cutaway_mesh;
    ARRAY<OPENGL_COLOR>* color_map;
protected:
    bool smooth_normals;
    ARRAY<VECTOR<T,3> >* vertex_normals;
    OPENGL_SELECTION* current_selection;
public:

    OPENGL_TETRAHEDRALIZED_VOLUME(const OPENGL_MATERIAL& material_input)
        :material(material_input),mesh(0),particles(0),current_tetrahedron(1),current_node(1),current_boundary_triangle(1),boundary_only(true),
        draw_subsets(false),cutaway_mode(0),cutaway_fraction((T).5),color_map(0),smooth_normals(false),vertex_normals(0),current_selection(0)
    {
        Initialize();
    }

    OPENGL_TETRAHEDRALIZED_VOLUME(TETRAHEDRON_MESH* mesh_input,GEOMETRY_PARTICLES<VECTOR<T,3> >* particles_input,const OPENGL_MATERIAL& material_input,bool initialize=true,
        ARRAY<OPENGL_COLOR>* color_map_input=0)
        :material(material_input),mesh(mesh_input),particles(particles_input),current_tetrahedron(1),current_node(1),current_boundary_triangle(1),boundary_only(true),
        draw_subsets(false),cutaway_mode(0),cutaway_fraction((T).5),color_map(color_map_input),smooth_normals(false),vertex_normals(0),current_selection(0)
    {
        if(initialize)Initialize();
    }

    ~OPENGL_TETRAHEDRALIZED_VOLUME()
    {delete vertex_normals;}

    void Initialize()
    {if(!mesh->boundary_mesh) mesh->Initialize_Boundary_Mesh(); // Neighboring nodes is no longer initialized here to conserve memory.
    if(!mesh->node_on_boundary) mesh->Initialize_Node_On_Boundary();minimum_valence=mesh->Minimum_Valence();mesh->Initialize_Boundary_Nodes();}

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION* Get_Selection(GLuint* buffer,int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION* selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;

    OPENGL_SELECTION* Get_Vertex_Selection(int index);
    OPENGL_SELECTION* Get_Tetrahedron_Selection(int index);

    void Set_Boundary_Only(bool boundary_only_input)
    {boundary_only=boundary_only_input;}

    void Draw_Boundary_Only()
    {boundary_only=true;}

    void Draw_Interior()
    {boundary_only=false;}

    void Set_Color_Map(ARRAY<OPENGL_COLOR>* color_map_input)
    {color_map=color_map_input;}

    bool Toggle_Boundary_Only()
    {boundary_only=!boundary_only;return boundary_only;}

    bool Is_Smooth_Normals()
    {return smooth_normals;}

    void Cycle_Cutaway_Mode()
    {cutaway_mode=(cutaway_mode+1)%8;if(cutaway_mode)Update_Cutaway_Plane();}

    void Decrease_Cutaway_Fraction()
    {cutaway_fraction-=(T).05;if(cutaway_fraction<0)cutaway_fraction+=1;if(cutaway_mode)Update_Cutaway_Plane();}

    void Increase_Cutaway_Fraction()
    {cutaway_fraction+=(T).05;if(cutaway_fraction>1)cutaway_fraction-=1;if(cutaway_mode)Update_Cutaway_Plane();}

//#####################################################################
    void Draw_Boundary_Triangles(const TETRAHEDRON_MESH& tetrahedron_mesh) const;
    void Draw_Current_Tetrahedron() const;
    void Draw_Subset() const;
    void Draw_Subset_Triangles() const;
    void Draw_Subset_Particles() const;
    void Draw_In_Color_From_Spectrum() const;
    void Draw_In_Color_From_Color_Map() const;
    void Draw_Wireframe_Mesh(const TETRAHEDRON_MESH& tetrahedron_mesh) const;
    void Highlight_Boundary_Nodes_Of_Current_Tetrahedron() const;
    void Highlight_Current_Node() const;
    void Highlight_Nodes_Of_Minimum_Valence() const;
    void Highlight_Boundary_Normal_Vectors_Of_Current_Tetrahedron() const;
    void Highlight_Current_Boundary_Triangle() const;
    void Turn_Smooth_Shading_On() PHYSBAM_OVERRIDE;
    void Turn_Smooth_Shading_Off() PHYSBAM_OVERRIDE;
    void Display_Subset();
    void Update_Cutaway_Plane();
    void Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION *selection) const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION* selection,MATRIX<T,4>* transform) const;
    void Initialize_Vertex_Normals();
protected:
    void Draw_Vertices_For_Selection() const;
    void Draw_Tetrahedra_For_Selection() const;
    int Find_Shortest_Spring(const VECTOR<int,4>& nodes,T& distance,TV& minimum_normal,TV& weights) const;
//#####################################################################
};

template<class T>
class OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_VERTEX:public OPENGL_SELECTION
{
public:
    int index;
    OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_VERTEX(OPENGL_OBJECT* object,int index=0) 
        :OPENGL_SELECTION(OPENGL_SELECTION::TETRAHEDRALIZED_VOLUME_VERTEX,object),index(index)
    {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_TETRAHEDRON:public OPENGL_SELECTION
{
public:
    int index;
    OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_TETRAHEDRON(OPENGL_OBJECT* object,int index=0) 
        :OPENGL_SELECTION(OPENGL_SELECTION::TETRAHEDRALIZED_VOLUME_TETRAHEDRON,object),index(index)
    {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};
}
#endif

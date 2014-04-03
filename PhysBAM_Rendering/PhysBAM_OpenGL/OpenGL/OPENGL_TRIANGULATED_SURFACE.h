//#####################################################################
// Copyright 2002-2004, Eran Guendelman, Eilene Hao, Neil Molino, Robert Bridson, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_TRIANGULATED_SURFACE
//##################################################################### 
#ifndef __OPENGL_TRIANGULATED_SURFACE__
#define __OPENGL_TRIANGULATED_SURFACE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T> class TRIANGULATED_SURFACE;
template<class T> class OPENGL_SELECTION_TRIANGULATED_SURFACE_VERTEX;
template<class T> class OPENGL_SELECTION_TRIANGULATED_SURFACE_SEGMENT;
template<class T> class OPENGL_SELECTION_TRIANGULATED_SURFACE_TRIANGLE;

template<class T>
class OPENGL_TRIANGULATED_SURFACE:public OPENGL_OBJECT
{
public:
    OPENGL_TRIANGULATED_SURFACE(TRIANGULATED_SURFACE<T>& surface_input,bool smooth_normals_input=true,
                                const OPENGL_MATERIAL& material_input=OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Cyan()));
    OPENGL_TRIANGULATED_SURFACE(TRIANGULATED_SURFACE<T>& surface_input,bool smooth_normals_input,
                                const OPENGL_MATERIAL& front_material_input,const OPENGL_MATERIAL& back_material_input);
    ~OPENGL_TRIANGULATED_SURFACE();

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION *selection) const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION* selection,MATRIX<T,4>* transform) const;

    OPENGL_SELECTION *Get_Vertex_Selection(int index);
    OPENGL_SELECTION *Get_Segment_Selection(int index);
    OPENGL_SELECTION *Get_Triangle_Selection(int index);

    void Turn_Smooth_Shading_On() PHYSBAM_OVERRIDE;
    void Turn_Smooth_Shading_Off() PHYSBAM_OVERRIDE;
    bool Is_Smooth_Normals() {return smooth_normals;}
    void Set_Two_Sided(bool two_sided_input=true) {two_sided=two_sided_input;}
    void Set_Front_Material(const OPENGL_MATERIAL &material_input);
    void Set_Back_Material(const OPENGL_MATERIAL &material_input);
    void Initialize_Vertex_Normals();
    void Delete_Vertex_Normals();
    void Print_Triangles_Incident_On_Current_Node();
    void Draw_Triangles_Incident_On_Current_Node() const;
    void Print_Neighbor_Nodes_Of_Current_Node();

    // Create_Display_Lists creates a display list which this object will own.
    // Use_Display_List makes this object use a display list owned by another object.
    int Create_Display_List();
    int Get_Display_List_Id() { return display_list_id; }
    void Use_Display_List(int input_display_list_id);

    void Highlight_Current_Node() const;

    void Rescale(const T scaling_x,const T scaling_y,const T scaling_z)
    {
        surface.Rescale(scaling_x,scaling_y,scaling_z);        
        if (owns_display_list) Reinitialize_Display_List();
    }
    
public:
    void Use_Vertex_Colors();
    void Set_Vertex_Color(const int i,const OPENGL_COLOR color);
    void Reinitialize_Display_List();
    void Draw() const;
    void Draw_Subsets() const;
    void Draw_Vertices_For_Selection() const;
    void Draw_Segments_For_Selection() const;
    void Draw_Triangles_For_Selection() const;

    TRIANGULATED_SURFACE<T>& surface;
    bool two_sided;
public:
    OPENGL_MATERIAL front_material, front_material_gray;
    OPENGL_MATERIAL back_material, back_material_gray;

    ARRAY<VECTOR<T,3> >* vertex_normals;
    ARRAY<OPENGL_COLOR>* vertex_colors;
protected:
    bool smooth_normals;

    bool use_display_list;
    bool owns_display_list;
    int display_list_id;

public:
    OPENGL_SELECTION *current_selection;
    int current_node;
    bool highlight_current_node;
    bool highlight_neighbors_of_current_node;
    bool highlight_boundary;
    bool wireframe_only;
    bool draw_subsets;
    bool draw_velocities;
    bool draw_particles;
    T velocity_scale;

    ARRAY<int> subset;
    ARRAY<int> subset_particles;

    friend class OPENGL_SELECTION_TRIANGULATED_SURFACE_VERTEX<T>;
    friend class OPENGL_SELECTION_TRIANGULATED_SURFACE_SEGMENT<T>;
    friend class OPENGL_SELECTION_TRIANGULATED_SURFACE_TRIANGLE<T>;
};

template<class T>
class OPENGL_SELECTION_TRIANGULATED_SURFACE_VERTEX : public OPENGL_SELECTION
{
public:
    int index;
    OPENGL_SELECTION_TRIANGULATED_SURFACE_VERTEX(OPENGL_OBJECT *object, int index=0) 
        : OPENGL_SELECTION(OPENGL_SELECTION::TRIANGULATED_SURFACE_VERTEX, object), index(index) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_TRIANGULATED_SURFACE_SEGMENT : public OPENGL_SELECTION
{
public:
    int index;
    OPENGL_SELECTION_TRIANGULATED_SURFACE_SEGMENT(OPENGL_OBJECT *object, int index=0) 
        : OPENGL_SELECTION(OPENGL_SELECTION::TRIANGULATED_SURFACE_SEGMENT, object), index(index) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_TRIANGULATED_SURFACE_TRIANGLE : public OPENGL_SELECTION
{
public:
    int index;
    OPENGL_SELECTION_TRIANGULATED_SURFACE_TRIANGLE(OPENGL_OBJECT *object, int index=0) 
        : OPENGL_SELECTION(OPENGL_SELECTION::TRIANGULATED_SURFACE_TRIANGLE, object), index(index) {}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}
#endif

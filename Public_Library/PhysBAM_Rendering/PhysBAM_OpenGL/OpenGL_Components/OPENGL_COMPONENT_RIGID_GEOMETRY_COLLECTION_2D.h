//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Tamar Shinar, Rachel Weinstein, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D__
#define __OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T> class OPENGL_SEGMENTED_CURVE_2D;
template<class T> class OPENGL_TRIANGULATED_AREA;
template<class T> class OPENGL_LEVELSET_2D;
template<class T> class OPENGL_AXES;

template<class T,class RW=T>
class OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D : public OPENGL_COMPONENT
{
protected:
    typedef VECTOR<T,2> TV;

    std::string basedir;
    int frame_loaded;
    bool valid;
    bool show_object_names;
    bool output_positions;
    bool draw_velocity_vectors;
    bool draw_individual_axes;
    bool draw_node_velocity_vectors;
    bool draw_segmented_curve;
    bool draw_triangulated_area;
    bool draw_implicit_curve;
    bool has_init_destroy_information;
    ARRAY<int> needs_init,needs_destroy;

public:
    RIGID_GEOMETRY_COLLECTION<TV>* rigid_geometry_collection;
    ARRAY<int> colors;
protected:
    ARRAY<OPENGL_SEGMENTED_CURVE_2D<T>*,int> opengl_segmented_curve;
    ARRAY<OPENGL_TRIANGULATED_AREA<T>*,int> opengl_triangulated_area;
    ARRAY<OPENGL_LEVELSET_2D<T>*,int> opengl_levelset;
    ARRAY<ARRAY<OPENGL_COMPONENT*>,int> extra_components;
    ARRAY<OPENGL_AXES<T>*,int> opengl_axes;
    ARRAY<bool,int> draw_object;
    ARRAY<bool,int> use_object_bounding_box;
    ARRAY<VECTOR<T,2> > positions;
    ARRAY<VECTOR<T,2> > velocity_vectors;
    ARRAY<VECTOR<T,2> > node_positions;
    ARRAY<VECTOR<T,2> > node_velocity_vectors;
    OPENGL_VECTOR_FIELD_2D<ARRAY<TV> > velocity_field;
    OPENGL_VECTOR_FIELD_2D<ARRAY<TV> > node_velocity_field;
    OPENGL_SELECTION* current_selection;
    bool need_destroy_rigid_geometry_collection;

public:
    OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D(RIGID_GEOMETRY_PARTICLES<TV>* particles_input,const std::string& basedir);
    OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D(RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection,const std::string& basedir);
    ~OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D();
    
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION *selection) const PHYSBAM_OVERRIDE;

    void Read_Hints(const std::string& filename);

    void Set_Draw_Object(int i, bool draw_it);  // Need to call Reinitialize after changing draw objects
    bool Get_Draw_Object(int i) const;
    void Set_Object_Color(int i, const OPENGL_COLOR &color);
    void Set_Use_Object_Bounding_Box(int i, bool use_it);
    void Set_Vector_Size(double size);

    void Toggle_Velocity_Vectors();
    void Toggle_Individual_Axes();
    void Toggle_Output_Positions();
    void Toggle_Show_Object_Names();
    void Toggle_Node_Velocity_Vectors();
    void Toggle_Draw_Mode();
    void Increase_Vector_Size();
    void Decrease_Vector_Size();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D, Toggle_Velocity_Vectors, "Toggle velocity vectors");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D, Toggle_Individual_Axes, "Toggle individual axes");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D, Toggle_Output_Positions, "Toggle output positions");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D, Toggle_Show_Object_Names, "Toggle show object names");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D, Toggle_Node_Velocity_Vectors, "Toggle node velocity vectors");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D, Toggle_Draw_Mode, "Toggle draw mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D, Increase_Vector_Size, "Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D, Decrease_Vector_Size, "Decrease vector size");

public:
    virtual void Reinitialize(const bool force=false,const bool read_geometry=true);    // Needs to be called after some state changes
protected:
    void Set_Draw_Mode(const int mode);
    int Get_Draw_Mode() const;
    void Create_Geometry(const int id);
    void Update_Geometry(const int id);
    void Destroy_Geometry(const int id);
    virtual void Update_Object_Labels();
};

template<class T>
class OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D : public OPENGL_SELECTION
{
public:
    int body_id;
    OPENGL_SELECTION *body_selection;

    OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D(OPENGL_OBJECT *object,const int body_id,OPENGL_SELECTION* body_selection)
        :OPENGL_SELECTION(OPENGL_SELECTION::COMPONENT_RIGID_BODIES_2D,object),body_id(body_id),body_selection(body_selection)
    {}

    virtual OPENGL_SELECTION::TYPE Actual_Type() const 
    {return body_selection->Actual_Type();}

    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}

#endif

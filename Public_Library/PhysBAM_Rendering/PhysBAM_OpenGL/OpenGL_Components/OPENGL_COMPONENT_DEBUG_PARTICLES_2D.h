//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_DEBUG_PARTICLES_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_DEBUG_PARTICLES_2D__
#define __OPENGL_COMPONENT_DEBUG_PARTICLES_2D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_DEBUG_PARTICLES_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{

template<class T,class RW=T>
class OPENGL_COMPONENT_DEBUG_PARTICLES_2D:public OPENGL_COMPONENT
{
    typedef VECTOR<T,2> TV;
public:
    OPENGL_COMPONENT_DEBUG_PARTICLES_2D(const std::string &filename);
    ~OPENGL_COMPONENT_DEBUG_PARTICLES_2D();

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
    OPENGL_SELECTION* Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection) PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Selection_Bounding_Box(OPENGL_SELECTION *selection) const PHYSBAM_OVERRIDE;
    void Toggle_Draw_Velocities();
    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Arrowhead();
    void Command_Prompt();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEBUG_PARTICLES_2D,Toggle_Draw_Velocities,"Toggle draw velocities");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEBUG_PARTICLES_2D,Increase_Vector_Size,"Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEBUG_PARTICLES_2D,Decrease_Vector_Size,"Decrease vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEBUG_PARTICLES_2D,Toggle_Arrowhead,"Toggle arrow heads");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEBUG_PARTICLES_2D,Command_Prompt,"Command prompt");

private:
    void Reinitialize(bool force=false);

    void Command_Prompt_Response();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEBUG_PARTICLES_2D, Command_Prompt_Response, "");

public:
    GEOMETRY_PARTICLES<TV>& particles;
    OPENGL_DEBUG_PARTICLES_2D<T>& opengl_particles;

private:
    std::string filename;
    int frame_loaded;
    int set;
    int set_loaded;
    bool valid;
    bool draw_multiple_particle_sets;
};

template<class T>
class OPENGL_SELECTION_COMPONENT_DEBUG_PARTICLES_2D : public OPENGL_SELECTION
{
public:
    int index;  // index into particles array
    VECTOR<T,2> location;

    OPENGL_SELECTION_COMPONENT_DEBUG_PARTICLES_2D(OPENGL_OBJECT *object) : OPENGL_SELECTION(OPENGL_SELECTION::DEBUG_PARTICLES_2D, object) {}
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}

#endif

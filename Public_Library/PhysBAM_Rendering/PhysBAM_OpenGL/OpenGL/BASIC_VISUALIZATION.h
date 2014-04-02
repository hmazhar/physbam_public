//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BASIC_VISUALIZATION
//#####################################################################
#ifndef __BASIC_VISUALIZATION__
#define __BASIC_VISUALIZATION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_BASIC_CALLBACKS.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
class PARSE_ARGS;
class OPENGL_COMPONENT;

class BASIC_VISUALIZATION
{
public:
    static const int OWNED;
    static const int SELECTABLE;
    static const int START_HIDDEN;
    
             BASIC_VISUALIZATION();
    virtual ~BASIC_VISUALIZATION();

    void    Initialize(int argc, char *argv[]);
    void    Run();
    void    Initialize_And_Run(int argc, char *argv[]);

    virtual void Process_Hits(GLint hits, GLuint buffer[]);

    DEFINE_CALLBACK_CREATOR(BASIC_VISUALIZATION, Draw_All_Objects);
protected:
    virtual void Add_Arguments(PARSE_ARGS &parse_args);
    virtual void Parse_Arguments(PARSE_ARGS &parse_args);
    virtual void Initialize_Components_And_Key_Bindings();
    virtual void Add_OpenGL_Initialization();
    virtual void Initialize_Scene();
    virtual void Update_OpenGL_Strings();
    virtual void Goto_Start_Frame() {}
    virtual void Render_Offscreen() {}
    virtual void Selection_Callback();

    void Add_Component(OPENGL_COMPONENT *component, const std::string &name,const char toggle_draw_key,const int flags);
    const OPENGL_COMPONENT* Find_Component(const std::string& name) const;
    OPENGL_COMPONENT* Find_Component(const std::string& name);
    void Set_Current_Selection(OPENGL_SELECTION *selection);
    int &Selection_Priority(OPENGL_SELECTION::TYPE selection_type);

private:
    void Parse_Args(int argc, char *argv[]);
    void PreInitialize_OpenGL_World();
    void PostInitialize_OpenGL_World();
    void Reset_Objects_In_World();

private:
    void Reset_View();
    void Reset_Up();
    void Toggle_Axes();
    void Draw_All_Objects();

    DEFINE_CALLBACK_CREATOR(BASIC_VISUALIZATION, Reset_View);
    DEFINE_CALLBACK_CREATOR(BASIC_VISUALIZATION, Reset_Up);
    DEFINE_CALLBACK_CREATOR(BASIC_VISUALIZATION, Toggle_Axes);

protected:
    ARRAY<OPENGL_COMPONENT*> component_list;
    ARRAY<OPENGL_COMPONENT*> owned_components;
    HASHTABLE<std::string,OPENGL_COMPONENT*> component_by_name;

    OPENGL_WORLD            opengl_world;
    int                     width, height;
    bool                    set_window_position;
    VECTOR<int,2>          window_position;
    float                   fovy;
    std::string             opengl_window_title;
    bool                    add_axes;
    bool                    render_offscreen;
    std::string             camera_script_filename;
    std::string             initialization_key_sequence;

    // Selection stuff
    bool                                selection_enabled;
    OPENGL_SELECTION                    *current_selection;
    ARRAY<int>                     selection_priority;     // higher priority takes precedence; priority=0 is unselectable
};

}

#endif

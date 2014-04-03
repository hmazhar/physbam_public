//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT
//#####################################################################
#ifndef __OPENGL_COMPONENT__
#define __OPENGL_COMPONENT__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <string>

#define DEFINE_COMPONENT_CALLBACK(classname, func, str) \
    class classname##func : public OPENGL_CALLBACK { \
        public: \
            classname##func(classname *obj) : obj(obj) {} \
            void operator() () { obj->func(); } \
            void Print(std::ostream &out) { if (!obj->component_name.empty()) out << obj->component_name << ": " << str; else out << str; } \
        private: \
             classname *obj;}; \
    OPENGL_CALLBACK *func##_CB() \
    { return new classname##func(this); }

namespace PhysBAM
{

class OPENGL_COMPONENT : public OPENGL_OBJECT
{
public:
    OPENGL_COMPONENT(const std::string &name = "");
    virtual ~OPENGL_COMPONENT();

    void Set_Name(const std::string &name) { component_name = name; }

    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual bool Valid_Frame(int frame_input) const;
    virtual bool Is_Up_To_Date(int frame) const;

    virtual void Set_Frame(int frame_input);
    virtual void Set_Draw(bool draw_input = true);
    virtual void Draw_All_Objects();

    bool Is_Animated() { return is_animation; }
    
    virtual void Next_Frame();
    virtual void Prev_Frame();
    void Toggle_Draw() { Set_Draw(!draw); }

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT, Next_Frame, "Next frame");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT, Prev_Frame, "Prev frame");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT, Toggle_Draw, "Toggle draw");

protected:
    int     frame;
    bool    draw;
    bool    is_animation;

public:
    std::string component_name;
};

}

#endif

//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_BASIC
//#####################################################################
#ifndef __OPENGL_COMPONENT_BASIC__
#define __OPENGL_COMPONENT_BASIC__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
template<class T>
class OPENGL_COMPONENT_BASIC : public OPENGL_COMPONENT
{
public:
    T  &object;
    bool own_object;
    bool use_clip_planes;

    OPENGL_COMPONENT_BASIC(T &object)
        :object(object),own_object(true)
    {
        Use_Clip_Planes(false);
    }

    ~OPENGL_COMPONENT_BASIC()
    {
        if(own_object) delete &object;
    }

    void Use_Clip_Planes(bool use=true) 
    {use_clip_planes=use;}

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE
    {if(draw){
        if(use_clip_planes){glEnable(GL_CLIP_PLANE0);glEnable(GL_CLIP_PLANE1);}
        object.Display(in_color);
        if(use_clip_planes){glDisable(GL_CLIP_PLANE0);glDisable(GL_CLIP_PLANE1);}}}

    bool Use_Bounding_Box() const  PHYSBAM_OVERRIDE
    {if(draw) return object.Use_Bounding_Box();else return false;}

    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE
    {if(draw) return object.Bounding_Box();else return RANGE<VECTOR<float,3> >::Centered_Box();}

    bool Is_Transparent() const PHYSBAM_OVERRIDE
    {return object.Is_Transparent();}

    virtual void Turn_Smooth_Shading_Off() PHYSBAM_OVERRIDE
    {object.Turn_Smooth_Shading_Off();}

    virtual void Turn_Smooth_Shading_On() PHYSBAM_OVERRIDE
    {object.Turn_Smooth_Shading_On();}

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer,int buffer_size)
    {return object.Get_Selection(buffer,buffer_size);}

    virtual void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE
    {object.Highlight_Selection(selection);}

    virtual void Clear_Highlight() PHYSBAM_OVERRIDE
    {object.Clear_Highlight();}

    virtual void Set_Slice(OPENGL_SLICE *slice_input) PHYSBAM_OVERRIDE
    {slice=slice_input;object.Set_Slice(slice_input);}

    virtual void Slice_Has_Changed() PHYSBAM_OVERRIDE
    {object.Slice_Has_Changed();}

    void Print_Selection_Info(std::ostream& ostream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE
    {object.Print_Selection_Info(ostream,selection);}
};
}

#endif

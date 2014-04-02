//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Eilene Hao, Geoffrey Irving, Michael Lentine, Neil Molino, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_OBJECT
//#####################################################################
#ifndef __OPENGL_OBJECT__
#define __OPENGL_OBJECT__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/Convert_1d_To_3d.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/Convert_2d_To_3d.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h> // just so we get gl in the right order
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SLICE.h>
namespace PhysBAM{

class OPENGL_SELECTION;

class OPENGL_OBJECT:public NONCOPYABLE
{
public:
    FRAME<VECTOR<float,3> >* frame; // pointer so you can enslave this body to another's motion
    std::string name;
    bool selectable;
    bool visible;
    bool show_name;
    OPENGL_SLICE *slice;
private:
    FRAME<VECTOR<float,3> > default_frame;

public:
    OPENGL_OBJECT();
    virtual ~OPENGL_OBJECT();

    void Set_Name(const std::string& name_input)
    {name=name_input;}

    void Enslave_Transform_To(OPENGL_OBJECT& object)
    {frame=object.frame;}

    virtual void Display(const int in_color=1) const;
    virtual bool Use_Bounding_Box() const;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const;
    virtual bool Is_Transparent() const;
    virtual void Turn_Smooth_Shading_On();
    virtual void Turn_Smooth_Shading_Off();

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer,int buffer_size);
    virtual void Set_Selection(OPENGL_SELECTION *selection);
    virtual void Highlight_Selection(OPENGL_SELECTION *selection);
    virtual void Clear_Highlight();
    virtual void Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION *selection) const;
    virtual RANGE<VECTOR<float,3> > Selection_Bounding_Box(OPENGL_SELECTION *selection) const;
    virtual OPENGL_SELECTION* Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection);

    virtual void Set_Slice(OPENGL_SLICE *slice_input);
    virtual void Slice_Has_Changed();

    VECTOR<float,3> World_Space_Point(const VECTOR<float,1>& object_space_point) const
    {return World_Space_Point(Convert_1d_To_3d(object_space_point));}

    VECTOR<float,3> World_Space_Point(const VECTOR<float,2>& object_space_point) const
    {return World_Space_Point(Convert_2d_To_3d(object_space_point));}

    VECTOR<float,3> World_Space_Point(const VECTOR<float,3>& object_space_point) const
    {return *frame*object_space_point;}

    RANGE<VECTOR<float,3> > World_Space_Box(const RANGE<VECTOR<float,1> >& object_space_box) const
    {return World_Space_Box(RANGE<VECTOR<float,3> >(object_space_box.min_corner.x,object_space_box.max_corner.x,0,0,0,0));}

    RANGE<VECTOR<float,3> > World_Space_Box(const RANGE<VECTOR<float,2> >& object_space_box) const
    {return World_Space_Box(Convert_2d_To_3d(object_space_box));}

    RANGE<VECTOR<float,3> > World_Space_Box(const RANGE<VECTOR<float,3> >& object_space_box) const
    {return ORIENTED_BOX<VECTOR<float,3> >(object_space_box,*frame).Axis_Aligned_Bounding_Box();}

    RANGE<VECTOR<float,3> > World_Space_Box(const RANGE<VECTOR<double,1> >& object_space_box) const
    {return World_Space_Box(RANGE<VECTOR<float,1> >(object_space_box));}

    RANGE<VECTOR<float,3> > World_Space_Box(const RANGE<VECTOR<double,2> >& object_space_box) const
    {return World_Space_Box(RANGE<VECTOR<float,2> >(object_space_box));}

    RANGE<VECTOR<float,3> > World_Space_Box(const RANGE<VECTOR<double,3> >& object_space_box) const
    {return World_Space_Box(RANGE<VECTOR<float,3> >(object_space_box));}

    void Send_Transform_To_GL_Pipeline() const
    {OpenGL_Translate(frame->t);OpenGL_Rotate(frame->r);}
};
}
#endif

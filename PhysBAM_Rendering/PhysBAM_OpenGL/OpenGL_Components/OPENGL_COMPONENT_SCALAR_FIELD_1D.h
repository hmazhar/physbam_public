//#####################################################################
// Copyright 2006-2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_SCALAR_FIELD_1D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_SCALAR_FIELD_1D__
#define __OPENGL_COMPONENT_SCALAR_FIELD_1D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <string>
namespace PhysBAM{

template<class T,class T2=T,class RW=T>
class OPENGL_COMPONENT_SCALAR_FIELD_1D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,1> TV;
private:
    std::string scalar_field_filename;
    int frame_loaded;
    bool valid;
public:
    OPENGL_SCALAR_FIELD_1D<T,T2>  opengl_scalar_field;

    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE
    {return valid && frame_loaded==frame;}

    bool Use_Bounding_Box() const  PHYSBAM_OVERRIDE
    {return draw && valid;}

//##################################################################### 
    OPENGL_COMPONENT_SCALAR_FIELD_1D(const GRID<TV> &grid_input,const std::string &scalar_field_filename_input,OPENGL_COLOR point_color,OPENGL_COLOR line_color);
    ~OPENGL_COMPONENT_SCALAR_FIELD_1D();
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Scale(const T scale);
    void Increase_Scale();
    void Decrease_Scale();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SCALAR_FIELD_1D,Increase_Scale,"Increase scale");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SCALAR_FIELD_1D,Decrease_Scale,"Decrease scale");
private:
    void Reinitialize();
//##################################################################### 
};
}
#endif

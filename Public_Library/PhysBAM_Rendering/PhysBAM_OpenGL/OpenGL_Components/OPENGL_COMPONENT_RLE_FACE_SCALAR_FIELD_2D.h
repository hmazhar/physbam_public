#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_2D__
#define __OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_2D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_FACE_SCALAR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T,class T2=T,class RW=T>
class OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_2D:public OPENGL_COMPONENT
{
public:
    OPENGL_RLE_FACE_SCALAR_FIELD_2D<T,T2> opengl_rle_face_scalar_field;
private:
    std::string scalar_field_filename;
    int frame_loaded;
    bool valid;
public:

    OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_2D(RLE_GRID_2D<T>& grid,const std::string& scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map,bool draw_points);
    ~OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_2D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE {return valid && frame_loaded==frame;}

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE {return draw && valid;}
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;

    void Increase_Point_Size();
    void Decrease_Point_Size();

    void Increase_Vector_Size()
    {opengl_rle_face_scalar_field.Scale_Vector_Size(1.1);}

    void Decrease_Vector_Size()
    {opengl_rle_face_scalar_field.Scale_Vector_Size(1/1.1);}

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_2D, Increase_Point_Size, "Increase point size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_2D, Decrease_Point_Size, "Decrease point size");

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_2D, Increase_Vector_Size, "Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_2D, Decrease_Vector_Size, "Decrease vector size");
private:
    void Reinitialize();

};
}
#endif
#endif

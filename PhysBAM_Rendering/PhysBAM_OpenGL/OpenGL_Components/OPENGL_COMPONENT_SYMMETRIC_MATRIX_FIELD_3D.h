//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D__
#define __OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D__

#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SYMMETRIC_MATRIX_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T,class RW=T>
class OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D:public OPENGL_COMPONENT
{
    typedef VECTOR<T,3> TV;
public:
    ARRAY<SYMMETRIC_MATRIX<T,3>,VECTOR<int,3> > field;
    OPENGL_SYMMETRIC_MATRIX_FIELD_3D<T> opengl_symmetric_matrix_field;
    std::string field_filename;
    int frame_loaded;
    bool valid;

    OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D(const GRID<TV>& grid,const std::string& field_filename_input)
        :OPENGL_COMPONENT("Symmetric Matrix Field 3D"),opengl_symmetric_matrix_field(grid,field),
        field_filename(field_filename_input),frame_loaded(-1),valid(false)
    {
        is_animation=FILE_UTILITIES::Is_Animated(field_filename);Reinitialize();
    }

    ~OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D()
    {}

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE
    {return FILE_UTILITIES::Frame_File_Exists(field_filename,frame_input);}

    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE
    {return valid && frame_loaded==frame;}

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE
    {OPENGL_COMPONENT::Set_Frame(frame_input);Reinitialize();}

    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE
    {OPENGL_COMPONENT::Set_Draw(draw_input);Reinitialize();}

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE
    {if(valid&&draw)opengl_symmetric_matrix_field.Display(in_color);}
    
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE
    {return draw&&valid;}

    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE
    {if(valid&&draw)return opengl_symmetric_matrix_field.Bounding_Box();return RANGE<VECTOR<float,3> >::Centered_Box();}

    void Set_Slice(OPENGL_SLICE* slice_input) PHYSBAM_OVERRIDE
    {slice=slice_input;opengl_symmetric_matrix_field.Set_Slice(slice_input);}

    void Slice_Has_Changed() PHYSBAM_OVERRIDE
    {opengl_symmetric_matrix_field.Slice_Has_Changed();} 

    void Increase_Size()
    {opengl_symmetric_matrix_field.size*=(T)1.1;}

    void Decrease_Size()
    {opengl_symmetric_matrix_field.size*=1/(T)1.1;}

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D,Increase_Size,"Increase symmetric matrix size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D,Decrease_Size,"Decrease symmetric matrix size");

//#####################################################################
    void Reinitialize(bool force_load_even_if_not_drawn=false);
//#####################################################################
};
}

#endif

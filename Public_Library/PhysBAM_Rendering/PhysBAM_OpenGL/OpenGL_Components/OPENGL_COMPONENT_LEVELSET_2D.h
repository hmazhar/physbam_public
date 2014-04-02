//#####################################################################
// Copyright 2004, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_LEVELSET_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_LEVELSET_2D__
#define __OPENGL_COMPONENT_LEVELSET_2D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_LEVELSET_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T,class RW=T>
class OPENGL_COMPONENT_LEVELSET_2D:public OPENGL_COMPONENT
{
    typedef VECTOR<T,2> TV;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
public:
    OPENGL_COMPONENT_LEVELSET_2D(const std::string& levelset_filename_input,const std::string filename_set_input="");
    ~OPENGL_COMPONENT_LEVELSET_2D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE {return valid && frame_loaded == frame;}

    void Toggle_Color_Mode();
    void Toggle_Smooth();
    void Toggle_Normals();
    void Toggle_Draw_Mode();
    void Toggle_Draw_Sign();
    void Next_Set();
    void Previous_Set();
    void Toggle_Draw_Multiple_Levelsets();
    void Toggle_Draw_Ghost_Values();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_LEVELSET_2D,Toggle_Color_Mode,"Toggle color mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_LEVELSET_2D,Toggle_Smooth,"Toggle smooth levelset draw");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_LEVELSET_2D,Toggle_Normals,"Toggle levelset normals draw");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_LEVELSET_2D,Toggle_Draw_Mode,"Toggle levelset contour/cellview");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_LEVELSET_2D,Toggle_Draw_Sign,"Toggle levelset sign direction");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_LEVELSET_2D,Next_Set,"Switch to next set");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_LEVELSET_2D,Previous_Set,"Switch to previous set");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_LEVELSET_2D,Toggle_Draw_Multiple_Levelsets,"Toggle mutliple/single levelset draw");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_LEVELSET_2D, Toggle_Draw_Ghost_Values, "Toggle draw ghost values");
    
private:
    void Reinitialize(const bool force_even_if_not_drawn=false);
    template<class,class> friend class OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D;

public:
    OPENGL_LEVELSET_2D<T>* opengl_levelset;
    ARRAY<OPENGL_LEVELSET_2D<T>*> opengl_levelsets;

private:
    std::string levelset_filename;
    std::string filename_set;
    int frame_loaded;
    int set;
    bool use_sets;
    int set_loaded;
    bool valid;
    bool draw_multiple_levelsets;
};

}

#endif

//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_SCALAR_FIELD_3D
//#####################################################################
#ifndef __OPENGL_SCALAR_FIELD_3D__
#define __OPENGL_SCALAR_FIELD_3D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>

namespace PhysBAM
{
class OPENGL_TEXTURED_RECT;
template<class T,class T_ARRAY> class OPENGL_POINTS_3D;

template<class T,class T2=T>
class OPENGL_SCALAR_FIELD_3D : public OPENGL_OBJECT
{
    typedef VECTOR<T,3> TV;
public:
    GRID<TV> grid;
    ARRAY<T2,VECTOR<int,3> > &values;

    enum DRAW_MODE { DRAW_TEXTURE, DRAW_POINTS };
    DRAW_MODE draw_mode;
    ARRAY<OPENGL_COLOR_MAP<T2>*> color_maps; // all owned by us
    int current_color_map;

    OPENGL_SCALAR_FIELD_3D(const GRID<TV> &grid_input,ARRAY<T2,VECTOR<int,3> > &values_input,OPENGL_COLOR_MAP<T2> *color_map_input,DRAW_MODE draw_mode_input=DRAW_TEXTURE);
    ~OPENGL_SCALAR_FIELD_3D();
    void Initialize_Color_Maps(OPENGL_COLOR_MAP<T2>* color_map);

    void Set_Scale_Range(const T2 range_min,const T2 range_max);
    void Reset_Scale_Range();
    T2 Pre_Map_Value(const T2 value) const;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Slice_Has_Changed() PHYSBAM_OVERRIDE;    

    void Set_Draw_Mode(DRAW_MODE draw_mode);
    virtual void Update();  // Call when values or other attributes have changed
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;

    // convenience functions
    void Toggle_Draw_Mode();
    void Set_Smooth_Slice_Texture(bool smooth_slice_texture_input=true);
    void Toggle_Smooth_Slice_Texture();
    void Toggle_Color_Map();

private:
    void Display_3D() const;
    void Update_Slice();
    void Update_Points();
    void Delete_Textured_Rect();
    void Delete_Points();

public:
    OPENGL_TEXTURED_RECT *opengl_textured_rect;
    OPENGL_POINTS_3D<T,ARRAY<VECTOR<T,3> > > *opengl_points;
    bool smooth_slice_texture;
    bool scale_range;
    T2 scale_range_min,scale_range_dx;
};

}

#endif

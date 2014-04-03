//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_SCALAR_FIELD_2D
//#####################################################################
#ifndef __OPENGL_SCALAR_FIELD_2D__
#define __OPENGL_SCALAR_FIELD_2D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>

namespace PhysBAM
{
class OPENGL_TEXTURED_RECT;
template<class T,class T_ARRAY> class OPENGL_POINTS_2D;
template<class T> class OPENGL_SEGMENTED_CURVE_2D;

template<class T,class T2=T>
class OPENGL_SCALAR_FIELD_2D : public OPENGL_OBJECT
{
    typedef VECTOR<T,2> TV;
public:
    GRID<TV>& grid;
    ARRAY<T2,VECTOR<int,2> > &values;
    ARRAY<bool,VECTOR<int,2> > *active_cells;

    enum DRAW_MODE { DRAW_TEXTURE, DRAW_POINTS, DRAW_CONTOURS };
    DRAW_MODE draw_mode;
    bool draw_ghost_values;
    ARRAY<OPENGL_COLOR_MAP<T2>*> color_maps; // all owned by us
    int current_color_map;
private: 
    OPENGL_TEXTURED_RECT *opengl_textured_rect;
    OPENGL_POINTS_2D<T,ARRAY<VECTOR<T,2> > > *opengl_points;
    ARRAY<T2> contour_values;
    ARRAY<OPENGL_SEGMENTED_CURVE_2D<T>*> contour_curves;
    bool scale_range;
    T2 scale_range_min,scale_range_dx;

public:
    OPENGL_SCALAR_FIELD_2D(GRID<TV> &grid_input,ARRAY<T2,VECTOR<int,2> > &values_input,OPENGL_COLOR_MAP<T2>* color_map_input,DRAW_MODE draw_mode_input=DRAW_TEXTURE);
    OPENGL_SCALAR_FIELD_2D(GRID<TV> &grid_input,ARRAY<T2,VECTOR<int,2> > &values_input,OPENGL_COLOR_MAP<T2>* color_map_input,ARRAY<bool,VECTOR<int,2> >* active_cells_input,DRAW_MODE draw_mode_input=DRAW_TEXTURE);
    ~OPENGL_SCALAR_FIELD_2D();

    void Set_Scale_Range(const T2 range_min,const T2 range_max);
    void Reset_Scale_Range();
    T2 Pre_Map_Value(const T2 value) const;

    void Set_Uniform_Contour_Values(const T2 min_avlue,const T2 max_value,const T2 increment);

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    void Set_Draw_Mode(DRAW_MODE draw_mode);
    virtual void Update();  // Call when values or other attributes have changed

    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;

    // convenience functions
    void Toggle_Draw_Mode();
    void Toggle_Smooth_Texture();
    void Toggle_Draw_Ghost_Values();
    void Toggle_Color_Map();

    DEFINE_CALLBACK_CREATOR(OPENGL_SCALAR_FIELD_2D, Toggle_Draw_Ghost_Values);

private:
    void Initialize_Color_Maps(OPENGL_COLOR_MAP<T2>* color_map_input);
    void Update_Texture(const VECTOR<int,2>& start_index,const VECTOR<int,2>& end_index);
    void Update_Points(const VECTOR<int,2>& start_index,const VECTOR<int,2>& end_index);
    void Update_Contour_Curves();
};

}

#endif

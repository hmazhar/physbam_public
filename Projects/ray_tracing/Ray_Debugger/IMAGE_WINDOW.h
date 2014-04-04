//#####################################################################
// Copyright 2002, 2003, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __IMAGE_WINDOW__
#define __IMAGE_WINDOW__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h> // To get gl, glu in a portable manner

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <FL/Fl_Gl_Window.h>
namespace PhysBAM{

template<class T> 
class IMAGE_WINDOW : public Fl_Gl_Window
{
private:
    // camera setup
    T view_x_offset,view_y_offset,view_scale,view_scale_power;
    // callback data
    void (*user_callback)(void*,int,int);
    void *user_callback_data;
    bool have_last_click;VECTOR<int,2> last_click;
    // texture data
    bool update_image;
    int texture_dimension;
    GLuint image_texture;
    ARRAY<VECTOR<T,3> ,VECTOR<int,2> > image;
    bool texture_initialized;
public:
    IMAGE_WINDOW(int x,int y,int width,int height,int picture_width,int picture_height);
    void Set_Color(int i,int j,const VECTOR<T,3>& color);    
    void Update_Image();
    void Set_Callback(void (*)(void*,int,int),void *data);
    void Save_RGB(const char *filename);
    void Save_PhysBAM(const char *filename);
private:
    // fltk window draw/event handle interface functions
    int handle(int event);
    void draw();
    // helper functions
    void Init_Gl();
    void Print_Gl_Error();
    void Load_Camera();
};

}
#endif

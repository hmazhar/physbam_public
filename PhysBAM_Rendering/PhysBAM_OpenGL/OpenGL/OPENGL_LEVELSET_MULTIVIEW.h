//#####################################################################
// Copyright 2004-2008, Eran Guendelman, Michael Lentine, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_LEVELSET_MULTIVIEW
//#####################################################################
#ifndef __OPENGL_LEVELSET_MULTIVIEW__
#define __OPENGL_LEVELSET_MULTIVIEW__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>

namespace PhysBAM
{
template<class TV> class RANGE;
template<class TV> class GRID;
template<class T_GRID> class LEVELSET_3D;
template<class TV> class LEVELSET_IMPLICIT_OBJECT;
template<class T> class TRIANGULATED_SURFACE;
template<class T> class OPENGL_TRIANGULATED_SURFACE;
template<class T> class OPENGL_COLOR_MAP;
template<class T,class T2> class OPENGL_SCALAR_FIELD_3D;

template<class T,class RW=T>
class OPENGL_LEVELSET_MULTIVIEW:public OPENGL_OBJECT
{
    typedef VECTOR<T,3> TV;typedef LEVELSET_3D<GRID<TV> > T_LEVELSET;
public:
    enum COLOR_MODE {COLOR_SOLID,COLOR_GRADIENT}; // currently only used for slice mode
    COLOR_MODE color_mode;
    OPENGL_MATERIAL front_surface_material, back_surface_material;
    OPENGL_MATERIAL overlayed_surface_material;
    OPENGL_COLOR inside_slice_color, outside_slice_color;
    const FRAME<TV>* implicit_object_transform;
private:
    bool display_overlay;    // in slice mode
    bool smooth_shading;
    bool smooth_slice_texture;
    std::string levelset_filename;
    std::string triangulated_surface_filename;
    bool generate_triangulated_surface;
    bool write_generated_triangulated_surface;
    bool use_marching_tetrahedra;
    T_LEVELSET* levelset;
    LEVELSET_IMPLICIT_OBJECT<TV>* levelset_implicit_surface;
    TRIANGULATED_SURFACE<T>* triangulated_surface;
    bool i_own_levelset; 
    bool i_own_triangulated_surface;
    OPENGL_TRIANGULATED_SURFACE<T>* opengl_triangulated_surface;
    OPENGL_SCALAR_FIELD_3D<T,T>* opengl_scalar_field;
    bool two_sided;

public:
    OPENGL_LEVELSET_MULTIVIEW()
        :color_mode(COLOR_SOLID),implicit_object_transform(0),display_overlay(false),smooth_shading(true),smooth_slice_texture(false),
        use_marching_tetrahedra(false),levelset(0),levelset_implicit_surface(0),triangulated_surface(0),i_own_levelset(true),i_own_triangulated_surface(true),
        opengl_triangulated_surface(0),opengl_scalar_field(0),two_sided(true)
    {
        back_surface_material=OPENGL_MATERIAL::Plastic(OPENGL_COLOR((float)0.8,(float)0,0));
        front_surface_material=OPENGL_MATERIAL::Plastic(OPENGL_COLOR((float)0,0.6,(float)0.9));
        overlayed_surface_material=OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Gray((float)0.6,(float)0.5));
        outside_slice_color=OPENGL_COLOR::Red(0.5);
        inside_slice_color=OPENGL_COLOR::Blue();
        Reset();
    }

    ~OPENGL_LEVELSET_MULTIVIEW()
    {Reset();}

//#####################################################################
    // First call one of these
    void Set_Levelset(T_LEVELSET &levelset);
    void Read_Levelset(const std::string& levelset_filename);
    // Then (optionally) call one of these
    const TRIANGULATED_SURFACE<T>* Get_Triangulated_Surface() const;
    void Set_Triangulated_Surface(TRIANGULATED_SURFACE<T> &triangulated_surface);
    void Read_Triangulated_Surface(const std::string& triangulated_surface_filename);
    void Generate_Triangulated_Surface(bool write_generated_triangulatd_surface = false,const std::string& triangulated_surface_filename="");
    // Must call this next, and must be after you've created an OpenGL context
    void Initialize();
    void Reset();
    void Reset_Surface();
    void Set_Surface_Material(const OPENGL_MATERIAL &front_surface_mat,const OPENGL_MATERIAL &back_surface_mat);
    void Set_Overlayed_Surface_Material(const OPENGL_MATERIAL &overlayed_surface_mat);
    void Set_Slice_Color(const OPENGL_COLOR &inside_slice_color, const OPENGL_COLOR &outside_slice_color);
    void Set_Slice_Color_Mode(COLOR_MODE color_mode_input);
    void Set_Two_Sided(bool two_sided_input);
    void Toggle_Slice_Color_Mode();
    void Toggle_Display_Overlay();
    void Toggle_Smooth_Slice_Texture();
    void Toggle_Marching_Tetrahedra();
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Turn_Smooth_Shading_On() PHYSBAM_OVERRIDE;
    void Turn_Smooth_Shading_Off() PHYSBAM_OVERRIDE;
    void Slice_Has_Changed() PHYSBAM_OVERRIDE;
    const T_LEVELSET *Levelset() const;
    void Update();
private:
    void Initialize_Levelset();
    void Initialize_Triangulated_Surface();
    void Initialize_OpenGL_Triangulated_Surface();
    void Initialize_OpenGL_Scalar_Field();
//#####################################################################
};
}

#endif

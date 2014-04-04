//#####################################################################
// Copyright 2002-2005, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SCENE_WINDOW__
#define __SCENE_WINDOW__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_ARCBALL.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h> // To get gl, glu in a portable manner
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY_DEBUG.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE.h>
#include "FL/Fl_Gl_Window.h"

namespace PhysBAM{
template<class T> class TRIANGULATED_SURFACE;
template<class T> 
class SCENE_WINDOW : public Fl_Gl_Window
{
public:
    RENDER_WORLD<T>& world;
    RENDERING_RAY_DEBUG<T> *debug_ray,*selected_ray;
    RENDERING_OBJECT<T>* selected_object;
    ARRAY<TRIANGULATED_SURFACE<T>*> object_triangles;
    OPENGL_ARCBALL<GLfloat> arcball; 
    bool display_global_photon_map,display_caustic_photon_map,display_volume_photon_map,display_irradiance_cache,display_used_photons,display_used_irradiance_estimates,
        display_subsurface_scattering_samples;

    SCENE_WINDOW(int x,int y,int size_x,int size_y,RENDER_WORLD<T>& world)
        :Fl_Gl_Window(x,y,size_x,size_y,"Scene Window"),world(world),debug_ray(0),selected_ray(0),selected_object(0),arcball(*new OPENGL_WORLD()), 
         display_global_photon_map(false),display_caustic_photon_map(false),display_volume_photon_map(false),display_irradiance_cache(false),display_used_photons(false),
        display_used_irradiance_estimates(false),camera_rotation(false),camera_zoom(false),camera_translate(false)
    {
        Set_Look_At(VECTOR<GLfloat,3>(-4,2,0),VECTOR<GLfloat,3>(0,0,0),VECTOR<GLfloat,3>(0,1,0));
    }

    void Set_Debug_Ray(RENDERING_RAY_DEBUG<T>* debug_ray);
    void Select_Ray(RENDERING_RAY_DEBUG<T>* selected_debug_ray);
    void Select_Object(RENDERING_OBJECT<T>* select_object);
    void View_Global_Photon_Map(bool value){display_global_photon_map=value;damage(1);}
    void View_Caustic_Photon_Map(bool value){display_caustic_photon_map=value;damage(1);}
    void View_Volume_Photon_Map(bool value){display_volume_photon_map=value;printf("volume %s\n",display_volume_photon_map?"true":"false");damage(1);}
    void View_Irradiance_Cache(bool value){display_irradiance_cache=value;damage(1);}
    void View_Used_Photons(bool value){display_used_photons=value;damage(1);}
    void View_Used_Irradiance_Samples(bool value){display_used_irradiance_estimates=value;damage(1);}
    void View_Subsurface_Scattering_Samples(bool value){display_subsurface_scattering_samples=value;damage(1);}
    void Initialize_Objects();

private:
    // fltk window draw/event handle interface functions
    T camera_distance;
    VECTOR<GLfloat,3> camera_target;
    MATRIX<GLfloat,4> rotation_matrix;
    MATRIX<GLfloat,4> arcball_matrix;
    bool camera_rotation,camera_zoom,camera_translate;
    VECTOR<GLfloat,2> old_mouse_vector;
    VECTOR<GLfloat,2> Convert_Mouse_Coordinates(int x,int y);
    void Set_Look_At(const VECTOR<GLfloat,3> &camera, const VECTOR<GLfloat,3> &target, const VECTOR<GLfloat,3> &up);
    int handle(int event);
    void draw();
    void Init_Gl();
    void Draw_Ray(RENDERING_RAY_DEBUG<T>* debug_ray);
    void Draw_Photons(const PHOTON_MAP<T>& photon_map);
    void Draw_Subsurface_Scattering_Samples(SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<T>* bssrdf_tree);
    void Draw_Subsurface_Candidate_List();
    void Draw_Camera();
    void Draw_Photons();
    void Draw_Irradiance_Cache();
    VECTOR<GLfloat,3> Translate_In_Camera_Plane(VECTOR<GLfloat,2> old_mouse_vector,VECTOR<GLfloat,2> new_mouse_vector);
};

}
#endif

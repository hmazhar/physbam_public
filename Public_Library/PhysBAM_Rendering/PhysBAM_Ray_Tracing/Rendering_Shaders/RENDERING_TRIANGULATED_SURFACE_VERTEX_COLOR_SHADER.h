//#####################################################################
// Copyright 2005, Geoffrey Irving, Robert Travis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RENDERING_TRIANGULATED_SURFACE_VERTEX_COLOR_SHADER__
#define __RENDERING_TRIANGULATED_SURFACE_VERTEX_COLOR_SHADER__

#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/MATERIAL_SHADER.h>
namespace PhysBAM{

template<class T,class RW>
class RENDERING_TRIANGULATED_SURFACE_VERTEX_COLOR_SHADER:public MATERIAL_SHADER<T>
{
public:
    ARRAY<VECTOR<T,3> > vertex_colors;

    RENDERING_TRIANGULATED_SURFACE_VERTEX_COLOR_SHADER(RENDER_WORLD<T>& world_input)  
        :MATERIAL_SHADER<T>(world_input)
    {}

    void Initialize(const std::string& filename)
    {FILE_UTILITIES::Read_From_File<RW>(filename,vertex_colors);}

    VECTOR<T,3> Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
        const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
    {const TRIANGULATED_SURFACE<T>& triangulated_surface=((const RENDERING_TRIANGULATED_SURFACE<T>&) intersection_object).triangulated_surface; 
    int i,j,k;triangulated_surface.mesh.elements(ray.ray.aggregate_id).Get(i,j,k);
    VECTOR<T,3> barycenter=TRIANGLE_3D<T>::Barycentric_Coordinates(intersection_object.Object_Space_Point(intersection_point),triangulated_surface.particles.X(i),triangulated_surface.particles.X(j),triangulated_surface.particles.X(k));      
    return barycenter.x*vertex_colors(i)+barycenter.y*vertex_colors(j)+barycenter.z*vertex_colors(k);}

//#####################################################################
};
}
#endif

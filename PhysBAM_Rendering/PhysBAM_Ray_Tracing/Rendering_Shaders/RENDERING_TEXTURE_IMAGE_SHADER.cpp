//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Igor Neverov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_TEXTURE_IMAGE_SHADER  
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_MAC_GRID_PERIODIC.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/RENDERING_TEXTURE_IMAGE_SHADER.h>
using namespace PhysBAM;
//#####################################################################
// Function Setup_Interpolation
//#####################################################################
template<class T> void RENDERING_TEXTURE_IMAGE_SHADER<T>::
Setup_Interpolation(const int m,const int n)
{
    int ghost_cells=3;
    grid.Initialize(m,n,RANGE<VECTOR<T,2> >::Unit_Box(),true);pixels.Resize(grid.Domain_Indices(ghost_cells));
    if(wrap_s||wrap_t) BOUNDARY_MAC_GRID_PERIODIC<GRID<VECTOR<T,2> >,VECTOR<T,3> >().Fill_Ghost_Cells(grid,pixels,pixels,0,0,ghost_cells);
    else BOUNDARY_UNIFORM<GRID<VECTOR<T,2> >,VECTOR<T,3> >().Fill_Ghost_Cells(grid,pixels,pixels,0,0,ghost_cells);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void RENDERING_TEXTURE_IMAGE_SHADER<T>::
Initialize(const std::string& filename)
{
    IMAGE<T>::Read(filename,pixels);
    Setup_Interpolation(pixels.counts.x,pixels.counts.y);
}
//#####################################################################
// Function Shade_Surface_Using_Direct_Illumination
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_TEXTURE_IMAGE_SHADER<T>::
Shade_Surface_Using_Direct_Illumination(const RENDERING_RAY<T>& ray,const RENDERING_OBJECT<T>& exiting_object,const RENDERING_OBJECT<T>& entering_object,
     const RENDERING_OBJECT<T>& intersection_object,const VECTOR<T,3>& intersection_point,const VECTOR<T,3>& same_side_normal) const
{
    T s,t;intersection_object.Get_Texture_Coordinates(intersection_object.Object_Space_Point(intersection_point),ray.ray.aggregate_id,s,t);
    if(wrap_s) s-=floor(s);else clamp(s,(T)0,(T)1);
    if(wrap_t) t-=floor(t);else clamp(t,(T)0,(T)1);
    return interpolation->Clamped_To_Array(grid,pixels,VECTOR<T,2>(s,t));
}
//#####################################################################
template class RENDERING_TEXTURE_IMAGE_SHADER<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RENDERING_TEXTURE_IMAGE_SHADER<double>;
#endif

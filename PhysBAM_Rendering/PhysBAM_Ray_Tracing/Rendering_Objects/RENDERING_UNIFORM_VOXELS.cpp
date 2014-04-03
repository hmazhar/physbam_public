//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_UNIFORM_VOXELS
//#####################################################################
#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Computations/SMOOTH_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_UNIFORM_VOXELS.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/RENDERING_VOXEL_SHADER.h>
using namespace PhysBAM;
//#####################################################################
// Function Postprocess_Light_Field
//#####################################################################
template<class T> void RENDERING_UNIFORM_VOXELS<T>::
Postprocess_Light_Field()
{
    if(number_of_smoothing_steps) for(int light=1;light<=precomputed_light.m;light++){
        LOG::cout<<"Smoothing light "<<light<<" "<<number_of_smoothing_steps<<" steps"<<std::endl;
        SMOOTH::Smooth<GRID<TV> >(*precomputed_light(light),number_of_smoothing_steps,0);}
}
//#####################################################################
template class RENDERING_UNIFORM_VOXELS<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RENDERING_UNIFORM_VOXELS<double>;
#endif

//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T>::
OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD(TRIANGULATED_AREA<T>& triangulated_area,ARRAY<TV>& V)
    :OPENGL_VECTOR_FIELD_2D<ARRAY<TV> >(vector_field,vector_locations),triangulated_area(triangulated_area),V(V)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T>::
~OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD()
{}
//#####################################################################
// Function Update
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T>::
Update()
{
    V.Resize(triangulated_area.particles.array_collection->Size());
    vector_field.Resize(V.m);
    vector_locations.Resize(V.m);
    for(int i=1;i<=triangulated_area.particles.array_collection->Size();i++){
        vector_field(i)=V(i);vector_locations(i)=triangulated_area.particles.X(i);}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T>::
Bounding_Box() const
{
    return RANGE<VECTOR<float,3> >(Convert_2d_To_3d(triangulated_area.bounding_box?*triangulated_area.bounding_box:RANGE<TV>::Centered_Box()));
}
//#####################################################################
template class OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_TRIANGULATED_AREA_BASED_VECTOR_FIELD<double>;
#endif

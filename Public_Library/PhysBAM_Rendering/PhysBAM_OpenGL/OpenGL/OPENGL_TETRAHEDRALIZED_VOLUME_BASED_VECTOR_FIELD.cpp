//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T>::
OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume, ARRAY<VECTOR<T,3> >& V)
    :OPENGL_VECTOR_FIELD_3D<T>(*(new ARRAY<VECTOR<T,3> >),*(new ARRAY<VECTOR<T,3> >)),tetrahedralized_volume(tetrahedralized_volume),V(V)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T>::
~OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD()
{
    delete &vector_field;delete &vector_locations;
}
//#####################################################################
// Function Update
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T>::
Update()
{
    V.Resize(tetrahedralized_volume.particles.array_collection->Size());
    vector_field.Resize(V.m);
    vector_locations.Resize(V.m);
    for(int i=1;i<=tetrahedralized_volume.particles.array_collection->Size();i++){
        vector_field(i)=V(i);vector_locations(i)=tetrahedralized_volume.particles.X(i);}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T>::
Bounding_Box() const
{
    if(tetrahedralized_volume.bounding_box)
        return RANGE<VECTOR<float,3> >(RANGE<VECTOR<T,3> >(*tetrahedralized_volume.bounding_box));
    else
        return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
template class OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<double>;
#endif

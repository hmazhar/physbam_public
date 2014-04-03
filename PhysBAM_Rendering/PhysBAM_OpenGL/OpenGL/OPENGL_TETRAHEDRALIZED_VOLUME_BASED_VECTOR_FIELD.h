//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD
//#####################################################################
#ifndef __OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD__
#define __OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
namespace PhysBAM{

template<class T> class TETRAHEDRALIZED_VOLUME;
    
template<class T_input>
class OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD:public OPENGL_VECTOR_FIELD_3D<T_input>
{
    typedef T_input T;
public:
    using OPENGL_VECTOR_FIELD_3D<T>::vector_field;using OPENGL_VECTOR_FIELD_3D<T>::vector_locations;
    using OPENGL_VECTOR_FIELD_3D<T>::size;

    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume;
    ARRAY<VECTOR<T,3> >& V;

    OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,ARRAY<VECTOR<T,3> >& V);
    ~OPENGL_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD();

    void Update();  // Call when tetrahedralized volume/V change

    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};
}
#endif

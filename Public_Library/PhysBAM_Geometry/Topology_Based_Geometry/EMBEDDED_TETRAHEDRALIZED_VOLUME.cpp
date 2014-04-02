//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_TETRAHEDRALIZED_VOLUME
//#####################################################################
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
namespace PhysBAM{
//#####################################################################
// Register this class as read-write
//#####################################################################
namespace {
bool Register_Embedded_Tetrahedralized_Volume(){
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<EMBEDDED_TETRAHEDRALIZED_VOLUME<float> >();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<EMBEDDED_TETRAHEDRALIZED_VOLUME<double> >();
#endif
    return true;
}
static bool registered=Register_Embedded_Tetrahedralized_Volume();
}
//#####################################################################
// Register this class as read-write
//#####################################################################
template<class T> EMBEDDED_TETRAHEDRALIZED_VOLUME<T>::
EMBEDDED_TETRAHEDRALIZED_VOLUME(TETRAHEDRALIZED_VOLUME<T>& simplicial_object_input)
    :EMBEDDED_OBJECT<TV,3>(simplicial_object_input)
{
    PHYSBAM_ASSERT(registered);
}
//#####################################################################
template class EMBEDDED_TETRAHEDRALIZED_VOLUME<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EMBEDDED_TETRAHEDRALIZED_VOLUME<double>;
#endif
}

//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_TRIANGULATED_OBJECT
//#####################################################################
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TRIANGULATED_OBJECT.h>
namespace PhysBAM{
//#####################################################################
// Register this class as read-write
//#####################################################################
namespace {
bool Register_Embedded_Triangulated_Object(){
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<EMBEDDED_TRIANGULATED_OBJECT<VECTOR<float,2> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<EMBEDDED_TRIANGULATED_OBJECT<VECTOR<float,3> > >();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<EMBEDDED_TRIANGULATED_OBJECT<VECTOR<double,2> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<EMBEDDED_TRIANGULATED_OBJECT<VECTOR<double,3> > >();
#endif
    return true;
}
static bool registered=Register_Embedded_Triangulated_Object();
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> EMBEDDED_TRIANGULATED_OBJECT<TV>::
EMBEDDED_TRIANGULATED_OBJECT(T_TRIANGULATED_OBJECT& simplicial_object_input)
    :EMBEDDED_OBJECT<TV,2>(simplicial_object_input)
{
    PHYSBAM_ASSERT(registered);
}
//#####################################################################
template class EMBEDDED_TRIANGULATED_OBJECT<VECTOR<float,2> >;
template class EMBEDDED_TRIANGULATED_OBJECT<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EMBEDDED_TRIANGULATED_OBJECT<VECTOR<double,2> >;
template class EMBEDDED_TRIANGULATED_OBJECT<VECTOR<double,3> >;
#endif
}

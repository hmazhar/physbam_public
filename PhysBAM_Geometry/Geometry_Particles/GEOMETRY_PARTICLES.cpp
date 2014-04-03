//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GEOMETRY_PARTICLES<TV>::
GEOMETRY_PARTICLES(ARRAY_COLLECTION* array_collection_input)
    :V(0,0),store_velocity(false)
{
    delete array_collection;array_collection=array_collection_input;
    Initialize_Array_Collection();
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GEOMETRY_PARTICLES<TV>::
GEOMETRY_PARTICLES()
    :V(0,0),store_velocity(false)
{}
//#####################################################################
// Destructor 
//#####################################################################
template<class TV> GEOMETRY_PARTICLES<TV>::
~GEOMETRY_PARTICLES()
{}
//#####################################################################
// Function Initialize_Array_Collection
//#####################################################################
template<class TV> void GEOMETRY_PARTICLES<TV>::
Initialize_Array_Collection()
{
    POINT_CLOUD<TV>::Initialize_Array_Collection();
}
//#####################################################################
template class GEOMETRY_PARTICLES<VECTOR<float,1> >;
template class GEOMETRY_PARTICLES<VECTOR<float,2> >;
template class GEOMETRY_PARTICLES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GEOMETRY_PARTICLES<VECTOR<double,1> >;
template class GEOMETRY_PARTICLES<VECTOR<double,2> >;
template class GEOMETRY_PARTICLES<VECTOR<double,3> >;
#endif 
}

//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Geoffrey Irving, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> STRUCTURE_REGISTRY<TV>::
STRUCTURE_REGISTRY()
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> STRUCTURE_REGISTRY<TV>::
~STRUCTURE_REGISTRY()
{}
template class STRUCTURE_REGISTRY<VECTOR<float,1> >;
template class STRUCTURE_REGISTRY<VECTOR<float,2> >;
template class STRUCTURE_REGISTRY<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class STRUCTURE_REGISTRY<VECTOR<double,1> >;
template class STRUCTURE_REGISTRY<VECTOR<double,2> >;
template class STRUCTURE_REGISTRY<VECTOR<double,3> >;
#endif

//#####################################################################
// Copyright 2008, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDED_HORIZONTAL_PLANE
//#####################################################################
#include <PhysBAM_Geometry/Basic_Geometry/BOUNDED_HORIZONTAL_PLANE.h>
using namespace PhysBAM;
//#####################################################################
template class BOUNDED_HORIZONTAL_PLANE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BOUNDED_HORIZONTAL_PLANE<VECTOR<double,3> >;
#endif

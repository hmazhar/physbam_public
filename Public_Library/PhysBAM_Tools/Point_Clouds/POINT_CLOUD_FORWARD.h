//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header POINT_CLOUD_FORWARD
//#####################################################################
#ifndef __POINT_CLOUD_FORWARD__
#define __POINT_CLOUD_FORWARD__
#include <PhysBAM_Tools/Arrays/ATTRIBUTE_ID.h>

namespace PhysBAM{

template<class TV> class POINT_CLOUD;

template<class T_POINT_CLOUD> class POINT_CLOUD_POOL;

template<class TV> class POINT_CLOUD_CONNECTIVITY;
template<class TV,class T_POINT_CLOUD> class POINT_CLOUD_SUBSET;
const ATTRIBUTE_ID ATTRIBUTE_ID_X(1);
const ATTRIBUTE_ID ATTRIBUTE_ID_ID(20);
}
#endif

//#####################################################################
// Copyright 2006-2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header SPATIAL_ACCELERATION_FORWARD
//#####################################################################
#ifndef __SPATIAL_ACCELERATION_FORWARD__
#define __SPATIAL_ACCELERATION_FORWARD__

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>

namespace PhysBAM{
template<class TV> class BOX_HIERARCHY;
template<class TV> class SEGMENT_HIERARCHY;
template<class T> class TRIANGLE_HIERARCHY;
template<class T> class TRIANGLE_HIERARCHY_2D;
template<class T> class TETRAHEDRON_HIERARCHY;
template<class T> class PARTICLE_PARTITION;
template<class TV> class GEOMETRY_PARTICLES;
template<class TV,class T_ARRAY=ARRAY_VIEW<TV> > class PARTICLE_HIERARCHY;
}
#endif

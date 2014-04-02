//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Jon Gretarsson, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RAY
//#####################################################################
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
namespace PhysBAM{
//#####################################################################
// Function Compute_Lazy_Box_Intersection_Acceleration_Data
//#####################################################################
template<class TV> void RAY<TV>::
Compute_Lazy_Box_Intersection_Acceleration_Data()
{
    STATIC_ASSERT(TV::dimension==3);
    if(computed_lazy_box_intersection_acceleration_data) return;
    const T tolerance=(T)1e-10,large_number=(T)1e10;
    if(direction.x<0){direction_is_negative.x=1;inverse_direction.x=(direction.x<-tolerance)?1/direction.x:-large_number;}
    else{direction_is_negative.x=0;inverse_direction.x=(direction.x>tolerance)?1/direction.x:large_number;}
    if(direction.y<0){direction_is_negative.y=1;inverse_direction.y=(direction.y<-tolerance)?1/direction.y:-large_number;}
    else{direction_is_negative.y=0;inverse_direction.y=(direction.y>tolerance)?1/direction.y:large_number;}
    if(direction.z<0){direction_is_negative.z=1;inverse_direction.z=(direction.z<-tolerance)?1/direction.z:-large_number;}
    else{direction_is_negative.z=0;inverse_direction.z=(direction.z>tolerance)?1/direction.z:large_number;}
    computed_lazy_box_intersection_acceleration_data=true;
}
//#####################################################################
template void RAY<VECTOR<float,3> >::Compute_Lazy_Box_Intersection_Acceleration_Data();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void RAY<VECTOR<double,3> >::Compute_Lazy_Box_Intersection_Acceleration_Data();
#endif
}

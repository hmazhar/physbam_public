//####################################################################
// Copyright 2004, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __IRRADIANCE_SAMPLE__
#define __IRRADIANCE_SAMPLE__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T>
class IRRADIANCE_SAMPLE
{
public:
    VECTOR<T,3> location;
    VECTOR<T,3> direction;
    VECTOR<T,3> irradiance;
    T harmonic_mean_distance;
    
    IRRADIANCE_SAMPLE()
    {}

    IRRADIANCE_SAMPLE(const VECTOR<T,3>& location_input,const VECTOR<T,3>& direction_input,const VECTOR<T,3>& irradiance_input,const T harmonic_mean_distance_input)
        :location(location_input),direction(direction_input),irradiance(irradiance_input),harmonic_mean_distance(harmonic_mean_distance_input)
    {}
    
//#####################################################################
};
}
#endif

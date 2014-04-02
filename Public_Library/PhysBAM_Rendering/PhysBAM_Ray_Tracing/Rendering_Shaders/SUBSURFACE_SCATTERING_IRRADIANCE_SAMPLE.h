//#####################################################################
// Copyright 2004-2006, Jiayi Chong, Igor Neverov, Andrew Selle, Mike Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE
//#####################################################################
#ifndef __SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE__
#define __SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T>
class SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE
{
public:
    VECTOR<T,3> position;
    VECTOR<T,3> normal;
    VECTOR<T,3> transmitted_irradiance;
    VECTOR<T,3> transmitted_irradiance_product;
    T area;

    SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE()
        :area(0)
    {}

    SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE(const VECTOR<T,3>& position_input,const VECTOR<T,3>& normal_input,const T area_input)
        :position(position_input),normal(normal_input),area(area_input)
    {}

//#####################################################################
};
}
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Read_Write/Rendering_Shaders/READ_WRITE_SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE.h>
#endif

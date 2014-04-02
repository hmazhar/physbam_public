//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE
//#####################################################################
#ifndef __READ_WRITE_SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE__
#define __READ_WRITE_SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE<T>,RW>
{
public:
    static void Read(std::istream& input,SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE<T>& object)
    {Read_Binary<RW>(input,object.position,object.normal,object.transmitted_irradiance,object.transmitted_irradiance_product,object.area);}

    static void Write(std::ostream& output,const SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE<T>& object)
    {Write_Binary<RW>(output,object.position,object.normal,object.transmitted_irradiance,object.transmitted_irradiance_product,object.area);}
};
}
#endif

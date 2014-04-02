#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE
//#####################################################################
#ifndef __READ_WRITE_SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE__
#define __READ_WRITE_SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE__

#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Math_Tools/READ_WRITE_RANGE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<T>,RW>
{
public:
    static void Read(std::istream& input,SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<T>& object)
    {Read_Binary<RW>(input,object.samples,object.bounding_box);}

    static void Write(std::ostream& output,const SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<T>& object)
    {Write_Binary<RW>(output,object.samples,object.bounding_box);}
};
}
#endif
#endif

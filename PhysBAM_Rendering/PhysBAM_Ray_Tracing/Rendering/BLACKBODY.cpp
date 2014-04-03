//#####################################################################
// Copyright 2002, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//##################################################################### 
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/BLACKBODY.h>
#include <cmath>
using namespace PhysBAM;
//#####################################################################
// Function Calculate_Blackbody_Spectrum
//#####################################################################
// radiance_spectrum is in units of Watts/(steradian*m^2)/m 
// input grid of wavelengths in nanometers
template<class T> void BLACKBODY<T>::
Calculate_Radiance_Spectrum(const T temperature,const GRID<VECTOR<T,1> >& grid,ARRAY<T,VECTOR<int,1> >& radiance_spectrum) const
{
    T constant_1=T(plancks_constant*sqr(speed_of_light)/2),constant_2=T(plancks_constant*speed_of_light/boltzmanns_constant);
    for(int i=1;i<=grid.counts.x;i++){
        T lambda=grid.Axis_X(i,1);
        radiance_spectrum(i)=constant_1/(cube(lambda)*sqr(lambda)*(exp(constant_2/(lambda*temperature))-1));}
}
//#####################################################################
// Function Calculate_XYZ
//#####################################################################
// assumes f is defined on the grid
template<class T> VECTOR<T,3> BLACKBODY<T>::
Calculate_XYZ(const T temperature) const
{
    ARRAY<T,VECTOR<int,1> > radiance_spectrum(1,cie.grid.counts.x);
    Calculate_Radiance_Spectrum(temperature,cie.grid,radiance_spectrum);
    return cie.Calculate_XYZ(radiance_spectrum);
}
//#####################################################################
template class BLACKBODY<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BLACKBODY<double>;
#endif

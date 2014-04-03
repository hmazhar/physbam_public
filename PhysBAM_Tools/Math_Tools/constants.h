//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Avi Robinson-Mosher, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __constants__
#define __constants__

#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <cmath>
namespace PhysBAM{

using ::std::sqrt;

const double pi=4*atan(1.);
const double pi_squared=pi*pi;
const double two_pi=2*pi;
const double half_pi=.5*pi;
const double one_fourth_pi=.25*pi;
const double one_sixth_pi=1./6*pi;
const double two_thirds_pi=2./3*pi;
const double four_thirds_pi=4./3*pi;
const double one_over_pi=1./pi;
const double one_over_two_pi=.5/pi;
const double one_over_four_pi=.25/pi;

const double one_third=1./3;
const double two_thirds=2./3;
const double four_thirds=4./3;
const double five_thirds=5./3;
const double one_sixth=1./6;
const double one_ninth=1./9;
const double one_twelfth=1./12;
const double one_twenty_fourth=1./24;
const double one_twenty_seventh=1./27;
const double one_sixtieth=1./60;
const double thirteen_over_twelve=13./12;
const double root_two=sqrt(2.);
const double root_three=sqrt(3.);
const double root_six=sqrt(6.);
const double root_two_thirds=sqrt(2./3);
const double one_over_root_two=1./sqrt(2.);
const double one_over_root_three=1./sqrt(3.);

const double speed_of_light=2.99792458e8; // m/s
const double plancks_constant=6.6260755e-34; // J*s
const double boltzmanns_constant=1.380658e-23; // J/K

template<int d> struct unit_sphere_size{STATIC_ASSERT(d<4);static const double value;};
template<int d> const double unit_sphere_size<d>::value=d==0?(double)0:d==1?(double)2:d==2?pi:four_thirds_pi;

}
#endif

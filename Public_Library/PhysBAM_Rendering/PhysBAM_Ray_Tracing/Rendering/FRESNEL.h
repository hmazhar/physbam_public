//#####################################################################
// Copyright 2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Reference: PBRT Book, Matt Pharr and Greg Humphrey
//#####################################################################
#ifndef __FRESNEL__
#define __FRESNEL__

namespace PhysBAM{

template<class T>
class FRESNEL
{
public:
    static T Dielectric_Reflection(T eta_incident,T eta_transmitted,T cos_incident)
    {cos_incident=clamp(cos_incident,(T)-1,(T)1);
    bool entering=cos_incident>0;
    if(!entering){T temp=eta_transmitted;eta_transmitted=eta_incident;eta_incident=temp;}
    T sin_transmitted=eta_incident/eta_transmitted*sqrt(max((T)0,(T)1-cos_incident*cos_incident));
    if(sin_transmitted>1) return 1; // total internal reflection
    T cos_transmitted=sqrt(max((T)0,1-sqr(sin_transmitted)));
    T reflection_parallel=((eta_transmitted*cos_incident)-(eta_incident*cos_transmitted))/((eta_transmitted*cos_incident)+(eta_incident*cos_transmitted));
    T reflection_perpendicular=((eta_incident*cos_incident)-(eta_transmitted*cos_transmitted))/((eta_incident*cos_incident)+(eta_transmitted*cos_transmitted));
    return (T).5*(sqr(reflection_parallel)+sqr(reflection_perpendicular));}

    static T Conductor_Reflection(T index_of_refraction,T cos_incident,T absorption_coefficient)
    {T sum=sqr(index_of_refraction)+sqr(absorption_coefficient),cos_incident_squared=sqr(cos_incident),product=2*index_of_refraction*cos_incident;
    T reflection_parallel_sqauared=(sum*cos_incident_squared-product+1)/(sum*cos_incident_squared+product+1);
    T reflection_perp_squared=(sum-product+cos_incident_squared)/(sum+product+cos_incident_squared);
    return (T).5*(reflection_perp_squared+reflection_parallel_sqauared);}

    static T Conductor_Eta(const T reflection_in_normal_direction) // assumes absorption is zero
    {return (1+sqrt(reflection_in_normal_direction))/(1-sqrt(reflection_in_normal_direction));}

    static T Conductor_Absorption(const T reflectance_in_normal_direction) // assumes eta is one
    {T clamped_reflectance=clamp(reflectance_in_normal_direction,0,(T).9999);
    return (T)2*clamped_reflectance/sqrt(1-clamped_reflectance);}

//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Neil Molino, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Geometry/Basic_Geometry/ELLIPSOID.h>
using namespace PhysBAM;
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> ELLIPSOID<T>::
Normal(const TV& location) const   
{
    return orientation.Rotate(sqr(radii).Solve_Linear_System(orientation.Inverse_Rotate(location-center))).Normalized();
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool ELLIPSOID<T>::
Inside(const TV& location,const T thickness_over_two) const
{
    TV scaled_offset=radii.Solve_Linear_System(orientation.Inverse_Rotate(location-center));
    return scaled_offset.Magnitude_Squared() <= sqr(1-thickness_over_two);
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool ELLIPSOID<T>::
Outside(const TV& location,const T thickness_over_two) const
{
    TV scaled_offset=radii.Solve_Linear_System(orientation.Inverse_Rotate(location-center));
    return scaled_offset.Magnitude_Squared() >= sqr(1+thickness_over_two);
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class T> bool ELLIPSOID<T>::
Boundary(const TV& location,const T thickness_over_two) const
{
    return !Inside(location,thickness_over_two) && !Outside(location,thickness_over_two);
}
//#####################################################################
// Function Approximate_Surface
//#####################################################################
// not necessarily the closest point on the surface, but approximates it in some sense
template<class T> VECTOR<T,3> ELLIPSOID<T>::
Approximate_Surface(const TV& location) const 
{
    TV scaled_offset=radii.Solve_Linear_System(orientation.Inverse_Rotate(location-center));
    return orientation.Rotate(radii*scaled_offset.Normalized())+center;
}
//#####################################################################
// Function Approximate_Surface
//#####################################################################
template<class T> T ELLIPSOID<T>::
Approximate_Signed_Distance(const TV& location) const       
{
    TV offset=orientation.Inverse_Rotate(location-center),scaled_offset=radii.Solve_Linear_System(offset);
    T scaled_magnitude=scaled_offset.Normalize();
    TV closest_offset=radii*scaled_offset;
    T distance=(closest_offset-offset).Magnitude();
    return scaled_magnitude<1?-distance:distance;
}
//#####################################################################
// Function Covariance_Ellipsoid
//#####################################################################
template<class T> template<class T_ARRAY_TV> ELLIPSOID<T> ELLIPSOID<T>::
Covariance_Ellipsoid(const T_ARRAY_TV& points)
{
    TV average=ARRAYS_COMPUTATIONS::Average(points);
    SYMMETRIC_MATRIX<T,3> covariance;DIAGONAL_MATRIX<T,3> eigenvalues;MATRIX<T,3> eigenvectors;
    for(int p=1;p<=points.m;p++){
        TV variance=points(p)-average;
        covariance.x11+=variance.x*variance.x;covariance.x21+=variance.y*variance.x;covariance.x31+=variance.z*variance.x;
        covariance.x22+=variance.y*variance.y;covariance.x32+=variance.z*variance.y;covariance.x33+=variance.z*variance.z;}
    covariance/=points.m-(T)1;covariance.Fast_Solve_Eigenproblem(eigenvalues,eigenvectors);
    return ELLIPSOID<T>(average,eigenvalues.Sqrt(),eigenvectors);
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template class ELLIPSOID<T>; \
    template ELLIPSOID<T> ELLIPSOID<T>::Covariance_Ellipsoid(const ARRAY<VECTOR<T,3> >& points);
INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif

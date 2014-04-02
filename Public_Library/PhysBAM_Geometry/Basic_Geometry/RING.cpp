//#####################################################################
// Copyright 2003-2008, Eran Guendelman, Geoffrey Irving, Sergey Koltakov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RING
//##################################################################### 
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/RING.h>
namespace PhysBAM{
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > RING<T>::
Bounding_Box() const
{
    return CYLINDER<T>(plane1.x1,plane2.x1,outer_radius).Bounding_Box();
}
//#####################################################################
// Function Signed_Distance
//#####################################################################
template<class T> T RING<T>::
Signed_Distance(const TV& X) const
{
    TV v=X-plane1.x1;
    T plane1_phi=TV::Dot_Product(v,plane1.normal),plane_phi=max(plane1_phi,-height-plane1_phi);
    T radius=(v-plane1_phi*plane1.normal).Magnitude();
    T radial_phi=max(radius-outer_radius,inner_radius-radius);
    return radial_phi>0 && plane_phi>0?Vector(radial_phi,plane_phi).Magnitude():max(radial_phi,plane_phi);
}
//#####################################################################
// Function Suface
//#####################################################################
template<class T> VECTOR<T,3> RING<T>::
Surface(const TV& X) const 
{
    TV v=X-plane1.x1;
    T plane1_phi=TV::Dot_Product(v,plane1.normal),plane2_phi=-height-plane1_phi,plane_phi=max(plane1_phi,plane2_phi);
    TV radial_direction=v-plane1_phi*plane1.normal;
    T radius=radial_direction.Normalize();
    T outer_phi=radius-outer_radius,inner_phi=inner_radius-radius;
    T radial_phi=max(inner_phi,outer_phi);
    if(radial_phi>0 || radial_phi>plane_phi){ // closest point is on one of the cylinders
        T nearest_radius=outer_phi>inner_phi?outer_radius:inner_radius;
        return plane1.x1+clamp(plane1_phi,-height,(T)0)*plane1.normal+nearest_radius*radial_direction;}
    return plane1_phi > plane2_phi?X-plane1_phi*plane1.normal:X+plane2_phi*plane1.normal;
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> RING<T>::
Normal(const TV& X) const 
{
    // compute plane normal
    TV v=X-plane1.x1;
    T plane1_phi=TV::Dot_Product(v,plane1.normal),plane2_phi=-height-plane1_phi;
    TV plane_normal=plane1_phi>plane2_phi?plane1.normal:plane2.normal;
    T plane_phi=max(plane1_phi,plane2_phi);
    // compute radial normal
    TV radial_normal=v-plane1_phi*plane1.normal;
    T radius=radial_normal.Normalize();
    T outer_phi=radius-outer_radius,inner_phi=inner_radius-radius;
    if(inner_phi>outer_phi) radial_normal=-radial_normal;
    T radial_phi=max(inner_phi,outer_phi);
    // combine normals
    if(radial_phi>0 && plane_phi>0){
        VECTOR<T,2> weights=Vector(radial_phi,plane_phi).Normalized();
        return weights.x*radial_normal+weights.y*plane_normal;}
    return radial_phi>plane_phi?radial_normal:plane_normal;
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> RING<T>::
Normal(const TV& X,const int aggregate) const 
{
    assert(aggregate >= 1 && aggregate <= 4);
    if(aggregate==1){ // outer cylinder
        TV normal=X-plane1.x1;
        normal-=TV::Dot_Product(normal,plane1.normal)*plane1.normal;
        return normal.Normalized();}
    else if(aggregate==2){ // inner cylinder
        TV normal=X-plane1.x1;
        normal-=TV::Dot_Product(normal,plane1.normal)*plane1.normal;
        return -normal.Normalized();}
    else if(aggregate==3)
        return plane1.Normal();
    else 
        return plane2.Normal();
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class T> VECTOR<T,2> RING<T>::
Principal_Curvatures(const TV& X) const
{
    TV v=X-plane1.x1;
    T plane1_phi=TV::Dot_Product(v,plane1.normal),plane_phi=max(plane1_phi,-height-plane1_phi);
    T radius=(v-plane1_phi*plane1.normal).Magnitude();
    T outer_phi=radius-outer_radius,inner_phi=inner_radius-radius;
    if(plane_phi>max(inner_phi,outer_phi)) return VECTOR<T,2>();
    return VECTOR<T,2>(1/(outer_phi>inner_phi?outer_radius:-inner_radius),0);
}
//#####################################################################
// Function Lazy_Inside
//#####################################################################
template<class T> bool RING<T>::
Lazy_Inside(const TV& X) const 
{
    if(!plane1.Lazy_Inside(X) || !plane2.Lazy_Inside(X)) return false;
    TV v=X-plane1.x1;
    v-=TV::Dot_Product(v,plane1.normal)*plane1.normal;
    T radius_sqr=v.Magnitude_Squared();
    return radius_sqr<=sqr(outer_radius) && radius_sqr>=sqr(inner_radius);
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class T> bool RING<T>::
Lazy_Outside(const TV& X) const  
{
    return !Lazy_Inside(X);
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool RING<T>::
Inside(const TV& X,const T thickness_over_two) const 
{
    if(plane1.Inside(X,thickness_over_two) && plane2.Inside(X,thickness_over_two)){ // inside both planes, check the cylinder
        TV distance=X-plane1.x1;
        distance-=TV::Dot_Product(distance,plane1.normal)*plane1.normal;
        T distance_squared = TV::Dot_Product(distance,distance);
        if(distance_squared <= sqr(outer_radius-thickness_over_two) && distance_squared >= sqr(inner_radius+thickness_over_two)) return true;}
    return false;
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool RING<T>::
Outside(const TV& X,const T thickness_over_two) const  
{
    if(plane1.Outside(X,thickness_over_two) || plane2.Outside(X,thickness_over_two)) return true;
    else{ // check the cylinder
        TV distance=X-plane1.x1;
        distance-=TV::Dot_Product(distance,plane1.normal)*plane1.normal;
        T distance_squared = TV::Dot_Product(distance,distance);
        if(distance_squared >= sqr(outer_radius+thickness_over_two) || distance_squared <= sqr(inner_radius-thickness_over_two)) return true;
        else return false;}
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class T> bool RING<T>::
Boundary(const TV& X,const T thickness_over_two) const
{
    return !Inside(X,thickness_over_two) && !Outside(X,thickness_over_two);
}
//#####################################################################
// Function Name
//#####################################################################
template<class T> std::string RING<T>::
Name()
{
    return "RING<T>";
}
//#####################################################################
template class RING<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RING<double>;
#endif
}

//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CYLINDER
//##################################################################### 
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
using namespace PhysBAM;
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> CYLINDER<T>::
Normal(const TV& location,const int aggregate) const 
{
    assert(aggregate >= 1 && aggregate <= 3);
    if(aggregate == 1) return (location-plane1.x1).Projected_Orthogonal_To_Unit_Direction(plane1.normal).Normalized(); // cylinder
    else if(aggregate == 2) return plane1.Normal();
    else return plane2.Normal();
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool CYLINDER<T>::
Inside(const TV& location,const T thickness_over_two) const 
{
    if(!plane1.Inside(location,thickness_over_two) || !plane2.Inside(location,thickness_over_two)) return false;
    return (location-plane1.x1).Projected_Orthogonal_To_Unit_Direction(plane1.normal).Magnitude_Squared() <= sqr(radius-thickness_over_two);
}
//#####################################################################
// Function Lazy_Inside
//#####################################################################
template<class T> bool CYLINDER<T>::
Lazy_Inside(const TV& location) const 
{
    if(!plane1.Lazy_Inside(location) || !plane2.Lazy_Inside(location)) return false;
    return (location-plane1.x1).Projected_Orthogonal_To_Unit_Direction(plane1.normal).Magnitude_Squared() <= sqr(radius);
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool CYLINDER<T>::
Outside(const TV& location,const T thickness_over_two) const  
{
    if(plane1.Outside(location,thickness_over_two) || plane2.Outside(location,thickness_over_two)) return true;
    return (location-plane1.x1).Projected_Orthogonal_To_Unit_Direction(plane1.normal).Magnitude_Squared() >= sqr(radius+thickness_over_two);
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class T> bool CYLINDER<T>::
Lazy_Outside(const TV& location) const  
{
    if(plane1.Lazy_Outside(location) || plane2.Lazy_Outside(location)) return true;
    return (location-plane1.x1).Projected_Orthogonal_To_Unit_Direction(plane1.normal).Magnitude_Squared() >= sqr(radius);
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class T> bool CYLINDER<T>::
Boundary(const TV& location,const T thickness_over_two) const
{
    return !Inside(location,thickness_over_two) && !Outside(location,thickness_over_two);
}
//#####################################################################
// Function Surface
//#####################################################################
template<class T> VECTOR<T,3> CYLINDER<T>::
Surface(const TV& location) const
{
    TV v=location-plane1.x1;T axial_distance=-TV::Dot_Product(v,plane1.normal);
    TV radial_direction=v+axial_distance*plane1.normal;T radial_distance=radial_direction.Normalize();
    T radial_distance_minus_radius=radial_distance-radius;
    if(radial_distance_minus_radius>0 || (radial_distance_minus_radius>-axial_distance && radial_distance_minus_radius>axial_distance-height)) // closest point is on infinite cylinder
        return plane1.x1-clamp(axial_distance,(T)0,height)*plane1.normal+radius*radial_direction;
    if(axial_distance < height-axial_distance) return location+axial_distance*plane1.normal; // closest point is on plane1
    else return location-(height-axial_distance)*plane1.normal; // closest point is on plane2
}
//#####################################################################
// Function Signed_Distance
//#####################################################################
template<class T> T CYLINDER<T>::
Signed_Distance(const TV& location) const
{  
    TV v=location-plane1.x1;
    T plane1_distance=TV::Dot_Product(v,plane1.normal),plane_distance=max(plane1_distance,-height-plane1_distance);
    T cylinder_distance=(v-plane1_distance*plane1.normal).Magnitude()-radius;
    return cylinder_distance>0 && plane_distance>0?sqrt(sqr(cylinder_distance)+sqr(plane_distance)):max(cylinder_distance,plane_distance);
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> CYLINDER<T>::
Normal(const TV& location) const
{
    TV v=location-plane1.x1;
    T plane1_distance=TV::Dot_Product(v,plane1.normal),plane2_distance=-height-plane1_distance;
    TV infinite_cylinder_normal=v-plane1_distance*plane1.normal;
    T cylinder_distance=infinite_cylinder_normal.Normalize()-radius;
    if(plane1_distance>=plane2_distance){
        if(cylinder_distance>0 && plane1_distance>0){
            T magnitude=sqrt(sqr(cylinder_distance)+sqr(plane1_distance));
            return cylinder_distance/magnitude*infinite_cylinder_normal+plane1_distance/magnitude*plane1.normal;}
        else if(cylinder_distance>plane1_distance) return infinite_cylinder_normal;
        return plane1.normal;}
    if(cylinder_distance>0 && plane2_distance>0){
        T magnitude=sqrt(sqr(cylinder_distance)+sqr(plane2_distance));
        return cylinder_distance/magnitude*infinite_cylinder_normal+plane2_distance/magnitude*plane2.normal;}
    else if(cylinder_distance>plane2_distance) return infinite_cylinder_normal;
    return plane2.normal;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > CYLINDER<T>::
Bounding_Box() const
{
    TV edge1=(T)2*radius*plane1.normal.Unit_Orthogonal_Vector();
    TV edge2=TV::Cross_Product(plane1.normal,edge1);
    TV corner=plane1.x1-(T).5*(edge1+edge2);
    ORIENTED_BOX<TV> oriented_box(corner,-height*plane1.normal,edge1,edge2); // TODO(jontg): ...
    return oriented_box.Axis_Aligned_Bounding_Box();
}
//#####################################################################
template class CYLINDER<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CYLINDER<double>;
#endif

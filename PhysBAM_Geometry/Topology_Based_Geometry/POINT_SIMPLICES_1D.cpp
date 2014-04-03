//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
namespace PhysBAM{
//#####################################################################
// Register this class as read-write
//#####################################################################
namespace{
bool Register_Point_Simplices_1d(){
    STRUCTURE_REGISTRY<VECTOR<float,1> >::Register<POINT_SIMPLICES_1D<float> >();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    STRUCTURE_REGISTRY<VECTOR<double,1> >::Register<POINT_SIMPLICES_1D<double> >();
#endif
    return true;
}
static bool registered=Register_Point_Simplices_1d();
}
//#####################################################################
// Function Closest_Point_On_Boundary
//#####################################################################
template<class T> POINT_SIMPLICES_1D<T>::
POINT_SIMPLICES_1D(POINT_SIMPLEX_MESH& point_simplex_mesh_input,GEOMETRY_PARTICLES<VECTOR<T,1> >& particles_input)
    :MESH_OBJECT<TV,POINT_SIMPLEX_MESH>(point_simplex_mesh_input,particles_input),point_simplex_list(0),particle_partition(0),hierarchy(0),number_point_simplices(0)
{
    PHYSBAM_ASSERT(registered);
}
//#####################################################################
// Function Closest_Point_On_Boundary
//#####################################################################
template<class T> VECTOR<T,1> POINT_SIMPLICES_1D<T>::
Closest_Point_On_Boundary(const TV& location,const T max_depth,const T thickness_over_two,int* closest_point_simplex,T* distance) const 
{
    T min_distance_squared=FLT_MAX;TV closest_point;
    for(int i=1;i<=mesh.elements.m;i++){
        POINT_SIMPLEX_1D<T> point_simplex=point_simplex_list?(*point_simplex_list)(i):Get_Element(i);
        TV new_point=point_simplex.x1;
        T distance_squared=(new_point-location).Magnitude_Squared();
        if(distance_squared < min_distance_squared){
            min_distance_squared=distance_squared;closest_point=new_point;
            if(closest_point_simplex) *closest_point_simplex=i;}}
    if(distance) *distance=sqrt(min_distance_squared);
    return closest_point;
}
//#####################################################################
// Function Calculate_Signed_Distance
//#####################################################################
template<class T> T POINT_SIMPLICES_1D<T>:: 
Calculate_Signed_Distance(const TV& location,T thickness_over_two) const
{
    T distance;
    Closest_Point_On_Boundary(location,0,thickness_over_two,0,&distance);
    if(Inside(location,thickness_over_two)) distance*=-1;
    return distance;
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool POINT_SIMPLICES_1D<T>::
Inside(const TV& location,T thickness_over_two) const
{
    for(int i=1;i<=mesh.elements.m;i++){
        POINT_SIMPLEX_1D<T> point_simplex=point_simplex_list?(*point_simplex_list)(i):Get_Element(i);
        T direction=mesh.directions(i)?1:-1;
        T robust_point_simplex_location=point_simplex.x1.x+direction*thickness_over_two;
        if(direction*(location.x-robust_point_simplex_location)>0) return false;}
    return true;
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool POINT_SIMPLICES_1D<T>::
Outside(const TV& location,T thickness_over_two) const
{
    return !Inside(location,thickness_over_two);
}
//#####################################################################
// Function Inside_Any_Simplex
//#####################################################################
template<class T> bool POINT_SIMPLICES_1D<T>::
Inside_Any_Simplex(const TV& point,int& point_simplex_id,const T thickness_over_two) const
{
    for(int i=1;i<=mesh.elements.m;++i){
        POINT_SIMPLEX_1D<T> point_simplex=point_simplex_list?(*point_simplex_list)(i):Get_Element(i);
        if(point_simplex.x1.x-thickness_over_two <= point.x && point_simplex.x1.x+thickness_over_two >= point.x) return true;}
    return false;
}
//#####################################################################
template class POINT_SIMPLICES_1D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class POINT_SIMPLICES_1D<double>;
#endif 
}

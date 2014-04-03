//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_PLANE_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(const RANGE<VECTOR<T,3> >& box,const PLANE<T>& plane,const T thickness_over_two)
{
    typedef VECTOR<T,3> TV;
    bool points_on_positive_side=false,points_on_negative_side=false;
    for(int i=0;i<=1;i++) for(int j=0;j<=1;j++) for(int k=0;k<=1;k++){
        TV test_point(i?box.min_corner.x:box.max_corner.x,j?box.min_corner.y:box.max_corner.y,k?box.min_corner.z:box.max_corner.z);
        T distance=TV::Dot_Product(plane.normal,test_point-plane.x1);
        if(distance>-thickness_over_two)points_on_positive_side=true;
        if(distance<thickness_over_two)points_on_negative_side=true;}
    return points_on_positive_side && points_on_negative_side;
}
//#####################################################################
// Function Halfspace_Intersection_Size
//#####################################################################
// The case enumeration relies on z00>=z10>=z01>=z11, which follows from normal.z>=normal.y>=normal.x>=0.  Cases are distinguished by which of
// these four values are negative or greater than 1.  Note that normal.z>=normal.y implies z00-z01=z10-z11<=1.
// 1>0>z00>z10>z01>z11 - case 1 - no intersection
// 1>z00>0>z10>z01>z11 - case 3 - intersection contains (0,0,0)
// 1>z00>z10>0>z01>z11 - case 4 - intersection contains (0,0,0) and (1,0,0)
// 1>z00>z10>z01>0>z11 - case 5 with d=0 - intersection contains (0,0,0), (1,0,0), and (0,1,0)
// 1>z00>z10>z01>z11>0 - case 2 - intersection contains (0,0,0), (1,0,0), (0,1,0), and (1,1,0)
// z00>1>z10>z01>0>z11 - case 5 with d>0 - intersection contains (0,0,0), (1,0,0), (0,1,0), and (0,0,1)
// Impossible since normal.z>=normal.y: z00>1>z10>0>z01>z11, z00>z10>1>0>z01>z11, z00>z10>1>z01>0>z11, z00>z10>z01>1>0>z11, z00>1>0>z10>z01>z11
// Impossible since cube center is in intersection: z00>1>z10>z01>z11>0, z00>z10>1>z01>z11>0, z00>z10>z01>1>z11>0, z00>z10>z01>z11>1>0
// Let E(w) be the expression obtained by integrating w over the intersection volume.  Then, volume=E(1) and centroid=TV(E(x),E(y),E(z))/E(1).
template<class T> T Halfspace_Intersection_Size(const RANGE<VECTOR<T,3> >& box,const PLANE<T>& halfspace,VECTOR<T,3>* centroid)
{
    // Normalize input and remember how.  After normalization: box is a unit box, normal.z>=normal.y>=normal.x>=0.
    bool flip[3]={false},exchange_xz=false,exchange_yz=false,exchange_xy=false,complement=false;
    T z00=VECTOR<T,3>::Dot_Product(halfspace.x1-box.min_corner,halfspace.normal);VECTOR<T,3> scaling=box.Edge_Lengths(),normal=halfspace.normal*scaling;
    for(int i=1;i<=3;i++) if(normal(i)<0){normal(i)=-normal(i);z00+=normal(i);flip[i-1]=true;}
    if(normal.z<normal.x){exchange(normal.z,normal.x);exchange_xz=true;}
    if(normal.z<normal.y){exchange(normal.z,normal.y);exchange_yz=true;}
    if(normal.y<normal.x){exchange(normal.y,normal.x);exchange_xy=true;}
    // z00, z10, z01, z11 = intersections with lines x={0,1}, y={0,1}.  Normalize to cube center outside or on the boundary.
    z00/=normal.z;normal/=normal.z;T nx=normal.x,ny=normal.y;if(2*z00-nx-ny>1){z00=1-z00+nx+ny;complement=true;}
    T z01=z00-ny,z10=z00-nx,z11=z10-ny,volume=0;
    // At this point, there are six case.  Compute them analytically.  See above for details on the cases.
    if(z00<=0) {if(centroid) *centroid=box.Center();return complement?box.Size():0;}
    else if(z11>=0){volume=(T).5*(z00+z11);if(centroid) *centroid=VECTOR<T,3>(6*volume-nx,6*volume-ny,6*sqr(volume)+(T).5*(sqr(nx)+sqr(ny)))/(12*volume);}
    else if(z10<=0){volume=cube(z00)/(6*nx*ny);if(centroid) *centroid=(T).25*z00/normal;}
    else if(z01<=0){volume=(3*z00*z10+sqr(nx))/(6*ny);if(centroid){T z=(z00+z10)*(sqr(nx)+2*z00*z10);*centroid=VECTOR<T,3>(3*z10*(z00+z10)+z00*nx,z/ny,z)/(24*ny*volume);}}
    else{T a=cube(z00),b=cube(z10),c=cube(z01),d=max((T)0,z00-1),e=cube(d);volume=a-b-c-e;
        if(centroid) *centroid=((T).25*(a*z00-b*z10-c*z01-e*d)/normal-VECTOR<T,3>(b,c,e))/volume;volume/=(6*nx*ny);}
    // Undo the normalizations to convert normalized centroid and volume into world space.
    if(complement){volume=1-volume;if(centroid) *centroid=-((volume-1)**centroid+(T).5)/volume+1;}
    if(centroid){
        if(exchange_xy) exchange(centroid->y,centroid->x);if(exchange_yz) exchange(centroid->z,centroid->y);if(exchange_xz) exchange(centroid->z,centroid->x);
        for(int i=1;i<=3;i++) if(flip[i-1]) (*centroid)(i)=1-(*centroid)(i);
        *centroid=*centroid*scaling+box.min_corner;}
    return volume*box.Size();
}
//#####################################################################
template bool Intersects(const RANGE<VECTOR<float,3> >&,const PLANE<float>&,const float);
template float Halfspace_Intersection_Size(const RANGE<VECTOR<float,3> >&,const PLANE<float>&,VECTOR<float,3>*);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(const RANGE<VECTOR<double,3> >&,const PLANE<double>&,const double);
template double Halfspace_Intersection_Size(const RANGE<VECTOR<double,3> >&,const PLANE<double>&,VECTOR<double,3>*);
#endif
};
};

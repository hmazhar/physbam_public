//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_LINE_2D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(const RANGE<VECTOR<T,2> >& box,const LINE_2D<T>& line,const T thickness_over_two)
{
    bool points_on_positive_side=false,points_on_negative_side=false;
    for(int i=0;i<=1;i++)for(int j=0;j<=1;j++){
        VECTOR<T,2> test_point(i?box.min_corner.x:box.max_corner.x,j?box.min_corner.y:box.max_corner.y);
        T distance=VECTOR<T,2>::Dot_Product(line.normal,test_point-line.x1);
        if(distance>-thickness_over_two) points_on_positive_side=true;
        if(distance<thickness_over_two) points_on_negative_side=true;}
    return points_on_positive_side&&points_on_negative_side;
}
//#####################################################################
// Function Halfspace_Intersection_Size
//#####################################################################
template<class T> T Halfspace_Intersection_Size(const RANGE<VECTOR<T,2> >& box,const LINE_2D<T>& halfspace,VECTOR<T,2>* centroid)
{
    // Normalize input and remember how.  After normalization: box is a unit box, normal.z>=normal.y>=normal.x>=0.  Normalize to cube center outside or on the boundary.
    bool flip[2]={false},exchange_xy=false,compliment=false;
    T y0=VECTOR<T,2>::Dot_Product(halfspace.x1-box.min_corner,halfspace.normal);VECTOR<T,2> scaling=box.Edge_Lengths(),normal=halfspace.normal*scaling;
    for(int i=1;i<=2;i++) if(normal(i)<0){normal(i)=-normal(i);y0+=normal(i);flip[i-1]=true;}
    if(normal.y<normal.x){exchange(normal.y,normal.x);exchange_xy=true;}
    y0/=normal.y;normal/=normal.y;T nx=normal.x,volume=0;if(2*y0-nx>1){y0=1-y0+nx;compliment=true;}
    // At this point, there are three case.  Compute them analytically.  Same approach as 3D, but much simpler.
    if(y0<=0){if(centroid) *centroid=box.Center();return compliment?box.Size():0;}
    else if(y0-nx<=0){volume=sqr(y0)/(2*nx);if(centroid) *centroid=(T)one_third*y0/normal;}
    else{volume=y0-(T).5*nx;if(centroid) *centroid=VECTOR<T,2>(y0-(T)two_thirds*nx,sqr(volume)+sqr(nx)/12)/(volume*2);}
    // Undo the normalizations to convert normalized centroid and volume into world space.
    if(compliment){volume=1-volume;if(centroid) *centroid=-((volume-1)**centroid+(T).5)/volume+1;}
    if(centroid){if(exchange_xy) exchange(centroid->y,centroid->x);for(int i=1;i<=2;i++) if(flip[i-1]) (*centroid)(i)=1-(*centroid)(i);*centroid=*centroid*scaling+box.min_corner;}
    return volume*box.Size();
}
//#####################################################################
template bool Intersects(const RANGE<VECTOR<float,2> >&,const LINE_2D<float>&,const float);
template float Halfspace_Intersection_Size(const RANGE<VECTOR<float,2> >&,const LINE_2D<float>&,VECTOR<float,2>*);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Intersects(const RANGE<VECTOR<double,2> >&,const LINE_2D<double>&,const double);
template double Halfspace_Intersection_Size(const RANGE<VECTOR<double,2> >&,const LINE_2D<double>&,VECTOR<double,2>*);
#endif
};
};

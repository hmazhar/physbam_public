//#####################################################################
// Copyright 2002-2005, Robert Bridson, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POLYGON  
//##################################################################### 
#include <PhysBAM_Geometry/Basic_Geometry/POLYGON.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> POLYGON<TV>::
POLYGON()
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> POLYGON<TV>::
POLYGON(const int number_of_vertices) 
    :X(number_of_vertices)
{
}
//#####################################################################
// Constructor
//#####################################################################
namespace {
template<class T>
void Constructor_Helper(POLYGON<VECTOR<T,1> >& p, const RANGE<VECTOR<T,1> >& range)
{p.X(1)=VECTOR<T,1>(range.min_corner.x);p.X(2)=VECTOR<T,1>(range.max_corner.x);}

template<class T>
void Constructor_Helper(POLYGON<VECTOR<T,2> >& p, const RANGE<VECTOR<T,2> >& range)
{
    p.X(1)=VECTOR<T,2>(range.min_corner.x,range.min_corner.y); p.X(2)=VECTOR<T,2>(range.max_corner.x,range.min_corner.y);
    p.X(3)=VECTOR<T,2>(range.max_corner.x,range.max_corner.y); p.X(4)=VECTOR<T,2>(range.min_corner.x,range.max_corner.y);
}

template<class T>
void Constructor_Helper(POLYGON<VECTOR<T,3> >& p, const RANGE<VECTOR<T,3> >& range)
{
    p.X(1)=VECTOR<T,3>(range.min_corner.x,range.min_corner.y, range.min_corner.z); p.X(2)=VECTOR<T,3>(range.max_corner.x,range.min_corner.y, range.min_corner.z);
    p.X(3)=VECTOR<T,3>(range.max_corner.x,range.max_corner.y, range.min_corner.z); p.X(4)=VECTOR<T,3>(range.min_corner.x,range.max_corner.y, range.min_corner.z);
    p.X(5)=VECTOR<T,3>(range.min_corner.x,range.max_corner.y, range.max_corner.z); p.X(6)=VECTOR<T,3>(range.max_corner.x,range.max_corner.y, range.max_corner.z);
    p.X(7)=VECTOR<T,3>(range.max_corner.x,range.min_corner.y, range.max_corner.z); p.X(8)=VECTOR<T,3>(range.min_corner.x,range.min_corner.y, range.max_corner.z);
}
};
template<class TV> POLYGON<TV>::
POLYGON(const RANGE<TV>& range)
    :X(1<<TV::dimension)
{Constructor_Helper<T>(*this,range);}
//#####################################################################
// Function Area
//#####################################################################
// doesn't work if the polygon is self intersecting
namespace {
template <class T> T Area_Helper(const POLYGON<VECTOR<T,1> >& p)
{return abs(p.X(2).x-p.X(1).x);}

template <class T> T Area_Helper(const POLYGON<VECTOR<T,2> >& p)
{
    T area=0;
    for(int k=1;k<p.X.m;k++) area+=(p.X(k).y+p.X(k+1).y)*(p.X(k+1).x-p.X(k).x);
    area+=(p.X(p.X.m).y+p.X(1).y)*(p.X(1).x-p.X(p.X.m).x); // last edge
    return abs(area)/2; // should have done ((*y)(k)+(*y)(k+1))/2 above
}

template <class T> T Area_Helper(const POLYGON<VECTOR<T,3> >& p)
{PHYSBAM_FATAL_ERROR("Area_Helper<VECTOR<T,3> >(...) NOT IMPLEMENTED");}
};
template<class TV> typename TV::SCALAR POLYGON<TV>::
Area() const
{       
    return Area_Helper<T>(*this);
}
//#####################################################################
// Function Find_Closest_Point_On_Polygon
//#####################################################################
namespace {
template<class T> VECTOR<T,1> Find_Closest_Point_On_Polygon_Helper(const POLYGON<VECTOR<T,1> >& p,const VECTOR<T,1>& X_point,int& side)
{
    assert(p.X.Size()==2);
    if((p.X(1)-X_point).Magnitude_Squared() < (p.X(2)-X_point).Magnitude_Squared()){side=1;return p.X(1);}
    else{side=2;return p.X(2);}
}

template<class T> VECTOR<T,2> Find_Closest_Point_On_Polygon_Helper(const POLYGON<VECTOR<T,2> >& p,const VECTOR<T,2>& X_point,int& side)
{
    T distance_squared=FLT_MAX;side=0;VECTOR<T,2> result;

    // all sides except for the last one
    for(int k=1;k<p.X.m;k++){
        VECTOR<T,2> closest=SEGMENT_2D<T>(p.X(k),p.X(k+1)).Closest_Point_On_Segment(X_point);T d=(closest-X_point).Magnitude_Squared();
        if(distance_squared > d){distance_squared=d;result=closest;side=k;}}

    // last side - if the polygon is closed
    VECTOR<T,2> closest=SEGMENT_2D<T>(p.X(p.X.m),p.X(1)).Closest_Point_On_Segment(X_point);T d=(closest-X_point).Magnitude_Squared();
    if(distance_squared > d){result=closest;side=p.X.m;}
    return result;
}

template<class T> VECTOR<T,3> Find_Closest_Point_On_Polygon_Helper(const POLYGON<VECTOR<T,3> >& p,const VECTOR<T,3>& X_point,int& side)
{PHYSBAM_FATAL_ERROR("Find_Closest_Point_On_Polygon<VECTOR<T,3> >(...) NOT IMPLEMENTED");}
};
template<class TV> TV POLYGON<TV>::
Find_Closest_Point_On_Polygon(const TV& X_point,int& side) const
{return Find_Closest_Point_On_Polygon_Helper<T>(*this,X_point,side);}                  
//#####################################################################
// Function Distance_From_Polygon_To_Point
//#####################################################################
template<class TV> typename TV::SCALAR POLYGON<TV>::
Distance_From_Polygon_To_Point(const TV& X_point) const
{
    int side;return (X_point-Find_Closest_Point_On_Polygon(X_point,side)).Magnitude();
}
//#####################################################################
// Function Inside_Polygon
//#####################################################################
namespace {
template<class T> bool Inside_Polygon_Helper(const POLYGON<VECTOR<T,1> >& p, const VECTOR<T,1>& X_point)
{
    assert(p.X.Size()==2);
    return (p.X(1).x <= p.X(2).x && p.X(1).x <= X_point.x && X_point.x <= p.X(2).x) ||
           (p.X(1).x >= p.X(2).x && p.X(1).x >= X_point.x && X_point.x >= p.X(2).x);
}

template<class T> bool Inside_Polygon_Helper(const POLYGON<VECTOR<T,2> >& p, const VECTOR<T,2>& X_point)
{
    T theta_total=0;

    // all sides except for the last one
    for(int k=1;k<p.X.m;k++){
        VECTOR<T,2> X1=p.X(k)-X_point,X2=p.X(k+1)-X_point;
        if(X1==VECTOR<T,2>() || X2==VECTOR<T,2>()) return true; // (x,y) lies on the polygon
        T theta1=atan2(X1.y,X1.x),theta2=atan2(X2.y,X2.x); // atan2 returns values between -pi and pi, if (x,y) != (0,0)
        T theta=theta2-theta1;
        if(theta == (T)pi || theta == -(T)pi) return true; // (x,y) lies on the polygon
        if(theta > (T)pi) theta-=2*(T)pi;          // make sure the smaller angle is swept out
        else if(theta < -(T)pi) theta+=2*(T)pi; // make sure the smaller angle is swept out
        theta_total+=theta;}

    // last side
    VECTOR<T,2> X1=p.X(p.X.m)-X_point,X2=p.X(1)-X_point;
    if(X1==VECTOR<T,2>() || X2==VECTOR<T,2>()) return true; // (x,y) lies on the polygon
    T theta1=atan2(X1.y,X1.x),theta2=atan2(X2.y,X2.x); // atan2 returns values between -pi and pi, if (x,y) != (0,0)
    T theta=theta2-theta1;
    if(theta == (T)pi || theta == -(T)pi) return true; // (x,y) lies on the polygon
    if(theta > (T)pi) theta-=2*(T)pi;          // make sure the smaller angle is swept out
    else if(theta < -(T)pi) theta+=2*(T)pi; // make sure the smaller angle is swept out
    theta_total+=theta;

    // decide on inside or outside
    return abs(theta_total) >= (T)pi; // theta_total = +2*pi or -2*pi for a point inside the polygon
}

template<class T> bool Inside_Polygon_Helper(const POLYGON<VECTOR<T,3> >& p, const VECTOR<T,3>& X_point)
{PHYSBAM_FATAL_ERROR("Inside_Polygon_Helper<VECTOR<T,3> >(...) NOT IMPLEMENTED");}
};
template<class TV> bool POLYGON<TV>::
Inside_Polygon(const TV& X_point) const
{return Inside_Polygon_Helper<T>(*this,X_point);}
//#####################################################################
// Function Signed_Distance
//#####################################################################
template<class TV> typename TV::SCALAR POLYGON<TV>::
Signed_Distance(const TV& location) const
{
    T distance=Distance_From_Polygon_To_Point(location);
    return Inside_Polygon(location)?-distance:distance;
}
//#####################################################################
template class POLYGON<VECTOR<float,1> >;
template class POLYGON<VECTOR<float,2> >;
template class POLYGON<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class POLYGON<VECTOR<double,1> >;
template class POLYGON<VECTOR<double,2> >;
template class POLYGON<VECTOR<double,3> >;
#endif

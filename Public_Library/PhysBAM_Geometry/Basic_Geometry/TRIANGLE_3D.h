//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Neil Molino, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_3D
//##################################################################### 
#ifndef __TRIANGLE_3D__
#define __TRIANGLE_3D__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.h>
namespace PhysBAM{

template<class T>
class TRIANGLE_3D:public PLANE<T>
{
    typedef VECTOR<T,3> TV;
public:
    using PLANE<T>::x1;using PLANE<T>::normal;

    TV x2,x3; // x1 (in PLANE), x2 and x3 - clockwise order when looking at the plane

    TRIANGLE_3D()
    {
        Specify_Three_Points(TV(0,0,0),TV(0,1,0),TV(1,0,0));
    }

    TRIANGLE_3D(const TV& x1_input,const TV& x2_input,const TV& x3_input)
    {
        Specify_Three_Points(x1_input,x2_input,x3_input); 
    }

    template<class T_ARRAY>
    TRIANGLE_3D(const T_ARRAY& X_input)
    {
        STATIC_ASSERT(T_ARRAY::m==3);
        Specify_Three_Points(X_input(1),X_input(2),X_input(3));
    }

    void Specify_Three_Points(const TV& x1_input,const TV& x2_input,const TV& x3_input)
    {PLANE<T>::Specify_Three_Points(x1_input,x2_input,x3_input);x2=x2_input;x3=x3_input;}

    T Area() const 
    {return Area(x1,x2,x3);}
    
    static T Area(const TV& x1,const TV& x2,const TV& x3) // always positive for clockwise vertices: x1, x2, x3 
    {return (T).5*TV::Cross_Product(x2-x1,x3-x1).Magnitude();}
    
    static T Area_Squared(const TV& x1,const TV& x2,const TV& x3) // always positive for clockwise vertices: x1, x2, x3 
    {return (T).25*TV::Cross_Product(x2-x1,x3-x1).Magnitude_Squared();}

    T Size() const
    {return Area();}

    template<class T_ARRAY>
    static T Size(const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==3);return Area(X(1),X(2),X(3));}

    template<class T_ARRAY>
    static T Signed_Size(const T_ARRAY& X)
    {return Size(X);}

    T Aspect_Ratio() const
    {return Aspect_Ratio(x1,x2,x3);}

    static T Aspect_Ratio(const TV& x1_input,const TV& x2_input,const TV& x3_input)
    {TV u=x1_input-x2_input,v=x2_input-x3_input,w=x3_input-x1_input;
    T u2=TV::Dot_Product(u,u),v2=TV::Dot_Product(v,v),w2=TV::Dot_Product(w,w);
    return max(u2,v2,w2)/sqrt(TV::Cross_Product(u,v).Magnitude_Squared());}

    static T Minimum_Edge_Length(const TV& x1,const TV& x2,const TV& x3)
    {return sqrt(min((x2-x1).Magnitude_Squared(),(x3-x1).Magnitude_Squared(),(x3-x2).Magnitude_Squared()));}

    static T Maximum_Edge_Length(const TV& x1,const TV& x2,const TV& x3)
    {return sqrt(max((x2-x1).Magnitude_Squared(),(x3-x1).Magnitude_Squared(),(x3-x2).Magnitude_Squared()));}

    static T Minimum_Altitude(const TV& x1,const TV& x2,const TV& x3)
    {return 2*Area(x1,x2,x3)/Maximum_Edge_Length(x1,x2,x3);}
    
    static TV Barycentric_Coordinates(const TV& location,const TV& x1,const TV& x2,const TV& x3) // clockwise vertices
    {TV u=x2-x1,v=x3-x1,w=location-x1;
    T u_dot_u=TV::Dot_Product(u,u),v_dot_v=TV::Dot_Product(v,v),u_dot_v=TV::Dot_Product(u,v),
       u_dot_w=TV::Dot_Product(u,w),v_dot_w=TV::Dot_Product(v,w);
    T denominator=u_dot_u*v_dot_v-sqr(u_dot_v),one_over_denominator;
    if(abs(denominator)>(T)1e-16) one_over_denominator=1/denominator;else one_over_denominator=(T)1e16;
    T a=(v_dot_v*u_dot_w-u_dot_v*v_dot_w)*one_over_denominator,b=(u_dot_u*v_dot_w-u_dot_v*u_dot_w)*one_over_denominator;
    return TV(1-a-b,a,b);}
    
    static TV Clamped_Barycentric_Coordinates(const TV& location,const TV& x1,const TV& x2,const TV& x3,const T tolerance=1e-7) // clockwise vertices
    {TV u=x2-x1,v=x3-x1,w=location-x1;
    T u_dot_u=TV::Dot_Product(u,u),v_dot_v=TV::Dot_Product(v,v),u_dot_v=TV::Dot_Product(u,v),
       u_dot_w=TV::Dot_Product(u,w),v_dot_w=TV::Dot_Product(v,w);
    if(abs(u_dot_u)<tolerance){
        if(abs(v_dot_v)<tolerance) return TV((T)one_third,(T)one_third,(T)one_third); // single point
        T c=clamp(v_dot_w/v_dot_v,(T)0,(T)1);T a_and_b=(T).5*(1-c);return TV(a_and_b,a_and_b,c);} // x1 and x2 are a single point
    else if(abs(v_dot_v)<tolerance){
        T b=clamp(u_dot_w/u_dot_u,(T)0,(T)1);T a_and_c=(T).5*(1-b);return TV(a_and_c,b,a_and_c);} // x1 and x3 are a single point
    else{
        T denominator=u_dot_u*v_dot_v-sqr(u_dot_v); 
        if(abs(denominator)<tolerance){
            if(u_dot_v>0){ // u and v point in the same direction
                if(u_dot_u>u_dot_v){T b=clamp(u_dot_w/u_dot_u,(T)0,(T)1);return TV(1-b,b,0);}
                else{T c=clamp(v_dot_w/v_dot_v,(T)0,(T)1);return TV(1-c,0,c);}}
            else if(u_dot_w>0){T b=clamp(u_dot_w/u_dot_u,(T)0,(T)1);return TV(1-b,b,0);} // u and v point in opposite directions, and w is on the u segment
            else{T c=clamp(v_dot_w/v_dot_v,(T)0,(T)1);return TV(1-c,0,c);}} // u and v point in opposite directions, and w is on the v segment
        T one_over_denominator=1/denominator;
        T a=clamp((v_dot_v*u_dot_w-u_dot_v*v_dot_w)*one_over_denominator,(T)0,(T)1),b=clamp((u_dot_u*v_dot_w-u_dot_v*u_dot_w)*one_over_denominator,(T)0,(T)1);
        return TV(1-a-b,a,b);}}

    template<class T_ARRAY>
    static TV Barycentric_Coordinates(const TV& location,const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==3);return Barycentric_Coordinates(location,X(1),X(2),X(3));}

    template<class T_ARRAY>
    static TV Clamped_Barycentric_Coordinates(const TV& location,const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==3);return Clamped_Barycentric_Coordinates(location,X(1),X(2),X(3));}

    TV Sum_Barycentric_Coordinates(const TRIANGLE_3D<T>& embedded_triangle) const
    {return Barycentric_Coordinates(embedded_triangle.x1)+Barycentric_Coordinates(embedded_triangle.x2)+Barycentric_Coordinates(embedded_triangle.x3);}

    TV Barycentric_Coordinates(const TV& location) const 
    {return Barycentric_Coordinates(location,x1,x2,x3);}

    static TV Point_From_Barycentric_Coordinates(const TV& weights,const TV& x1,const TV& x2,const TV& x3) // clockwise vertices
    {return weights.x*x1+weights.y*x2+weights.z*x3;}

    template<class T_ARRAY>
    static TV Point_From_Barycentric_Coordinates(const TV& weights,const T_ARRAY& X) // clockwise vertices
    {STATIC_ASSERT(T_ARRAY::m==3);return weights.x*X(1)+weights.y*X(2)+weights.z*X(3);}

    TV Point_From_Barycentric_Coordinates(const TV& weights) const
    {return Point_From_Barycentric_Coordinates(weights,x1,x2,x3);}

    static TV Center(const TV& x1,const TV& x2,const TV& x3) // centroid
    {return (T)one_third*(x1+x2+x3);}

    TV Center() const // centroid
    {return Center(x1,x2,x3);}

    TV Incenter() const // intersection of angle bisectors
    {TV edge_lengths((x3-x2).Magnitude(),(x1-x3).Magnitude(),(x2-x1).Magnitude());T perimeter=edge_lengths.x+edge_lengths.y+edge_lengths.z;assert(perimeter>0);
    return Point_From_Barycentric_Coordinates(edge_lengths/perimeter);}

    bool Point_Face_Interaction(const TV& x,const TV& v,const INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,3>&> V_face,const T interaction_distance,T& distance,
        TV& interaction_normal,TV& weights,T& relative_speed,const bool allow_negative_weights,const bool exit_early) const
    {return Point_Face_Interaction(x,v,V_face(1),V_face(2),V_face(3),interaction_distance,distance,interaction_normal,weights,relative_speed,allow_negative_weights,exit_early);}

    bool Point_Face_Collision(const TV& x,const TV& v,const INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,3>&> V_face,const T dt,const T collision_thickness,T& collision_time,TV& normal,
        TV& weights,T& relative_speed,const bool exit_early) const
    {return Point_Face_Collision(x,v,V_face(1),V_face(2),V_face(3),dt,collision_thickness,collision_time,normal,weights,relative_speed,exit_early);}

    RANGE<TV> Bounding_Box() const
    {return RANGE<TV>::Bounding_Box(x1,x2,x3);}

    const TV& X(const int i) const
    {assert(1<=i && i<=3);
    switch(i){
        case 1: return x1;
        case 2: return x2;
        case 3: return x3;}
    PHYSBAM_FATAL_ERROR();}

    TV& X(const int i)
    {assert(1<=i && i<=3);
    switch(i){
        case 1: return x1;
        case 2: return x2;
        case 3: return x3;}
    PHYSBAM_FATAL_ERROR();}

//#####################################################################
    void Change_Size(const T delta);
    bool Inside(const TV& point,const T thickness_over_two=0) const;
    bool Point_Inside_Triangle(const TV& point,const T thickness_over_2=0) const;
    bool Planar_Point_Inside_Triangle(const TV& point,const T thickness_over_2=0) const;
    bool Lazy_Planar_Point_Inside_Triangle(const TV& point) const;
    T Minimum_Edge_Length() const;
    T Maximum_Edge_Length() const;
    int Region(const TV& location,int& region_id,const T tolerance) const;
    TV Closest_Point(const TV& location,TV& weights) const;
    T Distance_To_Triangle(const TV& location) const;
    T Minimum_Angle() const;
    T Maximum_Angle() const;
    T Signed_Solid_Angle(const TV& center) const;
    bool Point_Face_Interaction(const TV& x,const T interaction_distance,const bool allow_negative_weights,T& distance) const;
    void Point_Face_Interaction_Data(const TV& x,T& distance,TV& interaction_normal,TV& weights,const bool perform_attractions) const;    
    bool Point_Face_Interaction(const TV& x,const TV& v,const TV& v1,const TV& v2,const TV& v3,const T interaction_distance,
        T& distance,TV& interaction_normal,TV& weights,T& relative_speed,const bool allow_negative_weights,const bool exit_early) const;
    static POINT_SIMPLEX_COLLISION_TYPE Robust_Point_Triangle_Collision(const TRIANGLE_3D<T>& initial_triangle,const TRIANGLE_3D<T>& final_triangle,const TV& x,
        const TV& final_x,const T dt,const T collision_thickness,T& collision_time,TV& normal,TV& weights,T& relative_speed);
    bool Point_Face_Collision(const TV& x,const TV& v,const TV& v1,const TV& v2,const TV& v3,const T dt,const T collision_thickness,
        T& collision_time,TV& normal,TV& weights,T& relative_speed,const bool exit_early) const;
    void Clip_To_Box(const RANGE<TV>& box,ARRAY<TRIANGLE_3D<T> >& clipped_simplices) const;
    static void Cut_With_Hyperplane_And_Discard_Outside_Simplices(const TRIANGLE_3D<T>& triangle,const PLANE<T>& cutting_plane,ARRAY<TRIANGLE_3D<T> >& negative_triangles);
//#####################################################################
};

template<class T> std::ostream& operator<<(std::ostream& output,const TRIANGLE_3D<T>& triangle)
{output<<triangle.x1<<", "<<triangle.x2<<", "<<triangle.x3;return output;}

}
#endif

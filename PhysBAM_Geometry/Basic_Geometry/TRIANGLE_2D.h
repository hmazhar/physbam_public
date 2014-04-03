//#####################################################################
// Copyright 2003-2006, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Duc Nguyen, Avi Robinson-Mosher, Andrew Selle, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_2D
//#####################################################################
#ifndef __TRIANGLE_2D__
#define __TRIANGLE_2D__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class T>
class TRIANGLE_2D
{
    typedef VECTOR<T,2> TV;
public:
    VECTOR<TV,3> X;

    TRIANGLE_2D()
    {
        Specify_Three_Points(TV(0,0),TV(0,1),TV(1,0));
    }

    TRIANGLE_2D(const TV& x1_input,const TV& x2_input,const TV& x3_input)
    {
        Specify_Three_Points(x1_input,x2_input,x3_input);
    }

    template<class T_ARRAY>
    TRIANGLE_2D(const T_ARRAY& X_input)
        :X(X_input)
    {
        STATIC_ASSERT(T_ARRAY::m==3);
    }

    void Specify_Three_Points(const TV& x1_input,const TV& x2_input,const TV& x3_input)
    {X[1]=x1_input;X[2]=x2_input;X[3]=x3_input;}

    static T Signed_Area(const TV& x1,const TV& x2,const TV& x3)
    {return (T).5*TV::Cross_Product(x2-x1,x3-x1).x;}

    T Signed_Area() const
    {return Signed_Area(X[1],X[2],X[3]);}

    static T Area(const TV& x1,const TV& x2,const TV& x3)
    {return abs(Signed_Area(x1,x2,x3));}

    T Area() const
    {return abs(Signed_Area());}

    T Size() const
    {return Area();}

    T Signed_Size() const
    {return Signed_Area();}

    template<class T_ARRAY>
    static T Signed_Size(const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==3);return Signed_Area(X(1),X(2),X(3));}

    template<class T_ARRAY>
    static T Size(const T_ARRAY& X)
    {return abs(Signed_Size(X));}

    static bool Check_Orientation(const TV& x1,const TV& x2,const TV& x3)
    {return Signed_Area(x1,x2,x3)>=0;}

    bool Check_Orientation() const
    {return Signed_Area()>=0;}

    bool Fix_Orientation()
    {if(Check_Orientation()) return false;exchange(X[2],X[3]);return true;}

    template<class T_ARRAY>
    static T Half_Boundary_Measure(const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==3);return (T).5*((X(1)-X(2)).Magnitude()+(X(2)-X(3)).Magnitude()+(X(3)-X(1)).Magnitude());}

    T Aspect_Ratio() const
    {return Aspect_Ratio(X[1],X[2],X[3]);}

    static T Aspect_Ratio(const TV& x1_input,const TV& x2_input,const TV& x3_input)
    {TV u=x1_input-x2_input,v=x2_input-x3_input,w=x3_input-x1_input;
    T u2=TV::Dot_Product(u,u),v2=TV::Dot_Product(v,v),w2=TV::Dot_Product(w,w);
    return max(u2,v2,w2)/abs(TV::Cross_Product(u,v).x);}

    static T Minimum_Edge_Length(const TV& x1,const TV& x2,const TV& x3)
    {return sqrt(Minimum_Edge_Length_Squared(x1,x2,x3));}

    static T Minimum_Edge_Length_Squared(const TV& x1,const TV& x2,const TV& x3)
    {return min((x2-x1).Magnitude_Squared(),(x3-x1).Magnitude_Squared(),(x3-x2).Magnitude_Squared());}

    static T Maximum_Edge_Length(const TV& x1,const TV& x2,const TV& x3)
    {return sqrt(Maximum_Edge_Length_Squared(x1,x2,x3));}

    static T Maximum_Edge_Length_Squared(const TV& x1,const TV& x2,const TV& x3)
    {return max((x2-x1).Magnitude_Squared(),(x3-x1).Magnitude_Squared(),(x3-x2).Magnitude_Squared());}

    T Minimum_Altitude() const
    {return Minimum_Altitude(X[1],X[2],X[3]);}

    static T Minimum_Altitude(const TV& x1,const TV& x2,const TV& x3)
    {return 2*Area(x1,x2,x3)/Maximum_Edge_Length(x1,x2,x3);}

    static TV First_Two_Barycentric_Coordinates(const TV& location,const TV& x1,const TV& x2,const TV& x3)
    {return MATRIX<T,2>(x1-x3,x2-x3).Robust_Solve_Linear_System(location-x3);}

    static VECTOR<T,3> Barycentric_Coordinates(const TV& location,const TV& x1,const TV& x2,const TV& x3)
    {TV w=First_Two_Barycentric_Coordinates(location,x1,x2,x3);return VECTOR<T,3>(w.x,w.y,1-w.x-w.y);}

    template<class T_ARRAY>
    static VECTOR<T,3> Barycentric_Coordinates(const TV& location,const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==3);return Barycentric_Coordinates(location,X(1),X(2),X(3));}

    static TV Point_From_Barycentric_Coordinates(const VECTOR<T,3>& weights,const TV& x1,const TV& x2,const TV& x3)
    {return weights.x*x1+weights.y*x2+weights.z*x3;}

    static TV Point_From_Barycentric_Coordinates(const TV& weights,const TV& x1,const TV& x2,const TV& x3)
    {return weights.x*x1+weights.y*x2+(1-weights.x-weights.y)*x3;}

    TV Point_From_Barycentric_Coordinates(const VECTOR<T,3>& weights) const
    {return Point_From_Barycentric_Coordinates(weights,X[1],X[2],X[3]);}

    template<class T_ARRAY>
    static TV Point_From_Barycentric_Coordinates(const VECTOR<T,3>& weights,const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==3);return Point_From_Barycentric_Coordinates(weights,X(1),X(2),X(3));}

    VECTOR<T,3> Sum_Barycentric_Coordinates(const TRIANGLE_2D<T>& embedded_triangle) const
    {return Barycentric_Coordinates(embedded_triangle.X[1],X)+Barycentric_Coordinates(embedded_triangle.X[2],X)+Barycentric_Coordinates(embedded_triangle.X[3],X);}

    static TV Center(const TV& x1,const TV& x2,const TV& x3) // centroid
    {return (T)one_third*(x1+x2+x3);}

    TV Center() const // centroid
    {return Center(X[1],X[2],X[3]);}

    TV Incenter() const // intersection of angle bisectors
    {VECTOR<T,3> edge_lengths((X[3]-X[2]).Magnitude(),(X[1]-X[3]).Magnitude(),(X[2]-X[1]).Magnitude());T perimeter=edge_lengths.x+edge_lengths.y+edge_lengths.z;assert(perimeter>0);
    return Point_From_Barycentric_Coordinates(edge_lengths/perimeter);}

    static TV Circumcenter(const TV& x1,const TV& x2,const TV& x3)
    {TV x1x2=x2-x1,x1x3=x3-x1,m1=(T).5*(x1+x2),m2=(T).5*(x1+x3),m1m2=m2-m1,x1x2_perp(-x1x2.y,x1x2.x);
    return m1+x1x2_perp*(TV::Dot_Product(m1m2,x1x3)/TV::Dot_Product(x1x2_perp,x1x3));}

    static VECTOR<T,3> Circumcenter_Barycentric_Coordinates(const TV& x1,const TV& x2,const TV& x3)
    {TV a=x3-x2,b=x3-x1,c=x2-x1;T aa=a.Magnitude_Squared(),bb=b.Magnitude_Squared(),cc=c.Magnitude_Squared();
    VECTOR<T,3> weights(aa*(bb+cc-aa),bb*(cc+aa-bb),cc*(aa+bb-cc));return weights/(weights.x+weights.y+weights.z);}

    T Minimum_Angle() const
    {TV s1=(X[1]-X[2]).Normalized(),s2=(X[2]-X[3]).Normalized(),s3=(X[3]-X[1]).Normalized();
    return acos(max(TV::Dot_Product(s1,-s2),TV::Dot_Product(-s1,s3),TV::Dot_Product(s2,-s3)));}

    T Maximum_Angle() const
    {TV s1=(X[1]-X[2]).Normalized(),s2=(X[2]-X[3]).Normalized(),s3=(X[3]-X[1]).Normalized();
    return acos(min(TV::Dot_Product(s1,-s2),TV::Dot_Product(-s1,s3),TV::Dot_Product(s2,-s3)));}

    bool Outside(const TV& location,const T thickness_over_2=0) const
    {return Outside(location,X[1],X[2],X[3],thickness_over_2);}

    static bool Outside(const TV& location,const TV& x1,const TV& x2,const TV& x3,const T thickness_over_2=0)
    {assert(Check_Orientation(x1,x2,x3));TV location_minus_x1=location-x1;
    TV edge1=x2-x1;if(TV::Cross_Product(location_minus_x1,edge1).x>thickness_over_2*edge1.Magnitude()) return true;
    TV edge3=x1-x3;if(TV::Cross_Product(location-x3,edge3).x>thickness_over_2*edge3.Magnitude()) return true;
    TV edge2=x3-x2;if(TV::Cross_Product(location-x2,edge2).x>thickness_over_2*edge2.Magnitude()) return true;
    return false;}

    bool Inside(const TV& location,const T thickness_over_2=0) const
    {return !Outside(location,thickness_over_2);}

    RANGE<TV> Bounding_Box() const
    {return RANGE<TV>::Bounding_Box(X[1],X[2],X[3]);}

    static bool Check_Delaunay_Criterion(TV a,TV b,TV c,TV d)
    {assert(Check_Orientation(a,b,c) && Check_Orientation(d,c,b));b-=a;c-=a;d-=a;
    return VECTOR<T,3>::Triple_Product(b.Append(b.Magnitude_Squared()),c.Append(c.Magnitude_Squared()),d.Append(d.Magnitude_Squared()))>=0;}

//#####################################################################
    static T Negative_Material(const ARRAY<TV>& X,const ARRAY<T>& phis,const VECTOR<int,3>& indices);
    static void Cut_With_Hyperplane(ARRAY<TV>& X,const LINE_2D<T>& cutting_plane,const VECTOR<int,3>& indices,ARRAY<VECTOR<int,3> >& left_tris,
        ARRAY<VECTOR<int,3> >& right_tris);
    void Clip_To_Box(const RANGE<TV>& box,ARRAY<TRIANGLE_2D<T> >& clipped_simplices) const;
    void Cut_With_Hyperplane_And_Discard_Outside_Simplices(const TRIANGLE_2D<T>& triangle,const LINE_2D<T>& cutting_plane,ARRAY<TRIANGLE_2D<T> >& negative_triangles) const;
private:
    static void Cut_Simplex(ARRAY<TV>& X,const VECTOR<int,3>& indices,const VECTOR<TV,3>& X_nodes,const VECTOR<T,3>& phi_nodes,
        ARRAY<VECTOR<int,3> >& left_tris,ARRAY<VECTOR<int,3> >& right_tris);
//#####################################################################
};

template<class T> std::ostream& operator<<(std::ostream& output,const TRIANGLE_2D<T>& triangle)
{return output<<triangle.X;}

}
#endif


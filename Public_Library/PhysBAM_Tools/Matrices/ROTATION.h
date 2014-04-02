//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Igor Neverov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ROTATION
//#####################################################################
#ifndef __ROTATION__
#define __ROTATION__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Tools/Vectors/COMPLEX.h>
namespace PhysBAM{

using ::std::atan2;

template<class TV> class ROTATION;

template<class TV> struct IS_SCALAR_BLOCK<ROTATION<TV> > {static const bool value=(TV::m>1) && IS_SCALAR_BLOCK<TV>::value;};
template<class TV,class RW> struct IS_BINARY_IO_SAFE<ROTATION<TV>,RW> {static const bool value=(TV::m>1) && IS_BINARY_IO_SAFE<TV,RW>::value;};
template<class TV> struct HAS_CHEAP_COPY<ROTATION<TV> > {static const bool value=true;};

//#####################################################################
// 1D
//#####################################################################
template<class T> 
class ROTATION<VECTOR<T,1> >
{
    typedef VECTOR<T,1> TV;
    typedef VECTOR<T,0> T_SPIN;
public:
    typedef T SCALAR;

    bool operator==(const ROTATION<TV>& r) const
    {return true;}

    bool operator!=(const ROTATION<TV>& r) const
    {return false;}

    ROTATION<TV>& operator*=(const ROTATION<TV>& r)
    {return *this;}

    ROTATION<TV> operator*(const ROTATION<TV>& r) const
    {return *this;}

    T Normalize()
    {return 1;}

    ROTATION<TV> Normalized() const
    {return *this;}

    ROTATION<TV> Inverse() const
    {return *this;}

    const TV& Inverse_Rotate(const TV& x) const
    {return x;}

    const TV& Rotate(const TV& x) const
    {return x;}

    const T_SPIN& Rotate_Spin(const T_SPIN& spin) const
    {return spin;}

    const TWIST<TV>& Rotate(const TWIST<TV>& twist) const
    {return twist;}

    VECTOR<T,0> Euler_Angles() const
    {return VECTOR<T,0>();}

    static ROTATION From_Euler_Angles(const VECTOR<T,0>&)
    {return ROTATION();}

    VECTOR<T,0> Rotation_Vector() const
    {return VECTOR<T,0>();}

    MATRIX<T,1> Rotation_Matrix() const
    {return MATRIX<T,1>(1);}

    TV Rotated_X_Axis() const // Q*(1)
    {return TV(1);}

    TV Rotated_Axis(const int axis) const
    {assert(axis==1);return TV(1);}

    static ROTATION<TV> From_Rotation_Vector(const VECTOR<T,0>)
    {return ROTATION<TV>();}

    static ROTATION<TV> From_Rotated_Vector(const TV&,const TV&)
    {return ROTATION<TV>();}

    static ROTATION<TV> Spherical_Linear_Interpolation(const ROTATION<TV>,const ROTATION<TV>,const T)
    {return ROTATION<TV>();}

    bool Is_Normalized() const
    {return true;}

    T Angle() const
    {return 0;}

    ROTATION Scale_Angle() const
    {return *this;}

    static ROTATION Average_Rotation(const ARRAY<ROTATION>&)
    {return ROTATION();}

};

//#####################################################################
// 2D
//#####################################################################
template<class T> 
class ROTATION<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<T,1> T_SPIN;

public:
    COMPLEX<T> c;

private:
    ROTATION(const COMPLEX<T>& c2)
        :c(c2)
    {}

public:
    typedef T SCALAR;

    ROTATION()
        :c(1,0)
    {}

    explicit ROTATION(const MATRIX<T,2>& A)
        :c(A.Column(1).Normalized())
    {}

    const COMPLEX<T>& Complex() const
    {return c;}

    static ROTATION<TV> From_Complex(const COMPLEX<T>& c2)
    {return ROTATION<TV>(c2).Normalized();}

    bool operator==(const ROTATION<TV>& r) const
    {return c==r.c;}

    bool operator!=(const ROTATION<TV>& r) const
    {return c!=r.c;}

    ROTATION<TV>& operator*=(const ROTATION<TV>& r)
    {c*=r.c;return *this;}

    ROTATION<TV> operator*(const ROTATION<TV>& r) const
    {return ROTATION<TV>(c*r.c);}

    ROTATION<TV> Inverse() const
    {return ROTATION<TV>(c.Conjugated());}

    TV Rotate(const TV& v) const
    {return TV(c.re*v.x-c.im*v.y,c.im*v.x+c.re*v.y);}

    TV Inverse_Rotate(const TV& v) const
    {return TV(c.re*v.x+c.im*v.y,c.re*v.y-c.im*v.x);}

    const T_SPIN& Rotate_Spin(const T_SPIN& spin) const
    {return spin;}

    TWIST<TV> Rotate(const TWIST<TV>& twist) const
    {return TWIST<TV>(Rotate(twist.linear),twist.angular);}

    T Normalize()
    {return c.Normalize();}

    ROTATION<TV> Normalized() const
    {ROTATION<TV> r(*this);r.Normalize();return r;}

    bool Is_Normalized(const T tolerance=(T)1e-3) const
    {return abs(c.Magnitude_Squared()-(T)1)<=tolerance;}

    void Get_Rotated_Frame(TV& x_axis,TV& y_axis) const
    {assert(Is_Normalized());x_axis=TV(c.re,c.im);y_axis=TV(-c.im,c.re);}

    T Angle() const
    {return atan2(c.im,c.re);}

    VECTOR<T,1> Euler_Angles() const
    {return Rotation_Vector();}

    VECTOR<T,1> Rotation_Vector() const
    {return VECTOR<T,1>(Angle());}

    MATRIX<T,2> Rotation_Matrix() const
    {return MATRIX<T,2>(c.re,c.im,-c.im,c.re);}

    TV Rotated_X_Axis() const // Q*(1,0)
    {return TV(c.re,c.im);}

    TV Rotated_Y_Axis() const // Q*(0,1)
    {return TV(-c.im,c.re);}

    TV Rotated_Axis(const int axis) const
    {assert(axis>=1 && axis<=2);if(axis==1) return Rotated_X_Axis();return Rotated_Y_Axis();}

    static ROTATION<TV> From_Angle(const T& a)
    {return ROTATION<TV>(COMPLEX<T>::Unit_Polar(a));}

    static ROTATION<TV> From_Rotation_Vector(const VECTOR<T,1>& v)
    {return From_Angle(v.x);}

    static ROTATION<TV> From_Euler_Angles(const VECTOR<T,1>& angle)
    {return From_Rotation_Vector(angle);}

//#####################################################################
    static ROTATION<TV> From_Rotated_Vector(const TV& v1,const TV& v2);
    ROTATION<TV> Scale_Angle(const T a) const;
    static ROTATION<TV> Spherical_Linear_Interpolation(const ROTATION<TV>& r1,const ROTATION<TV>& r2,const T t);
    static ROTATION<TV> Average_Rotation(const ARRAY<ROTATION<TV> >& rotations);
};

//#####################################################################
// 3D
//#####################################################################
template<class T> 
class ROTATION<VECTOR<T,3> >
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<T,3> T_SPIN;

public:
    QUATERNION<T> q;

private:
    ROTATION(const QUATERNION<T>& q2)
        :q(q2)
    {}

    ROTATION(const T s,const T x,const T y,const T z)
        :q(s,x,y,z)
    {}

    class UNUSABLE {};
public:
    typedef T SCALAR;

    ROTATION()
        :q(1,0,0,0)
    {}

    template<class T2> explicit ROTATION(const ROTATION<T2>& r)
        :q(r.Quaternion())
    {
        STATIC_ASSERT(!IS_SAME<T,T2>::value);
        if(!IS_SAME<T,int>::value) Normalize();
    }

    ROTATION(const T angle,const TV& direction);
    explicit ROTATION(const MATRIX<T,3>& A); // matches A with a quaternion

private:
    template<class T2>
    static ROTATION<TV> From_Components_Helper(const T2 s,const T x,const T y,const T z)
    {return ROTATION<TV>(s,x,y,z).Normalized();}

    static ROTATION<TV> From_Components_Helper(const int s,const int x,const int y,const int z)
    {return ROTATION<TV>((T)s,(T)x,(T)y,(T)z);}
public:

    static ROTATION<TV> From_Components(const T s,const T x,const T y,const T z)
    {return From_Components_Helper(s,x,y,z);}

    const QUATERNION<T>& Quaternion() const
    {return q;}

    static ROTATION<TV> From_Quaternion(const QUATERNION<T>& q)
    {return ROTATION<TV>(q).Normalized();}

    bool operator==(const ROTATION<TV>& r) const
    {return q==r.q || q==-r.q;}

    bool operator!=(const ROTATION<TV>& r) const
    {return !(*this==r);}

    ROTATION<TV>& operator*=(const ROTATION<TV>& r)
    {q*=r.q;return *this;}

    QUATERNION<T> operator*(const QUATERNION<T>& w) const
    {return q*w;}

    ROTATION<TV> operator*(const ROTATION<TV>& r) const
    {return ROTATION<TV>(q*r.q);}

    T Normalize()
    {return q.Normalize();}

    ROTATION<TV> Normalized() const
    {ROTATION<TV> r(*this);r.Normalize();return r;}

    bool Is_Normalized(const T tolerance=(T)1e-4) const
    {return q.Is_Normalized(tolerance);}

    ROTATION<TV> Inverse() const
    {return ROTATION<TV>(q.Conjugate());}

    static ROTATION<TV> From_Euler_Angles(const TV& euler_angles)
    {return From_Euler_Angles(euler_angles.x,euler_angles.y,euler_angles.z);}

    TV Rotate(const TV& v) const // 20 mult and 13 add/sub
    {assert(Is_Normalized());T two_s=q.s+q.s;return two_s*TV::Cross_Product(q.v,v)+(two_s*q.s-(T)1)*v+(T)2*TV::Dot_Product(q.v,v)*q.v;}

    TV Inverse_Rotate(const TV& v) const // 20 mult and 13 add/sub
    {assert(Is_Normalized());T two_s=q.s+q.s;return two_s*TV::Cross_Product(v,q.v)+(two_s*q.s-(T)1)*v+(T)2*TV::Dot_Product(q.v,v)*q.v;}

    T_SPIN Rotate_Spin(const T_SPIN& spin) const
    {return Rotate(spin);}

    TWIST<TV> Rotate(const TWIST<TV>& twist) const
    {return TWIST<TV>(Rotate(twist.linear),Rotate(twist.angular));}

    TV Euler_Angles() const
    {TV euler_angles;Euler_Angles(euler_angles.x,euler_angles.y,euler_angles.z);return euler_angles;}

    TV Rotated_Axis(const int axis) const
    {assert(axis>=1 && axis<=3);if(axis==1) return Rotated_X_Axis();if(axis==2) return Rotated_Y_Axis();return Rotated_Z_Axis();}

//#####################################################################
    TV Rotation_Vector() const;
    static ROTATION<TV> From_Rotation_Vector(const TV& v);
    static ROTATION<TV> From_Euler_Angles(const T euler_angle_x,const T euler_angle_y,const T euler_angle_z); // rotation about fixed axes x, then y, then z
    void Euler_Angles(T& euler_angle_x,T& euler_angle_y,T& euler_angle_z) const;
    TV Rotated_X_Axis() const; // Q*(1,0,0)
    TV Rotated_Y_Axis() const; // Q*(0,1,0)
    TV Rotated_Z_Axis() const; // Q*(0,0,1)
    void Get_Rotated_Frame(TV& x_axis,TV& y_axis,TV& z_axis) const;
    void Get_Angle_Axis(T& angle,TV& axis) const;
    T Angle() const;
    TV Get_Axis() const;
    MATRIX<T,3> Rotation_Matrix() const; // 18 mult and 12 add/sub
    ROTATION<TV> Scale_Angle(const T a) const;
    static ROTATION<TV> Spherical_Linear_Interpolation(const ROTATION<TV>& r1,const ROTATION<TV>& r2,const T t);
    static ROTATION<TV> Average_Rotation(const ARRAY<ROTATION<TV> >& rotations);
    static ROTATION<TV> From_Rotated_Vector(const TV& initial_vector,const TV& final_vector);
};
//#####################################################################

template<> inline VECTOR<int,3> ROTATION<VECTOR<int,3> >::
Rotate(const VECTOR<int,3>& v) const // homogenous of degree 2 in q, since we can't usefully assume normalization for integer case
{return 2*q.s*VECTOR<int,3>::Cross_Product(q.v,v)+(q.s*q.s-q.v.Magnitude_Squared())*v+2*VECTOR<int,3>::Dot_Product(q.v,v)*q.v;}

template<> inline VECTOR<int,3> ROTATION<VECTOR<int,3> >::
Inverse_Rotate(const VECTOR<int,3>& v) const // homogenous of degree 2 in q, since we can't usefully assume normalization for integer case
{return 2*q.s*VECTOR<int,3>::Cross_Product(v,q.v)+(q.s*q.s-q.v.Magnitude_Squared())*v+2*VECTOR<int,3>::Dot_Product(q.v,v)*q.v;}
}
#endif

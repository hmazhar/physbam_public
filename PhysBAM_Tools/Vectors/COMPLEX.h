//#####################################################################
// Copyright 2003-2007, Eran Guendelman, Geoffrey Irving, Craig Schroeder, Andrew Selle, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMPLEX
//#####################################################################
#ifndef __COMPLEX__
#define __COMPLEX__

#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
namespace PhysBAM{

using ::std::sin;
using ::std::cos;

template<class T> struct IS_SCALAR_BLOCK<COMPLEX<T> >:public IS_SCALAR_BLOCK<T>{};
template<class T> struct IS_SCALAR_VECTOR_SPACE<COMPLEX<T> >:public IS_SCALAR_VECTOR_SPACE<T>{};
template<class T,class RW> struct IS_BINARY_IO_SAFE<COMPLEX<T>,RW>:public IS_BINARY_IO_SAFE<T,RW>{};

template<class T>
class COMPLEX
{
public:
    typedef T SCALAR;

    T re,im;

    COMPLEX()
        :re(0),im(0)
    {}

    COMPLEX(T re_input,T im_input)
        :re(re_input),im(im_input)
    {}

    template<class T2> explicit COMPLEX(const COMPLEX<T2>& complex_input)
        :re((T)complex_input.re),im((T)complex_input.im)
    {}

    explicit COMPLEX(const VECTOR<T,2>& input)
        :re(input.x),im(input.y)
    {}

    VECTOR<T,2> Vector() const
    {return VECTOR<T,2>(re,im);}

    static COMPLEX<T> One()
    {return COMPLEX(1,0);}

    bool operator==(const COMPLEX<T>& c) const
    {return re==c.re && im==c.im;}

    bool operator!=(const COMPLEX<T>& c) const
    {return re!=c.re || im!=c.im;}

    COMPLEX<T>& operator*=(const COMPLEX<T>& c)
    {T old_re=re;re=re*c.re-im*c.im;im=old_re*c.im+im*c.re;return *this;}

    COMPLEX<T>& operator*=(const T a)
    {re*=a;im*=a;return *this;}

    COMPLEX<T> operator*(const COMPLEX<T>& c) const
    {return COMPLEX<T>(re*c.re-im*c.im,re*c.im+im*c.re);}

    COMPLEX<T> operator*(const T a) const
    {return COMPLEX<T>(a*re,a*im);}

    COMPLEX<T>& operator+=(const COMPLEX<T>& c)
    {re+=c.re;im+=c.im;return *this;}

    COMPLEX<T> operator+(const COMPLEX<T>& c) const
    {return COMPLEX<T>(re+c.re,im+c.im);}

    COMPLEX<T> operator+(const T& a) const
    {return COMPLEX<T>(re+a,im);}

    COMPLEX<T> operator-=(const COMPLEX<T>& c)
    {re-=c.re;im-=c.im;return *this;}

    COMPLEX<T> operator-(const COMPLEX<T>& c) const
    {return COMPLEX<T>(re-c.re,im-c.im);}

    COMPLEX<T> operator-(const T a) const
    {return COMPLEX<T>(re-a,im);}

    T Magnitude_Squared() const
    {return sqr(re)+sqr(im);}

    T Magnitude() const
    {return sqrt(sqr(re)+sqr(im));}

    void Conjugate()
    {im*=-1;}

    COMPLEX<T> Conjugated() const
    {return COMPLEX<T>(re,-im);}

    COMPLEX<T> Inverse() const
    {assert(re!=0 || im!=0);T denominator=(T)1/(re*re+im*im);return COMPLEX<T>(re*denominator,-im*denominator);}

    COMPLEX<T> Sqrt() const
    {T magnitude=Magnitude();return COMPLEX<T>(sqrt((T).5*(magnitude+re)),sign(im)*sqrt((T).5*(magnitude-re)));}

    static T Dot_Product(const COMPLEX<T>& c1,const COMPLEX<T>& c2)
    {return c1.re*c2.re+c1.im*c2.im;}

    T Normalize()
    {T magnitude=Magnitude();if(magnitude) *this*=1/magnitude;else{re=1;im=0;}return magnitude;}

    COMPLEX<T> Normalized() const
    {COMPLEX<T> c(*this);c.Normalize();return c;}

    bool Is_Normalized(const T tolerance=(T)1e-3) const
    {return abs(Magnitude_Squared()-(T)1)<=tolerance;}

    COMPLEX<T> Rotated_Counter_Clockwise_90() const
    {return COMPLEX<T>(-im,re);}

    COMPLEX<T> Rotated_Clockwise_90() const
    {return COMPLEX<T>(im,-re);}

    static COMPLEX<T> Polar(const T r,const T theta)   // r*e^(i*theta) = r(cos(theta)+i*sin(theta))
    {return COMPLEX<T>(r*cos(theta),r*sin(theta));}

    static COMPLEX<T> Unit_Polar(const T theta)
    {return COMPLEX<T>(cos(theta),sin(theta));}

//#####################################################################
};
template<class T>
inline COMPLEX<T> operator*(const T a,const COMPLEX<T>& c)
{return COMPLEX<T>(a*c.re,a*c.im);}
}
#endif

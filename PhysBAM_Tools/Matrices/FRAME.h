//#####################################################################
// Copyright 2002-2008, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Neil Molino, Igor Neverov, Andrew Selle, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRAME
//#####################################################################
#ifndef __FRAME__
#define __FRAME__

#include <PhysBAM_Tools/Math_Tools/ONE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/COMPLEX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV> struct HAS_CHEAP_COPY<FRAME<TV> > {static const bool value=true;};
template<class TV> struct IS_SCALAR_BLOCK<FRAME<TV> >:public IS_SCALAR_BLOCK<TV>{};
template<class TV,class RW> struct IS_BINARY_IO_SAFE<FRAME<TV>,RW>:public IS_BINARY_IO_SAFE<TV,RW>{};

template<class TV>
class FRAME
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
public:
    typedef T SCALAR;

    TV t; // defaults to 0
    ROTATION<TV> r; // defaults to 1

    FRAME()
    {}

    explicit FRAME(const TV& t)
        :t(t)
    {}

    explicit FRAME(const ROTATION<TV>& r)
        :r(r)
    {}

    FRAME(const TV& t,const ROTATION<TV>& r)
        :t(t),r(r)
    {}

    explicit FRAME(const MATRIX<T,d+1>& m_input)
        :t(m_input.Translation()),r(m_input.Extract_Rotation())
    {}

    template<class T2> explicit FRAME(const FRAME<VECTOR<T2,d> >& f)
        :t(f.t),r(f.r)
    {}

    bool operator==(const FRAME& f) const
    {return t==f.t && r==f.r;}

    bool operator!=(const FRAME& f) const
    {return t!=f.t || r!=f.r;}

    FRAME& operator*=(const FRAME& f)
    {t+=r.Rotate(f.t);r*=f.r;return *this;}

    FRAME operator*(const FRAME& f) const
    {return FRAME(t+r.Rotate(f.t),r*f.r);}

    TV operator*(const TV& v) const
    {return t+r.Rotate(v);}

    void Invert()
    {*this=Inverse();}

    FRAME Inverse() const
    {ROTATION<TV> r_inverse=r.Inverse();return FRAME(-r_inverse.Rotate(t),r_inverse);}

    TV Inverse_Times(const TV& v) const
    {return r.Inverse_Rotate(v-t);}

    FRAME Inverse_Times(const FRAME& f) const
    {return FRAME(r.Inverse_Rotate(f.t-t),r.Inverse()*f.r);}

    static FRAME Interpolation(const FRAME& f1,const FRAME& f2,const T s)
    {return FRAME((1-s)*f1.t+s*f2.t,ROTATION<TV>::Spherical_Linear_Interpolation(f1.r,f2.r,s));}

    MATRIX<T,d+1> Matrix() const
    {MATRIX<T,d+1> matrix=MATRIX<T,d+1>::From_Linear(r.Rotation_Matrix());matrix.Translation()=t;return matrix;}

    static FRAME<TV> From_Delta(const TWIST<TV>& twist)
    {return FRAME<TV>(twist.linear,ROTATION<TV>::From_Rotation_Vector(twist.angular));}

    TWIST<TV> Delta() const
    {return TWIST<TV>(t,r.Rotation_Vector());}

    const FRAME& Frame() const
    {return *this;}
    
    const TV& X() const
    {return t;}

    const ROTATION<TV>& Rotation() const
    {return r;}

    std::string Name() const {return Static_Name();}
    static std::string Static_Name()
    {return STRING_UTILITIES::string_sprintf("FRAME<VECTOR<T,%d> >",TV::m);}

//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ZERO
//#####################################################################
#ifndef __ZERO__
#define __ZERO__

namespace PhysBAM{

struct ZERO
{
    bool operator!() const
    {return true;}

    ZERO operator-() const
    {return ZERO();}

    ZERO operator-(const ZERO) const
    {return ZERO();}

    ZERO operator+(const ZERO) const
    {return ZERO();}

    ZERO operator*(const ZERO) const
    {return ZERO();}

    template<class T>
    ZERO operator*=(const T&)
    {return ZERO();}

    template<class T>
    ZERO operator/=(const T&)
    {return ZERO();}

    ZERO operator+=(const ZERO)
    {return ZERO();}

    ZERO operator-=(const ZERO)
    {return ZERO();}


    bool operator==(const ZERO) const
    {return true;}

    bool operator!=(const ZERO) const
    {return false;}
    
//#####################################################################
};

inline bool operator<(const float x,const ZERO)
{return x<0;}

inline bool operator<(const double x,const ZERO)
{return x<0;}

inline bool operator>(const float x,const ZERO)
{return x>0;}

inline bool operator>(const double x,const ZERO)
{return x>0;}

template<class T> inline const T& operator+(const T& x,const ZERO)
{return x;}

template<class T> inline const T& operator+(const ZERO,const T& x)
{return x;}

template<class T> inline const T& operator-(const T& x,const ZERO)
{return x;}

template<class T> inline T operator-(const ZERO,const T& x)
{return -x;}

template<class T> inline ZERO operator*(const ZERO,const T&)
{return ZERO();}

template<class T> inline ZERO operator*(const T&,const ZERO)
{return ZERO();}

}
#endif

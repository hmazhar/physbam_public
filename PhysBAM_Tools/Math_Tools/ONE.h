//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ONE
//#####################################################################
#ifndef __ONE__
#define __ONE__

namespace PhysBAM{

struct ONE
{
    bool operator!() const
    {return false;}

    ONE operator*(const ONE) const
    {return ONE();}

    bool operator==(const ONE) const
    {return true;}

    ONE Inverse() const
    {return ONE();}

    static ONE One()
    {return ONE();}

//#####################################################################
};

template<class T> inline const T& operator*(const T& x,const ONE)
{return x;}

template<class T> inline const T& operator*(const ONE,const T& x)
{return x;}

template<class T> inline const T& operator/(const T& x,const ONE)
{return x;}

template<class T> inline T& operator*=(T& x,const ONE)
{return x;}

template<class T> inline T& operator/=(T& x,const ONE)
{return x;}

}
#endif

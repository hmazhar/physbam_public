//#####################################################################
// Copyright 2008, Eftychios Sifakis, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXACT_FLOAT
//#####################################################################
#ifndef __EXACT_FLOAT__
#define __EXACT_FLOAT__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include <PhysBAM_Tools/Math_Tools/ZERO.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>

namespace PhysBAM{

template<class T> class EXACT_RATIONAL;

template<class T>
class EXACT_FLOAT
{
    STATIC_ASSERT(IS_FLOATING_POINT<T>::value);
    friend class EXACT_RATIONAL<T>;

    ARRAY<T> expansion;

public:
    typedef T FP_PRIMITIVE;

    explicit EXACT_FLOAT(T base_float=0)
    {expansion.Append(base_float);}

    EXACT_FLOAT operator-(const EXACT_FLOAT& exact_float) const
    {return *this+(-exact_float);}

    EXACT_FLOAT& operator+=(const EXACT_FLOAT& exact_float)
    {EXACT_FLOAT result;result=*this+exact_float;expansion.Exchange(result.expansion);return *this;}

    EXACT_FLOAT& operator-=(const EXACT_FLOAT& exact_float)
    {return *this+=-exact_float;}

    EXACT_FLOAT& operator*=(const EXACT_FLOAT& exact_float)
    {EXACT_FLOAT result;result=*this*exact_float;expansion.Exchange(result.expansion);return *this;}

    EXACT_FLOAT operator-() const
    {EXACT_FLOAT result(*this);result.expansion*=-1;return result;}

    bool operator<(ZERO zero) const
    {return sign(*this)<0;}

    bool operator<(const EXACT_FLOAT& exact_float) const
    {return sign(*this-exact_float)<0;}

    bool operator<=(ZERO zero) const
    {return sign(*this)<=0;}

    bool operator<=(const EXACT_FLOAT& exact_float) const
    {return sign(*this-exact_float)<=0;}

    bool operator>=(ZERO zero) const
    {return sign(*this)>=0;}

    bool operator>=(const EXACT_FLOAT& exact_float) const
    {return sign(*this-exact_float)>=0;}

    bool operator>(ZERO zero) const
    {return sign(*this)>0;}

    bool operator>(const EXACT_FLOAT& exact_float) const
    {return sign(*this-exact_float)>0;}

    bool operator==(ZERO zero) const
    {return sign(*this)==0;}

    bool operator==(const EXACT_FLOAT& exact_float) const
    {return sign(*this-exact_float)==0;}

    bool operator!=(ZERO zero) const
    {return sign(*this)!=0;}

    bool operator!=(const EXACT_FLOAT& exact_float) const
    {return sign(*this-exact_float)!=0;}

    friend inline int sign(const EXACT_FLOAT& exact_float)
    {const T& msf=exact_float.expansion.Last();if(msf>0)return 1;if(msf<0)return -1;return 0;}

    T Compress_And_Estimate() const
    {Compress();return expansion.Last();}

private:
    EXACT_FLOAT& Negate()
    {expansion*=-1;return *this;}

//#####################################################################
public:
    EXACT_FLOAT operator+(const EXACT_FLOAT& exact_float) const;
    EXACT_FLOAT operator*(const EXACT_FLOAT& exact_float) const;
private:
    EXACT_FLOAT operator*(T base_float) const;
    const EXACT_FLOAT& Compress() const;
//#####################################################################
};

template<class T> EXACT_RATIONAL<T>
operator/(const EXACT_FLOAT<T>& exact_float1,const EXACT_FLOAT<T>& exact_float2)
{return EXACT_RATIONAL<T>(exact_float1,exact_float2);}

template<class T> std::ostream&
operator<<(std::ostream& output,const EXACT_FLOAT<T>& exact_float)
{return output<<exact_float.Estimate();}

}
#endif

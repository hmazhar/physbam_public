//#####################################################################
// Copyright 2008, Eftychios Sifakis, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXACT_RATIONAL
//#####################################################################
#ifndef __EXACT_RATIONAL__
#define __EXACT_RATIONAL__

#include <PhysBAM_Tools/Adaptive_Arithmetic/EXACT_FLOAT.h>

namespace PhysBAM{

template<class T>
class EXACT_RATIONAL
{
    EXACT_FLOAT<T> numerator;
    EXACT_FLOAT<T> denominator; // Assumed to be positive

public:
    typedef T FP_PRIMITIVE;
    struct DIVISION_BY_ZERO{};

    explicit EXACT_RATIONAL(T base_float=0)
        :numerator(base_float),denominator(1)
    {}

    EXACT_RATIONAL(const EXACT_FLOAT<T>& numerator_input,const EXACT_FLOAT<T>& denominator_input)
        : numerator(numerator_input),denominator(denominator_input)
    {
        if(denominator==ZERO()) throw DIVISION_BY_ZERO();
        if(denominator<ZERO()){numerator.Negate();denominator.Negate();}
        Normalize();
    }

    EXACT_RATIONAL operator+(const EXACT_RATIONAL& exact_rational) const
    {EXACT_RATIONAL result(*this);result+=exact_rational;return result;}
    
    EXACT_RATIONAL operator-(const EXACT_RATIONAL& exact_rational) const
    {EXACT_RATIONAL result(*this);result-=exact_rational;return result;}
    
    EXACT_RATIONAL operator*(const EXACT_RATIONAL& exact_rational) const
    {EXACT_RATIONAL result(*this);result*=exact_rational;return result;}
    
    EXACT_RATIONAL operator/(const EXACT_RATIONAL& exact_rational) const
    {EXACT_RATIONAL result(*this);result/=exact_rational;return result;}
    
    bool operator<(ZERO zero) const
    {return sign(numerator)<0;}

    bool operator==(ZERO zero) const
    {return sign(numerator)==0;}

    friend inline int sign(const EXACT_RATIONAL& exact_rational)
    {return sign(exact_rational.numerator);}

//#####################################################################
    EXACT_RATIONAL& operator+=(const EXACT_RATIONAL& exact_rational);
    EXACT_RATIONAL& operator-=(const EXACT_RATIONAL& exact_rational);
    EXACT_RATIONAL& operator*=(const EXACT_RATIONAL& exact_rational);
    EXACT_RATIONAL& operator/=(const EXACT_RATIONAL& exact_rational);
    T Compress_And_Estimate() const;
    void Normalize() const;
//#####################################################################
};

template<class T> inline std::ostream&
operator<<(std::ostream& output, const EXACT_RATIONAL<T>& exact_rational)
{return output<<exact_rational.Estimate();}

}
#endif

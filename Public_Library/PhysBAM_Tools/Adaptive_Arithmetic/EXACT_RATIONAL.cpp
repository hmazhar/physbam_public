//#####################################################################
// Copyright 2008, Eftychios Sifakis, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Adaptive_Arithmetic/EXACT_RATIONAL.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <limits>
#undef _DEBUG_EXACT_RATIONAL_
using namespace PhysBAM;
//#####################################################################
// operator +=
//#####################################################################
template<class T> inline EXACT_RATIONAL<T>& EXACT_RATIONAL<T>::
operator+=(const EXACT_RATIONAL& exact_rational)
{
    if(denominator==exact_rational.denominator) numerator+=exact_rational.numerator;
    else{numerator*=exact_rational.denominator;numerator+=exact_rational.numerator*denominator;denominator*=exact_rational.denominator;Normalize();}
    return *this;
}
//#####################################################################
// operator -=
//#####################################################################
template<class T> inline EXACT_RATIONAL<T>& EXACT_RATIONAL<T>::
operator-=(const EXACT_RATIONAL& exact_rational)
{
    if(denominator==exact_rational.denominator) numerator-=exact_rational.numerator;
    else{numerator*=exact_rational.denominator;numerator-=exact_rational.numerator*denominator;denominator*=exact_rational.denominator;Normalize();}
    return *this;
}
//#####################################################################
// operator *=
//#####################################################################
template<class T> inline EXACT_RATIONAL<T>& EXACT_RATIONAL<T>::
operator*=(const EXACT_RATIONAL& exact_rational)
{
    numerator*=exact_rational.numerator;denominator*=exact_rational.denominator;Normalize();
    return *this;
}
//#####################################################################
//operator /=
//#####################################################################
template<class T> inline EXACT_RATIONAL<T>& EXACT_RATIONAL<T>::
operator/=(const EXACT_RATIONAL& exact_rational)
{
    if(exact_rational==ZERO()) throw DIVISION_BY_ZERO();
    if(denominator==exact_rational.denominator) denominator=exact_rational.numerator;else{numerator*=exact_rational.denominator;denominator*=exact_rational.numerator;}
    if(exact_rational<ZERO()){numerator.Negate();denominator.Negate();}
    Normalize();
    return *this;
}
//#####################################################################
// Function Compress_And_Estimate
//#####################################################################
template<class T> inline T EXACT_RATIONAL<T>::
Compress_And_Estimate() const
{
    return numerator.Compress_And_Estimate()/denominator.Compress_And_Estimate();
}
//#####################################################################
// Function Normalize
//#####################################################################
// TODO: clean up and determine sources of overflow

namespace{
    template<class T> struct FLOATING_POINT_HELPER;
    
    template<>
    struct FLOATING_POINT_HELPER<float>
    {
        inline static float exp(int exp)
        {
            STATIC_ASSERT((sizeof(float)==sizeof(unsigned int)&&std::numeric_limits<float>::is_iec559));
            PHYSBAM_ASSERT(-126<=exp&&exp<=127);
            union {unsigned int i;float f;} u;
            u.i=(static_cast<unsigned int>(exp+127)<<23);
            return u.f;
        }
    
        inline static int log(float x)
        {
            STATIC_ASSERT((sizeof(float)==sizeof(unsigned int)&&std::numeric_limits<float>::is_iec559));
            PHYSBAM_ASSERT(x!=0);
            union {unsigned int i;float f;} u;
            u.f=x;
            return (static_cast<int>((u.i>>23)&0xFF)-127);
        }
    };
    
    template<>
    struct FLOATING_POINT_HELPER<double>
    {
        inline static double exp(int exp)
        {
            STATIC_ASSERT((sizeof(double)==sizeof(unsigned long long)&&std::numeric_limits<double>::is_iec559));
            PHYSBAM_ASSERT(-1022<=exp&&exp<=1023);
            union {unsigned long long i;double f;} u;
            u.i=(static_cast<unsigned long long>(exp+1023)<<52);
            return u.f;
        }
    
        inline static int log(double x)
        {
            STATIC_ASSERT((sizeof(double)==sizeof(unsigned long long)&&std::numeric_limits<double>::is_iec559));
            PHYSBAM_ASSERT(x!=0);
            union {unsigned long long i;double f;} u;
            u.f=x;
            return (static_cast<int>((u.i>>52)&0x7FF)-1023);
        }
    };
}

template<class T> inline void EXACT_RATIONAL<T>::
Normalize() const
{
    if(numerator==ZERO()){const_cast<EXACT_FLOAT<T>&>(denominator)=EXACT_FLOAT<T>(1);return;}
    numerator.Compress();denominator.Compress();
#ifdef _DEBUG_EXACT_RATIONAL_
    // TODO: determine whether normalization is necessary.
    int max_exp=floating_point_exponent(denominator.expansion.Last());
    int min_exp=floating_point_exponent(denominator.expansion(1));
    T scale=FLOATING_POINT_HELPER<T>::exp(-(max_exp+min_exp)/2);
    const_cast<EXACT_FLOAT<T>&>(numerator).expansion*=scale;
    const_cast<EXACT_FLOAT<T>&>(denominator).expansion*=scale;

    max_exp=max(FLOATING_POINT_HELPER<T>::log(numerator.expansion.Last()),FLOATING_POINT_HELPER<T>::log(denominator.expansion.Last()));
    min_exp=min(FLOATING_POINT_HELPER<T>::log(numerator.expansion(1)),FLOATING_POINT_HELPER<T>::log(denominator.expansion(1)));
    if(max_exp>=std::numeric_limits<T>::max_exponent/4||min_exp<=std::numeric_limits<T>::min_exponent/4)
        LOG::cerr<<"*** DEBUG WARNING ***\n"<<"Large exponents encountered in EXACT_RATIONAL<T>::Normalize()\n"<<"    this = "<<this<<"\n"<<"    this.numerator.expansion = "<<numerator.expansion
            <<"    this.denominator.expansion = "<<denominator.expansion<<std::endl;
#endif
}
//#####################################################################
template class EXACT_RATIONAL<float>;
template class EXACT_RATIONAL<double>;


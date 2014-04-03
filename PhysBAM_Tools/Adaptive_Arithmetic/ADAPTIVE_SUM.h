//#####################################################################
// Copyright 2008, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADAPTIVE_SUM__
#define __ADAPTIVE_SUM__
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_BASE.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_FORWARD.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_WRAPPER_BASE.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/EXACT_ARITHMETIC_POLICY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <limits>

namespace PhysBAM{
namespace ADAPTIVE_DETAIL{

using PhysBAM::EXACT_ARITHMETIC_POLICY;

template<class T_ADAPTIVE1,class T_ADAPTIVE2>
struct IS_ADAPTIVE<ADAPTIVE_SUM<T_ADAPTIVE1,T_ADAPTIVE2> > {static const bool value=true;};

template<class T_ADAPTIVE1,class T_ADAPTIVE2>
class ADAPTIVE_SUM:public ADAPTIVE_BASE<ADAPTIVE_SUM<T_ADAPTIVE1,T_ADAPTIVE2>,typename EXACT_ARITHMETIC_POLICY<typename T_ADAPTIVE1::EXACT_TYPE,typename T_ADAPTIVE2::EXACT_TYPE>::SUM_TYPE>
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE1>::value));
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE2>::value));

    const T_ADAPTIVE1 a1;
    const T_ADAPTIVE2 a2;

public:
    typedef ADAPTIVE_SUM<T_ADAPTIVE1,T_ADAPTIVE2> THIS;
    typedef ADAPTIVE_BASE<THIS,typename EXACT_ARITHMETIC_POLICY<typename T_ADAPTIVE1::EXACT_TYPE,typename T_ADAPTIVE2::EXACT_TYPE>::SUM_TYPE > BASE;
    typedef typename BASE::EXACT_TYPE EXACT_TYPE;
    typedef typename BASE::FP_TYPE FP_TYPE;

    ADAPTIVE_SUM(const T_ADAPTIVE1& a1_input,const T_ADAPTIVE2& a2_input)
        :a1(a1_input),a2(a2_input)
    {}

private:
    ADAPTIVE_SUM();
    THIS& operator=(const THIS&);
    
public:
    PAIR<FP_TYPE,FP_TYPE> Estimate_And_Error_Implementation() const
    {const FP_TYPE eps=std::numeric_limits<FP_TYPE>::epsilon()/2;
    FP_TYPE a1_estimate,a1_error;a1.Estimate_And_Error().Get(a1_estimate,a1_error);
    FP_TYPE a2_estimate,a2_error;a2.Estimate_And_Error().Get(a2_estimate,a2_error);
    FP_TYPE estimate=a1_estimate+a2_estimate;
    FP_TYPE error=eps*fabs(estimate)+a1_error+a2_error;
    error+=2*eps*error;
    return PAIR<FP_TYPE,FP_TYPE>(estimate,error);}

    EXACT_TYPE Exact_Implementation() const
    {return a1.Exact()+a2.Exact();}

    int Quick_Sign_Implementation() const
    {int a1_sign=a1.Quick_Sign();
    if(a1_sign==ADAPTIVE_SIGN_UNKNOWN) return ADAPTIVE_SIGN_UNKNOWN;
    int a2_sign=a2.Quick_Sign();
    if(a2_sign==ADAPTIVE_SIGN_UNKNOWN)return ADAPTIVE_SIGN_UNKNOWN;
    int sign=a1_sign+a2_sign;
    if(sign>0) return ADAPTIVE_SIGN_POSITIVE;
    else if (sign<0) return ADAPTIVE_SIGN_NEGATIVE;
    else if(a1_sign==0) return ADAPTIVE_SIGN_ZERO;
    else return ADAPTIVE_SIGN_UNKNOWN;}
//#####################################################################
};
}
}
#endif

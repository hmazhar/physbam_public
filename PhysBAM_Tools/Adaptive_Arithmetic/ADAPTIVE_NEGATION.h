//#####################################################################
// Copyright 2008, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADAPTIVE_DETAIL_ADAPTIVE_NEGATION_H__
#define __ADAPTIVE_DETAIL_ADAPTIVE_NEGATION_H__

#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_BASE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>

#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>

namespace PhysBAM{
namespace ADAPTIVE_DETAIL{

template<class T_ADAPTIVE>
struct IS_ADAPTIVE<ADAPTIVE_NEGATION<T_ADAPTIVE> > {static const bool value=true;};

template<class T_ADAPTIVE>
class ADAPTIVE_NEGATION:public ADAPTIVE_BASE<ADAPTIVE_NEGATION<T_ADAPTIVE>,typename T_ADAPTIVE::EXACT_TYPE>
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));

public:
    typedef ADAPTIVE_NEGATION<T_ADAPTIVE> THIS;
    typedef ADAPTIVE_BASE<THIS,typename T_ADAPTIVE::EXACT_TYPE> BASE;

    typedef typename BASE::EXACT_TYPE EXACT_TYPE;
    typedef typename BASE::FP_TYPE FP_TYPE;

    ADAPTIVE_NEGATION(const T_ADAPTIVE& a0);

private:
    const T_ADAPTIVE a0;

    ADAPTIVE_NEGATION();
    THIS& operator=(const THIS&);

public:
    PAIR<FP_TYPE,FP_TYPE> Estimate_And_Error_Implementation() const
    {PAIR<FP_TYPE,FP_TYPE> estimate_and_error=a0.Estimate_And_Error();
    estimate_and_error.x*=-1;
    return estimate_and_error;}

    FP_TYPE Refined_Estimate() const
    {return -a0.Refined_Estimate();}

    EXACT_TYPE Exact_Implementation() const
    {return -a0.Exact();}

    int Sign_Implementation() const
    {return -a0.Sign();}

    int Quick_Sign_Implementation() const
    {int sign=a0.Quick_Sign();
    return (sign==ADAPTIVE_SIGN_UNKNOWN?ADAPTIVE_SIGN_UNKNOWN:-sign);}
//#####################################################################
};
}
using ADAPTIVE_DETAIL::ADAPTIVE_NEGATION;
}
#endif

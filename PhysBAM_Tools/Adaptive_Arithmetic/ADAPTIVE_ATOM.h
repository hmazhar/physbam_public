//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADAPTIVE_ATOM__
#define __ADAPTIVE_ATOM__

#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_BASE.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_FORWARD.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <limits>

namespace PhysBAM{
namespace ADAPTIVE_DETAIL{

template<class T,class EXACT_TYPE>
struct IS_ADAPTIVE<ADAPTIVE_ATOM<T,EXACT_TYPE> > {static const bool value=true;};

template<class T, class EXACT_TYPE> struct IS_ATOMIZABLE
    :public AND<IS_CONVERTIBLE<T,typename EXACT_TYPE::FP_PRIMITIVE>::value,
    std::numeric_limits<T>::digits<=std::numeric_limits<typename EXACT_TYPE::FP_PRIMITIVE>::digits>
{};

template<class T,class EXACT_TYPE_T>
class ADAPTIVE_ATOM
    :public ADAPTIVE_BASE<ADAPTIVE_ATOM<T,EXACT_TYPE_T>,EXACT_TYPE_T>
{
    STATIC_ASSERT((IS_ATOMIZABLE<T,EXACT_TYPE_T>::value));

public:
    typedef EXACT_TYPE_T EXACT_TYPE;
    typedef ADAPTIVE_BASE<ADAPTIVE_ATOM,EXACT_TYPE> BASE;

    using BASE::Sign;
    typedef typename BASE::FP_TYPE FP_TYPE;

    ADAPTIVE_ATOM()
        :value(0)
    {}

    ADAPTIVE_ATOM(T value_input)
        :value(value_input)
    {}

private:
    T value;

public:
    PAIR<FP_TYPE,FP_TYPE > Estimate_And_Error_Implementation() const
    {return PAIR<FP_TYPE,FP_TYPE >(static_cast<FP_TYPE>(value),static_cast<FP_TYPE>(0));}

    FP_TYPE Refined_Estimate_Implementation() const
    {return static_cast<FP_TYPE>(value);}

    EXACT_TYPE Exact_Implementation() const
    {return EXACT_TYPE(static_cast<FP_TYPE>(value));}

    int Sign_Implementation() const
    {return (value>0?1:value<0?-1:0);}

    int Quick_Sign_Implementation() const
    {return Sign();}

//#####################################################################
};
}
using ADAPTIVE_DETAIL::IS_ATOMIZABLE;
using ADAPTIVE_DETAIL::ADAPTIVE_ATOM;
}
#endif

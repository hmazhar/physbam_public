//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADAPTIVE_WRAPPER_BASE__
#define __ADAPTIVE_WRAPPER_BASE__

#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_BASE.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_FORWARD.h>
#include <limits>

namespace PhysBAM{
namespace ADAPTIVE_DETAIL{

template<class DERIVED,class T_WRAPPED,CACHED_TYPE cached_type>
class ADAPTIVE_WRAPPER_BASE
    :public ADAPTIVE_BASE<DERIVED,typename T_WRAPPED::EXACT_TYPE>
{
    STATIC_ASSERT((IS_ADAPTIVE<T_WRAPPED>::value));

public:
    typedef T_WRAPPED WRAPPED_TYPE;
    typedef ADAPTIVE_BASE<DERIVED,typename WRAPPED_TYPE::EXACT_TYPE> BASE;

    using BASE::Derived;
    typedef typename BASE::EXACT_TYPE EXACT_TYPE;
    typedef typename BASE::FP_TYPE FP_TYPE;

    WRAPPED_TYPE Wrapped() const
    {return Derived().Wrapped_Implementation();}

protected:
    ADAPTIVE_WRAPPER_BASE()
    {}

    ADAPTIVE_WRAPPER_BASE(const ADAPTIVE_WRAPPER_BASE& other)
    {}

    ~ADAPTIVE_WRAPPER_BASE()
    {}

private:
    ADAPTIVE_WRAPPER_BASE& operator=(const ADAPTIVE_WRAPPER_BASE&);

public:
    PAIR<FP_TYPE,FP_TYPE> Estimate_And_Error_Implementation() const
    {return Wrapped().Estimate_And_Error();}

    FP_TYPE Refined_Estimate_Implementation() const
    {return Wrapped().Refined_Estimate();}

    EXACT_TYPE Exact_Implementation() const
    {return Wrapped().Exact();}

    int Sign_Implementation() const
    {return Wrapped().Sign();}

    int Quick_Sign_Implementation() const 
    {return Wrapped().Quick_Sign();}

//#####################################################################
};

template<class DERIVED,class T_WRAPPED>
class ADAPTIVE_WRAPPER_BASE<DERIVED,T_WRAPPED,CACHED>
    :public ADAPTIVE_BASE<DERIVED,typename T_WRAPPED::EXACT_TYPE>
{
    STATIC_ASSERT((IS_ADAPTIVE<T_WRAPPED>::value));

public:
    typedef T_WRAPPED WRAPPED_TYPE;
    typedef ADAPTIVE_BASE<DERIVED,typename WRAPPED_TYPE::EXACT_TYPE> BASE;

    using BASE::Derived;
    using BASE::Sign;
    typedef typename BASE::EXACT_TYPE EXACT_TYPE;
    typedef typename BASE::FP_TYPE FP_TYPE;

    const WRAPPED_TYPE& Wrapped() const
    {return Derived().Wrapped_Implementation();}

protected:
    ADAPTIVE_WRAPPER_BASE()
        :cached_error(-1),cached_exact(0),cached_sign(ADAPTIVE_SIGN_UNKNOWN)
    {}

    ~ADAPTIVE_WRAPPER_BASE()
    {delete cached_exact;}

private:
    mutable FP_TYPE cached_estimate;
    mutable FP_TYPE cached_error;
    mutable EXACT_TYPE* cached_exact;
    mutable int cached_sign;

    void Compute_Exact_If_Necessary() const
    {if(!cached_exact){
        cached_exact=new EXACT_TYPE(Wrapped().Exact());
        cached_estimate=cached_exact->Compress_And_Estimate();
        const FP_TYPE epsilon=(FP_TYPE).5*std::numeric_limits<FP_TYPE>::epsilon();
        cached_error=epsilon*cached_estimate;}}

    ADAPTIVE_WRAPPER_BASE& operator=(const ADAPTIVE_WRAPPER_BASE&);

public:
    PAIR<FP_TYPE,FP_TYPE> Estimate_And_Error_Implementation() const
    {if(cached_error<0) Wrapped().Estimate_And_Error().Get(cached_estimate,cached_error);
    return PAIR<FP_TYPE,FP_TYPE>(cached_estimate,cached_error);}

    EXACT_TYPE Exact_Implementation() const
    {Compute_Exact_If_Necessary();return *cached_exact;}

    FP_TYPE Refined_Estimate_Implementation() const
    {Compute_Exact_If_Necessary();return cached_estimate;}

    int Sign_Implementation() const
    {if(cached_sign==ADAPTIVE_SIGN_UNKNOWN && (cached_sign=Wrapped().Quick_Sign())==ADAPTIVE_SIGN_UNKNOWN)
        cached_sign=BASE::Sign_Implementation();
    return cached_sign;}

    int Quick_Sign_Implementation() const
    {return Sign();}

//#####################################################################
};
}
}
#endif

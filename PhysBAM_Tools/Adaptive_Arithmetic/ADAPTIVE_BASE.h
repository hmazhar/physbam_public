//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADAPTIVE_BASE__
#define __ADAPTIVE_BASE__

#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>

namespace PhysBAM{
namespace ADAPTIVE_DETAIL{


template<class T,class T_EXACT> class ADAPTIVE_ATOM;

template<class DERIVED,class EXACT_TYPE_T>
class ADAPTIVE_BASE
{
public:
    typedef EXACT_TYPE_T EXACT_TYPE;
    typedef typename EXACT_TYPE::FP_PRIMITIVE FP_TYPE;

    DERIVED& Derived()
    {return *static_cast<DERIVED*>(this);}

    const DERIVED& Derived() const
    {return *static_cast<const DERIVED*>(this);}

    PAIR<FP_TYPE,FP_TYPE> Estimate_And_Error() const
    {return Derived().Estimate_And_Error_Implementation();}

    FP_TYPE Refined_Estimate() const
    {return Derived().Refined_Estimate_Implementation();}

    EXACT_TYPE Exact() const
    {return Derived().Exact_Implementation();}

    int Sign() const
    {return Derived().Sign_Implementation();}

    int Quick_Sign() const
    {return Derived().Quick_Sign_Implementation();}

protected:
    ADAPTIVE_BASE()
    {}

    ADAPTIVE_BASE(const ADAPTIVE_BASE&)
    {}

    inline FP_TYPE Refined_Estimate_Implementation() const
    {return Exact().Compress_And_Estimate();}
    
    int Sign_Implementation() const
    {FP_TYPE estimate, error;
    Estimate_And_Error().Get(estimate,error);
    if(estimate>error) return ADAPTIVE_SIGN_POSITIVE;if(-estimate>error) return ADAPTIVE_SIGN_NEGATIVE;
    if(error==0) return ADAPTIVE_SIGN_ZERO; // estimate has to be zero, too, since the previous conditionals failed
    estimate=Refined_Estimate();
    return (estimate>0?ADAPTIVE_SIGN_POSITIVE:estimate<0?ADAPTIVE_SIGN_NEGATIVE:ADAPTIVE_SIGN_ZERO);}
    
    int Quick_Sign_Implementation() const
    {return ADAPTIVE_SIGN_UNKNOWN;}

//#####################################################################
};
}
using ADAPTIVE_DETAIL::ADAPTIVE_BASE;
using ADAPTIVE_DETAIL::IS_ADAPTIVE;
using ADAPTIVE_DETAIL::IS_NOT_ADAPTIVE;
}
#endif

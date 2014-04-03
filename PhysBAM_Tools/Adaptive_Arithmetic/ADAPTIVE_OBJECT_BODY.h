//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADAPTIVE_OBJECT_BODY__
#define __ADAPTIVE_OBJECT_BODY__

#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_WRAPPER_BASE.h>
#include <cassert>

namespace PhysBAM{
namespace ADAPTIVE_DETAIL{

template<class T_ADAPTIVE,class T_EXACT>
class ADAPTIVE_OBJECT_BODY
    :public ADAPTIVE_DYNAMIC_BASE<T_EXACT>,public ADAPTIVE_WRAPPER_BASE<ADAPTIVE_OBJECT_BODY<T_ADAPTIVE,T_EXACT>,T_ADAPTIVE,CACHED>
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE>::value));
    STATIC_ASSERT((IS_SAME<T_EXACT,typename T_ADAPTIVE::EXACT_TYPE>::value));

public:
    typedef T_EXACT EXACT_TYPE;
    typedef ADAPTIVE_DYNAMIC_BASE<EXACT_TYPE> BASE_DYNAMIC;
    typedef ADAPTIVE_WRAPPER_BASE<ADAPTIVE_OBJECT_BODY<T_ADAPTIVE,EXACT_TYPE>,T_ADAPTIVE,CACHED> BASE_WRAPPER;

    using BASE_DYNAMIC::reference_count;
    typedef typename BASE_WRAPPER::FP_TYPE FP_TYPE;
    typedef typename BASE_WRAPPER::WRAPPED_TYPE WRAPPED_TYPE;

    ADAPTIVE_OBJECT_BODY(const WRAPPED_TYPE& wrapped)
        :wrapped(wrapped)
    {}

    ~ADAPTIVE_OBJECT_BODY()
    {
        assert(reference_count==0);
    }

    PAIR<FP_TYPE,FP_TYPE> Estimate_And_Error() const
    {return BASE_WRAPPER::Estimate_And_Error();}

    FP_TYPE Refined_Estimate() const
    {return BASE_WRAPPER::Refined_Estimate();}

    EXACT_TYPE Exact() const
    {return BASE_WRAPPER::Exact();}

    int Sign() const
    {return BASE_WRAPPER::Sign();}

    int Quick_Sign() const
    {return BASE_WRAPPER::Quick_Sign();}

private:
    WRAPPED_TYPE wrapped;

    ADAPTIVE_OBJECT_BODY();
    ADAPTIVE_OBJECT_BODY(const ADAPTIVE_OBJECT_BODY&);
    ADAPTIVE_OBJECT_BODY& operator=(const ADAPTIVE_OBJECT_BODY&);

public:
    const WRAPPED_TYPE& Wrapped_Implementation() const {return wrapped;}
//#####################################################################
};
}
}
#endif

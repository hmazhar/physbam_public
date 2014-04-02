//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADAPTIVE_OBJECT__
#define __ADAPTIVE_OBJECT__

#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_DYNAMIC_BASE.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_OBJECT_BODY.h>

namespace PhysBAM{
namespace ADAPTIVE_DETAIL{

template<class T_EXACT>
class ADAPTIVE_OBJECT
    :public ADAPTIVE_BASE<ADAPTIVE_OBJECT<T_EXACT>,T_EXACT>
{
public:
    typedef T_EXACT EXACT_TYPE;
    typedef ADAPTIVE_BASE<ADAPTIVE_OBJECT,EXACT_TYPE> BASE;

    typedef typename BASE::FP_TYPE FP_TYPE;

    ADAPTIVE_OBJECT()
        :body(0)
    {}

    template<class DERIVED>
    ADAPTIVE_OBJECT(const ADAPTIVE_BASE<DERIVED,EXACT_TYPE>& expr)
        :body(new ADAPTIVE_OBJECT_BODY<DERIVED,EXACT_TYPE>(expr.Derived()))
    {}

    ADAPTIVE_OBJECT(const ADAPTIVE_OBJECT& other)
        :body(other.body)
    {
        if(body) body->reference_count++;
    }

    ~ADAPTIVE_OBJECT()
    {
        Release();
    }

    ADAPTIVE_OBJECT&
    Swap(ADAPTIVE_OBJECT& other)
    {exchange(body,other.body);return *this;}

    ADAPTIVE_OBJECT&
    operator=(const ADAPTIVE_OBJECT& other)
    {ADAPTIVE_OBJECT temp(other);return Swap(temp);}

    template<class T_ADAPTIVE>
    ADAPTIVE_OBJECT&
    operator=(const ADAPTIVE_BASE<T_ADAPTIVE,EXACT_TYPE>& expr)
    {ADAPTIVE_OBJECT temp(expr.Derived());return Swap(temp);}

    template<class T_ADAPTIVE>
    ADAPTIVE_OBJECT&
    operator+=(const ADAPTIVE_BASE<T_ADAPTIVE,EXACT_TYPE>& expr)
    {return (*this=*this+expr);}

    void Release()
    {if(body && --(body->reference_count)==0) delete body;body=0;}

private:
    const ADAPTIVE_DYNAMIC_BASE<EXACT_TYPE>* body;

public:
    PAIR<FP_TYPE,FP_TYPE> Estimate_And_Error_Implementation() const
    {assert(body);return body->Estimate_And_Error();}

    FP_TYPE Refined_Estimate_Implementation() const
    {assert(body);return body->Refined_Estimate();}

    EXACT_TYPE Exact_Implementation() const
    {assert(body);return body->Exact();}

    int Sign_Implementation() const
    {assert(body);return body->Sign();}

    int Quick_Sign_Implementation() const
    {assert(body);return body->Quick_Sign();}

//#####################################################################
};

template<class EXACT_TYPE>
struct IS_ADAPTIVE<ADAPTIVE_OBJECT<EXACT_TYPE> > {static const bool value=true;};

}
using ADAPTIVE_DETAIL::ADAPTIVE_OBJECT;
}
#endif

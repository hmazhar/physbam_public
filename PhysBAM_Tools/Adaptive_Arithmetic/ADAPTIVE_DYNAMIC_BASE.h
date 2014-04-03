//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADAPTIVE_DYNAMIC_BASE__
#define __ADAPTIVE_DYNAMIC_BASE__
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>

namespace PhysBAM{
namespace ADAPTIVE_DETAIL{

template<class EXACT_TYPE_T>
class ADAPTIVE_DYNAMIC_BASE:public NONCOPYABLE
{
public:
    typedef EXACT_TYPE_T EXACT_TYPE;
    typedef typename EXACT_TYPE::FP_PRIMITIVE FP_TYPE;

protected:
    ADAPTIVE_DYNAMIC_BASE()
        :reference_count(1)
    {}

public:
    virtual ~ADAPTIVE_DYNAMIC_BASE()
    {}

    mutable int reference_count;

//#####################################################################
    virtual PAIR<FP_TYPE,FP_TYPE> Estimate_And_Error() const=0;
    virtual FP_TYPE Refined_Estimate() const=0;
    virtual EXACT_TYPE Exact() const=0;
    virtual int Sign() const=0;
    virtual int Quick_Sign() const=0;
//#####################################################################
};
}
}
#endif

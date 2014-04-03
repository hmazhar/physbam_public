//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOBOL
//#####################################################################
#ifndef __SOBOL__
#define __SOBOL__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <limits>
namespace PhysBAM{

template<class TV> class RANGE;
template<class TV>
class SOBOL:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    STATIC_ASSERT(std::numeric_limits<T>::radix==2);
    static const int max_bits=std::numeric_limits<T>::digits;
    typedef typename IF<(max_bits<=30),unsigned long,unsigned long long>::TYPE TI; // pick an integer large enough to hold T's mantissa
    STATIC_ASSERT(std::numeric_limits<TI>::digits>max_bits);
    enum WORKAROUND {d=TV::m};
private:
    const TV offset,scales;
    ARRAY<VECTOR<TI,d> > v; // direction numbers
    VECTOR<TI,d> x; // last result
    TI n;
public:

    SOBOL(const RANGE<TV>& box);
    ~SOBOL();

//#####################################################################
    TV Get_Vector();
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACTORIAL
//#####################################################################
#ifndef __FACTORIAL__
#define __FACTORIAL__

#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
namespace PhysBAM{

template<unsigned d> struct FACTORIAL;

template<> struct FACTORIAL<0>{enum WORKAROUND {value=1};};

template<unsigned d> struct FACTORIAL
{
    STATIC_ASSERT((d<=12));
    enum WORKAROUND {value=d*FACTORIAL<d-1>::value};
};
}
#endif

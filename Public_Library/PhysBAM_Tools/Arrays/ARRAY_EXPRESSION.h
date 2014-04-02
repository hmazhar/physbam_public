//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_EXPRESSION
//#####################################################################
#ifndef __ARRAY_EXPRESSION__
#define __ARRAY_EXPRESSION__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
namespace PhysBAM{

template<class T,class T_ARRAY,class ID> // ID = int
class ARRAY_EXPRESSION:public ARRAY_BASE<T,T_ARRAY,ID>
{
public:
    typedef const T RESULT_TYPE;
    typedef const T CONST_RESULT_TYPE;
//#####################################################################
};
}
#endif

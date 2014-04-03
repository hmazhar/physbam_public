//#####################################################################
// Copyright 2006-2009, Geoffrey Irving, Nipun Kwatra, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header ARRAYS_FORWARD
//#####################################################################
#ifndef __ARRAYS_FORWARD__
#define __ARRAYS_FORWARD__

#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
namespace PhysBAM{

template<class T,class ID=int> class ARRAY_VIEW;
template<class T,class T_ARRAY,class ID=int> class ARRAY_BASE;
template<class T,class T_ARRAY,class ID=int> class ARRAY_EXPRESSION;
template<class T,class T_ARRAY,class ID=int> class INDEX_ITERATED_ARRAY;

template<class T,class ID=int> class ARRAY;

template<class ID=int> class IDENTITY_ARRAY;
template<class T,class ID=int> class CONSTANT_ARRAY;

template<class T_ARRAY,class T_INDICES=ARRAY<int>&> class INDIRECT_ARRAY;

template<class T_ARRAY,class T_PROJECTOR> class PROJECTED_ARRAY;
template<class T_STRUCT,class T_FIELD,T_FIELD T_STRUCT::* field> struct FIELD_PROJECTOR;
struct INDEX_PROJECTOR;

template<class T_ARRAY,class T_NEW> struct REBIND;
template<class T,class T_NEW> struct REBIND<ARRAY<T>,T_NEW>{typedef ARRAY<T_NEW> TYPE;};
template<class T,class T_NEW,class T_ARRAY,class ID> struct REBIND<ARRAY_BASE<T,T_ARRAY,ID>,T_NEW>{typedef ARRAY_BASE<T_NEW,typename T_ARRAY::template REBIND<T_NEW>::TYPE,ID> TYPE;};

template<class T,int d> class VECTOR;
template<class T_ARRAY,int length_new> struct REBIND_LENGTH;

template<class T_ARRAY> struct ARRAY_RESULT_TYPE{typedef typename T_ARRAY::RESULT_TYPE TYPE;};
template<class T_ARRAY> struct ARRAY_RESULT_TYPE<const T_ARRAY>{typedef typename T_ARRAY::CONST_RESULT_TYPE TYPE;};

template<class T,class ID,class SCALAR> struct REPLACE_FLOATING_POINT<ARRAY<T,ID>,SCALAR>{typedef ARRAY<typename REPLACE_FLOATING_POINT<T,SCALAR>::TYPE,ID> TYPE;};

}
#endif

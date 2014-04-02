//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADAPTIVE_FORWARD__
#define __ADAPTIVE_FORWARD__

#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>

namespace PhysBAM{
namespace ADAPTIVE_DETAIL{

enum ADAPTIVE_SIGN {ADAPTIVE_SIGN_NEGATIVE=-1,ADAPTIVE_SIGN_ZERO=0,ADAPTIVE_SIGN_POSITIVE=1,ADAPTIVE_SIGN_UNKNOWN=2};
struct NIL {};
template<class T>
struct IS_NIL {static const bool value=false;};
template<>
struct IS_NIL<NIL> {static const bool value=true;};

template<class T>
struct IS_ADAPTIVE {static const bool value=false;};

template<class T>
struct IS_NOT_ADAPTIVE {static const bool value=!IS_ADAPTIVE<T>::value;};

template<class DERIVED,class EXACT_TYPE_T> class ADAPTIVE_BASE;
template<class EXACT_TYPE_T> class ADAPTIVE_OBJECT;
template<class EXACT_TYPE_T> class ADAPTIVE_DYNAMIC_BASE;
template<class T_ADAPTIVE,class T_EXACT> class ADAPTIVE_OBJECT_BODY;

enum CACHED_TYPE {CACHED,NOT_CACHED};
template<class DERIVED,class T_WRAPPED,CACHED_TYPE cached_type=NOT_CACHED> class ADAPTIVE_WRAPPER_BASE;

template<class T,class EXACT_TYPE_T> class ADAPTIVE_ATOM;
template<class T_ADAPTIVE1> class ADAPTIVE_NEGATION;
template<class T_ADAPTIVE1,class T_ADAPTIVE2> class ADAPTIVE_SUM;
template<class T_ADAPTIVE1,class T_ADAPTIVE2> class ADAPTIVE_DIFFERENCE;
template<class T_ADAPTIVE1,class T_ADAPTIVE2> class ADAPTIVE_PRODUCT;
template<class T_ADAPTIVE1,class T_ADAPTIVE2> class ADAPTIVE_QUOTIENT;
template<class T_ADAPTIVE,int d> class ADAPTIVE_DETERMINANT;
template<class T_ADAPTIVE,int d> class ADAPTIVE_SIGNED_VOLUME;

struct GEOMETRIC_DEGENERACY{}; // exception

template<class T_EXACT,class T> struct ADAPTIVE_POLICY
{typedef typename IF<IS_ADAPTIVE<T>::value,T,ADAPTIVE_ATOM<T,T_EXACT> >::TYPE ADAPTIVE;};

}
using ADAPTIVE_DETAIL::ADAPTIVE_SUM;
using ADAPTIVE_DETAIL::ADAPTIVE_DIFFERENCE;
using ADAPTIVE_DETAIL::ADAPTIVE_PRODUCT;
using ADAPTIVE_DETAIL::GEOMETRIC_DEGENERACY;
}
#endif

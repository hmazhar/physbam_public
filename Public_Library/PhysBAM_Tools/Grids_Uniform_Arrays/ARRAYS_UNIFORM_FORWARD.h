//#####################################################################
// Copyright 2006-2009, Geoffrey Irving, Nipun Kwatra, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header ARRAYS_UNIFORM_FORWARD
//#####################################################################
#ifndef __ARRAYS_UNIFORM_FORWARD__
#define __ARRAYS_UNIFORM_FORWARD__

#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T,class T_ARRAY,class ID> class ARRAY_BASE;
template<class TV> class ARRAYS_ND_BASE;
template<class T,class ID> class ARRAY;
template<int d> class FACE_INDEX;
template<int d> class SIDED_FACE_INDEX;

class FLOOD_FILL_1D;
class FLOOD_FILL_2D;
class FLOOD_FILL_3D;

template<class T_ARRAY,class T_NEW> struct REBIND;
template<class T,class T_NEW,int d> struct REBIND<ARRAY<T,VECTOR<int,d> >,T_NEW>{typedef ARRAY<T_NEW,VECTOR<int,d> > TYPE;};
template<class T,class T_NEW,int d> struct REBIND<ARRAY<T,FACE_INDEX<d> >,T_NEW>{typedef ARRAY<T_NEW,FACE_INDEX<d> > TYPE;};
template<class T,class T_NEW,int d> struct REBIND<ARRAY<T,SIDED_FACE_INDEX<d> >,T_NEW>{typedef ARRAY<T_NEW,SIDED_FACE_INDEX<d> > TYPE;};
//template<class T,class T_NEW,int d> struct REBIND<ARRAYS_BASE<T,ARRAYS_ND_BASE,VECTOR<int,d> >,T_NEW>{typedef ARRAY_BASE<T_NEW,ARRAYS_ND_BASE,VECTOR<int,d> > TYPE;};
template<class T,class T_NEW,int d> struct REBIND<ARRAYS_ND_BASE<VECTOR<T,d> >,T_NEW>{typedef ARRAYS_ND_BASE<VECTOR<T_NEW,d> > TYPE;};

template<class T,int d> class VECTOR;
template<class T_ARRAY,int length_new> struct REBIND_LENGTH;
template<class T,int m,int d,int length_new> struct REBIND_LENGTH<ARRAY<VECTOR<T,d>,VECTOR<int,m> >,length_new>{typedef ARRAY<VECTOR<T,length_new>,VECTOR<int,m> > TYPE;};

}
#endif

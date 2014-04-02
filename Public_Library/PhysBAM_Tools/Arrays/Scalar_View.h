//#####################################################################
// Copyright 2008-2010, Geoffrey Irving, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function Scalar_View 
//#####################################################################
//
// Flattens any applicable type into a ARRAY_VIEW of the underlying scalar type.
//
//#####################################################################
#ifndef __Scalar_View__
#define __Scalar_View__

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UNIFORM_FORWARD.h>
#include <PhysBAM_Tools/Vectors/SCALAR_POLICY.h>
namespace PhysBAM{

// check if it's safe to use vector space operations on the result of Scalar_View<TV>
template<class TV> struct SCALAR_VIEW_IS_VECTOR_SPACE:public IS_SCALAR_VECTOR_SPACE<TV>{};
template<class TV> struct SCALAR_VIEW_IS_VECTOR_SPACE<ARRAY<TV> >:public IS_SCALAR_VECTOR_SPACE<TV>{};
template<class TV> struct SCALAR_VIEW_IS_VECTOR_SPACE<ARRAY_VIEW<TV> >:public IS_SCALAR_VECTOR_SPACE<TV>{};
//template<class TV,int d> struct SCALAR_VIEW_IS_VECTOR_SPACE<ARRAY<VECTOR<TV,d> ,VECTOR<int,1> > >:public IS_SCALAR_VECTOR_SPACE<TV>{};
//template<class TV,int d> struct SCALAR_VIEW_IS_VECTOR_SPACE<ARRAY<VECTOR<TV,d> ,VECTOR<int,2> > >:public IS_SCALAR_VECTOR_SPACE<TV>{};
//template<class TV,int d> struct SCALAR_VIEW_IS_VECTOR_SPACE<ARRAY<VECTOR<TV,d> ,VECTOR<int,3> > >:public IS_SCALAR_VECTOR_SPACE<TV>{};
template<class TV,int d> struct SCALAR_VIEW_IS_VECTOR_SPACE<ARRAY<TV,VECTOR<int,d> > >:public IS_SCALAR_VECTOR_SPACE<TV>{};
template<class T,int d> struct SCALAR_VIEW_IS_VECTOR_SPACE<ARRAYS_ND_BASE<VECTOR<T,d> > >:public IS_SCALAR_VECTOR_SPACE<T>{};

template<class TV> typename ENABLE_IF<IS_SCALAR_BLOCK<TV>::value,ARRAY_VIEW<typename SCALAR_POLICY<TV>::TYPE> >::TYPE
Scalar_View(TV& block)
{
    typedef typename SCALAR_POLICY<TV>::TYPE T;
    return ARRAY_VIEW<T>(sizeof(TV)/sizeof(T),reinterpret_cast<T*>(&block));
}

template<class TV,class T_ARRAY> ARRAY_VIEW<typename SCALAR_POLICY<TV>::TYPE>
Scalar_View(ARRAY_BASE<TV,T_ARRAY>& array)
{
    typedef typename SCALAR_POLICY<TV>::TYPE T;
    STATIC_ASSERT(IS_SCALAR_BLOCK<TV>::value);
    T_ARRAY& array_=array.Derived();
    return ARRAY_VIEW<T>(sizeof(TV)/sizeof(T)*array_.Size(),reinterpret_cast<T*>(array_.Get_Array_Pointer()));
}

template<class T,int d> ARRAY_VIEW<typename SCALAR_POLICY<T>::TYPE>
Scalar_View(ARRAYS_ND_BASE<VECTOR<T,d> >& array)
{
    return Scalar_View(array.array);
}

template<class T,int d> ARRAY_VIEW<T>
Scalar_View(ARRAY<T,FACE_INDEX<d> >& array)
{
    return ARRAY_VIEW<T>(array.buffer_size,array.base_pointer);
}
    
}
#endif

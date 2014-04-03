//#####################################################################
// Copyright 2002-2007, Ron Fedkiw, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FFTW
//#####################################################################
#ifndef __FFTW__
#define __FFTW__

#ifdef USE_FFTW

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRIDS_UNIFORM_ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/COMPLEX.h>
#include <fftw3.h>
namespace PhysBAM{

template<class T> struct FFTW_POLICY;
template<> struct FFTW_POLICY<float>{typedef fftwf_plan PLAN;};
template<> struct FFTW_POLICY<double>{typedef fftw_plan PLAN;};

template<class T,int d>
class FFTW
{
public:
    typedef typename ARRAYS_ND<T,d>::TYPE T_ARRAYS_T;
    typedef typename T_ARRAYS_T::template REBIND<COMPLEX<T> >::TYPE T_ARRAYS_COMPLEX;
    typedef typename FFTW_POLICY<T>::PLAN T_PLAN;

    static const unsigned int plan_flags=FFTW_ESTIMATE | FFTW_UNALIGNED;

//#####################################################################
    static void Destroy_Plan(const T_PLAN plan);
    static T_PLAN Plan_R2C(const VECTOR<int,d>& counts,const T_ARRAYS_T& u,const T_ARRAYS_COMPLEX& u_hat);
    static T_PLAN Plan_C2R(const VECTOR<int,d>& counts,const T_ARRAYS_COMPLEX& u_hat,const T_ARRAYS_T& u);
    static void Execute_R2C(const T_PLAN plan,const T_ARRAYS_T& u,T_ARRAYS_COMPLEX& u_hat);
    static void Execute_C2R(const T_PLAN plan,T_ARRAYS_COMPLEX& u_hat,T_ARRAYS_T& u); // destroys input
//#####################################################################
};

#ifndef PLATFORM_WINDOWS
#define SPECIALIZE_FFTW_d(T_FFTW,T,d) \
    template<> inline void FFTW<T,d>::Destroy_Plan(const T_PLAN plan) \
    {if(plan) T_FFTW##_destroy_plan(plan);} \
    \
    template<> inline FFTW_POLICY<T>::PLAN FFTW<T,d>::Plan_R2C(const VECTOR<int,d>& counts,const T_ARRAYS_T& u,const T_ARRAYS_COMPLEX& u_hat) \
    {return T_FFTW##_plan_dft_r2c(d,&counts[1],const_cast<T_ARRAYS_T&>(u).array.Get_Array_Pointer(),(T(*)[2])const_cast<T_ARRAYS_COMPLEX&>(u_hat).array.Get_Array_Pointer(),plan_flags);} \
    \
    template<> inline FFTW_POLICY<T>::PLAN FFTW<T,d>::Plan_C2R(const VECTOR<int,d>& counts,const T_ARRAYS_COMPLEX& u_hat,const T_ARRAYS_T& u) \
    {return T_FFTW##_plan_dft_c2r(d,&counts[1],(T(*)[2])const_cast<T_ARRAYS_COMPLEX&>(u_hat).array.Get_Array_Pointer(),const_cast<T_ARRAYS_T&>(u).array.Get_Array_Pointer(),plan_flags);} \
    \
    template<> inline void FFTW<T,d>::Execute_R2C(const T_PLAN plan,const T_ARRAYS_T& u,T_ARRAYS_COMPLEX& u_hat) \
    {return T_FFTW##_execute_dft_r2c(plan,const_cast<T_ARRAYS_T&>(u).array.Get_Array_Pointer(),(T(*)[2])u_hat.array.Get_Array_Pointer());} \
    \
    template<> inline void FFTW<T,d>::Execute_C2R(const T_PLAN plan,T_ARRAYS_COMPLEX& u_hat,T_ARRAYS_T& u) \
    {return T_FFTW##_execute_dft_c2r(plan,(T(*)[2])u_hat.array.Get_Array_Pointer(),u.array.Get_Array_Pointer());}
#else

// MSVC does not understand "float(*)[2]" as a legitimate type, so we must express it in terms
//   of derived classes (in this case: "fftwf_complex *")
#define SPECIALIZE_FFTW_d(T_FFTW,T,d) \
    template<> inline void FFTW<T,d>::Destroy_Plan(const T_PLAN plan) \
    {if(plan) T_FFTW##_destroy_plan(plan);} \
    \
    template<> inline FFTW_POLICY<T>::PLAN FFTW<T,d>::Plan_R2C(const VECTOR<int,d>& counts,const T_ARRAYS_T& u,const T_ARRAYS_COMPLEX& u_hat) \
    {return T_FFTW##_plan_dft_r2c(d,&counts[1],const_cast<T_ARRAYS_T&>(u).array.Get_Array_Pointer(),(T_FFTW##_complex *)const_cast<T_ARRAYS_COMPLEX&>(u_hat).array.Get_Array_Pointer(),plan_flags);} \
    \
    template<> inline FFTW_POLICY<T>::PLAN FFTW<T,d>::Plan_C2R(const VECTOR<int,d>& counts,const T_ARRAYS_COMPLEX& u_hat,const T_ARRAYS_T& u) \
    {return T_FFTW##_plan_dft_c2r(d,&counts[1],(T_FFTW##_complex *)const_cast<T_ARRAYS_COMPLEX&>(u_hat).array.Get_Array_Pointer(),const_cast<T_ARRAYS_T&>(u).array.Get_Array_Pointer(),plan_flags);} \
    \
    template<> inline void FFTW<T,d>::Execute_R2C(const T_PLAN plan,const T_ARRAYS_T& u,T_ARRAYS_COMPLEX& u_hat) \
    {return T_FFTW##_execute_dft_r2c(plan,const_cast<T_ARRAYS_T&>(u).array.Get_Array_Pointer(),(T_FFTW##_complex *)u_hat.array.Get_Array_Pointer());} \
    \
    template<> inline void FFTW<T,d>::Execute_C2R(const T_PLAN plan,T_ARRAYS_COMPLEX& u_hat,T_ARRAYS_T& u) \
        {return T_FFTW##_execute_dft_c2r(plan,(T_FFTW##_complex *)u_hat.array.Get_Array_Pointer(),u.array.Get_Array_Pointer());}

#endif


#define SPECIALIZE_FFTW(T_FFTW,T) \
    SPECIALIZE_FFTW_d(T_FFTW,T,1) \
    SPECIALIZE_FFTW_d(T_FFTW,T,2) \
    SPECIALIZE_FFTW_d(T_FFTW,T,3)

SPECIALIZE_FFTW(fftwf,float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
SPECIALIZE_FFTW(fftw,double)
#endif

#undef SPECIALIZE_FFTW_d
#undef SPECIALIZE_FFTW

}
#endif
#endif

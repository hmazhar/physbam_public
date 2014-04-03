//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Neil Molino, Igor Neverov, Duc Nguyen, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RANDOM_NUMBERS
//#####################################################################
#ifndef __RANDOM_NUMBERS__
#define __RANDOM_NUMBERS__

#include <PhysBAM_Tools/Random_Numbers/MT19937.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <time.h>
namespace PhysBAM{

template<class TV> class RANGE;
template<class TV> class ROTATION;
template<class TV> class FRAME;
template<class TV> class TWIST;
template<class T,int d> class VECTOR;
template<class T,int d> class DIAGONAL_MATRIX;
template<class T,int d> class SYMMETRIC_MATRIX;
template<class T,int d> class UPPER_TRIANGULAR_MATRIX;
template<class T,class T_VECTOR> class VECTOR_BASE;
template<class T,class T_MATRIX> class MATRIX_BASE;
template<class T,class T_ARRAY,class ID> class ARRAY_BASE;
class TYPED_ISTREAM;
class TYPED_OSTREAM;

template<class T,class GENERATOR=MT19937<T> >
class RANDOM_NUMBERS:public NONCOPYABLE
{
public:
    int gaussian_iset; // Used to force Get_Gaussian to reset
    T gset; // used internally by Get_Gaussian
    GENERATOR random_number_generator;
public:

    T Get_Number()
    {return random_number_generator();} // in [0,1)

    explicit RANDOM_NUMBERS(const unsigned int seed=time(0));
    virtual ~RANDOM_NUMBERS();

    template<class T2,class T_ARRAY,class ID> void Fill_Uniform(ARRAY_BASE<T2,T_ARRAY,ID>& array,const T a,const T b)
    {T_ARRAY& derived=array.Derived();for(ID i(1);i<=derived.Size();i++) Fill_Uniform(derived(i),a,b);}

//#####################################################################
    void Set_Seed(const unsigned int seed_input=time(0));
    int Get_Uniform_Integer(const int a,const int b);
    T Get_Uniform_Number(const T a,const T b);
    template<int d> VECTOR<T,d> Get_Uniform_Vector(const VECTOR<T,d>& v0,const VECTOR<T,d>& v1);
    template<int d> VECTOR<T,d> Get_Uniform_Vector(const T a,const T b);
    template<class TV> TV Get_Uniform_Vector(const RANGE<TV>& box);
    void Fill_Uniform(T& x,const T a,const T b);
    template<class T_VECTOR> void Fill_Uniform(VECTOR_BASE<T,T_VECTOR>& v,const T a,const T b);
    template<class T_MATRIX> void Fill_Uniform(MATRIX_BASE<T,T_MATRIX>& m,const T a,const T b);
    template<int d> void Fill_Uniform(DIAGONAL_MATRIX<T,d>& m,const T a,const T b);
    template<int d> void Fill_Uniform(SYMMETRIC_MATRIX<T,d>& m,const T a,const T b);
    template<int d> void Fill_Uniform(UPPER_TRIANGULAR_MATRIX<T,d>& m,const T a,const T b);
    template<class TV> void Fill_Uniform(TWIST<TV>& m,const T a,const T b);
    T Get_Gaussian();
    template<class TV> TV Get_Vector_In_Unit_Sphere();
    template<class TV> TV Get_Direction();
    template<class TV> ROTATION<TV> Get_Rotation();
    template<class TV> FRAME<TV> Get_Frame(const TV& v0,const TV& v1);
    template<class TV> TWIST<TV> Get_Twist(const T& a);
//#####################################################################
};
}
#endif

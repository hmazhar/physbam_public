//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BANDED_SYMMETRIC_MATRIX
//#####################################################################
#ifndef __BANDED_SYMMETRIC_MATRIX__
#define __BANDED_SYMMETRIC_MATRIX__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T,int bandwidth>
class BANDED_SYMMETRIC_MATRIX
{
    STATIC_ASSERT(bandwidth&1);
public:
    ARRAY<VECTOR<T,1+bandwidth/2> > A; // A(i)[j] = a_i,i+j-1

    BANDED_SYMMETRIC_MATRIX(const int m=0)
        :A(m)
    {}

    int Size() const
    {return A.Size();}

    void Resize(const int m_new)
    {A.Resize(m_new);}

    // eigenvalue is within tolerance of a real eigenvalue if return true
    bool Power_Iterate(VECTOR_ND<T>& x,T& eigenvalue,const T tolerance,const int max_iterations) const
    {return Power_Iterate_Shifted(x,0,eigenvalue,tolerance,max_iterations);}

    // see Golub and Van Loan, 5.1.8, p. 216
    static void Givens(const T a,const T b,T& c,T& s)
    {if(b==0){c=1;s=0;}
    else if(abs(b)>abs(a)){T t=-a/b;s=1/sqrt(1+sqr(t));c=s*t;}
    else{T t=-b/a;c=1/sqrt(1+sqr(t));s=c*t;}}

//#####################################################################
    void Multiply(const VECTOR_ND<T>& x,VECTOR_ND<T>& result) const;
    bool Power_Iterate_Shifted(VECTOR_ND<T>& x,const T shift,T& eigenvalue,const T tolerance,const int max_iterations) const;
    bool Eigenvalue_Range(VECTOR_ND<T>& x,INTERVAL<T>& eigenvalue_range,const T tolerance,const int max_iterations) const;
    void Diagonalize();
    void Eigenvalues(ARRAY<T>& eigenvalues) const;
    void Print_Spectral_Information() const;
private:
    void Implicit_QR_Step(const int start,const T shift);
//#####################################################################
};
}
#endif

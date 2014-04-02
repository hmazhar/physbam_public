//#####################################################################
// Copyright 2007, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRANSPOSE_MATRIX
//#####################################################################
#ifndef __TRANSPOSE_MATRIX__
#define __TRANSPOSE_MATRIX__

#include <PhysBAM_Tools/Matrices/MATRIX_ARITHMETIC_POLICY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_BASE.h>
namespace PhysBAM{

template<class T_MATRIX>
class TRANSPOSE_MATRIX:public MATRIX_BASE<typename T_MATRIX::SCALAR,typename TRANSPOSE<T_MATRIX>::TYPE>
{
    typedef typename T_MATRIX::SCALAR T;
    typedef typename TRANSPOSE<T_MATRIX>::TYPE DERIVED;
public:
    typedef T SCALAR;
    enum WORKAROUND1 {m=T_MATRIX::n,n=T_MATRIX::m};

    T_MATRIX transpose;

protected:
    TRANSPOSE_MATRIX()
    {
        STATIC_ASSERT((IS_BASE_OF<TRANSPOSE_MATRIX,DERIVED>::value));
    }

    TRANSPOSE_MATRIX(INITIAL_SIZE mm,INITIAL_SIZE nn)
        :transpose(nn,mm)
    {}

    TRANSPOSE_MATRIX(const TRANSPOSE_MATRIX& M)
        :transpose(M.transpose)
    {}

    template<class T_MATRIX2>
    TRANSPOSE_MATRIX(const TRANSPOSE_MATRIX<T_MATRIX2>& M)
        :transpose(M.transpose)
    {}
public:

    DERIVED& operator=(const TRANSPOSE_MATRIX& M)
    {
        transpose=M.transpose;return Derived();
    }

    DERIVED& Derived()
    {return static_cast<DERIVED&>(*this);}

    int Rows() const
    {return transpose.Columns();}

    int Columns() const
    {return transpose.Rows();}

    T& operator()(const int i,const int j)
    {return transpose(j,i);}

    const T& operator()(const int i,const int j) const
    {return transpose(j,i);}

    bool Valid_Index(const int i,const int j) const
    {return transpose.Valid_Index(j,i);}

    bool operator==(const DERIVED& A) const
    {return transpose==A.transpose;}

    bool operator!=(const DERIVED& A) const
    {return transpose!=A.transpose;}

    DERIVED operator-() const
    {DERIVED r;r.transpose=-transpose;return r;}

    DERIVED operator*(const T a) const
    {DERIVED r;r.transpose=transpose*a;return r;}

    DERIVED& operator*=(const T a)
    {transpose*=a;return Derived();}

    DERIVED operator/(const T a) const
    {DERIVED r;r.transpose=transpose/a;return r;}

    DERIVED& operator/=(const T a)
    {transpose/=a;return Derived();}

    DERIVED operator+(const DERIVED& A) const
    {DERIVED r;r.transpose=transpose+A.transpose;return r;}

    DERIVED operator+=(const DERIVED& A)
    {transpose+=A.transpose;return Derived();}

    DERIVED operator-(const DERIVED& A) const
    {DERIVED r;r.transpose=transpose-A.transpose;return r;}

    DERIVED operator-=(const DERIVED& A)
    {transpose-=A.transpose;return Derived();}

    template<class TM>
    typename PRODUCT<DERIVED,TM>::TYPE operator*(const TM& A) const
    {return this->transpose.Transpose_Times(A);}

    static DERIVED Transposed(const T_MATRIX& A)
    {DERIVED transpose;transpose.transpose=A;return transpose;} 

    const T_MATRIX& Transposed() const
    {return transpose;}

    DERIVED Cofactor_Matrix() const
    {return Transposed(transpose.Cofactor_Matrix());}

    static DERIVED Outer_Product(const VECTOR<T,m>& u,const VECTOR<T,n>& v)
    {DERIVED result;result.transpose=T_MATRIX::Outer_Product(v,u);return result;}

    T Max_Abs() const
    {return transpose.Max_Abs();}

    T Frobenius_Norm_Squared() const
    {return transpose.Frobenius_Norm_Squared();}

    template<class TM> typename PRODUCT<T_MATRIX,TM>::TYPE Transpose_Times(const TM& A) const
    {return transpose*A;}

    template<class TM> typename PRODUCT<DERIVED,typename TRANSPOSE<TM>::TYPE>::TYPE Times_Transpose(const TM& A) const
    {return (A*transpose).Transposed();}

//#####################################################################
};

template<class TM> typename TRANSPOSE<TM>::TYPE
operator*(const typename TM::SCALAR& a,const TRANSPOSE_MATRIX<TM>& A)
{
    typename TRANSPOSE<TM>::TYPE result;result.transpose=A.transpose*a;return result;
}

template<class TM1,class TM2> typename TRANSPOSE<typename PRODUCT<TM2,typename TRANSPOSE<TM1>::TYPE>::TYPE>::TYPE
operator*(const TM1& A,const TRANSPOSE_MATRIX<TM2>& B)
{
    return A.Times_Transpose(B.transpose);
}

}
#endif

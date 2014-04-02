//#####################################################################
// Copyright 2007, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MATRIX_0X0__
#define __MATRIX_0X0__

#include <PhysBAM_Tools/Matrices/MATRIX_BASE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_0D.h>
#include <cfloat>
namespace PhysBAM{

template<class T>
class MATRIX<T,0>:public MATRIX_BASE<T,MATRIX<T,0> >
{
public:
    typedef T SCALAR;
    enum WORKAROUND1 {m=0,n=0};

    explicit MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(),INITIAL_SIZE nn=INITIAL_SIZE())
    {assert(mm==INITIAL_SIZE() && nn==INITIAL_SIZE());}

    MATRIX(const MATRIX& matrix)
    {}

    template<class T2>
    explicit MATRIX(const MATRIX<T2,0>& matrix)
    {}

    template<class T_MATRIX>
    explicit MATRIX(const MATRIX_BASE<T,T_MATRIX>& A)
    {
        assert(A.Rows()==0 && A.Columns()==0);
    }

    MATRIX& operator=(const MATRIX& matrix)
    {return *this;}

    int Rows() const
    {return 0;}

    int Columns() const
    {return 0;}

    T& operator()(const int i,const int j=0)
    {PHYSBAM_FATAL_ERROR();}

    const T& operator()(const int i,const int j=0) const
    {PHYSBAM_FATAL_ERROR();}

    VECTOR<T,0>& Column(const int i)
    {PHYSBAM_FATAL_ERROR();}

    const VECTOR<T,0>& Column(const int i) const
    {PHYSBAM_FATAL_ERROR();}

    bool Valid_Index(const int i,const int j) const
    {return false;}

    bool operator==(const MATRIX& A) const
    {return true;}

    bool operator!=(const MATRIX& A) const
    {return false;}

    MATRIX Transposed() const
    {return *this;}

    void Transpose()
    {}

    MATRIX Inverse() const
    {return *this;}

    VECTOR<T,0> Solve_Linear_System(const VECTOR<T,0>& b) const
    {return b;}

    VECTOR<T,0> Robust_Solve_Linear_System(const VECTOR<T,0>& b) const
    {return b;}

    MATRIX Cofactor_Matrix() const
    {return *this;}

    MATRIX Normal_Equations_Matrix() const
    {return *this;}

    MATRIX operator-() const
    {return *this;}

    MATRIX operator+(const T a) const
    {return *this;}

    MATRIX operator+(const MATRIX& A) const
    {return *this;}

    MATRIX operator*(const T a) const
    {return *this;}

    MATRIX operator/(const T a) const
    {return *this;}

    VECTOR<T,0> operator*(const VECTOR<T,0>& v) const
    {return v;}

    MATRIX operator*(const MATRIX& A) const
    {return *this;}

    MATRIX& operator+=(const T a)
    {return *this;}

    MATRIX& operator+=(const MATRIX& A)
    {return *this;}

    MATRIX operator-(const T a) const
    {return *this;}

    MATRIX operator-(const MATRIX& A) const
    {return *this;}

    MATRIX& operator-=(const MATRIX& A)
    {return *this;}

    MATRIX& operator-=(const T& a)
    {return *this;}

    MATRIX& operator*=(const MATRIX& A)
    {return *this;}

    MATRIX& operator*=(const T a)
    {return *this;}

    MATRIX& operator/=(const T a)
    {return *this;}

    MATRIX_MXN<T> operator*(const MATRIX_MXN<T>& A) const
    {assert(A.m==0);return A;}

    static MATRIX Outer_Product(const VECTOR<T,0>& v)
    {return MATRIX();}

    static MATRIX Outer_Product(const VECTOR<T,0>& u,const VECTOR<T,0>& v)
    {return MATRIX();}

    T Max() const
    {return -FLT_MAX;}

    T Min() const
    {return FLT_MAX;}

    T Frobenius_Norm() const
    {return 0;}

    T Frobenius_Norm_Squared() const
    {return 0;}

    T Inner_Product(const VECTOR<T,0>& u,const VECTOR<T,0>& v) const
    {return 0;}

    bool Positive_Definite() const
    {return true;}

    bool Positive_Semidefinite() const
    {return true;}

    MATRIX Positive_Definite_Part() const
    {return *this;}

    MATRIX Diagonal_Part() const
    {return *this;}

    MATRIX Sqrt() const
    {return *this;}

    MATRIX Clamp_Min(const T a) const
    {return *this;}

    MATRIX Clamp_Max(const T a) const
    {return *this;}

    MATRIX Abs() const
    {return *this;}

    MATRIX Sign() const
    {return *this;}

    static MATRIX Identity_Matrix()
    {return MATRIX();}

    T Trace() const
    {return 0;}

    T Determinant() const
    {return 1;}

    MATRIX Symmetric_Part() const
    {return MATRIX();}

    void Fast_Solve_Eigenproblem(MATRIX&,MATRIX&) const
    {}

    void Fast_Singular_Value_Decomposition(MATRIX&,MATRIX&,MATRIX&) const
    {}

    MATRIX<T,0,1> Times_Cross_Product_Matrix(const VECTOR<T,1>& v)
    {return MATRIX<T,0,1>();}

    MATRIX<T,1,0> Cross_Product_Matrix_Transpose_Times(const VECTOR<T,1>& v)
    {return MATRIX<T,1,0>();}

//#####################################################################
};

template<class T>
inline MATRIX<T,0> operator*(const T a,const MATRIX<T,0>& A)
{return A;}

template<class T>
inline MATRIX<T,0> operator+(const T a,const MATRIX<T,0>& A)
{return A;}

template<class T>
inline MATRIX<T,0> operator-(const T a,const MATRIX<T,0>& A)
{return A;}

template<class T>
inline MATRIX<T,0> log(const MATRIX<T,0>& A)
{return MATRIX<T,0>();}

template<class T>
inline MATRIX<T,0> exp(const MATRIX<T,0>& A)
{return MATRIX<T,0>();}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
template<class T>
inline std::istream& operator>>(std::istream& input,MATRIX<T,0>& A)
{FILE_UTILITIES::Ignore(input,'[');FILE_UTILITIES::Ignore(input,']');return input;}
#endif
}
#endif

//#####################################################################
// Copyright 2007, Nipun Kwatra, Craig Schroeder, Andrew Selle, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MATRIX_1X1__
#define __MATRIX_1X1__

#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Matrices/MATRIX_BASE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
namespace PhysBAM{

using std::log;
using std::exp;
using std::sqrt;

template<class T>
class MATRIX<T,1>:public MATRIX_BASE<T,MATRIX<T,1> >
{
    struct UNUSABLE{};
public:
    typedef T SCALAR;
    enum WORKAROUND1 {m=1,n=1};
    typedef MATRIX_BASE<T,MATRIX<T,1> > BASE;using BASE::operator*;

    T x11;

    explicit MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(1),INITIAL_SIZE nn=INITIAL_SIZE(1))
        :x11(0)
    {
        assert(mm==INITIAL_SIZE(1) && nn==INITIAL_SIZE(1));
    }

    MATRIX(const MATRIX& matrix)
        :x11(matrix.x11)
    {}

    template<class T2>
    explicit MATRIX(const MATRIX<T2,1>& matrix)
        :x11(matrix.x11)
    {}

    template<class T_MATRIX>
    explicit MATRIX(const MATRIX_BASE<T,T_MATRIX>& A)
        :x11(A(1,1))
    {
        assert(A.Rows()==1 && A.Columns()==1);
    }

    explicit MATRIX(const T x11)
        :x11(x11)
    {}

    explicit MATRIX(const VECTOR<T,1>& v)
        :x11(v.x)
    {}

    MATRIX& operator=(const MATRIX& matrix)
    {x11=matrix.x11;return *this;}

    int Rows() const
    {return 1;}

    int Columns() const
    {return 1;}

    T& operator()(const int i,const int j=1)
    {assert(i==1 && j==1);return x11;}

    const T& operator()(const int i,const int j=1) const
    {assert(i==1 && j==1);return x11;}

    bool Valid_Index(const int i,const int j) const
    {return i==1 && j==1;}

    VECTOR<T,1>& Column(const int j)
    {assert(j==1);return *(VECTOR<T,1>*)this;}

    const VECTOR<T,1>& Column(const int j) const
    {assert(j==1);return *(const VECTOR<T,1>*)this;}

    bool operator==(const MATRIX& A) const
    {return x11==A.x11;}

    bool operator!=(const MATRIX& A) const
    {return !(*this==A);}

    VECTOR<T,1> Column_Sum() const
    {return Column(1);}

    VECTOR<T,1> Column_Magnitudes() const
    {return VECTOR<T,1>(Column(1).Magnitude());}

    MATRIX Inverse() const
    {assert(x11);return MATRIX(1/x11);}

    VECTOR<T,1> Solve_Linear_System(const VECTOR<T,1>& b) const
    {return VECTOR<T,1>(b.x/x11);}

    VECTOR<T,1> Robust_Solve_Linear_System(const VECTOR<T,1>& b) const
    {return VECTOR<T,1>(Robust_Divide(b.x,x11));}

    MATRIX Cofactor_Matrix() const
    {return MATRIX(1);}

    MATRIX Normal_Equations_Matrix() const // 1 mult
    {return MATRIX(sqr(x11));}

    MATRIX operator-() const
    {return MATRIX(-x11);}

    MATRIX operator+(const T a) const
    {return MATRIX(x11+a);}

    MATRIX operator+(const MATRIX& A) const
    {return MATRIX(x11+A.x11);}

    MATRIX operator-(const T a) const
    {return MATRIX(x11-a);}

    MATRIX operator-(const MATRIX& A) const
    {return MATRIX(x11-A.x11);}

    MATRIX operator*(const MATRIX& A) const
    {return MATRIX(x11*A.x11);}

    MATRIX_MXN<T> operator*(const MATRIX_MXN<T>& A) const
    {assert(A.m==1);return A*x11;}

    MATRIX operator*(const T a) const
    {return MATRIX(a*x11);}

    MATRIX operator/(const T a) const
    {return MATRIX(x11/a);}

    VECTOR<T,1> operator*(const VECTOR<T,1>& v) const
    {return VECTOR<T,1>(x11*v.x);}

    MATRIX& operator+=(const MATRIX& A)
    {x11+=A.x11;return *this;}

    MATRIX& operator-=(const MATRIX& A)
    {x11-=A.x11;return *this;}

    MATRIX& operator+=(const T& a)
    {x11+=a;return *this;}

    MATRIX& operator-=(const T& a)
    {x11-=a;return *this;}

    MATRIX& operator*=(const T a)
    {x11*=a;return *this;}

    MATRIX& operator*=(const MATRIX& A)
    {x11*=A.x11;return *this;}

    MATRIX& operator/=(const T a)
    {x11/=a;return *this;}

    VECTOR<T,1> Transpose_Times(const VECTOR<T,1>& v) const
    {return VECTOR<T,1>(x11*v.x);}

    template<class T_MATRIX>
    T_MATRIX Transpose_Times(const MATRIX_BASE<T,T_MATRIX>& A) const
    {return *this*A.Derived();}

    template<class T_MATRIX>
    typename TRANSPOSE<T_MATRIX>::TYPE Times_Transpose(const MATRIX_BASE<T,T_MATRIX>& A) const
    {return *this*A.Derived().Transposed();}

    static MATRIX Outer_Product(const VECTOR<T,1>& u)
    {return MATRIX(u.x*u.x);}

    static MATRIX Outer_Product(const VECTOR<T,1>& u,const VECTOR<T,1>& v)
    {return MATRIX(u.x*v.x);}

    void Transpose()
    {}

    MATRIX Transposed() const
    {return *this;}

    T Trace() const
    {return x11;}

    T Determinant() const
    {return x11;}

    MATRIX<T,1> Fast_Eigenvalues() const
    {return *this;}

    T Max() const
    {return x11;}

    T Min() const
    {return x11;}

    T Frobenius_Norm() const
    {return abs(x11);}

    T Frobenius_Norm_Squared() const
    {return sqr(x11);}

    SYMMETRIC_MATRIX<T,2> Conjugate_With_Cross_Product_Matrix(const VECTOR<T,2>& v) const
    {return x11*SYMMETRIC_MATRIX<T,2>(sqr(v.y),-v.x*v.y,sqr(v.x));}

    T Inner_Product(const VECTOR<T,1>& u,const VECTOR<T,1>& v) const
    {return x11*u.x*v.x;}

    T Inverse_Inner_Product(const VECTOR<T,1>& u,const VECTOR<T,1>& v) const
    {return u.x*v.x/x11;}

    MATRIX<T,1,2> Times_Cross_Product_Matrix(const VECTOR<T,2>& v) const
    {return x11*MATRIX<T,1,2>::Cross_Product_Matrix(v);}

    MATRIX<T,0,1> Cross_Product_Matrix_Times(const VECTOR<T,1>& v)
    {return MATRIX<T,0,1>();}

    bool Positive_Definite() const
    {return x11>0;}

    bool Positive_Semidefinite() const
    {return x11>=0;}

    MATRIX Positive_Definite_Part() const
    {return Clamp_Min(0);}

    MATRIX Diagonal_Part() const
    {return *this;}

    MATRIX Sqrt() const
    {return MATRIX(sqrt(x11));}

    MATRIX Clamp_Min(const T a) const
    {return MATRIX(clamp_min(x11,a));}

    MATRIX Clamp_Max(const T a) const
    {return MATRIX(clamp_max(x11,a));}

    MATRIX Abs() const
    {return MATRIX(abs(x11));}

    MATRIX Sign() const
    {return MATRIX(sign(x11));}

    static MATRIX Identity_Matrix()
    {return MATRIX((T)1);}

    MATRIX Symmetric_Part() const
    {return *this;}

    VECTOR<T,1> To_Vector() const
    {return VECTOR<T,1>(x11);}

    static MATRIX Conjugate(const MATRIX& A,const MATRIX& B)
    {return A*B*A;}

    void Solve_Eigenproblem(MATRIX& D,MATRIX& V) const
    {Fast_Solve_Eigenproblem(D,V);}

    void Fast_Solve_Eigenproblem(MATRIX& D,MATRIX& V) const
    {V.x11=1;D=*this;}
    
    void Fast_Singular_Value_Decomposition(MATRIX& U,MATRIX& D,MATRIX& V) const
    {U.x11=V.x11=1;D=*this;}

//#####################################################################
};

template<class T>
inline MATRIX<T,1> operator*(const T a,const MATRIX<T,1>& A)
{return A*a;}

template<class T>
inline MATRIX<T,1> operator+(const T a,const MATRIX<T,1>& A)
{return A+a;}

template<class T>
inline MATRIX<T,1> operator-(const T a,const MATRIX<T,1>& A)
{return MATRIX<T,1>(a-A.x11);}

template<class T>
inline MATRIX<T,1> clamp(const MATRIX<T,1>& x,const MATRIX<T,1>& xmin,const MATRIX<T,1>& xmax)
{return MATRIX<T,1>(clamp(x.x11,xmin.x11,xmax.x11));}

template<class T>
inline MATRIX<T,1> clamp_min(const MATRIX<T,1>& x,const MATRIX<T,1>& xmin)
{return MATRIX<T,1>(clamp_min(x.x11,xmin.x11));}

template<class T>
inline MATRIX<T,1> clamp_max(const MATRIX<T,1>& x,const MATRIX<T,1>& xmax)
{return MATRIX<T,1>(clamp_max(x.x11,xmax.x11));}

template<class T>
inline MATRIX<T,1> log(const MATRIX<T,1>& A)
{return MATRIX<T,1>(log(A.x11));}

template<class T>
inline MATRIX<T,1> exp(const MATRIX<T,1>& A)
{return MATRIX<T,1>(exp(A.x11));}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
template<class T>
inline std::istream& operator>>(std::istream& input,MATRIX<T,1>& A)
{FILE_UTILITIES::Ignore(input,'[');input>>A.x11;FILE_UTILITIES::Ignore(input,']');return input;}
#endif
}
#endif

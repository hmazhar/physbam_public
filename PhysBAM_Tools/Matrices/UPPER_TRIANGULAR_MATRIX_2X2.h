//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UPPER_TRIANGULAR_MATRIX_2X2
//#####################################################################
#ifndef __UPPER_TRIANGULAR_MATRIX_2X2__
#define __UPPER_TRIANGULAR_MATRIX_2X2__

#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <cassert>
#include <cmath>
#include <ostream>
namespace PhysBAM{

using ::std::sqrt;

template<class T> struct IS_SCALAR_BLOCK<UPPER_TRIANGULAR_MATRIX<T,2> >:public IS_SCALAR_BLOCK<T>{};
template<class T> struct IS_SCALAR_VECTOR_SPACE<UPPER_TRIANGULAR_MATRIX<T,2> >:public IS_SCALAR_VECTOR_SPACE<T>{};
template<class T,class RW> struct IS_BINARY_IO_SAFE<UPPER_TRIANGULAR_MATRIX<T,2>,RW>:public IS_BINARY_IO_SAFE<T,RW>{};

template<class T>
class UPPER_TRIANGULAR_MATRIX<T,2>
{
public:
    typedef T SCALAR;
    enum WORKAROUND1 {m=2,n=2};

    T x11,x12,x22;

    UPPER_TRIANGULAR_MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(2),INITIAL_SIZE nn=INITIAL_SIZE(2))
        :x11(0),x12(0),x22(0)
    {
        STATIC_ASSERT(sizeof(UPPER_TRIANGULAR_MATRIX)==3*sizeof(T));assert(mm==INITIAL_SIZE(2) && nn==INITIAL_SIZE(2));
    }

    template<class T2> explicit
    UPPER_TRIANGULAR_MATRIX(const UPPER_TRIANGULAR_MATRIX<T2,2>& matrix_input)
        :x11(matrix_input.x11),x12(matrix_input.x12),x22(matrix_input.x22)
    {}

    UPPER_TRIANGULAR_MATRIX(const T x11_input,const T x12_input,const T x22_input)
        :x11(x11_input),x12(x12_input),x22(x22_input)
    {}

    int Rows() const
    {return 2;}

    int Columns() const
    {return 2;}

    T& operator()(const int i,const int j)
    {assert(1<=i && i<=j && j<=2);return ((T*)this)[((j*(j-1))>>1)+i-1];}

    const T& operator()(const int i,const int j) const
    {assert(1<=i && i<=j && j<=2);return ((const T*)this)[((j*(j-1))>>1)+i-1];}

    bool Valid_Index(const int i,const int j) const
    {return 1<=i && i<=j && j<=2;}

    bool operator==(const UPPER_TRIANGULAR_MATRIX& A) const
    {return x11==A.x11 && x12==A.x12 && x22==A.x22;}

    bool operator!=(const UPPER_TRIANGULAR_MATRIX& A) const
    {return !(*this==A);}

    UPPER_TRIANGULAR_MATRIX operator-() const
    {return UPPER_TRIANGULAR_MATRIX(-x11,-x12,-x22);}

    UPPER_TRIANGULAR_MATRIX& operator+=(const UPPER_TRIANGULAR_MATRIX& A)
    {x11+=A.x11;x12+=A.x12;x22+=A.x22;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator+=(const T& a)
    {x11+=a;x22+=a;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator-=(const UPPER_TRIANGULAR_MATRIX& A)
    {x11-=A.x11;x12-=A.x12;x22-=A.x22;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator-=(const T& a)
    {x11-=a;x22-=a;return *this;}

    UPPER_TRIANGULAR_MATRIX operator+(const UPPER_TRIANGULAR_MATRIX& A) const
    {return UPPER_TRIANGULAR_MATRIX(x11+A.x11,x12+A.x12,x22+A.x22);}

    UPPER_TRIANGULAR_MATRIX operator+(const DIAGONAL_MATRIX<T,2>& A) const
    {return UPPER_TRIANGULAR_MATRIX(x11+A.x11,x12,x22+A.x22);}

    UPPER_TRIANGULAR_MATRIX operator+(const T a) const
    {return UPPER_TRIANGULAR_MATRIX(x11+a,x12,x22+a);}

    UPPER_TRIANGULAR_MATRIX operator-(const UPPER_TRIANGULAR_MATRIX& A) const
    {return UPPER_TRIANGULAR_MATRIX(x11-A.x11,x12-A.x12,x22-A.x22);}

    UPPER_TRIANGULAR_MATRIX operator-(const DIAGONAL_MATRIX<T,2>& A) const
    {return UPPER_TRIANGULAR_MATRIX(x11-A.x11,x12,x22-A.x22);}

    UPPER_TRIANGULAR_MATRIX operator-(const T a) const
    {return UPPER_TRIANGULAR_MATRIX(x11-a,x12,x22-a);}

    UPPER_TRIANGULAR_MATRIX& operator*=(const UPPER_TRIANGULAR_MATRIX& A)
    {return *this=*this*A;}

    UPPER_TRIANGULAR_MATRIX& operator*=(const T a)
    {x11*=a;x12*=a;x22*=a;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator/=(const T a)
    {assert(a!=0);T s=1/a;x11*=s;x12*=s;x22*=s;return *this;}

    UPPER_TRIANGULAR_MATRIX operator*(const T a) const
    {return UPPER_TRIANGULAR_MATRIX(a*x11,a*x12,a*x22);}

    UPPER_TRIANGULAR_MATRIX operator/(const T a) const
    {assert(a!=0);return *this*(1/a);}

    VECTOR<T,2> operator*(const VECTOR<T,2>& v) const
    {return VECTOR<T,2>(x11*v.x+x12*v.y,x22*v.y);}

    UPPER_TRIANGULAR_MATRIX operator*(const UPPER_TRIANGULAR_MATRIX& A) const // 4 mults, 1 add
    {return UPPER_TRIANGULAR_MATRIX(x11*A.x11,x11*A.x12+x12*A.x22,x22*A.x22);}

    UPPER_TRIANGULAR_MATRIX operator*(const DIAGONAL_MATRIX<T,2>& A) const // 3 mults
    {return UPPER_TRIANGULAR_MATRIX(x11*A.x11,x12*A.x22,x22*A.x22);}

    MATRIX_MXN<T> operator*(const MATRIX_MXN<T>& A) const
    {assert(A.Rows()==2);MATRIX_MXN<T> M(2,A.Columns());
    for(int j=1;j<=A.Columns();j++) for(int k=1;k<=2;k++) for(int i=1;i<=k;i++) M(i,j)+=(*this)(i,k)*A(k,j);return M;}

    MATRIX<T,2,3> Times_Transpose(const MATRIX<T,3,2>& A) const
    {return MATRIX<T,2,3>(x11*A.x[0]+x12*A.x[3],x22*A.x[3],x11*A.x[1]+x12*A.x[4],x22*A.x[4],x11*A.x[2]+x12*A.x[5],x22*A.x[5]);}

    SYMMETRIC_MATRIX<T,2> Outer_Product_Matrix() const // 4 mults, 1 add
    {return SYMMETRIC_MATRIX<T,2>(x11*x11+x12*x12,x12*x22,x22*x22);}

    T Determinant() const
    {return x11*x22;}

    T Trace() const
    {return x11+x22;}

    UPPER_TRIANGULAR_MATRIX Inverse() const
    {assert(x11!=0 && x22!=0);T one_over_x11=1/x11,one_over_x22=1/x22;
    return UPPER_TRIANGULAR_MATRIX(one_over_x11,-x12*one_over_x11*one_over_x22,one_over_x22);}

    VECTOR<T,2> Solve_Linear_System(const VECTOR<T,2>& b) const
    {return Inverse()*b;}

    UPPER_TRIANGULAR_MATRIX Cofactor_Matrix() const
    {return UPPER_TRIANGULAR_MATRIX(x22,-x12,x11);}

    static UPPER_TRIANGULAR_MATRIX Identity_Matrix()
    {return UPPER_TRIANGULAR_MATRIX(1,0,1);}

    T Max_Abs() const
    {return maxabs(x11,x12,x22);}

    T Frobenius_Norm() const
    {return sqrt(Frobenius_Norm_Squared());}

    T Frobenius_Norm_Squared() const
    {return sqr(x11)+sqr(x12)+sqr(x22);}

    T Simplex_Minimum_Altitude() const
    {return Determinant()/max(abs(x11),sqrt(sqr(x22)+max(sqr(x12),sqr(x12-x11))));}

//#####################################################################
};
// global functions
template<class T>
inline UPPER_TRIANGULAR_MATRIX<T,2> operator*(const T a,const UPPER_TRIANGULAR_MATRIX<T,2>& A) // 3 mults
{return A*a;}

template<class T>
inline UPPER_TRIANGULAR_MATRIX<T,2> operator+(const T a,const UPPER_TRIANGULAR_MATRIX<T,2>& A)
{return A+a;}

template<class T>
inline UPPER_TRIANGULAR_MATRIX<T,2> operator-(const T a,const UPPER_TRIANGULAR_MATRIX<T,2>& A)
{return -A+a;}

template<class T>
inline UPPER_TRIANGULAR_MATRIX<T,2> operator*(const DIAGONAL_MATRIX<T,2>& A,const UPPER_TRIANGULAR_MATRIX<T,2>& B) // 3 mults
{return UPPER_TRIANGULAR_MATRIX<T,2>(A.x11*B.x11,A.x11*B.x12,A.x22*B.x22);}

template<class T>
inline UPPER_TRIANGULAR_MATRIX<T,2> operator+(const DIAGONAL_MATRIX<T,2>& A,const UPPER_TRIANGULAR_MATRIX<T,2>& B)
{return B+A;}

template<class T>
inline UPPER_TRIANGULAR_MATRIX<T,2> operator-(const DIAGONAL_MATRIX<T,2>& A,const UPPER_TRIANGULAR_MATRIX<T,2>& B)
{return -B+A;}

template<class T>
inline std::ostream& operator<<(std::ostream& output_stream,const UPPER_TRIANGULAR_MATRIX<T,2>& A)
{return output_stream<<"["<<A.x11<<" "<<A.x12<<" ; 0 "<<A.x22<<"]";}
//#####################################################################
}
#endif


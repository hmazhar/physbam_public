//#####################################################################
// Copyright 2003-2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UPPER_TRIANGULAR_MATRIX_3X3
//#####################################################################
#ifndef __UPPER_TRIANGULAR_MATRIX_3X3__
#define __UPPER_TRIANGULAR_MATRIX_3X3__

#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <cassert>
#include <ostream>
namespace PhysBAM{

template<class T> struct IS_SCALAR_BLOCK<UPPER_TRIANGULAR_MATRIX<T,3> >:public IS_SCALAR_BLOCK<T>{};
template<class T> struct IS_SCALAR_VECTOR_SPACE<UPPER_TRIANGULAR_MATRIX<T,3> >:public IS_SCALAR_VECTOR_SPACE<T>{};
template<class T,class RW> struct IS_BINARY_IO_SAFE<UPPER_TRIANGULAR_MATRIX<T,3>,RW>:public IS_BINARY_IO_SAFE<T,RW>{};

template<class T>
class UPPER_TRIANGULAR_MATRIX<T,3>
{
public:
    typedef T SCALAR;
    enum WORKAROUND1 {m=3,n=3};

    T x11,x12,x22,x13,x23,x33;

    UPPER_TRIANGULAR_MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(3),INITIAL_SIZE nn=INITIAL_SIZE(3))
        :x11(0),x12(0),x22(0),x13(0),x23(0),x33(0)
    {
        STATIC_ASSERT(sizeof(UPPER_TRIANGULAR_MATRIX)==6*sizeof(T));assert(mm==INITIAL_SIZE(3) && nn==INITIAL_SIZE(3));
    }

    template<class T2> explicit
    UPPER_TRIANGULAR_MATRIX(const UPPER_TRIANGULAR_MATRIX<T2,3>& matrix_input)
        :x11(matrix_input.x11),x12(matrix_input.x12),x22(matrix_input.x22),x13(matrix_input.x13),x23(matrix_input.x23),x33(matrix_input.x33)
    {}

    UPPER_TRIANGULAR_MATRIX(const T x11_input,const T x12_input,const T x22_input,const T x13_input,const T x23_input,const T x33_input)
        :x11(x11_input),x12(x12_input),x22(x22_input),x13(x13_input),x23(x23_input),x33(x33_input)
    {}

    int Rows() const
    {return 3;}

    int Columns() const
    {return 3;}

    T& operator()(const int i,const int j)
    {assert(1<=i && i<=j && j<=3);return ((T*)this)[((j*(j-1))>>1)+i-1];}

    const T& operator()(const int i,const int j) const
    {assert(1<=i && i<=j && j<=3);return ((const T*)this)[((j*(j-1))>>1)+i-1];}

    bool Valid_Index(const int i,const int j) const
    {return 1<=i && i<=j && j<=3;}

    bool operator==(const UPPER_TRIANGULAR_MATRIX& A) const
    {return x11==A.x11 && x12==A.x12 && x22==A.x22 && x13==A.x13 && x23==A.x23 && x33==A.x33;}

    bool operator!=(const UPPER_TRIANGULAR_MATRIX& A) const
    {return !(*this==A);}

    UPPER_TRIANGULAR_MATRIX operator-() const
    {return UPPER_TRIANGULAR_MATRIX(-x11,-x12,-x22,-x13,-x23,-x33);}

    UPPER_TRIANGULAR_MATRIX& operator+=(const UPPER_TRIANGULAR_MATRIX& A)
    {x11+=A.x11;x12+=A.x12;x22+=A.x22;x13+=A.x13;x23+=A.x23;x33+=A.x33;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator+=(const T& a)
    {x11+=a;x22+=a;x33+=a;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator-=(const UPPER_TRIANGULAR_MATRIX& A)
    {x11-=A.x11;x12-=A.x12;x22-=A.x22;x13-=A.x13;x23-=A.x23;x33-=A.x33;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator-=(const T& a)
    {x11-=a;x22-=a;x33-=a;return *this;}

    UPPER_TRIANGULAR_MATRIX operator+(const UPPER_TRIANGULAR_MATRIX& A) const
    {return UPPER_TRIANGULAR_MATRIX(x11+A.x11,x12+A.x12,x22+A.x22,x13+A.x13,x23+A.x23,x33+A.x33);}

    UPPER_TRIANGULAR_MATRIX operator+(const DIAGONAL_MATRIX<T,3>& A) const
    {return UPPER_TRIANGULAR_MATRIX(x11+A.x11,x12,x22+A.x22,x13,x23,x33+A.x33);}

    UPPER_TRIANGULAR_MATRIX operator+(const T a) const
    {return UPPER_TRIANGULAR_MATRIX(x11+a,x12,x22+a,x13,x23,x33+a);}

    UPPER_TRIANGULAR_MATRIX operator-(const UPPER_TRIANGULAR_MATRIX& A) const
    {return UPPER_TRIANGULAR_MATRIX(x11-A.x11,x12-A.x12,x22-A.x22,x13-A.x13,x23-A.x23,x33-A.x33);}

    UPPER_TRIANGULAR_MATRIX operator-(const DIAGONAL_MATRIX<T,3>& A) const
    {return UPPER_TRIANGULAR_MATRIX(x11-A.x11,x12,x22-A.x22,x13,x23,x33-A.x33);}

    UPPER_TRIANGULAR_MATRIX operator-(const T a) const
    {return UPPER_TRIANGULAR_MATRIX(x11-a,x12,x22-a,x13,x23,x33-a);}

    UPPER_TRIANGULAR_MATRIX& operator*=(const UPPER_TRIANGULAR_MATRIX& A)
    {return *this=*this*A;}

    UPPER_TRIANGULAR_MATRIX& operator*=(const T a)
    {x11*=a;x12*=a;x22*=a;x13*=a;x23*=a;x33*=a;return *this;}

    UPPER_TRIANGULAR_MATRIX& operator/=(const T a)
    {assert(a!=0);T s=1/a;x11*=s;x12*=s;x22*=s;x13*=s;x23*=s;x33*=s;return *this;}

    UPPER_TRIANGULAR_MATRIX operator*(const T a) const
    {return UPPER_TRIANGULAR_MATRIX(a*x11,a*x12,a*x22,a*x13,a*x23,a*x33);}

    UPPER_TRIANGULAR_MATRIX operator/(const T a) const
    {assert(a!=0);T s=1/a;return UPPER_TRIANGULAR_MATRIX(s*x11,s*x12,s*x22,s*x13,s*x23,s*x33);}

    VECTOR<T,3> operator*(const VECTOR<T,3>& v) const // 6 mults, 3 adds
    {return VECTOR<T,3>(x11*v.x+x12*v.y+x13*v.z,x22*v.y+x23*v.z,x33*v.z);}

    UPPER_TRIANGULAR_MATRIX operator*(const UPPER_TRIANGULAR_MATRIX& A) const
    {return UPPER_TRIANGULAR_MATRIX(x11*A.x11,x11*A.x12+x12*A.x22,x22*A.x22,x11*A.x13+x12*A.x23+x13*A.x33,x22*A.x23+x23*A.x33,x33*A.x33);}

    UPPER_TRIANGULAR_MATRIX operator*(const DIAGONAL_MATRIX<T,3>& A) const
    {return UPPER_TRIANGULAR_MATRIX(x11*A.x11,x12*A.x22,x22*A.x22,x13*A.x33,x23*A.x33,x33*A.x33);}

    MATRIX_MXN<T> operator*(const MATRIX_MXN<T>& A) const
    {assert(A.Rows()==3);MATRIX_MXN<T> M(3,A.Columns());
    for(int j=1;j<=A.Columns();j++) for(int k=1;k<=3;k++) for(int i=1;i<=k;i++) M(i,j)+=(*this)(i,k)*A(k,j);return M;}

    T Determinant() const
    {return x11*x22*x33;}

    T Trace() const
    {return x11+x22+x33;}

    UPPER_TRIANGULAR_MATRIX Inverse() const
    {T determinant=x11*x22*x33;assert(determinant!=0);T s=1/determinant;
    return s*UPPER_TRIANGULAR_MATRIX(x22*x33,-x12*x33,x11*x33,x12*x23-x22*x13,-x11*x23,x11*x22);}

    VECTOR<T,3> Solve_Linear_System(const VECTOR<T,3>& b) const
    {return Cofactor_Matrix()*(b/Determinant());}

    UPPER_TRIANGULAR_MATRIX Cofactor_Matrix() const
    {return UPPER_TRIANGULAR_MATRIX(x22*x33,-x12*x33,x11*x33,x12*x23-x22*x13,-x11*x23,x11*x22);}

    static UPPER_TRIANGULAR_MATRIX Identity_Matrix()
    {return UPPER_TRIANGULAR_MATRIX(1,0,1,0,0,1);}

    T Max_Abs() const
    {return maxabs(x11,x12,x22,x13,x23,x33);}

    T Frobenius_Norm() const
    {return sqrt(Frobenius_Norm_Squared());}

    T Frobenius_Norm_Squared() const
    {return sqr(x11)+sqr(x12)+sqr(x22)+sqr(x13)+sqr(x23)+sqr(x33);}

    T Simplex_Minimum_Altitude() const
    {return MATRIX<T,3>(*this).Simplex_Minimum_Altitude();}

//#####################################################################
};
// global functions
template<class T>
inline UPPER_TRIANGULAR_MATRIX<T,3> operator*(const T a,const UPPER_TRIANGULAR_MATRIX<T,3>& A)
{return A*a;}

template<class T>
inline UPPER_TRIANGULAR_MATRIX<T,3> operator+(const T a,const UPPER_TRIANGULAR_MATRIX<T,3>& A)
{return A+a;}

template<class T>
inline UPPER_TRIANGULAR_MATRIX<T,3> operator-(const T a,const UPPER_TRIANGULAR_MATRIX<T,3>& A)
{return -A+a;}

template<class T>
inline UPPER_TRIANGULAR_MATRIX<T,3> operator*(const DIAGONAL_MATRIX<T,3>& A,const UPPER_TRIANGULAR_MATRIX<T,3>& B)
{return UPPER_TRIANGULAR_MATRIX<T,3>(A.x11*B.x11,A.x11*B.x12,A.x22*B.x22,A.x11*B.x13,A.x22*B.x23,A.x33*B.x33);}

template<class T>
inline UPPER_TRIANGULAR_MATRIX<T,3> operator+(const DIAGONAL_MATRIX<T,3>& A,const UPPER_TRIANGULAR_MATRIX<T,3>& B)
{return B+A;}

template<class T>
inline UPPER_TRIANGULAR_MATRIX<T,3> operator-(const DIAGONAL_MATRIX<T,3>& A,const UPPER_TRIANGULAR_MATRIX<T,3>& B)
{return -B+A;}

template<class T>
inline std::ostream& operator<<(std::ostream& output_stream,const UPPER_TRIANGULAR_MATRIX<T,3>& A)
{return output_stream<<"["<<A.x11<<" "<<A.x12<<" "<<A.x13<<" ; 0 "<<A.x22<<" "<<A.x23<<" ; 0 0 "<<A.x33<<"]";}
//#####################################################################
}
#endif

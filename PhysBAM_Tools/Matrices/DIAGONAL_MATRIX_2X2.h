//#####################################################################
// Copyright 2003-2007, Geoffrey Irving, Igor Neverov, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIAGONAL_MATRIX_2X2
//#####################################################################
#ifndef __DIAGONAL_MATRIX_2X2__
#define __DIAGONAL_MATRIX_2X2__

#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include <PhysBAM_Tools/Matrices/MATRIX_BASE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
namespace PhysBAM{

using ::std::log;

template<class T> struct IS_SCALAR_BLOCK<DIAGONAL_MATRIX<T,2> >:public IS_SCALAR_BLOCK<T>{};
template<class T> struct IS_SCALAR_VECTOR_SPACE<DIAGONAL_MATRIX<T,2> >:public IS_SCALAR_VECTOR_SPACE<T>{};
template<class T,class RW> struct IS_BINARY_IO_SAFE<DIAGONAL_MATRIX<T,2>,RW>:public IS_BINARY_IO_SAFE<T,RW>{};

template<class T>
class DIAGONAL_MATRIX<T,2>
{
public:
    typedef T SCALAR;
    enum WORKAROUND1 {m=2,n=2};

    T x11,x22;

    DIAGONAL_MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(2),INITIAL_SIZE nn=INITIAL_SIZE(2))
        :x11(0),x22(0)
    {
        STATIC_ASSERT(sizeof(DIAGONAL_MATRIX)==2*sizeof(T));assert(mm==INITIAL_SIZE(2) && nn==INITIAL_SIZE(2));
    }

    template<class T2> explicit
    DIAGONAL_MATRIX(const DIAGONAL_MATRIX<T2,2>& matrix_input)
        :x11((T)matrix_input.x11),x22((T)matrix_input.x22)
    {}

    DIAGONAL_MATRIX(const T y11,const T y22)
        :x11(y11),x22(y22)
    {}

    explicit DIAGONAL_MATRIX(const VECTOR<T,2>& v)
        :x11(v.x),x22(v.y)
    {}

    int Rows() const
    {return 2;}

    int Columns() const
    {return 2;}

    T& operator()(const int i)
    {assert(i>=1 && i<=2);return ((T*)this)[i-1];}

    const T& operator()(const int i) const
    {assert(i>=1 && i<=2);return ((T*)this)[i-1];}

    T& operator()(const int i,const int j)
    {assert(i>=1 && i<=2 && i==j);return ((T*)this)[i-1];}

    const T& operator()(const int i,const int j) const
    {assert(i>=1 && i<=2 && i==j);return ((T*)this)[i-1];}

    bool Valid_Index(const int i,const int j) const
    {return 1<=i && i<=2 && i==j;}

    T First() const
    {return x11;}

    T Last() const
    {return x22;}

    bool operator==(const DIAGONAL_MATRIX& A) const
    {return x11==A.x11 && x22==A.x22;}

    bool operator!=(const DIAGONAL_MATRIX& A) const
    {return x11!=A.x11 || x22!=A.x22;}

    DIAGONAL_MATRIX operator-() const
    {return DIAGONAL_MATRIX(-x11,-x22);}

    DIAGONAL_MATRIX& operator+=(const DIAGONAL_MATRIX& A)
    {x11+=A.x11;x22+=A.x22;return *this;}

    DIAGONAL_MATRIX& operator+=(const T& a)
    {x11+=a;x22+=a;return *this;}

    DIAGONAL_MATRIX& operator-=(const DIAGONAL_MATRIX& A)
    {x11-=A.x11;x22-=A.x22;return *this;}

    DIAGONAL_MATRIX& operator-=(const T& a)
    {x11-=a;x22-=a;return *this;}

    DIAGONAL_MATRIX& operator*=(const T a)
    {x11*=a;x22*=a;return *this;}

    DIAGONAL_MATRIX& operator/=(const T a)
    {assert(a!=0);T s=1/a;x11*=s;x22*=s;return *this;}

    DIAGONAL_MATRIX operator+(const DIAGONAL_MATRIX& A) const
    {return DIAGONAL_MATRIX(x11+A.x11,x22+A.x22);}

    MATRIX<T,2> operator+(const MATRIX<T,2>& A) const
    {return MATRIX<T,2>(x11+A.x[0],A.x[1],A.x[2],x22+A.x[3]);}

    MATRIX<T,2> operator-(const MATRIX<T,2>& A) const
    {return MATRIX<T,2>(x11-A.x[0],-A.x[1],-A.x[2],x22-A.x[3]);}

    DIAGONAL_MATRIX operator+(const T a) const
    {return DIAGONAL_MATRIX(x11+a,x22+a);}

    DIAGONAL_MATRIX operator-(const DIAGONAL_MATRIX& A) const
    {return DIAGONAL_MATRIX(x11-A.x11,x22-A.x22);}

    DIAGONAL_MATRIX operator-(const T a) const
    {return DIAGONAL_MATRIX(x11-a,x22-a);}

    DIAGONAL_MATRIX operator*(const T a) const
    {return DIAGONAL_MATRIX(a*x11,a*x22);}

    DIAGONAL_MATRIX operator/(const T a) const
    {assert(a!=0);return *this*(1/a);}

    VECTOR<T,2> operator*(const VECTOR<T,2>& v) const
    {return VECTOR<T,2>(x11*v.x,x22*v.y);}

    DIAGONAL_MATRIX operator*(const DIAGONAL_MATRIX& A) const
    {return DIAGONAL_MATRIX(x11*A.x11,x22*A.x22);}

    DIAGONAL_MATRIX& operator*=(const DIAGONAL_MATRIX& A)
    {return *this=*this*A;}

    DIAGONAL_MATRIX operator/(const DIAGONAL_MATRIX& A) const
    {return DIAGONAL_MATRIX(x11/A.x11,x22/A.x22);}

    T Determinant() const
    {return x11*x22;}

    DIAGONAL_MATRIX Inverse() const
    {assert(x11!=0 && x22!=0);return DIAGONAL_MATRIX(1/x11,1/x22);}

    VECTOR<T,2> Solve_Linear_System(const VECTOR<T,2>& v) const
    {assert(x11!=0 && x22!=0);return VECTOR<T,2>(v.x/x11,v.y/x22);}

    VECTOR<T,2> Robust_Solve_Linear_System(const VECTOR<T,2>& v) const
    {return VECTOR<T,2>(Robust_Divide(v.x,x11),Robust_Divide(v.y,x22));}

    DIAGONAL_MATRIX Transposed() const
    {return *this;}

    void Transpose()
    {}

    DIAGONAL_MATRIX Cofactor_Matrix() const
    {return DIAGONAL_MATRIX(x22,x11);}

    T Trace() const
    {return x11+x22;}

    T Dilational() const
    {return (T).5*Trace();}

    T Min() const
    {return min(x11,x22);}

    T Max() const
    {return max(x11,x22);}

    T Max_Abs() const
    {return maxabs(x11,x22);}

    static T Inner_Product(const DIAGONAL_MATRIX& A,const DIAGONAL_MATRIX& B)
    {return A.x11*B.x11+A.x22*B.x22;}

    T Inner_Product(const VECTOR<T,2>& a,const VECTOR<T,2>& b) const // inner product with respect to this matrix
    {return a.x*x11*b.x+a.y*x22*b.y;}

    T Inverse_Inner_Product(const VECTOR<T,2>& a,const VECTOR<T,2>& b) const // inner product with respect to the inverse of this matrix
    {assert(x11!=0 && x22!=0);return a.x/x11*b.x+a.y/x22*b.y;}

    T Frobenius_Norm_Squared() const
    {return sqr(x11)+sqr(x22);}

    T Frobenius_Norm() const
    {return sqrt(Frobenius_Norm_Squared());}

    bool Positive_Definite() const
    {return x11>0 && x22>0;}

    bool Positive_Semidefinite() const
    {return x11>=0 && x22>=0;}

    DIAGONAL_MATRIX Positive_Definite_Part() const
    {return Clamp_Min(0);}

    DIAGONAL_MATRIX Sqrt() const
    {return DIAGONAL_MATRIX(sqrt(x11),sqrt(x22));}

    DIAGONAL_MATRIX Clamp_Min(const T a) const
    {return DIAGONAL_MATRIX(clamp_min(x11,a),clamp_min(x22,a));}

    DIAGONAL_MATRIX Clamp_Max(const T a) const
    {return DIAGONAL_MATRIX(clamp_max(x11,a),clamp_max(x22,a));}

    DIAGONAL_MATRIX Abs() const
    {return DIAGONAL_MATRIX(abs(x11),abs(x22));}

    DIAGONAL_MATRIX Sign() const
    {return DIAGONAL_MATRIX(sign(x11),sign(x22));}

    static DIAGONAL_MATRIX Identity_Matrix()
    {return DIAGONAL_MATRIX(1,1);}

    VECTOR<T,2> To_Vector() const
    {return VECTOR<T,2>(x11,x22);}

    template<class T_MATRIX>
    typename PRODUCT<DIAGONAL_MATRIX,T_MATRIX>::TYPE Transpose_Times(const T_MATRIX& M) const
    {return *this*M;}

    template<class T_MATRIX>
    typename PRODUCT_TRANSPOSE<DIAGONAL_MATRIX,T_MATRIX>::TYPE Times_Transpose(const MATRIX_BASE<T,T_MATRIX>& B) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);assert(B.Columns()==2);typename PRODUCT_TRANSPOSE<DIAGONAL_MATRIX,T_MATRIX>::TYPE M((INITIAL_SIZE)B.Columns(),(INITIAL_SIZE)B.Rows());
    for(int k=1;k<=B.Rows();k++) for(int i=1;i<=B.Columns();i++) M(i,k)=(*this)(i,i)*B(k,i);return M;}

    DIAGONAL_MATRIX Times_Transpose(const DIAGONAL_MATRIX& M) const
    {return *this*M;}

//#####################################################################
    MATRIX<T,2> Times_Transpose(const MATRIX<T,2>& A) const;
    MATRIX<T,2,3> Times_Transpose(const MATRIX<T,3,2>& A) const;
    static T Inner_Product_Conjugate(const DIAGONAL_MATRIX& A,const MATRIX<T,2>& Q,const DIAGONAL_MATRIX B);
//#####################################################################
};
// global functions
template<class T>
inline DIAGONAL_MATRIX<T,2> operator*(const T a,const DIAGONAL_MATRIX<T,2>& A)
{return A*a;}

template<class T>
inline DIAGONAL_MATRIX<T,2> operator+(const T a,const DIAGONAL_MATRIX<T,2>& A)
{return A+a;}

template<class T>
inline DIAGONAL_MATRIX<T,2> operator-(const T a,const DIAGONAL_MATRIX<T,2>& A)
{return -A+a;}

template<class T>
inline MATRIX<T,2> operator+(const MATRIX<T,2>& A,const DIAGONAL_MATRIX<T,2>& B)
{return B+A;}

template<class T>
inline MATRIX<T,2> operator-(const MATRIX<T,2>& A,const DIAGONAL_MATRIX<T,2>& B)
{return -B+A;}

template<class T>
inline std::ostream& operator<<(std::ostream& output_stream,const DIAGONAL_MATRIX<T,2>& A)
{return output_stream<<"["<<A.x11<<" 0 ; 0 "<<A.x22<<"]";}

template<class T>
inline DIAGONAL_MATRIX<T,2> log(const DIAGONAL_MATRIX<T,2>& A)
{return DIAGONAL_MATRIX<T,2>(log(A.x11),log(A.x22));}

template<class T>
inline DIAGONAL_MATRIX<T,2> exp(const DIAGONAL_MATRIX<T,2>& A)
{return DIAGONAL_MATRIX<T,2>(exp(A.x11),exp(A.x22));}

//#####################################################################
// Function Times_Transpose
//#####################################################################
template<class T> inline MATRIX<T,2,3> DIAGONAL_MATRIX<T,2>::
Times_Transpose(const MATRIX<T,3,2>& A) const
{
    return MATRIX<T,2,3>(x11*A.x[0],x22*A.x[3],x11*A.x[1],x22*A.x[4],x11*A.x[2],x22*A.x[5]);
}
//#####################################################################
// Function Times_Transpose
//#####################################################################
template<class T> inline MATRIX<T,2> DIAGONAL_MATRIX<T,2>::
Times_Transpose(const MATRIX<T,2>& A) const
{
    return MATRIX<T,2>(x11*A.x[0],x22*A.x[2],x11*A.x[1],x22*A.x[3]);
}
//#####################################################################
// Function Inner_Product_Conjugate
//#####################################################################
template<class T> inline T DIAGONAL_MATRIX<T,2>::
Inner_Product_Conjugate(const DIAGONAL_MATRIX<T,2>& A,const MATRIX<T,2>& Q,const DIAGONAL_MATRIX<T,2> B)
{
    MATRIX<T,2> BQ=B*Q.Transposed();
    return A.x11*(Q.x[0]*BQ.x[0]+Q.x[2]*BQ.x[1])+A.x22*(Q.x[1]*BQ.x[2]+Q.x[3]*BQ.x[3]);
}
//#####################################################################
}
#endif

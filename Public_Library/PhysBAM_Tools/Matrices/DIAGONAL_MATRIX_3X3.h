//#####################################################################
// Copyright 2003-2007, Geoffrey Irving, Neil Molino, Craig Schroeder, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIAGONAL_MATRIX_3X3
//#####################################################################
#ifndef __DIAGONAL_MATRIX_3X3__
#define __DIAGONAL_MATRIX_3X3__

#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Tools/Matrices/MATRIX_BASE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

using ::std::log;
using ::std::exp;

template<class T> struct IS_SCALAR_BLOCK<DIAGONAL_MATRIX<T,3> >:public IS_SCALAR_BLOCK<T>{};
template<class T> struct IS_SCALAR_VECTOR_SPACE<DIAGONAL_MATRIX<T,3> >:public IS_SCALAR_VECTOR_SPACE<T>{};
template<class T,class RW> struct IS_BINARY_IO_SAFE<DIAGONAL_MATRIX<T,3>,RW>:public IS_BINARY_IO_SAFE<T,RW>{};

template<class T>
class DIAGONAL_MATRIX<T,3>
{
public:
    typedef T SCALAR;
    enum WORKAROUND1 {m=3,n=3};

    T x11,x22,x33;

    DIAGONAL_MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(3),INITIAL_SIZE nn=INITIAL_SIZE(3))
        :x11(0),x22(0),x33(0)
    {
        STATIC_ASSERT(sizeof(DIAGONAL_MATRIX)==3*sizeof(T));assert(mm==INITIAL_SIZE(3) && nn==INITIAL_SIZE(3));
    }

    template<class T2> explicit
    DIAGONAL_MATRIX(const DIAGONAL_MATRIX<T2,3>& matrix_input)
        :x11((T)matrix_input.x11),x22((T)matrix_input.x22),x33((T)matrix_input.x33)
    {}

    DIAGONAL_MATRIX(const T y11,const T y22,const T y33)
        :x11(y11),x22(y22),x33(y33)
    {}

    explicit DIAGONAL_MATRIX(const VECTOR<T,3>& v)
        :x11(v.x),x22(v.y),x33(v.z)
    {}

    int Rows() const
    {return 3;}

    int Columns() const
    {return 3;}

    T& operator()(const int i)
    {assert(i>=1 && i<=3);return ((T*)this)[i-1];}

    const T& operator()(const int i) const
    {assert(i>=1 && i<=3);return ((T*)this)[i-1];}

    T& operator()(const int i,const int j)
    {assert(i>=1 && i<=3 && i==j);return ((T*)this)[i-1];}

    const T& operator()(const int i,const int j) const
    {assert(i>=1 && i<=3 && i==j);return ((T*)this)[i-1];}

    bool Valid_Index(const int i,const int j) const
    {return 1<=i && i<=3 && i==j;}

    T First() const
    {return x11;}

    T Last() const
    {return x33;}

    bool operator==(const DIAGONAL_MATRIX& A) const
    {return x11==A.x11 && x22==A.x22 && x33==A.x33;}

    bool operator!=(const DIAGONAL_MATRIX& A) const
    {return x11!=A.x11 || x22!=A.x22 || x33!=A.x33;}

    DIAGONAL_MATRIX operator-() const
    {return DIAGONAL_MATRIX(-x11,-x22,-x33);}

    DIAGONAL_MATRIX& operator+=(const DIAGONAL_MATRIX& A)
    {x11+=A.x11;x22+=A.x22;x33+=A.x33;return *this;}

    DIAGONAL_MATRIX& operator+=(const T& a)
    {x11+=a;x22+=a;x33+=a;return *this;}

    DIAGONAL_MATRIX& operator-=(const DIAGONAL_MATRIX& A)
    {x11-=A.x11;x22-=A.x22;x33-=A.x33;return *this;}

    DIAGONAL_MATRIX& operator-=(const T& a)
    {x11-=a;x22-=a;x33-=a;return *this;}

    DIAGONAL_MATRIX& operator*=(const T a)
    {x11*=a;x22*=a;x33*=a;return *this;}

    DIAGONAL_MATRIX& operator/=(const T a)
    {assert(a!=0);T s=1/a;x11*=s;x22*=s;x33*=s;return *this;}

    DIAGONAL_MATRIX operator+(const DIAGONAL_MATRIX& A) const
    {return DIAGONAL_MATRIX(x11+A.x11,x22+A.x22,x33+A.x33);}

    MATRIX<T,3> operator+(const MATRIX<T,3>& A) const
    {return MATRIX<T,3>(x11+A.x[0],A.x[1],A.x[2],A.x[3],x22+A.x[4],A.x[5],A.x[6],A.x[7],x33+A.x[8]);}

    MATRIX<T,3> operator-(const MATRIX<T,3>& A) const
    {return MATRIX<T,3>(x11-A.x[0],-A.x[1],-A.x[2],-A.x[3],x22-A.x[4],-A.x[5],-A.x[6],-A.x[7],x33-A.x[8]);}

    DIAGONAL_MATRIX operator+(const T a) const
    {return DIAGONAL_MATRIX(x11+a,x22+a,x33+a);}

    DIAGONAL_MATRIX operator-(const DIAGONAL_MATRIX& A) const
    {return DIAGONAL_MATRIX(x11-A.x11,x22-A.x22,x33-A.x33);}

    DIAGONAL_MATRIX operator-(const T a) const
    {return DIAGONAL_MATRIX(x11-a,x22-a,x33-a);}

    DIAGONAL_MATRIX operator*(const T a) const
    {return DIAGONAL_MATRIX(a*x11,a*x22,a*x33);}

    DIAGONAL_MATRIX operator/(const T a) const
    {assert(a!=0);return *this*(1/a);}

    VECTOR<T,3> operator*(const VECTOR<T,3>& v) const
    {return VECTOR<T,3>(x11*v.x,x22*v.y,x33*v.z);}

    DIAGONAL_MATRIX operator*(const DIAGONAL_MATRIX& A) const
    {return DIAGONAL_MATRIX(x11*A.x11,x22*A.x22,x33*A.x33);}

    DIAGONAL_MATRIX& operator*=(const DIAGONAL_MATRIX& A)
    {return *this=*this*A;}

    DIAGONAL_MATRIX operator/(const DIAGONAL_MATRIX& A) const
    {return DIAGONAL_MATRIX(x11/A.x11,x22/A.x22,x33/A.x33);}

    T Determinant() const
    {return x11*x22*x33;}

    DIAGONAL_MATRIX Inverse() const
    {assert(x11!=0 && x22!=0 && x33!=0);return DIAGONAL_MATRIX(1/x11,1/x22,1/x33);}

    VECTOR<T,3> Solve_Linear_System(const VECTOR<T,3>& v) const
    {assert(x11!=0 && x22!=0 && x33!=0);return VECTOR<T,3>(v.x/x11,v.y/x22,v.z/x33);}

    VECTOR<T,3> Robust_Solve_Linear_System(const VECTOR<T,3>& v) const
    {return VECTOR<T,3>(Robust_Divide(v.x,x11),Robust_Divide(v.y,x22),Robust_Divide(v.z,x33));}

    template<class T_MATRIX>
    typename PRODUCT_TRANSPOSE<DIAGONAL_MATRIX,T_MATRIX>::TYPE Times_Transpose(const MATRIX_BASE<T,T_MATRIX>& B) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);assert(B.Columns()==3);typename PRODUCT_TRANSPOSE<DIAGONAL_MATRIX,T_MATRIX>::TYPE M((INITIAL_SIZE)B.Columns(),(INITIAL_SIZE)B.Rows());
    for(int k=1;k<=B.Rows();k++) for(int i=1;i<=B.Columns();i++) M(i,k)=(*this)(i,i)*B(k,i);return M;}

    DIAGONAL_MATRIX Times_Transpose(const DIAGONAL_MATRIX& M) const
    {return *this*M;}

    DIAGONAL_MATRIX Transposed() const
    {return *this;}

    void Transpose()
    {}

    template<class T_MATRIX>
    typename PRODUCT<DIAGONAL_MATRIX,T_MATRIX>::TYPE Transpose_Times(const T_MATRIX& M) const
    {return *this*M;}

    SYMMETRIC_MATRIX<T,3> Conjugate_With_Cross_Product_Matrix(const VECTOR<T,3>& v) const
    {T yy=sqr(v.y),zz=sqr(v.z),bx=x22*v.x,cx=x33*v.x;return SYMMETRIC_MATRIX<T,3>(x22*zz+x33*yy,-cx*v.y,-bx*v.z,x11*zz+cx*v.x,-v.y*v.z*x11,x11*yy+bx*v.x);}

    DIAGONAL_MATRIX Cofactor_Matrix() const
    {return DIAGONAL_MATRIX(x22*x33,x11*x33,x11*x22);}

    T Trace() const
    {return x11+x22+x33;}

    T Dilational() const
    {return (T)one_third*Trace();}

    T Min() const
    {return min(x11,x22,x33);}

    T Max() const
    {return max(x11,x22,x33);}

    T Max_Abs() const
    {return maxabs(x11,x22,x33);}

    static T Inner_Product(const DIAGONAL_MATRIX& A,const DIAGONAL_MATRIX& B)
    {return A.x11*B.x11+A.x22*B.x22+A.x33*B.x33;}

    T Inner_Product(const VECTOR<T,3>& a,const VECTOR<T,3>& b) const // inner product with respect to this matrix
    {return a.x*x11*b.x+a.y*x22*b.y+a.z*x33*b.z;}

    T Inverse_Inner_Product(const VECTOR<T,3>& a,const VECTOR<T,3>& b) const // inner product with respect to the inverse of this matrix
    {assert(x11!=0 && x22!=0 && x33!=0);return a.x/x11*b.x+a.y/x22*b.y+a.z/x33*b.z;}

    T Frobenius_Norm_Squared() const
    {return sqr(x11)+sqr(x22)+sqr(x33);}

    T Frobenius_Norm() const
    {return sqrt(Frobenius_Norm_Squared());}

    bool Positive_Definite() const
    {return x11>0 && x22>0 && x33>0;}

    bool Positive_Semidefinite() const
    {return x11>=0 && x22>=0 && x33>=0;}

    DIAGONAL_MATRIX Positive_Definite_Part() const
    {return Clamp_Min(0);}

    DIAGONAL_MATRIX Sqrt() const
    {return DIAGONAL_MATRIX(sqrt(x11),sqrt(x22),sqrt(x33));}

    DIAGONAL_MATRIX Clamp_Min(const T a) const
    {return DIAGONAL_MATRIX(clamp_min(x11,a),clamp_min(x22,a),clamp_min(x33,a));}

    DIAGONAL_MATRIX Clamp_Max(const T a) const
    {return DIAGONAL_MATRIX(clamp_max(x11,a),clamp_max(x22,a),clamp_max(x33,a));}

    DIAGONAL_MATRIX Abs() const
    {return DIAGONAL_MATRIX(abs(x11),abs(x22),abs(x33));}

    DIAGONAL_MATRIX Sign() const
    {return DIAGONAL_MATRIX(sign(x11),sign(x22),sign(x33));}

    static DIAGONAL_MATRIX Identity_Matrix()
    {return DIAGONAL_MATRIX(1,1,1);}

    VECTOR<T,3> To_Vector() const
    {return VECTOR<T,3>(x11,x22,x33);}

//#####################################################################
    MATRIX<T,3> Times_Transpose(const MATRIX<T,3>& A) const;
    static T Inner_Product_Conjugate(const DIAGONAL_MATRIX& A,const MATRIX<T,3>& Q,const DIAGONAL_MATRIX B);
//#####################################################################
};
// global functions
template<class T>
inline DIAGONAL_MATRIX<T,3> operator*(const T a,const DIAGONAL_MATRIX<T,3>& A)
{return A*a;}

template<class T>
inline DIAGONAL_MATRIX<T,3> operator+(const T a,const DIAGONAL_MATRIX<T,3>& A)
{return A+a;}

template<class T>
inline DIAGONAL_MATRIX<T,3> operator-(const T a,const DIAGONAL_MATRIX<T,3>& A)
{return -A+a;}

template<class T>
inline MATRIX<T,3> operator+(const MATRIX<T,3>& A,const DIAGONAL_MATRIX<T,3>& B)
{return B+A;}

template<class T>
inline MATRIX<T,3> operator-(const MATRIX<T,3>& A,const DIAGONAL_MATRIX<T,3>& B)
{return -B+A;}

template<class T>
inline std::istream& operator>>(std::istream& input_stream,DIAGONAL_MATRIX<T,3>& A)
{return input_stream>>A.x11>>A.x22>>A.x33;}

template<class T>
inline std::ostream& operator<<(std::ostream& output_stream,const DIAGONAL_MATRIX<T,3>& A)
{return output_stream<<"["<<A.x11<<" 0 0 ; 0 "<<A.x22<<" 0 ; 0 0 "<<A.x33<<"]";}

template<class T>
inline DIAGONAL_MATRIX<T,3> log(const DIAGONAL_MATRIX<T,3>& A)
{return DIAGONAL_MATRIX<T,3>(log(A.x11),log(A.x22),log(A.x33));}

template<class T>
inline DIAGONAL_MATRIX<T,3> exp(const DIAGONAL_MATRIX<T,3>& A)
{return DIAGONAL_MATRIX<T,3>(exp(A.x11),exp(A.x22),exp(A.x33));}

//#####################################################################
// Function Times_Transpose
//#####################################################################
template<class T> inline MATRIX<T,3> DIAGONAL_MATRIX<T,3>::
Times_Transpose(const MATRIX<T,3>& A) const
{
    return MATRIX<T,3>(x11*A.x[0],x22*A.x[3],x33*A.x[6],x11*A.x[1],x22*A.x[4],x33*A.x[7],x11*A.x[2],x22*A.x[5],x33*A.x[8]);
}
//#####################################################################
// Function Inner_Product_Conjugate
//#####################################################################
template<class T> inline T DIAGONAL_MATRIX<T,3>::
Inner_Product_Conjugate(const DIAGONAL_MATRIX<T,3>& A,const MATRIX<T,3>& Q,const DIAGONAL_MATRIX<T,3> B)
{
    MATRIX<T,3> BQ=B*Q.Transposed();
    return A.x11*(Q.x[0]*BQ.x[0]+Q.x[3]*BQ.x[1]+Q.x[6]*BQ.x[2])+A.x22*(Q.x[1]*BQ.x[3]+Q.x[4]*BQ.x[4]+Q.x[7]*BQ.x[5])+A.x33*(Q.x[2]*BQ.x[6]+Q.x[5]*BQ.x[7]+Q.x[8]*BQ.x[8]);
}
//#####################################################################
}
#endif


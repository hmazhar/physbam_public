//#####################################################################
// Copyright 2003-2007, Geoffrey Irving, Craig Schroeder, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMMETRIC_MATRIX_2X2
//#####################################################################
#ifndef __SYMMETRIC_MATRIX_2X2__
#define __SYMMETRIC_MATRIX_2X2__

#include <PhysBAM_Tools/Math_Tools/exchange_sort.h>
#include <PhysBAM_Tools/Matrices/MATRIX_ARITHMETIC_POLICY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
namespace PhysBAM{

template<class T> struct IS_SCALAR_BLOCK<SYMMETRIC_MATRIX<T,2> >:public IS_SCALAR_BLOCK<T>{};
template<class T> struct IS_SCALAR_VECTOR_SPACE<SYMMETRIC_MATRIX<T,2> >:public IS_SCALAR_VECTOR_SPACE<T>{};
template<class T,class RW> struct IS_BINARY_IO_SAFE<SYMMETRIC_MATRIX<T,2>,RW>:public IS_BINARY_IO_SAFE<T,RW>{};

template<class T>
class SYMMETRIC_MATRIX<T,2>
{
public:
    typedef T SCALAR;
    enum WORKAROUND1 {m=2,n=2};

    T x11,x21,x22;

    SYMMETRIC_MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(2),INITIAL_SIZE nn=INITIAL_SIZE(2))
        :x11(0),x21(0),x22(0)
    {
        assert(mm==INITIAL_SIZE(2) && nn==INITIAL_SIZE(2));
    }

    template<class T2> explicit
    SYMMETRIC_MATRIX(const SYMMETRIC_MATRIX<T2,2>& matrix_input)
        :x11((T)matrix_input.x11),x21((T)matrix_input.x21),x22((T)matrix_input.x22)
    {}

    SYMMETRIC_MATRIX(const DIAGONAL_MATRIX<T,2>& matrix_input)
        :x11(matrix_input.x11),x21(0),x22(matrix_input.x22)
    {}

    SYMMETRIC_MATRIX(const T y11,const T y21,const T y22)
        :x11(y11),x21(y21),x22(y22)
    {}

    int Rows() const
    {return 2;}

    int Columns() const
    {return 2;}

    VECTOR<T,2> Column(const int axis) const
    {assert(1<=axis && axis<=2);return axis==1?VECTOR<T,2>(x11,x21):VECTOR<T,2>(x21,x22);}

    T& operator()(int i,int j)
    {return i<j?Element_Upper(i,j):Element_Lower(i,j);}

    const T& operator()(int i,int j) const
    {return i<j?Element_Upper(i,j):Element_Lower(i,j);}

    bool Valid_Index(const int i,const int j) const
    {return 1<=i && i<=2 && 1<=j && j<=2;}

    T& Element_Upper(int i,int j)
    {return Element_Lower(j,i);}

    const T& Element_Upper(int i,int j) const
    {return Element_Lower(j,i);}

    T& Element_Lower(int i,int j)
    {assert(i<=2 && j>=1 && j<=i);return ((T*)this)[((4-j)*(j-1)>>1)+i-1];}

    const T& Element_Lower(int i,int j) const
    {assert(i<=2 && j>=1 && j<=i);return ((const T*)this)[((4-j)*(j-1)>>1)+i-1];}

    bool operator==(const SYMMETRIC_MATRIX& A) const
    {return x11==A.x11 && x21==A.x21 && x22==A.x22;}

    bool operator!=(const SYMMETRIC_MATRIX& A) const
    {return !(*this==A);}

    static SYMMETRIC_MATRIX Componentwise_Min(const SYMMETRIC_MATRIX& v1,const SYMMETRIC_MATRIX& v2)
    {return SYMMETRIC_MATRIX(min(v1.x11,v2.x11),min(v1.x21,v2.x21),min(v1.x22,v2.x22));}

    static SYMMETRIC_MATRIX Componentwise_Max(const SYMMETRIC_MATRIX& v1,const SYMMETRIC_MATRIX& v2)
    {return SYMMETRIC_MATRIX(max(v1.x11,v2.x11),max(v1.x21,v2.x21),max(v1.x22,v2.x22));}

    SYMMETRIC_MATRIX operator-() const
    {return SYMMETRIC_MATRIX(-x11,-x21,-x22);}

    SYMMETRIC_MATRIX& operator+=(const SYMMETRIC_MATRIX& A)
    {x11+=A.x11;x21+=A.x21;x22+=A.x22;return *this;}

    SYMMETRIC_MATRIX& operator+=(const T& a)
    {x11+=a;x22+=a;return *this;}

    SYMMETRIC_MATRIX& operator-=(const SYMMETRIC_MATRIX& A)
    {x11-=A.x11;x21-=A.x21;x22-=A.x22;return *this;}

    SYMMETRIC_MATRIX& operator-=(const T& a)
    {x11-=a;x22-=a;return *this;}

    SYMMETRIC_MATRIX& operator*=(const T a)
    {x11*=a;x21*=a;x22*=a;return *this;}

    SYMMETRIC_MATRIX& operator/=(const T a)
    {assert(a!=0);T s=1/a;x11*=s;x21*=s;x22*=s;return *this;}

    SYMMETRIC_MATRIX operator+(const SYMMETRIC_MATRIX& A) const
    {return SYMMETRIC_MATRIX(x11+A.x11,x21+A.x21,x22+A.x22);}

    SYMMETRIC_MATRIX operator+(const T a) const
    {return SYMMETRIC_MATRIX(x11+a,x21,x22+a);}

    SYMMETRIC_MATRIX operator-(const SYMMETRIC_MATRIX& A) const
    {return SYMMETRIC_MATRIX(x11-A.x11,x21-A.x21,x22-A.x22);}

    SYMMETRIC_MATRIX operator-(const T a) const
    {return SYMMETRIC_MATRIX(x11-a,x21,x22-a);}

    SYMMETRIC_MATRIX operator*(const T a) const
    {return SYMMETRIC_MATRIX(a*x11,a*x21,a*x22);}

    SYMMETRIC_MATRIX operator/(const T a) const
    {assert(a!=0);return *this*(1/a);}

    VECTOR<T,2> operator*(const VECTOR<T,2>& v) const
    {return VECTOR<T,2>(x11*v.x+x21*v.y,x21*v.x+x22*v.y);}

    template<class T_MATRIX>
    typename PRODUCT<SYMMETRIC_MATRIX,T_MATRIX>::TYPE Transpose_Times(const T_MATRIX& M) const
    {return *this*M;}

    template<class T_MATRIX>
    typename PRODUCT_TRANSPOSE<SYMMETRIC_MATRIX,T_MATRIX>::TYPE Times_Transpose(const MATRIX_BASE<T,T_MATRIX>& A) const
    {typename PRODUCT_TRANSPOSE<SYMMETRIC_MATRIX,T_MATRIX>::TYPE M((INITIAL_SIZE)A.Columns(),(INITIAL_SIZE)A.Rows());A.Add_Times_Transpose(*this,A.Derived(),M);return M;}

    T Determinant() const
    {return x11*x22-x21*x21;}

    SYMMETRIC_MATRIX Inverse() const
    {return SYMMETRIC_MATRIX(x22,-x21,x11)/Determinant();}

    SYMMETRIC_MATRIX Transposed() const
    {return *this;}

    void Transpose()
    {}

    MATRIX<T,2> Times_Transpose(const MATRIX<T,2>& A) const
    {return *this*A.Transposed();}

    MATRIX<T,2,3> Times_Transpose(const MATRIX<T,3,2>& A) const
    {return MATRIX<T,2,3>(x11*A.x[0]+x21*A.x[3],x21*A.x[0]+x22*A.x[3],x11*A.x[1]+x21*A.x[4],x21*A.x[1]+x22*A.x[4],x11*A.x[2]+x21*A.x[5],x21*A.x[2]+x22*A.x[5]);}

    MATRIX<T,2> Times_Transpose(const DIAGONAL_MATRIX<T,2>& A) const
    {return *this*A;}

    MATRIX<T,2> Times_Transpose(const SYMMETRIC_MATRIX<T,2>& A) const
    {return *this*A;}

    MATRIX<T,1,2> Cross_Product_Matrix_Times(const VECTOR<T,2>& v) const
    {return MATRIX<T,1,2>::Cross_Product_Matrix(v)*(*this);}

    VECTOR<T,2> Solve_Linear_System(const VECTOR<T,2>& b) const
    {return SYMMETRIC_MATRIX(x22,-x21,x11)*b/Determinant();}

    VECTOR<T,2> Robust_Solve_Linear_System(const VECTOR<T,2>& b) const
    {T determinant=Determinant();
    VECTOR<T,2> unscaled_result=SYMMETRIC_MATRIX(x22,-x21,x11)*b;
    T relative_tolerance=(T)FLT_MIN*unscaled_result.Max_Abs();
    if(abs(determinant)<=relative_tolerance){relative_tolerance=max(relative_tolerance,(T)FLT_MIN);determinant=determinant>=0?relative_tolerance:-relative_tolerance;}
    return unscaled_result/determinant;}

    T Trace() const
    {return x11+x22;}

    static T Inner_Product(const SYMMETRIC_MATRIX& A,const SYMMETRIC_MATRIX& B)
    {return A.x11*B.x11+A.x22*B.x22+2*A.x21*B.x21;}

    T Frobenius_Norm_Squared() const
    {return x11*x11+x22*x22+2*x21*x21;}

    T Frobenius_Norm() const
    {return sqrt(Frobenius_Norm_Squared());}

    SYMMETRIC_MATRIX Cofactor_Matrix()
    {return SYMMETRIC_MATRIX(x22,-x21,x11);}

    VECTOR<T,2> Largest_Column() const
    {return abs(x11)>abs(x22)?VECTOR<T,2>(x11,x21):VECTOR<T,2>(x21,x22);}

    VECTOR<T,2> Largest_Column_Normalized() const // 5 mults, 2 adds, 1 div, 1 sqrt
    {T sqr11=sqr(x11),sqr12=sqr(x21),sqr22=sqr(x22);
    T scale1=sqr11+sqr12,scale2=sqr12+sqr22;
    if(scale1>scale2) return VECTOR<T,2>(x11,x21)/sqrt(scale1);
    else if(scale2>0) return VECTOR<T,2>(x21,x22)/sqrt(scale2);
    else return VECTOR<T,2>(1,0);}

    T Max_Abs() const
    {return maxabs(x11,x21,x22);}

    static SYMMETRIC_MATRIX Outer_Product(const VECTOR<T,2>& u)
    {return SYMMETRIC_MATRIX(u.x*u.x,u.x*u.y,u.y*u.y);}

    static SYMMETRIC_MATRIX Identity_Matrix()
    {return SYMMETRIC_MATRIX(1,0,1);}

    static SYMMETRIC_MATRIX Unit_Matrix(const T scale=1)
    {return SYMMETRIC_MATRIX(scale,scale,scale);}

    bool Positive_Definite() const
    {return x11>0 && x11*x22>x21*x21;}

    bool Positive_Semidefinite(const T tolerance=(T)1e-7) const
    {T scale=Max_Abs();return !scale || (*this+tolerance*scale).Positive_Definite();}

    VECTOR<T,2> First_Eigenvector_From_Ordered_Eigenvalues(const DIAGONAL_MATRIX<T,2>& eigenvalues) const
    {return (*this-eigenvalues.x11).Cofactor_Matrix().Largest_Column_Normalized();}

    VECTOR<T,2> Last_Eigenvector_From_Ordered_Eigenvalues(const DIAGONAL_MATRIX<T,2>& eigenvalues) const
    {return (*this-eigenvalues.x22).Cofactor_Matrix().Largest_Column_Normalized();}

    DIAGONAL_MATRIX<T,2> Fast_Eigenvalues() const
    {T da;
    if(x21==0) da=0;
    else{T theta=(T).5*(x22-x11)/x21,t=1/(abs(theta)+sqrt(1+sqr(theta)));if(theta<0) t=-t;da=t*x21;}
    DIAGONAL_MATRIX<T,2> eigenvalues(x11-da,x22+da);
    exchange_sort(eigenvalues.x11,eigenvalues.x22);return eigenvalues;}

    SYMMETRIC_MATRIX Positive_Definite_Part() const
    {DIAGONAL_MATRIX<T,2> D;MATRIX<T,2> V;Solve_Eigenproblem(D,V);D=D.Clamp_Min(0);return Conjugate(V,D);}

    DIAGONAL_MATRIX<T,2> Diagonal_Part() const
    {return DIAGONAL_MATRIX<T,2>(x11,x22);}

    void Fast_Solve_Eigenproblem(DIAGONAL_MATRIX<T,2>& eigenvalues,MATRIX<T,2>& eigenvectors) const
    {Solve_Eigenproblem(eigenvalues,eigenvectors);}

    static SYMMETRIC_MATRIX Transpose_Times_With_Symmetric_Result(const MATRIX<T,2>& A,const MATRIX<T,2>& B) // A^t*B and assume symmetric result, 6 mults, 3 adds
    {return SYMMETRIC_MATRIX(A.x[0]*B.x[0]+A.x[1]*B.x[1],A.x[2]*B.x[0]+A.x[3]*B.x[1],A.x[2]*B.x[2]+A.x[3]*B.x[3]);}

    static SYMMETRIC_MATRIX Transpose_Times_With_Symmetric_Result(const MATRIX<T,3,2>& A,const MATRIX<T,3,2>& B) // A^t*B and assume symmetric result, 9 mults, 6 adds
    {return SYMMETRIC_MATRIX(A.x[0]*B.x[0]+A.x[1]*B.x[1]+A.x[2]*B.x[2],A.x[3]*B.x[0]+A.x[4]*B.x[1]+A.x[5]*B.x[2],A.x[3]*B.x[3]+A.x[4]*B.x[4]+A.x[5]*B.x[5]);}

    static SYMMETRIC_MATRIX Transpose_Times_With_Symmetric_Result(const MATRIX<T,2>& A,const UPPER_TRIANGULAR_MATRIX<T,2>& B) // A^t*B and assume symmetric result, 4 mults, 1 adds
    {return SYMMETRIC_MATRIX(A.x[0]*B.x11,A.x[2]*B.x11,A.x[2]*B.x12+A.x[3]*B.x22);}

//#####################################################################
    void Solve_Eigenproblem(DIAGONAL_MATRIX<T,2>& eigenvalues,MATRIX<T,2>& eigenvectors) const;
    MATRIX<T,2> operator*(const DIAGONAL_MATRIX<T,2>& A) const;
    MATRIX<T,2> operator*(const UPPER_TRIANGULAR_MATRIX<T,2>& A) const;
    SYMMETRIC_MATRIX operator+(const DIAGONAL_MATRIX<T,2>& A) const;
    MATRIX<T,2> Times_Transpose(const UPPER_TRIANGULAR_MATRIX<T,2>& A) const;
    static SYMMETRIC_MATRIX Conjugate(const MATRIX<T,2>& A,const DIAGONAL_MATRIX<T,2>& B);
    static SYMMETRIC_MATRIX Conjugate(const MATRIX<T,2>& A,const SYMMETRIC_MATRIX& B);
    static SYMMETRIC_MATRIX Conjugate(const MATRIX<T,2,3>& A,const DIAGONAL_MATRIX<T,3>& B);
    static SYMMETRIC_MATRIX Conjugate(const MATRIX<T,2,3>& A,const SYMMETRIC_MATRIX<T,3>& B);
    static SYMMETRIC_MATRIX Conjugate(const UPPER_TRIANGULAR_MATRIX<T,2>& A,const SYMMETRIC_MATRIX& B);
    static SYMMETRIC_MATRIX Conjugate_With_Transpose(const MATRIX<T,2>& A,const DIAGONAL_MATRIX<T,2>& B);
    static SYMMETRIC_MATRIX Conjugate_With_Transpose(const MATRIX<T,2>& A,const SYMMETRIC_MATRIX& B);
    static SYMMETRIC_MATRIX Conjugate_With_Transpose(const UPPER_TRIANGULAR_MATRIX<T,2>& A,const SYMMETRIC_MATRIX& B);
    static SYMMETRIC_MATRIX Conjugate_With_Transpose(const MATRIX<T,3,2>& A,const SYMMETRIC_MATRIX<T,3>& B);
//#####################################################################
};
// global functions
template<class T>
inline SYMMETRIC_MATRIX<T,2> operator*(const T a,const SYMMETRIC_MATRIX<T,2>& A) // 4 mults
{return A*a;}

template<class T>
inline SYMMETRIC_MATRIX<T,2> operator+(const T a,const SYMMETRIC_MATRIX<T,2>& A) // 2 adds
{return A+a;}

template<class T>
inline SYMMETRIC_MATRIX<T,2> operator-(const T a,const SYMMETRIC_MATRIX<T,2>& A) // 2 adds
{return -A+a;}

template<class T>
inline SYMMETRIC_MATRIX<T,2> clamp(const SYMMETRIC_MATRIX<T,2>& x,const SYMMETRIC_MATRIX<T,2>& xmin,const SYMMETRIC_MATRIX<T,2>& xmax)
{return SYMMETRIC_MATRIX<T,2>(clamp(x.x11,xmin.x11,xmax.x11),clamp(x.x21,xmin.x21,xmax.x21),clamp(x.x22,xmin.x22,xmax.x22));}

template<class T>
inline SYMMETRIC_MATRIX<T,2> clamp_min(const SYMMETRIC_MATRIX<T,2>& x,const SYMMETRIC_MATRIX<T,2>& xmin)
{return SYMMETRIC_MATRIX<T,2>(clamp_min(x.x11,xmin.x11),clamp_min(x.x21,xmin.x21),clamp_min(x.x22,xmin.x22));}

template<class T>
inline SYMMETRIC_MATRIX<T,2> clamp_max(const SYMMETRIC_MATRIX<T,2>& x,const SYMMETRIC_MATRIX<T,2>& xmax)
{return SYMMETRIC_MATRIX<T,2>(clamp_max(x.x11,xmax.x11),clamp_max(x.x21,xmax.x21),clamp_max(x.x22,xmax.x22));}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
template<class T>
inline std::ostream& operator<< (std::ostream& output_stream,const SYMMETRIC_MATRIX<T,2>& A)
{output_stream<<"["<<A.x11<<" "<<A.x21<<" ; "<<A.x21<<" "<<A.x22<<"]";return output_stream;}
#endif

template<class T>
inline SYMMETRIC_MATRIX<T,2> log(const SYMMETRIC_MATRIX<T,2>& A)
{DIAGONAL_MATRIX<T,2> D;MATRIX<T,2> Q;A.Solve_Eigenproblem(D,Q);return SYMMETRIC_MATRIX<T,2>::Conjugate(Q,log(D));}

template<class T>
inline SYMMETRIC_MATRIX<T,2> exp(const SYMMETRIC_MATRIX<T,2>& A)
{DIAGONAL_MATRIX<T,2> D;MATRIX<T,2> Q;A.Solve_Eigenproblem(D,Q);return SYMMETRIC_MATRIX<T,2>::Conjugate(Q,exp(D));}

//#####################################################################
// Function Solve_Eigenproblem
//#####################################################################
template<class T> void SYMMETRIC_MATRIX<T,2>::
Solve_Eigenproblem(DIAGONAL_MATRIX<T,2>& eigenvalues,MATRIX<T,2>& eigenvectors) const
{
    typedef VECTOR<T,2> TV;
    T a=(T).5*(x11+x22),b=(T).5*(x11-x22),c=x21;
    T c_squared=sqr(c),m=sqrt(sqr(b)+c_squared),k=x11*x22-c_squared;
    if(a>=0){eigenvalues.x11=a+m;eigenvalues.x22=eigenvalues.x11?k/eigenvalues.x11:0;}
    else{eigenvalues.x22=a-m;eigenvalues.x11=eigenvalues.x22?k/eigenvalues.x22:0;}
    exchange_sort(eigenvalues.x22,eigenvalues.x11); // if order is wrong, matrix is nearly scalar
    eigenvectors.Column(1)=(b>=0?TV(m+b,c):TV(-c,b-m)).Normalized();
    eigenvectors.Column(2)=eigenvectors.Column(1).Perpendicular();
}
//#####################################################################
// Function Conjugate
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,2> SYMMETRIC_MATRIX<T,2>::
Conjugate(const MATRIX<T,2>& A,const DIAGONAL_MATRIX<T,2>& B) // 10 mults, 3 adds
{
    MATRIX<T,2> BA=B*A.Transposed();
    return SYMMETRIC_MATRIX<T,2>(A.x[0]*BA.x[0]+A.x[2]*BA.x[1],A.x[1]*BA.x[0]+A.x[3]*BA.x[1],A.x[1]*BA.x[2]+A.x[3]*BA.x[3]);
}
//#####################################################################
// Function Conjugate
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,2> SYMMETRIC_MATRIX<T,2>::
Conjugate(const MATRIX<T,2>& A,const SYMMETRIC_MATRIX<T,2>& B) // 12 mults, 7 adds
{
    MATRIX<T,2> BA=(A*B).Transposed();
    return SYMMETRIC_MATRIX<T,2>(A.x[0]*BA.x[0]+A.x[2]*BA.x[1],A.x[1]*BA.x[0]+A.x[3]*BA.x[1],A.x[1]*BA.x[2]+A.x[3]*BA.x[3]);
}
//#####################################################################
// Function Conjugate
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,2> SYMMETRIC_MATRIX<T,2>::
Conjugate(const MATRIX<T,2,3>& A,const DIAGONAL_MATRIX<T,3>& B)
{
    return Transpose_Times_With_Symmetric_Result(A.transpose,B*A.transpose);
}
//#####################################################################
// Function Conjugate
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,2> SYMMETRIC_MATRIX<T,2>::
Conjugate(const MATRIX<T,2,3>& A,const SYMMETRIC_MATRIX<T,3>& B)
{
    return Transpose_Times_With_Symmetric_Result(A.transpose,B*A.transpose);
}
//#####################################################################
// Function Conjugate
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,2> SYMMETRIC_MATRIX<T,2>::
Conjugate(const UPPER_TRIANGULAR_MATRIX<T,2>& A,const SYMMETRIC_MATRIX<T,2>& B) // 
{
    MATRIX<T,2> BA=B.Times_Transpose(A);
    return SYMMETRIC_MATRIX<T,2>(A.x11*BA.x[0]+A.x12*BA.x[1],A.x22*BA.x[1],A.x22*BA.x[3]);
}
//#####################################################################
// Function Conjugate_With_Transpose
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,2> SYMMETRIC_MATRIX<T,2>::
Conjugate_With_Transpose(const MATRIX<T,2>& A,const DIAGONAL_MATRIX<T,2>& B) // 10 mults, 3 adds
{
    return Transpose_Times_With_Symmetric_Result(B*A,A);
}
//#####################################################################
// Function Conjugate_With_Transpose
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,2> SYMMETRIC_MATRIX<T,2>::
Conjugate_With_Transpose(const MATRIX<T,2>& A,const SYMMETRIC_MATRIX<T,2>& B) // 12 mults, 7 adds
{
    return Transpose_Times_With_Symmetric_Result(B*A,A);
}
//#####################################################################
// Function Conjugate_With_Transpose
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,2> SYMMETRIC_MATRIX<T,2>::
Conjugate_With_Transpose(const UPPER_TRIANGULAR_MATRIX<T,2>& A,const SYMMETRIC_MATRIX<T,2>& B) // 10 mults, 3 adds
{
    return Transpose_Times_With_Symmetric_Result(B*A,A);
}
//#####################################################################
// Function Conjugate_With_Transpose
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,2> SYMMETRIC_MATRIX<T,2>::
Conjugate_With_Transpose(const MATRIX<T,3,2>& A,const SYMMETRIC_MATRIX<T,3>& B) // 21 mults, 12 adds
{
    return Transpose_Times_With_Symmetric_Result(B*A,A);
}
//#####################################################################
// Function Times_Transpose
//#####################################################################
template<class T> inline MATRIX<T,2> SYMMETRIC_MATRIX<T,2>::
Times_Transpose(const UPPER_TRIANGULAR_MATRIX<T,2>& A) const // 6 mults, 2 adds
{
    return (A**this).Transposed();
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T> inline MATRIX<T,2> SYMMETRIC_MATRIX<T,2>::
operator*(const DIAGONAL_MATRIX<T,2>& A) const // 4 mults
{
    return MATRIX<T,2>(x11*A.x11,x21*A.x11,x21*A.x22,x22*A.x22);
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T> inline MATRIX<T,2> SYMMETRIC_MATRIX<T,2>::
operator*(const UPPER_TRIANGULAR_MATRIX<T,2>& A) const // 6 mults, 2 adds
{
    return MATRIX<T,2>(x11*A.x11,x21*A.x11,x11*A.x12+x21*A.x22,x21*A.x12+x22*A.x22);
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T> inline MATRIX<T,2> operator*(const DIAGONAL_MATRIX<T,2>& D,const SYMMETRIC_MATRIX<T,2>& A) // 4 mults
{
    return MATRIX<T,2>(D.x11*A.x11,D.x22*A.x21,D.x11*A.x21,D.x22*A.x22);
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T>
inline MATRIX<T,2> operator*(const UPPER_TRIANGULAR_MATRIX<T,2>& A,const SYMMETRIC_MATRIX<T,2>& B) // 6 mults, 2 adds
{
    return MATRIX<T,2>(A.x11*B.x11+A.x12*B.x21,A.x22*B.x21,A.x11*B.x21+A.x12*B.x22,A.x22*B.x22);
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T>
inline MATRIX<T,2> operator*(const SYMMETRIC_MATRIX<T,2>& A,const MATRIX<T,2>& B) // 8 mults, 4 mults
{
    return MATRIX<T,2>(A.x11*B.x[0]+A.x21*B.x[1],A.x21*B.x[0]+A.x22*B.x[1],A.x11*B.x[2]+A.x21*B.x[3],A.x21*B.x[2]+A.x22*B.x[3]);
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T>
inline MATRIX<T,2> operator*(const SYMMETRIC_MATRIX<T,2>& A,const SYMMETRIC_MATRIX<T,2>& B) // 8 mults, 4 adds
{
    return MATRIX<T,2>(A.x11*B.x11+A.x21*B.x21,A.x21*B.x11+A.x22*B.x21,A.x11*B.x21+A.x21*B.x22,A.x21*B.x21+A.x22*B.x22);
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,2> SYMMETRIC_MATRIX<T,2>::
operator+(const DIAGONAL_MATRIX<T,2>& A) const // 2 adds
{
    return SYMMETRIC_MATRIX<T,2>(x11+A.x11,x21,x22+A.x22);
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,2>
operator+(const DIAGONAL_MATRIX<T,2>& A,const SYMMETRIC_MATRIX<T,2>& B) // 2 adds
{
    return B+A;
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline MATRIX<T,2>
operator+(const SYMMETRIC_MATRIX<T,2>& A,const UPPER_TRIANGULAR_MATRIX<T,2>& B)
{
    return MATRIX<T,2>(A.x11+B.x11,A.x21,A.x21+B.x12,A.x22+B.x22);
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline MATRIX<T,2>
operator+(const UPPER_TRIANGULAR_MATRIX<T,2>& A,const SYMMETRIC_MATRIX<T,2>& B)
{
    return B+A;
}
//#####################################################################
// Function operator-
//#####################################################################
template<class T> inline MATRIX<T,2>
operator-(const SYMMETRIC_MATRIX<T,2>& A,const UPPER_TRIANGULAR_MATRIX<T,2>& B)
{
    return MATRIX<T,2>(A.x11-B.x11,A.x21,A.x21-B.x12,A.x22-B.x22);
}
//#####################################################################
// Function operator-
//#####################################################################
template<class T> inline MATRIX<T,2>
operator-(const UPPER_TRIANGULAR_MATRIX<T,2>& A,const SYMMETRIC_MATRIX<T,2>& B)
{
    return -B+A;
}
//#####################################################################
// Function operator-
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,2>
operator-(const DIAGONAL_MATRIX<T,2>& A,const SYMMETRIC_MATRIX<T,2>& B) // 2 adds
{
    return SYMMETRIC_MATRIX<T,2>(A.x11-B.x11,-B.x21,A.x22-B.x22);
}
//#####################################################################
}
#endif

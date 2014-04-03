//#####################################################################
// Copyright 2003-2007, Ronald Fedkiw, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMMETRIC_MATRIX_3X3
//#####################################################################
#ifndef __SYMMETRIC_MATRIX_3X3__
#define __SYMMETRIC_MATRIX_3X3__

#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Matrices/MATRIX_ARITHMETIC_POLICY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T> struct HAS_CHEAP_COPY<SYMMETRIC_MATRIX<T,3> > {static const bool value=true;};
template<class T> struct IS_SCALAR_BLOCK<SYMMETRIC_MATRIX<T,3> >:public IS_SCALAR_BLOCK<T>{};
template<class T> struct IS_SCALAR_VECTOR_SPACE<SYMMETRIC_MATRIX<T,3> >:public IS_SCALAR_VECTOR_SPACE<T>{};
template<class T,class RW> struct IS_BINARY_IO_SAFE<SYMMETRIC_MATRIX<T,3>,RW>:public IS_BINARY_IO_SAFE<T,RW>{};

template<class T>
class SYMMETRIC_MATRIX<T,3>
{
public:
    typedef T SCALAR;
    enum WORKAROUND1 {m=3,n=3};

    T x11,x21,x31,x22,x32,x33;

    SYMMETRIC_MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(3),INITIAL_SIZE nn=INITIAL_SIZE(3))
        :x11(0),x21(0),x31(0),x22(0),x32(0),x33(0)
    {
        STATIC_ASSERT(sizeof(SYMMETRIC_MATRIX)==6*sizeof(T));assert(mm==INITIAL_SIZE(3) && nn==INITIAL_SIZE(3));
    }

    template<class T2> explicit
    SYMMETRIC_MATRIX(const SYMMETRIC_MATRIX<T2,3>& matrix_input)
        :x11((T)matrix_input.x11),x21((T)matrix_input.x21),x31((T)matrix_input.x31),x22((T)matrix_input.x22),x32((T)matrix_input.x32),x33((T)matrix_input.x33)
    {}

    SYMMETRIC_MATRIX(const DIAGONAL_MATRIX<T,3>& matrix_input)
        :x11(matrix_input.x11),x21(0),x31(0),x22(matrix_input.x22),x32(0),x33(matrix_input.x33)
    {}

    SYMMETRIC_MATRIX(const T y11,const T y21,const T y31,const T y22,const T y32,const T y33)
        :x11(y11),x21(y21),x31(y31),x22(y22),x32(y32),x33(y33)
    {}

    void From_Matrix(const MATRIX<T,3>& matrix_input)
    {
        x11=matrix_input(1,1);
        x21=matrix_input(2,1);
        x31=matrix_input(3,1);
        x22=matrix_input(2,2);
        x32=matrix_input(3,2);
        x33=matrix_input(3,3);
    }

    int Rows() const
    {return 3;}

    int Columns() const
    {return 3;}

    T& operator()(int i,int j)
    {return i<j?Element_Upper(i,j):Element_Lower(i,j);}

    const T& operator()(int i,int j) const
    {return i<j?Element_Upper(i,j):Element_Lower(i,j);}

    bool Valid_Index(const int i,const int j) const
    {return 1<=i && i<=3 && 1<=j && j<=3;}

    T& Element_Upper(int i,int j)
    {return Element_Lower(j,i);}

    const T& Element_Upper(int i,int j) const
    {return Element_Lower(j,i);}

    T& Element_Lower(int i,int j)
    {assert(i<=3 && j>=1 && j<=i);return ((T*)this)[((6-j)*(j-1)>>1)+i-1];}

    const T& Element_Lower(int i,int j) const
    {assert(i<=3 && j>=1 && j<=i);return ((const T*)this)[((6-j)*(j-1)>>1)+i-1];}

    VECTOR<T,3> Column(const int axis) const
    {assert(1<=axis && axis<=3);return axis==1?VECTOR<T,3>(x11,x21,x31):axis==2?VECTOR<T,3>(x21,x22,x32):VECTOR<T,3>(x31,x32,x33);}

    bool operator==(const SYMMETRIC_MATRIX& A) const
    {return x11==A.x11 && x21==A.x21 && x31==A.x31 && x22==A.x22 && x32==A.x32 && x33==A.x33;}

    bool operator!=(const SYMMETRIC_MATRIX& A) const
    {return !(*this==A);}

    static SYMMETRIC_MATRIX Componentwise_Min(const SYMMETRIC_MATRIX& v1,const SYMMETRIC_MATRIX& v2)
    {return SYMMETRIC_MATRIX(min(v1.x11,v2.x11),min(v1.x21,v2.x21),min(v1.x31,v2.x31),min(v1.x22,v2.x22),min(v1.x32,v2.x32),min(v1.x33,v2.x33));}

    static SYMMETRIC_MATRIX Componentwise_Max(const SYMMETRIC_MATRIX& v1,const SYMMETRIC_MATRIX& v2)
    {return SYMMETRIC_MATRIX(max(v1.x11,v2.x11),max(v1.x21,v2.x21),max(v1.x31,v2.x31),max(v1.x22,v2.x22),max(v1.x32,v2.x32),max(v1.x33,v2.x33));}

    SYMMETRIC_MATRIX operator-() const
    {return SYMMETRIC_MATRIX(-x11,-x21,-x31,-x22,-x32,-x33);}

    SYMMETRIC_MATRIX& operator+=(const SYMMETRIC_MATRIX& A)
    {x11+=A.x11;x21+=A.x21;x31+=A.x31;x22+=A.x22;x32+=A.x32;x33+=A.x33;return *this;}

    SYMMETRIC_MATRIX& operator+=(const T& a)
    {x11+=a;x22+=a;x33+=a;return *this;}

    SYMMETRIC_MATRIX& operator-=(const SYMMETRIC_MATRIX& A)
    {x11-=A.x11;x21-=A.x21;x31-=A.x31;x22-=A.x22;x32-=A.x32;x33-=A.x33;return *this;}

    SYMMETRIC_MATRIX& operator-=(const T& a)
    {x11-=a;x22-=a;x33-=a;return *this;}

    SYMMETRIC_MATRIX& operator*=(const T a)
    {x11*=a;x21*=a;x31*=a;x22*=a;x32*=a;x33*=a;return *this;}

    SYMMETRIC_MATRIX& operator/=(const T a)
    {assert(a!=0);T s=1/a;x11*=s;x21*=s;x31*=s;x22*=s;x32*=s;x33*=s;return *this;}

    SYMMETRIC_MATRIX operator+(const SYMMETRIC_MATRIX& A) const
    {return SYMMETRIC_MATRIX(x11+A.x11,x21+A.x21,x31+A.x31,x22+A.x22,x32+A.x32,x33+A.x33);}

    SYMMETRIC_MATRIX operator+(const T a) const
    {return SYMMETRIC_MATRIX(x11+a,x21,x31,x22+a,x32,x33+a);}

    SYMMETRIC_MATRIX operator-(const SYMMETRIC_MATRIX& A) const
    {return SYMMETRIC_MATRIX(x11-A.x11,x21-A.x21,x31-A.x31,x22-A.x22,x32-A.x32,x33-A.x33);}

    SYMMETRIC_MATRIX operator-(const T a) const
    {return SYMMETRIC_MATRIX(x11-a,x21,x31,x22-a,x32,x33-a);}

    SYMMETRIC_MATRIX operator*(const T a) const
    {return SYMMETRIC_MATRIX(a*x11,a*x21,a*x31,a*x22,a*x32,a*x33);}

    MATRIX<T,3> operator*(const SYMMETRIC_MATRIX& A) const // 27 mults, 18 adds
    {return MATRIX<T,3>(x11*A.x11+x21*A.x21+x31*A.x31,x21*A.x11+x22*A.x21+x32*A.x31,x31*A.x11+x32*A.x21+x33*A.x31,
                        x11*A.x21+x21*A.x22+x31*A.x32,x21*A.x21+x22*A.x22+x32*A.x32,x31*A.x21+x32*A.x22+x33*A.x32,
                        x11*A.x31+x21*A.x32+x31*A.x33,x21*A.x31+x22*A.x32+x32*A.x33,x31*A.x31+x32*A.x32+x33*A.x33);}

    MATRIX_MXN<T> operator*(const MATRIX_MXN<T>& A) const
    {assert(3==A.m);MATRIX_MXN<T> matrix(3,A.n);for(int j=1;j<=A.n;j++) matrix.Set_Column(j,(*this)*VECTOR<T,3>(A(1,j),A(2,j),A(3,j)));return matrix;}

    template<class T_MATRIX>
    typename PRODUCT<SYMMETRIC_MATRIX,T_MATRIX>::TYPE Transpose_Times(const T_MATRIX& M) const
    {return *this*M;}

    template<class T_MATRIX>
    typename PRODUCT_TRANSPOSE<SYMMETRIC_MATRIX,T_MATRIX>::TYPE Times_Transpose(const MATRIX_BASE<T,T_MATRIX>& A) const
    {typename PRODUCT_TRANSPOSE<SYMMETRIC_MATRIX,T_MATRIX>::TYPE M((INITIAL_SIZE)A.Columns(),(INITIAL_SIZE)A.Rows());A.Add_Times_Transpose(*this,A.Derived(),M);return M;}

    MATRIX<T,3> Times_Transpose(const SYMMETRIC_MATRIX& M) const // 27 mults, 18 adds
    {return *this*M;}

    MATRIX<T,3> Cross_Product_Matrix_Times(const VECTOR<T,3>& v) const // (v*) * (*this)
    {return MATRIX<T,3>(-v.z*x21+v.y*x31,v.z*x11-v.x*x31,-v.y*x11+v.x*x21,-v.z*x22+v.y*x32,v.z*x21-v.x*x32,-v.y*x21+v.x*x22,-v.z*x32+v.y*x33,v.z*x31-v.x*x33,-v.y*x31+v.x*x32);}

    MATRIX<T,3> Cross_Product_Matrix_Transpose_Times(const VECTOR<T,3>& v) const // (v*)^T * (*this)
    {return MATRIX<T,3>(v.z*x21-v.y*x31,-v.z*x11+v.x*x31,v.y*x11-v.x*x21,v.z*x22-v.y*x32,-v.z*x21+v.x*x32,v.y*x21-v.x*x22,v.z*x32-v.y*x33,-v.z*x31+v.x*x33,v.y*x31-v.x*x32);}

    MATRIX<T,3> Times_Cross_Product_Matrix(const VECTOR<T,3>& v) const // (*this) * (v*)
    {return MATRIX<T,3>(x21*v.z-x31*v.y,x22*v.z-x32*v.y,x32*v.z-x33*v.y,-x11*v.z+x31*v.x,-x21*v.z+x32*v.x,-x31*v.z+x33*v.x,x11*v.y-x21*v.x,x21*v.y-x22*v.x,x31*v.y-x32*v.x);}

    SYMMETRIC_MATRIX operator/(const T a) const
    {assert(a!=0);T s=1/a;return SYMMETRIC_MATRIX(s*x11,s*x21,s*x31,s*x22,s*x32,s*x33);}

    VECTOR<T,3> operator*(const VECTOR<T,3>& v) const
    {return VECTOR<T,3>(x11*v.x+x21*v.y+x31*v.z,x21*v.x+x22*v.y+x32*v.z,x31*v.x+x32*v.y+x33*v.z);}

    T Determinant() const
    {return x11*(x22*x33-x32*x32)+x21*(2*x32*x31-x21*x33)-x31*x22*x31;}

    SYMMETRIC_MATRIX Inverse() const
    {T cofactor11=x22*x33-x32*x32,cofactor12=x32*x31-x21*x33,cofactor13=x21*x32-x22*x31;
    return SYMMETRIC_MATRIX(cofactor11,cofactor12,cofactor13,x11*x33-x31*x31,x21*x31-x11*x32,x11*x22-x21*x21)/(x11*cofactor11+x21*cofactor12+x31*cofactor13);}

    SYMMETRIC_MATRIX Transposed() const
    {return *this;}

    void Transpose()
    {}

    T Dilational() const
    {return (T)one_third*Trace();}

    SYMMETRIC_MATRIX Deviatoric() const
    {return *this-Dilational();}

    VECTOR<T,3> Solve_Linear_System(const VECTOR<T,3>& b) const // 18 mults, 8 adds
    {T cofactor11=x22*x33-x32*x32,cofactor12=x32*x31-x21*x33,cofactor13=x21*x32-x22*x31;
    return SYMMETRIC_MATRIX(cofactor11,cofactor12,cofactor13,x11*x33-x31*x31,x21*x31-x11*x32,x11*x22-x21*x21)*b/(x11*cofactor11+x21*cofactor12+x31*cofactor13);}

    VECTOR<T,3> Robust_Solve_Linear_System(const VECTOR<T,3>& b) const
    {T cofactor11=x22*x33-x32*x32,cofactor12=x32*x31-x21*x33,cofactor13=x21*x32-x22*x31;
    T determinant=x11*cofactor11+x21*cofactor12+x31*cofactor13;
    VECTOR<T,3> unscaled_result=SYMMETRIC_MATRIX(cofactor11,cofactor12,cofactor13,x11*x33-x31*x31,x21*x31-x11*x32,x11*x22-x21*x21)*b;
    T relative_tolerance=(T)FLT_MIN*unscaled_result.Max_Abs();
    if(abs(determinant)<=relative_tolerance){relative_tolerance=max(relative_tolerance,(T)FLT_MIN);determinant=determinant>=0?relative_tolerance:-relative_tolerance;}
    return unscaled_result/determinant;}

    SYMMETRIC_MATRIX Squared() const
    {return SYMMETRIC_MATRIX(x11*x11+x21*x21+x31*x31,x21*x11+x22*x21+x32*x31,x31*x11+x32*x21+x33*x31,x21*x21+x22*x22+x32*x32,x31*x21+x32*x22+x33*x32,x31*x31+x32*x32+x33*x33);}

    T Trace() const
    {return x11+x22+x33;}

    static T Inner_Product(const SYMMETRIC_MATRIX& A,const SYMMETRIC_MATRIX& B)
    {return A.x11*B.x11+A.x22*B.x22+A.x33*B.x33+2*(A.x21*B.x21+A.x31*B.x31+A.x32*B.x32);}

    T Frobenius_Norm_Squared() const
    {return x11*x11+x22*x22+x33*x33+2*(x21*x21+x31*x31+x32*x32);}

    T Frobenius_Norm() const
    {return sqrt(Frobenius_Norm_Squared());}

    SYMMETRIC_MATRIX Cofactor_Matrix() // 12 mults, 6 adds
    {return SYMMETRIC_MATRIX(x22*x33-x32*x32,x32*x31-x21*x33,x21*x32-x22*x31,x11*x33-x31*x31,x21*x31-x11*x32,x11*x22-x21*x21);}

    VECTOR<T,3> Largest_Column() const
    {T sqr11=sqr(x11),sqr12=sqr(x21),sqr13=sqr(x31),sqr22=sqr(x22),sqr23=sqr(x32),sqr33=sqr(x33);
    T scale1=sqr11+sqr12+sqr13,scale2=sqr12+sqr22+sqr23,scale3=sqr13+sqr23+sqr33;
    return scale1>scale2?(scale1>scale3?VECTOR<T,3>(x11,x21,x31):VECTOR<T,3>(x31,x32,x33)):(scale2>scale3?VECTOR<T,3>(x21,x22,x32):VECTOR<T,3>(x31,x32,x33));}

    VECTOR<T,3> Largest_Column_Normalized() const // 9 mults, 6 adds, 1 div, 1 sqrt
    {T sqr11=sqr(x11),sqr12=sqr(x21),sqr13=sqr(x31),sqr22=sqr(x22),sqr23=sqr(x32),sqr33=sqr(x33);
    T scale1=sqr11+sqr12+sqr13,scale2=sqr12+sqr22+sqr23,scale3=sqr13+sqr23+sqr33;
    if(scale1>scale2){if(scale1>scale3) return VECTOR<T,3>(x11,x21,x31)/sqrt(scale1);}
    else if(scale2>scale3) return VECTOR<T,3>(x21,x22,x32)/sqrt(scale2);
    if(scale3>0) return VECTOR<T,3>(x31,x32,x33)/sqrt(scale3);else return VECTOR<T,3>(1,0,0);}

    T Max_Abs() const
    {return maxabs(x11,x21,x31,x22,x32,x33);}

    static SYMMETRIC_MATRIX Outer_Product(const VECTOR<T,3>& u) // 6 mults
    {return SYMMETRIC_MATRIX(u.x*u.x,u.x*u.y,u.x*u.z,u.y*u.y,u.y*u.z,u.z*u.z);}

    static SYMMETRIC_MATRIX Identity_Matrix()
    {return SYMMETRIC_MATRIX(1,0,0,1,0,1);}

    static SYMMETRIC_MATRIX Unit_Matrix(const T scale=1)
    {return SYMMETRIC_MATRIX(scale,scale,scale,scale,scale,scale);}

    bool Positive_Definite() const
    {return x11>0 && x11*x22>x21*x21 && Determinant()>0;}

    bool Positive_Semidefinite(const T tolerance=(T)1e-7) const
    {T scale=Max_Abs();return !scale || (*this+tolerance*scale).Positive_Definite();}

    VECTOR<T,3> First_Eigenvector_From_Ordered_Eigenvalues(const DIAGONAL_MATRIX<T,3>& eigenvalues,const T tolerance=1e-5) const
    {T scale=maxabs(eigenvalues.x11,eigenvalues.x33),scale_inverse=Robust_Inverse(scale),tiny=tolerance*scale;
    if(eigenvalues.x11-eigenvalues.x22>tiny) return ((*this-eigenvalues.x11)*scale_inverse).Cofactor_Matrix().Largest_Column_Normalized();
    return ((*this-eigenvalues.x33)*scale_inverse).Cofactor_Matrix().Largest_Column().Unit_Orthogonal_Vector();}

    VECTOR<T,3> Last_Eigenvector_From_Ordered_Eigenvalues(const DIAGONAL_MATRIX<T,3>& eigenvalues,const T tolerance=1e-5) const
    {T scale=maxabs(eigenvalues.x11,eigenvalues.x33),scale_inverse=Robust_Inverse(scale),tiny=tolerance*scale;
    if(eigenvalues.x22-eigenvalues.x33>tiny) return ((*this-eigenvalues.x33)*scale_inverse).Cofactor_Matrix().Largest_Column_Normalized();
    return ((*this-eigenvalues.x11)*scale_inverse).Cofactor_Matrix().Largest_Column().Unit_Orthogonal_Vector();}

    SYMMETRIC_MATRIX Positive_Definite_Part() const
    {DIAGONAL_MATRIX<T,3> D;MATRIX<T,3> V;Fast_Solve_Eigenproblem(D,V);D=D.Clamp_Min(0);return Conjugate(V,D);}

    DIAGONAL_MATRIX<T,3> Diagonal_Part() const
    {return DIAGONAL_MATRIX<T,3>(x11,x22,x33);}

    static SYMMETRIC_MATRIX Multiply_With_Symmetric_Result(const MATRIX<T,3>& A,const MATRIX<T,3>& B) // A*B and assume symmetric result, 18 mults, 12 adds
    {return SYMMETRIC_MATRIX(A.x[0]*B.x[0]+A.x[3]*B.x[1]+A.x[6]*B.x[2],A.x[1]*B.x[0]+A.x[4]*B.x[1]+A.x[7]*B.x[2],
                             A.x[2]*B.x[0]+A.x[5]*B.x[1]+A.x[8]*B.x[2],A.x[1]*B.x[3]+A.x[4]*B.x[4]+A.x[7]*B.x[5],
                             A.x[2]*B.x[3]+A.x[5]*B.x[4]+A.x[8]*B.x[5],A.x[2]*B.x[6]+A.x[5]*B.x[7]+A.x[8]*B.x[8]);}

    static SYMMETRIC_MATRIX Times_Transpose_With_Symmetric_Result(const MATRIX<T,3>& A,const MATRIX<T,3>& B) // A*B^t and assume symmetric result, 18 mults, 12 adds
    {return SYMMETRIC_MATRIX(A.x[0]*B.x[0]+A.x[3]*B.x[3]+A.x[6]*B.x[6],A.x[1]*B.x[0]+A.x[4]*B.x[3]+A.x[7]*B.x[6],
                             A.x[2]*B.x[0]+A.x[5]*B.x[3]+A.x[8]*B.x[6],A.x[1]*B.x[1]+A.x[4]*B.x[4]+A.x[7]*B.x[7],
                             A.x[2]*B.x[1]+A.x[5]*B.x[4]+A.x[8]*B.x[7],A.x[2]*B.x[2]+A.x[5]*B.x[5]+A.x[8]*B.x[8]);}

    static SYMMETRIC_MATRIX Times_Transpose_With_Symmetric_Result(const MATRIX<T,3,2>& A,const MATRIX<T,3,2>& B) // A*B^t and assume symmetric result, 12 mults, 6 adds
    {return SYMMETRIC_MATRIX(A.x[0]*B.x[0]+A.x[3]*B.x[3],A.x[1]*B.x[0]+A.x[4]*B.x[3],A.x[2]*B.x[0]+A.x[5]*B.x[3],
                             A.x[1]*B.x[1]+A.x[4]*B.x[4],A.x[2]*B.x[1]+A.x[5]*B.x[4],A.x[2]*B.x[2]+A.x[5]*B.x[5]);}

    static SYMMETRIC_MATRIX Transpose_Times_With_Symmetric_Result(const MATRIX<T,3>& A,const MATRIX<T,3>& B) // A^t*B and assume symmetric result, 18 mults, 12 adds
    {return SYMMETRIC_MATRIX(A.x[0]*B.x[0]+A.x[1]*B.x[1]+A.x[2]*B.x[2],A.x[3]*B.x[0]+A.x[4]*B.x[1]+A.x[5]*B.x[2],A.x[6]*B.x[0]+A.x[7]*B.x[1]+A.x[8]*B.x[2],
                             A.x[3]*B.x[3]+A.x[4]*B.x[4]+A.x[5]*B.x[5],A.x[6]*B.x[3]+A.x[7]*B.x[4]+A.x[8]*B.x[5],A.x[6]*B.x[6]+A.x[7]*B.x[7]+A.x[8]*B.x[8]);}

    static SYMMETRIC_MATRIX Transpose_Times_With_Symmetric_Result(const MATRIX<T,3>& A,const UPPER_TRIANGULAR_MATRIX<T,3>& B) // A^t*B and assume symmetric result, 10 mults, 4 adds
    {return SYMMETRIC_MATRIX(A.x[0]*B.x11,A.x[3]*B.x11,A.x[6]*B.x11,A.x[3]*B.x12+A.x[4]*B.x22,A.x[6]*B.x12+A.x[7]*B.x22,A.x[6]*B.x13+A.x[7]*B.x23+A.x[8]*B.x33);}

    SYMMETRIC_MATRIX<T,3> Conjugate_With_Cross_Product_Matrix(const VECTOR<T,3>& v) const
    {return Cross_Product_Matrix_Times(v).Times_Cross_Product_Matrix_With_Symmetric_Result(-v);}

//#####################################################################
    MATRIX<T,3> operator*(const DIAGONAL_MATRIX<T,3>& A) const;
    MATRIX<T,3> operator*(const UPPER_TRIANGULAR_MATRIX<T,3>& A) const;
    SYMMETRIC_MATRIX operator+(const DIAGONAL_MATRIX<T,3>& A) const;
    static SYMMETRIC_MATRIX Conjugate(const MATRIX<T,3>& A,const DIAGONAL_MATRIX<T,3>& B);
    static SYMMETRIC_MATRIX Conjugate(const MATRIX<T,3>& A,const SYMMETRIC_MATRIX& B);
    static SYMMETRIC_MATRIX Conjugate(const SYMMETRIC_MATRIX<T,3>& A,const SYMMETRIC_MATRIX& B);
    static SYMMETRIC_MATRIX Conjugate(const MATRIX<T,3,2>& A,const DIAGONAL_MATRIX<T,2>& B);
    static SYMMETRIC_MATRIX Conjugate(const MATRIX<T,3,2>& A,const SYMMETRIC_MATRIX<T,2>& B);
    static SYMMETRIC_MATRIX Conjugate_With_Transpose(const MATRIX<T,3>& A,const DIAGONAL_MATRIX<T,3>& B);
    static SYMMETRIC_MATRIX Conjugate_With_Transpose(const MATRIX<T,3>& A,const SYMMETRIC_MATRIX& B);
    static SYMMETRIC_MATRIX Conjugate_With_Transpose(const UPPER_TRIANGULAR_MATRIX<T,3>& A,const SYMMETRIC_MATRIX& B);
    DIAGONAL_MATRIX<T,3> Fast_Eigenvalues() const;
    void Fast_Solve_Eigenproblem(DIAGONAL_MATRIX<T,3>& eigenvalues,MATRIX<T,3>& eigenvectors) const;
    void Solve_Eigenproblem(DIAGONAL_MATRIX<T,3>& eigenvalues,MATRIX<T,3>& eigenvectors) const;
private:
    static void Jacobi_Transform(const int sweep,const T threshold,T& app,T& apq,T& aqq,T& arp,T& arq,T& v1p,T& v1q,T& v2p,T& v2q,T& v3p,T& v3q);
//#####################################################################
};
// global functions
template<class T>
inline SYMMETRIC_MATRIX<T,3> operator*(const T a,const SYMMETRIC_MATRIX<T,3>& A)
{return A*a;}

template<class T>
inline SYMMETRIC_MATRIX<T,3> operator+(const T a,const SYMMETRIC_MATRIX<T,3>& A)
{return A+a;}

template<class T>
inline SYMMETRIC_MATRIX<T,3> operator-(const T a,const SYMMETRIC_MATRIX<T,3>& A)
{return -A+a;}

template<class T>
inline SYMMETRIC_MATRIX<T,3> clamp(const SYMMETRIC_MATRIX<T,3>& x,const SYMMETRIC_MATRIX<T,3>& xmin,const SYMMETRIC_MATRIX<T,3>& xmax)
{return SYMMETRIC_MATRIX<T,3>(clamp(x.x11,xmin.x11,xmax.x11),clamp(x.x21,xmin.x21,xmax.x21),clamp(x.x31,xmin.x31,xmax.x31),clamp(x.x22,xmin.x22,xmax.x22),clamp(x.x32,xmin.x32,xmax.x32),clamp(x.x33,xmin.x33,xmax.x33));}

template<class T>
inline SYMMETRIC_MATRIX<T,3> clamp_min(const SYMMETRIC_MATRIX<T,3>& x,const SYMMETRIC_MATRIX<T,3>& xmin)
{return SYMMETRIC_MATRIX<T,3>(clamp_min(x.x11,xmin.x11),clamp_min(x.x21,xmin.x21),clamp_min(x.x31,xmin.x31),clamp_min(x.x22,xmin.x22),clamp_min(x.x32,xmin.x32),clamp_min(x.x33,xmin.x33));}

template<class T>
inline SYMMETRIC_MATRIX<T,3> clamp_max(const SYMMETRIC_MATRIX<T,3>& x,const SYMMETRIC_MATRIX<T,3>& xmax)
{return SYMMETRIC_MATRIX<T,3>(clamp_max(x.x11,xmax.x11),clamp_max(x.x21,xmax.x21),clamp_max(x.x31,xmax.x31),clamp_max(x.x22,xmax.x22),clamp_max(x.x32,xmax.x32),clamp_max(x.x33,xmax.x33));}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
template<class T>
inline std::ostream& operator<< (std::ostream& output_stream,const SYMMETRIC_MATRIX<T,3>& A)
{output_stream<<"["<<A.x11<<" "<<A.x21<<" "<<A.x31<<" ; "<<A.x21<<" "<<A.x22<<" "<<A.x32<<" ; "<<A.x31<<" "<<A.x32<<" "<<A.x33<<"]";return output_stream;}
#endif

template<class T>
inline SYMMETRIC_MATRIX<T,3> log(const SYMMETRIC_MATRIX<T,3>& A)
{DIAGONAL_MATRIX<T,3> D;MATRIX<T,3> Q;A.Fast_Solve_Eigenproblem(D,Q);return A.Conjugate(Q,log(D));}

template<class T>
inline SYMMETRIC_MATRIX<T,3> exp(const SYMMETRIC_MATRIX<T,3>& A)
{DIAGONAL_MATRIX<T,3> D;MATRIX<T,3> Q;A.Fast_Solve_Eigenproblem(D,Q);return A.Conjugate(Q,exp(D));}

//#####################################################################
// Function Conjugate
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Conjugate(const MATRIX<T,3>& A,const DIAGONAL_MATRIX<T,3>& B) // 27 mults, 12 adds
{
    return Times_Transpose_With_Symmetric_Result(A*B,A);
}
//#####################################################################
// Function Conjugate
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Conjugate(const MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B)
{
    return Times_Transpose_With_Symmetric_Result(A*B,A);
}
//#####################################################################
// Function Conjugate
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Conjugate(const SYMMETRIC_MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B)
{
    return Times_Transpose_With_Symmetric_Result(A*B,A);
}
//#####################################################################
// Function Conjugate
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Conjugate(const MATRIX<T,3,2>& A,const DIAGONAL_MATRIX<T,2>& B)
{
    return Times_Transpose_With_Symmetric_Result(A*B,A);
}
//#####################################################################
// Function Conjugate
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Conjugate(const MATRIX<T,3,2>& A,const SYMMETRIC_MATRIX<T,2>& B)
{
    return Times_Transpose_With_Symmetric_Result(A*B,A);
}
//#####################################################################
// Function Conjugate_With_Transpose
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Conjugate_With_Transpose(const MATRIX<T,3>& A,const DIAGONAL_MATRIX<T,3>& B)
{
    return Transpose_Times_With_Symmetric_Result(B*A,A);
}
//#####################################################################
// Function Conjugate_With_Transpose
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Conjugate_With_Transpose(const MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B)
{
    return Transpose_Times_With_Symmetric_Result(B*A,A);
}
//#####################################################################
// Function Conjugate_With_Transpose
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
Conjugate_With_Transpose(const UPPER_TRIANGULAR_MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B)
{
    return Transpose_Times_With_Symmetric_Result(B*A,A);
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T> inline MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
operator*(const DIAGONAL_MATRIX<T,3>& A) const // 9 mults
{
    return MATRIX<T,3>(x11*A.x11,x21*A.x11,x31*A.x11,x21*A.x22,x22*A.x22,x32*A.x22,x31*A.x33,x32*A.x33,x33*A.x33);
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T> inline MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
operator*(const UPPER_TRIANGULAR_MATRIX<T,3>& A) const // 18 mults, 9 adds
{
    return MATRIX<T,3>(x11*A.x11,x21*A.x11,x31*A.x11,x11*A.x12+x21*A.x22,x21*A.x12+x22*A.x22,x31*A.x12+x32*A.x22,
                       x11*A.x13+x21*A.x23+x31*A.x33,x21*A.x13+x22*A.x23+x32*A.x33,x31*A.x13+x32*A.x23+x33*A.x33);
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T> inline MATRIX<T,3>
operator*(const DIAGONAL_MATRIX<T,3>& D,const SYMMETRIC_MATRIX<T,3>& A) // 9 mults, 
{
    return MATRIX<T,3>(D.x11*A.x11,D.x22*A.x21,D.x33*A.x31,D.x11*A.x21,D.x22*A.x22,D.x33*A.x32,D.x11*A.x31,D.x22*A.x32,D.x33*A.x33);
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T> inline MATRIX<T,3>
operator*(const UPPER_TRIANGULAR_MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B) // 18 mults, 9 adds
{
    return MATRIX<T,3>(A.x11*B.x11+A.x12*B.x21+A.x13*B.x31,A.x22*B.x21+A.x23*B.x31,A.x33*B.x31,A.x11*B.x21+A.x12*B.x22+A.x13*B.x32,
                       A.x22*B.x22+A.x23*B.x32,A.x33*B.x32,A.x11*B.x31+A.x12*B.x32+A.x13*B.x33,A.x22*B.x32+A.x23*B.x33,A.x33*B.x33);
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3> SYMMETRIC_MATRIX<T,3>::
operator+(const DIAGONAL_MATRIX<T,3>& A) const // 3 adds
{
    return SYMMETRIC_MATRIX<T,3>(x11+A.x11,x21,x31,x22+A.x22,x32,x33+A.x33);
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3>
operator+(const DIAGONAL_MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B) // 3 adds
{
    return B+A;
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline MATRIX<T,3>
operator+(const SYMMETRIC_MATRIX<T,3>& A,const UPPER_TRIANGULAR_MATRIX<T,3>& B)
{
    return MATRIX<T,3>(A.x11+B.x11,A.x21,A.x31,A.x21+B.x12,A.x22+B.x22,A.x32,A.x31+B.x13,A.x32+B.x23,A.x33+B.x33);
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> inline MATRIX<T,3>
operator+(const UPPER_TRIANGULAR_MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B)
{
    return B+A;
}
//#####################################################################
// Function operator-
//#####################################################################
template<class T> inline MATRIX<T,3>
operator-(const SYMMETRIC_MATRIX<T,3>& A,const UPPER_TRIANGULAR_MATRIX<T,3>& B)
{
    return MATRIX<T,3>(A.x11-B.x11,A.x21,A.x31,A.x21-B.x12,A.x22-B.x22,A.x32,A.x31-B.x13,A.x32-B.x23,A.x33-B.x33);
}
//#####################################################################
// Function operator-
//#####################################################################
template<class T> inline MATRIX<T,3>
operator-(const UPPER_TRIANGULAR_MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B)
{
    return -B+A;
}
//#####################################################################
// Function operator-
//#####################################################################
template<class T> inline SYMMETRIC_MATRIX<T,3>
operator-(const DIAGONAL_MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& B) // 3 adds
{
    return SYMMETRIC_MATRIX<T,3>(A.x11-B.x11,-B.x21,-B.x31,A.x22-B.x22,-B.x32,A.x33-B.x33);
}
//#####################################################################
}
#endif

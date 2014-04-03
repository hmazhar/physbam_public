//#####################################################################
// Copyright 2007, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX
//#####################################################################
#ifndef __MATRIX__
#define __MATRIX__

#include <PhysBAM_Tools/Matrices/MATRIX_0X0.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T,int m_input,int n_input> // n_input=m_input
class MATRIX:public MATRIX_BASE<T,MATRIX<T,m_input,n_input> >
{
public:
    enum WORKAROUND1 {m=m_input,n=n_input,size=n_input*m_input};
    STATIC_ASSERT((!((m>=2 && m<=3 && n>=2 && n<=3) || (m==4 && n==4) || (m==0 && n==0)))); // 0x0, 1x1, 2x2, 2x3, 3x2, 3x3, and 4x4 are specialized
    typedef T SCALAR;typedef MATRIX_BASE<T,MATRIX<T,m_input,n_input> > BASE;
    using BASE::Frobenius_Norm_Squared;using BASE::operator*;using BASE::Transpose_Times;using BASE::Times_Transpose;

    T x[m*n+(m*n==0)]; // pointer to the one dimensional data

    explicit MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(m),INITIAL_SIZE nn=INITIAL_SIZE(n))
    {
        STATIC_ASSERT(sizeof(MATRIX)==(size+!size)*sizeof(T));assert(mm==INITIAL_SIZE(m) && nn==INITIAL_SIZE(n));
        for(int i=0;i<size;i++) x[i]=T();
    }

    MATRIX(const MATRIX& A)
        :BASE()
    {
        for(int i=0;i<size;i++) x[i]=A.x[i];
    }

    explicit MATRIX(const VECTOR<T,size>& column1)
    {
        STATIC_ASSERT(m==1 || n==1);
        for(int i=0;i<column1.Size();i++) x[i]=column1(i+1);
    }

    template<class T_MATRIX>
    explicit MATRIX(const MATRIX_BASE<T,T_MATRIX>& A)
    {
        assert(m==A.Rows() && n==A.Columns());for(int j=1;j<=n;j++) for(int i=1;i<=m;i++) (*this)(i,j)=A(i,j);
    }

    int Rows() const
    {return m;}

    int Columns() const
    {return n;}

    T& operator()(const int i,const int j)
    {assert(i>=1 && i<=m);assert(j>=1 && j<=n);return x[(j-1)*m+(i-1)];}

    const T& operator()(const int i,const int j) const
    {assert(i>=1 && i<=m);assert(j>=1 && j<=n);return x[(j-1)*m+(i-1)];}

    bool Valid_Index(const int i,const int j) const
    {return 1<=i && i<=m && 1<=j && j<=n;}

    VECTOR<T,m>& Column(const int j)
    {assert(1<=j && j<=n);return *(VECTOR<T,m>*)(x+m*(j-1));}

    const VECTOR<T,m>& Column(const int j) const
    {assert(1<=j && j<=n);return *(const VECTOR<T,n>*)(x+m*(j-1));}

    bool operator==(const MATRIX& A) const
    {for(int i=0;i<size;i++) if(x[i]!=A.x[i]) return false;return true;}

    bool operator!=(const MATRIX& A) const
    {return !(*this==A);}

    MATRIX& operator=(const MATRIX& A)
    {for(int i=0;i<size;i++) x[i]=A.x[i];return *this;}

    template<class T_MATRIX>
    MATRIX& operator=(const MATRIX_BASE<T,T_MATRIX>& A)
    {assert(Rows()==A.Rows() && Columns()==A.Columns());for(int j=1;j<=Columns();j++) for(int i=1;i<=Rows();i++) (*this)(i,j)=A(i,j);return *this;}

    MATRIX& operator*=(const T a)
    {for(int i=0;i<size;i++) x[i]*=a;return *this;}

    MATRIX& operator+=(const MATRIX& A)
    {for(int i=0;i<size;i++) x[i]+=A.x[i];return *this;}

    MATRIX& operator-=(const MATRIX& A)
    {for(int i=0;i<size;i++) x[i]-=A.x[i];return *this;}

    MATRIX operator+(const MATRIX& A) const
    {assert(n==A.n && m==A.m);MATRIX matrix;for(int i=0;i<size;i++) matrix.x[i]=x[i]+A.x[i];return matrix;}

    MATRIX operator-(const MATRIX& A) const
    {assert(n==A.n && m==A.m);MATRIX matrix;for(int i=0;i<size;i++) matrix.x[i]=x[i]-A.x[i];return matrix;}

    MATRIX operator-() const
    {MATRIX matrix;for(int i=0;i<size;i++) matrix.x[i]=-x[i];return matrix;}

    MATRIX operator*(const T a) const
    {MATRIX matrix;for(int i=0;i<size;i++) matrix.x[i]=x[i]*a;return matrix;}

    VECTOR<T,m> operator*(const VECTOR<T,n>& y) const
    {VECTOR<T,m> result;for(int j=1;j<=n;j++) for(int i=1;i<=m;i++) result(i)+=(*this)(i,j)*y(j);return result;}

    template<int p>
    MATRIX<T,m,p> operator*(const MATRIX<T,n,p>& A) const
    {MATRIX<T,m,p> matrix;for(int j=1;j<=p;j++) for(int k=1;k<=n;k++) for(int i=1;i<=m;i++) matrix(i,j)+=(*this)(i,k)*A(k,j);return matrix;}

    MATRIX<T,m,n> operator*(const SYMMETRIC_MATRIX<T,n>& A) const
    {MATRIX<T,m,n> matrix;for(int j=1;j<=n;j++) for(int k=1;k<=n;k++) for(int i=1;i<=m;i++) matrix(i,j)+=(*this)(i,k)*A(k,j);return matrix;}

    MATRIX<T,m,n> operator*(const DIAGONAL_MATRIX<T,n>& A) const
    {MATRIX<T,m,n> matrix;for(int j=1;j<=n;j++) for(int i=1;i<=m;i++) matrix(i,j)=(*this)(i,j)*A(j,j);return matrix;}

    MATRIX_MXN<T> operator*(const MATRIX_MXN<T>& A) const
    {assert(n==A.m);MATRIX_MXN<T> matrix(m,A.n);for(int j=1;j<=A.n;j++) for(int i=1;i<=m;i++) for(int k=1;k<=n;k++) matrix(i,j)+=(*this)(i,k)*A(k,j);return matrix;}

    MATRIX<T,n,m> Transposed() const
    {MATRIX<T,n,m> matrix;for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) matrix(j,i)=(*this)(i,j);return matrix;}

    VECTOR<T,n> Transpose_Times(const VECTOR<T,m>& y) const
    {VECTOR<T,n> result;for(int j=1;j<=n;j++) for(int i=1;i<=m;i++) result(j)+=(*this)(i,j)*y(i);return result;}

    template<int p>
    MATRIX<T,n,p> Transpose_Times(const MATRIX<T,m,p>& A) const
    {MATRIX<T,n,p> matrix;for(int j=1;j<=p;j++) for(int i=1;i<=n;i++) for(int k=1;k<=m;k++) matrix(i,j)+=(*this)(k,i)*A(k,j);return matrix;}

    template<int p>
    MATRIX<T,m,p> Times_Transpose(const MATRIX<T,p,n>& A) const
    {MATRIX<T,m,p> matrix;for(int j=1;j<=p;j++) for(int i=1;i<=m;i++) for(int k=1;k<=n;k++) matrix(i,j)+=(*this)(i,k)*A(j,k);return matrix;}

    MATRIX Times_Cross_Product_Matrix(const VECTOR<T,3>& v) const
    {STATIC_ASSERT(n==3);MATRIX matrix;for(int i=1;i<=m;i++) matrix.Set_Row(i,VECTOR<T,3>::Cross_Product(VECTOR<T,3>((*this)(i,1),(*this)(i,2),(*this)(i,3)),v));return matrix;}

    MATRIX Times_Cross_Product_Matrix_Transpose(const VECTOR<T,3>& v) const
    {STATIC_ASSERT(n==3);MATRIX matrix;for(int i=1;i<=m;i++) matrix.Set_Row(i,VECTOR<T,3>::Cross_Product(v,VECTOR<T,3>((*this)(i,1),(*this)(i,2),(*this)(i,3))));return matrix;}

    MATRIX Cross_Product_Matrix_Times(const VECTOR<T,3>& v) const
    {STATIC_ASSERT(m==3);MATRIX matrix;for(int i=1;i<=n;i++) matrix.Set_Column(i,VECTOR<T,3>::Cross_Product(v,VECTOR<T,3>((*this)(1,i),(*this)(2,i),(*this)(3,i))));return matrix;}

    MATRIX Cross_Product_Matrix_Transpose_Times(const VECTOR<T,3>& v) const
    {STATIC_ASSERT(m==3);MATRIX matrix;for(int i=1;i<=n;i++) matrix.Set_Column(i,VECTOR<T,3>::Cross_Product(VECTOR<T,3>((*this)(1,i),(*this)(2,i),(*this)(3,i)),v));return matrix;}

    MATRIX<T,m,2> Times_Cross_Product_Matrix(const VECTOR<T,2>& v) const
    {STATIC_ASSERT(n==1);return (*this)*MATRIX<T,1,2>::Cross_Product_Matrix(v);}

    MATRIX<T,m,1> Times_Cross_Product_Matrix_Transpose(const VECTOR<T,2>& v) const
    {STATIC_ASSERT(n==2);return Times_Transpose(MATRIX<T,1,2>::Cross_Product_Matrix(v));}

    MATRIX<T,1,n> Cross_Product_Matrix_Times(const VECTOR<T,2>& v) const
    {STATIC_ASSERT(m==2);return MATRIX<T,1,2>::Cross_Product_Matrix(v)*(*this);}

    MATRIX<T,2,n> Cross_Product_Matrix_Transpose_Times(const VECTOR<T,2>& v) const
    {STATIC_ASSERT(m==1);return MATRIX<T,1,2>::Cross_Product_Matrix(v).Transpose_Times(*this);}

    static MATRIX<T,1,2> Cross_Product_Matrix(const VECTOR<T,2>& v)
    {STATIC_ASSERT(m==1 && n==2);MATRIX<T,1,2> M;M(1,1)=-v.y;M(1,2)=v.x;return M;}

    MATRIX<T,1> Times_Cross_Product_Matrix_Transpose_With_Symmetric_Result(const VECTOR<T,2>& v) const
    {STATIC_ASSERT((m==1 && n==2));return Times_Cross_Product_Matrix_Transpose(v);}

    MATRIX<T,m,1> Times_Cross_Product_Matrix(const VECTOR<T,1>& v) const
    {STATIC_ASSERT(n==0);return MATRIX<T,m,1>();}

    MATRIX<T,m,0> Times_Cross_Product_Matrix_Transpose(const VECTOR<T,1>& v) const
    {STATIC_ASSERT(n==1);return MATRIX<T,m,0>();}

    MATRIX<T,0,n> Cross_Product_Matrix_Times(const VECTOR<T,1>& v) const
    {STATIC_ASSERT(m==1);return MATRIX<T,0,n>();}

    MATRIX<T,1,n> Cross_Product_Matrix_Transpose_Times(const VECTOR<T,1>& v) const
    {STATIC_ASSERT(m==0);return MATRIX<T,1,n>();}

    static MATRIX<T,0,1> Cross_Product_Matrix(const VECTOR<T,1>& v)
    {STATIC_ASSERT(m==0 && n==1);return MATRIX<T,0,1>();}

    MATRIX Permute_Columns(const VECTOR<int,n>& p) const
    {MATRIX x;for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) x(i,j)=(*this)(i,p(j));return x;}

    MATRIX Unpermute_Columns(const VECTOR<int,n>& p) const
    {MATRIX x;for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) x(i,p(j))=(*this)(i,j);return x;}

    static MATRIX Outer_Product(const VECTOR<T,m>& u,const VECTOR<T,n>& v)
    {MATRIX result;for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) result(i,j)=u(i)*v(j);return result;}

    MATRIX<T,n> Normal_Equations_Matrix() const
    {MATRIX<T,n> result;for(int j=1;j<=n;j++) for(int i=1;i<=n;i++) for(int k=1;k<=m;k++) result(i,j)+=(*this)(k,i)*(*this)(k,j);return result;}

    VECTOR<T,n> Normal_Equations_Solve(const VECTOR<T,m>& b) const
    {MATRIX<T,n> A_transpose_A(Normal_Equations_Matrix());VECTOR<T,n> A_transpose_b(Transpose_Times(b));return A_transpose_A.Cholesky_Solve(A_transpose_b);}

    template<class T_VECTOR>
    VECTOR<T,n> Solve_Linear_System(const VECTOR_BASE<T,T_VECTOR>& b)
    {return PLU_Solve(b);}

    T Parallelepiped_Measure() const
    {STATIC_ASSERT(n==1);return sqrt(Frobenius_Norm_Squared());}
};

template<class T,int m,int n>
inline MATRIX<T,m,n> operator*(const T a,const MATRIX<T,m,n>& A)
{return A*a;}
}
#endif

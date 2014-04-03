//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Huamin Wang, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_MXN
//#####################################################################
#ifndef __MATRIX_MXN__
#define __MATRIX_MXN__

#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Tools/Matrices/MATRIX_BASE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
namespace PhysBAM{

template<class T>
class MATRIX_MXN:public MATRIX_BASE<T,MATRIX_MXN<T> >
{
public:
    typedef MATRIX_BASE<T,MATRIX_MXN<T> > BASE;
    using BASE::operator+;using BASE::operator+=;using BASE::operator-;using BASE::operator-=;using BASE::operator*;using BASE::operator*=;
    using BASE::Times_Transpose;using BASE::Transpose_Times;
    typedef T SCALAR;
    int m,n; // size of the m by n matrix
    T *x; // pointer to the one dimensional data

    MATRIX_MXN()
        :m(0),n(0),x(0)
    {}

    MATRIX_MXN(const int m_input,const int n_input)
        :m(m_input),n(n_input)
    {
        x=new T[m*n];
        for(int k=0;k<m*n;k++) x[k]=0;
    }

    explicit MATRIX_MXN(const int n_input)
        :m(n_input),n(n_input)
    {
        x=new T[m*n];
        for(int k=0;k<m*n;k++) x[k]=0;
    }

    MATRIX_MXN(const INITIAL_SIZE m_input,const INITIAL_SIZE n_input)
        :m(Value(m_input)),n(Value(n_input))
    {
        x=new T[m*n];
        for(int k=0;k<m*n;k++) x[k]=0;
    }

    explicit MATRIX_MXN(const INITIAL_SIZE n_input)
        :m(Value(n_input)),n(Value(n_input))
    {
        x=new T[m*n];
        for(int k=0;k<m*n;k++) x[k]=0;
    }

    MATRIX_MXN(const MATRIX_MXN<T>& A)
        :m(A.m),n(A.n)
    {
        x=new T[m*n];
        for(int k=0;k<m*n;k++) x[k]=A.x[k];
    }

    template<class T2>
    explicit MATRIX_MXN(const MATRIX_MXN<T2>& A)
        :m(A.m),n(A.n)
    {
        x=new T[m*n];
        for(int k=0;k<m*n;k++) x[k]=(T)A.x[k];
    }

    template<class T_MATRIX>
    explicit MATRIX_MXN(const MATRIX_BASE<T,T_MATRIX>& A)
        :m(A.Rows()),n(A.Columns())
    {x=new T[m*n];for(int j=1;j<=n;j++) for(int i=1;i<=m;i++) (*this)(i,j)=A(i,j);}

    template<int d>
    explicit MATRIX_MXN(const DIAGONAL_MATRIX<T,d>& A)
        :m(d),n(d)
    {x=new T[m*n];for(int k=0;k<m*n;k++) x[k]=0;for(int i=1;i<=d;i++) (*this)(i,i)=A(i,i);}

    template<int d>
    explicit MATRIX_MXN(const SYMMETRIC_MATRIX<T,d>& A)
        :m(A.Rows()),n(A.Columns())
    {x=new T[m*n];for(int j=1;j<=n;j++) for(int i=1;i<=j;i++) (*this)(i,j)=(*this)(j,i)=A.Element_Lower(j,i);}

    template<int d>
    explicit MATRIX_MXN(const UPPER_TRIANGULAR_MATRIX<T,d>& A)
        :m(A.Rows()),n(A.Columns())
    {x=new T[m*n];
    for(int j=1;j<=n;j++){
        for(int i=1;i<=j;i++) (*this)(i,j)=A(i,j);
        for(int i=j+1;i<=m;i++) (*this)(i,j)=0;}}

    ~MATRIX_MXN()
    {delete[] x;}

    int Rows() const
    {return m;}

    int Columns() const
    {return n;}

    void Resize(const int m_new,const int n_new)
    {if(m_new==m && n_new==n) return;
    T* x_new=new T[m_new*n_new];for(int t=0;t<m_new*n_new;t++) x_new[t]=(T)0;
    int m1=min(m,m_new),n1=min(n,n_new);for(int i=1;i<=m1;i++) for(int j=1;j<=n1;j++) x_new[(j-1)*m_new+(i-1)]=(*this)(i,j);
    delete[] x;x=x_new;m=m_new;n=n_new;}

    T& operator()(const int i,const int j)
    {assert(i>=1 && i<=m);assert(j>=1 && j<=n);return x[(j-1)*m+(i-1)];}

    const T& operator()(const int i,const int j) const
    {assert(i>=1 && i<=m);assert(j>=1 && j<=n);return x[(j-1)*m+(i-1)];}

    bool Valid_Index(const int i,const int j) const
    {return 1<=i && i<=m && 1<=j && j<=n;}

    bool operator==(const MATRIX_MXN& A) const
    {if(m!=A.m || n!=A.n) return false;
    for(int i=0;i<m*n;i++) if(x[i]!=A.x[i]) return false;return true;}

    bool operator!=(const MATRIX_MXN& A) const
    {return !(*this==A);}

    MATRIX_MXN<T>& operator=(const MATRIX_MXN<T>& A)
    {if(&A==this) return *this;
    if(!x || m*n!=A.m*A.n){delete[] x;x=new T[A.m*A.n];};
    m=A.m;n=A.n;for(int k=0;k<m*n;k++) x[k]=A.x[k];return *this;}

    template<class T_MATRIX>
    MATRIX_MXN<T>& operator=(const MATRIX_BASE<T,T_MATRIX>& A)
    {if((void*)&A==(void*)this) return *this;
    if(!x || m*n!=A.Rows()*A.Columns()){delete[] x;x=new T[A.Rows()*A.Columns()];}
    m=A.Rows();n=A.Columns();for(int j=1;j<=n;j++) for(int i=1;i<=m;i++) (*this)(i,j)=A(i,j);return *this;}

    template<int d>
    MATRIX_MXN<T>& operator=(const DIAGONAL_MATRIX<T,d>& A)
    {if(!x || m*n!=A.Rows()*A.Columns()){delete[] x;x=new T[d*d];}
    m=n=d;for(int k=0;k<m*n;k++) x[k]=0;for(int i=1;i<=d;i++) (*this)(i,i)=A(i,i);return *this;}

    template<int d>
    MATRIX_MXN<T>& operator=(const SYMMETRIC_MATRIX<T,d>& A)
    {if(!x || m*n!=A.Rows()*A.Columns()){delete[] x;x=new T[d*d];}
    m=n=d;for(int j=1;j<=n;j++) for(int i=1;i<=j;i++) (*this)(i,j)=(*this)(j,i)=A.Element_Lower(j,i);return *this;}

    T Trace() const
    {assert(m==n);T trace=0;for(int i=1;i<=n;i++) trace+=(*this)(i,i);return trace;}

    void Transpose()
    {*this=Transposed();}

    MATRIX_MXN<T> Transposed() const
    {MATRIX_MXN<T> matrix(n,m);for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) matrix(j,i)=(*this)(i,j);return matrix;}

    MATRIX_MXN<T> Times_Cross_Product_Matrix_Transpose(const VECTOR<T,2>& v) const
    {assert(n==2);return (*this)*MATRIX<T,1,2>::Cross_Product_Matrix(v).Transposed();}

    MATRIX_MXN<T> Cross_Product_Matrix_Times(const VECTOR<T,2>& v) const
    {assert(m==2);return MATRIX<T,1,2>::Cross_Product_Matrix(v)*(*this);}

    MATRIX_MXN<T> Times_Cross_Product_Matrix_Transpose(const VECTOR<T,3>& v) const
    {assert(n==3);MATRIX_MXN<T> matrix(m,3);for(int i=1;i<=m;i++) matrix.Set_Row(i,VECTOR<T,3>::Cross_Product(v,VECTOR<T,3>((*this)(i,1),(*this)(i,2),(*this)(i,3))));return matrix;}

    MATRIX_MXN<T> Cross_Product_Matrix_Times(const VECTOR<T,3>& v) const
    {assert(m==3);MATRIX_MXN<T> matrix(3,n);for(int i=1;i<=n;i++) matrix.Set_Column(i,VECTOR<T,3>::Cross_Product(v,VECTOR<T,3>((*this)(1,i),(*this)(2,i),(*this)(3,i))));return matrix;}

    MATRIX_MXN<T> Permute_Columns(const VECTOR_ND<int>& p) const
    {assert(n==p.n);MATRIX_MXN<T> x(m,n);for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) x(i,j)=(*this)(i,p(j));return x;}

    MATRIX_MXN<T> Unpermute_Columns(const VECTOR_ND<int>& p) const
    {assert(n==p.n);MATRIX_MXN<T> x(m,n);for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) x(i,p(j))=(*this)(i,j);return x;}

    static MATRIX_MXN<T> Outer_Product(const VECTOR_ND<T> u,const VECTOR_ND<T> v)
    {MATRIX_MXN<T> result(u.n,v.n);for(int i=1;i<=u.n;i++) for(int j=1;j<=v.n;j++) result(i,j)=u(i)*v(j);return result;}

    MATRIX_MXN<T> Normal_Equations_Matrix() const
    {MATRIX_MXN<T> result(n);for(int j=1;j<=n;j++) for(int i=j;i<=n;i++){T a=0;for(int k=1;k<=m;k++) a+=(*this)(k,i)*(*this)(k,j);result(i,j)=result(j,i)=a;}return result;}

    VECTOR_ND<T> Normal_Equations_Solve(const VECTOR_ND<T>& b) const
    {MATRIX_MXN<T> A_transpose_A(Normal_Equations_Matrix());VECTOR_ND<T> A_transpose_b(Transpose_Times(b));return A_transpose_A.Cholesky_Solve(A_transpose_b);}

    void Gauss_Seidel_Single_Iteration(VECTOR_ND<T>& x,const VECTOR_ND<T>& b) const
    {assert(m==n && x.n==b.n && x.n==n);
    for(int i=1;i<=n;i++){
        T rho=0;
        for(int j=1;j<i;j++) rho+=(*this)(i,j)*x(j);
        for(int j=i+1;j<=n;j++) rho+=(*this)(i,j)*x(j);
        x(i)=(b(i)-rho)/(*this)(i,i);}}

    void Left_Givens_Rotation(const int i,const int j,const T c,const T s)
    {assert(1<=i && i<j && j<=m);for(int k=1;k<=n;k++){T x=(*this)(i,k);(*this)(i,k)=c*(*this)(i,k)-s*(*this)(j,k);(*this)(j,k)=s*x+c*(*this)(j,k);}}
    
    void Right_Givens_Rotation(const int i,const int j,const T c,const T s)
    {assert(1<=i && i<j && j<=n);for(int k=1;k<=m;k++){T x=(*this)(k,i);(*this)(k,i)=c*(*this)(k,i)-s*(*this)(k,j);(*this)(k,j)=s*x+c*(*this)(k,j);}}

//#####################################################################
    void Jacobi_Singular_Value_Decomposition(ARRAY<VECTOR<int,2> >& left_givens_pairs,ARRAY<VECTOR<T,2> >& left_givens_coefficients,
        ARRAY<VECTOR<int,2> >& right_givens_pairs,ARRAY<VECTOR<T,2> >& right_givens_coefficients,const T tolerance=(T)1e-10,
        const int max_iterations=1000000);
//#####################################################################
};
template<class T>
inline MATRIX_MXN<T> operator*(const T a,const MATRIX_MXN<T>& A)
{return A*a;}
}
#endif

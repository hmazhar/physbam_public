//#####################################################################
// Copyright 2007, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Huamin Wang, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_BASE
//#####################################################################
#ifndef __MATRIX_BASE__
#define __MATRIX_BASE__

#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Matrices/MATRIX_ARITHMETIC_POLICY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_BASE.h>
#include <iomanip>
#include <iostream>
namespace PhysBAM{

template<class T_MATRIX> struct MATRIX_INFO;
template<class T,int m,int n> struct MATRIX_INFO<MATRIX<T,m,n> >{typedef VECTOR<T,m> LEFT_VECTOR;typedef VECTOR<T,n> RIGHT_VECTOR;typedef MATRIX<T,n,m> TRANSPOSE;};
template<class T,int d> struct MATRIX_INFO<DIAGONAL_MATRIX<T,d> >{typedef VECTOR<T,d> LEFT_VECTOR;typedef VECTOR<T,d> RIGHT_VECTOR;typedef DIAGONAL_MATRIX<T,d> TRANSPOSE;};
template<class T,int d> struct MATRIX_INFO<SYMMETRIC_MATRIX<T,d> >{typedef VECTOR<T,d> LEFT_VECTOR;typedef VECTOR<T,d> RIGHT_VECTOR;typedef SYMMETRIC_MATRIX<T,d> TRANSPOSE;};
template<class T,int d> struct MATRIX_INFO<UPPER_TRIANGULAR_MATRIX<T,d> >{typedef VECTOR<T,d> LEFT_VECTOR;typedef VECTOR<T,d> RIGHT_VECTOR;};
template<class T> struct MATRIX_INFO<MATRIX_MXN<T> >{typedef VECTOR_ND<T> LEFT_VECTOR;typedef VECTOR_ND<T> RIGHT_VECTOR;typedef MATRIX_MXN<T> TRANSPOSE;};
template<class T,class T_MATRIX> struct MATRIX_INFO<MATRIX_BASE<T,T_MATRIX> >{typedef typename MATRIX_INFO<T_MATRIX>::LEFT_VECTOR LEFT_VECTOR;
    typedef typename MATRIX_INFO<T_MATRIX>::RIGHT_VECTOR RIGHT_VECTOR;typedef typename MATRIX_INFO<T_MATRIX>::TRANSPOSE TRANSPOSE;};

template<class T_MATRIX> struct EFFICIENT_MATRIX {static const bool value=false;};
template<class T,int m,int n> struct EFFICIENT_MATRIX<MATRIX<T,m,n> > {static const bool value=((m>=2 && m<=3 && n>=2 && n<=3) || (m==4 && n==4) || (m==0 && n==0));};
template<class T,int d> struct EFFICIENT_MATRIX<DIAGONAL_MATRIX<T,d> > {static const bool value=(d==2 || d==3);};
template<class T,int d> struct EFFICIENT_MATRIX<SYMMETRIC_MATRIX<T,d> > {static const bool value=(d==2 || d==3);};
template<class T,class T_MATRIX> struct EFFICIENT_MATRIX<MATRIX_BASE<T,T_MATRIX> >:public EFFICIENT_MATRIX<T_MATRIX>{};
template<class T,int d> struct EFFICIENT_MATRIX<VECTOR<T,d> > {static const bool value=(d<=3);};
template<class T,class T_VECTOR> struct EFFICIENT_MATRIX<VECTOR_BASE<T,T_VECTOR> >:public EFFICIENT_MATRIX<T_VECTOR>{};

namespace{
template<int line,class A,class B=void> struct ASSERT_EFFICIENT
{
    template<class T> struct EFFICIENT_OR_VOID:public OR<EFFICIENT_MATRIX<T>::value,IS_SAME<T,void>::value>{};
    static const bool efficient=EFFICIENT_OR_VOID<A>::value && EFFICIENT_OR_VOID<B>::value;
    struct UNUSABLE{};

    ASSERT_EFFICIENT(typename IF<efficient,UNUSABLE,const char*>::TYPE str)
    {}

    ASSERT_EFFICIENT(typename IF<efficient,const char*,UNUSABLE>::TYPE str)
    {PHYSBAM_WARNING(std::string("Base implementation used for: ")+str+".");}
};

#ifdef NDEBUG
#define WARN_IF_NOT_EFFICIENT(...)
#else
#ifdef WIN32
#define WARN_IF_NOT_EFFICIENT(...) ASSERT_EFFICIENT<__LINE__,__VA_ARGS__>(__FUNCTION__)
#else
#define WARN_IF_NOT_EFFICIENT(...) ASSERT_EFFICIENT<__LINE__,__VA_ARGS__>(__PRETTY_FUNCTION__)
#endif
#endif
}

template<class T,class T_MATRIX>
class MATRIX_BASE
{
    template<class T_MATRIX2>
    T_MATRIX& operator=(const T_MATRIX&) const
    {STATIC_ASSERT((T)false);}

    MATRIX_BASE& operator=(const MATRIX_BASE&) const
    {STATIC_ASSERT((T)false);}

    MATRIX_BASE(const MATRIX_BASE&)
    {STATIC_ASSERT((T)false);}

public:
    typedef T SCALAR;
    typedef typename MATRIX_INFO<T_MATRIX>::RIGHT_VECTOR RIGHT_VECTOR;typedef typename MATRIX_INFO<T_MATRIX>::LEFT_VECTOR LEFT_VECTOR;
    typedef typename RIGHT_VECTOR::template REBIND<int>::TYPE COLUMN_PERMUTATION;

    MATRIX_BASE()
    {}

    ~MATRIX_BASE()
    {}

    T_MATRIX& Derived()
    {return static_cast<T_MATRIX&>(*this);}

    const T_MATRIX& Derived() const
    {return static_cast<const T_MATRIX&>(*this);}

    T operator()(const int i,const int j) const
    {return Derived()(i,j);}

    T& operator()(const int i,const int j)
    {return Derived()(i,j);}

    int Rows() const
    {return Derived().Rows();}

    int Columns() const
    {return Derived().Columns();}

protected:
    // Logic only; no type or bounds checking; requires all arguments to be distinct
    template<class T_MATRIX1,class T_MATRIX2,class T_MATRIX3>
    static void Add_Times_Matrix_Helper(const T_MATRIX1& A,const T_MATRIX2& B,T_MATRIX3& C)
    {for(int j=1;j<=B.Columns();j++) for(int k=1;k<=A.Columns();k++) for(int i=1;i<=A.Rows();i++) C(i,j)+=A(i,k)*B(k,j);}

    template<class T_MATRIX1,class T_MATRIX2,class T_MATRIX3>
    static void Add_Transpose_Times_Matrix_Helper(const T_MATRIX1& A,const T_MATRIX2& B,T_MATRIX3& C)
    {for(int j=1;j<=B.Columns();j++) for(int i=1;i<=A.Columns();i++) for(int k=1;k<=A.Rows();k++) C(i,j)+=A(k,i)*B(k,j);}

    template<class T_MATRIX1,class T_MATRIX2,class T_MATRIX3>
    static void Add_Times_Transpose_Matrix_Helper(const T_MATRIX1& A,const T_MATRIX2& B,T_MATRIX3& C)
    {for(int k=1;k<=A.Columns();k++) for(int j=1;j<=B.Rows();j++) for(int i=1;i<=A.Rows();i++) C(i,j)+=A(i,k)*B(j,k);}

    template<class T_MATRIX1,class T_VECTOR,class T_VECTOR2>
    static void Add_Transpose_Times_Vector_Helper(const T_MATRIX1& A,const VECTOR_BASE<T,T_VECTOR>& v,VECTOR_BASE<T,T_VECTOR2>& y)
    {for(int i=1;i<=y.Size();i++) for(int j=1;j<=v.Size();j++) y(i)+=A(j,i)*v(j);}

    template<class T_MATRIX1,class T_VECTOR,class T_VECTOR2>
    static void Add_Times_Vector_Helper(const T_MATRIX1& A,const VECTOR_BASE<T,T_VECTOR>& v,VECTOR_BASE<T,T_VECTOR2>& y)
    {for(int j=1;j<=v.Size();j++) for(int i=1;i<=y.Size();i++) y(i)+=A(i,j)*v(j);}

    template<class T_MATRIX1,class T_VECTOR,class T_VECTOR2>
    static void Subtract_Times_Vector_Helper(const T_MATRIX1& A,const VECTOR_BASE<T,T_VECTOR>& v,VECTOR_BASE<T,T_VECTOR2>& y)
    {for(int j=1;j<=v.Size();j++) for(int i=1;i<=y.Size();i++) y(i)-=A(i,j)*v(j);}

    template<class T_MATRIX1,class T_MATRIX2>
    static bool Test_Aliased(const T_MATRIX1& A,const T_MATRIX2& B)
    {assert((void*)&A!=(void*)&B);return false;}

    template<class T_MATRIX1>
    static bool Test_Aliased(const T_MATRIX1& A,const T_MATRIX1& B)
    {return &A==&B;}
public:

    // With bounds checking; arguments may be the same.
    template<class T_MATRIX1,class T_MATRIX2,class T_MATRIX3>
    static void Add_Times(const T_MATRIX1& A,const T_MATRIX2& B,MATRIX_BASE<T,T_MATRIX3>& C)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX1,T_MATRIX2);assert(A.Columns()==B.Rows() && C.Rows()==A.Rows() && C.Columns()==B.Columns());
    if(Test_Aliased(A,C) || Test_Aliased(B,C)){T_MATRIX3 M(INITIAL_SIZE(C.Rows()),INITIAL_SIZE(C.Columns()));
    Add_Times_Matrix_Helper(A,B,M);C+=M;}else Add_Times_Matrix_Helper(A,B,C);}

    template<class T_MATRIX1,class T_MATRIX2,class T_MATRIX3>
    static void Add_Transpose_Times(const T_MATRIX1& A,const T_MATRIX2& B,MATRIX_BASE<T,T_MATRIX3>& C)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX1,T_MATRIX2);assert(A.Rows()==B.Rows() && C.Rows()==A.Columns() && C.Columns()==B.Columns());
    if(Test_Aliased(A,C) || Test_Aliased(B,C)){T_MATRIX3 M(INITIAL_SIZE(C.Rows()),INITIAL_SIZE(C.Columns()));
    Add_Transpose_Times_Matrix_Helper(A,B,M);C+=M;}else Add_Transpose_Times_Matrix_Helper(A,B,C);}

    template<class T_MATRIX1,class T_MATRIX2,class T_MATRIX3>
    static void Add_Times_Transpose(const T_MATRIX1& A,const T_MATRIX2& B,MATRIX_BASE<T,T_MATRIX3>& C)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX1,T_MATRIX2);assert(A.Columns()==B.Columns() && C.Rows()==A.Rows() && C.Columns()==B.Rows());
    if(Test_Aliased(A,C) || Test_Aliased(B,C)){T_MATRIX3 M(INITIAL_SIZE(C.Rows()),INITIAL_SIZE(C.Columns()));
    Add_Times_Transpose_Matrix_Helper(A,B,M);C+=M;}else Add_Times_Transpose_Matrix_Helper(A,B,C);}

    template<class T_MATRIX1,class T_VECTOR,class T_VECTOR2>
    static void Add_Transpose_Times(const MATRIX_BASE<T,T_MATRIX1>& A,const VECTOR_BASE<T,T_VECTOR>& v,VECTOR_BASE<T,T_VECTOR2>& y)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX1,T_VECTOR2);assert(A.Rows()==v.Size() && A.Columns()==y.Size());
    if(Test_Aliased(v,y)){T_VECTOR u(v);Add_Transpose_Times_Vector_Helper(A,u,y);}else Add_Transpose_Times_Vector_Helper(A,v,y);}

    template<class T_MATRIX1,class T_VECTOR,class T_VECTOR2>
    static void Add_Times(const MATRIX_BASE<T,T_MATRIX1>& A,const VECTOR_BASE<T,T_VECTOR>& v,VECTOR_BASE<T,T_VECTOR2>& y)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX1,T_VECTOR);assert(A.Columns()==v.Size() && A.Rows()==y.Size());
    if(Test_Aliased(v,y)){T_VECTOR u(v);Add_Times_Vector_Helper(A,u,y);}else Add_Times_Vector_Helper(A,v,y);}

    template<class T_MATRIX1,class T_VECTOR,class T_VECTOR2>
    static void Subtract_Times(const MATRIX_BASE<T,T_MATRIX1>& A,const VECTOR_BASE<T,T_VECTOR>& v,VECTOR_BASE<T,T_VECTOR2>& y)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX1,T_VECTOR);assert(A.Columns()==v.Size() && A.Rows()==y.Size());
    if(Test_Aliased(v,y)){T_VECTOR u(v);Subtract_Times_Vector_Helper(A,u,y);}else Subtract_Times_Vector_Helper(A,v,y);}

    template<class T_MATRIX1,class T_VECTOR,class T_VECTOR2>
    static void Add_Transpose_Times(const MATRIX_BASE<T,T_MATRIX1>& A,const VECTOR_EXPRESSION<T,T_VECTOR>& v,VECTOR_BASE<T,T_VECTOR2>& y)
    {assert(A.Rows()==v.Size() && A.Columns()==y.Size());Add_Transpose_Times_Vector_Helper(A,typename VECTOR_TYPE<T_VECTOR>::TYPE(v),y);}

    template<class T_MATRIX1,class T_VECTOR,class T_VECTOR2>
    static void Add_Times(const MATRIX_BASE<T,T_MATRIX1>& A,const VECTOR_EXPRESSION<T,T_VECTOR>& v,VECTOR_BASE<T,T_VECTOR2>& y)
    {assert(A.Columns()==v.Size() && A.Rows()==y.Size());Add_Times_Vector_Helper(A,typename VECTOR_TYPE<T_VECTOR>::TYPE(v),y);}

    template<class T_MATRIX1,class T_VECTOR,class T_VECTOR2>
    static void Subtract_Times(const MATRIX_BASE<T,T_MATRIX1>& A,const VECTOR_EXPRESSION<T,T_VECTOR>& v,VECTOR_BASE<T,T_VECTOR2>& y)
    {assert(A.Columns()==v.Size() && A.Rows()==y.Size());Subtract_Times_Vector_Helper(A,typename VECTOR_TYPE<T_VECTOR>::TYPE(v),y);}

    template<class T_MATRIX1>
    typename TRANSPOSE_PRODUCT<T_MATRIX,T_MATRIX1>::TYPE Transpose_Times(const MATRIX_BASE<T,T_MATRIX1>& A) const
    {typename TRANSPOSE_PRODUCT<T_MATRIX,T_MATRIX1>::TYPE M((INITIAL_SIZE)Columns(),(INITIAL_SIZE)A.Columns());Add_Transpose_Times(Derived(),A,M);return M;}

    template<int d>
    typename TRANSPOSE_PRODUCT<T_MATRIX,SYMMETRIC_MATRIX<T,d> >::TYPE Transpose_Times(const SYMMETRIC_MATRIX<T,d>& A) const
    {typename TRANSPOSE_PRODUCT<T_MATRIX,SYMMETRIC_MATRIX<T,d> >::TYPE M((INITIAL_SIZE)Columns(),(INITIAL_SIZE)d);Add_Transpose_Times(Derived(),A,M);return M;}

    template<int d>
    typename TRANSPOSE_PRODUCT<T_MATRIX,DIAGONAL_MATRIX<T,d> >::TYPE Transpose_Times(const DIAGONAL_MATRIX<T,d>& A) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);typename TRANSPOSE_PRODUCT<T_MATRIX,DIAGONAL_MATRIX<T,d> >::TYPE M((INITIAL_SIZE)Columns(),(INITIAL_SIZE)d);
    for(int i=1;i<=Columns();i++) for(int k=1;k<=Rows();k++) M(i,k)=(*this)(k,i)*A(k,k);return M;}

    template<class T_MATRIX1>
    typename PRODUCT_TRANSPOSE<T_MATRIX,T_MATRIX1>::TYPE Times_Transpose(const MATRIX_BASE<T,T_MATRIX1>& A) const
    {typename PRODUCT_TRANSPOSE<T_MATRIX,T_MATRIX1>::TYPE M((INITIAL_SIZE)Rows(),(INITIAL_SIZE)A.Rows());Add_Times_Transpose(Derived(),A,M);return M;}

    template<int d>
    typename PRODUCT<T_MATRIX,DIAGONAL_MATRIX<T,d> >::TYPE Times_Transpose(const DIAGONAL_MATRIX<T,d>& A) const
    {return Derived()*A;}

    template<int d>
    typename PRODUCT<T_MATRIX,SYMMETRIC_MATRIX<T,d> >::TYPE Times_Transpose(const SYMMETRIC_MATRIX<T,d>& A) const
    {return Derived()*A;}

    template<class T_VECTOR>
    typename TRANSPOSE_PRODUCT<T_MATRIX,typename VECTOR_TYPE<T_VECTOR>::TYPE>::TYPE Transpose_Times(const VECTOR_BASE<T,T_VECTOR>& y) const
    {assert(y.Size()==Rows());typename TRANSPOSE_PRODUCT<T_MATRIX,typename VECTOR_TYPE<T_VECTOR>::TYPE>::TYPE result((INITIAL_SIZE)Columns());
    Add_Transpose_Times(Derived(),y.Derived(),result);return result;}

    template<class T_MATRIX1>
    typename PRODUCT<T_MATRIX,T_MATRIX1>::TYPE operator*(const MATRIX_BASE<T,T_MATRIX1>& A) const
    {assert(Columns()==A.Rows());typename PRODUCT<T_MATRIX,T_MATRIX1>::TYPE M((INITIAL_SIZE)Rows(),(INITIAL_SIZE)A.Columns());Add_Times(Derived(),A,M);return M;}

    template<class T_VECTOR>
    typename PRODUCT<T_MATRIX,typename VECTOR_TYPE<T_VECTOR>::TYPE>::TYPE operator*(const VECTOR_BASE<T,T_VECTOR>& y) const
    {assert(y.Size()==Columns());typename PRODUCT<T_MATRIX,typename VECTOR_TYPE<T_VECTOR>::TYPE>::TYPE result((INITIAL_SIZE)Rows());
    Add_Times(Derived(),y.Derived(),result);return result;}

    template<int d>
    typename PRODUCT<T_MATRIX,DIAGONAL_MATRIX<T,d> >::TYPE operator*(const DIAGONAL_MATRIX<T,d>& A) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);typename PRODUCT<T_MATRIX,DIAGONAL_MATRIX<T,d> >::TYPE M((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    for(int k=1;k<=Columns();k++) for(int i=1;i<=Rows();i++) M(i,k)=(*this)(i,k)*A(k,k);return M;}

    template<int d>
    typename PRODUCT<T_MATRIX,SYMMETRIC_MATRIX<T,d> >::TYPE operator*(const SYMMETRIC_MATRIX<T,d>& A) const
    {assert(Columns()==A.Rows());typename PRODUCT<T_MATRIX,SYMMETRIC_MATRIX<T,d> >::TYPE M((INITIAL_SIZE)Rows(),(INITIAL_SIZE)A.Columns());Add_Times(Derived(),A,M);return M;}

    template<int d>
    T_MATRIX operator*(const UPPER_TRIANGULAR_MATRIX<T,d>& A) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);assert(Columns()==d);
    T_MATRIX M((INITIAL_SIZE)Rows(),(INITIAL_SIZE)d);
    for(int j=1;j<=d;j++) for(int k=1;k<=j;k++) for(int i=1;i<=Rows();i++) M(i,j)+=(*this)(i,k)*A(k,j);return M;}

    T_MATRIX& operator*=(const T a)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);for(int j=1;j<=Columns();j++) for(int i=1;i<=Rows();i++) (*this)(i,j)*=a;return Derived();}

    T_MATRIX& operator/=(const T a)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);return Derived()*=(1/a);}

    T_MATRIX operator*(const T a) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);T_MATRIX matrix((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    for(int j=1;j<=Columns();j++) for(int i=1;i<=Rows();i++) matrix(i,j)=(*this)(i,j)*a;return matrix;}

    T_MATRIX operator/(const T a) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);return Derived()*(1/a);}

    T_MATRIX& operator+=(const T a)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);assert(Rows()==Columns());for(int i=1;i<=Rows();i++) (*this)(i,i)+=a;return Derived();}

    T_MATRIX operator+(const T a) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);return T_MATRIX(Derived())+=a;}

    template<class T_MATRIX1>
    T_MATRIX& operator+=(const MATRIX_BASE<T,T_MATRIX1>& A)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX,T_MATRIX1);assert(Rows()==A.Rows() && Columns()==A.Columns());
    for(int j=1;j<=Columns();j++) for(int i=1;i<=Rows();i++) (*this)(i,j)+=A(i,j);return Derived();}

    template<int d>
    T_MATRIX& operator+=(const SYMMETRIC_MATRIX<T,d>& A)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX,SYMMETRIC_MATRIX<T,d>);assert(Rows()==A.Rows() && Columns()==A.Columns());
    for(int j=1;j<=Columns();j++){for(int i=j+1;i<=Rows();i++){T element=A.Element_Lower(i,j);(*this)(i,j)+=element;(*this)(j,i)+=element;}(*this)(j,j)+=A.Element_Lower(j,j);}
    return Derived();}

    template<class T_MATRIX1>
    T_MATRIX& operator-=(const MATRIX_BASE<T,T_MATRIX1>& A)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX,T_MATRIX1);assert(Rows()==A.Rows() && Columns()==A.Columns());
    for(int j=1;j<=Columns();j++) for(int i=1;i<=Rows();i++) (*this)(i,j)-=A(i,j);return Derived();}

    template<class T_MATRIX1>
    T_MATRIX operator+(const MATRIX_BASE<T,T_MATRIX1>& A) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX,T_MATRIX1);assert(Rows()==A.Rows() && Columns()==A.Columns());T_MATRIX matrix((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    for(int j=1;j<=Columns();j++) for(int i=1;i<=Rows();i++) matrix(i,j)=(*this)(i,j)+A(i,j);return matrix;}

    template<class T_MATRIX1>
    T_MATRIX operator-(const MATRIX_BASE<T,T_MATRIX1>& A) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX,T_MATRIX1);assert(Rows()==A.Rows() && Columns()==A.Columns());T_MATRIX matrix((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    for(int j=1;j<=Columns();j++) for(int i=1;i<=Rows();i++) matrix(i,j)=(*this)(i,j)-A(i,j);return matrix;}

    T_MATRIX operator-() const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);T_MATRIX matrix((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    for(int j=1;j<=Columns();j++) for(int i=1;i<=Rows();i++) matrix(i,j)=-(*this)(i,j);return matrix;}

    T Max_Abs() const
    {T max_abs=0;for(int j=1;j<=Columns();j++) for(int i=1;i<=Rows();i++) max_abs=max(max_abs,abs((*this)(i,j)));return max_abs;}

    T Infinity_Norm() const
    {T max_sum=0;for(int j=1;j<=Columns();j++){T sum=0;for(int i=1;i<=Rows();i++) sum+=abs((*this)(i,j));max_sum=max(sum,max_sum);}return max_sum;}

    T Frobenius_Norm() const
    {return sqrt(Derived().Frobenius_Norm_Squared());}

    T Frobenius_Norm_Squared() const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);T sum=0;for(int j=1;j<=Columns();j++) for(int i=1;i<=Rows();i++) sum+=sqr((*this)(i,j));return sum;}

    template<class T_MATRIX2>
    void Set_Submatrix(const int istart,const int jstart,const MATRIX_BASE<T,T_MATRIX2>& a)
    {for(int i=1;i<=a.Rows();i++) for(int j=1;j<=a.Columns();j++) (*this)(istart+i-1,jstart+j-1)=a(i,j);}

    void Set_Submatrix(const int istart,const int jstart,const SYMMETRIC_MATRIX<T,3>& a)
    {(*this)(istart,jstart)=a.x11;(*this)(istart,jstart+1)=a.x21;(*this)(istart+1,jstart)=a.x21;(*this)(istart,jstart+2)=a.x31;(*this)(istart+2,jstart)=a.x31;
    (*this)(istart+1,jstart+1)=a.x22;(*this)(istart+1,jstart+2)=a.x32;(*this)(istart+2,jstart+1)=a.x32;(*this)(istart+2,jstart+2)=a.x33;}

    void Set_Submatrix(const int istart,const int jstart,const DIAGONAL_MATRIX<T,3>& a)
    {(*this)(istart,jstart)=a.x11;(*this)(istart+1,jstart+1)=a.x22;(*this)(istart+2,jstart+2)=a.x33;}

    template<class T_VECTOR>
    void Set_Submatrix(const int istart,const int jstart,const VECTOR_BASE<T,T_VECTOR>& a)
    {for(int i=1;i<=a.Size();i++) (*this)(istart+i-1,jstart)=a(i);}

    void Set_Submatrix(const int istart,const int jstart,const SYMMETRIC_MATRIX<T,2>& a)
    {(*this)(istart,jstart)=a.x11;(*this)(istart,jstart+1)=a.x21;(*this)(istart+1,jstart)=a.x21;(*this)(istart+1,jstart+1)=a.x22;}

    template<class T_MATRIX2>
    void Add_To_Submatrix(const int istart,const int jstart,const MATRIX_BASE<T,T_MATRIX2>& a)
    {for(int i=1;i<=a.Rows();i++) for(int j=1;j<=a.Columns();j++) (*this)(istart+i-1,jstart+j-1)+=a(i,j);}

    void Add_To_Submatrix(const int istart,const int jstart,const SYMMETRIC_MATRIX<T,3>& a)
    {(*this)(istart,jstart)+=a.x11;(*this)(istart,jstart+1)+=a.x21;(*this)(istart+1,jstart)+=a.x21;(*this)(istart,jstart+2)+=a.x31;(*this)(istart+2,jstart)+=a.x31;
    (*this)(istart+1,jstart+1)+=a.x22;(*this)(istart+1,jstart+2)+=a.x32;(*this)(istart+2,jstart+1)+=a.x32;(*this)(istart+2,jstart+2)+=a.x33;}

    void Add_To_Submatrix(const int istart,const int jstart,const DIAGONAL_MATRIX<T,3>& a)
    {(*this)(istart,jstart)+=a.x11;(*this)(istart+1,jstart+1)+=a.x22;(*this)(istart+2,jstart+2)+=a.x33;}

    void Add_To_Submatrix(const int istart,const int jstart,const DIAGONAL_MATRIX<T,2>& a)
    {(*this)(istart,jstart)+=a.x11;(*this)(istart+1,jstart+1)+=a.x22;}

    template<class T_VECTOR>
    void Add_To_Submatrix(const int istart,const int jstart,const VECTOR_BASE<T,T_VECTOR>& a)
    {for(int i=1;i<=a.Size();i++) (*this)(istart+i-1,jstart)+=a(i);}

    template<class T_MATRIX2>
    void Subtract_From_Submatrix(const int istart,const int jstart,const MATRIX_BASE<T,T_MATRIX2>& a)
    {for(int i=1;i<=a.Rows();i++) for(int j=1;j<=a.Columns();j++) (*this)(istart+i-1,jstart+j-1)-=a(i,j);}

    template<class T_MATRIX2>
    void Get_Submatrix(const int istart,const int jstart,MATRIX_BASE<T,T_MATRIX2>& a) const
    {for(int i=1;i<=a.Rows();i++) for(int j=1;j<=a.Columns();j++) a(i,j)=(*this)(istart+i-1,jstart+j-1);}

    void Get_Submatrix(const int istart,const int jstart,SYMMETRIC_MATRIX<T,3>& a) const
    {a.x11=(*this)(istart,jstart);a.x21=(*this)(istart,jstart+1);a.x21=(*this)(istart+1,jstart);a.x31=(*this)(istart,jstart+2);a.x31=(*this)(istart+2,jstart);
    a.x22=(*this)(istart+1,jstart+1);a.x32=(*this)(istart+1,jstart+2);a.x32=(*this)(istart+2,jstart+1);a.x33=(*this)(istart+2,jstart+2);}

    void Get_Submatrix(const int istart,const int jstart,DIAGONAL_MATRIX<T,3>& a) const
    {a.x11=(*this)(istart,jstart);a.x22=(*this)(istart+1,jstart+1);a.x33=(*this)(istart+2,jstart+2);}

    template<class T_VECTOR>
    void Set_Column(const int j,const VECTOR_BASE<T,T_VECTOR>& a)
    {for(int i=1;i<=Rows();i++) (*this)(i,j)=a(i);}

    template<class T_VECTOR>
    void Get_Column(const int j,VECTOR_BASE<T,T_VECTOR>& a) const
    {for(int i=1;i<=Rows();i++) a(i)=(*this)(i,j);}

    template<class T_VECTOR>
    void Set_Row(const int i,const VECTOR_BASE<T,T_VECTOR>& a)
    {for(int j=1;j<=Columns();j++) (*this)(i,j)=a(j);}

    template<class T_VECTOR>
    void Get_Row(const int i,VECTOR_BASE<T,T_VECTOR>& a) const
    {for(int j=1;j<=Columns();j++) a(j)=(*this)(i,j);}

    void Set_Zero_Matrix()
    {for(int i=1;i<=Rows();i++) for(int j=1;j<=Columns();j++) (*this)(i,j)=0;}

    void Add_Identity_Matrix()
    {assert(Rows()==Columns());for(int i=1;i<=Columns();i++) (*this)(i,i)+=1;}

    void Set_Identity_Matrix()
    {Set_Zero_Matrix();Add_Identity_Matrix();}

    static T_MATRIX Identity_Matrix(const int n)
    {T_MATRIX A(n);A.Add_Identity_Matrix();return A;}

    template<class T_VECTOR>
    T Symmetric_Conjugate(const VECTOR_BASE<T,T_VECTOR>& v) const
    {assert(Rows()==Columns());T r=0;for(int j=1;j<=Columns();j++){T a=0;for(int i=1;i<j;i++) a+=v(i)*(*this)(i,j);r+=(a+a+v(j)*(*this)(j,j))*v(j);}return r;}

    template<class T_VECTOR>
    RIGHT_VECTOR Lower_Triangular_Solve(const VECTOR_BASE<T,T_VECTOR>& b) const
    {assert(Rows()==Columns() && Columns()==b.Size());RIGHT_VECTOR x(INITIAL_SIZE(b.Size()));
    for(int i=1;i<=Columns();i++){x(i)=b(i);for(int j=1;j<=i-1;j++) x(i)-=(*this)(i,j)*x(j);x(i)/=(*this)(i,i);}
    return x;}

    template<class T_VECTOR>
    LEFT_VECTOR Transpose_Lower_Triangular_Solve(const VECTOR_BASE<T,T_VECTOR>& b) const
    {assert(Rows()==Columns() && Columns()==b.Size());LEFT_VECTOR x(b);
    for(int i=1;i<=Columns();i++){x(i)/=(*this)(i,i);for(int j=i+1;j<=Columns();j++) x(j)-=(*this)(i,j)*x(i);}
    return x;}

    template<class T_VECTOR>
    RIGHT_VECTOR Upper_Triangular_Solve(const VECTOR_BASE<T,T_VECTOR>& b) const
    {assert(Rows()==Columns() && Columns()==b.Size());RIGHT_VECTOR x(INITIAL_SIZE(b.Size()));
    for(int i=Columns();i>=1;i--){x(i)=b(i);for(int j=Columns();j>=i+1;j--) x(i)-=(*this)(i,j)*x(j);x(i)/=(*this)(i,i);}
    return x;}

    template<class T_VECTOR>
    LEFT_VECTOR Transpose_Upper_Triangular_Solve(const VECTOR_BASE<T,T_VECTOR>& b) const
    {assert(Rows()==Columns() && Columns()==b.Size());LEFT_VECTOR x(b);
    for(int i=Columns();i>=1;i--){x(i)/=(*this)(i,i);for(int j=i-1;j>=1;j--) x(j)-=(*this)(i,j)*x(i);}
    return x;}

    template<class T_MATRIX2>
    T_MATRIX2 Upper_Triangular_Solve(const MATRIX_BASE<T,T_MATRIX2>& b) const
    {assert(Rows()==Columns() && Columns()==b.Rows());T_MATRIX2 x(INITIAL_SIZE(b.Rows()),INITIAL_SIZE(b.Columns()));
    for(int bcol=1;bcol<=b.Columns();bcol++) for(int i=Columns();i>=1;i--){
        x(i,bcol)=b(i,bcol);for(int j=Columns();j>=i+1;j--) x(i,bcol)-=(*this)(i,j)*x(j,bcol);x(i,bcol)/=(*this)(i,i);}
    return x;}

    template<class T_VECTOR>
    RIGHT_VECTOR In_Place_Cholesky_Solve(const VECTOR_BASE<T,T_VECTOR>& b)
    {assert(Rows()==Columns());In_Place_Cholesky_Factorization();return Transpose_Upper_Triangular_Solve(Lower_Triangular_Solve(b));}

    template<class T_VECTOR>
    RIGHT_VECTOR Cholesky_Solve(const VECTOR_BASE<T,T_VECTOR>& b) const
    {return T_MATRIX(Derived()).In_Place_Cholesky_Solve(b);}

    void In_Place_Cholesky_Inverse()
    {return T_MATRIX(Derived()).In_Place_Cholesky_Inverse(*this);}

    template<class T_MATRIX2>
    void In_Place_Cholesky_Inverse(MATRIX_BASE<T,T_MATRIX2>& inverse)
    {assert(Rows()==Columns());inverse.Derived()=T_MATRIX2((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());LEFT_VECTOR b((INITIAL_SIZE)Rows()); // holds piece of the identity matrix
    In_Place_Cholesky_Factorization();
    for(int j=1;j<=Columns();j++){b(j)=1;inverse.Set_Column(j,Transpose_Upper_Triangular_Solve(Lower_Triangular_Solve(b)));b(j)=0;}}

    template<class T_MATRIX2>
    void Cholesky_Inverse(MATRIX_BASE<T,T_MATRIX2>& inverse) const
    {return T_MATRIX(Derived()).In_Place_Cholesky_Inverse(inverse);}

    template<class T_VECTOR>
    RIGHT_VECTOR In_Place_Gram_Schmidt_QR_Solve(const VECTOR_BASE<T,T_VECTOR>& b)
    {T_MATRIX R;In_Place_Gram_Schmidt_QR_Factorization(R);return R.Upper_Triangular_Solve(Derived().Transpose_Times(b.Derived()));}

    template<class T_VECTOR>
    RIGHT_VECTOR Gram_Schmidt_QR_Solve(const VECTOR_BASE<T,T_VECTOR>& b) const
    {return T_MATRIX(Derived()).In_Place_Gram_Schmidt_QR_Solve(b);}

    template<class T_VECTOR,class T_MATRIX2>
    static VECTOR_ND<T> Householder_Transform(const VECTOR_BASE<T,T_VECTOR>& b,const MATRIX_BASE<T,T_MATRIX2>& V) // TODO: don't assume VECTOR_ND
    {assert(V.Rows()==b.Size());VECTOR_ND<T> result(b);
    for(int j=1;j<=V.Columns();j++){VECTOR_ND<T> v(V.Rows());for(int i=1;i<=V.Rows();i++) v(i)=V(i,j);result=result.Householder_Transform(v);}
    return result;}

    template<class T_VECTOR>
    RIGHT_VECTOR Householder_QR_Solve(const VECTOR_BASE<T,T_VECTOR>& b)
    {T_MATRIX V,R;Householder_QR_Factorization(V,R);VECTOR_ND<T> c=Householder_Transform(b,V),c_short(Columns());for(int i=1;i<=Columns();i++) c_short(i)=c(i); // TODO: don't assume VECTOR_ND
    return R.Upper_Triangular_Solve(c_short);}

    T Condition_Number() const
    {assert(Rows()==Columns());T_MATRIX inverse;PLU_Inverse(inverse);return Infinity_Norm()*inverse.Infinity_Norm();}

    template<class T_VECTOR>
    RIGHT_VECTOR In_Place_PLU_Solve(const VECTOR_BASE<T,T_VECTOR>& b)
    {assert(Rows()==Columns());COLUMN_PERMUTATION p;T_MATRIX L;In_Place_PLU_Factorization(L,p);
    return Upper_Triangular_Solve(L.Lower_Triangular_Solve(typename VECTOR_TYPE<T_VECTOR>::TYPE(b).Permute(p)));}

    template<class T_VECTOR>
    RIGHT_VECTOR PLU_Solve(const VECTOR_BASE<T,T_VECTOR>& b) const
    {return T_MATRIX(Derived()).In_Place_PLU_Solve(b);}

    template<class T_MATRIX2>
    void In_Place_PLU_Inverse(MATRIX_BASE<T,T_MATRIX2>& inverse)
    {assert(Rows()==Columns());inverse.Derived()=T_MATRIX2((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());COLUMN_PERMUTATION p;T_MATRIX L;
    In_Place_PLU_Factorization(L,p);
    RIGHT_VECTOR b((INITIAL_SIZE)Columns()); // used for piece of the identity matrix
    for(int j=1;j<=Columns();j++){b(j)=0;inverse.Set_Column(j,Upper_Triangular_Solve(L.Lower_Triangular_Solve(b.Permute(p))));b(j)=0;}}

    template<class T_MATRIX2>
    void PLU_Inverse(MATRIX_BASE<T,T_MATRIX2>& inverse) const
    {(T_MATRIX2(Derived())).In_Place_PLU_Inverse(inverse);}

    template<class T_VECTOR>
    RIGHT_VECTOR In_Place_LU_Solve(const VECTOR_BASE<T,T_VECTOR>& b)
    {assert(Rows()==Columns());T_MATRIX L;In_Place_LU_Factorization(L);return Upper_Triangular_Solve(L.Lower_Triangular_Solve(b));}

    template<class T_VECTOR>
    RIGHT_VECTOR LU_Solve(const VECTOR_BASE<T,T_VECTOR>& b) const
    {return T_MATRIX(Derived()).In_Place_LU_Inverse(b);}

    template<class T_MATRIX2>
    void In_Place_LU_Inverse(MATRIX_BASE<T,T_MATRIX2>& inverse) // don't assume VECTOR_ND
    {assert(Rows()==Columns());inverse.Derived()=T_MATRIX2((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());T_MATRIX L;VECTOR_ND<T> b(Columns()); // used forpiece of the identity matrix
    In_Place_LU_Factorization(L);
    for(int j=1;j<=Columns();j++){b(j)=1;inverse.Set_Column(j,Upper_Triangular_Solve(L.Lower_Triangular_Solve(b)));b(j)=0;}}

    template<class T_MATRIX2>
    void LU_Inverse(MATRIX_BASE<T,T_MATRIX2>& inverse) const
    {(T_MATRIX2(Derived())).In_Place_LU_Inverse(inverse);}

//#####################################################################
    template<class T_MATRIX2> void In_Place_Gram_Schmidt_QR_Factorization(MATRIX_BASE<T,T_MATRIX2>& R); // this=Q
    template<class T_MATRIX2,class T_MATRIX3> void Householder_QR_Factorization(MATRIX_BASE<T,T_MATRIX2>& V,MATRIX_BASE<T,T_MATRIX3>& R);
    template<class T_VECTOR1,class T_VECTOR2> void In_Place_Robust_Householder_QR_Solve(VECTOR_BASE<T,T_VECTOR1>& b,VECTOR_BASE<int,T_VECTOR2>& p); // this=Q
    template<class T_MATRIX2> void In_Place_PLU_Factorization(MATRIX_BASE<T,T_MATRIX2>& L,COLUMN_PERMUTATION& p); // this=U
    void In_Place_Cholesky_Factorization(); // this=L
    template<class T_MATRIX2> void In_Place_LU_Factorization(MATRIX_BASE<T,T_MATRIX2>& L); // this=U
    int Number_Of_Nonzero_Rows(const T threshold) const;
//#####################################################################
};
template<class T,class T_MATRIX> inline std::ostream& operator<<(std::ostream& output_stream,const MATRIX_BASE<T,T_MATRIX>& A)
{output_stream<<"[";for(int i=1;i<=A.Rows();i++){for(int j=1;j<=A.Columns();j++){output_stream<<A(i,j);if(j<A.Columns()) output_stream<<" ";}if(i<A.Rows()) output_stream<<"; ";}output_stream<<"]";return output_stream;}

template<class T,class T_MATRIX,int d>
typename PRODUCT<DIAGONAL_MATRIX<T,d>,T_MATRIX>::TYPE operator*(const DIAGONAL_MATRIX<T,d>& A,const MATRIX_BASE<T,T_MATRIX>& B)
{assert(d==B.Rows());typename PRODUCT<DIAGONAL_MATRIX<T,d>,T_MATRIX>::TYPE M((INITIAL_SIZE)B.Rows(),(INITIAL_SIZE)B.Columns());
for(int i=1;i<=B.Rows();i++){T a=A(i,i);for(int k=1;k<=B.Columns();k++) M(i,k)=a*B(i,k);}return M;}

template<class T,class T_MATRIX,int d>
typename PRODUCT<SYMMETRIC_MATRIX<T,d>,T_MATRIX>::TYPE operator*(const SYMMETRIC_MATRIX<T,d>& A,const MATRIX_BASE<T,T_MATRIX>& B)
{typename PRODUCT<SYMMETRIC_MATRIX<T,d>,T_MATRIX>::TYPE M((INITIAL_SIZE)B.Rows(),(INITIAL_SIZE)B.Columns());B.Add_Times(A,B.Derived(),M);return M;}
}
#endif

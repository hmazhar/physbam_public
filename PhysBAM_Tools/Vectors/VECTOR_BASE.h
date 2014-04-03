//#####################################################################
// Copyright 2007, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VECTOR_BASE
//#####################################################################
#ifndef __VECTOR_BASE__
#define __VECTOR_BASE__

#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/exchange.h>
#include <PhysBAM_Tools/Math_Tools/Inverse.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Vectors/ARITHMETIC_POLICY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T_VECTOR> struct VECTOR_TYPE;
template<class T> struct VECTOR_TYPE<VECTOR_ND<T> > {typedef VECTOR_ND<T> TYPE;};
template<class T,int d> struct VECTOR_TYPE<VECTOR<T,d> > {typedef VECTOR<T,d> TYPE;};
template<class T,class T_VECTOR> struct VECTOR_TYPE<VECTOR_EXPRESSION<T,T_VECTOR> > {typedef typename VECTOR_TYPE<T_VECTOR>::TYPE TYPE;};
template<class T_VECTOR,class T_VECTOR2> struct VECTOR_TYPE<VECTOR_SUM<T_VECTOR,T_VECTOR2> > {typedef typename VECTOR_TYPE<T_VECTOR>::TYPE TYPE;};
template<class T_VECTOR,class T_VECTOR2> struct VECTOR_TYPE<VECTOR_DIFFERENCE<T_VECTOR,T_VECTOR2> > {typedef typename VECTOR_TYPE<T_VECTOR>::TYPE TYPE;};
template<class T_VECTOR> struct VECTOR_TYPE<VECTOR_NEGATION<T_VECTOR> > {typedef typename VECTOR_TYPE<T_VECTOR>::TYPE TYPE;};
template<class T2,class T_VECTOR> struct VECTOR_TYPE<VECTOR_SCALE<T2,T_VECTOR> > {typedef typename VECTOR_TYPE<T_VECTOR>::TYPE TYPE;};

template<class TV> struct EFFICIENT_VECTOR{static const bool value=false;};
template<class T> struct EFFICIENT_VECTOR<VECTOR<T,0> >{static const bool value=true;};
template<class T> struct EFFICIENT_VECTOR<VECTOR<T,1> >{static const bool value=true;};
template<class T> struct EFFICIENT_VECTOR<VECTOR<T,2> >{static const bool value=true;};
template<class T> struct EFFICIENT_VECTOR<VECTOR<T,3> >{static const bool value=true;};
template<class TV> struct INEFFICIENT_VECTOR{static const bool value=IS_VECTOR<TV>::value && !EFFICIENT_VECTOR<TV>::value;};

template<class T,class T_VECTOR>
class VECTOR_BASE
{
    VECTOR_BASE& operator=(const VECTOR_BASE&) const
    {STATIC_ASSERT((T)false);}

    VECTOR_BASE(const VECTOR_BASE&)
    {STATIC_ASSERT((T)false);}

public:

    VECTOR_BASE()
    {}

    ~VECTOR_BASE()
    {}

    T_VECTOR& Derived()
    {return static_cast<T_VECTOR&>(*this);}

    const T_VECTOR& Derived() const
    {return static_cast<const T_VECTOR&>(*this);}

    T operator()(const int i) const
    {return Derived()(i);}

    T& operator()(const int i)
    {return Derived()(i);}

    int Size() const
    {return Derived().Size();}

    template<int d>
    static void Static_Assert_Not_Small(const VECTOR<T,d>*)
    {STATIC_ASSERT(d>3);}

    template<class T_VECTOR2>
    static void Static_Assert_Not_Small(const T_VECTOR2*)
    {}

    void Static_Assert_Not_Small() const
    {Static_Assert_Not_Small((T_VECTOR*)0);}

    template<class T2,class T3,int p,int q>
    static void Assert_Same_Size_Helper(const VECTOR<T2,p>&,const VECTOR<T3,q>&)
    {STATIC_ASSERT(p==q);}

    template<class T2,class T3,class T_VECTOR2,class T_VECTOR3>
    static void Assert_Same_Size_Helper(const VECTOR_BASE<T2,T_VECTOR2>& u,const VECTOR_BASE<T3,T_VECTOR3>& v)
    {assert(u.Size()==v.Size());}

    template<class T2,class T3,class T_VECTOR2,class T_VECTOR3>
    static void Assert_Same_Size(const VECTOR_BASE<T2,T_VECTOR2>& u,const VECTOR_BASE<T3,T_VECTOR3>& v)
    {Assert_Same_Size_Helper(u.Derived(),v.Derived());}

    void Set_Zero()
    {for(int i=1;i<=Size();i++) (*this)(i)=T();}

    void Fill(const T& constant)
    {Static_Assert_Not_Small();for(int i=1;i<=Size();i++) (*this)(i)=constant;}

    void Negate()
    {for(int i=1;i<=Size();i++) (*this)(i)=-(*this)(i);}

    T Magnitude_Squared() const
    {Static_Assert_Not_Small();T norm_squared=0;for(int i=1;i<=Size();i++) norm_squared+=sqr((*this)(i));return norm_squared;}

    T Magnitude() const
    {Static_Assert_Not_Small();return sqrt(Magnitude_Squared());}

    T Lp_Norm(const T& p) const
    {T sum=0;for(int i=1;i<=Size();i++) sum+=pow((*this)(i),p);return pow(sum,1/p);}

    T L1_Norm() const
    {T sum=0;for(int i=1;i<=Size();i++) sum+=abs((*this)(i));return sum;}

    template<class T_VECTOR1,class T_VECTOR2>
    static T Dot_Product(const VECTOR_BASE<T,T_VECTOR1>& v1,const VECTOR_BASE<T,T_VECTOR2>& v2)
    {v1.Static_Assert_Not_Small();Assert_Same_Size(v1,v2);T sum=0;for(int i=1;i<=v1.Size();i++) sum+=v1(i)*v2(i);return sum;}

    template<class T_VECTOR1,class T_VECTOR2>
    static double Dot_Product_Double_Precision(const VECTOR_BASE<T,T_VECTOR1>& v1,const VECTOR_BASE<T,T_VECTOR2>& v2)
    {v1.Static_Assert_Not_Small();Assert_Same_Size(v1,v2);double sum=0;for(int i=1;i<=v1.Size();i++) sum+=(double)v1(i)*(double)v2(i);return sum;}

    template<class T_VECTOR1,class T_VECTOR2>
    static double Dot_Product_Double_Precision(const VECTOR_BASE<T,T_VECTOR1>& v1,const VECTOR_BASE<T,T_VECTOR2>& v2,const int start_index,const int end_index)
    {v1.Static_Assert_Not_Small();Assert_Same_Size(v1,v2);double sum=0;for(int i=start_index;i<=end_index;i++) sum+=(double)v1(i)*(double)v2(i);return sum;}

    T Sum() const
    {Static_Assert_Not_Small();T result=0;for(int i=1;i<=Size();i++) result+=(*this)(i);return result;}

    T Average() const
    {Static_Assert_Not_Small();return Sum()/Size();}

    T Product() const
    {Static_Assert_Not_Small();T result=1;for(int i=1;i<=Size();i++) result*=(*this)(i);return result;}

    double Sum_Double_Precision() const
    {Static_Assert_Not_Small();double result=0;for(int i=1;i<=Size();i++) result+=(*this)(i);return result;}
    
    double Sum_Double_Precision(int start_index,int end_index) const
    {Static_Assert_Not_Small();double result=0;for(int i=start_index;i<=end_index;i++) result+=(*this)(i);return result;}

    T Maximum_Magnitude() const
    {T result=0;for(int i=1;i<=Size();i++) result=PhysBAM::max(result,abs((*this)(i)));return result;}

    template<class T2,class T_VECTOR1,class T_VECTOR2>
    static void Copy(const T2 c,const T_VECTOR1& v,T_VECTOR2& result)
    {result=c*v;}

    template<class T2,class T_VECTOR1,class T_VECTOR2,class T_VECTOR3>
    static void Copy(const T2 c1,const T_VECTOR1& v1,const T_VECTOR2& v2,T_VECTOR3& result)
    {result=c1*v1+v2;}

    template<class T_VECTOR2>
    void Set_Subvector(const int istart,const VECTOR_BASE<T,T_VECTOR2>& v)
    {for(int i=1;i<=v.Size();i++) (*this)(istart+i-1)=v(i);}

    template<class T_VECTOR2>
    void Add_Subvector(const int istart,const VECTOR_BASE<T,T_VECTOR2>& v)
    {for(int i=1;i<=v.Size();i++) (*this)(istart+i-1)+=v(i);}
    
    template<class T_VECTOR2>
    void Get_Subvector(const int istart,VECTOR_BASE<T,T_VECTOR2>& v) const
    {for(int i=1;i<=v.Size();i++) v(i)=(*this)(istart+i-1);}

    T Sup_Norm() const
    {return Max_Abs();}

    T Max_Abs() const
    {T result=0;for(int i=1;i<=Size();i++) result=PhysBAM::max(result,abs((*this)(i)));return result;}

    T Max_Abs(int start_index,int end_index) const
    {T result=0;for(int i=start_index;i<=end_index;i++) result=PhysBAM::max(result,abs((*this)(i)));return result;}

    T Normalize()
    {Static_Assert_Not_Small();T magnitude=Magnitude();if(magnitude) Derived()*=1/magnitude;else{Set_Zero();(*this)(1)=(T)1;}return magnitude;}

    T_VECTOR Normalized() const
    {Static_Assert_Not_Small();T magnitude=Magnitude();if(magnitude) return Derived()*(1/magnitude);T_VECTOR v((INITIAL_SIZE(Size())));v(1)=(T)1;return v;}

    template<class T_VECTOR1,class T_VECTOR2>
    static T_VECTOR Interpolate(const T_VECTOR1& p1,const T_VECTOR2& p2,const T alpha)
    {Assert_Same_Size(p1,p2);T_VECTOR x((INITIAL_SIZE)p1.Size());for(int i=1;i<=p1.Size();i++) x(i)=(1-alpha)*p1(i)+alpha*p2(i);return x;}

    template<class T_VECTOR2>
    T_VECTOR Permute(const VECTOR_BASE<int,T_VECTOR2>& p) const
    {Assert_Same_Size(*this,p);T_VECTOR x((INITIAL_SIZE)Size());for(int i=1;i<=Size();i++) x(i)=(*this)(p(i));return x;}

    template<class T_VECTOR2>
    T_VECTOR Unpermute(const VECTOR_BASE<int,T_VECTOR2>& p) const
    {Assert_Same_Size(*this,p);T_VECTOR x((INITIAL_SIZE)Size());for(int i=1;i<=Size();i++) x(p(i))=(*this)(i);return x;}

    T_VECTOR Householder_Vector(const int k) const
    {T_VECTOR v((INITIAL_SIZE)Size());T v_dot_v=0;for(int i=k;i<=Size();i++){v(i)=(*this)(i);v_dot_v+=sqr(v(i));}
    if((*this)(k)>=0) v(k)+=sqrt(v_dot_v);else v(k)-=sqrt(v_dot_v);
    return v;}

    template<class T_VECTOR2>
    VECTOR_DIFFERENCE<T_VECTOR,VECTOR_SCALE<T,T_VECTOR2> > Householder_Transform(const VECTOR_BASE<T,T_VECTOR2>& v) const
    {Assert_Same_Size(*this,v);
    T v_dot_a=0,v_dot_v=0;for(int i=1;i<=Size();i++){v_dot_a+=v(i)*(*this)(i);v_dot_v+=sqr(v(i));}
    return *this-2*v_dot_a/v_dot_v*v;}

    void Givens_Rotate(const int i,const int j,const T c,const T s)
    {if(i>j){Givens_Rotate(j,i,c,-s);return;}assert(1<=i && i<j && j<=Size());T u=(*this)(i),v=(*this)(j);(*this)(i)=c*u-s*v;(*this)(j)=s*u+c*v;}

    void Set_To_Orthogonal_Vector() // result isn't normalized
    {assert(Size()>=2);int m1=1,m2=2;
    if(abs((*this)(m1))>abs((*this)(m2))) exchange(m1,m2);
    for(int i=3;i<=Size();i++) if(abs((*this)(i))>abs((*this)(m1))){
        m1=i;if(abs((*this)(m1))>abs((*this)(m2))) exchange(m1,m2);}
    T x1=(*this)(m1),x2=(*this)(m2);
    Set_Zero();
    (*this)(m2)=x1;(*this)(m1)=-x2;}

    T Min() const
    {T result=(*this)(1);for(int i=2;i<=Size();i++) result=min(result,(*this)(i));return result;}

    T Max() const
    {T result=(*this)(1);for(int i=2;i<=Size();i++) result=max(result,(*this)(i));return result;}

    template<class T_VECTOR1,class T_VECTOR2>
    static T Angle_Between(const VECTOR_BASE<T,T_VECTOR1>& u,const VECTOR_BASE<T,T_VECTOR2>& v) // 0 .. pi
    {Assert_Same_Size(u,v);T u2=0,u1=u(1),v1=v(1),uv=0;for(int i=2;i<=u.Size();i++){T ui=u(i);u2+=sqr(ui);uv+=ui*v(i);}u1+=sign_nonzero(u1)*sqrt(u2+sqr(u1));
    T factor=2*(uv+u1*v1)/(u2+sqr(u1)),R12=v1-factor*u1,R22=0;for(int i=2;i<=u.Size();i++) R22+=sqr(v(i)-factor*u(i));return atan2(sqrt(R22),-R12*sign_nonzero(u1));}

    template<class T_VECTOR1,class T_VECTOR2>
    static typename T_VECTOR1::template REBIND<bool>::TYPE Componentwise_Greater_Equal(const VECTOR_BASE<T,T_VECTOR1>& u,const VECTOR_BASE<T,T_VECTOR2>& v)
    {Assert_Same_Size(u,v);typename T_VECTOR1::template REBIND<bool>::TYPE result(INITIAL_SIZE(u.Size()));for(int i=1;i<=u.Size();i++) result(i)=u(i)>=v(i);return result;}

    template<class T_VECTOR1,class T_VECTOR2>
    static T_VECTOR1 Componentwise_And(const VECTOR_BASE<bool,T_VECTOR1>& u,const VECTOR_BASE<bool,T_VECTOR2>& v)
    {Assert_Same_Size(u,v);T_VECTOR1 result(INITIAL_SIZE(u.Size()));for(int i=1;i<=u.Size();i++) result(i)=(u(i) && v(i));return result;}

//#####################################################################
};

template<class T,class T_VECTOR,class T_VECTOR2> T_VECTOR& operator+=(VECTOR_BASE<T,T_VECTOR>& v,const VECTOR_BASE<T,T_VECTOR2>& w)
{v.Static_Assert_Not_Small();v.Assert_Same_Size(v,w);for(int i=1;i<=v.Size();i++) v(i)+=w(i);return v.Derived();}

template<class T,class T_VECTOR,class T_VECTOR2> T_VECTOR& operator-=(VECTOR_BASE<T,T_VECTOR>& v,const VECTOR_BASE<T,T_VECTOR2>& w)
{v.Static_Assert_Not_Small();v.Assert_Same_Size(v,w);for(int i=1;i<=v.Size();i++) v(i)-=w(i);return v.Derived();}

template<class T,class T_VECTOR,class T_VECTOR2> T_VECTOR& operator*=(VECTOR_BASE<T,T_VECTOR>& v,const VECTOR_BASE<T,T_VECTOR2>& w)
{v.Static_Assert_Not_Small();v.Assert_Same_Size(v,w);for(int i=1;i<=v.Size();i++) v(i)*=w(i);return v.Derived();}

template<class T,class T_VECTOR,class T_VECTOR2> T_VECTOR& operator/=(VECTOR_BASE<T,T_VECTOR>& v,const VECTOR_BASE<T,T_VECTOR2>& w)
{v.Static_Assert_Not_Small();v.Assert_Same_Size(v,w);for(int i=1;i<=v.Size();i++) v(i)/=w(i);return v.Derived();}

template<class T,class T_VECTOR> T_VECTOR& operator*=(VECTOR_BASE<T,T_VECTOR>& v,const T& a)
{v.Static_Assert_Not_Small();for(int i=1;i<=v.Size();i++) v(i)*=a;return v.Derived();}

template<class T,class T_VECTOR> T_VECTOR& operator*=(VECTOR_BASE<T,T_VECTOR>& v,const INT_INVERSE& a)
{v.Static_Assert_Not_Small();for(int i=1;i<=v.Size();i++) v(i)*=a;return v.Derived();}

template<class T,class T_VECTOR> T_VECTOR& operator+=(VECTOR_BASE<T,T_VECTOR>& v,const T& a)
{v.Static_Assert_Not_Small();for(int i=1;i<=v.Size();i++) v(i)+=a;return v.Derived();}

template<class T,class T_VECTOR> T_VECTOR& operator-=(VECTOR_BASE<T,T_VECTOR>& v,const T& a)
{v.Static_Assert_Not_Small();for(int i=1;i<=v.Size();i++) v(i)-=a;return v.Derived();}

template<class T,class T_VECTOR> T_VECTOR& operator/=(VECTOR_BASE<T,T_VECTOR>& v,const T& a)
{return v*=Inverse(a);}
//#####################################################################
template<class T_VECTOR1,class T_VECTOR2> struct CAN_ASSIGN<T_VECTOR1,T_VECTOR2,typename ENABLE_IF<IS_VECTOR<T_VECTOR1>::value && IS_VECTOR<T_VECTOR2>::value && CAN_ASSIGN<typename T_VECTOR1::ELEMENT,typename T_VECTOR2::ELEMENT>::value && !IS_SAME<T_VECTOR1,T_VECTOR2>::value>::TYPE>
{static const bool value=true;};
//#####################################################################
}
#include <PhysBAM_Tools/Vectors/VECTOR_EXPRESSION.h>
#endif

//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VECTOR
//#####################################################################
#ifndef __VECTOR__
#define __VECTOR__

#include <PhysBAM_Tools/Math_Tools/Inverse.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <PhysBAM_Tools/Vectors/VECTOR_0D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#ifdef DIFFERENCE // Windows workaround.
#undef DIFFERENCE
#endif
namespace PhysBAM{

template<class T_ARRAY,class T_INDICES> class INDIRECT_ARRAY;

template<class T,int d>
class VECTOR:public VECTOR_BASE<T,VECTOR<T,d> >
{
    STATIC_ASSERT(d>3);
    struct UNUSABLE{};
public:
    typedef VECTOR_BASE<T,VECTOR<T,d> > BASE;
    template<class T2> struct REBIND{typedef VECTOR<T2,d> TYPE;};
    typedef typename IF<IS_SCALAR<T>::value,T,UNUSABLE>::TYPE SCALAR;
    typedef T ELEMENT;
    typedef UNUSABLE SPIN;
    typedef int INDEX;
    typedef T value_type; // for stl
    typedef T* iterator; // for stl
    enum WORKAROUND1 {dimension=d};
    enum WORKAROUND2 {m=d};

    T array[d];

    explicit VECTOR(INITIAL_SIZE n=INITIAL_SIZE(d))
    {
        STATIC_ASSERT(sizeof(VECTOR)==d*sizeof(T));assert(n==INITIAL_SIZE(d));
        for(int i=0;i<d;i++) array[i]=T();
    }

    VECTOR(const T& x1,const T& x2,const T& x3,const T& x4)
    {
        STATIC_ASSERT(d==4);array[0]=x1;array[1]=x2;array[2]=x3;array[3]=x4;
    }

    VECTOR(const T& x1,const T& x2,const T& x3,const T& x4,const T& x5)
    {
        STATIC_ASSERT(d==5);array[0]=x1;array[1]=x2;array[2]=x3;array[3]=x4;array[4]=x5;
    }

    VECTOR(const T& x1,const T& x2,const T& x3,const T& x4,const T& x5,const T& x6)
    {
        STATIC_ASSERT(d==6);array[0]=x1;array[1]=x2;array[2]=x3;array[3]=x4;array[4]=x5;array[5]=x6;
    }

    VECTOR(const T& x1,const T& x2,const T& x3,const T& x4,const T& x5,const T& x6,const T& x7)
    {
        STATIC_ASSERT(d==7);array[0]=x1;array[1]=x2;array[2]=x3;array[3]=x4;array[4]=x5;array[5]=x6;array[6]=x7;
    }

    VECTOR(const T& x1,const T& x2,const T& x3,const T& x4,const T& x5,const T& x6,const T& x7,const T& x8)
    {
        STATIC_ASSERT(d==8);array[0]=x1;array[1]=x2;array[2]=x3;array[3]=x4;array[4]=x5;array[5]=x6;array[6]=x7;array[7]=x8;
    }

    template<class T2,int d2>
    explicit VECTOR(const VECTOR<T2,d2>& v)
    {
        STATIC_ASSERT(d2<=d);
        for(int i=0;i<d2;i++) array[i]=(T)v[i+1];
        for(int i=d2;i<d;i++) array[i]=T();
    }

    template<class T_VECTOR,class T_INDICES>
    explicit VECTOR(const INDIRECT_ARRAY<T_VECTOR,T_INDICES>& v)
    {
        STATIC_ASSERT((IS_SAME<T,typename INDIRECT_ARRAY<T_VECTOR,T_INDICES>::ELEMENT>::value && INDIRECT_ARRAY<T_VECTOR,T_INDICES>::m==d));
        for(int i=0;i<d;i++) array[i]=v(i+1);
    }

    VECTOR(const VECTOR& v)
        :BASE()
    {
        for(int i=0;i<d;i++) array[i]=v.array[i];
    }

    template<class T_VECTOR>
    explicit VECTOR(const VECTOR_BASE<T,T_VECTOR>& v)
    {
        assert(d==v.Size());for(int i=0;i<d;i++) array[i]=v(i+1);
    }

    template<class T_VECTOR> // TODO: This constructor should go away.
    VECTOR(const VECTOR_EXPRESSION<T,T_VECTOR>& v,typename ENABLE_IF<IS_SAME<typename VECTOR_TYPE<T_VECTOR>::TYPE,VECTOR>::value,UNUSABLE>::TYPE unusable=UNUSABLE())
    {
        assert(d==v.Size());for(int i=0;i<d;i++) array[i]=v(i+1);
    }

    template<int n>
    VECTOR(const VECTOR<T,n>& v1,const VECTOR<T,d-n>& v2)
    {
        for(int i=1;i<=n;i++) (*this)(i)=v1(i);for(int i=n+1;i<=d;i++) (*this)(i)=v2(i-n);
    }

    template<class T_VECTOR> typename ENABLE_IF<AND<IS_SAME<T,typename T_VECTOR::ELEMENT>::value,INTS_EQUAL<T_VECTOR::m,d>::value>::value,VECTOR&>::TYPE
    operator=(const T_VECTOR& v)
    {
        for(int i=0;i<d;i++) array[i]=v(i+1);return *this;
    }

    VECTOR& operator=(const VECTOR& v)
    {
        for(int i=0;i<d;i++) array[i]=v.array[i];return *this;
    }

    template<class T_VECTOR>
    VECTOR& operator=(const VECTOR_BASE<T,T_VECTOR>& v)
    {
        assert(d==v.Size());for(int i=0;i<d;i++) array[i]=v(i+1);return *this;
    }

    int Size() const
    {return m;}

    const T& operator[](const int i) const
    {assert(1<=i && i<=d);return array[i-1];}

    T& operator[](const int i)
    {assert(1<=i && i<=d);return array[i-1];}

    const T& operator()(const int i) const
    {assert(1<=i && i<=d);return array[i-1];}

    T& operator()(const int i)
    {assert(1<=i && i<=d);return array[i-1];}

    bool operator==(const VECTOR& v) const
    {for(int i=0;i<d;i++) if(array[i]!=v.array[i]) return false;return true;}

    bool operator!=(const VECTOR& v) const
    {return !((*this)==v);}

    VECTOR operator+(const T& a) const
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=array[i]+a;return r;}

    VECTOR operator-(const T& a) const
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=array[i]-a;return r;}

    VECTOR operator*(const T& a) const
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=array[i]*a;return r;}

    VECTOR operator*(const VECTOR& v) const
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=array[i]*v.array[i];return r;}

    VECTOR operator/(const VECTOR& v) const
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=array[i]/v.array[i];return r;}

    VECTOR operator*(const INT_INVERSE a) const
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=array[i]*a;return r;}

    bool Elements_Equal() const
    {bool equal=true;for(int i=1;i<d;i++) equal&=(array[i]==array[0]);return equal;}

    static VECTOR Componentwise_Min(const VECTOR& v1,const VECTOR& v2)
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=min(v1.array[i],v2.array[i]);return r;}

    static VECTOR Componentwise_Max(const VECTOR& v1,const VECTOR& v2)
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=max(v1.array[i],v2.array[i]);return r;}

    VECTOR Projected_On_Unit_Direction(const VECTOR& direction) const
    {return Dot_Product(*this,direction)*direction;}

    VECTOR Projected(const VECTOR& direction) const // un-normalized direction
    {return Dot_Product(*this,direction)/direction.Magnitude_Squared()*direction;}

    void Project_On_Unit_Direction(const VECTOR& direction)
    {*this=Dot_Product(*this,direction)*direction;}

    void Project(const VECTOR& direction) // un-normalized direction
    {*this=Dot_Product(*this,direction)/direction.Magnitude_Squared()*direction;}

    VECTOR Projected_Orthogonal_To_Unit_Direction(const VECTOR& direction) const
    {return *this-Dot_Product(*this,direction)*direction;}

    void Project_Orthogonal_To_Unit_Direction(const VECTOR& direction)
    {*this-=Dot_Product(*this,direction)*direction;}

    const VECTOR& Column_Sum() const
    {return *this;}

    int Number_True() const
    {STATIC_ASSERT((IS_SAME<T,bool>::value));int count=0;for(int i=0;i<d;i++)if(array[i]) count++;return count;}

    static VECTOR Axis_Vector(const int axis)
    {VECTOR r;r(axis)=(T)1;return r;}

    static VECTOR Constant_Vector(const T& constant)
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=constant;return r;}

    static VECTOR All_Ones_Vector()
    {return Constant_Vector(1);}

    void Fill(const T& constant)
    {for(int i=0;i<d;i++) array[i]=constant;}

    void Get(T& element1,T& element2,T& element3,T& element4) const
    {STATIC_ASSERT(d==4);element1=array[0];element2=array[1];element3=array[2];element4=array[3];}

    void Get(T& element1,T& element2,T& element3,T& element4,T& element5) const
    {STATIC_ASSERT(d==5);element1=array[0];element2=array[1];element3=array[2];element4=array[3];element5=array[4];}

    void Get(T& element1,T& element2,T& element3,T& element4,T& element5,T& element6) const
    {STATIC_ASSERT(d==6);element1=array[0];element2=array[1];element3=array[2];element4=array[3];element5=array[4];element6=array[5];}

    void Get(T& element1,T& element2,T& element3,T& element4,T& element5,T& element6,T& element7) const
    {STATIC_ASSERT(d==7);element1=array[0];element2=array[1];element3=array[2];element4=array[3];element5=array[4];element6=array[5];element7=array[6];}

    void Get(T& element1,T& element2,T& element3,T& element4,T& element5,T& element6,T& element7,T& element8) const
    {STATIC_ASSERT(d==8);element1=array[0];element2=array[1];element3=array[2];element4=array[3];element5=array[4];element6=array[5];element7=array[6];element8=array[7];}

    void Set(const T& element1,const T& element2,const T& element3,const T& element4)
    {STATIC_ASSERT(d==4);array[0]=element1;array[1]=element2;array[2]=element3;array[3]=element4;}

    void Set(const T& element1,const T& element2,const T& element3,const T& element4,const T& element5)
    {STATIC_ASSERT(d==5);array[0]=element1;array[1]=element2;array[2]=element3;array[3]=element4;array[4]=element5;}

    void Set(const T& element1,const T& element2,const T& element3,const T& element4,const T& element5,const T& element6)
    {STATIC_ASSERT(d==6);array[0]=element1;array[1]=element2;array[2]=element3;array[3]=element4;array[4]=element5;array[5]=element6;}

    void Set(const T& element1,const T& element2,const T& element3,const T& element4,const T& element5,const T& element6,const T& element7)
    {STATIC_ASSERT(d==7);array[0]=element1;array[1]=element2;array[2]=element3;array[3]=element4;array[4]=element5;array[5]=element6;array[6]=element7;}

    void Set(const T& element1,const T& element2,const T& element3,const T& element4,const T& element5,const T& element6,const T& element7,const T& element8)
    {STATIC_ASSERT(d==8);array[0]=element1;array[1]=element2;array[2]=element3;array[3]=element4;array[4]=element5;array[5]=element6;array[6]=element7;array[7]=element8;}

    template<class T_FUNCTION>
    static VECTOR Map(const T_FUNCTION& f,const VECTOR& v)
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=f(v.array[i]);return r;}

    int Find(const T& element) const
    {for(int i=0;i<d;i++) if(array[i]==element) return i+1;return 0;}

    bool Contains(const T& element) const
    {for(int i=0;i<d;i++) if(array[i]==element) return true;return false;}

    template<class T_ARRAY>
    bool Contains_All(const T_ARRAY& elements) const
    {STATIC_ASSERT((IS_SAME<typename T_ARRAY::ELEMENT,T>::value));
    for(int i=1;i<=elements.Size();i++) if(!Contains(elements(i))) return false;
    return true;}

    template<class T_ARRAY>
    bool Contains_Any(const T_ARRAY& elements) const
    {STATIC_ASSERT((IS_SAME<typename T_ARRAY::ELEMENT,T>::value));
    for(int i=1;i<=elements.Size();i++) if(Contains(elements(i))) return true;
    return false;}

    VECTOR<T,d-1> Remove_Index(const int index) const
    {assert(1<=index && index<=d);VECTOR<T,d-1> r;for(int i=1;i<=d-1;i++) r[i]=(*this)[i+(i>=index)];return r;}

    VECTOR<T,d+1> Insert(const T& element,const int index) const
    {VECTOR<T,d+1> r;r[index]=element;for(int i=1;i<=d;i++) r[i+(i>=index)]=(*this)[i];return r;}

    VECTOR<T,d+1> Append(const T& element) const
    {VECTOR<T,d+1> r;for(int i=1;i<=d;i++) r[i]=(*this)[i];r[d+1]=element;return r;}

    template<int d2> VECTOR<T,d+d2> Append_Elements(const VECTOR<T,d2>& elements) const
    {VECTOR<T,d+d2> r;
    for(int i=1;i<=d;i++) r[i]=(*this)[i];
    for(int i=1;i<=d2;i++) r[i+d]=elements[i];
    return r;}

    VECTOR<T,4> Sorted() const
    {STATIC_ASSERT(d==4);VECTOR<T,4> r(*this);exchange_sort(r[1],r[2],r[3],r[4]);return r;}

    VECTOR Reversed() const
    {VECTOR r;for(int i=0;i<d;i++) r.array[d-i-1]=array[i];return r;}

    template<int d1,int d2> VECTOR<T,d2-d1+1> Slice() const
    {STATIC_ASSERT(((1<=d1) && (d2<=d)));
    VECTOR<T,d2-d1+1> r;for(int i=d1;i<=d2;i++) r[i-d1+1]=(*this)[i];return r;}

    template<int n> void Split(VECTOR<T,n>& v1,VECTOR<T,d-n>& v2) const
    {for(int i=1;i<=n;i++) v1(i)=(*this)(i);
    for(int i=n+1;i<=d;i++) v2(i-n)=(*this)(i);}

    T* begin() // for stl
    {return array;}

    const T* begin() const // for stl
    {return array;}

    T* end() // for stl
    {return array+d;}

    const T* end() const // for stl
    {return array+d;}

//#####################################################################
    int Dominant_Axis() const;
//#####################################################################
};

//#####################################################################
// Miscellaneous free operators and functions
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
operator+(const T& a,const VECTOR<T,d>& v)
{VECTOR<T,d> r;for(int i=0;i<d;i++) r.array[i]=a+v.array[i];return r;}

template<class T,int d> inline VECTOR<T,d>
operator-(const T& a,const VECTOR<T,d>& v)
{VECTOR<T,d> r;for(int i=0;i<d;i++) r.array[i]=a-v.array[i];return r;}

template<class T,int d> inline VECTOR<T,d>
operator/(const T& a,const VECTOR<T,d>& v)
{VECTOR<T,d> r;for(int i=0;i<d;i++) r.array[i]=a/v.array[i];return r;}

template<class T,int d> inline VECTOR<T,d>
Inverse(const VECTOR<T,d>& v)
{VECTOR<T,d> r;for(int i=0;i<d;i++) r.array[i]=1/v.array[i];return r;}

template<class T,int d> VECTOR<T,d> abs(const VECTOR<T,d>& v);
template<class T,int d> VECTOR<T,d> floor(const VECTOR<T,d>& v);
template<class T,int d> VECTOR<T,d> ceil(const VECTOR<T,d>& v);
template<class T,int d> VECTOR<T,d> rint(const VECTOR<T,d>& v);
template<class T,int d> VECTOR<T,d> exp(const VECTOR<T,d>& v);
template<class T,int d> VECTOR<T,d> sin(const VECTOR<T,d>& v);
template<class T,int d> VECTOR<T,d> cos(const VECTOR<T,d>& v);
template<class T,int d> VECTOR<T,d> sqrt(const VECTOR<T,d>& v);

//#####################################################################
// Functions clamp, clamp_min, clamp_max, in_bounds
//#####################################################################
template<class T,int d> VECTOR<T,d> clamp(const VECTOR<T,d>& v,const VECTOR<T,d>& vmin,const VECTOR<T,d>& vmax);
template<class T,int d> VECTOR<T,d> clamp(const VECTOR<T,d>& v,const T& min,const T& max);
template<class T,int d> VECTOR<T,d> clamp_min(const VECTOR<T,d>& v,const VECTOR<T,d>& vmin);
template<class T,int d> VECTOR<T,d> clamp_min(const VECTOR<T,d>& v,const T& min);
template<class T,int d> VECTOR<T,d> clamp_max(const VECTOR<T,d>& v,const VECTOR<T,d>& vmax);
template<class T,int d> VECTOR<T,d> clamp_max(const VECTOR<T,d>& v,const T& max);
template<class T,int d> VECTOR<T,d> in_bounds(const VECTOR<T,d>& v,const VECTOR<T,d>& vmin,const VECTOR<T,d>& vmax);

//#####################################################################
// Vector construction
//#####################################################################
template<class T> VECTOR<T,1>
Vector(const T& d1)
{return VECTOR<T,1>(d1);}

template<class T> VECTOR<T,2>
Vector(const T& d1,const T& d2)
{return VECTOR<T,2>(d1,d2);}

template<class T> VECTOR<T,3>
Vector(const T& d1,const T& d2,const T& d3)
{return VECTOR<T,3>(d1,d2,d3);}

template<class T> VECTOR<T,4>
Vector(const T& d1,const T& d2,const T& d3,const T& d4)
{return VECTOR<T,4>(d1,d2,d3,d4);}

template<class T> VECTOR<T,5>
Vector(const T& d1,const T& d2,const T& d3,const T& d4,const T& d5)
{return VECTOR<T,5>(d1,d2,d3,d4,d5);}

template<class T> VECTOR<T,6>
Vector(const T& d1,const T& d2,const T& d3,const T& d4,const T& d5,const T& d6)
{return VECTOR<T,6>(d1,d2,d3,d4,d5,d6);}

template<class T> VECTOR<T,7>
Vector(const T& d1,const T& d2,const T& d3,const T& d4,const T& d5,const T& d6,const T& d7)
{return VECTOR<T,7>(d1,d2,d3,d4,d5,d6,d7);}

template<class T> VECTOR<T,8>
Vector(const T& d1,const T& d2,const T& d3,const T& d4,const T& d5,const T& d6,const T& d7,const T& d8)
{return VECTOR<T,8>(d1,d2,d3,d4,d5,d6,d7,d8);}
//#####################################################################
template<class T,int d> struct SUM<VECTOR<T,d>,T,typename ENABLE_IF<(d>3)>::TYPE>{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct SUM<T,VECTOR<T,d>,typename ENABLE_IF<(d>3)>::TYPE>{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct DIFFERENCE<VECTOR<T,d>,T,typename ENABLE_IF<(d>3)>::TYPE>{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct DIFFERENCE<T,VECTOR<T,d>,typename ENABLE_IF<(d>3)>::TYPE>{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct PRODUCT<VECTOR<T,d>,VECTOR<T,d>,typename ENABLE_IF<(d>3)>::TYPE>{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct QUOTIENT<VECTOR<T,d>,VECTOR<T,d>,typename ENABLE_IF<(d>3)>::TYPE>{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct QUOTIENT<T,VECTOR<T,d>,typename ENABLE_IF<(d>3)>::TYPE>{typedef VECTOR<T,d> TYPE;};
//#####################################################################
}
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR.h>
#endif
#endif

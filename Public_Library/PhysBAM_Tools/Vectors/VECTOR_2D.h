//#####################################################################
// Copyright 2002-2008, Robert Bridson, Kevin Der, Doug Enright, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Igor Neverov, Duc Nguyen, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VECTOR_2D
//#####################################################################
#ifndef __VECTOR_2D__
#define __VECTOR_2D__

#include <PhysBAM_Tools/Math_Tools/argmax.h>
#include <PhysBAM_Tools/Math_Tools/argmin.h>
#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/exchange_sort.h>
#include <PhysBAM_Tools/Math_Tools/Inverse.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Tools/Math_Tools/wrap.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
#include <cmath>
#ifdef DIFFERENCE // Windows workaround.
#undef DIFFERENCE
#endif
namespace PhysBAM{

using ::std::floor;
using ::std::ceil;
using ::std::sin;
using ::std::cos;

template<class T>
class VECTOR<T,2>:public VECTOR_BASE<T,VECTOR<T,2> >
{
    struct UNUSABLE{};
public:
    template<class T2> struct REBIND{typedef VECTOR<T2,2> TYPE;};
    typedef typename IF<IS_SCALAR<T>::value,T,UNUSABLE>::TYPE SCALAR;
    typedef T ELEMENT;
    typedef VECTOR<T,1> SPIN;
    typedef int INDEX;
    typedef T value_type; // for stl
    typedef T* iterator; // for stl
    enum WORKAROUND1 {dimension=2};
    enum WORKAROUND2 {m=2};

    T x,y;

    explicit VECTOR(INITIAL_SIZE n=INITIAL_SIZE(2))
        :x(),y()
    {
        STATIC_ASSERT(sizeof(VECTOR)==2*sizeof(T));assert(n==INITIAL_SIZE(2));
    }

    VECTOR(const T& x_input,const T& y_input)
        :x(x_input),y(y_input)
    {}

    VECTOR(const VECTOR& vector_input)
        :x(vector_input.x),y(vector_input.y)
    {}

    template<class T2> explicit VECTOR(const VECTOR<T2,2>& vector_input)
        :x((T)vector_input.x),y((T)vector_input.y)
    {}

    template<class T2> explicit VECTOR(const VECTOR<T2,1>& vector_input)
        :x((T)vector_input.x),y(T())
    {}

    template<class T_VECTOR>
    explicit VECTOR(const T_VECTOR& v)
        :x(v(1)),y(v(2))
    {
        STATIC_ASSERT((IS_SAME<T,typename T_VECTOR::ELEMENT>::value && T_VECTOR::m==2));
    }

    explicit VECTOR(const ARRAY<T>& v)
        :x(v(1)),y(v(2))
    {
        assert(2==v.Size());
    }

    template<class T_ARRAY>
    explicit VECTOR(const ARRAY_BASE<T,T_ARRAY>& v)
    {
        const T_ARRAY& v_=v.Derived();
        assert(2==v_.Size());
        x=v_(1);y=v_(2); // doing this here instead of as initializers dodges a bug in 4.1.1
    }

    template<class T_VECTOR>
    explicit VECTOR(const VECTOR_BASE<T,T_VECTOR>& v)
        :x(v(1)),y(v(2))
    {
        assert(2==v.Size());
    }

    VECTOR(const VECTOR_ND<T>& vector_input)
        :x(vector_input(1)),y(vector_input(2))
    {
        assert(vector_input.n==2);
    }

    template<class T_VECTOR> typename ENABLE_IF<AND<IS_SAME<T,typename T_VECTOR::ELEMENT>::value,INTS_EQUAL<T_VECTOR::m,2>::value>::value,VECTOR&>::TYPE
    operator=(const T_VECTOR& v)
    {
        x=v(1);y=v(2);return *this;
    }

    VECTOR& operator=(const VECTOR& v)
    {
        x=v(1);y=v(2);return *this;
    }

    int Size() const
    {return 2;}

    const T& operator[](const int i) const
    {assert(1<=i && i<=2);return *((const T*)(this)+i-1);}

    T& operator[](const int i)
    {assert(1<=i && i<=2);return *((T*)(this)+i-1);}

    const T& operator()(const int i) const
    {assert(1<=i && i<=2);return *((const T*)(this)+i-1);}

    T& operator()(const int i)
    {assert(1<=i && i<=2);return *((T*)(this)+i-1);}

    bool operator==(const VECTOR& v) const
    {return x==v.x && y==v.y;}

    bool operator!=(const VECTOR& v) const
    {return x!=v.x || y!=v.y;}

    VECTOR operator-() const
    {return VECTOR(-x,-y);}

    VECTOR& operator+=(const VECTOR& v)
    {x+=v.x;y+=v.y;return *this;}

    VECTOR& operator-=(const VECTOR& v)
    {x-=v.x;y-=v.y;return *this;}

    VECTOR& operator*=(const VECTOR& v)
    {x*=v.x;y*=v.y;return *this;}

    VECTOR& operator+=(const T& a)
    {x+=a;y+=a;return *this;}

    VECTOR& operator-=(const T& a)
    {x-=a;y-=a;return *this;}

    VECTOR& operator*=(const T& a)
    {x*=a;y*=a;return *this;}

    VECTOR& operator*=(const INT_INVERSE a)
    {x*=a;y*=a;return *this;}

    VECTOR& operator/=(const T& a)
    {return *this*=Inverse(a);}

    VECTOR& operator/=(const VECTOR& v)
    {x/=v.x;y/=v.y;return *this;}

    VECTOR operator+(const VECTOR& v) const
    {return VECTOR(x+v.x,y+v.y);}

    VECTOR operator-(const VECTOR& v) const
    {return VECTOR(x-v.x,y-v.y);}

    VECTOR operator*(const VECTOR& v) const
    {return VECTOR(x*v.x,y*v.y);}

    VECTOR operator/(const VECTOR& v) const
    {return VECTOR(x/v.x,y/v.y);}

    VECTOR operator+(const T& a) const
    {return VECTOR(x+a,y+a);}

    VECTOR operator-(const T& a) const
    {return VECTOR(x-a,y-a);}

    VECTOR operator*(const T& a) const
    {return VECTOR(x*a,y*a);}

    VECTOR operator*(const INT_INVERSE a) const
    {return VECTOR(x*a,y*a);}

    VECTOR operator/(const T& a) const
    {return *this*Inverse(a);}

    T Magnitude_Squared() const
    {return sqr(x)+sqr(y);}

    T Magnitude() const
    {return sqrt(sqr(x)+sqr(y));}

    T Lp_Norm(const T& p) const
    {return pow(pow(abs(x),p)+pow(abs(y),p),1/p);}

    T L1_Norm() const
    {return abs(x)+abs(y);}

    T Normalize()
    {T magnitude=Magnitude();if(magnitude) *this*=1/magnitude;else *this=VECTOR(1,0);return magnitude;}

    VECTOR Normalized() const
    {T magnitude=Magnitude();if(magnitude) return *this*(1/magnitude);else return VECTOR(1,0);}

    VECTOR Rotate_Clockwise_90() const
    {return VECTOR(y,-x);}

    VECTOR Rotate_Counterclockwise_90() const
    {return VECTOR(-y,x);}

    VECTOR Rotate_Counterclockwise_Multiple_90(const int n) const
    {VECTOR r(*this);if(n&2) r=-r;
    return n&1?r.Rotate_Counterclockwise_90():r;}

    VECTOR Perpendicular() const
    {return VECTOR(-y,x);}

    VECTOR Orthogonal_Vector() const
    {return VECTOR(-y,x);}

    VECTOR Unit_Orthogonal_Vector() const
    {return Orthogonal_Vector().Normalized();}

    T Min() const
    {return min(x,y);}

    T Max() const
    {return max(x,y);}

    T Max_Abs() const
    {return maxabs(x,y);}

    int Arg_Min() const
    {return argmin(x,y);}

    int Arg_Max() const
    {return argmax(x,y);}

    bool Elements_Equal() const
    {return x==y;}

    int Dominant_Axis() const
    {return (abs(x)>abs(y))?1:2;}

    static T Dot_Product(const VECTOR& v1,const VECTOR& v2)
    {return v1.x*v2.x+v1.y*v2.y;}

    static VECTOR Componentwise_Min(const VECTOR& v1,const VECTOR& v2)
    {return VECTOR(min(v1.x,v2.x),min(v1.y,v2.y));}

    static VECTOR Componentwise_Max(const VECTOR& v1,const VECTOR& v2)
    {return VECTOR(max(v1.x,v2.x),max(v1.y,v2.y));}

    bool All_Greater(const VECTOR& v) const
    {return x>v.x && y>v.y;}

    bool All_Less(const VECTOR& v) const
    {return x<v.x && y<v.y;}

    bool All_Greater_Equal(const VECTOR& v) const
    {return x>=v.x && y>=v.y;}

    bool All_Less_Equal(const VECTOR& v) const
    {return x<=v.x && y<=v.y;}

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

    static VECTOR<T,1> Cross_Product(const VECTOR& v1,const VECTOR& v2)
    {return VECTOR<T,1>(v1.x*v2.y-v1.y*v2.x);}

    static VECTOR Cross_Product(const VECTOR& v1,const VECTOR<T,1>& v2) // v2 is out of plane
    {return VECTOR(v1.y*v2.x,-v1.x*v2.x);}

    static VECTOR Cross_Product(const VECTOR<T,1>& v1,const VECTOR& v2) // v1 is out of plane
    {return VECTOR(-v1.x*v2.y,v1.x*v2.x);}

    static T Angle_Between(const VECTOR& u,const VECTOR& v)
    {T s=Cross_Product(u,v).Magnitude(),c=Dot_Product(u,v);return atan2(s,c);}

    static T Oriented_Angle_Between(const VECTOR& u,const VECTOR& v)
    {T s=Cross_Product(u,v).x,c=Dot_Product(u,v);return atan2(s,c);}

    T Sum() const
    {return x+y;}

    T Average() const
    {return (T).5*Sum();}

    T Product() const
    {return x*y;}

    const VECTOR& Column_Sum() const
    {return *this;}

    int Number_True() const
    {STATIC_ASSERT((IS_SAME<T,bool>::value));return x+y;}

    static VECTOR Axis_Vector(const int axis)
    {VECTOR vec;vec[axis]=(T)1;return vec;}

    static VECTOR Constant_Vector(const T& constant)
    {return VECTOR(constant,constant);}

    static VECTOR All_Ones_Vector()
    {return Constant_Vector(1);}

    VECTOR<T,1> Horizontal_Vector() const
    {return VECTOR<T,1>(x);}

    VECTOR<T,1> Vertical_Vector() const
    {return VECTOR<T,1>(y);}

    void Fill(const T& constant)
    {x=y=constant;}

    void Get(T& element1,T& element2) const
    {element1=x;element2=y;}

    void Set(const T& element1,const T& element2)
    {x=element1;y=element2;}

    template<class T_FUNCTION>
    static VECTOR Map(const T_FUNCTION& f,const VECTOR& v)
    {return VECTOR(f(v.x),f(v.y));}

    int Find(const T& element) const
    {return x==element?1:y==element?2:0;}

    bool Contains(const T& element) const
    {return x==element || y==element;}

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

    VECTOR<T,1> Remove_Index(const int index) const
    {assert(1<=index && index<=2);return VECTOR<T,1>((*this)[3-index]);}

    VECTOR<T,3> Insert(const T& element,const int index) const
    {VECTOR<T,3> r;r[index]=element;for(int i=1;i<=2;i++) r[i+(i>=index)]=(*this)[i];return r;}

    VECTOR<T,3> Append(const T& element) const
    {return VECTOR<T,3>(x,y,element);}

    template<int d2> VECTOR<T,2+d2> Append_Elements(const VECTOR<T,d2>& elements) const
    {VECTOR<T,2+d2> r;r[1]=x;r[2]=y;for(int i=1;i<=d2;i++) r[i+2]=elements[i];return r;}

    VECTOR<T,2> Sorted() const
    {VECTOR<T,2> r(*this);exchange_sort(r.x,r.y);return r;}

    VECTOR Reversed() const
    {return VECTOR(y,x);}

    template<int d1,int d2> VECTOR<T,d2-d1+1> Slice() const
    {STATIC_ASSERT(((1<=d1) && (d2<=2)));
    VECTOR<T,d2-d1+1> r;for(int i=d1;i<=d2;i++) r[i-d1+1]=(*this)[i];return r;}

    template<int n> void Split(VECTOR<T,n>& v1,VECTOR<T,2-n>& v2) const
    {for(int i=1;i<=n;i++) v1(i)=(*this)(i);
    for(int i=n+1;i<=2;i++) v2(i-n)=(*this)(i);}

    template<class T_VECTOR>
    void Set_Subvector(const int istart,const T_VECTOR& v)
    {for(int i=1;i<=v.Size();i++) (*this)(istart+i-1)=v(i);}

    template<class T_VECTOR>
    void Add_Subvector(const int istart,const T_VECTOR& v)
    {for(int i=1;i<=v.Size();i++) (*this)(istart+i-1)+=v(i);}
    
    template<class T_VECTOR>
    void Get_Subvector(const int istart,T_VECTOR& v) const
    {for(int i=1;i<=v.Size();i++) v(i)=(*this)(istart+i-1);}

    T* begin() // for stl
    {return &x;}

    const T* begin() const // for stl
    {return &x;}

    T* end() // for stl
    {return &y+1;}

    const T* end() const // for stl
    {return &y+1;}

//#####################################################################
};

//#####################################################################
// Miscellaneous free operators and functions
//#####################################################################
template<class T> inline VECTOR<T,2>
operator+(const T& a,const VECTOR<T,2>& v)
{return VECTOR<T,2>(a+v.x,a+v.y);}

template<class T> inline VECTOR<T,2>
operator-(const T& a,const VECTOR<T,2>& v)
{return VECTOR<T,2>(a-v.x,a-v.y);}

template<class T> inline VECTOR<T,2>
operator*(const T& a,const VECTOR<T,2>& v)
{return VECTOR<T,2>(a*v.x,a*v.y);}

template<class T> inline VECTOR<T,2>
operator/(const T& a,const VECTOR<T,2>& v)
{return VECTOR<T,2>(a/v.x,a/v.y);}

template<class T> inline VECTOR<T,2>
abs(const VECTOR<T,2>& v)
{return VECTOR<T,2>(abs(v.x),abs(v.y));}

template<class T> inline VECTOR<T,2>
floor(const VECTOR<T,2>& v)
{return VECTOR<T,2>(floor(v.x),floor(v.y));}

template<class T> inline VECTOR<T,2>
ceil(const VECTOR<T,2>& v)
{return VECTOR<T,2>(ceil(v.x),ceil(v.y));}

template<class T> inline VECTOR<T,2>
rint(const VECTOR<T,2>& v)
{return VECTOR<T,2>(rint(v.x),rint(v.y));}

template<class T> inline VECTOR<T,2>
exp(const VECTOR<T,2>& v)
{return VECTOR<T,2>(exp(v.x),exp(v.y));}

template<class T> inline VECTOR<T,2>
sin(const VECTOR<T,2>& v)
{return VECTOR<T,2>(sin(v.x),sin(v.y));}

template<class T> inline VECTOR<T,2>
cos(const VECTOR<T,2>& v)
{return VECTOR<T,2>(cos(v.x),cos(v.y));}

template<class T> inline VECTOR<T,2>
sqrt(const VECTOR<T,2>& v)
{return VECTOR<T,2>(sqrt(v.x),sqrt(v.y));}

template<class T> inline VECTOR<T,2>
Inverse(const VECTOR<T,2>& v)
{return VECTOR<T,2>(1/v.x,1/v.y);}

//#####################################################################
// Functions clamp, clamp_min, clamp_max, in_bounds
//#####################################################################
template<class T> inline VECTOR<T,2>
clamp(const VECTOR<T,2>& v,const VECTOR<T,2>& vmin,const VECTOR<T,2>& vmax)
{return VECTOR<T,2>(clamp(v.x,vmin.x,vmax.x),clamp(v.y,vmin.y,vmax.y));}

template<class T> inline VECTOR<T,2>
clamp(const VECTOR<T,2>& v,T min,T max)
{return VECTOR<T,2>(clamp(v.x,min,max),clamp(v.y,min,max));}

template<class T> inline VECTOR<T,2>
clamp_min(const VECTOR<T,2>& v,const VECTOR<T,2>& vmin)
{return VECTOR<T,2>(clamp_min(v.x,vmin.x),clamp_min(v.y,vmin.y));}

template<class T> inline VECTOR<T,2>
clamp_min(const VECTOR<T,2>& v,const T& min)
{return VECTOR<T,2>(clamp_min(v.x,min),clamp_min(v.y,min));}

template<class T> inline VECTOR<T,2>
clamp_max(const VECTOR<T,2>& v,const VECTOR<T,2>& vmax)
{return VECTOR<T,2>(clamp_max(v.x,vmax.x),clamp_max(v.y,vmax.y));}

template<class T> inline VECTOR<T,2>
clamp_max(const VECTOR<T,2>& v,const T& max)
{return VECTOR<T,2>(clamp_max(v.x,max),clamp_max(v.y,max));}

template<class T> inline bool
in_bounds(const VECTOR<T,2>& v,const VECTOR<T,2>& vmin,const VECTOR<T,2>& vmax)
{return in_bounds(v.x,vmin.x,vmax.x) && in_bounds(v.y,vmin.y,vmax.y);}

template<class T> inline VECTOR<T,2>
wrap(const VECTOR<T,2>& v,const VECTOR<T,2>& vmin,const VECTOR<T,2>& vmax)
{return VECTOR<T,2>(wrap(v.x,vmin.x,vmax.x),wrap(v.y,vmin.y,vmax.y));}
//#####################################################################
template<class T> struct SUM<VECTOR<T,2>,VECTOR<T,2> >{typedef VECTOR<T,2> TYPE;};
template<class T> struct SUM<VECTOR<T,2>,T>{typedef VECTOR<T,2> TYPE;};
template<class T> struct SUM<T,VECTOR<T,2> >{typedef VECTOR<T,2> TYPE;};
template<class T> struct DIFFERENCE<VECTOR<T,2>,VECTOR<T,2> >{typedef VECTOR<T,2> TYPE;};
template<class T> struct DIFFERENCE<VECTOR<T,2>,T>{typedef VECTOR<T,2> TYPE;};
template<class T> struct DIFFERENCE<T,VECTOR<T,2> >{typedef VECTOR<T,2> TYPE;};
template<class T> struct PRODUCT<VECTOR<T,2>,VECTOR<T,2> >{typedef VECTOR<T,2> TYPE;};
template<class T> struct PRODUCT<VECTOR<T,2>,T>{typedef VECTOR<T,2> TYPE;};
template<class T> struct PRODUCT<T,VECTOR<T,2> >{typedef VECTOR<T,2> TYPE;};
template<class T> struct QUOTIENT<VECTOR<T,2>,VECTOR<T,2> >{typedef VECTOR<T,2> TYPE;};
template<class T> struct QUOTIENT<VECTOR<T,2>,T>{typedef VECTOR<T,2> TYPE;};
template<class T> struct QUOTIENT<T,VECTOR<T,2> >{typedef VECTOR<T,2> TYPE;};
template<class T> struct NEGATION<VECTOR<T,2> >{typedef VECTOR<T,2> TYPE;};
//#####################################################################
}
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_2D.h>
#endif
#endif

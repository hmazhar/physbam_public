//#####################################################################
// Copyright 2002-2008, Robert Bridson, Doug Enright, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Igor Neverov, Duc Nguyen, Craig Schroeder, Tamar Shinar, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VECTOR_3D
//#####################################################################
#ifndef __VECTOR_3D__
#define __VECTOR_3D__
#include <PhysBAM_Tools/Math_Tools/argmax.h>
#include <PhysBAM_Tools/Math_Tools/argmin.h>
#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/Inverse.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Tools/Math_Tools/wrap.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <cmath>
#ifdef DIFFERENCE // Windows workaround.
#undef DIFFERENCE
#endif
namespace PhysBAM{

using ::std::exp;
using ::std::sin;
using ::std::cos;

template<class T_ARRAY,class T_INDICES> class INDIRECT_ARRAY;

template<class T>
class VECTOR<T,3>:public VECTOR_BASE<T,VECTOR<T,3> >
{
    struct UNUSABLE{};
public:
    template<class T2> struct REBIND{typedef VECTOR<T2,3> TYPE;};
    typedef typename IF<IS_SCALAR<T>::value,T,UNUSABLE>::TYPE SCALAR;
    typedef T ELEMENT;
    typedef VECTOR<T,3> SPIN;
    typedef int INDEX;
    typedef T value_type; // for stl
    typedef T* iterator; // for stl
    enum WORKAROUND1 {dimension=3};
    enum WORKAROUND2 {m=3};

    T x,y,z;

    explicit VECTOR(INITIAL_SIZE n=INITIAL_SIZE(3))
        :x(),y(),z()
    {
        STATIC_ASSERT(sizeof(VECTOR)==3*sizeof(T));assert(n==INITIAL_SIZE(3));
    }

    VECTOR(const T& x_input,const T& y_input,const T& z_input)
        :x(x_input),y(y_input),z(z_input)
    {}

    VECTOR(const VECTOR& vector_input)
        :x(vector_input.x),y(vector_input.y),z(vector_input.z)
    {}

    template<class T2> explicit VECTOR(const VECTOR<T2,3>& vector_input)
        :x((T)vector_input.x),y((T)vector_input.y),z((T)vector_input.z)
    {}

    explicit VECTOR(const VECTOR<T,2>& vector_input)
        :x(vector_input.x),y(vector_input.y),z(0)
    {}

    template<class T_VECTOR,class T_INDICES>
    explicit VECTOR(const INDIRECT_ARRAY<T_VECTOR,T_INDICES>& v)
        :x(v(1)),y(v(2)),z(v(3))
    {
        STATIC_ASSERT((IS_SAME<T,typename INDIRECT_ARRAY<T_VECTOR,T_INDICES>::ELEMENT>::value && INDIRECT_ARRAY<T_VECTOR,T_INDICES>::m==3));
    }

    explicit VECTOR(const ARRAY<T>& v)
        :x(v(1)),y(v(2)),z(v(3))
    {
        assert(3==v.Size());
    }

    template<class T_VECTOR>
    explicit VECTOR(const VECTOR_BASE<T,T_VECTOR>& v)
        :x(v(1)),y(v(2)),z(v(3))
    {
        assert(3==v.Size());
    }

    VECTOR(const VECTOR_ND<T>& vector_input)
        :x(vector_input(1)),y(vector_input(2)),z(vector_input(3))
    {
        assert(vector_input.n==3);
    }

    template<int n>
    VECTOR(const VECTOR<T,n>& v1,const VECTOR<T,3-n>& v2)
    {
        for(int i=1;i<=n;i++) (*this)(i)=v1(i);for(int i=n+1;i<=3;i++) (*this)(i)=v2(i-n);
    }

    template<class T_VECTOR> typename ENABLE_IF<AND<IS_SAME<T,typename T_VECTOR::ELEMENT>::value,INTS_EQUAL<T_VECTOR::m,3>::value>::value,VECTOR&>::TYPE
    operator=(const T_VECTOR& v)
    {
        x=v(1);y=v(2);z=v(3);return *this;
    }

    VECTOR& operator=(const VECTOR& v)
    {
        x=v(1);y=v(2);z=v(3);return *this;
    }

    int Size() const
    {return 3;}

    const T& operator[](const int i) const
    {assert(1<=i && i<=3);return *((const T*)(this)+i-1);}

    T& operator[](const int i)
    {assert(1<=i && i<=3);return *((T*)(this)+i-1);}

    const T& operator()(const int i) const
    {assert(1<=i && i<=3);return *((const T*)(this)+i-1);}

    T& operator()(const int i)
    {assert(1<=i && i<=3);return *((T*)(this)+i-1);}

    bool operator==(const VECTOR& v) const
    {return x==v.x && y==v.y && z==v.z;}

    bool operator!=(const VECTOR& v) const
    {return x!=v.x || y!=v.y || z!=v.z;}

    VECTOR operator-() const
    {return VECTOR(-x,-y,-z);}

    VECTOR& operator+=(const VECTOR& v)
    {x+=v.x;y+=v.y;z+=v.z;return *this;}

    VECTOR& operator-=(const VECTOR& v)
    {x-=v.x;y-=v.y;z-=v.z;return *this;}

    VECTOR& operator*=(const VECTOR& v)
    {x*=v.x;y*=v.y;z*=v.z;return *this;}

    VECTOR& operator+=(const T& a)
    {x+=a;y+=a;z+=a;return *this;}

    VECTOR& operator-=(const T& a)
    {x-=a;y-=a;z-=a;return *this;}

    VECTOR& operator*=(const T& a)
    {x*=a;y*=a;z*=a;return *this;}

    VECTOR& operator*=(const INT_INVERSE a)
    {x*=a;y*=a;z*=a;return *this;}

    VECTOR& operator/=(const T& a)
    {return *this*=Inverse(a);}

    VECTOR& operator/=(const VECTOR& v)
    {x/=v.x;y/=v.y;z/=v.z;return *this;}

    VECTOR operator+(const VECTOR& v) const
    {return VECTOR(x+v.x,y+v.y,z+v.z);}

    VECTOR operator-(const VECTOR& v) const
    {return VECTOR(x-v.x,y-v.y,z-v.z);}

    VECTOR operator*(const VECTOR& v) const
    {return VECTOR(x*v.x,y*v.y,z*v.z);}

    VECTOR operator/(const VECTOR& v) const
    {return VECTOR(x/v.x,y/v.y,z/v.z);}

    VECTOR operator+(const T& a) const
    {return VECTOR(x+a,y+a,z+a);}

    VECTOR operator-(const T& a) const
    {return VECTOR(x-a,y-a,z-a);}

    VECTOR operator*(const T& a) const
    {return VECTOR(x*a,y*a,z*a);}

    VECTOR operator*(const INT_INVERSE a) const
    {return VECTOR(x*a,y*a,z*a);}

    VECTOR operator/(const T& a) const
    {return *this*Inverse(a);}

    T Magnitude_Squared() const
    {return x*x+y*y+z*z;}

    T Magnitude() const
    {return sqrt(x*x+y*y+z*z);}

    T Lp_Norm(const T& p) const
    {return pow(pow(abs(x),p)+pow(abs(y),p)+pow(abs(z),p),1/p);}

    T L1_Norm() const
    {return abs(x)+abs(y)+abs(z);}

    T Normalize()
    {T magnitude=Magnitude();if(magnitude) *this*=1/magnitude;else *this=VECTOR(1,0,0);return magnitude;}

    VECTOR Normalized() const // 6 mults, 2 adds, 1 div, 1 sqrt
    {T magnitude=Magnitude();if(magnitude) return *this*(1/magnitude);else return VECTOR(1,0,0);}

    VECTOR Orthogonal_Vector() const
    {T abs_x=abs(x),abs_y=abs(y),abs_z=abs(z);
    if(abs_x<abs_y) return abs_x<abs_z?VECTOR(0,z,-y):VECTOR(y,-x,0);
    else return abs_y<abs_z?VECTOR(-z,0,x):VECTOR(y,-x,0);}

    VECTOR Unit_Orthogonal_Vector() const // roughly 6 mults, 2 adds, 1 div, 1 sqrt
    {return Orthogonal_Vector().Normalized();}

    VECTOR Distinct_Vector() const
    {return Axis_Vector(Arg_Abs_Min());}

    T Min() const
    {return min(x,y,z);}

    T Max() const
    {return max(x,y,z);}

    T Max_Abs() const
    {return maxabs(x,y,z);}

    int Arg_Min() const
    {return argmin(x,y,z);}

    int Arg_Max() const
    {return argmax(x,y,z);}

    int Arg_Abs_Min() const
    {return argmin(abs(x),abs(y),abs(z));}

    int Arg_Abs_Max() const
    {return argmax(abs(x),abs(y),abs(z));}

    bool Elements_Equal() const
    {return x==y && x==z;}

    bool All_Greater(const VECTOR& v) const
    {return x>v.x && y>v.y && z>v.z;}

    bool All_Less(const VECTOR& v) const
    {return x<v.x && y<v.y && z<v.z;}

    bool All_Greater_Equal(const VECTOR& v) const
    {return x>=v.x && y>=v.y && z>=v.z;}

    bool All_Less_Equal(const VECTOR& v) const
    {return x<=v.x && y<=v.y && z<=v.z;}

    int Dominant_Axis() const
    {T abs_x=abs(x),abs_y=abs(y),abs_z=abs(z);return (abs_x>abs_y && abs_x>abs_z) ? 1:((abs_y>abs_z) ? 2:3);}

    static T Dot_Product(const VECTOR& v1,const VECTOR& v2)
    {return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;}

    static VECTOR Componentwise_Min(const VECTOR& v1,const VECTOR& v2)
    {return VECTOR(min(v1.x,v2.x),min(v1.y,v2.y),min(v1.z,v2.z));}

    static VECTOR Componentwise_Max(const VECTOR& v1,const VECTOR& v2)
    {return VECTOR(max(v1.x,v2.x),max(v1.y,v2.y),max(v1.z,v2.z));}

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

    static VECTOR Cross_Product(const VECTOR& v1,const VECTOR& v2) // 6 mults, 3 adds
    {return VECTOR(v1.y*v2.z-v1.z*v2.y,v1.z*v2.x-v1.x*v2.z,v1.x*v2.y-v1.y*v2.x);}

    static T Angle_Between(const VECTOR& u,const VECTOR& v) // 0 .. pi
    {T s=Cross_Product(u,v).Magnitude(),c=Dot_Product(u,v);return atan2(s,c);}

    static T Triple_Product(const VECTOR& u,const VECTOR& v,const VECTOR& w)
    {return Dot_Product(u,Cross_Product(v,w));}

    T Sum() const
    {return x+y+z;}

    T Average() const
    {return (T)one_third*Sum();}

    T Product() const
    {return x*y*z;}

    const VECTOR& Column_Sum() const
    {return *this;}

    int Number_True() const
    {STATIC_ASSERT((IS_SAME<T,bool>::value));return x+y+z;}

    static VECTOR Axis_Vector(const int axis)
    {VECTOR vec;vec[axis]=(T)1;return vec;}

    static VECTOR Constant_Vector(const T& constant)
    {return VECTOR(constant,constant,constant);}

    static VECTOR All_Ones_Vector()
    {return Constant_Vector(1);}

    VECTOR<T,2> Horizontal_Vector() const
    {return VECTOR<T,2>(x,z);}

    void Fill(const T& constant)
    {x=y=z=constant;}

    void Get(T& element1,T& element2,T& element3) const
    {element1=x;element2=y;element3=z;}

    void Set(const T& element1,const T& element2,const T& element3)
    {x=element1;y=element2;z=element3;}

    template<class T_FUNCTION>
    static VECTOR Map(const T_FUNCTION& f,const VECTOR& v)
    {return VECTOR(f(v.x),f(v.y),f(v.z));}

    int Find(const T& element) const
    {return x==element?1:y==element?2:z==element?3:0;}

    bool Contains(const T& element) const
    {return x==element || y==element || z==element;}

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

    VECTOR<T,2> Remove_Index(const int index) const
    {assert(1<=index && index<=3);return VECTOR<T,2>(index>1?x:y,index<3?z:y);}

    VECTOR<T,4> Insert(const T& element,const int index) const
    {VECTOR<T,4> r;r[index]=element;for(int i=1;i<=3;i++) r[i+(i>=index)]=(*this)[i];return r;}

    VECTOR<T,4> Append(const T& element) const
    {return VECTOR<T,4>(x,y,z,element);}

    template<int d2> VECTOR<T,3+d2> Append_Elements(const VECTOR<T,d2>& elements) const
    {VECTOR<T,3+d2> r;r[1]=x;r[2]=y;r[3]=z;for(int i=1;i<=d2;i++) r[i+3]=elements[i];return r;}

    VECTOR<T,3> Sorted() const
    {VECTOR<T,3> r(*this);exchange_sort(r.x,r.y,r.z);return r;}

    VECTOR Reversed() const
    {return VECTOR(z,y,x);}

    template<int d1,int d2> VECTOR<T,d2-d1+1> Slice() const
    {STATIC_ASSERT(((1<=d1) && (d2<=3)));
    VECTOR<T,d2-d1+1> r;for(int i=d1;i<=d2;i++) r[i-d1+1]=(*this)[i];return r;}

    template<int n> void Split(VECTOR<T,n>& v1,VECTOR<T,3-n>& v2) const
    {for(int i=1;i<=n;i++) v1(i)=(*this)(i);
    for(int i=n+1;i<=3;i++) v2(i-n)=(*this)(i);}

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
    {return &z+1;}

    const T* end() const // for stl
    {return &z+1;}

//#####################################################################
};

//#####################################################################
// Miscellaneous free operators and functions
//#####################################################################
template<class T> inline VECTOR<T,3>
operator+(const T& a,const VECTOR<T,3>& v)
{return VECTOR<T,3>(a+v.x,a+v.y,a+v.z);}

template<class T> inline VECTOR<T,3>
operator-(const T& a,const VECTOR<T,3>& v)
{return VECTOR<T,3>(a-v.x,a-v.y,a-v.z);}

template<class T> inline VECTOR<T,3>
operator*(const T& a,const VECTOR<T,3>& v)
{return VECTOR<T,3>(a*v.x,a*v.y,a*v.z);}

template<class T> inline VECTOR<T,3>
operator/(const T& a,const VECTOR<T,3>& v)
{return VECTOR<T,3>(a/v.x,a/v.y,a/v.z);}

template<class T> inline VECTOR<T,3>
abs(const VECTOR<T,3>& v)
{return VECTOR<T,3>(abs(v.x),abs(v.y),abs(v.z));}

template<class T> inline VECTOR<T,3>
floor(const VECTOR<T,3>& v)
{return VECTOR<T,3>(floor(v.x),floor(v.y),floor(v.z));}

template<class T> inline VECTOR<T,3>
ceil(const VECTOR<T,3>& v)
{return VECTOR<T,3>(ceil(v.x),ceil(v.y),ceil(v.z));}

template<class T> inline VECTOR<T,3>
rint(const VECTOR<T,3>& v)
{return VECTOR<T,3>(rint(v.x),rint(v.y),rint(v.z));}

template<class T> inline VECTOR<T,3>
exp(const VECTOR<T,3>& v)
{return VECTOR<T,3>(exp(v.x),exp(v.y),exp(v.z));}

template<class T> inline VECTOR<T,3>
sin(const VECTOR<T,3>& v)
{return VECTOR<T,3>(sin(v.x),sin(v.y),sin(v.z));}

template<class T> inline VECTOR<T,3>
cos(const VECTOR<T,3>& v)
{return VECTOR<T,3>(cos(v.x),cos(v.y),cos(v.z));}

template<class T> inline VECTOR<T,3>
sqrt(const VECTOR<T,3>& v)
{return VECTOR<T,3>(sqrt(v.x),sqrt(v.y),sqrt(v.z));}

template<class T> inline VECTOR<T,3>
Inverse(const VECTOR<T,3>& v)
{return VECTOR<T,3>(1/v.x,1/v.y,1/v.z);}

//#####################################################################
// Functions clamp, clamp_min, clamp_max, in_bounds
//#####################################################################
template<class T> inline VECTOR<T,3>
clamp(const VECTOR<T,3>& v,const VECTOR<T,3>& vmin,const VECTOR<T,3>& vmax)
{return VECTOR<T,3>(clamp(v.x,vmin.x,vmax.x),clamp(v.y,vmin.y,vmax.y),clamp(v.z,vmin.z,vmax.z));}

template<class T> inline VECTOR<T,3>
clamp(const VECTOR<T,3>& v,T min,T max)
{return VECTOR<T,3>(clamp(v.x,min,max),clamp(v.y,min,max),clamp(v.z,min,max));}

template<class T> inline VECTOR<T,3>
clamp_min(const VECTOR<T,3>& v,const VECTOR<T,3>& vmin)
{return VECTOR<T,3>(clamp_min(v.x,vmin.x),clamp_min(v.y,vmin.y),clamp_min(v.z,vmin.z));}

template<class T> inline VECTOR<T,3>
clamp_min(const VECTOR<T,3>& v,const T& min)
{return VECTOR<T,3>(clamp_min(v.x,min),clamp_min(v.y,min),clamp_min(v.z,min));}

template<class T> inline VECTOR<T,3>
clamp_max(const VECTOR<T,3>& v,const VECTOR<T,3>& vmax)
{return VECTOR<T,3>(clamp_max(v.x,vmax.x),clamp_max(v.y,vmax.y),clamp_max(v.z,vmax.z));}

template<class T> inline VECTOR<T,3>
clamp_max(const VECTOR<T,3>& v,const T& max)
{return VECTOR<T,3>(clamp_max(v.x,max),clamp_max(v.y,max),clamp_max(v.z,max));}

template<class T> inline bool
in_bounds(const VECTOR<T,3>& v,const VECTOR<T,3>& vmin,const VECTOR<T,3>& vmax)
{return in_bounds(v.x,vmin.x,vmax.x) && in_bounds(v.y,vmin.y,vmax.y) && in_bounds(v.z,vmin.z,vmax.z);}

template<class T> inline VECTOR<T,3>
wrap(const VECTOR<T,3>& v,const VECTOR<T,3>& vmin,const VECTOR<T,3>& vmax)
{return VECTOR<T,3>(wrap(v.x,vmin.x,vmax.x),wrap(v.y,vmin.y,vmax.y),wrap(v.z,vmin.z,vmax.z));}

//#####################################################################
template<class T> struct SUM<VECTOR<T,3>,VECTOR<T,3> >{typedef VECTOR<T,3> TYPE;};
template<class T> struct SUM<VECTOR<T,3>,T>{typedef VECTOR<T,3> TYPE;};
template<class T> struct SUM<T,VECTOR<T,3> >{typedef VECTOR<T,3> TYPE;};
template<class T> struct DIFFERENCE<VECTOR<T,3>,VECTOR<T,3> >{typedef VECTOR<T,3> TYPE;};
template<class T> struct DIFFERENCE<VECTOR<T,3>,T>{typedef VECTOR<T,3> TYPE;};
template<class T> struct DIFFERENCE<T,VECTOR<T,3> >{typedef VECTOR<T,3> TYPE;};
template<class T> struct PRODUCT<VECTOR<T,3>,VECTOR<T,3> >{typedef VECTOR<T,3> TYPE;};
template<class T> struct PRODUCT<VECTOR<T,3>,T>{typedef VECTOR<T,3> TYPE;};
template<class T> struct PRODUCT<T,VECTOR<T,3> >{typedef VECTOR<T,3> TYPE;};
template<class T> struct QUOTIENT<VECTOR<T,3>,VECTOR<T,3> >{typedef VECTOR<T,3> TYPE;};
template<class T> struct QUOTIENT<VECTOR<T,3>,T>{typedef VECTOR<T,3> TYPE;};
template<class T> struct QUOTIENT<T,VECTOR<T,3> >{typedef VECTOR<T,3> TYPE;};
template<class T> struct NEGATION<VECTOR<T,3> >{typedef VECTOR<T,3> TYPE;};
//#####################################################################
}
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_3D.h>
#endif
#endif

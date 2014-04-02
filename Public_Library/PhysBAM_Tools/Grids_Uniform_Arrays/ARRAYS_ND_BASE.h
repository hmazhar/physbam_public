//#####################################################################
// Copyright 2007-2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAYS_ND_BASE
//#####################################################################
#ifndef __ARRAYS_ND_BASE__
#define __ARRAYS_ND_BASE__

#include <PhysBAM_Tools/Arrays/ARRAY_BASE.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/UNIFORM_ARRAY_ITERATOR.h>
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Log/LOG.h>
#endif
#include <PhysBAM_Tools/Math_Tools/ONE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T,class T_ARRAY,int dimension>
class ARRAY_BASE<T,T_ARRAY,VECTOR<int,dimension> >
{
    typedef VECTOR<int,dimension> TV_INT;
public:

    T_ARRAY& Derived()
    {return static_cast<T_ARRAY&>(*this);}

    const T_ARRAY& Derived() const
    {return static_cast<const T_ARRAY&>(*this);}

    template<class T_ARRAY1> static bool Same_Array(const T_ARRAY1& array1,const T_ARRAY1& array2)
    {return &array1==&array2;}

    template<class T_ARRAY1,class T_ARRAY2> static bool
    Same_Array(const T_ARRAY1& array1,const T_ARRAY2& array2)
    {return SAME_ARRAY<T_ARRAY1,T_ARRAY2>::Same_Array(array1,array2);}

    template<class T_INDICES>
    INDIRECT_ARRAY<T_ARRAY,T_INDICES&> Subset(const T_INDICES& indices)
    {return INDIRECT_ARRAY<T_ARRAY,T_INDICES&>(Derived(),indices);}

    template<class T_INDICES>
    INDIRECT_ARRAY<const T_ARRAY,T_INDICES&> Subset(const T_INDICES& indices) const
    {return INDIRECT_ARRAY<const T_ARRAY,T_INDICES&>(Derived(),indices);}

protected:
    T_ARRAY& operator=(const ARRAY_BASE& source)
    {T_ARRAY& self=Derived();RANGE<TV_INT> domain_indices=self.Domain_Indices();const T_ARRAY& source_=source.Derived();assert(domain_indices==source_.Domain_Indices());
    if(!T_ARRAY::Same_Array(self,source_)) for(UNIFORM_ARRAY_ITERATOR<dimension> iterator(domain_indices);iterator.Valid();iterator.Next()) self(iterator.Index())=source_(iterator.Index());
    return self;}

    template<class T_ARRAY1>
    T_ARRAY& operator=(const T_ARRAY1& source)
    {STATIC_ASSERT(CAN_ASSIGN<T,typename T_ARRAY1::ELEMENT>::value);
    T_ARRAY& self=Derived();RANGE<TV_INT> domain_indices=self.Domain_Indices();assert(domain_indices==source.Domain_Indices());
    if(!T_ARRAY::Same_Array(self,source)) for(UNIFORM_ARRAY_ITERATOR<dimension> iterator(domain_indices);iterator.Valid();iterator.Next()) self(iterator.Index())=source(iterator.Index());
    return self;}
};

template<class TV>
class ARRAYS_ND_BASE:public ARRAY_BASE<typename TV::ELEMENT,ARRAYS_ND_BASE<TV>,VECTOR<int,TV::dimension> >
{
    struct UNUSABLE{};
public:
    typedef typename TV::ELEMENT T;
    enum WORKAROUND {d=TV::dimension};
    typedef T ELEMENT;typedef T& RESULT_TYPE;typedef const T& CONST_RESULT_TYPE;
    typedef typename SCALAR_POLICY<T>::TYPE SCALAR;
    template<class T2> struct REBIND{typedef ARRAYS_ND_BASE<VECTOR<T2,d> > TYPE;};
    typedef VECTOR<int,d> TV_INT;
    typedef ARRAY_BASE<T,ARRAYS_ND_BASE<TV>,TV_INT> BASE;
    typedef TV_INT INDEX;

    RANGE<TV_INT> domain;
    TV_INT counts;
    ARRAY_VIEW<T> array; // one-dimensional data storage
private:
    template<class S> struct ELEMENT_OF{typedef typename S::ELEMENT TYPE;};
    typedef typename IF<IS_VECTOR<T>::value,ELEMENT_OF<T>,FIRST<UNUSABLE> >::TYPE::TYPE ELEMENT_OF_T;

protected:
    T* base_pointer;

    ARRAYS_ND_BASE()
        :domain(TV_INT::All_Ones_Vector(),TV_INT()),array(0,0)
    {}

    ARRAYS_ND_BASE(const RANGE<TV_INT>& domain_input)
        :domain(domain_input),counts(domain.Edge_Lengths()+1),array(0,0)
    {}

public:
    void Calculate_Acceleration_Constants()
    {T* array_pointer=array.Get_Array_Pointer();
    base_pointer=array_pointer?array_pointer-Compute_Index(domain.min_corner):0;}

    template<class T_ARRAY2>
    ARRAYS_ND_BASE& operator=(const T_ARRAY2& source)
    {assert(Equal_Dimensions(*this,source));array=source.array;return *this;}

    const RANGE<TV_INT>& Domain_Indices() const
    {return domain;}

    const TV_INT Size() const
    {return counts;}

    bool operator==(const ARRAYS_ND_BASE& v) const
    {return Equal_Dimensions(*this,v) && array==v.array;}

    bool operator!=(const ARRAYS_ND_BASE& v) const
    {return !(*this==v);}

    ARRAYS_ND_BASE& operator+=(const ARRAYS_ND_BASE& v)
    {assert(Equal_Dimensions(*this,v));array+=v.array;return *this;}

    ARRAYS_ND_BASE& operator+=(const T& a)
    {array+=a;return *this;}

    ARRAYS_ND_BASE& operator-=(const ARRAYS_ND_BASE& v)
    {assert(Equal_Dimensions(*this,v));array-=v.array;return *this;}

    ARRAYS_ND_BASE& operator-=(const T& a)
    {array-=a;return *this;}

    template<class T2>
    ARRAYS_ND_BASE& operator*=(const ARRAYS_ND_BASE<VECTOR<T2,d> >& v)
    {assert(Equal_Dimensions(*this,v));array*=v.array;return *this;}

    template<class T2> typename ENABLE_IF<OR<IS_SCALAR<T2>::value,IS_SAME<T,T2>::value>::value,ARRAYS_ND_BASE&>::TYPE
    operator*=(const T2 a)
    {array*=a;return *this;}

    template<class T2>
    ARRAYS_ND_BASE& operator/=(const T2 a)
    {return *this*=Inverse(a);}

    int Number_True() const
    {return array.Number_True();}

    int Pack_Size() const
    {return sizeof(T);}

    template<class T_ARRAY1>
    void Pack(T_ARRAY1& buffer,typename T_ARRAY1::INDEX& position,const TV_INT& p)
    {*(T*)(&buffer(position+1))=(*this)(p);position+=sizeof(T);}
    
    template<class T_ARRAY1>
    void Unpack(const T_ARRAY1& buffer,typename T_ARRAY1::INDEX& position,const TV_INT& p)
    {(*this)(p)=*(const T*)(&buffer(position+1));position+=sizeof(T);}

    void Fill(const T& constant)
    {ARRAYS_COMPUTATIONS::Fill(array,constant);}

    static void Copy(const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy)
    {assert(old_copy.Domain_Indices()==new_copy.Domain_Indices());
    ARRAY_VIEW<T>::Copy(old_copy.array,new_copy.array);}

    template<class T2>
    static void Copy(const T2 constant,const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy)
    {assert(old_copy.Domain_Indices()==new_copy.Domain_Indices());
    new_copy.array=constant*old_copy.array;}

    template<class T2>
    static void Copy(const T2 c1,const ARRAYS_ND_BASE& v1,const ARRAYS_ND_BASE& v2,ARRAYS_ND_BASE& result)
    {assert(Equal_Dimensions(v1,v2)&&Equal_Dimensions(v2,result));
    result.array=c1*v1.array+v2.array;}

    template<class T2>
    static void Copy(const T2 c1,const ARRAYS_ND_BASE& v1,const T2 c2,const ARRAYS_ND_BASE& v2,ARRAYS_ND_BASE& result)
    {assert(Equal_Dimensions(v1,v2)&&Equal_Dimensions(v2,result));
    result.array=c1*v1.array+c2*v2.array;}

    template<class T2>
    static void Copy(const T2 c1,const ARRAYS_ND_BASE& v1,const T2 c2,const ARRAYS_ND_BASE& v2,const T2 c3,const ARRAYS_ND_BASE& v3,ARRAYS_ND_BASE& result)
    {assert(Equal_Dimensions(v1,v2)&&Equal_Dimensions(v2,v3)&&Equal_Dimensions(v3,result));
    result.array=c1*v1.array+c2*v2.array+c3*v3.array;}

    void Clamp_Below(const T& value)
    {array.Clamp_Below(value);}

    T Average() const
    {return ARRAYS_COMPUTATIONS::Average(array);}

    T Max() const
    {return ARRAYS_COMPUTATIONS::Max(array);}

    T Maxabs() const
    {return ARRAYS_COMPUTATIONS::Maxabs(array);}

    T Maxmag() const
    {return ARRAYS_COMPUTATIONS::Maxmag(array);}

    T Min() const
    {return ARRAYS_COMPUTATIONS::Min(array);}

    T Minmag() const
    {return ARRAYS_COMPUTATIONS::Minmag(array);}

    T Sum() const
    {return ARRAYS_COMPUTATIONS::Sum(array);}

    T Sumabs() const
    {return ARRAYS_COMPUTATIONS::Sumabs(array);}

    T Componentwise_Maxabs() const
    {return ARRAYS_COMPUTATIONS::Componentwise_Maxabs(array);}

    static T Dot_Product(const ARRAYS_ND_BASE& a1,const ARRAYS_ND_BASE& a2)
    {assert(Equal_Dimensions(a1,a2));
    return ARRAY_VIEW<T>::Dot_Product(a1.array,a2.array);}

    template<class TV2>
    static typename SCALAR_POLICY<TV>::TYPE Maximum_Magnitude(const ARRAYS_ND_BASE<VECTOR<TV2,d> >& a)
    {STATIC_ASSERT(IS_SAME<T,TV2>::value);typedef typename TV2::SCALAR TS;
    return ARRAY_VIEW<T>::Maximum_Magnitude(a.array);}

private:
    int Compute_Index(const VECTOR<int,3>& index) const
    {return (index.x*counts.y+index.y)*counts.z+index.z;}

    int Compute_Index(const VECTOR<int,2>& index) const
    {return index.x*counts.y+index.y;}

    int Compute_Index(const VECTOR<int,1>& index) const
    {return index.x;}
public:

    T& operator()(const int i,const int j,const int ij)
    {STATIC_ASSERT(d==3);assert(domain.Lazy_Inside(VECTOR<int,3>(i,j,ij)));return base_pointer[(i*counts.y+j)*counts.z+ij];}

    const T& operator()(const int i,const int j,const int ij) const
    {STATIC_ASSERT(d==3);assert(domain.Lazy_Inside(VECTOR<int,3>(i,j,ij)));return base_pointer[(i*counts.y+j)*counts.z+ij];}

    T& operator()(const int i,const int j)
    {STATIC_ASSERT(d==2);assert(domain.Lazy_Inside(VECTOR<int,2>(i,j)));return base_pointer[i*counts.y+j];}

    const T& operator()(const int i,const int j) const
    {STATIC_ASSERT(d==2);assert(domain.Lazy_Inside(VECTOR<int,2>(i,j)));return base_pointer[i*counts.y+j];}

    T& operator()(const int i)
    {STATIC_ASSERT(d==1);assert(domain.Lazy_Inside(VECTOR<int,1>(i)));return base_pointer[i];}

    const T& operator()(const int i) const
    {STATIC_ASSERT(d==1);assert(domain.Lazy_Inside(VECTOR<int,1>(i)));return base_pointer[i];}

    T& operator()(const TV_INT& index)
    {assert(domain.Lazy_Inside(index));return base_pointer[Compute_Index(index)];}

    const T& operator()(const TV_INT& index) const
    {assert(domain.Lazy_Inside(index));return base_pointer[Compute_Index(index)];}

    bool Valid_Index(const TV_INT& index) const
    {return domain.Lazy_Inside(index);}

    bool Valid_Index(const int i,const int j,const int ij) const
    {STATIC_ASSERT(d==3);return domain.Lazy_Inside(TV_INT(i,j,ij));}

    bool Valid_Index(const int i,const int j) const
    {STATIC_ASSERT(d==2);return domain.Lazy_Inside(TV_INT(i,j));}

    bool Valid_Index(const int i) const
    {STATIC_ASSERT(d==1);return domain.Lazy_Inside(TV_INT(i));}

    int Standard_Index(const TV_INT& index) const
    {assert(Valid_Index(index));return Compute_Index(index-domain.min_corner)+1;}

    TV_INT Clamp(const TV_INT& i) const
    {return domain.Clamp(i);}

    TV_INT Clamp_End_Minus_One(const TV_INT& i) const
    {return clamp(i,domain.min_corner,domain.max_corner-1);}

    TV_INT Clamp_End_Minus_Two(const TV_INT& i) const
    {return clamp(i,domain.min_corner,domain.max_corner-2);}

    TV_INT Clamp_End_Minus_Three(const TV_INT& i) const
    {return clamp(i,domain.min_corner,domain.max_corner-3);}

    TV_INT Clamp_Interior(const TV_INT& i) const
    {return clamp(i,domain.min_corner+1,domain.max_corner-1);}

    TV_INT Clamp_Interior_End_Minus_One(const TV_INT& i) const
    {return clamp(i,domain.min_corner+1,domain.max_corner-2);}

    void Clamp(int& i,int& j,int& ij) const
    {STATIC_ASSERT(d==3);TV_INT index=Clamp(TV_INT(i,j,ij));i=index.x;j=index.y;ij=index.z;}

    void Clamp_End_Minus_One(int& i,int& j,int& ij) const
    {STATIC_ASSERT(d==3);TV_INT index=Clamp_End_Minus_One(TV_INT(i,j,ij));i=index.x;j=index.y;ij=index.z;}

    void Clamp_End_Minus_Two(int& i,int& j,int& ij) const
    {STATIC_ASSERT(d==3);TV_INT index=Clamp_End_Minus_Two(TV_INT(i,j,ij));i=index.x;j=index.y;ij=index.z;}

    void Clamp_End_Minus_Three(int& i,int& j,int& ij) const
    {STATIC_ASSERT(d==3);TV_INT index=Clamp_End_Minus_Three(TV_INT(i,j,ij));i=index.x;j=index.y;ij=index.z;}

    void Clamp_Interior(int& i,int& j,int& ij) const
    {STATIC_ASSERT(d==3);TV_INT index=Clamp_Interior(TV_INT(i,j,ij));i=index.x;j=index.y;ij=index.z;}

    void Clamp_Interior_End_Minus_One(int& i,int& j,int& ij) const
    {STATIC_ASSERT(d==3);TV_INT index=Clamp_Interior_End_Minus_One(TV_INT(i,j,ij));i=index.x;j=index.y;ij=index.z;}

    template<class T2>
    static bool Equal_Dimensions(const ARRAYS_ND_BASE& a,const ARRAYS_ND_BASE<VECTOR<T2,d> >& b)
    {return a.domain==b.domain;}

    static bool Equal_Dimensions(const ARRAYS_ND_BASE& a,const int m_start,const int m_end,const int n_start,const int n_end,const int mn_start,const int mn_end)
    {STATIC_ASSERT(d==3);return a.domain==RANGE<TV_INT>(m_start,m_end,n_start,n_end,mn_start,mn_end);}

    // note that these functions move the *contents* of the grid, not the grid itself, being careful about the order of grid traversal.
    static bool Equal_Dimensions(const ARRAYS_ND_BASE& a,const int m_start,const int m_end,const int n_start,const int n_end)
    {STATIC_ASSERT(d==2);return a.domain==RANGE<TV_INT>(m_start,m_end,n_start,n_end);}

    // note that these functions move the *contents* of the grid, not the grid itself, being careful about the order of grid traversal.
    static bool Equal_Dimensions(const ARRAYS_ND_BASE& a,const int m_start,const int m_end)
    {STATIC_ASSERT(d==1);return a.domain==RANGE<TV_INT>(m_start,m_end);}

    static void Extract_Dimension(const ARRAYS_ND_BASE& old_array,ARRAYS_ND_BASE<VECTOR<ELEMENT_OF_T,d> >& extracted_array,int dim)
    {STATIC_ASSERT(IS_VECTOR<T>::value);assert(Equal_Dimensions(old_array,extracted_array));//extracted_array.Resize(old_array.domain,false,false);
    for(int i=1;i<=old_array.array.m;i++) extracted_array.array(i)=old_array.array(i)(dim);}

    static void Get(ARRAYS_ND_BASE& new_copy,const ARRAYS_ND_BASE& old_copy)
    {if(&old_copy!=&new_copy) Put(old_copy,new_copy);}

    static void Shifted_Get(ARRAYS_ND_BASE& new_copy,const ARRAYS_ND_BASE& old_copy,const TV_INT& shift)
    {Shifted_Put(old_copy,new_copy,shift);}

    static void Limited_Shifted_Get(ARRAYS_ND_BASE& new_copy,const ARRAYS_ND_BASE& old_copy,const VECTOR<int,3>& shift)
    {STATIC_ASSERT(d==3);
    RANGE<TV_INT> new_domain=new_copy.Domain_Indices(),old_domain=old_copy.Domain_Indices();
    RANGE<TV_INT> box(TV_INT::Componentwise_Max(new_domain.min_corner,old_domain.min_corner-shift),TV_INT::Componentwise_Min(new_domain.max_corner,old_domain.max_corner-shift));
    TV_INT i;
    for(i.x=box.min_corner.x;i.x<=box.max_corner.x;i.x++)
        for(i.y=box.min_corner.y;i.y<=box.max_corner.y;i.y++)
            for(i.z=box.min_corner.z;i.z<=box.max_corner.z;i.z++)
                new_copy(i)=old_copy(i+shift);}

    static void Limited_Shifted_Get(ARRAYS_ND_BASE& new_copy,const ARRAYS_ND_BASE& old_copy,const VECTOR<int,2>& shift)
    {STATIC_ASSERT(d==2);
    RANGE<TV_INT> new_domain=new_copy.Domain_Indices(),old_domain=old_copy.Domain_Indices();
    RANGE<TV_INT> box(TV_INT::Componentwise_Max(new_domain.min_corner,old_domain.min_corner-shift),TV_INT::Componentwise_Min(new_domain.max_corner,old_domain.max_corner-shift));
    TV_INT i;
    for(i.x=box.min_corner.x;i.x<=box.max_corner.x;i.x++)
        for(i.y=box.min_corner.y;i.y<=box.max_corner.y;i.y++)
            new_copy(i)=old_copy(i+shift);}

    static void Limited_Shifted_Get(ARRAYS_ND_BASE& new_copy,const ARRAYS_ND_BASE& old_copy,const VECTOR<int,1>& shift)
    {STATIC_ASSERT(d==1);
    RANGE<TV_INT> new_domain=new_copy.Domain_Indices(),old_domain=old_copy.Domain_Indices();
    RANGE<TV_INT> box(TV_INT::Componentwise_Max(new_domain.min_corner,old_domain.min_corner-shift),TV_INT::Componentwise_Min(new_domain.max_corner,old_domain.max_corner-shift));
    TV_INT i;
    for(i.x=box.min_corner.x;i.x<=box.max_corner.x;i.x++)
        new_copy(i)=old_copy(i+shift);}

    static void Put(const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy)
    {if(&old_copy!=&new_copy) Put(ONE(),old_copy,new_copy,RANGE<TV_INT>::Intersect(old_copy.Domain_Indices(),new_copy.Domain_Indices()));}

    void Put_With_Range(RANGE<TV_INT>& range,const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy)
    {if(&old_copy!=&new_copy) Put(ONE(),old_copy,new_copy,range);}

    static void Shifted_Put(const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy,const VECTOR<int,3>& shift)
    {if(shift==TV_INT()) Put(old_copy,new_copy);
    else{
        TV_INT i;
        RANGE<TV_INT> new_domain=new_copy.Domain_Indices(),old_domain=old_copy.Domain_Indices();
        for(i.x=new_domain.min_corner.x;i.x<=new_domain.max_corner.x;i.x++)
            for(i.y=new_domain.min_corner.y;i.y<=new_domain.max_corner.y;i.y++)
                for(i.z=new_domain.min_corner.z;i.z<=new_domain.max_corner.z;i.z++)
                    new_copy(i)=old_copy(i+shift);}}

    static void Shifted_Put(const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy,const VECTOR<int,2>& shift)
    {if(shift==TV_INT()) Put(old_copy,new_copy);
    else{
        TV_INT i;
        RANGE<TV_INT> new_domain=new_copy.Domain_Indices(),old_domain=old_copy.Domain_Indices();
        for(i.x=new_domain.min_corner.x;i.x<=new_domain.max_corner.x;i.x++)
            for(i.y=new_domain.min_corner.y;i.y<=new_domain.max_corner.y;i.y++)
                new_copy(i)=old_copy(i+shift);}}

    static void Shifted_Put(const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy,const VECTOR<int,1>& shift)
    {if(shift==TV_INT()) Put(old_copy,new_copy);
    else{
        TV_INT i;
        RANGE<TV_INT> new_domain=new_copy.Domain_Indices(),old_domain=old_copy.Domain_Indices();
        for(i.x=new_domain.min_corner.x;i.x<=new_domain.max_corner.x;i.x++)
            new_copy(i)=old_copy(i+shift);}}

    template<class T2>
    static void Put(const T2 constant,const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy)
    {Put(constant,old_copy,new_copy,old_copy.Domain_Indices());}

private:
    template<class T2>
    static void Put(const T2 constant,const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy,const RANGE<VECTOR<int,3> >& box)
    {assert(old_copy.Domain_Indices().Contains(box));assert(new_copy.Domain_Indices().Contains(box));
    TV_INT i;
    for(i.x=box.min_corner.x;i.x<=box.max_corner.x;i.x++)
        for(i.y=box.min_corner.y;i.y<=box.max_corner.y;i.y++)
            for(i.z=box.min_corner.z;i.z<=box.max_corner.z;i.z++)
                new_copy(i)=constant*old_copy(i);}

    template<class T2>
    static void Put(const T2 constant,const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy,const RANGE<VECTOR<int,2> >& box)
    {assert(old_copy.Domain_Indices().Contains(box));assert(new_copy.Domain_Indices().Contains(box));
    TV_INT i;
    for(i.x=box.min_corner.x;i.x<=box.max_corner.x;i.x++)
        for(i.y=box.min_corner.y;i.y<=box.max_corner.y;i.y++)
            new_copy(i)=constant*old_copy(i);}

    template<class T2>
    static void Put(const T2 constant,const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy,const RANGE<VECTOR<int,1> >& box)
    {assert(old_copy.Domain_Indices().Contains(box));assert(new_copy.Domain_Indices().Contains(box));
    TV_INT i;
    for(i.x=box.min_corner.x;i.x<=box.max_corner.x;i.x++)
        new_copy(i)=constant*old_copy(i);}

    static void Put(const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy,const RANGE<TV_INT>& box)
    {Put(ONE(),old_copy,new_copy,box);}

    template<class T2>
    static void Put(const T2 constant,const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy,const int m_start,const int m_end,const int n_start,const int n_end,const int mn_start,const int mn_end)
    {STATIC_ASSERT(d==3);Put(constant,old_copy,new_copy,RANGE<TV_INT>(m_start,m_end,n_start,n_end,mn_start,mn_end));}

    static void Put(const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy,const int m_start,const int m_end,const int n_start,const int n_end,const int mn_start,const int mn_end)
    {STATIC_ASSERT(d==3);Put(old_copy,new_copy,RANGE<TV_INT>(m_start,m_end,n_start,n_end,mn_start,mn_end));}
public:
    // note that these functions move the *contents* of the grid, not the grid itself, being careful about the order of grid traversal.
    void Move_Contents_By_Offset(const VECTOR<int,3>& offset)
    {STATIC_ASSERT(d==3);TV_INT i,s(TV_INT::Componentwise_Greater_Equal(offset,TV_INT())),c(s*2-1),e(s*domain.Edge_Lengths()),a(domain.max_corner-e),b(domain.min_corner+e);
    for(i.x=a.x;i.x<=b.x;i.x+=c.x) for(i.y=a.y;i.y<=b.y;i.y+=c.y) for(i.z=a.z;i.z<=b.z;i.z+=c.z) (*this)(i)=(*this)(i+offset);}

    void Move_Contents_By_Offset(const VECTOR<int,2>& offset)
    {STATIC_ASSERT(d==2);TV_INT i,s(TV_INT::Componentwise_Greater_Equal(offset,TV_INT())),c(s*2-1),e(s*domain.Edge_Lengths()),a(domain.max_corner-e),b(domain.min_corner+e);
    for(i.x=a.x;i.x<=b.x;i.x+=c.x) for(i.y=a.y;i.y<=b.y;i.y+=c.y) (*this)(i)=(*this)(i+offset);}

    void Move_Contents_By_Offset(const VECTOR<int,1>& offset)
    {STATIC_ASSERT(d==1);TV_INT i,s(TV_INT::Componentwise_Greater_Equal(offset,TV_INT())),c(s*2-1),e(s*domain.Edge_Lengths()),a(domain.max_corner-e),b(domain.min_corner+e);
    for(i.x=a.x;i.x<=b.x;i.x+=c.x) (*this)(i)=(*this)(i+offset);}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    void Print_Grid_Array(VECTOR<int,1> counts,int ghost_cells=0)
    {for(int i=1-ghost_cells;i<=counts.x+ghost_cells;i++) LOG::cout<<(operator()(VECTOR<int,1>(i)));LOG::cout<<std::endl;}

    void Print_Grid_Array(VECTOR<int,2> counts,int ghost_cells=0)
    {for(int i=1-ghost_cells;i<=counts.x+ghost_cells;i++){
        for(int j=1-ghost_cells;j<=counts.y+ghost_cells;j++) LOG::cout<<(operator()(VECTOR<int,2>(i,j)));
        LOG::cout<<std::endl;}
    LOG::cout<<std::endl;}

    void Print_Grid_Array(VECTOR<int,3> counts,int ghost_cells=0)
    {for(int i=1-ghost_cells;i<=counts.x+ghost_cells;i++){
        for(int j=1-ghost_cells;j<=counts.y+ghost_cells;j++){
            for(int k=1-ghost_cells;k<=counts.z+ghost_cells;k++) LOG::cout<<(operator()(VECTOR<int,3>(i,j,k)));
            LOG::cout<<std::endl;}
        LOG::cout<<std::endl;}
    LOG::cout<<std::endl;}
#endif

//#####################################################################
};
}
#endif


//#####################################################################
// Copyright 2004-2009, Kevin Der, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jerry Talton, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_BASE
//#####################################################################
#ifndef __ARRAY_BASE__
#define __ARRAY_BASE__

#include <PhysBAM_Tools/Arrays/ARRAY_DIFFERENCE.h>
#include <PhysBAM_Tools/Arrays/ARRAY_LEFT_MULTIPLE.h>
#include <PhysBAM_Tools/Arrays/ARRAY_NEGATION.h>
#include <PhysBAM_Tools/Arrays/ARRAY_PLUS_SCALAR.h>
#include <PhysBAM_Tools/Arrays/ARRAY_PRODUCT.h>
#include <PhysBAM_Tools/Arrays/ARRAY_SUM.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_ARITHMETIC.h>
#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/SCALAR_POLICY.h>
#include <iostream>
namespace PhysBAM{

template<class T,class T_ARRAY,class ID> class ARRAY_BASE;

template<class T_ARRAY,class ENABLER=void> struct CANONICALIZE_CONST_ARRAY:public FIRST<T_ARRAY>{};

template<class T_ARRAY1,class T_ARRAY2> struct SAME_ARRAY_CANONICAL{static bool Same_Array(const T_ARRAY1& array1,const T_ARRAY2& array2)
{STATIC_ASSERT(!IS_SAME<T_ARRAY1,T_ARRAY2>::value);return false;}};

template<class T_ARRAY> struct SAME_ARRAY_CANONICAL<T_ARRAY,T_ARRAY>{static bool Same_Array(const T_ARRAY& array1,const T_ARRAY& array2)
{return T_ARRAY::Same_Array(array1,array2);}};

template<class TA1,class TA2> struct SAME_ARRAY:public SAME_ARRAY_CANONICAL<typename CANONICALIZE_CONST_ARRAY<TA1>::TYPE,typename CANONICALIZE_CONST_ARRAY<TA2>::TYPE>{};

template<class T,class T_ARRAY,class ID>
class ARRAY_BASE
{
    struct UNUSABLE{};
public:
    typedef T ELEMENT;
    typedef ID INDEX;
    typedef T value_type; // for stl
    typedef T* iterator; // for stl
    typedef const T* const_iterator; // for stl
    typedef int difference_type; // for stl

    typedef T& RESULT_TYPE;
    typedef const T& CONST_RESULT_TYPE;
    typedef typename SCALAR_POLICY<T>::TYPE SCALAR;
    template<class T2> struct REBIND{typedef ARRAY_BASE<T2,typename T_ARRAY::template REBIND<T2>::TYPE,ID> TYPE;};

protected:
    ARRAY_BASE(){}
    ARRAY_BASE(const ARRAY_BASE&){}
    ~ARRAY_BASE(){}
public:

    T_ARRAY& Derived()
    {return static_cast<T_ARRAY&>(*this);}

    const T_ARRAY& Derived() const
    {return static_cast<const T_ARRAY&>(*this);}

protected:
    template<class T_ARRAY2> struct IS_ARRAY_BASE {static const bool value=false;};
    template<class T2,class T_ARRAY2> struct IS_ARRAY_BASE<ARRAY_BASE<T2,T_ARRAY2,ID> > {static const bool value=true;};
public:

    template<class T_ARRAY1> static bool Same_Array(const T_ARRAY1& array1,const T_ARRAY1& array2)
    {return &array1==&array2;}

    template<class T_ARRAY1,class T_ARRAY2> static bool
    Same_Array(const T_ARRAY1& array1,const T_ARRAY2& array2)
    {return SAME_ARRAY<T_ARRAY1,T_ARRAY2>::Same_Array(array1,array2);}

protected:
    T_ARRAY& operator=(const ARRAY_BASE& source)
    {T_ARRAY& self=Derived();ID m=self.Size();const T_ARRAY& source_=source.Derived();assert(m==source_.Size());
    if(!T_ARRAY::Same_Array(self,source_)) for(ID i(1);i<=m;i++) self(i)=source_(i);
    return self;}

    template<class T_ARRAY1>
    T_ARRAY& operator=(const T_ARRAY1& source)
    {STATIC_ASSERT(CAN_ASSIGN<T,typename T_ARRAY1::ELEMENT>::value);
    T_ARRAY& self=Derived();ID m=self.Size();assert(m==source.Size());
    if(!T_ARRAY::Same_Array(self,source)) for(ID i(1);i<=m;i++) self(i)=source(i);
    return self;}
public:

    template<class T_INDICES>
    INDIRECT_ARRAY<T_ARRAY,T_INDICES&> Subset(const T_INDICES& indices)
    {return INDIRECT_ARRAY<T_ARRAY,T_INDICES&>(Derived(),indices);}

    template<class T_INDICES>
    INDIRECT_ARRAY<const T_ARRAY,T_INDICES&> Subset(const T_INDICES& indices) const
    {return INDIRECT_ARRAY<const T_ARRAY,T_INDICES&>(Derived(),indices);}

    INDIRECT_ARRAY<T_ARRAY,IDENTITY_ARRAY<> > Prefix(const ID prefix_size)
    {assert(prefix_size<=Derived().Size());return INDIRECT_ARRAY<T_ARRAY,IDENTITY_ARRAY<> >(Derived(),IDENTITY_ARRAY<>(prefix_size));}

    INDIRECT_ARRAY<const T_ARRAY,IDENTITY_ARRAY<> > Prefix(const ID prefix_size) const
    {assert(prefix_size<=Derived().Size());return INDIRECT_ARRAY<const T_ARRAY,IDENTITY_ARRAY<> >(Derived(),IDENTITY_ARRAY<>(prefix_size));}

private:
    typedef typename IF<IS_CLASS<T>::value,T,UNUSABLE>::TYPE T_IF_CLASS;
public:

    template<class T_FIELD,T_FIELD T_IF_CLASS::* field>
    PROJECTED_ARRAY<T_ARRAY,FIELD_PROJECTOR<T_IF_CLASS,T_FIELD,field> > Project()
    {return PROJECTED_ARRAY<T_ARRAY,FIELD_PROJECTOR<ELEMENT,T_FIELD,field> >(Derived());}

    template<class T_FIELD,T_FIELD T_IF_CLASS::* field>
    PROJECTED_ARRAY<const T_ARRAY,FIELD_PROJECTOR<T_IF_CLASS,T_FIELD,field> > Project() const
    {return PROJECTED_ARRAY<const T_ARRAY,FIELD_PROJECTOR<ELEMENT,T_FIELD,field> >(Derived());}

    PROJECTED_ARRAY<T_ARRAY,INDEX_PROJECTOR> Project(const ID index)
    {return PROJECTED_ARRAY<T_ARRAY,INDEX_PROJECTOR>(Derived(),INDEX_PROJECTOR(index));}

    PROJECTED_ARRAY<const T_ARRAY,INDEX_PROJECTOR> Project(const ID index) const
    {return PROJECTED_ARRAY<const T_ARRAY,INDEX_PROJECTOR>(Derived(),INDEX_PROJECTOR(index));}

private:
    template<class S> struct ELEMENT_OF{typedef typename S::ELEMENT TYPE;};
    typedef typename IF<IS_VECTOR<T>::value,ELEMENT_OF<T>,FIRST<UNUSABLE> >::TYPE::TYPE ELEMENT_OF_T;
public:

    T& operator()(const ID i)
    {return Derived()(i);}

    const T& operator()(const ID i) const
    {return Derived()(i);}

    bool Valid_Index(const ID i) const
    {return ID(1)<=i && i<=Size();}

    ARRAY_VIEW<ELEMENT_OF_T> Flattened() // valid only for contiguous arrays of VECTOR<T,d>
    {T_ARRAY& self=Derived();return ARRAY_VIEW<typename T::ELEMENT>(T::m*self.Size(),self.Get_Array_Pointer()->begin());}

    ARRAY_VIEW<const ELEMENT_OF_T> Flattened() const // valid only for contiguous arrays of VECTOR<T,d>
    {const T_ARRAY& self=Derived();return ARRAY_VIEW<const typename T::ELEMENT>(T::m*self.Size(),self.Get_Array_Pointer()->begin());}

    template<class ID2>
    ARRAY_VIEW<const T,ID2> Array_View(const ID first,const ID2 length) const
    {const T_ARRAY& self=Derived();assert(Value(first)+Value(length)-1<=Value(self.Size()) && Value(length)>=0);return ARRAY_VIEW<const T,ID2>(length,self.Get_Array_Pointer()+Value(first-1));}

    template<class ID2>
    ARRAY_VIEW<T,ID2> Array_View(const ID first,const ID2 length)
    {T_ARRAY& self=Derived();assert(Value(first)+Value(length)-1<=Value(self.Size()) && Value(length)>=0);return ARRAY_VIEW<T,ID2>(length,self.Get_Array_Pointer()+Value(first-1));}

    template<class T_ARRAY1>
    bool operator==(const T_ARRAY1& v) const
    {STATIC_ASSERT_SAME(T,typename T_ARRAY1::ELEMENT);
    const T_ARRAY& self=Derived();ID m=self.Size();
    if(m!=v.Size()) return false;for(ID i(1);i<=m;i++) if(self(i)!=v(i)) return false;return true;}

    template<class T_ARRAY1>
    bool operator!=(const T_ARRAY1& v) const
    {return !(*this==v);}

    template<class T_ARRAY1>
    T_ARRAY& operator+=(const ARRAY_BASE<T,T_ARRAY1,ID>& v)
    {return ARRAYS_COMPUTATIONS::Plus_Equals(*this,v);}

    T_ARRAY& operator+=(const T& a)
    {return ARRAYS_COMPUTATIONS::Plus_Equals(*this,a);}

    template<class T_ARRAY1>
    T_ARRAY& operator-=(const ARRAY_BASE<T,T_ARRAY1,ID>& v)
    {return ARRAYS_COMPUTATIONS::Minus_Equals(*this,v);}

    T_ARRAY& operator-=(const T& a)
    {return ARRAYS_COMPUTATIONS::Minus_Equals(*this,a);}

    template<class T2,class T_ARRAY_T2>
    T_ARRAY& operator*=(const ARRAY_BASE<T2,T_ARRAY_T2,ID>& v)
    {return ARRAYS_COMPUTATIONS::Times_Equals(*this,v);}

    T_ARRAY& operator*=(const SCALAR& a)
    {return ARRAYS_COMPUTATIONS::Times_Equals(*this,a);}

    template<class T2,class T_ARRAY_T2>
    T_ARRAY& operator/=(const ARRAY_BASE<T2,T_ARRAY_T2,ID>& v)
    {return ARRAYS_COMPUTATIONS::Divide_Equals(*this,v);}

    T_ARRAY& operator/=(const SCALAR& a)
    {return ARRAYS_COMPUTATIONS::Divide_Equals(*this,a);}

    T& Last()
    {T_ARRAY& self=Derived();return self(self.Size());}

    const T& Last() const
    {const T_ARRAY& self=Derived();return self(self.Size());}

    ID Size() const
    {return Derived().Size();}

    ID Find(const T& element) const
    {const T_ARRAY& self=Derived();ID m=self.Size();
    for(ID i(1);i<=m;i++) if(self(i)==element) return i;return 0;}

    bool Find(const T& element,ID& index) const // returns the first occurence of an element in an array
    {return Find(element,1,index);}

    bool Find(const T& element,const ID start_index,ID& index) const // returns the first occurence after start_index of an element in an array
    {const T_ARRAY& self=Derived();ID m=self.Size();
    for(ID i=start_index;i<=m;i++) if(self(i)==element){index=i;return true;}return false;}

    bool Contains(const T& element) const
    {const T_ARRAY& self=Derived();ID m=self.Size();
    for(ID i(1);i<=m;i++) if(self(i)==element) return true;return false;}

    bool Contains_Only(const T& element) const
    {const T_ARRAY& self=Derived();ID m=self.Size();
    for(ID i(1);i<=m;i++) if(self(i)!=element) return false;return true;}

    int Count_Matches(const T& value) const
    {const T_ARRAY& self=Derived();ID m=self.Size();
    int count=0;for(ID i(1);i<=m;i++) if(self(i)==value) count++;return count;}

    int Number_True() const
    {STATIC_ASSERT_SAME(T,bool);return Count_Matches(true);}

    int Number_False() const
    {STATIC_ASSERT_SAME(T,bool);return Count_Matches(false);}

    void Get_Unique(ARRAY<T>& array) const
    {const T_ARRAY& self=Derived();HASHTABLE<T> hash(Value(self.Size())*3/2);hash.Set_All(self);hash.Get_Keys(array);}

    void Prune_Duplicates()
    {T_ARRAY& self=Derived();HASHTABLE<T> hash(Value(self.Size())*3/2);int j=0;for(int i=1;i<=self.Size();i++) if(hash.Set(self(i))) self(++j)=self(i);self.Resize(j);}

    void Coalesce()
    {Sort(*this);T_ARRAY& self=Derived();int j=0;if(self.Size()>0) j=1;for(int i=2;i<=self.Size();i++){if(!(self(j)<self(i))) self(j).Merge(self(i));else self(++j)=self(i);}self.Resize(j);}

    template<class T_ARRAY1,class T_ARRAY2>
    static void Copy(const T_ARRAY1& old_copy,T_ARRAY2& new_copy)
    {new_copy=old_copy;}

    template<class T2,class T_ARRAY1,class T_ARRAY2>
    static void Copy(const T2 constant,const T_ARRAY1& array,T_ARRAY2& result)
    {result=constant*array;}

    template<class T2,class T_ARRAY1,class T_ARRAY2,class T_ARRAY3>
    static void Copy(const T2 c1,const T_ARRAY1& v1,const T_ARRAY2& v2,T_ARRAY3& result)
    {result=c1*v1+v2;}

    template<class T2,class T_ARRAY1,class T_ARRAY2,class T_ARRAY3>
    static void Copy(const T2 c1,const T_ARRAY1& v1,const T2 c2,const T_ARRAY2& v2,T_ARRAY3& result)
    {result=c1*v1+c2*v2;}

    template<class T2,class T_ARRAY1,class T_ARRAY2,class T_ARRAY3,class T_ARRAY4>
    static void Copy(const T2 c1,const T_ARRAY1& v1,const T2 c2,const T_ARRAY2& v2,const T2 c3,const T_ARRAY3& v3,T_ARRAY4& result)
    {result=c1*v1+c2*v2+c3*v3;}

    static void Copy_Element(const ID from,const ID to)
    {T_ARRAY& self=Derived();self(to)=self(from);}
    
    template<class T_ARRAY1>
    static void Copy_Element(const T_ARRAY1& from_array,const ID from,const ID to)
    {T_ARRAY& self=Derived();self(to)=static_cast<const T_ARRAY&>(from_array)(from);}

    static void Get(T_ARRAY& new_copy,const T_ARRAY& old_copy)
    {if(&old_copy!=&new_copy) new_copy=old_copy.Prefix(new_copy.Size());}

    static void Put(const T_ARRAY& old_copy,T_ARRAY& new_copy)
    {if(&old_copy!=&new_copy) new_copy.Prefix(old_copy.Size())=old_copy;}

    template<class T2>
    static void Put(const T2 constant,const T_ARRAY& old_copy,T_ARRAY& new_copy)
    {new_copy.Prefix(old_copy.Size())=constant*old_copy;}

    void Clamp_Below(const T& value)
    {T_ARRAY& self=Derived();ID m=self.Size();for(ID i(1);i<=m;i++) self(i)=clamp_min(self(i),value);}

    static void Find_Common_Elements(const T_ARRAY& a,const T_ARRAY& b,T_ARRAY& result)
    {assert(&a!=&result);assert(&b!=&result);result.Remove_All();
    ID m=a.Size();for(ID i(1);i<=m;i++) if(b.Contains(a(i))) result.Append(a(i));}

    template<class T_ARRAY1,class T_ARRAY2>
    static bool Equal_Dimensions(const T_ARRAY1& a,const T_ARRAY2& b)
    {return a.Size()==b.Size();}

    template<class T_ARRAY1,class T_ARRAY_INT>
    static void Permute(const T_ARRAY1& source,T_ARRAY1& destination,const T_ARRAY_INT& permutation)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY1::ELEMENT);
    STATIC_ASSERT_SAME(ID,typename T_ARRAY_INT::ELEMENT);
    ID m=permutation.Size();for(ID i(1);i<=m;i++) destination(i)=source(permutation(i));}

    template<class T_ARRAY1,class T_ARRAY_INT>
    static void Unpermute(const T_ARRAY1& source,T_ARRAY1& destination,const T_ARRAY_INT& permutation)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY1::ELEMENT);
    STATIC_ASSERT_SAME(ID,typename T_ARRAY_INT::ELEMENT);
    ID m=permutation.Size();for(ID i(1);i<=m;i++) destination(permutation(i))=source(i);}

    template<class T_ARRAY1>
    void Remove_Sorted_Indices(const T_ARRAY1& index)
    {STATIC_ASSERT_SAME(ID,typename T_ARRAY1::ELEMENT);
    T_ARRAY& self=Derived();ID m=self.Size(),index_m=index.Size();
    if(index_m==0) return;
    for(ID kk(1);kk<=index_m-1;kk++){
        assert(1<=index(kk) && index(kk)<=m);
        for(ID i=index(kk)+1-kk;i<=index(kk+1)-1-kk;i++) self(i)=self(i+kk);}
    for(ID i=index(index_m)+1-index_m;i<=m-index_m;i++) self(i)=self(i+index_m);
    self.Resize(m-index_m);}

    template<class T_ARRAY1>
    void Remove_Sorted_Indices_Lazy(const T_ARRAY1& index)
    {STATIC_ASSERT_SAME(ID,typename T_ARRAY1::ELEMENT);
    T_ARRAY& self=Derived();ID index_m=index.Size();
    if(index_m==0) return;
    ID curr=0;
    for(ID k=index_m;k>=ID(1);k--)if(index(k)!=curr){curr=index(k);self.Remove_Index_Lazy(curr);}
    self.Compact();}

    int Pack_Size() const
    {return sizeof(T);}

    template<class T_ARRAY1>
    void Pack(T_ARRAY1& buffer,typename T_ARRAY1::INDEX& position,const ID p)
    {T_ARRAY& self=Derived();*(T*)(&buffer(position+1))=self(p);position+=sizeof(T);}
    
    template<class T_ARRAY1>
    void Unpack(const T_ARRAY1& buffer,typename T_ARRAY1::INDEX& position,const ID p)
    {T_ARRAY& self=Derived();self(p)=*(const T*)(&buffer(position+1));position+=sizeof(T);}

    T* begin() // for stl
    {return Derived().Get_Array_Pointer();}

    const T* begin() const // for stl
    {return Derived().Get_Array_Pointer();}

    T* end() // for stl
    {return Derived().Get_Array_Pointer()+Derived().Size();}

    const T* end() const // for stl
    {return Derived().Get_Array_Pointer()+Derived().Size();}

    ID Binary_Search(const T& value) const// lower_bound binary search
    {const T_ARRAY& self=Derived();
    ID first=1,last=self.Size()+1,length=last-first,half,middle;    
    while(length>ID(0)){
        half=length>>1;
        middle=first+half;
        if(self(middle)<value){
            first=middle+1;
            length-=half+1;}
        else length=half;}
    return first;}

//#####################################################################
};
template<class T,class T_ARRAY,class ID>
inline std::ostream& operator<<(std::ostream& output,const ARRAY_BASE<T,T_ARRAY,ID>& a)
{output<<"(";
const T_ARRAY& a_=a.Derived();
ID m=a_.Size();
for(ID i(1);i<=m;i++){
    output<<a_(i);
    if(i<m) output<<" ";}
output<<")";
return output;}
//#####################################################################
template<class T_ARRAY1,class T_ARRAY2> struct CAN_ASSIGN<T_ARRAY1,T_ARRAY2,typename ENABLE_IF<IS_ARRAY<T_ARRAY1>::value && IS_ARRAY<T_ARRAY2>::value && IS_SAME<typename T_ARRAY1::ELEMENT,typename T_ARRAY2::ELEMENT>::value && !IS_SAME<T_ARRAY1,T_ARRAY2>::value>::TYPE>
{static const bool value=true;};
}
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#endif

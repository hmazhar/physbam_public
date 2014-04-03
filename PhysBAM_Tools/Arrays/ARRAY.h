//#####################################################################
// Copyright 2004-2007, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY
//#####################################################################
#ifndef __ARRAY__
#define __ARRAY__

#include <PhysBAM_Tools/Arrays/ARRAY_BASE.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
namespace PhysBAM{

template<class T,class ID> class ARRAY;
template<class T,class ID> struct IS_ARRAY<ARRAY<T,ID> > {static const bool value=true;};
template<class T,class ID> struct IS_ARRAY_VIEW<ARRAY<T,ID> > {static const bool value=false;};

template<class T,class ID> struct CANONICALIZE_CONST_ARRAY<ARRAY<T,ID> >:public FIRST<ARRAY_VIEW<const T,ID> >{};

template<class T,class ID>
class ARRAY:public ARRAY_BASE<T,ARRAY<T,ID>,ID>
{
public:
    template<class T2> struct REBIND{typedef ARRAY<T2,ID> TYPE;};
    template<int length> struct REBIND_LENGTH:public PhysBAM::REBIND_LENGTH<ARRAY,length>{};
    typedef T ELEMENT;typedef ID INDEX;
private:
    struct UNUSABLE{};
public:
    T* base_pointer;
    ID buffer_size;

    ID m; // the current size of the array (buffer_size may be larger for elbow room)

    ARRAY()
        :base_pointer(0),buffer_size(0),m(0)
    {}

    explicit ARRAY(const ID m_input,const bool initialize_using_default_constructor=true)
        :base_pointer(0),buffer_size(m_input),m(m_input)
    {
        assert(m>=ID());base_pointer=new T[Value(m)];
        if(!IS_CLASS<T>::value && initialize_using_default_constructor){ID m=Size();for(ID i(1);i<=m;i++) (*this)(i)=T();}
    }

    ARRAY(const ARRAY& array)
        :base_pointer(0),buffer_size(array.m),m(array.m)
    {
        base_pointer=new T[Value(m)];
        for(int i=0;i<Value(m);i++) base_pointer[i]=array.base_pointer[i];
    }

    template<class T_ARRAY>
    explicit ARRAY(const T_ARRAY& array,typename ENABLE_IF<IS_SAME<T,typename T_ARRAY::ELEMENT>::value,UNUSABLE>::TYPE unused=UNUSABLE())
        :base_pointer(0),buffer_size(array.Size()),m(array.Size())
    {
        base_pointer=new T[Value(m)];
        for(ID i(1);i<=m;i++) (*this)(i)=array(i);
    }

    ~ARRAY()
    {
        delete[] base_pointer;
    }

    ARRAY& operator=(const ARRAY& source)
    {ID source_m=source.m;
    if(buffer_size<source_m) Resize_Helper(source_m,false,false);
    else if(Same_Array(*this,source)) return *this;
    m=source_m;
    for(ID i(1);i<=source_m;i++) (*this)(i)=source(i);
    return *this;}

    template<class T_ARRAY>
    ARRAY& operator=(const T_ARRAY& source)
    {STATIC_ASSERT(CAN_ASSIGN<T,typename T_ARRAY::ELEMENT>::value);
    ID source_m=source.Size();
    if(buffer_size<source_m) Resize_Helper(source_m,false,false);
    else if(Same_Array(*this,source)) return *this;
    m=source_m;
    for(ID i(1);i<=source_m;i++) (*this)(i)=source(i);
    return *this;}

    ID Size() const
    {return m;}

    T& operator()(const ID i)
    {assert(i>=ID(1) && i<=m);return base_pointer[Value(i)-1];}

    const T& operator()(const ID i) const
    {assert(i>=ID(1) && i<=m);return base_pointer[Value(i)-1];}

    bool Valid_Index(const ID i) const
    {return ID(1)<=i && i<=m;}

    T* Get_Array_Pointer()
    {return base_pointer;}

    const T* Get_Array_Pointer() const
    {return base_pointer;}

    ID Max_Size() const
    {return buffer_size;}

    void Compact()
    {if(m<buffer_size) Exact_Resize(m);}

    void Fill(T value)
    {ARRAYS_COMPUTATIONS::Fill(*this,value);}

private:
    void Resize_Helper(const ID buffer_new,const bool initialize_new_elements=true,const bool copy_existing_elements=true)
    {if(buffer_size==buffer_new) return;assert(m<=buffer_size);
    T* p=new T[Value(buffer_new)];
    int m_end=Value(PhysBAM::min(m,buffer_new));
    if(copy_existing_elements) for(int i=0;i<m_end;i++) p[i]=base_pointer[i];
    if(!IS_CLASS<T>::value && initialize_new_elements) for(int i=m_end;i<Value(buffer_new);i++) p[i]=T();
    delete[] base_pointer;
    base_pointer=p;
    buffer_size=buffer_new;}

    void Resize_Helper(const ID buffer_new,const bool initialize_new_elements,const bool copy_existing_elements,const T& initialization_value)
    {if(buffer_size==buffer_new) return;assert(m<=buffer_size);
    T* p=new T[Value(buffer_new)];
    int m_end=Value(PhysBAM::min(m,buffer_new));
    if(copy_existing_elements) for(int i=0;i<m_end;i++) p[i]=base_pointer[i];
    if(!IS_CLASS<T>::value && initialize_new_elements) for(int i=m_end;i<Value(buffer_new);i++) p[i]=initialization_value;
    delete[] base_pointer;
    base_pointer=p;
    buffer_size=buffer_new;}

    void Ensure_Enough_Space(const ID m_new,const bool copy_existing_elements=true) PHYSBAM_ALWAYS_INLINE
    {if(buffer_size<m_new) Resize_Helper(ID(4*Value(m_new)/3+2),false,copy_existing_elements);}

public:
    void Preallocate(const ID max_size)
    {if(buffer_size<max_size) Resize_Helper(max_size,false);}

    void Resize(const ID m_new,const bool initialize_new_elements=true,const bool copy_existing_elements=true)
    {Ensure_Enough_Space(m_new,copy_existing_elements);if(initialize_new_elements && m_new>m) for(int i=Value(m);i<Value(m_new);i++) base_pointer[i]=T();
    m=m_new;}

    void Resize(const ID m_new,const bool initialize_new_elements,const bool copy_existing_elements,const T& initialization_value)
    {Ensure_Enough_Space(m_new,copy_existing_elements);if(initialize_new_elements && m_new>m) for(int i=Value(m);i<Value(m_new);i++) base_pointer[i]=initialization_value;
    m=m_new;}

    void Exact_Resize(const ID m_new,const bool initialize_new_elements=true) // zero elbow room
    {Resize_Helper(m_new,initialize_new_elements);m=buffer_size;}

    ID Append(const T& element) PHYSBAM_ALWAYS_INLINE
    {Ensure_Enough_Space(m+1);m++;(*this)(m)=element;return m;}

    template<class T_ARRAY>
    void Append_Elements(const T_ARRAY& append_array)
    {STATIC_ASSERT_SAME(ELEMENT,typename T_ARRAY::ELEMENT);ID m_new=m+Value(append_array.Size());Ensure_Enough_Space(m_new);m=m_new;
    for(typename T_ARRAY::INDEX i(1);i<=append_array.Size();i++) (*this)(m-Value(append_array.Size())+Value(i))=append_array(i);}

    void Append_Unique(const T& element)
    {for(ID i(1);i<=m;i++) if((*this)(i)==element) return;Append(element);}

    template<class T_ARRAY>
    void Append_Unique_Elements(const T_ARRAY& append_array)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY::ELEMENT);
    typename T_ARRAY::INDEX append_m=append_array.Size();for(typename T_ARRAY::INDEX i(1);i<=append_m;i++) Append_Unique(append_array(i));}

    void Remove_End()
    {assert(m>ID());m--;}

    void Remove_Index(const ID index) // preserves ordering of remaining elements
    {assert(ID(1)<=index && index<=m);for(ID i=index;i<m;i++) (*this)(i)=(*this)(i+1);Remove_End();}

    void Remove_Index_Lazy(const ID index)
    {assert(ID(1)<=index && index<=m);
    if(index<m) (*this)(index)=(*this)(m);
    Remove_End();}

    void Remove_All() // if elements are non-primitive this may waste memory
    {m=ID();}

    void Clean_Memory()
    {Exact_Resize(ID());}

    void Delete_Pointers_And_Clean_Memory() // only valid if T is a pointer type
    {for(ID i(1);i<=m;i++) delete (*this)(i);Clean_Memory();}

    void Insert(const T& element,const ID index)
    {Ensure_Enough_Space(m+1);m++;for(ID i=m;i>index;i--) (*this)(i)=(*this)(i-1);(*this)(index)=element;}

    T Pop()
    {Remove_End();return base_pointer[m];}

    ARRAY_VIEW<const T> Pop_Elements(const int count) // return value should be copied immediately, not kept around
    {static const bool has_trivial_destructor=HAS_TRIVIAL_DESTRUCTOR<T>::value;
    STATIC_ASSERT(has_trivial_destructor);
    assert(m-count>=ID());m-=count;
    return ARRAY_VIEW<const T>(count,base_pointer+m);}

    template<class T2>
    static void Compact_Array_Using_Compaction_Array(ARRAY<T2,ID>& array,const ARRAY<ID,ID>& compaction_array,ARRAY<T2,ID>* temporary_array=0)
    {ID compaction_array_m=compaction_array.Size();
    bool temporary_array_defined=temporary_array!=0;if(!temporary_array_defined) temporary_array=new ARRAY<T2,ID>(compaction_array_m,false);
    ARRAY<T2,ID>::Put(array,*temporary_array);for(ID i(1);i<=compaction_array_m;i++) if(compaction_array(i)>0) array(compaction_array(i))=(*temporary_array)(i);
    if(!temporary_array_defined){delete temporary_array;temporary_array=0;}}

private:
    template<class T2>
    static void exchange(T2& a,T2& b)
    {T2 tmp=a;a=b;b=tmp;}

public:
    void Exchange(ARRAY<T,ID>& other)
    {STATIC_ASSERT(!IS_CONST<T>::value); // make ARRAY_VIEW<const T> equivalent to const ARRAY_VIEW<const T>
    exchange(m,other.m);exchange(base_pointer,other.base_pointer);exchange(buffer_size,other.buffer_size);}

//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2005-2007, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_ARRAY
//#####################################################################
#ifndef __FACE_ARRAY__
#define __FACE_ARRAY__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND_VIEW.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/SCALAR_POLICY.h>
namespace PhysBAM{

template<class T,int d>
class ARRAY<T,FACE_INDEX<d> >
{
    typedef VECTOR<T,d> TV;typedef VECTOR<int,d> TV_INT;
    struct UNUSABLE{};
    template<class S> struct ELEMENT_OF{typedef typename S::ELEMENT TYPE;};
    typedef typename IF<IS_VECTOR<T>::value,ELEMENT_OF<T>,FIRST<UNUSABLE> >::TYPE::TYPE ELEMENT_OF_T;    
public:
    T* base_pointer;
    int buffer_size;

    enum {dimension=TV::dimension};
    template<class T2> struct REBIND{typedef ARRAY<T2,FACE_INDEX<dimension> > TYPE;};
    typedef ARRAY_VIEW<T,TV_INT> T_ARRAY_VIEW;
    typedef T ELEMENT;typedef FACE_INDEX<d> INDEX;typedef T& RESULT_TYPE;typedef const T& CONST_RESULT_TYPE;
    typedef typename SCALAR_POLICY<T>::TYPE SCALAR;

    RANGE<TV_INT> domain_indices;
    VECTOR<T_ARRAY_VIEW,dimension> data;
   
    ARRAY()
        :base_pointer(0),buffer_size(0)
    {}

    template<class T2>
    ARRAY(const GRID<VECTOR<T2,dimension> >& grid,const int ghost_cells=0,const bool initialize_using_default_constructor=true)
        :base_pointer(0),buffer_size(0)
    {Resize(grid,ghost_cells,initialize_using_default_constructor);}

    ARRAY(const RANGE<TV_INT>& domain_indices_input,const bool initialize_using_default_constructor=true)
        :base_pointer(0),buffer_size(0)
    {Resize(domain_indices_input,initialize_using_default_constructor);}

    ARRAY(const ARRAY& old_array)
        :base_pointer(0),buffer_size(0)
    {Resize(old_array.Domain_Indices(),false,false);Copy(old_array,*this);}

    virtual ~ARRAY()
    {delete[] base_pointer;}

    template<class T2>
    ARRAY& operator*=(const T2 a)
    {for(int i=1;i<=dimension;i++) data(i)*=a;return *this;}

    ARRAY& operator+=(const ARRAY& a)
    {data+=a.data;return *this;}

    const RANGE<TV_INT>& Domain_Indices() const
    {return domain_indices;}

    ARRAY& operator=(const ARRAY& source)
    {if(source.Domain_Indices()!=Domain_Indices()) Resize(source.Domain_Indices(),false,false);Copy(source,*this);return *this;}

    template<class T2>
    void Resize(const GRID<VECTOR<T2,dimension> >& grid,const int ghost_cells=0,const bool initialize_new_elements=true,const bool copy_existing_elements=true,const T& initialization_value=T())
    {Resize(grid.Domain_Indices(ghost_cells),initialize_new_elements,copy_existing_elements,initialization_value);}

    virtual void Resize(const RANGE<TV_INT>& domain,const bool initialize_new_elements=true,const bool copy_existing_elements=true,const T& initialization_value=T())
    {TV_INT sizes_new;bool early_return=true;
    VECTOR<RANGE<TV_INT>,dimension> domains;
    domain_indices=domain;
    for(int i=1;i<=dimension;i++){domains(i)=domain;domains(i).max_corner(i)++;sizes_new(i)=(domains(i).Edge_Lengths()+1).Product();if(domains(i)!=data(i).Domain_Indices()) early_return=false;}
    if(early_return) return;
    buffer_size=sizes_new.Sum();
    T* p=new T[buffer_size];
    if(initialize_new_elements){T* p_start=p;
        for(int i=1;i<=dimension;i++){
            T_ARRAY_VIEW array_new(domains(i),p_start);
            array_new.Fill(initialization_value);
            p_start+=sizes_new(i);}}
    if(copy_existing_elements){T* p_start=p;
        for(int i=1;i<=dimension;i++){
            T_ARRAY_VIEW array_new(domains(i),p_start);
            T_ARRAY_VIEW::Limited_Shifted_Get(array_new,data(i),TV_INT());
            p_start+=sizes_new(i);}}
    T* p_start=p;
    for(int i=1;i<=dimension;i++){
        T_ARRAY_VIEW array_new(domains(i),p_start);
        T_ARRAY_VIEW::Exchange_Arrays(array_new,data(i));
        p_start+=sizes_new(i);}
    delete[] base_pointer;
    base_pointer=p;}

    void Clean_Memory()
    {Resize(RANGE<TV_INT>(TV_INT::All_Ones_Vector(),TV_INT()),false,false);}

    void Delete_Pointers_And_Clean_Memory()
    {for(int i=0;i<buffer_size;i++) delete base_pointer[i];Clean_Memory();}

    int Number_Of_Ghost_Cells() const
    {return 1-domain_indices.min_corner.x;}

    T& operator()(const int axis,const TV_INT& face)
    {return Component(axis)(face);}

    const T& operator()(const int axis,const TV_INT& face) const
    {return Component(axis)(face);}

    T& operator()(const FACE_INDEX<dimension>& index)
    {return Component(index.axis)(index.index);}

    const T& operator()(const FACE_INDEX<dimension>& index) const
    {return Component(index.axis)(index.index);}

    bool Valid_Index(const FACE_INDEX<dimension>& index) const
    {return (1<=index.axis && index.axis<=dimension) && Component(index.axis).Valid_Index(index.index);}

    T_ARRAY_VIEW& Component(const int axis)
    {assert(1<=axis&&axis<=dimension);return *(&data.x+(axis-1));}

    const T_ARRAY_VIEW& Component(const int axis) const
    {assert(1<=axis&&axis<=dimension);return *(&data.x+(axis-1));}

    TV Cell_Centered_Average(const TV_INT& cell_index) const
    {TV average;for(int i=1;i<=dimension;i++) average(i)=(T).5*(data(i)(cell_index)+data(i)(cell_index+TV_INT::Axis_Vector(i)));return average;}

    void Set_All_Faces(const T& value,const TV_INT& cell_index)
    {for(int i=1;i<=dimension;i++) data(i)(cell_index)=data(i)(cell_index+TV_INT::Axis_Vector(i))=value;}

    static void Extract_Dimension(const ARRAY& old_array,ARRAY<ELEMENT_OF_T,FACE_INDEX<dimension> >& extracted_array,const TV_INT& dimensions_to_extract)
    {STATIC_ASSERT(IS_VECTOR<T>::value);for(int i=1;i<=dimension;i++) T_ARRAY_VIEW::Extract_Dimension(old_array.data(i),extracted_array.data(i),dimensions_to_extract(i));}

    void Fill(const T& constant)
    {for(int i=1;i<=dimension;i++) data(i).Fill(constant);}

    static void Copy(const ARRAY& old_copy,ARRAY& new_copy)
    {for(int i=1;i<=dimension;i++) T_ARRAY_VIEW::Copy(old_copy.data(i),new_copy.data(i));}

    template<class T2>
    static void Copy(const T2 c,const ARRAY& old,ARRAY& result)
    {for(int i=1;i<=dimension;i++) T_ARRAY_VIEW::Copy(c,old.data(i),result.data(i));}

    void Fill(const TV& value)
    {for(int i=1;i<=dimension;i++) data(i).Fill(value.x);}

    template<class T2>
    static void Copy(const T2 c1,const ARRAY& v1,const T2 c2,const ARRAY& v2,ARRAY& result)
    {for(int i=1;i<=dimension;i++) T_ARRAY_VIEW::Copy(c1,v1.data(i),c2,v2.data(i),result.data(i));}

    static void Put(const ARRAY& old_copy,ARRAY& new_copy)
    {for(int i=1;i<=dimension;i++) T_ARRAY_VIEW::Put(old_copy.data(i),new_copy.data(i));}

    TV Maxabs() const
    {TV maxabs_vector;for(int i=1;i<=dimension;i++) maxabs_vector(i)=data(i).Maxabs();return maxabs_vector;}

    static void Exchange_Arrays(ARRAY& a,ARRAY& b)
    {exchange(a.domain_indices,b.domain_indices);exchange(a.base_pointer,b.base_pointer);exchange(a.buffer_size,b.buffer_size);
    for(int i=1;i<=dimension;i++) T_ARRAY_VIEW::Exchange_Arrays(a.data(i),b.data(i));}

    template<class T_INDICES>
    INDIRECT_ARRAY<ARRAY,T_INDICES&> Subset(const T_INDICES& indices)
    {return INDIRECT_ARRAY<ARRAY,T_INDICES&>(*this,indices);}

    template<class T_INDICES>
    INDIRECT_ARRAY<const ARRAY,T_INDICES&> Subset(const T_INDICES& indices) const
    {return INDIRECT_ARRAY<const ARRAY,T_INDICES&>(*this,indices);}
    
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2005-2008, Geoffrey Irving, Frank Losasso, Avi Robinson-Mosher, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_PACKAGE
//#####################################################################
#ifndef __MPI_PACKAGE__
#define __MPI_PACKAGE__

#ifdef USE_MPI

#include <PhysBAM_Tools/Arrays/CONSTANT_ARRAY.h>
#include <PhysBAM_Tools/Arrays/PROJECTED_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
namespace PhysBAM{

class MPI_PACKAGE
{
public:
    void* data;
    int count;
    MPI::Datatype type;
    bool type_needs_free;

    MPI_PACKAGE()
        :data(0),count(0),type(MPI::DATATYPE_NULL),type_needs_free(false)
    {}

    template<class T>
    MPI_PACKAGE(ARRAY_VIEW<T> array,const RANGE<VECTOR<int,1> >& indices)
        :data(&array(indices.min_corner.x)),count(indices.max_corner.x-indices.min_corner.x+1),type(MPI_UTILITIES::Datatype<T>()),type_needs_free(false)
    {
        assert(array.Valid_Index(indices.max_corner.x));
    }

    template<class T>
    MPI_PACKAGE(ARRAY<T>& array,const RANGE<VECTOR<int,1> >& indices)
        :data(&array(indices.min_corner.x)),count(indices.max_corner.x-indices.min_corner.x+1),type(MPI_UTILITIES::Datatype<T>()),type_needs_free(false)
    {
        assert(array.Valid_Index(indices.max_corner.x));
    }

    template<class T>
    MPI_PACKAGE(ARRAY_VIEW<T> array,const MPI::Datatype& type_input)
        :data(array.Get_Array_Pointer()-1),count(1),type(type_input),type_needs_free(true)
    {
        type.Commit();
    }

    template<class T>
    MPI_PACKAGE(ARRAY<T>& array,const MPI::Datatype& type_input)
        :data(array.Get_Array_Pointer()-1),count(1),type(type_input),type_needs_free(true)
    {
        type.Commit();
    }

    template<class T,class T_ARRAY,class T_BOX_INT>
    MPI_PACKAGE(ARRAY_BASE<T,T_ARRAY,typename T_ARRAY::INDEX>& array,const T_BOX_INT& indices)
        :data(&array(indices.Minimum_Corner().x)),count(1),type_needs_free(true)
    {
        type=Make_Arrays_Type(array,indices);
        type.Commit();
    }

    template<class T,int dimension,class T_BOX_INT>
    MPI_PACKAGE(ARRAYS_ND_BASE<VECTOR<T,dimension> >& array,const T_BOX_INT& indices)
        :data(&array(indices.Minimum_Corner())),count(1),type_needs_free(true)
    {
        type=Make_Arrays_Type(array,indices);
        type.Commit();
    }

    template<class T>
    MPI_PACKAGE(ARRAY_VIEW<VECTOR<T,0> > array)
    {PHYSBAM_NOT_IMPLEMENTED();}

    template<class T>
    explicit MPI_PACKAGE(INDIRECT_ARRAY<ARRAY_VIEW<T> > array)
        :data(array.array.Get_Array_Pointer()-1),count(1),type_needs_free(true)
    {
        ARRAY<int> lengths(CONSTANT_ARRAY<int>(array.indices.m,1));
        type=MPI_UTILITIES::Datatype<T>().Create_indexed(array.indices.m,lengths.Get_Array_Pointer(),array.indices.Get_Array_Pointer());
        type.Commit();
    }

    template<class T>
    explicit MPI_PACKAGE(ARRAY_VIEW<T> array)
        :data(array.Get_Array_Pointer()),count(array.Size()),type(MPI_UTILITIES::Datatype<T>()),type_needs_free(false)
    {}

    template<class T>
    explicit MPI_PACKAGE(ARRAY<T>& array)
        :data(array.Get_Array_Pointer()),count(array.Size()),type(MPI_UTILITIES::Datatype<T>()),type_needs_free(false)
    {}

    explicit MPI_PACKAGE(const MPI::Datatype& type_input)
        :data((void*)MPI::BOTTOM),count(1),type(type_input),type_needs_free(true)
    {
        type.Commit();
    }

    ~MPI_PACKAGE()
    {}

    void Free()
    {if(type_needs_free) type.Free();type_needs_free=false;}

    int Size() const
    {return count*type.Get_size();}

    int Pack_Size(const MPI::Comm& comm) const
    {return type.Pack_size(count,comm);}

    void Pack(ARRAY<char>& buffer,const MPI::Comm& comm) const
    {int position=0;type.Pack(data,count,buffer.Get_Array_Pointer(),buffer.m,position,comm);}

    void Unpack(ARRAY<char>& buffer,const MPI::Comm& comm)
    {int position=0;type.Unpack(buffer.Get_Array_Pointer(),buffer.m,data,count,position,comm);}

    static void Free_All(ARRAY<MPI_PACKAGE>& packages)
    {for(int p=1;p<=packages.m;p++)packages(p).Free();packages.Clean_Memory();}

    MPI::Request Isend(const MPI::Comm& comm,const int rank,const int tag) const
    {return comm.Isend(data,count,type,rank,tag);}

    MPI::Request Irecv(const MPI::Comm& comm,const int rank,const int tag)
    {return comm.Irecv(data,count,type,rank,tag);}

    void Recv(const MPI::Comm& comm,const int rank,const int tag)
    {comm.Recv(data,count,type,rank,tag);}

    void Broadcast(const MPI::Intracomm& comm)
    {comm.Bcast(data,count,type,0);}

    template<class T2> static MPI::Datatype Make_Arrays_Type(const ARRAY_BASE<T2,ARRAY<T2> >& array,const RANGE<VECTOR<int,1> >& box)
    {assert(array.Valid_Index(box.Minimum_Corner().x) && array.Valid_Index(box.Maximum_Corner().x));
    return Make_Arrays_Type_1D(MPI_UTILITIES::Datatype<T2>(),1,box);}

    template<class T2> static MPI::Datatype Make_Arrays_Type(const ARRAYS_ND_BASE<VECTOR<T2,1> >& array,const RANGE<VECTOR<int,1> >& box)
    {assert(array.Valid_Index(box.Minimum_Corner()) && array.Valid_Index(box.Maximum_Corner()));
    return Make_Arrays_Type_1D(MPI_UTILITIES::Datatype<T2>(),1,box);}

    template<class T2> static MPI::Datatype Make_Arrays_Type(const ARRAYS_ND_BASE<VECTOR<T2,2> >& array,const RANGE<VECTOR<int,2> >& box)
    {assert(array.Valid_Index(box.Minimum_Corner()) && array.Valid_Index(box.Maximum_Corner()));
    return Make_Arrays_Type_2D(MPI_UTILITIES::Datatype<T2>(),1,array.counts.y,box);}

    template<class T2> static MPI::Datatype Make_Arrays_Type(const ARRAYS_ND_BASE<VECTOR<T2,3> >& array,const RANGE<VECTOR<int,3> >& box)
    {assert(array.Valid_Index(box.Minimum_Corner()) && array.Valid_Index(box.Maximum_Corner()));
    return Make_Arrays_Type_3D(MPI_UTILITIES::Datatype<T2>(),1,array.counts.y,array.counts.z,box);}

    template<class T_ARRAY>
    static MPI_PACKAGE Union(const T_ARRAY& packages)
    {STATIC_ASSERT((IS_SAME<MPI_PACKAGE,typename T_ARRAY::ELEMENT>::value));
    ARRAY<MPI::Datatype> datatypes(packages.template Project<MPI::Datatype,&MPI_PACKAGE::type>());
    ARRAY<int> lengths(CONSTANT_ARRAY<int>(packages.Size(),1));ARRAY<void*> pointers(packages.template Project<void*,&MPI_PACKAGE::data>());
    return MPI_PACKAGE(MPI::Datatype::Create_struct(packages.Size(),lengths.Get_Array_Pointer(),(MPI::Aint*)pointers.Get_Array_Pointer(),datatypes.Get_Array_Pointer()));}

    MPI_PACKAGE Subset(const ARRAY<int>& indices) const
    {MPI_PACKAGE package;package.data=(unsigned char*)data-type.Get_size();package.count=1;package.type_needs_free=true;
    package.type=type.Create_indexed(indices.m,ARRAY<int>(CONSTANT_ARRAY<int>(indices.m,1)).Get_Array_Pointer(),indices.Get_Array_Pointer());
    package.type.Commit();return package;}

private:
    static MPI::Datatype Make_Arrays_Type_1D(const MPI::Datatype base_datatype,const int array_length,const RANGE<VECTOR<int,1> >& box)
    {return base_datatype.Create_vector(box.max_corner.x-box.min_corner.x+1,array_length,array_length);}

    static MPI::Datatype Make_Arrays_Type_2D(const MPI::Datatype base_datatype,const int array_length,const int n,const RANGE<VECTOR<int,2> >& box)
    {MPI::Datatype type_1d=Make_Arrays_Type_1D(base_datatype,array_length,RANGE<VECTOR<int,1> >(box.min_corner.y,box.max_corner.y));
    MPI::Datatype type_2d=type_1d.Create_hvector(box.max_corner.x-box.min_corner.x+1,1,base_datatype.Get_size()*array_length*n);
    type_1d.Free();return type_2d;}

    static MPI::Datatype Make_Arrays_Type_3D(const MPI::Datatype base_datatype,const int array_length,const int n,const int mn,const RANGE<VECTOR<int,3> >& box)
    {MPI::Datatype type_2d=Make_Arrays_Type_2D(base_datatype,array_length,mn,RANGE<VECTOR<int,2> >(box.min_corner.y,box.max_corner.y,box.min_corner.z,box.max_corner.z));
    MPI::Datatype type_3d=type_2d.Create_hvector(box.max_corner.x-box.min_corner.x+1,1,base_datatype.Get_size()*array_length*mn*n);
    type_2d.Free();return type_3d;}

//#####################################################################
};
}
#endif
#endif

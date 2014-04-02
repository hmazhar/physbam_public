//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace MPI_UTILITIES
//#####################################################################
#ifdef USE_MPI
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE/RLE_RUN_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_RUN_3D.h>
#endif
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
namespace PhysBAM{
namespace MPI_UTILITIES{
//#####################################################################
// Function RLE_Run_Datatype
//#####################################################################
template<class T_RUN> MPI::Datatype
RLE_Run_Datatype()
{
    static MPI::Datatype datatype=MPI::DATATYPE_NULL;
    if(datatype==MPI::DATATYPE_NULL){
        static T_RUN zero_run;static int lengths[3]={1,1,1};
        static MPI::Aint displacements[3]={0,(char*)&zero_run.jmin-(char*)&zero_run,sizeof(T_RUN)};
        static MPI::Datatype old_types[3]={MPI::BYTE,MPI::SHORT,MPI::UB};
        datatype=MPI::Datatype::Create_struct(3,lengths,displacements,old_types);
        datatype.Commit();}
    return datatype;
}
//#####################################################################
// Datatype conversion for SPARSE_MATRIX_ENTRY
//#####################################################################
template<class T> MPI::Datatype
DATATYPE_HELPER<SPARSE_MATRIX_ENTRY<T> >::Datatype()
{
    static MPI::Datatype datatype=MPI::DATATYPE_NULL;
    if(datatype==MPI::DATATYPE_NULL){
        static SPARSE_MATRIX_ENTRY<T> zero_entry;static int lengths[2]={1,1};
        static MPI::Aint displacements[2]={0,(char*)&zero_entry.a-(char*)&zero_entry};
        static MPI::Datatype old_types[2]={MPI::INT,MPI_UTILITIES::Datatype<T>()};
        datatype=MPI::Datatype::Create_struct(2,lengths,displacements,old_types);
        datatype.Commit();}
    return datatype;
}
//#####################################################################
// Datatype conversion for PAIR<VECTOR<int,2>,VECTOR<T,2> >
//#####################################################################
template<class T>  MPI::Datatype
DATATYPE_HELPER<PAIR<VECTOR<int,2>,VECTOR<T,2> > >::Datatype()
{
    static MPI::Datatype datatype=MPI::DATATYPE_NULL;
    if(datatype==MPI::DATATYPE_NULL){
        static PAIR<VECTOR<int,2>,VECTOR<T,2> > zero_entry;static int lengths[3]={2,2,1};
        static MPI::Aint displacements[3]={0,(char*)&zero_entry.y-(char*)&zero_entry,sizeof(zero_entry)};
        static MPI::Datatype old_types[3]={MPI::INT,MPI_UTILITIES::Datatype<T>(),MPI::UB};
        datatype=MPI::Datatype::Create_struct(3,lengths,displacements,old_types);
        datatype.Commit();}
    return datatype;
}
//#####################################################################
// Function Scalar_Block_Datatype
//#####################################################################
template<class T,int d> MPI::Datatype Scalar_Block_Datatype()
{
    STATIC_ASSERT((AND<IS_SCALAR<T>::value,(d>1)>::value));
    static MPI::Datatype datatype=MPI::DATATYPE_NULL;
    if(datatype==MPI::DATATYPE_NULL){datatype=MPI_UTILITIES::Datatype<T>().Create_contiguous(d);datatype.Commit();}
    return datatype;
}
//#####################################################################
// Pack/Unpack for VECTOR_ND
//#####################################################################
template<class T> int Pack_Size(const VECTOR_ND<T>& data,const MPI::Comm& comm)
{return Pack_Size<int>(comm)+Datatype<T>().Pack_size(data.n,comm);}

template<class T> void Pack(const VECTOR_ND<T>& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm)
{assert(Pack_Size(data,comm)<=buffer.Size()-position);
Pack(data.n,buffer,position,comm);
Datatype<T>().Pack(data.Get_Array_Pointer(),data.n,&buffer(1),buffer.Size(),position,comm);}

template<class T> void Unpack(VECTOR_ND<T>& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm)
{int n;Unpack(n,buffer,position,comm);data.Resize(n);
Datatype<T>().Unpack(&buffer(1),buffer.Size(),data.Get_Array_Pointer(),data.n,position,comm);}
//#####################################################################
// Pack/Unpack for SPARSE_MATRIX_FLAT_NXN
//#####################################################################
template<class T> int Pack_Size(const SPARSE_MATRIX_FLAT_NXN<T>& data,const MPI::Comm& comm)
{return Pack_Size(data.n,data.offsets,data.A,comm);}

template<class T> void Pack(const SPARSE_MATRIX_FLAT_NXN<T>& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm)
{Pack(data.n,data.offsets,data.A,buffer,position,comm);}

template<class T> void Unpack(SPARSE_MATRIX_FLAT_NXN<T>& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm)
{Unpack(data.n,data.offsets,data.A,buffer,position,comm);}
//#####################################################################
// Pack/Unpack for SPARSE_MATRIX_FLAT_MXN
//#####################################################################
template<class T> int Pack_Size(const SPARSE_MATRIX_FLAT_MXN<T>& data,const MPI::Comm& comm)
{return Pack_Size(data.m,data.n,data.offsets,data.A,comm);}

template<class T> void Pack(const SPARSE_MATRIX_FLAT_MXN<T>& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm)
{Pack(data.m,data.n,data.offsets,data.A,buffer,position,comm);}

template<class T> void Unpack(SPARSE_MATRIX_FLAT_MXN<T>& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm)
{Unpack(data.m,data.n,data.offsets,data.A,buffer,position,comm);}
//#####################################################################
// Pack/Unpack for SPARSE_MATRIX_PARTITION
//#####################################################################
int Pack_Size(const SPARSE_MATRIX_PARTITION& data,const MPI::Comm& comm)
{int size=Pack_Size(data.number_of_sides,data.interior_indices,data.ghost_indices,data.neighbor_ranks,comm);
for(int s=1;s<=data.number_of_sides;s++)size+=Pack_Size(data.boundary_indices(s),comm);
return size;}

void Pack(const SPARSE_MATRIX_PARTITION& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm)
{Pack(data.number_of_sides,data.interior_indices,data.ghost_indices,data.neighbor_ranks,buffer,position,comm);
for(int s=1;s<=data.number_of_sides;s++)Pack(data.boundary_indices(s),buffer,position,comm);}

void Unpack(SPARSE_MATRIX_PARTITION& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm)
{Unpack(data.number_of_sides,data.interior_indices,data.ghost_indices,data.neighbor_ranks,buffer,position,comm);data.neighbors.Resize(data.number_of_sides);
data.boundary_indices.Resize(data.number_of_sides);for(int s=1;s<=data.number_of_sides;s++)Unpack(data.boundary_indices(s),buffer,position,comm);}
//#####################################################################
// Pack/Unpack for ARRAY<ARRAY<int> >
//#####################################################################
int Pack_Size(const ARRAY<ARRAY<int> >& data,const MPI::Comm& comm)
{int size=Pack_Size<int>(comm);for(int i=1;i<=data.m;i++) size+=Pack_Size(data(i),comm);return size;}

void Pack(const ARRAY<ARRAY<int> >& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm)
{assert(Pack_Size(data,comm)<=buffer.Size()-position);
Pack(data.m,buffer,position,comm);
for(int i=1;i<=data.m;i++) Pack(data(i),buffer,position,comm);}

void Unpack(ARRAY<ARRAY<int> >& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm)
{int m;Unpack(m,buffer,position,comm);data.Resize(m);
for(int i=1;i<=data.m;i++) Unpack(data(i),buffer,position,comm);}
//#####################################################################
// Pack/Unpack for UNION_FIND
//#####################################################################
int Pack_Size(const UNION_FIND<>& data,const MPI::Comm& comm)
{return Pack_Size(data.parents,comm)+Pack_Size(data.ranks,comm);}

void Pack(const UNION_FIND<>& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm)
{Pack(data.parents,data.ranks,buffer,position,comm);}

void Unpack(UNION_FIND<>& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm)
{Unpack(data.parents,data.ranks,buffer,position,comm);}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template class DATATYPE_HELPER<PAIR<VECTOR<int,2>,VECTOR<T,2> > >; \
    template MPI::Datatype Scalar_Block_Datatype<T,2>(); \
    template MPI::Datatype Scalar_Block_Datatype<T,3>(); \
    template MPI::Datatype Scalar_Block_Datatype<T,4>(); \
    template MPI::Datatype Scalar_Block_Datatype<T,5>(); \
    template MPI::Datatype Scalar_Block_Datatype<T,6>(); \
    template MPI::Datatype Scalar_Block_Datatype<T,7>(); \
    template MPI::Datatype Scalar_Block_Datatype<T,8>(); \
    template MPI::Datatype Scalar_Block_Datatype<T,9>(); \
    template MPI::Datatype Scalar_Block_Datatype<T,10>(); \
    template int Pack_Size<T>(const VECTOR_ND<T>&,const MPI::Comm&); \
    template void Pack<T>(const VECTOR_ND<T>&,ARRAY_VIEW<char>,int&,const MPI::Comm&); \
    template void Unpack<T>(VECTOR_ND<T>&,ARRAY_VIEW<const char>,int&,const MPI::Comm&); \
    template int Pack_Size<T>(const SPARSE_MATRIX_FLAT_NXN<T>&,const MPI::Comm&); \
    template void Pack<T>(const SPARSE_MATRIX_FLAT_NXN<T>&,ARRAY_VIEW<char>,int&,const MPI::Comm&); \
    template void Unpack<T>(SPARSE_MATRIX_FLAT_NXN<T>&,ARRAY_VIEW<const char>,int&,const MPI::Comm&); \
    template int Pack_Size<T>(const SPARSE_MATRIX_FLAT_MXN<T>&,const MPI::Comm&); \
    template void Pack<T>(const SPARSE_MATRIX_FLAT_MXN<T>&,ARRAY_VIEW<char>,int&,const MPI::Comm&); \
    template void Unpack<T>(SPARSE_MATRIX_FLAT_MXN<T>&,ARRAY_VIEW<const char>,int&,const MPI::Comm&);
template MPI::Datatype Scalar_Block_Datatype<int,2>();
template MPI::Datatype Scalar_Block_Datatype<int,3>();
template MPI::Datatype Scalar_Block_Datatype<int,4>();
template MPI::Datatype Scalar_Block_Datatype<int,5>();
template MPI::Datatype Scalar_Block_Datatype<int,6>();
template MPI::Datatype Scalar_Block_Datatype<int,7>();
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template MPI::Datatype RLE_Run_Datatype<RLE_RUN>();
template MPI::Datatype RLE_Run_Datatype<RLE_RUN_2D>();
template MPI::Datatype RLE_Run_Datatype<RLE_RUN_3D>();
#endif
INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
}
}
#endif

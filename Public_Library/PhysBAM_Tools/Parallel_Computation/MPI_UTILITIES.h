//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace MPI_UTILITIES
//#####################################################################
#ifndef __PhysBAM_Tools_MPI_UTILITIES__
#define __PhysBAM_Tools_MPI_UTILITIES__

#ifdef USE_MPI

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <mpi.h>
namespace PhysBAM{

class SPARSE_MATRIX_PARTITION;
template<class T> class SPARSE_MATRIX_ENTRY;
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
class RLE_RUN;
class RLE_RUN_2D;
class RLE_RUN_3D;
#endif

namespace MPI_UTILITIES{

//#####################################################################
// Datatype conversion
//#####################################################################
template<class T,class VOID=void> struct DATATYPE_HELPER{typedef int NO_DATATYPE_EXISTS;};
template<class T> struct DATATYPE_HELPER<const T>:public DATATYPE_HELPER<T>{};

template<> struct DATATYPE_HELPER<int>{static MPI::Datatype Datatype(){return MPI::INT;}};
template<> struct DATATYPE_HELPER<bool>{static MPI::Datatype Datatype(){return MPI::BOOL;}};
template<> struct DATATYPE_HELPER<unsigned char>{static MPI::Datatype Datatype(){return MPI::UNSIGNED_CHAR;}};
template<> struct DATATYPE_HELPER<unsigned short>{static MPI::Datatype Datatype(){return MPI::UNSIGNED_SHORT;}};
template<> struct DATATYPE_HELPER<float>{static MPI::Datatype Datatype(){return MPI::FLOAT;}};
template<> struct DATATYPE_HELPER<double>{static MPI::Datatype Datatype(){return MPI::DOUBLE;}};
template<class T> struct DATATYPE_HELPER<VECTOR<T,1> >:public DATATYPE_HELPER<T>{};
template<class T> struct DATATYPE_HELPER<MATRIX<T,1,1> >:public DATATYPE_HELPER<T>{};

template<class T,int d> MPI::Datatype Scalar_Block_Datatype();
template<class TV> struct DATATYPE_HELPER<TV,typename ENABLE_IF<AND<IS_SCALAR_BLOCK<TV>::value,(sizeof(TV)>sizeof(typename TV::SCALAR))>::value>::TYPE>{static MPI::Datatype Datatype()
{typedef typename TV::SCALAR T;return Scalar_Block_Datatype<T,sizeof(TV)/sizeof(T)>();}};

#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template<class T_RUN> MPI::Datatype RLE_Run_Datatype();
template<> struct DATATYPE_HELPER<RLE_RUN>{static MPI::Datatype Datatype(){return RLE_Run_Datatype<RLE_RUN>();}};
template<> struct DATATYPE_HELPER<RLE_RUN_2D>{static MPI::Datatype Datatype(){return RLE_Run_Datatype<RLE_RUN_2D>();}};
template<> struct DATATYPE_HELPER<RLE_RUN_3D>{static MPI::Datatype Datatype(){return RLE_Run_Datatype<RLE_RUN_3D>();}};
#endif

template<class T> struct DATATYPE_HELPER<SPARSE_MATRIX_ENTRY<T> >{static MPI::Datatype Datatype();};
template<class T> struct DATATYPE_HELPER<PAIR<VECTOR<int,2>,VECTOR<T,2> > >{static MPI::Datatype Datatype();};

template<class T> inline MPI::Datatype Datatype()
{return DATATYPE_HELPER<T>::Datatype();}

template<class T,class HAS=void> struct HAS_DATATYPE{enum {value=true};};
template<class T> struct HAS_DATATYPE<T,typename FIRST<void,typename DATATYPE_HELPER<T>::NO_DATATYPE_EXISTS>::TYPE>{enum {value=false};};
//#####################################################################
// Pack/Unpack for classes with datatypes
//#####################################################################
template<class T> inline int Pack_Size(const MPI::Comm& comm)
{return Datatype<T>().Pack_size(1,comm);}

template<class T> inline int Pack_Size(const T& data,typename ENABLE_IF<HAS_DATATYPE<T>::value,const MPI::Comm&>::TYPE comm)
{return Pack_Size<T>(comm);}

template<class T> inline void Pack(const T& data,ARRAY_VIEW<char> buffer,int& position,typename ENABLE_IF<HAS_DATATYPE<T>::value,const MPI::Comm&>::TYPE comm)
{assert(Pack_Size(data,comm)<=buffer.Size()-position);
Datatype<T>().Pack(&data,1,buffer.Get_Array_Pointer(),buffer.Size(),position,comm);}

template<class T> inline void Unpack(T& data,ARRAY_VIEW<const char> buffer,int& position,const typename ENABLE_IF<HAS_DATATYPE<T>::value,const MPI::Comm&>::TYPE comm)
{Datatype<T>().Unpack(buffer.Get_Array_Pointer(),buffer.Size(),&data,1,position,comm);}
//#####################################################################
// Pack/Unpack for arrays
//#####################################################################
template<class T> inline int Pack_Size(const ARRAY_VIEW<T>& data,const MPI::Comm& comm)
{return Pack_Size<int>(comm)+Datatype<T>().Pack_size(data.Size(),comm);}

template<class T_ARRAY,class T_INDICES> int Pack_Size(const INDIRECT_ARRAY<T_ARRAY,T_INDICES>& data,const MPI::Comm& comm)
{return Pack_Size<int>(comm)+Datatype<typename T_ARRAY::ELEMENT>().Pack_size(data.Size(),comm);}

template<class T> inline int Pack_Size(const ARRAY<T>& data,const MPI::Comm& comm)
{return Pack_Size(ARRAY_VIEW<const T>(data),comm);}

template<class T,class ID> inline int Pack_Size(const ARRAY<T,ID>& data,const MPI::Comm& comm)
{return Pack_Size(ARRAY_VIEW<const T>(data),comm);}

template<class T> inline void Pack(const ARRAY_VIEW<T>& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm)
{assert(Pack_Size(data,comm)<=buffer.Size()-position);
Pack(data.Size(),buffer,position,comm);
Datatype<T>().Pack(data.Get_Array_Pointer(),data.Size(),&buffer(1),buffer.Size(),position,comm);}

template<class T_ARRAY,class T_INDICES> inline void Pack(const INDIRECT_ARRAY<T_ARRAY,T_INDICES>& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm)
{assert(Pack_Size(data,comm)<=buffer.Size()-position);
Pack(data.Size(),buffer,position,comm);MPI::Datatype type=Datatype<typename T_ARRAY::ELEMENT>();
for(int i=1;i<=data.Size();i++) type.Pack(&data(i),1,&buffer(1),buffer.Size(),position,comm);}

template<class T> inline void Pack(const ARRAY<T>& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm)
{Pack(ARRAY_VIEW<const T>(data),buffer,position,comm);}

template<class T,class ID> inline void Pack(const ARRAY<T,ID>& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm)
{Pack(ARRAY_VIEW<const T>(data),buffer,position,comm);}

template<class T> inline void Unpack(ARRAY_VIEW<T> data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm)
{int m;Unpack(m,buffer,position,comm);PHYSBAM_ASSERT(m==data.Size());
Datatype<T>().Unpack(&buffer(1),buffer.Size(),data.Get_Array_Pointer(),data.Size(),position,comm);}

template<class T> inline void Unpack(ARRAY<T>& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm)
{int m;Unpack(m,buffer,position,comm);data.Resize(m);
Datatype<T>().Unpack(&buffer(1),buffer.Size(),data.Get_Array_Pointer(),data.Size(),position,comm);}

template<class T,class ID> inline void Unpack(ARRAY<T,ID>& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm)
{int m;Unpack(m,buffer,position,comm);data.Resize(ID(m));
Datatype<T>().Unpack(&buffer(1),buffer.Size(),data.Get_Array_Pointer(),Value(data.Size()),position,comm);}

template<class T_ARRAY,class T_INDICES> inline void Unpack(INDIRECT_ARRAY<T_ARRAY,T_INDICES>& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm)
{int m;Unpack(m,buffer,position,comm);PHYSBAM_ASSERT(m==data.Size());
MPI::Datatype type=Datatype<typename T_ARRAY::ELEMENT>();
for(int i=1;i<=data.Size();i++) type.Unpack(&buffer(1),buffer.Size(),&data(i),1,position,comm);}
//#####################################################################
// Pack/Unpack for particles
//#####################################################################
template<class T,int d> inline int Pack_Size(const POINT_CLOUD<VECTOR<T,d> >& data,const MPI::Comm& comm)
{int size=data.array_collection->Pack_Size();
PHYSBAM_ASSERT(size==MPI::UNSIGNED_CHAR.Pack_size(size,comm)); // assert that we can implement pack ourselves for particles
return size;}

template<class T_POINT_CLOUD> inline typename ENABLE_IF<IS_BASE_OF<POINT_CLOUD<typename T_POINT_CLOUD::VECTOR>,T_POINT_CLOUD>::value,int>::TYPE // work around compiler bug with enable_if
Pack_Size(const T_POINT_CLOUD& particles,const MPI::Comm& comm)
{return Pack_Size(particles,comm);}

template<class T_POINT_CLOUD> void Pack(const T_POINT_CLOUD& particles,int index,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm)
{particles.array_collection->Pack(buffer,position,index);} // valid as long as the assertion in Pack_Size succeeds

template<class T_POINT_CLOUD> void Unpack(T_POINT_CLOUD& particles,int index,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm)
{particles.array_collection->Unpack(buffer,position,index);} // valid as long as the assertion in Pack_Size succeeds
//#####################################################################
// Function Pack_Size
//#####################################################################
template<class T> int Pack_Size(const VECTOR_ND<T>& data,const MPI::Comm& comm);
template<class T> int Pack_Size(const SPARSE_MATRIX_FLAT_NXN<T>& data,const MPI::Comm& comm);
template<class T> int Pack_Size(const SPARSE_MATRIX_FLAT_MXN<T>& data,const MPI::Comm& comm);
int Pack_Size(const SPARSE_MATRIX_PARTITION& data,const MPI::Comm& comm);
int Pack_Size(const ARRAY<ARRAY<int> >& data,const MPI::Comm& comm);
template<class ID> int Pack_Size(const ARRAY<ARRAY<int>,ID>& data,const MPI::Comm& comm)
{return Pack_Size(reinterpret_cast<const ARRAY<ARRAY<int> >&>(data),comm);}
int Pack_Size(const UNION_FIND<>& data,const MPI::Comm& comm);

template<class T1,class T2> int Pack_Size(const T1& d1,const T2& d2,const MPI::Comm& comm)
{return Pack_Size(d1,comm)+Pack_Size(d2,comm);}
template<class T1,class T2,class T3> int Pack_Size(const T1& d1,const T2& d2,const T3& d3,const MPI::Comm& comm)
{return Pack_Size(d1,comm)+Pack_Size(d2,comm)+Pack_Size(d3,comm);}
template<class T1,class T2,class T3,class T4> int Pack_Size(const T1& d1,const T2& d2,const T3& d3,const T4& d4,const MPI::Comm& comm)
{return Pack_Size(d1,comm)+Pack_Size(d2,comm)+Pack_Size(d3,comm)+Pack_Size(d4,comm);}
template<class T1,class T2,class T3,class T4,class T5> int Pack_Size(const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,const MPI::Comm& comm)
{return Pack_Size(d1,comm)+Pack_Size(d2,comm)+Pack_Size(d3,comm)+Pack_Size(d4,comm)+Pack_Size(d5,comm);}
template<class T1,class T2,class T3,class T4,class T5,class T6> int Pack_Size(const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,const T6& d6,const MPI::Comm& comm)
{return Pack_Size(d1,comm)+Pack_Size(d2,comm)+Pack_Size(d3,comm)+Pack_Size(d4,comm)+Pack_Size(d5,comm)+Pack_Size(d6,comm);}
//#####################################################################
// Function Pack
//#####################################################################
template<class T> void Pack(const VECTOR_ND<T>& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm);
template<class T> void Pack(const SPARSE_MATRIX_FLAT_NXN<T>& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm);
template<class T> void Pack(const SPARSE_MATRIX_FLAT_MXN<T>& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm);
void Pack(const SPARSE_MATRIX_PARTITION& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm);
void Pack(const ARRAY<ARRAY<int> >& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm); // TODO: generalize ARRAY<T> to handle this
template<class ID> void Pack(const ARRAY<ARRAY<int>,ID>& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm) // TODO: generalize ARRAY<T> to handle this
{Pack(reinterpret_cast<const ARRAY<ARRAY<int> >&>(data),buffer,position,comm);}
void Pack(const UNION_FIND<>& data,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm);


template<class T1,class T2> inline void Pack(const T1& d1,const T2& d2,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm)
{Pack(d1,buffer,position,comm);Pack(d2,buffer,position,comm);}
template<class T1,class T2,class T3> inline void Pack(const T1& d1,const T2& d2,const T3& d3,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm)
{Pack(d1,buffer,position,comm);Pack(d2,buffer,position,comm);Pack(d3,buffer,position,comm);}
template<class T1,class T2,class T3,class T4> inline void Pack(const T1& d1,const T2& d2,const T3& d3,const T4& d4,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm)
{Pack(d1,buffer,position,comm);Pack(d2,buffer,position,comm);Pack(d3,buffer,position,comm);Pack(d4,buffer,position,comm);}
template<class T1,class T2,class T3,class T4,class T5> inline void Pack(const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm)
{Pack(d1,buffer,position,comm);Pack(d2,buffer,position,comm);Pack(d3,buffer,position,comm);Pack(d4,buffer,position,comm);Pack(d5,buffer,position,comm);}
template<class T1,class T2,class T3,class T4,class T5,class T6> inline void Pack(const T1& d1,const T2& d2,const T3& d3,const T4& d4,const T5& d5,const T6& d6,ARRAY_VIEW<char> buffer,int& position,const MPI::Comm& comm)
{Pack(d1,buffer,position,comm);Pack(d2,buffer,position,comm);Pack(d3,buffer,position,comm);Pack(d4,buffer,position,comm);Pack(d5,buffer,position,comm);Pack(d6,buffer,position,comm);}
//#####################################################################
// Function Unpack
//#####################################################################
template<class T> void Unpack(VECTOR_ND<T>& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm);
template<class T> void Unpack(SPARSE_MATRIX_FLAT_NXN<T>& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm);
template<class T> void Unpack(SPARSE_MATRIX_FLAT_MXN<T>& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm);
void Unpack(SPARSE_MATRIX_PARTITION& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm);
void Unpack(ARRAY<ARRAY<int> >& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm); // TODO: generalize ARRAY<T> to handle this
template<class ID> void Unpack(ARRAY<ARRAY<int>,ID>& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm) // TODO: generalize ARRAY<T> to handle this
{Unpack(reinterpret_cast<ARRAY<ARRAY<int> >&>(data),buffer,position,comm);}
void Unpack(UNION_FIND<>& data,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm);

template<class T1,class T2> inline void Unpack(T1& d1,T2& d2,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm)
{Unpack(d1,buffer,position,comm);Unpack(d2,buffer,position,comm);}
template<class T1,class T2,class T3> inline void Unpack(T1& d1,T2& d2,T3& d3,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm)
{Unpack(d1,buffer,position,comm);Unpack(d2,buffer,position,comm);Unpack(d3,buffer,position,comm);}
template<class T1,class T2,class T3,class T4> inline void Unpack(T1& d1,T2& d2,T3& d3,T4& d4,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm)
{Unpack(d1,buffer,position,comm);Unpack(d2,buffer,position,comm);Unpack(d3,buffer,position,comm);Unpack(d4,buffer,position,comm);}
template<class T1,class T2,class T3,class T4,class T5> inline void Unpack(T1& d1,T2& d2,T3& d3,T4& d4,T5& d5,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm)
{Unpack(d1,buffer,position,comm);Unpack(d2,buffer,position,comm);Unpack(d3,buffer,position,comm);Unpack(d4,buffer,position,comm);Unpack(d5,buffer,position,comm);}
template<class T1,class T2,class T3,class T4,class T5,class T6> inline void Unpack(T1& d1,T2& d2,T3& d3,T4& d4,T5& d5,T6& d6,ARRAY_VIEW<const char> buffer,int& position,const MPI::Comm& comm)
{Unpack(d1,buffer,position,comm);Unpack(d2,buffer,position,comm);Unpack(d3,buffer,position,comm);Unpack(d4,buffer,position,comm);Unpack(d5,buffer,position,comm);Unpack(d6,buffer,position,comm);}
//#####################################################################
// Function Wait_All
//#####################################################################
inline void Wait_All(ARRAY<MPI::Request>& requests)
{if(requests.m) MPI::Request::Waitall(requests.m,&requests(1));requests.Clean_Memory();}
//#####################################################################
// Function Wait_Any
//#####################################################################
inline bool Wait_Any(ARRAY<MPI::Request>& requests,MPI::Status& status)
{if(!requests.m) return false;int index=MPI::Request::Waitany(requests.m,&requests(1),status);requests.Remove_Index_Lazy(index+1);return true;}
//#####################################################################
// Function Free_Elements_And_Clean_Memory
//#####################################################################
template<class T> inline void Free_Elements_And_Clean_Memory(ARRAY<T>& array)
{for(int i=1;i<=array.m;i++) if(array(i)!=T()) array(i).Free();array.Clean_Memory();}
//#####################################################################
// Function Free_Elements_And_Clean_Memory
//#####################################################################
template<class T,class ID> inline void Free_Elements_And_Clean_Memory(ARRAY<T,ID>& array)
{for(ID i(1);i<=array.Size();i++) if(array(i)!=T()) array(i).Free();array.Clean_Memory();}
//#####################################################################
// Function Reduce
//#####################################################################
template<class T> inline void Reduce(const T& input,T& output,const MPI::Op& op,const MPI::Intracomm& comm)
{MPI::Datatype type=MPI_UTILITIES::Datatype<T>();assert(type==MPI::INT || type==MPI::FLOAT || type==MPI::DOUBLE);
comm.Allreduce(&input,&output,1,type,op);}
template<class T> inline void Reduce(const ARRAY<T>& input,ARRAY<T>& output,const MPI::Op& op,const MPI::Intracomm& comm)
{MPI::Datatype type=MPI_UTILITIES::Datatype<T>();assert(type==MPI::INT || type==MPI::FLOAT || type==MPI::DOUBLE);assert(input.m==output.m);
comm.Allreduce(input.Get_Array_Pointer(),output.Get_Array_Pointer(),output.m,type,op);}
//#####################################################################

}
}
#endif
#endif

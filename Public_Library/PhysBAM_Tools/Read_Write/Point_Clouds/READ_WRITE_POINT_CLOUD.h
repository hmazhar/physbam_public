//#####################################################################
// Copyright 2009, Michael Lentine, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_POINT_CLOUD
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_POINT_CLOUD__
#define __READ_WRITE_POINT_CLOUD__

#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY_COLLECTION.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<POINT_CLOUD<TV>,RW>
{
public:
    static void Read(std::istream& input,POINT_CLOUD<TV>& object)
    {int version;
    Read_Binary<RW>(input,version);
    if(version!=1) throw READ_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized particle version %d",(int)version));
    Read_Binary<RW>(input,*object.array_collection);}

    static void Write(std::ostream& output,const POINT_CLOUD<TV>& object)
    {Write_Binary<RW>(output,1,*object.array_collection);}

    static void Print(std::ostream& output,const POINT_CLOUD<TV>& object,const int p)
    {Read_Write<ARRAY_COLLECTION,RW>::Print(output,*object.array_collection,p);}
};

template<class RW,class T_POINT_CLOUD>
class Read_Write<T_POINT_CLOUD,RW,typename ENABLE_IF<AND<IS_BASE_OF<POINT_CLOUD<typename T_POINT_CLOUD::VECTOR_T> ,T_POINT_CLOUD>::value,NOT<IS_SAME<T_POINT_CLOUD,POINT_CLOUD<typename T_POINT_CLOUD::VECTOR_T> >::value>::value>::value>::TYPE>:public Read_Write<POINT_CLOUD<typename T_POINT_CLOUD::VECTOR_T>,RW>
{
};
}
#endif
#endif

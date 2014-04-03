//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_ELEMENT_ID
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_ELEMENT_ID__
#define __READ_WRITE_ELEMENT_ID__

#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class ID,class T,int flags>
class Read_Write<ELEMENT_ID<ID,T,flags>,RW>
{
public:
    static void Read(std::istream& input,ELEMENT_ID<ID,T,flags>& object)
    {T value;Read_Binary<RW>(input,value);object=ELEMENT_ID<ID,T,flags>(value);}

    static void Write(std::ostream& output,const ELEMENT_ID<ID,T,flags>& object)
    {Write_Binary<RW>(output,object.Value());}
};

template<class RW,class ID>
class Read_Write<ID,RW,typename ENABLE_IF<IS_BASE_OF<ELEMENT_ID<ID,typename ID::VALUE,ID::capability_flags>,ID>::value>::TYPE>:
    public Read_Write<ELEMENT_ID<ID,typename ID::VALUE,ID::capability_flags>,RW>
{};

template<class ID,class T,int flags> inline std::ostream& operator<<(std::ostream& output,const ELEMENT_ID<ID,T,flags> id)
{return output<<id.Value();}

template<class ID,class T,int flags> inline std::istream& operator>>(std::istream& input,ELEMENT_ID<ID,T,flags>& id)
{T i;input>>i;id=ID(i);return input;}
}
#endif
#endif

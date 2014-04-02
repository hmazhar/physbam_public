//#####################################################################
// Copyright 2009, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_SYMBOL
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_SYMBOL__
#define __READ_WRITE_SYMBOL__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Tools/Utilities/SYMBOL.h>
namespace PhysBAM{

template<class RW>
class Read_Write<SYMBOL,RW>:public SYMBOL
{
public:
    static void Read(std::istream& input,SYMBOL& object)
    {std::string name;Read_Binary<RW>(input,name);object.key=Mapping().Name_To_Key(name);}

    static void Write(std::ostream& output,const SYMBOL& object)
    {Write_Binary<RW>(output,object.Name());}
};
//#####################################################################
// Stream input and output
//#####################################################################
inline std::ostream& operator<<(std::ostream& output,const SYMBOL symbol)
{return output<<symbol.Name();}
}
#endif
#endif

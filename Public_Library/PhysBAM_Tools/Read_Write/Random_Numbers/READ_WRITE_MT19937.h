//#####################################################################
// Copyright 2009, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_MT19937
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_MT19937__
#define __READ_WRITE_MT19937__

#include <PhysBAM_Tools/Random_Numbers/MT19937.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<MT19937<T>,RW>
{
public:
    static void Read(std::istream& input,MT19937<T>& object)
    {Read_Binary<RW>(input,object.index);Read_Binary_Array<RW>(input,object.mt,object.n);}

    static void Write(std::ostream& output,const MT19937<T>& object)
    {Write_Binary<RW>(output,object.index);Write_Binary_Array<RW>(output,object.mt,object.n);}
};
//#####################################################################
// Stream input and output
//#####################################################################
template<class T> inline std::istream& 
operator>>(std::istream& input,MT19937<T>& g)
{input>>g.index;for(int i=0;i<g.n;i++) input>>g.mt[i];return input;}

template<class T> inline std::ostream&
operator<<(std::ostream& output,const MT19937<T>& g)
{output<<g.index;for(int i=0;i<g.n;i++) output<<" "<<g.mt[i];return output;}
}
#endif
#endif

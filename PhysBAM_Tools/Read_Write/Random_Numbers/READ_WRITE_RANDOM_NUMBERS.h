//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_RANDOM_NUMBERS
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_RANDOM_NUMBERS__
#define __READ_WRITE_RANDOM_NUMBERS__

#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Random_Numbers/READ_WRITE_MT19937.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T,class GENERATOR>
class Read_Write<RANDOM_NUMBERS<T,GENERATOR>,RW>
{
public:
    static void Read(std::istream& input,RANDOM_NUMBERS<T,GENERATOR>& object)
    {Read_Binary<RW>(input,object.gaussian_iset,object.gset);input>>object.random_number_generator;}

    static void Write(std::ostream& output,const RANDOM_NUMBERS<T,GENERATOR>& object)
    {Write_Binary<RW>(output,object.gaussian_iset,object.gset);output<<object.random_number_generator;}
};
}
#endif
#endif

//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_FRAME
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_FRAME__
#define __READ_WRITE_FRAME__

#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_ROTATION.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<FRAME<TV>,RW,typename DISABLE_IF<IS_BINARY_IO_SAFE<FRAME<TV>,RW>::value>::TYPE>
{
public:
    static void Read(std::istream& input,FRAME<TV>& object)
    {Read_Binary<RW>(input,object.t,object.r);}

    static void Write(std::ostream& output,const FRAME<TV>& object)
    {Write_Binary<RW>(output,object.t,object.r);}
};
// global functions
template<class TV> inline std::istream& operator>>(std::istream& input,FRAME<TV>& f)
{FILE_UTILITIES::Ignore(input,'(');input>>f.t>>f.r;FILE_UTILITIES::Ignore(input,')');return input;}

template<class TV> inline std::ostream& operator<<(std::ostream& output,const FRAME<TV>& f)
{output<<"("<<f.t<<"  "<<f.r<<")";return output;}
}
#endif
#endif

//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_RANGE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_RANGE__
#define __READ_WRITE_RANGE__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<RANGE<TV>,RW>
{
public:
    static void Read(std::istream& input,RANGE<TV>& object)
    {for(int i=1;i<=TV::dimension;i++) Read_Binary<RW>(input,object.min_corner(i),object.max_corner(i));}

    static void Write(std::ostream& output,const RANGE<TV>& object)
    {for(int i=1;i<=TV::dimension;i++) Write_Binary<RW>(output,object.min_corner(i),object.max_corner(i));}
};
template<class TV>
inline std::ostream& operator<<(std::ostream& output,const RANGE<TV>& box)
{output<<"("<<box.min_corner<<" "<<box.max_corner<<")";return output;}

template<class TV>
inline std::istream& operator>>(std::istream& input,RANGE<TV>& box)
{FILE_UTILITIES::Ignore(input,'(');input>>box.min_corner>>box.max_corner;FILE_UTILITIES::Ignore(input,')');return input;}
}
#endif
#endif

//#####################################################################
// Copyright 3009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_VECTOR_BASE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_VECTOR_BASE__
#define __READ_WRITE_VECTOR_BASE__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_BASE.h>
namespace PhysBAM{

template<class T,class T_VECTOR> void
Write_Raw(std::ostream& output,const VECTOR_BASE<T,T_VECTOR>& a)
{int m=a.Size();for(int i=1;i<=m;i++){output<<a(i);if(i<m) output<<" ";}}

template<class T,class T_VECTOR> inline std::ostream&
operator<<(std::ostream& output,const VECTOR_BASE<T,T_VECTOR>& a)
{output<<"[";Write_Raw(output,a);output<<"]";return output;}}
#endif
#endif

//#####################################################################
// Copyright 2007, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function argmin
//#####################################################################
#ifndef __argmin__
#define __argmin__

namespace PhysBAM{

template<class T>
inline int argmin(const T a,const T b)
{return a<=b?1:2;}

template<class T>
inline int argmin(const T a,const T b,const T c)
{if(a<=c) return a<=b?1:2;return b<=c?2:3;}

}
#endif

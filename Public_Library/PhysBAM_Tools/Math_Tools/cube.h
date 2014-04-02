//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function cube 
//#####################################################################
//
// Finds the cube. That is raises is to the cube power.  
//               
//#####################################################################
#ifndef __cube__
#define __cube__

namespace PhysBAM{

template<class T>
inline T cube(const T a)
{return a*a*a;}

}
#endif


//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Igor Neverov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function Minmod  
//#####################################################################
//
//  Returns the flux. 
//               
//#####################################################################
#ifndef __Minmod__
#define __Minmod__

#include <PhysBAM_Tools/Math_Tools/minmag.h>
namespace PhysBAM{

template<class T>
inline T Minmod(const T dx,const T D1,const T D2_left,const T D2_right)
{if(D2_left*D2_right>0) return D1+minmag(D2_left,D2_right)*dx;
else return D1;}

}
#endif


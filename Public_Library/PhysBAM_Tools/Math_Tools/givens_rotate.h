//#####################################################################
// Copyright 2006, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function givens_rotate
//#####################################################################
//
// Applies a givens rotation to a pair of scalars
//               
//#####################################################################
#ifndef __givens_rotate__
#define __givens_rotate__

namespace PhysBAM{

template<class T>
inline void givens_rotate(T& x,T& y,const T c,const T s)
{T w=c*x+s*y;y=c*y-s*x;x=w;}

}
#endif


//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function sign 
//#####################################################################
//
// finds the sign as +1, -1, or 0  
//               
//#####################################################################
#ifndef __sign__
#define __sign__

namespace PhysBAM{

template<class T>
inline T sign(const T a)
{if(a>0) return (T)1;
else if(a<0) return -(T)1;
else return (T)0;}

template<class T>
inline T sign_nonzero(const T a)
{return a>=0?(T)1:-(T)1;}

}
#endif


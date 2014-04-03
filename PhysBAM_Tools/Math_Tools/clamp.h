//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function clamp
//#####################################################################
#ifndef __clamp__
#define __clamp__

namespace PhysBAM{

template<class T> 
inline T clamp(const T x,const T xmin,const T xmax)
{if(x<=xmin) return xmin;
else if(x>=xmax) return xmax;
else return x;}

template<class T>
inline T clamp_min(const T x,const T xmin)
{if(x<=xmin) return xmin;else return x;}

template<class T>
inline T clamp_max(const T x,const T xmax)
{if(x>=xmax) return xmax;else return x;}

template<class T>
inline bool in_bounds(const T x,const T xmin,const T xmax)
{return xmin<=x && x<=xmax;}

}
#endif

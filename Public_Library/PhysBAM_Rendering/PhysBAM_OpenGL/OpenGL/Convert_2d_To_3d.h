//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function Convert_2d_To_3d
//#####################################################################
#ifndef __Convert_2d_To_3d__
#define __Convert_2d_To_3d__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Vectors/COMPLEX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T> inline VECTOR<T,3> Convert_2d_To_3d(const VECTOR<T,2>& v)
{return VECTOR<T,3>(v.x,v.y,0);}

template<class T> inline RANGE<VECTOR<T,3> > Convert_2d_To_3d(const RANGE<VECTOR<T,2> >& box)
{return RANGE<VECTOR<T,3> >(box.min_corner.x,box.max_corner.x,box.min_corner.y,box.max_corner.y,0,0);}

template<class T> inline ROTATION<VECTOR<T,3> > Convert_2d_To_3d(const ROTATION<VECTOR<T,2> >& c)
{COMPLEX<T> sqrt_c=c.Complex().Sqrt();return ROTATION<VECTOR<T,3> >::From_Components(sqrt_c.re,0,0,sqrt_c.im);}

template<class T> inline FRAME<VECTOR<T,3> > Convert_2d_To_3d(const FRAME<VECTOR<T,2> >& f)
{return FRAME<VECTOR<T,3> >(Convert_2d_To_3d(f.t),Convert_2d_To_3d(f.r));}

}
#endif

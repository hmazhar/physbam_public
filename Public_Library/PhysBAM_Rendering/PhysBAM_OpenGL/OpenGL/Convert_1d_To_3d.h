//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function Convert_1d_To_3d
//#####################################################################
#ifndef __Convert_1d_To_3d__
#define __Convert_1d_To_3d__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T> inline VECTOR<T,3> Convert_1d_To_3d(const VECTOR<T,1>& v)
{return VECTOR<T,3>(v.x,0,0);}

template<class T> inline RANGE<VECTOR<T,3> > Convert_1d_To_3d(const RANGE<VECTOR<T,1> >& box)
{return RANGE<VECTOR<T,3> >(box.min_corner.x,box.max_corner.x,0,0,0,0);}

template<class T> inline ROTATION<VECTOR<T,3> > Convert_1d_To_3d(const ROTATION<VECTOR<T,1> >& c)
{return ROTATION<VECTOR<T,3> >();}

template<class T> inline FRAME<VECTOR<T,3> > Convert_1d_To_3d(const FRAME<VECTOR<T,1> >& f)
{return FRAME<VECTOR<T,3> >(Convert_1d_To_3d(f.t),Convert_1d_To_3d(f.r));}

}
#endif

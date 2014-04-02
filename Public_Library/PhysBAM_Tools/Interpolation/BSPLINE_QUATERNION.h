//#####################################################################
// Copyright 2005, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BSPLINE_QUATERNION
//#####################################################################
#ifndef __BSPLINE_QUATERNION__
#define __BSPLINE_QUATERNION__

#include <PhysBAM_Tools/Interpolation/BSPLINE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
namespace PhysBAM{

template<class T>
class BSPLINE_QUATERNION:public BSPLINE<T,ROTATION<VECTOR<T,3> > >
{
public:
    typedef VECTOR<T,3> TV;typedef BSPLINE<T,ROTATION<TV> > BASE;
    using BASE::control_points_times;using BASE::control_points;using BASE::k;using BASE::closed;
    using BASE::Start_Time;using BASE::End_Time;using BASE::Range;using BASE::Create_Closed_Points;

    BSPLINE_QUATERNION(const ARRAY<T>& control_points_times,const ARRAY<ROTATION<TV> >& control_points,const int order=1);
    ~BSPLINE_QUATERNION();

    ROTATION<TV>  Evaluate(const T t) PHYSBAM_OVERRIDE;
    T Quaternion_Basis_Function(const int i,const int k,const T t);
    ROTATION<TV> Omega(const int i);
    void Create_Closed_Points() PHYSBAM_OVERRIDE;
    void Quaternion_Check();
//#####################################################################
};
}
#endif

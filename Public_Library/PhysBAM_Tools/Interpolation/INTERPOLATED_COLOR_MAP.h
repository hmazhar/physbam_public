//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Macro INTERPOLATED_COLOR_MAP
//#####################################################################
#ifndef __INTERPOLATED_COLOR_MAP__
#define __INTERPOLATED_COLOR_MAP__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
namespace PhysBAM{
template<class T>
class INTERPOLATED_COLOR_MAP
{
public:
    T mn,mx;
    bool use_log;
    INTERPOLATION_CURVE<T,VECTOR<T,3> >& colors;

    INTERPOLATED_COLOR_MAP();
    ~INTERPOLATED_COLOR_MAP();

    void Initialize_Colors(T min_value,T max_value,bool log_scale,bool reverse,bool two_sets);

    VECTOR<T,3> operator()(T x) const;
};
}

#endif

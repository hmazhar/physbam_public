//#####################################################################
// Copyright 2002-2003, Ronald Fedkiw, Eran Guendelman, Neil Molino, Igor Neverov, Robert Bridson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_NORMAL_COMPUTATION
//#####################################################################
#ifndef __LEVELSET_NORMAL_COMPUTATION__
#define __LEVELSET_NORMAL_COMPUTATION__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
namespace PhysBAM{

template<class T,int d> class VECTOR;
template<class T_GRID> class LEVELSET_3D;

template<class T_GRID>
class LEVELSET_NORMAL_COMPUTATION
{
    typedef typename T_GRID::SCALAR T;
public:
//#####################################################################
    virtual ~LEVELSET_NORMAL_COMPUTATION(){}
    virtual VECTOR<T,3> Compute_Normal(const LEVELSET_3D<T_GRID>& levelset,const VECTOR<T,3>& location) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2006-2009, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_POLICY 
//#####################################################################
#ifndef __LEVELSET_POLICY__
#define __LEVELSET_POLICY__
namespace PhysBAM{

template<class T_GRID> class LEVELSET_POLICY;
template<class T_GRID> class FAST_LEVELSET;
template<class T,int d> class VECTOR;

enum PARTICLE_LEVELSET_PARTICLE_TYPE {PARTICLE_LEVELSET_POSITIVE,PARTICLE_LEVELSET_NEGATIVE,PARTICLE_LEVELSET_REMOVED_POSITIVE,PARTICLE_LEVELSET_REMOVED_NEGATIVE};

}
#endif

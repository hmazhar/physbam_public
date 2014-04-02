//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Geoffrey Irving, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STRUCTURE_REGISTRY
//##################################################################### 
#ifndef __STRUCTURE_REGISTRY__
#define __STRUCTURE_REGISTRY__

#include <PhysBAM_Geometry/Registry/REGISTRY.h>
namespace PhysBAM{

template<class TV> class STRUCTURE;  // Needed before the declaration of STRUCTURE_REGISTRY
template<class TV> class STRUCTURE_REGISTRY:public REGISTRY<STRUCTURE<TV>,std::string,STRUCTURE_REGISTRY<TV> > {
public:
    STRUCTURE_REGISTRY();
    virtual ~STRUCTURE_REGISTRY();
};
}
#endif

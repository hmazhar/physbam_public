//#####################################################################
// Copyright 2005-2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_WORLD
//#####################################################################
#ifndef __MPI_WORLD__
#define __MPI_WORLD__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

class MPI_WORLD:public NONCOPYABLE
{
public:
    bool initialized;
    int rank;

    MPI_WORLD();
    MPI_WORLD(int& argc,char**& argv);
    ~MPI_WORLD();

//#####################################################################
    static bool Initialized();
private:
    void Initialize(bool force_mpi);
//#####################################################################
};
}
#endif

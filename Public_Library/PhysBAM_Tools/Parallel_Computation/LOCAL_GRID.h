//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LOCAL_GRID
//#####################################################################
#ifndef __LOCAL_GRID__
#define __LOCAL_GRID__

#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
namespace PhysBAM{

template<class T_GRID>
class LOCAL_GRID
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
public:
    T_GRID grid;
    MPI_UNIFORM_GRID<T_GRID> mpi_grid;
    const T_GRID& global_grid;
    TV_INT offset;
    ARRAY<bool> neighbor_overlaps;

    LOCAL_GRID(const T_GRID& global_grid_input);
    LOCAL_GRID(const T_GRID& global_grid_input,const T_GRID& local_grid_input);
    ~LOCAL_GRID();

    RANGE<TV_INT> Interior_Region(const RANGE<TV_INT>& sentinels) const
    {return grid.Get_MAC_Grid().Domain_Indices()+sentinels;}

//#####################################################################
    void Initialize();
    template<class T_ARRAYS> void Put(const T_ARRAYS& local_data,const RANGE<TV_INT>& region,T_ARRAYS& global_data) const;
    template<class T_ARRAYS> void Put(const T_ARRAYS& local_data,T_ARRAYS& global_data,const RANGE<TV_INT>& sentinels=RANGE<TV_INT>::Zero_Box()) const;
    template<class T_ARRAYS> void Get(const T_ARRAYS& global_data,T_ARRAYS& local_data) const;
    template<class T_FACE_ARRAYS> void Put_Faces(const T_FACE_ARRAYS& local_data,T_FACE_ARRAYS& global_data) const;
    template<class T_ARRAYS> T Maximum_Error(const T_ARRAYS& local_data,const T_ARRAYS& global_data,const int bandwidth,TV_INT& index,const RANGE<TV_INT>& sentinels=RANGE<TV_INT>::Zero_Box()) const;
    template<class T_FACE_ARRAYS> T Maximum_Error(const std::string& prefix,const T_FACE_ARRAYS& local_data,const T_FACE_ARRAYS& global_data,const int bandwidth,const T threshold=1e-7);
//#####################################################################
};
}
#endif

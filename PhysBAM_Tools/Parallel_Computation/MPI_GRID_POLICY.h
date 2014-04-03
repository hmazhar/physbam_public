//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPI_GRID_POLICY__
#define __MPI_GRID_POLICY__

namespace PhysBAM{

template<class T_GRID> class MPI_UNIFORM_GRID;

template<class T_GRID>
struct MPI_GRID_POLICY
{
    typedef MPI_UNIFORM_GRID<T_GRID> MPI_GRID;
    typedef T_GRID PARALLEL_GRID;
};

#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template<class T_GRID> class MPI_RLE_GRID;
template<class T> class RLE_GRID_2D;
template<class T> class RLE_GRID_3D;

template<class T> 
struct MPI_GRID_POLICY<RLE_GRID_2D<T> >
{
    typedef MPI_RLE_GRID<RLE_GRID_2D<T> > MPI_GRID;
    typedef typename RLE_GRID_2D<T>::HORIZONTAL_GRID PARALLEL_GRID;
};

template<class T> 
struct MPI_GRID_POLICY<RLE_GRID_3D<T> >
{
    typedef MPI_RLE_GRID<RLE_GRID_3D<T> > MPI_GRID;
    typedef typename RLE_GRID_3D<T>::HORIZONTAL_GRID PARALLEL_GRID;
};
#endif
}
#endif

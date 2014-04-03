//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __LAPLACE_COLLIDABLE_POLICY__
#define __LAPLACE_COLLIDABLE_POLICY__

namespace PhysBAM{
template<class T_GRID> class LAPLACE_COLLIDABLE_UNIFORM;
template<class T_GRID> class LAPLACE_COLLIDABLE_RLE;
template<class T> class RLE_GRID_2D;
template<class T> class RLE_GRID_3D;

template<class T_GRID>
struct LAPLACE_COLLIDABLE_POLICY
{
    typedef LAPLACE_COLLIDABLE_UNIFORM<T_GRID> LAPLACE;
};

template<class T>
struct LAPLACE_COLLIDABLE_POLICY<RLE_GRID_2D<T> >
{
    typedef LAPLACE_COLLIDABLE_RLE<RLE_GRID_2D<T> > LAPLACE;
};

template<class T>
struct LAPLACE_COLLIDABLE_POLICY<RLE_GRID_3D<T> >
{
    typedef LAPLACE_COLLIDABLE_RLE<RLE_GRID_3D<T> > LAPLACE;
};

}
#endif

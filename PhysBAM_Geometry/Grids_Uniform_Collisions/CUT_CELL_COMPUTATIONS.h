//#####################################################################
// Copyright 2011, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUT_CELL_COMPUTATIONS
//#####################################################################
#ifndef __CUT_CELL_COMPUTATIONS__
#define __CUT_CELL_COMPUTATIONS__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/CUT_CELL.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
namespace PhysBAM {
namespace CUT_CELL_COMPUTATIONS{
template<class T>
void Compute_Cut_Geometries(const GRID<VECTOR<T,1> >& grid,const int num_ghost_cells,typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<T,1> > >::GRID_BASED_COLLISION_GEOMETRY& collision_bodies_affecting_fluid,ARRAY<CUT_CELLS<T,1>*,VECTOR<int,1> >& cut_cells);

template<class T>
void Compute_Cut_Geometries(const GRID<VECTOR<T,2> >& grid,const int num_ghost_cells,typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<T,2> > >::GRID_BASED_COLLISION_GEOMETRY& collision_bodies_affecting_fluid,ARRAY<CUT_CELLS<T,2>*,VECTOR<int,2> >& cut_cells);

template<class T>
void Compute_Cut_Geometries(const GRID<VECTOR<T,3> >& grid,const int num_ghost_cells,typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<T,3> > >::GRID_BASED_COLLISION_GEOMETRY& collision_bodies_affecting_fluid,ARRAY<CUT_CELLS<T,3>*,VECTOR<int,3> >& cut_cells);
}
}
#endif

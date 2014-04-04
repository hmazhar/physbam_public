#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/INTERPOLATION_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_3D_HELPER.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <PhysBAM_Geometry/Grids_RLE_Computations/DUALCONTOUR_RLE_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include "RLE_IMPLICIT_SURFACE_TWO_LEVEL.h"
using namespace PhysBAM;
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* RLE_IMPLICIT_SURFACE_TWO_LEVEL<T>::
Generate_Triangles() const
{
    // TODO: get rid of const casts somehow
    TRIANGULATED_SURFACE<T> *surface=DUALCONTOUR_RLE_3D<T>::Create_Triangulated_Surface_From_Levelset(const_cast<LEVELSET_RLE<T_GRID>&>(levelset),levelset.grid.number_of_ghost_cells),
        *fine_surface=DUALCONTOUR_RLE_3D<T>::Create_Triangulated_Surface_From_Levelset(const_cast<LEVELSET_RLE<T_GRID>&>(fine_levelset),fine_levelset.grid.number_of_ghost_cells);
    int fine_particles=fine_surface->particles.array_collection->Size();
    fine_surface->particles.array_collection->Take(*surface->particles.array_collection);
    fine_surface->mesh.number_nodes=fine_surface->particles.array_collection->Size();
    int fine_triangles=fine_surface->mesh.elements.m;
    fine_surface->mesh.elements.Resize(fine_triangles+surface->mesh.elements.m);
    for(int t=1;t<=surface->mesh.elements.m;t++){
        int i,j,k;surface->mesh.elements(t).Get(i,j,k);
        fine_surface->mesh.elements(fine_triangles+t).Set(fine_particles+i,fine_particles+j,fine_particles+k);}
    delete surface;
    return fine_surface;
}
//#####################################################################
template class RLE_IMPLICIT_SURFACE_TWO_LEVEL<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RLE_IMPLICIT_SURFACE_TWO_LEVEL<double>;
#endif
#endif

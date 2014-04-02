//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_IMPLICIT_OBJECT
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/LINEAR_INTERPOLATION_RLE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Geometry/Implicit_Objects_RLE/RLE_IMPLICIT_OBJECT.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RLE_IMPLICIT_OBJECT<TV>::
RLE_IMPLICIT_OBJECT(T_GRID& grid_input,ARRAY<T>& phi_input)
    :levelset(grid_input,phi_input)
{
    Initialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RLE_IMPLICIT_OBJECT<TV>::
~RLE_IMPLICIT_OBJECT()
{}
//###########################################################################
template RLE_IMPLICIT_OBJECT<VECTOR<float,2> >::RLE_IMPLICIT_OBJECT(RLE_GRID_2D<float>&,ARRAY<float>&);
template RLE_IMPLICIT_OBJECT<VECTOR<float,3> >::RLE_IMPLICIT_OBJECT(RLE_GRID_3D<float>&,ARRAY<float>&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template RLE_IMPLICIT_OBJECT<VECTOR<double,2> >::RLE_IMPLICIT_OBJECT(RLE_GRID_2D<double>&,ARRAY<double>&);
template RLE_IMPLICIT_OBJECT<VECTOR<double,3> >::RLE_IMPLICIT_OBJECT(RLE_GRID_3D<double>&,ARRAY<double>&);
#endif
#endif

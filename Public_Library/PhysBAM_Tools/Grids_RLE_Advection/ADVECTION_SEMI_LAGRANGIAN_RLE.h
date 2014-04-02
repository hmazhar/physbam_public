//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_RLE
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __ADVECTION_SEMI_LAGRANGIAN_RLE__
#define __ADVECTION_SEMI_LAGRANGIAN_RLE__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <cassert>
namespace PhysBAM{

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> // T_AVERAGING=AVERAGING_RLE<T,T,T_GRID>, T_INTERPOLATION=LINEAR_INTERPOLATION_RLE<T_GRID,T2>
class ADVECTION_SEMI_LAGRANGIAN_RLE:public ADVECTION<T_GRID,T2,typename T_AVERAGING::FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::BLOCK T_BLOCK;typedef typename T_AVERAGING::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;typedef typename REBIND<T_BOUNDARY,T2>::TYPE T_BOUNDARY_T2;
public:
    template<class T3> struct REBIND{typedef ADVECTION_SEMI_LAGRANGIAN_RLE<T_GRID,T3,T_AVERAGING,typename T_INTERPOLATION::template REBIND<T3>::TYPE> TYPE;};
    template<class T_INTERPOLATION_2> struct REBIND_INTERPOLATION{typedef ADVECTION_SEMI_LAGRANGIAN_RLE<T_GRID,T2,T_AVERAGING,T_INTERPOLATION_2> TYPE;};

    ADVECTION_SEMI_LAGRANGIAN_RLE()
    {}

    void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,ARRAY<T2>& Z,const ARRAY<T2>& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const ARRAY<T2>* Z_min_ghost,const ARRAY<T2>* Z_max_ghost,ARRAY<T2>* Z_min,ARRAY<T2>* Z_max)
    {assert(!Z_min && !Z_max);
    T_AVERAGING averaging;T_INTERPOLATION interpolation;
    for(CELL_ITERATOR cell(grid,0);cell;cell++)if(cell.Short()){
        TV sample_location=cell.X()-dt*averaging.Face_To_Cell_Vector(cell,face_velocities);
        T_BLOCK block(grid,sample_location);if(!block) continue;
        Z(cell.Cell())=interpolation.From_Block_Cell(block,Z_ghost,sample_location);}}

    void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,ARRAY<T>& Z,const T_FACE_LOOKUP& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,ARRAY<T>* Z_min,ARRAY<T>* Z_max)
    {assert(!Z_min && !Z_max);
    T_GRID::template Face_Loop<Update_Advection_Equation_Face_Helper>(grid,Z,Z_ghost,face_velocities,dt);}

private:
    struct Update_Advection_Equation_Face_Helper{template<class T_FACE> static void
    Apply(const T_GRID& grid,ARRAY<T>& Z,const T_FACE_LOOKUP& Z_ghost,const T_FACE_LOOKUP& face_velocities,const T dt)
    {T_AVERAGING averaging;T_INTERPOLATION interpolation;int axis=T_FACE::Axis();
    for(T_FACE face(grid,0);face;face++)if(face.Both_Cells_Short()){
        TV sample_location=face.X()-dt*averaging.Face_To_Face_Vector(face,face_velocities);
        T_BLOCK block(grid,sample_location);if(!block) continue;
        Z(face.Face())=interpolation.From_Block_Face_Component(axis,block,Z_ghost.Starting_Point_Face(face),sample_location);}}};

//#####################################################################
};
}
#endif
#endif

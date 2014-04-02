//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_REFLECTION_ATTENUATION
//#####################################################################
// Attenuates to far field value for open walls (constant_extrapolation=true) and reflects for solid walls
//#####################################################################
#ifndef __BOUNDARY_REFLECTION_ATTENUATION__
#define __BOUNDARY_REFLECTION_ATTENUATION__

#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_REFLECTION_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T_GRID,class T2>
class BOUNDARY_REFLECTION_ATTENUATION:public BOUNDARY_REFLECTION_UNIFORM<T_GRID,T2>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef VECTOR<bool,T_GRID::dimension> TV_BOOL;
    typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,T_GRID::dimension> TV_SIDES;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;
    typedef typename T_ARRAYS_BASE::template REBIND<T2>::TYPE T_ARRAYS_T2;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
public:
    typedef BOUNDARY_REFLECTION_UNIFORM<T_GRID,T2> BASE;
    using BASE::Set_Constant_Extrapolation;using BASE::Constant_Extrapolation;using BASE::Fill_Single_Ghost_Region;
    using BASE::fixed_boundary_value;

    T linear_attenuation;

    BOUNDARY_REFLECTION_ATTENUATION(const TV_SIDES& constant_extrapolation,const T2 far_field_value,const T linear_attenuation_input)
    {
        linear_attenuation=linear_attenuation_input;
        Set_Constant_Extrapolation(constant_extrapolation);
        fixed_boundary_value=far_field_value;
    }

//#####################################################################
    virtual void Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) PHYSBAM_OVERRIDE;
    virtual void Fill_Single_Ghost_Region(const T_GRID& grid,T_ARRAYS_T2& u_ghost,const RANGE<TV_INT>& region,const int side,const T dt,const T time,const int number_of_ghost_cells=3) const PHYSBAM_OVERRIDE ;
    virtual void Apply_Boundary_Condition(const T_GRID& grid,T_ARRAYS_T2& u,const T time) PHYSBAM_OVERRIDE {} // do nothing
private:
    T2 Attenuate_To_Far_Field_Value(const T2 boundary_value,const T dt) const;
//#####################################################################
};
}
#endif

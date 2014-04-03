//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SCATTERED_INTERPOLATION
//#####################################################################
#ifndef __SCATTERED_INTERPOLATION__
#define __SCATTERED_INTERPOLATION__

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
namespace PhysBAM{

template<class T_GRID>
class SCATTERED_INTERPOLATION
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename REBIND<T_ARRAYS_SCALAR,ARRAY<int> >::TYPE T_ARRAYS_ARRAY_INT;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
public:
    T radius_of_influence;
    bool use_tent_weights;
    bool use_distance_averaged_weights;
private:
    T radius_of_influence_squared;
public:

    SCATTERED_INTERPOLATION();
    ~SCATTERED_INTERPOLATION();

    void Set_Radius_Of_Influence(const T radius_of_influence_input)
    {radius_of_influence=radius_of_influence_input;radius_of_influence_squared=sqr(radius_of_influence);}

    void Use_Tent_Weights(const bool use=true)
    {use_tent_weights=use;if(use_tent_weights) use_distance_averaged_weights=false;}

    void Use_Distance_Averaged_Weights(const bool use=true)
    {use_distance_averaged_weights=use;if(use_distance_averaged_weights) use_tent_weights=false;}

private:
    RANGE<TV_INT> Grid_Influence_Bounds(const T_GRID& grid,const TV& X) const
    {return RANGE<TV_INT>(grid.Clamp_To_Cell(X-radius_of_influence),grid.Clamp_To_Cell(X+radius_of_influence));}
public:

    template<class T_ARRAY_T2,class T_ARRAYS_T2> void Transfer_To_Grid(ARRAY_VIEW<const TV> domain,const T_ARRAY_T2& range,const T_GRID& grid,T_ARRAYS_T2& grid_data) const
    {Transfer_To_Grid_Helper<T_ARRAYS_T2>(domain,range,grid,grid_data);}

//#####################################################################
    template<class T_ARRAYS_T2> void Transfer_To_Grid_Helper(ARRAY_VIEW<const TV> domain,ARRAY_VIEW<const typename T_ARRAYS_T2::ELEMENT> range,const T_GRID& grid,T_ARRAYS_T2& grid_data) const;
private:
    static void Bin_Domain_Values(ARRAY_VIEW<const TV> domain_values,const T_GRID& grid,T_ARRAYS_ARRAY_INT& points_in_cell);
    template<class T_ARRAYS_T2> void Transfer_With_Distance_Averaged_Weights(ARRAY_VIEW<const TV> domain,ARRAY_VIEW<const typename T_ARRAYS_T2::ELEMENT> range,const T_GRID& grid,
        T_ARRAYS_T2& grid_data,const T_ARRAYS_ARRAY_INT& points_in_cell) const;
    template<class T_ARRAYS_T2> void Transfer_With_Tent_Weights(ARRAY_VIEW<const TV> domain,ARRAY_VIEW<const typename T_ARRAYS_T2::ELEMENT> range,const T_GRID& grid,T_ARRAYS_T2& grid_data,
        const T_ARRAYS_ARRAY_INT& points_in_cell) const;
    void Instantiate();
//#####################################################################
};
}
#endif

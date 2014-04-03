//#####################################################################
// Copyright 2005-2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_IMPLICIT_OBJECT
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_IMPLICIT_OBJECT
#define __RLE_IMPLICIT_OBJECT

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_POLICY.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_RLE.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class TV>
class RLE_IMPLICIT_OBJECT:public IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename RLE_GRID_POLICY<TV>::RLE_GRID T_GRID;
    typedef typename INTERPOLATION_POLICY<T_GRID>::INTERPOLATION_SCALAR T_INTERPOLATION_SCALAR;
public:
    typedef IMPLICIT_OBJECT<TV> BASE;
    using BASE::box;

    LEVELSET_RLE<T_GRID> levelset;
private:
    T minimum_cell_size;
public:

    RLE_IMPLICIT_OBJECT(T_GRID& grid_input,ARRAY<T>& phi_input);
    ~RLE_IMPLICIT_OBJECT();

    void Initialize()
    {Update_Box();Update_Minimum_Cell_Size();}

    void Set_Custom_Interpolation(T_INTERPOLATION_SCALAR& interpolation_input)
    {levelset.Set_Custom_Interpolation(interpolation_input);}

    void Update_Box() PHYSBAM_OVERRIDE
    {box=levelset.grid.domain;}

    void Update_Minimum_Cell_Size(const int maximum_depth=0) PHYSBAM_OVERRIDE
    {minimum_cell_size=levelset.grid.Minimum_Edge_Length();}

    T operator()(const TV& location) const PHYSBAM_OVERRIDE
    {return levelset.Phi_Far(location);}

    T Extended_Phi(const TV& location) const PHYSBAM_OVERRIDE
    {return levelset.Extended_Phi(location);}

    T Phi_Secondary(const TV& location) const PHYSBAM_OVERRIDE
    {return levelset.Phi_Secondary(location);}

    TV Normal(const TV& location,const int aggregate=-1) const PHYSBAM_OVERRIDE
    {assert((aggregate >= 1 && aggregate <= T_GRID::number_of_neighbors_per_cell) || aggregate == -1);
    if(aggregate != -1) return box.Normal(aggregate);else return levelset.Normal(location);}

    TV Extended_Normal(const TV& location,const int aggregate=-1) const PHYSBAM_OVERRIDE
    {assert((aggregate >= 1 && aggregate <= T_GRID::number_of_neighbors_per_cell) || aggregate == -1);
    if(aggregate != -1) return box.Normal(aggregate);else return levelset.Extended_Normal(location);}

    void Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists=true) PHYSBAM_OVERRIDE
    {levelset.Compute_Block_Minimum_And_Maximum(recompute_if_exists);}

    // TODO: an rle specific intersection routine could take much larger steps that this
    T Integration_Step(const T phi) const PHYSBAM_OVERRIDE
    {return max(abs(phi),(T).1*minimum_cell_size);}

    T Minimum_Cell_Size() const  PHYSBAM_OVERRIDE
    {return minimum_cell_size;}

    T Min_Phi() const PHYSBAM_OVERRIDE
    {return levelset.phi.Min();}

//###########################################################################
};
}
#endif
#endif

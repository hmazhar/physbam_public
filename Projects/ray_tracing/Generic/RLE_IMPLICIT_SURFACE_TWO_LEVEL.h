#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_IMPLICIT_SURFACE_TWO_LEVEL
//#####################################################################
#ifndef __RLE_IMPLICIT_SURFACE_TWO_LEVEL__
#define __RLE_IMPLICIT_SURFACE_TWO_LEVEL__

#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_RLE.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class T_input>
class RLE_IMPLICIT_SURFACE_TWO_LEVEL:public IMPLICIT_OBJECT<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef RLE_GRID_3D<T> T_GRID;typedef typename T_GRID::BLOCK T_BLOCK;

    typedef IMPLICIT_OBJECT<TV> BASE;
    using BASE::box;

    LEVELSET_RLE<T_GRID> levelset;
    LEVELSET_RLE<T_GRID> fine_levelset;
    T minimum_cell_size;

    RLE_IMPLICIT_SURFACE_TWO_LEVEL(T_GRID& grid_input,ARRAY<T>& phi_input,T_GRID& fine_grid_input,ARRAY<T>& fine_phi_input)
        :levelset(grid_input,phi_input),fine_levelset(fine_grid_input,fine_phi_input)
    {
        Update_Box();Update_Minimum_Cell_Size();
        fine_levelset.Compute_Normals();
    }

    void Update_Box() PHYSBAM_OVERRIDE
    {box=levelset.grid.domain;}

    void Update_Minimum_Cell_Size(const int maximum_depth=0) PHYSBAM_OVERRIDE
    {minimum_cell_size=fine_levelset.grid.Minimum_Edge_Length();}

    T operator()(const VECTOR<T,3>& X) const PHYSBAM_OVERRIDE
    {T_BLOCK block(fine_levelset.grid,X);if(block) return fine_levelset.Phi(block,X);
    return levelset.Phi_Far(X);}

    T Extended_Phi(const VECTOR<T,3>& X) const PHYSBAM_OVERRIDE
    {T_BLOCK block(fine_levelset.grid,X);if(block) return fine_levelset.Phi(block,X);
    return levelset.Extended_Phi(X);}

    T Phi_Secondary(const VECTOR<T,3>& X) const PHYSBAM_OVERRIDE
    {T_BLOCK block(fine_levelset.grid,X);if(block) return fine_levelset.Phi(block,X);
    return levelset.Phi_Secondary(X);}

    VECTOR<T,3> Normal(const VECTOR<T,3>& X,const int aggregate=-1) const PHYSBAM_OVERRIDE
    {assert((aggregate >= 1 && aggregate <= 6) || aggregate == -1);
    if(aggregate != -1) return box.Normal(aggregate);
    T_BLOCK block(fine_levelset.grid,X);if(block) return fine_levelset.Normal(block,X);
    return levelset.Normal(X);}

    VECTOR<T,3> Extended_Normal(const VECTOR<T,3>& X,const int aggregate=-1) const PHYSBAM_OVERRIDE
    {assert((aggregate >= 1 && aggregate <= 6) || aggregate == -1);
    if(aggregate != -1) return box.Normal(aggregate);
    T_BLOCK block(fine_levelset.grid,X);if(block) return fine_levelset.Normal(block,X);
    return levelset.Extended_Normal(X);}

    // TODO: an rle specific intersection routine could take much larger steps that this
    T Integration_Step(const T phi) const PHYSBAM_OVERRIDE
    {T distance=fabs(phi);
    if(distance > 3*minimum_cell_size) return (T).5*distance;
    else if(distance > minimum_cell_size) return (T).25*distance;
    else return (T).1*minimum_cell_size;}

    T Minimum_Cell_Size() const PHYSBAM_OVERRIDE
    {return minimum_cell_size;}

//###########################################################################
    TRIANGULATED_SURFACE<T>* Generate_Triangles() const PHYSBAM_OVERRIDE;
//###########################################################################
};
}
#endif
#endif

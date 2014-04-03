//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXTRAPOLATION_RLE
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __EXTRAPOLATION_RLE__
#define __EXTRAPOLATION_RLE__

#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_RLE.h>
#include <PhysBAM_Geometry/Level_Sets/EXTRAPOLATION.h>
namespace PhysBAM{

template<class T_GRID,class T2>
class EXTRAPOLATION_RLE:public EXTRAPOLATION<T_GRID,T2>
{
    typedef typename T_GRID::SCALAR T;
    typedef ARRAY<VECTOR<int,T_GRID::number_of_neighbors_per_cell> > T_FACE_NEIGHBORS;
public:
    typedef EXTRAPOLATION<T_GRID,T2> BASE;
    using BASE::band_width;using BASE::boundary;using BASE::isobaric_fix_width;

private:
    const LEVELSET_RLE<T_GRID> levelset;
    bool use_default_value;
    T2 default_value;
    const ARRAY<VECTOR<bool,T_GRID::dimension> >* neighbors_visible; // used for collision aware extrapolation
    T dz_squared_over_dy_squared,dz_squared_over_dx_squared,dy_squared_over_dx_squared,tolerance; // optimizations (some of these are unused in 2d)
public:

    EXTRAPOLATION_RLE(const T_GRID& grid_input,ARRAY<T>& phi_input);
    ~EXTRAPOLATION_RLE();

    void Set_Band_Width(const T number_of_cells=(T)5)
    {band_width=number_of_cells*levelset.grid.uniform_grid.dX.Max();}

    void Set_Band_Width_To_Entire_Domain()
    {band_width=max(levelset.grid.negative_bandwidth,levelset.grid.positive_bandwidth);}

    void Set_Isobaric_Fix_Width(const int number_of_cells=0)
    {isobaric_fix_width=number_of_cells*levelset.grid.uniform_grid.dX.Min();}

    void Set_Default_Value(const T2& default_value_input)
    {use_default_value=true;default_value=default_value_input;}

    void Set_Collision_Aware_Extrapolation(const ARRAY<VECTOR<bool,T_GRID::dimension> >& neighbors_visible_input)
    {neighbors_visible=&neighbors_visible_input;}

private:
    bool Neighbor_Visible(const int i,const int node,const int neighbor) const
    {assert(1<=i && i<=T_GRID::number_of_neighbors_per_cell);
    return neighbor && (!neighbors_visible || (*neighbors_visible)(min(node,neighbor))(((i-1)>>1)+1));}
public:

//#####################################################################
    void Extrapolate_Cells(ARRAY<T2>& u,const T time);
    void Extrapolate_Faces(ARRAY<T2>& u,const ARRAY<bool>& fixed) const;
private:
    void Extrapolate(const T_FACE_NEIGHBORS& neighbors,ARRAY<T2>& u,const ARRAY<T>& phi,ARRAY<bool,VECTOR<int,1> >& done) const;
    void Update_Close_Point(const ARRAY<VECTOR<int,4> >& neighbors,ARRAY<T2>& u,const ARRAY<T>& phi,const ARRAY<bool,VECTOR<int,1> >& done,const int node) const;
    void Update_Close_Point(const ARRAY<VECTOR<int,6> >& neighbors,ARRAY<T2>& u,const ARRAY<T>& phi,const ARRAY<bool,VECTOR<int,1> >& done,const int node) const;
//#####################################################################
};
}
#endif
#endif

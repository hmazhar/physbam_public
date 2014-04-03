//#####################################################################
// Copyright 2003-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_PARTITION
//#####################################################################
#ifndef __PARTICLE_PARTITION__
#define __PARTICLE_PARTITION__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class TV>
class PARTICLE_PARTITION
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_T;
    typedef typename REBIND<T_ARRAYS_T,ARRAY<int> >::TYPE T_ARRAYS_ARRAY_INT;
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
public:
    GRID<TV> grid;
    T_ARRAYS_ARRAY_INT partition;
private:
    bool use_radius;
    T_ARRAYS_T radius;
public:

    PARTICLE_PARTITION(const RANGE<TV>& box,const VECTOR<int,d>& counts,const GEOMETRY_PARTICLES<TV>& particles,const bool use_radius=true,const bool is_mac_grid=true)
        :grid(counts,box,is_mac_grid),partition(grid.Domain_Indices()),use_radius(use_radius)
    {
        if(use_radius) radius.Resize(grid.Domain_Indices());
        for(int p=1;p<=particles.array_collection->Size();p++) Add_To_Partition(particles.X(p),p);
    }

    void Add_To_Partition(const TV& location,const int particle_id)
    {assert(grid.Domain().Lazy_Inside(location));
    VECTOR<int,d> cell=grid.Clamp_To_Cell(location);partition(cell).Append(particle_id);
    if(use_radius) radius(cell)=max(radius(cell),(location-grid.Center(cell)).Magnitude());}

    RANGE<VECTOR<int,d> > Range(const RANGE<TV>& box) const
    {return RANGE<VECTOR<int,d> >(grid.Clamp_To_Cell(box.Minimum_Corner()),grid.Clamp_To_Cell(box.Maximum_Corner()));}

    void Intersection_List(const IMPLICIT_OBJECT<TV>& test_surface,const MATRIX<T,d>& rotation,const TV& translation,ARRAY<VECTOR<int,d> >& intersection_list,const T contour_value=0) const
    {PHYSBAM_ASSERT(use_radius);intersection_list.Remove_All();
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){VECTOR<int,d> cell=iterator.Cell_Index(); 
        if(partition(cell).m && !test_surface.Lazy_Outside_Extended_Levelset(rotation*grid.Center(cell)+translation,radius(cell)+contour_value)) intersection_list.Append(cell);}}

    void Proximity_List(const TV& location,const T proximity,ARRAY<int>& proximity_list)
    {
        proximity_list.Remove_All();
        for(CELL_ITERATOR iterator(grid,grid.Clamp_To_Cell(RANGE<TV>(location).Thickened(proximity),0));iterator.Valid();iterator.Next())
            proximity_list.Append_Unique_Elements(partition(iterator.Cell_Index()));
    }

//#####################################################################
};
}
#endif

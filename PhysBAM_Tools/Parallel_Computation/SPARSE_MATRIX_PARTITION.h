//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPARSE_MATRIX_PARTITION
//#####################################################################
#ifndef __SPARSE_MATRIX_PARTITION__
#define __SPARSE_MATRIX_PARTITION__

#include <PhysBAM_Tools/Math_Tools/INTERVAL.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
namespace PhysBAM{

class SPARSE_MATRIX_PARTITION
{
public:
    int number_of_sides;
    INTERVAL<int> interior_indices;
    ARRAY<INTERVAL<int> > ghost_indices;
    ARRAY<ARRAY<int> > boundary_indices;
    int interior_offset;
    ARRAY<int> neighbor_ranks;
    ARRAY<SPARSE_MATRIX_PARTITION*> neighbors;

    SPARSE_MATRIX_PARTITION()
        :number_of_sides(0)
    {}

    SPARSE_MATRIX_PARTITION(const int number_of_sides_input)
    {
        Set_Number_Of_Sides(number_of_sides_input);
    }

    void Set_Number_Of_Sides(const int number_of_sides_input)
    {number_of_sides=number_of_sides_input;
    ghost_indices.Resize(number_of_sides);boundary_indices.Resize(number_of_sides);
    neighbor_ranks.Resize(number_of_sides);neighbors.Resize(number_of_sides);}

    int Interior_Rows() const
    {return interior_indices.Size()+1;}

    template<class T>
    int Interior_Entries(const SPARSE_MATRIX_FLAT_NXN<T>& A) const
    {return A.offsets(interior_indices.max_corner+1)-A.offsets(interior_indices.min_corner);}

    void Set_Interior_Offset(const int previous_rows)
    {interior_offset=previous_rows-interior_indices.min_corner+1;}

    int Translate_Index(const int j) const
    {if(interior_indices.Lazy_Inside(j)) return Translate_Interior_Index(j);
    for(int r=1;r<number_of_sides;r++) if(ghost_indices(r).Lazy_Inside(j)) return Translate_Ghost_Index(j,r);
    return Translate_Ghost_Index(j,number_of_sides);}

    int Translate_Interior_Index(const int j) const
    {assert(interior_indices.Lazy_Inside(j));return j+interior_offset;}

    int Translate_Ghost_Index(const int j,const int region) const
    {assert(ghost_indices(region).Lazy_Inside(j));
    assert(ghost_indices(region).Size()+1==neighbors(region)->boundary_indices(((region-1)^1)+1).m);
    return neighbors(region)->boundary_indices(((region-1)^1)+1)(j-ghost_indices(region).min_corner+1)+neighbors(region)->interior_offset;}

//#####################################################################
};
}
#endif

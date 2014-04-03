//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OBJECTS_IN_CELL
//#####################################################################
#ifndef __OBJECTS_IN_CELL__    
#define __OBJECTS_IN_CELL__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_GRID.h>
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_1D.h>
#endif
namespace PhysBAM{

template<class T_GRID> struct GRID_ARRAYS_POLICY;

template <class T_GRID,class ID>
class OBJECTS_IN_CELL
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<ID>::TYPE T_ARRAYS_ID;typedef typename T_GRID::INDEX T_INDEX;typedef typename T_GRID::BLOCK T_BLOCK;
public:
    T_ARRAYS_ID object_in_cell;
    ARRAY<ARRAY<ID>,ID> object_list;
private:
    mutable ARRAY<ID> merge;
    mutable OPERATION_HASH<ID> operation_hash;
public:
    
    OBJECTS_IN_CELL();
    virtual ~OBJECTS_IN_CELL();

    void Reset(const T_GRID& grid,const int number_of_ghost_cells=3);
    void Add_Object_To_Cell(const T_INDEX& cell_index,const ID object_id);

    void Get_Objects_For_Cell(const T_INDEX& cell_index,ARRAY<ID>& objects) const;
    void Get_Objects_For_Cells(const T_INDEX& cell_index1,const T_INDEX& cell_index2,const ID number_of_collision_bodies,ARRAY<ID>& objects) const;
    void Get_Objects_For_Cells(const T_INDEX* cells,const int number_of_cells,const ID number_of_collision_bodies,ARRAY<ID>& objects) const;
    void Get_Objects_For_Cells(const T_BLOCK& block,const ID number_of_collision_bodies,ARRAY<ID>& objects) const;

private:
    void Get_Objects_For_Cells_Start(const ID number_of_collision_bodies) const;
    void Get_Objects_For_Cells_Cell(const T_INDEX& cell) const;
    void Get_Objects_For_Cells_End(ARRAY<ID>& objects) const
    {objects=merge;}

//#####################################################################
};
//template<class T> class OBJECTS_IN_CELL<RLE_GRID_1D<T> >{};
//template<class T> class OBJECTS_IN_CELL<BINTREE_GRID<T> >{};
}
#endif

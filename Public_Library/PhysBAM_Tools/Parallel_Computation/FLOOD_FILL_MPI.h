//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLOOD_FILL_MPI
//#####################################################################
#ifndef __FLOOD_FILL_MPI__
#define __FLOOD_FILL_MPI__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE_Arrays/GRID_ARRAYS_POLICY_RLE.h>
#endif
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_GRID.h>
namespace PhysBAM{

template<class T_GRID>
class FLOOD_FILL_MPI:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename MPI_GRID_POLICY<T_GRID>::MPI_GRID T_MPI_GRID;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::INDEX T_INDEX;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename REBIND<T_ARRAYS_SCALAR,int>::TYPE T_ARRAYS_INT;typedef typename REBIND<T_FACE_ARRAYS_SCALAR,bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename MPI_GRID_POLICY<T_GRID>::PARALLEL_GRID T_PARALLEL_GRID;
public:
    const T_MPI_GRID& mpi_grid;
    const T_GRID& local_grid;
    const T_FACE_ARRAYS_BOOL& psi_N;
    int& number_of_regions;
    T_ARRAYS_INT& colors;
    ARRAY<ARRAY<int> >& color_ranks;
    ARRAY<bool>* color_touches_uncolorable;

    FLOOD_FILL_MPI(const T_MPI_GRID& mpi_grid_input,const T_GRID& local_grid_input,const T_FACE_ARRAYS_BOOL& psi_N_input,int& number_of_regions_input,T_ARRAYS_INT& colors_input,
        ARRAY<ARRAY<int> >& color_ranks_input,ARRAY<bool>* color_touches_uncolorable);
    virtual ~FLOOD_FILL_MPI();

//#####################################################################
    int Synchronize_Colors();
    int Synchronize_Colors_Threaded();
protected:
    void Find_Global_Colors(ARRAY<bool,VECTOR<int,1> >& color_is_global,const RANGE<TV_INT>&) const;
    template<class T_BOX_HORIZONTAL_INT> void Find_Global_Colors(ARRAY<bool,VECTOR<int,1> >& color_is_global,const T_BOX_HORIZONTAL_INT&) const;
    struct Find_Global_Colors_Helper{template<class T_FACE> static void Apply(const T_MPI_GRID& mpi_grid,const T_GRID& local_grid,const ARRAY<bool>& psi_N,const ARRAY<int>& colors,
        ARRAY<bool,VECTOR<int,1> >& color_is_global);};
    void Translate_Local_Colors_To_Global_Colors(const ARRAY<int,VECTOR<int,1> >& color_map,T_ARRAYS_INT& colors_copy,const RANGE<TV_INT>& region,const int global_color_offset) const;
    template<class T_BOX_HORIZONTAL_INT> void Translate_Local_Colors_To_Global_Colors(const ARRAY<int,VECTOR<int,1> >& color_map,ARRAY<int>& colors_copy,const T_BOX_HORIZONTAL_INT& region,
        const int global_color_offset) const;
    void Find_Color_Matches(const ARRAY<int,VECTOR<int,1> >& color_map,UNION_FIND<>& union_find,T_ARRAYS_INT& colors_copy,const RANGE<TV_INT>& region,const int global_color_offset) const;
    template<class T_BOX_HORIZONTAL_INT> void Find_Color_Matches(const ARRAY<int,VECTOR<int,1> >& color_map,UNION_FIND<>& union_find,ARRAY<int>& colors_copy,const T_BOX_HORIZONTAL_INT& region,
        const int global_color_offset) const;
    void Remap_Colors(ARRAY<int,VECTOR<int,1> >& color_map,const RANGE<TV_INT>&);
    template<class T_BOX_HORIZONTAL_INT> void Remap_Colors(ARRAY<int,VECTOR<int,1> >& color_map,const T_BOX_HORIZONTAL_INT&);
//#####################################################################
};
}
#endif

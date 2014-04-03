//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_RLE_GRID
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_RLE_GRID__
#define __READ_WRITE_RLE_GRID__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
namespace PhysBAM{

template<class RW,class T,class POLICY>
class Read_Write<RLE_GRID<T,POLICY>,RW>
{
    typedef typename POLICY::VECTOR_T TV;typedef typename POLICY::VECTOR_INT TV_INT;typedef typename RLE_GRID_POLICY<TV>::RLE_GRID T_GRID;
public:
    static void Read(std::istream& input,RLE_GRID<T,POLICY>& object)
    {Read_Binary<RW>(input,object.uniform_grid,object.number_of_ghost_cells,object.minimum_vertical_space,object.minimum_long_run_length,object.long_run_cells,object.long_run_faces_horizontal);
    Read_Binary<RW>(input,object.negative_bandwidth,object.positive_bandwidth,object.columns,object.jmin,object.jmax);
    object.Topology_Changed();dynamic_cast<T_GRID&>(object).Compute_Auxiliary_Information();}

    static void Write(std::ostream& output,const RLE_GRID<T,POLICY>& object)
    {Write_Binary<RW>(output,object.uniform_grid,object.number_of_ghost_cells,object.minimum_vertical_space,object.minimum_long_run_length,object.long_run_cells,object.long_run_faces_horizontal);
    Write_Binary<RW>(output,object.negative_bandwidth,object.positive_bandwidth,object.columns,object.jmin,object.jmax);}
};
}

#endif
#endif
#endif

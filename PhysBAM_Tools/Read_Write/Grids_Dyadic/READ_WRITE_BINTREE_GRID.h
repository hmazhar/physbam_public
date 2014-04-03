//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_BINTREE_GRID
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_BINTREE_GRID__
#define __READ_WRITE_BINTREE_GRID__

#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_GRID.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Grids_Dyadic/READ_WRITE_BINTREE_CELL.h>
#include <PhysBAM_Tools/Read_Write/Grids_Dyadic/READ_WRITE_BINTREE_CHILDREN.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<BINTREE_GRID<T>,RW>
{
public:
    static void Read(std::istream& input,BINTREE_GRID<T>& object)
    {Read_Binary<RW>(input,object.uniform_grid,object.number_of_ghost_cells,object.maximum_depth);
    int number_of_cells=0,number_of_nodes=0,number_of_faces=0;Read_Binary<RW>(input,number_of_cells,number_of_nodes,number_of_faces);
    object.Initialize(object.uniform_grid,object.maximum_depth,object.number_of_ghost_cells,number_of_nodes!=0,number_of_faces!=0);
    for(int i=object.cells.domain.min_corner.x;i<=object.cells.domain.max_corner.x;i++) Read_Write<BINTREE_CELL<T>,RW>::Read(input,*object.cells(i));
    object.number_of_cells=number_of_cells;object.number_of_nodes=number_of_nodes;object.number_of_faces=number_of_faces;object.Tree_Topology_Changed();}

    static void Write(std::ostream& output,const BINTREE_GRID<T>& object)
    {Write_Binary<RW>(output,object.uniform_grid,object.number_of_ghost_cells,object.maximum_depth);
    Write_Binary<RW>(output,object.number_of_cells,object.number_of_nodes,object.number_of_faces);
    for(int i=object.cells.domain.min_corner.x;i<=object.cells.domain.max_corner.x;i++) Read_Write<BINTREE_CELL<T>,RW>::Write(output,*object.cells(i));}
};
}
#endif
#endif
#endif

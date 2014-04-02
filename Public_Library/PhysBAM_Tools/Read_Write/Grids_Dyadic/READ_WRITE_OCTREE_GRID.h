#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_OCTREE_GRID
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_OCTREE_GRID__
#define __READ_WRITE_OCTREE_GRID__

#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Dyadic/READ_WRITE_OCTREE_CHILDREN.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<OCTREE_GRID<T>,RW>
{
public:
    static void Read(std::istream& input,OCTREE_GRID<T>& object)
    {Read_Binary<RW>(input,object.uniform_grid,object.number_of_ghost_cells,object.maximum_depth);
    object.Initialize(object.uniform_grid,object.maximum_depth,object.number_of_ghost_cells,false,false);
    for(int i=object.cells.domain.min_corner.x;i<=object.cells.domain.max_corner.x;i++)for(int j=object.cells.domain.min_corner.y;j<=object.cells.domain.max_corner.y;j++)for(int ij=object.cells.domain.min_corner.z;ij<=object.cells.domain.max_corner.z;ij++)Read_Binary<RW>(input,object.cells(i,j,ij));
    Read_Binary<RW>(input,object.number_of_cells,object.number_of_nodes,object.number_of_faces);}

    static void Write(std::ostream& output,const OCTREE_GRID<T>& object)
    {Write_Binary<RW>(output,object.uniform_grid,object.number_of_ghost_cells,object.maximum_depth);
    for(int i=object.cells.domain.min_corner.x;i<=object.cells.domain.max_corner.x;i++)for(int j=object.cells.domain.min_corner.y;j<=object.cells.domain.max_corner.y;j++)for(int ij=object.cells.domain.min_corner.z;ij<=object.cells.domain.max_corner.z;ij++)Write_Binary<RW>(output,object.cells(i,j,ij));
    Write_Binary<RW>(output,object.number_of_cells,object.number_of_nodes,object.number_of_faces);}
};
}
#endif
#endif
#endif

//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_RED_GREEN_GRID_2D
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_RED_GREEN_GRID_2D__
#define __READ_WRITE_RED_GREEN_GRID_2D__

#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Red_Green/READ_WRITE_GREEN_CHILDREN_2D.h>
#include <PhysBAM_Geometry/Read_Write/Red_Green/READ_WRITE_RED_CHILDREN_2D.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_2D.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<RED_GREEN_GRID_2D<T>,RW>
{
public:
    static void Read(std::istream& input,RED_GREEN_GRID_2D<T>& object)
    {object.Clean_Memory();Read_Binary<RW>(input,object.number_of_cells,object.number_of_nodes,object.maximum_depth,object.uniform_grid);
    object.elements.Resize(object.uniform_grid.Domain_Indices());
    for(int i=1;i<=object.uniform_grid.counts.x;i++)for(int j=1;j<=object.uniform_grid.counts.y;j++){
        RED_CHILDREN_2D<T>* owner;Read_Binary<RW>(input,owner);
        if(owner){owner->parent=0;object.elements(i,j)=&owner->children[0];}}}

    static void Write(std::ostream& output,const RED_GREEN_GRID_2D<T>& object)
    {Write_Binary<RW>(output,object.number_of_cells,object.number_of_nodes,object.maximum_depth,object.uniform_grid);
    for(int i=1;i<=object.uniform_grid.counts.x;i++)for(int j=1;j<=object.uniform_grid.counts.y;j++)
        Write_Binary<RW>(output,object.elements(i,j)?object.elements(i,j)->owner:0);}
};
}
#endif
#endif
#endif

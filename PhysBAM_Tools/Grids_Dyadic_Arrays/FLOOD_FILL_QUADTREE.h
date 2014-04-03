//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLOOD_FILL_QUADTREE
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __FLOOD_FILL_QUADTREE__
#define __FLOOD_FILL_QUADTREE__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/STACK.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

template<class T> class QUADTREE_CELL;
template<class T> class QUADTREE_GRID;

template<class T>
class FLOOD_FILL_QUADTREE:public NONCOPYABLE
{
private:
    ARRAY<const QUADTREE_CELL<T>*> cell_index_to_pointer;
    STACK<const QUADTREE_CELL<T>*> flood_fill_stack;
public:

    FLOOD_FILL_QUADTREE()
    {}

//#####################################################################
    // colors should be initialized by the user with 0's where colors will be filled and -1 for uncolorable nodes (where no color will be filled)
    int Flood_Fill(const QUADTREE_GRID<T>& grid,ARRAY<int>& colors,const ARRAY<bool>& edge_is_blocked,ARRAY<bool>* color_touches_uncolorable_node=0); // returns number of colors
    bool Find_Uncolored_Node(const ARRAY<int>& colors,const QUADTREE_CELL<T>*& seed_cell);
    void Flood_Fill_From_Seed_Node(const QUADTREE_GRID<T>& grid,ARRAY<int>& colors,const int fill_color,const ARRAY<bool>& edge_is_blocked,bool& touches_uncolorable_node,const QUADTREE_CELL<T>* seed_node);
//#####################################################################
};
}
#endif
#endif

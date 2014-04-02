//#####################################################################
// Copyright 2005, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLOOD_FILL_GRAPH
//#####################################################################
#ifndef __FLOOD_FILL_GRAPH__
#define __FLOOD_FILL_GRAPH__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/STACK.h>
namespace PhysBAM{

class GRAPH;

class FLOOD_FILL_GRAPH
{
private:
    STACK<int> flood_fill_stack;
public:

    FLOOD_FILL_GRAPH()
    {}

//#####################################################################
    // colors should be initialized by the user with 0's where colors will be filled and -1 for uncolorable nodes (where no color will be filled)
    int Flood_Fill(const GRAPH& graph,ARRAY<int>& colors,ARRAY<bool>* color_touches_uncolorable_node=0); // returns number of colors
    bool Find_Uncolored_Node(const GRAPH& graph,const ARRAY<int>& colors,ARRAY<int>& uncolored_nodes,int& seed_node);
    void Flood_Fill_From_Seed_Node(const GRAPH& graph,ARRAY<int>& colors,const int fill_color,bool& touches_uncolorable_node,const int seed_node);
//#####################################################################
};
}
#endif

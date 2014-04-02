//#####################################################################
// Copyright 2005, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLOOD_FILL_GRAPH
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/FLOOD_FILL_GRAPH.h>
#include <PhysBAM_Tools/Data_Structures/GRAPH.h>
using namespace PhysBAM;
//#####################################################################
// Function Flood_Fill
//#####################################################################
// colors should be initialized by the user with 0's where colors will be filled and -1 for uncolorable nodes (where no color will be filled)
// returns number of colors
int FLOOD_FILL_GRAPH::
Flood_Fill(const GRAPH& graph,ARRAY<int>& colors,ARRAY<bool>* color_touches_uncolorable_node)
{
    int seed_node;int fill_color=0;
    flood_fill_stack.Preallocate(graph.edges.m);
    ARRAY<int> uncolored_nodes;uncolored_nodes.Preallocate(graph.valid_nodes.m);
    for(int i=1;i<=graph.valid_nodes.m;i++) if(graph.valid_nodes(i)) uncolored_nodes.Append(i);
    while(Find_Uncolored_Node(graph,colors,uncolored_nodes,seed_node)){
        bool touches_uncolorable_node;fill_color++;Flood_Fill_From_Seed_Node(graph,colors,fill_color,touches_uncolorable_node,seed_node);
        if(color_touches_uncolorable_node)color_touches_uncolorable_node->Append(touches_uncolorable_node);}
    flood_fill_stack.Clean_Memory();
    return fill_color;
}
//#####################################################################
// Function Find_Uncolored_Node
//#####################################################################
bool FLOOD_FILL_GRAPH::
Find_Uncolored_Node(const GRAPH& graph,const ARRAY<int>& colors,ARRAY<int>& uncolored_nodes,int& seed_node)
{
    while(uncolored_nodes.m>0){
        int i=uncolored_nodes(uncolored_nodes.m);
        if(graph.valid_nodes(i)&&colors(i)==0){seed_node=i;return true;}
        else uncolored_nodes.Remove_Index_Lazy(uncolored_nodes.m);}
    return false;
}
//#####################################################################
// Function Flood_Fill_From_Seed_Node
//#####################################################################
void FLOOD_FILL_GRAPH::
Flood_Fill_From_Seed_Node(const GRAPH& graph,ARRAY<int>& colors,const int fill_color,bool& touches_uncolorable_node,const int seed_node)
{
    assert(colors(seed_node)==0);touches_uncolorable_node=false;
    flood_fill_stack.Remove_All();
    flood_fill_stack.Push(seed_node);
    while(!flood_fill_stack.Empty()){
        const int node=flood_fill_stack.Pop();
        if(colors(node)==-1){touches_uncolorable_node=true;continue;}else if(colors(node)!=0)continue;colors(node)=fill_color;
        for(int i=1;i<=graph.edges(node).m;i++){
            int neighbor_node=graph.edges(node)(i);
            if(graph.valid_nodes(neighbor_node)&&colors(neighbor_node)<=0) flood_fill_stack.Push(neighbor_node);}}
}
//#####################################################################

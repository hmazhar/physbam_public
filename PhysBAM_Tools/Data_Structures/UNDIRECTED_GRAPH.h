//#####################################################################
// Copyright 2004-2008, Ronald Fedkiw, Geoffrey Irving, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNDIRECTED_GRAPH
//#####################################################################
#ifndef __UNDIRECTED_GRAPH__
#define __UNDIRECTED_GRAPH__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
namespace PhysBAM{

class DIRECTED_GRAPH_CORE;
class UNDIRECTED_GRAPH_CORE
{
public:
    ARRAY<PAIR<int,int> > edges; // edges(e) returns (parent,child) of edge e
    ARRAY<ARRAY<int> > adjacent_edges; // adjacent_edges(j) returns list of edges connected to node j

    explicit UNDIRECTED_GRAPH_CORE(const int number_nodes=0,const int number_edges=0);

//##################################################################### 
    void Reset();
    void Add_Edge(const int parent,const int child,const int edge);
    void Remove_Edge(const int edge);
    void Modify_Edge(const int parent,const int child,const int edge);
    void Breadth_First_Directed_Graph(const int root_node,DIRECTED_GRAPH_CORE& directed_graph);
    int Connected_Components(ARRAY<int>& component_id); // returns total number of components
    void Depth_First_Directed_Graph(const int root_node,DIRECTED_GRAPH_CORE& directed_graph);
private:
    void Generate_Next_Level_Of_Breadth_First_Directed_Graph(DIRECTED_GRAPH_CORE& directed_graph,ARRAY<bool>& marked_nodes,ARRAY<bool,int>& marked_edges,QUEUE<int>& queue);
    void Generate_Next_Level_Of_Depth_First_Directed_Graph(DIRECTED_GRAPH_CORE& directed_graph,ARRAY<bool>& marked,const int node);
//#####################################################################
};

template<class ID,class EID>
class UNDIRECTED_GRAPH
{
public:
    UNDIRECTED_GRAPH_CORE core;

    explicit UNDIRECTED_GRAPH(const ID number_nodes=ID(),const EID number_edges=EID())
        :core(Value(number_nodes),Value(number_edges))
    {}

    ID Last_Node() const
    {return ID(core.adjacent_edges.m);}

    EID Last_Edge() const
    {return EID(core.edges.m);}

    const PAIR<ID,ID> Edges(const EID e) const // edges(e) returns (parent,child) of edge e
    {return PAIR<ID,ID>(core.edges(Value(e)));}

    void Ensure_Number_Nodes(const ID n)
    {if(n>core.adjacent_edges.m) core.adjacent_edges.Resize(n);}

    const ARRAY<EID>& Adjacent_Edges(const ID j) const // adjacent_edges(j) returns list of edges connected to node j
    {return (const ARRAY<EID>&)core.adjacent_edges(j);}

    void Reset()
    {core.Reset();}

    void Add_Edge(const ID parent,const ID child,const EID edge)
    {core.Add_Edge(Value(parent),Value(child),Value(edge));}

    void Remove_Edge(const EID edge)
    {core.Remove_Edge(Value(edge));}

    void Modify_Edge(const ID parent,const ID child,const EID edge)
    {core.Modify_Edge(Value(parent),Value(child),Value(edge));}

    void Breadth_First_Directed_Graph(const ID root_node,DIRECTED_GRAPH<ID>& directed_graph)
    {core.Breadth_First_Directed_Graph(Value(root_node),directed_graph.core);}

    int Connected_Components(ARRAY<int,ID>& component_id) // returns total number of components
    {return core.Connected_Components((ARRAY<int>&)component_id);}

    void Depth_First_Directed_Graph(const ID root_node,DIRECTED_GRAPH<ID>& directed_graph)
    {core.Depth_First_Directed_Graph(Value(root_node,directed_graph.core));}

//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2002-2009, Kevin Der, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Bridget Vuong.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIRECTED_GRAPH
//#####################################################################
#ifndef __DIRECTED_GRAPH__
#define __DIRECTED_GRAPH__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/CONSTANT_ARRAY.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
namespace PhysBAM{

class DIRECTED_GRAPH_CORE
{
public:
    ARRAY<ARRAY<int> > parents,children;
    ARRAY<int> level_of_node; // all cycles are condensed into one level
    ARRAY<ARRAY<int> > nodes_in_level;

    DIRECTED_GRAPH_CORE(int number_nodes);

//#####################################################################
    void Reset();
    void Initialize(const int number_nodes);
    void Resize(const int number_nodes);
    void Add_Edge(const int parent,const int child,bool ensure_distinct=false);
    void Remove_Edge(const int parent,const int child);
    bool Topological_Sort(ARRAY<int>& finish_time,ARRAY<int>& node_index); // returns false if cycles are found
    void Topological_Sort_Assuming_Cycle_Free(ARRAY<int>& finish_time,ARRAY<int>& node_index);
    int Strongly_Connected_Components(ARRAY<int>& component_id); // returns total number of components
    void Generate_Levels(); // collapses all cycles into a single level
    void Maximal_Depth_On_Acyclic_Graph(ARRAY<int>& depths) const; // depth is 1-based
private:
    void Visit(const int index,ARRAY<int>& finish_time,ARRAY<int>& node_index,ARRAY<short >& visit_tag,int& time,bool& cycle_free);
    void Visit_Transpose(const int index,ARRAY<short>& visit_tag,ARRAY<int>& component);
    void Visit_Assuming_Cycle_Free(const int index,ARRAY<int>& finish_time,ARRAY<int>& node_index,int& time);
//#####################################################################
};

template<class ID>
class DIRECTED_GRAPH
{
public:
    DIRECTED_GRAPH_CORE core;
    template<class PID,class EID> friend class UNDIRECTED_GRAPH;

    DIRECTED_GRAPH(const ID number_nodes)
        :core(Value(number_nodes))
    {}

    int Number_Of_Levels() const
    {return core.nodes_in_level.m;}

    ID Number_Nodes() const
    {return core.parents.m;}

    const ARRAY<ID>& Nodes_In_Level(const int level) const
    {return (const ARRAY<ID>&)core.nodes_in_level(level);}

    ARRAY<ID>& Nodes_In_Level(const int level)
    {return (ARRAY<ID>&)core.nodes_in_level(level);}

    int Level_Of_Node(const ID level) const
    {return core.level_of_node(Value(level));}

    const ARRAY<ID>& Parents(const ID e) const
    {return (const ARRAY<ID>&)core.parents(Value(e));}

    const ARRAY<ID>& Children(const ID e) const
    {return (const ARRAY<ID>&)core.children(Value(e));}

    ARRAY<ID>& Parents(const ID e)
    {return (ARRAY<ID>&)core.parents(Value(e));}

    ARRAY<ID>& Children(const ID e)
    {return (ARRAY<ID>&)core.children(Value(e));}

    void Reset()
    {core.Reset();}

    void Initialize(const ID number_nodes)
    {core.Initialize(Value(number_nodes));}

    void Resize(const ID number_nodes)
    {core.Resize(Value(number_nodes));}

    void Add_Edge(const ID parent,const ID child,bool ensure_distinct=false)
    {core.Add_Edge(Value(parent),Value(child),ensure_distinct);}

    void Remove_Edge(const ID parent,const ID child)
    {core.Remove_Edge(Value(parent),Value(child));}

    bool Topological_Sort(ARRAY<ID,ID>& finish_time,ARRAY<ID,ID>& node_index) // returns false if cycles are found
    {return core.Topological_Sort((ARRAY<int>&)finish_time,(ARRAY<int>&)node_index);}

    void Topological_Sort_Assuming_Cycle_Free(ARRAY<ID,ID>& finish_time,ARRAY<ID,ID>& node_index)
    {core.Topological_Sort_Assuming_Cycle_Free((ARRAY<int>&)finish_time,(ARRAY<int>&)node_index);}

    int Strongly_Connected_Components(ARRAY<int,ID>& component_id) // returns total number of components
    {return core.Strongly_Connected_Components((ARRAY<int>&)component_id);}

    void Generate_Levels() // collapses all cycles into a single level
    {core.Generate_Levels();}

    void Maximal_Depth_On_Acyclic_Graph(ARRAY<int,ID>& depths) const // depth is 1-based
    {core.Maximal_Depth_On_Acyclic_Graph((ARRAY<int>&)depths);}
//#####################################################################
};
}
#endif

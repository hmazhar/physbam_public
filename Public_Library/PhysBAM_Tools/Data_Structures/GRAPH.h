//#####################################################################
// Copyright 2005, Ron Fedkiw, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRAPH
//#####################################################################
#ifndef __GRAPH__
#define __GRAPH__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
namespace PhysBAM{

class GRAPH
{
public:
    ARRAY<ARRAY<int> > edges;
    ARRAY<bool> valid_nodes;

    GRAPH(const int number_of_nodes)
    {
        edges.Resize(number_of_nodes);
        valid_nodes.Resize(number_of_nodes,false);ARRAYS_COMPUTATIONS::Fill(valid_nodes,true);
    }

    void Add_Directed_Edge(const int from,const int to)
    {edges(from).Append(to);}

    void Add_Undirected_Edge(const int node1,const int node2) // not very space efficient, but i'm not sure of a better way right now
    {Add_Directed_Edge(node1,node2);Add_Directed_Edge(node2,node1);}

    bool Edge_Present(const int from,const int to) const // not efficient for more than sparsely connected graphs (need to use sorted arrays and binary search in that case)
    {for(int i=1;i<=edges(from).m;i++)if(edges(from)(i)==to)return true;return false;}

    void Remove_Directed_Edge(const int from,const int to)
    {for(int i=1;i<=edges(from).m;i++)if(edges(from)(i)==to)edges(from).Remove_Index(i);}

//#####################################################################
};
}
#endif


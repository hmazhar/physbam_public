//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_UNDIRECTED_GRAPH
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_UNDIRECTED_GRAPH__
#define __READ_WRITE_UNDIRECTED_GRAPH__

#include <PhysBAM_Tools/Data_Structures/UNDIRECTED_GRAPH.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW>
class Read_Write<UNDIRECTED_GRAPH_CORE,RW>
{
public:
    static void Read(std::istream& input,UNDIRECTED_GRAPH_CORE& object)
    {ARRAY<int,int> parent_node_index,child_node_index;
    Read_Binary<RW>(input,parent_node_index,child_node_index,object.adjacent_edges);
    if(parent_node_index.Size()!=child_node_index.Size()) PHYSBAM_FATAL_ERROR();
    object.edges.Resize(parent_node_index.Size(),false,false);
    for(int e(1);e<=object.edges.m;e++){object.edges(e).x=parent_node_index(e);object.edges(e).y=child_node_index(e);}}

    static void Write(std::ostream& output,const UNDIRECTED_GRAPH_CORE& object)
    {ARRAY<int,int> parent_node_index(object.edges.m,false),child_node_index(object.edges.m,false);
    for(int e(1);e<=object.edges.m;e++){parent_node_index(e)=object.edges(e).x;child_node_index(e)=object.edges(e).y;}
    Write_Binary<RW>(output,parent_node_index,child_node_index,object.adjacent_edges);}
};

template<class RW,class ID,class EID>
class Read_Write<UNDIRECTED_GRAPH<ID,EID>,RW>
{
public:
    static void Read(std::istream& input,UNDIRECTED_GRAPH<ID,EID>& object)
    {Read_Binary<RW>(input,object.core);}

    static void Write(std::ostream& output,const UNDIRECTED_GRAPH<ID,EID>& object)
    {Write_Binary<RW>(output,object.core);}
};
}
#endif
#endif

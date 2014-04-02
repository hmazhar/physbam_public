//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_DIRECTED_GRAPH
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_DIRECTED_GRAPH__
#define __READ_WRITE_DIRECTED_GRAPH__

#include <PhysBAM_Tools/Data_Structures/DIRECTED_GRAPH.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW>
class Read_Write<DIRECTED_GRAPH_CORE,RW>
{
public:
    static void Read(std::istream& input,DIRECTED_GRAPH_CORE& object)
    {Read_Binary<RW>(input,object.parents,object.children,object.level_of_node,object.nodes_in_level);}
        
    static void Write(std::ostream& output,const DIRECTED_GRAPH_CORE& object)
    {Write_Binary<RW>(output,object.parents,object.children,object.level_of_node,object.nodes_in_level);}
};

template<class RW,class ID>
class Read_Write<DIRECTED_GRAPH<ID>,RW>
{
public:
    static void Read(std::istream& input,DIRECTED_GRAPH<ID>& object)
    {Read_Binary<RW>(input,object.core);}

    static void Write(std::ostream& output,const DIRECTED_GRAPH<ID>& object)
    {Write_Binary<RW>(output,object.core);}
};
}
#endif
#endif

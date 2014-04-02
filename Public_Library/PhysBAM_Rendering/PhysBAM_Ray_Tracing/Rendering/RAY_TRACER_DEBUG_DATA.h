//#####################################################################
// Copyright 2004, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RAY_TRACER_DEBUG_DATA
//#####################################################################
#ifndef __RAY_TRACER_DEBUG_DATA__
#define __RAY_TRACER_DEBUG_DATA__

#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDERING_RAY_DEBUG.h>
namespace PhysBAM{

template<class T>
class RAY_TRACER_DEBUG_DATA:public NONCOPYABLE
{
public:
    RENDERING_RAY_DEBUG<T>* ray_tree;
    ARRAY<std::string> messages;

    RAY_TRACER_DEBUG_DATA()
        :ray_tree(0)
    {}

    ~RAY_TRACER_DEBUG_DATA()
    {delete ray_tree;}

    void Add_Message(const std::string& message)
    {messages.Append(message);}

    void Set_Ray_Tree(RENDERING_RAY_DEBUG<T>* ray_tree_in)
    {if(ray_tree) Add_Message("RAY_TREE being set even though it's already set");
    ray_tree=ray_tree_in;}

//#####################################################################
};
}
#endif

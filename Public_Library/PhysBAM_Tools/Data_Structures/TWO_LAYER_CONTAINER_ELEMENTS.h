//#####################################################################
// Copyright 2007-2008, Michael Lentine, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TWO_LAYER_CONTAINER_ELEMENTS
//#####################################################################
#ifndef __TWO_LAYER_CONTAINER_ELEMENTS__
#define __TWO_LAYER_CONTAINER_ELEMENTS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/CONTAINER_ELEMENTS.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Log/LOG.h>

namespace PhysBAM{
template<class T_CONTAINER,class T_CONTAINER2>
class TWO_LAYER_CONTAINER_ELEMENTS:public CONTAINER_ELEMENTS<T_CONTAINER>
{
    typedef CONTAINER_ELEMENTS<T_CONTAINER> BASE;
    typedef typename BASE::ITERATOR BASE_ITERATOR;
public:
    using BASE::elements_of_container;using BASE::Update_Containers;
    typedef typename T_CONTAINER::ID ID1;
    typedef typename T_CONTAINER2::ID ID2;
    HASHTABLE<ID2,ARRAY<ID1> > containers_of_top_level_container;

    TWO_LAYER_CONTAINER_ELEMENTS()
    {}

    ~TWO_LAYER_CONTAINER_ELEMENTS()
    {}

    void Update_Top_Level_Containers(const ARRAY<T_CONTAINER2*,ID2>& top_level_containers,const ARRAY<PAIR<ID2,ID2> >& swap_pairs,const ARRAY<ID2>& rebuild,const ARRAY<ID2>& remove)
    {
        for(int i=1;i<=swap_pairs.m;i++) containers_of_top_level_container.Exchange(swap_pairs(i).x,swap_pairs(i).y);
        for(int i=1;i<=remove.m;i++) containers_of_top_level_container.Delete_If_Present(remove(i));
        for(int i=1;i<=rebuild.m;i++){ID2 id=rebuild(i);
            ARRAY<ID1>& list=containers_of_top_level_container.Get_Or_Insert(id);list.Remove_All();
            for(int f=1;f<=top_level_containers(id)->containers.m;f++) if(elements_of_container.Contains(top_level_containers(id)->containers(f))) list.Append(top_level_containers(id)->containers(f));}
    }

    void Print() const
    {LOG::cout<<"containers_of_top_level_container"<<std::endl<<containers_of_top_level_container<<std::endl;BASE::Print();}

    class ITERATOR:public BASE_ITERATOR
    {
    public:
        using BASE_ITERATOR::container_elements;using BASE_ITERATOR::containers;using BASE_ITERATOR::current_container_index;
        using BASE_ITERATOR::current_list;using BASE_ITERATOR::current_list_index;using BASE_ITERATOR::Next_Container_If_Exists;
        using BASE_ITERATOR::Next;using BASE_ITERATOR::Data;using BASE_ITERATOR::Valid;using BASE_ITERATOR::Pointer_To_Array;

    protected:
        ITERATOR(const TWO_LAYER_CONTAINER_ELEMENTS<T_CONTAINER,T_CONTAINER2>& container_elements_input,ARRAY_VIEW<const ID1> containers_input)
            :BASE_ITERATOR(container_elements_input,containers_input)
        {}
        
    public:
        ITERATOR(const TWO_LAYER_CONTAINER_ELEMENTS<T_CONTAINER,T_CONTAINER2>& container_elements_input,const ID2 id)
            :BASE_ITERATOR(container_elements_input,Pointer_To_Array(container_elements.containers_of_top_level_container.Get_Pointer(id)))
        {Next();}
    
        ITERATOR(const TWO_LAYER_CONTAINER_ELEMENTS<T_CONTAINER,T_CONTAINER2>& container_elements_input,int no_id) // Don't create constructor with only one input, since that would result in misuse
            :BASE_ITERATOR(container_elements_input,container_elements_input.containers)
        {assert(!no_id);Next();}
    
        virtual ~ITERATOR() {};
    };
};
}
#endif

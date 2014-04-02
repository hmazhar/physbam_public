//#####################################################################
// Copyright 2007-2009, Michael Lentine, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONTAINER_ELEMENTS
//#####################################################################
#ifndef __CONTAINER_ELEMENTS__
#define __CONTAINER_ELEMENTS__

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Log/LOG.h>

namespace PhysBAM{
template<class T_CONTAINER>
class CONTAINER_ELEMENTS
{
public:
    typedef typename T_CONTAINER::ID ID;
    HASHTABLE<ID,ARRAY<int> > elements_of_container;
    ARRAY<ID> containers;

    CONTAINER_ELEMENTS()
    {}

    ~CONTAINER_ELEMENTS()
    {}

    template<class T_ARRAY>
    void Update_Containers(const ARRAY_BASE<int,T_ARRAY>& base_elements,const ARRAY<ID>& element_to_container_id,const bool remove_duplicates=false)
    {
        const T_ARRAY& elements=base_elements.Derived();
        elements_of_container.Remove_All();
        for(int t=1;t<=elements.Size();t++) if(ID id=element_to_container_id(elements(t))) elements_of_container.Get_Or_Insert(id).Append(elements(t));
        if(remove_duplicates)
            for(HASHTABLE_ITERATOR<ID,ARRAY<int> > iterator(elements_of_container);iterator.Valid();iterator.Next())
                iterator.Data().Prune_Duplicates();
        elements_of_container.Get_Keys(containers);
    }

    void Print() const
    {LOG::cout<<"elements_of_container"<<std::endl<<elements_of_container<<std::endl<<std::endl;}

    class ITERATOR
    {
    public:
        const CONTAINER_ELEMENTS<T_CONTAINER>& container_elements;
        ARRAY_VIEW<const ID> containers;
        int current_container_index,current_list_index;
        ARRAY_VIEW<const int> current_list;
        
    protected:
        ITERATOR(const CONTAINER_ELEMENTS<T_CONTAINER>& container_elements_input,ARRAY_VIEW<const ID> containers_input)
            :container_elements(container_elements_input),containers(containers_input),current_container_index(0),current_list_index(0),current_list(0,0)
        {}

    public:
        ITERATOR(const CONTAINER_ELEMENTS<T_CONTAINER>& container_elements_input)
            :container_elements(container_elements_input),containers(container_elements.containers),current_container_index(0),current_list_index(0),current_list(0,0)
        {Next();}

        virtual ~ITERATOR() {};

    protected:
        template<class S>
        static ARRAY_VIEW<const S> Pointer_To_Array(const ARRAY<S>* array)
        {return array?*array:ARRAY_VIEW<const S>(0,0);}

    public:
        void Next() PHYSBAM_ALWAYS_INLINE
        {
            if(current_list_index<current_list.Size()){current_list_index++;return;}
            Next_Container(); // separate this into a separate function call so that Next itself will get inlined
        }

        bool Valid() const
        {return current_list_index<=current_list.Size();} // it's important for performance that this condition is the same as the first line of Next

        int Data() const
        {return current_list(current_list_index);}

        bool Next_Container_If_Exists()
        {
            current_list_index=1;
            while(current_container_index<containers.Size()){
                current_container_index++;
                ARRAY_VIEW<const int> next_list(Pointer_To_Array(container_elements.elements_of_container.Get_Pointer(containers(current_container_index))));
                current_list.Const_Cast().Exchange(next_list.Const_Cast());
                if(current_list.Size()) return true;}
            return false;
        }

        //TODO: Make this non-virtual
        virtual void Next_Container()
        {if(!Next_Container_If_Exists()){ARRAY_VIEW<int> empty(0,0);current_list.Const_Cast().Exchange(empty);}}
    };
};
}
#endif

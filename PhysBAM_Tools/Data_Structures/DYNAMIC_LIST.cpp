//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DYNAMIC_LIST
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/DYNAMIC_LIST.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
DYNAMIC_LIST_CORE::
DYNAMIC_LIST_CORE(void (*deleter)(void*))
    :deleter(deleter),last_unique_id(0)
{}
//#####################################################################
// Destructor
//#####################################################################
DYNAMIC_LIST_CORE::
~DYNAMIC_LIST_CORE()
{
    Delete_All();
}
//#####################################################################
// Function Delete_And_Clear
//#####################################################################
void DYNAMIC_LIST_CORE::
Delete_And_Clear(void* pointer)
{
    deleter(pointer);pointer=0;
}
//#####################################################################
// Function Delete_All
//#####################################################################
void DYNAMIC_LIST_CORE::
Delete_All()
{
    for(int i=1;i<=array.m;i++) deleter(array(i));
    array.Clean_Memory();
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
void DYNAMIC_LIST_CORE::
Clean_Memory()
{
    Delete_All();
    pointer_to_id_map.Clean_Memory();
    needs_write.Clean_Memory();index_to_id_map.Clean_Memory();id_to_index_map.Clean_Memory();last_unique_id=0;deletion_list.Clean_Memory();
}
//#####################################################################
// Function Remove_All
//#####################################################################
void DYNAMIC_LIST_CORE::
Remove_All()
{
    Delete_All();
    pointer_to_id_map.Clean_Memory();
    index_to_id_map.Remove_All();id_to_index_map.Remove_All();
    needs_write.Remove_All();deletion_list.Remove_All();last_unique_id=0;
}
//#####################################################################
// Function Add_Element
//#####################################################################
int DYNAMIC_LIST_CORE::
Add_Element(void* element)
{
    int id;
    if(pointer_to_id_map.Get(element,id)){assert(array(id_to_index_map(id))==element);return id;}
    array.Append(element);
    if(deletion_list.m){id=deletion_list.Pop();assert(id_to_index_map(id)==0);id_to_index_map(id)=array.Size();}
    else{id=++last_unique_id;id_to_index_map.Append(array.Size());assert(id_to_index_map.Size()==id);}
    index_to_id_map.Append(id);needs_write.Append(id);
    pointer_to_id_map.Set(element,id);
    return id;
}
//#####################################################################
// Function Reactivate_Element
//#####################################################################
int DYNAMIC_LIST_CORE::
Reactivate_Element(void* element,const int id_number)
{
    assert(id_to_index_map(id_number)==0);
    array.Append(element);needs_write.Append(id_number);
    pointer_to_id_map.Set(element,id_number);
    index_to_id_map.Append(id_number);id_to_index_map(id_number)=array.Size();PHYSBAM_ASSERT(index_to_id_map.Size()==array.Size());
    return id_number;
}
//#####################################################################
// Function Deactivate_Element
//#####################################################################
void DYNAMIC_LIST_CORE::
Deactivate_Element(const int id,const bool delete_element)
{
    int index=id_to_index_map(id);assert(index);
    pointer_to_id_map.Delete(array(index));
    if(delete_element) Delete_And_Clear(array(index));
    id_to_index_map(id)=0;array.Remove_Index_Lazy(index);index_to_id_map.Remove_Index_Lazy(index);
    if(index<=array.Size()) id_to_index_map(index_to_id_map(index))=index;
}
//#####################################################################
// Function Swap_Elements
//#####################################################################
int DYNAMIC_LIST_CORE::
Swap_Elements(void* element,const int id_number,const int id)
{
    int index=id_to_index_map(id);assert(index);
    pointer_to_id_map.Delete(array(index));
    id_to_index_map(id)=0;
    assert(id_to_index_map(id_number)==0);
    array(index)=element;needs_write.Append(id_number);
    pointer_to_id_map.Set(element,id_number);
    index_to_id_map(index)=id_number;id_to_index_map(id_number)=index;assert(index_to_id_map.Size()==array.Size());
    return id_number;
}
//#####################################################################
// Function Remove_Elements
//#####################################################################
void DYNAMIC_LIST_CORE::
Remove_Element(const int id,const bool delete_element=true,const bool allow_id_reuse=true)
{
    Deactivate_Element(id,delete_element);if(allow_id_reuse)deletion_list.Append(id);
}
//#####################################################################
// Function Remove_Elements
//#####################################################################
void DYNAMIC_LIST_CORE::
Purge_Element(const int id)
{
    assert(id_to_index_map(id)==0);
    deletion_list.Append(id);
}
//#####################################################################
// Function Fill_Needs_Write
//#####################################################################
void DYNAMIC_LIST_CORE::
Fill_Needs_Write()
{
    needs_write=index_to_id_map;
}
//#####################################################################
}

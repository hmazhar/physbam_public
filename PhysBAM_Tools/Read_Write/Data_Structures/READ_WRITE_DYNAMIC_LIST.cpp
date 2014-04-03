//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/DYNAMIC_LIST.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_DYNAMIC_LIST.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class RW> void Read_Write<DYNAMIC_LIST_CORE,RW>::
Read(DYNAMIC_LIST_CORE& object,const std::string& prefix,ARRAY<int>& needs_init)
{
    object.pointer_to_id_map.Clean_Memory();
    needs_init.Remove_All();
    object.needs_write.Remove_All();
    ARRAY<int> active_ids;char version;
    FILE_UTILITIES::Read_From_File<RW>(STRING_UTILITIES::string_sprintf("%sactive_ids",prefix.c_str()),version,object.last_unique_id,active_ids);
    PHYSBAM_ASSERT(version==1);
    ARRAY<void*> new_array;
    ARRAY<bool> element_copied(object.array.Size());
    object.id_to_index_map.Resize(object.last_unique_id);
    for(int i=1;i<=active_ids.Size();i++){
        int index=object.id_to_index_map(active_ids(i));
        if(index){new_array.Append(object.array(index));element_copied(index)=true;}
        else{new_array.Append((void*)(0));needs_init.Append(active_ids(i));}}
    for(int i=1;i<=object.array.Size();i++)if(!element_copied(i)) object.Delete_And_Clear(object.array(i));
    object.index_to_id_map.Resize(active_ids.Size());
    ARRAYS_COMPUTATIONS::Fill(object.id_to_index_map,0);
    for(int i=1;i<=new_array.Size();i++){
        object.pointer_to_id_map.Set(new_array(i),active_ids(i));
        object.index_to_id_map(i)=active_ids(i);
        object.id_to_index_map(active_ids(i))=i;}
    object.array.Exchange(new_array);
}
//#####################################################################
// Constructor
//#####################################################################
template<class RW> void Read_Write<DYNAMIC_LIST_CORE,RW>::
Write(const DYNAMIC_LIST_CORE& object,const std::string& prefix)
{
    const char version=1;
    FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%sactive_ids",prefix.c_str()),version,object.last_unique_id,object.index_to_id_map);
    for(int i=object.needs_write.m;i>=1;i--) if(!object.id_to_index_map(object.needs_write(i))) object.needs_write.Remove_Index_Lazy(i);
    // handle case of new element which was removed without being written
}
template class Read_Write<DYNAMIC_LIST_CORE,float,void>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class Read_Write<DYNAMIC_LIST_CORE,double,void>;
#endif

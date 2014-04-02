//#####################################################################
// Copyright 2009, Michael Lentine, Avi Robinson-Mosher
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_ARRAY_COLLECTION
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_ELEMENT_ID.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_HASHTABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{
//#####################################################################
// Function Read_Write_Array_Collection_Registry
//#####################################################################
template<class RW> HASHTABLE<ATTRIBUTE_ID,READ_WRITE_ARRAY_COLLECTION_FUNCTIONS>& Read_Write<ARRAY_COLLECTION,RW>::
Read_Write_Array_Collection_Registry()
{
    static HASHTABLE<ATTRIBUTE_ID,READ_WRITE_ARRAY_COLLECTION_FUNCTIONS> read_write_array_collection_registry;
    return read_write_array_collection_registry;
}
//#####################################################################
// Function Attribute_Names_Registry
//#####################################################################
static HASHTABLE<ATTRIBUTE_ID,const char*>& Attribute_Names_Registry()
{
    static HASHTABLE<ATTRIBUTE_ID,const char*> names_registry;
    return names_registry;
}
//#####################################################################
// Function Register_Attribute_Name
//#####################################################################
void Register_Attribute_Name(const ATTRIBUTE_ID id,const char* name)
{
    PHYSBAM_ASSERT(Attribute_Names_Registry().Set(id,name));
}
//#####################################################################
// Function Register_Attribute_Name
//#####################################################################
const char* Get_Attribute_Name(const ATTRIBUTE_ID id)
{
    if(const char** name=Attribute_Names_Registry().Get_Pointer(id)) return *name;
    return 0;
}
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
//#####################################################################
// Function Read_Arrays
//#####################################################################
template<class RW> void Read_Write<ARRAY_COLLECTION,RW>::
Read_Arrays(std::istream& input,ARRAY_COLLECTION& object)
{
    int size;
    ATTRIBUTE_INDEX num_attributes;
    Read_Binary<RW>(input,size,num_attributes);
    object.Resize(size);

    for(ATTRIBUTE_INDEX i(1);i<=num_attributes;i++){
        ATTRIBUTE_ID hashed_id;int read_size;
        Read_Binary<RW>(input,hashed_id,read_size);

        READ_WRITE_ARRAY_COLLECTION_FUNCTIONS* read_write_functions=Read_Write_Array_Collection_Registry().Get_Pointer(Type_Only(hashed_id));
        if(!read_write_functions){input.ignore(read_size);continue;}
        if(!read_write_functions->Read && read_write_functions->sample_attribute)
            PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("No read registered for id %i\n",Value(hashed_id)));

        ATTRIBUTE_INDEX index=object.Get_Attribute_Index(Id_Only(hashed_id));
        if(!index) index=object.Add_Array(Id_Only(hashed_id),read_write_functions->sample_attribute->Clone_Default());
        // TODO: this really ought to know whether we're running in float or double
        else read_write_functions=Read_Write_Array_Collection_Registry().Get_Pointer(Type_Only(object.arrays(index)->Hashed_Id()));
        read_write_functions->Read(input,*object.arrays(index));
        object.arrays(index)->id=Id_Only(hashed_id);}
}
//#####################################################################
// Function Write_Arrays
//#####################################################################
template<class RW> void Read_Write<ARRAY_COLLECTION,RW>::
Write_Arrays(std::ostream& output,const ARRAY_COLLECTION& object)
{
    Write_Binary<RW>(output,object.number,object.arrays.m);
    for(ATTRIBUTE_INDEX i(1);i<=object.arrays.m;i++){
        const ARRAY_COLLECTION_ELEMENT_BASE* entry=object.arrays(i);
        const READ_WRITE_ARRAY_COLLECTION_FUNCTIONS* read_write_functions=Read_Write_Array_Collection_Registry().Get_Pointer(Type_Only(entry->Hashed_Id()));
        if(!read_write_functions || !read_write_functions->Write || !read_write_functions->Write_Size)
            PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("No write registered for id %i (type %s)\n",Value(entry->id),typeid(*entry).name()));

        int calculated_write_size=read_write_functions->Write_Size(*entry);
        Write_Binary<RW>(output,entry->Typed_Hashed_Id(RW()),calculated_write_size);
        if(calculated_write_size) read_write_functions->Write(output,*entry);}
}
//#####################################################################
// Function Print
//#####################################################################
template<class RW> void Read_Write<ARRAY_COLLECTION,RW>::
Print(std::ostream& output,const ARRAY_COLLECTION& object,const int p)
{
    if(p<1 || p>object.number) throw INDEX_ERROR("Index out of range");
    for(ATTRIBUTE_INDEX i(1);i<=object.arrays.m;i++){
        const ARRAY_COLLECTION_ELEMENT_BASE* entry=object.arrays(i);
        const READ_WRITE_ARRAY_COLLECTION_FUNCTIONS* read_write_functions=Read_Write_Array_Collection_Registry().Get_Pointer(Type_Only(entry->Hashed_Id()));
        if(!read_write_functions || !read_write_functions->Print)
            PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("No print registered for id %i (type %s)\n",Value(entry->id),typeid(*entry).name()));
        read_write_functions->Print(output,*entry,p);}
}
//#####################################################################
template class Read_Write<ARRAY_COLLECTION,float>;
template class Read_Write<ARRAY_COLLECTION,double>;
}
#endif

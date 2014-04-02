//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_COLLECTION
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/EXTERNAL_ARRAY_COLLECTION.h>
#include <sstream>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
EXTERNAL_ARRAY_COLLECTION::
EXTERNAL_ARRAY_COLLECTION()
    :ARRAY_COLLECTION()
{
    delete_data=false;
}
//#####################################################################
// Destructor
//#####################################################################
EXTERNAL_ARRAY_COLLECTION::
~EXTERNAL_ARRAY_COLLECTION()
{
    for(ATTRIBUTE_INDEX i(1);i<=arrays.m;i++) if(Owns_Element_From_Index(i)) delete arrays(i);
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
void EXTERNAL_ARRAY_COLLECTION::
Clean_Memory()
{
    number=buffer_size=0;
    for(ATTRIBUTE_INDEX i(1);i<=arrays.m;i++) if(Owns_Element_From_Index(i)) arrays(i)->Clean_Memory();
}
//#####################################################################
// Function Add_Array
//#####################################################################
ATTRIBUTE_INDEX EXTERNAL_ARRAY_COLLECTION::
Add_Array(const ATTRIBUTE_ID attribute_id,ARRAY_COLLECTION_ELEMENT_BASE* array,bool owned_element)
{
    ATTRIBUTE_INDEX index=Find_Attribute_Index(attribute_id);
    if(index<=arrays.m && arrays(index)->id==attribute_id){
        PHYSBAM_ASSERT(array==arrays(index));
        return index;}
    array->id=attribute_id;
    if(owned_element){array->Reallocate(buffer_size);array->Set_Size(number);}
    else Set_External_From_Index(index);
    arrays.Insert(array,index);
    return index;
}
//#####################################################################
// Function Reallocate_Buffer
//#####################################################################
void EXTERNAL_ARRAY_COLLECTION::
Reallocate_Buffer(int new_size)
{
    buffer_size=new_size;
    for(ATTRIBUTE_INDEX i(1);i<=arrays.m;i++) if(Owns_Element_From_Index(i)) arrays(i)->Reallocate(buffer_size);
}
//#####################################################################
// Function Resize
//#####################################################################
void EXTERNAL_ARRAY_COLLECTION::
Resize(const int new_size)
{
    if(buffer_size<new_size) Reallocate_Buffer(max(4*number/3+2,new_size));
    number=new_size;
    for(ATTRIBUTE_INDEX i(1);i<=arrays.m;i++) if(Owns_Element_From_Index(i)) arrays(i)->Set_Size(number);
}
//#####################################################################
}

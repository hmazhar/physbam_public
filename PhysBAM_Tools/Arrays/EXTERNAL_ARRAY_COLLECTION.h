//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXTERNAL_ARRAY_COLLECTION
//#####################################################################
#ifndef __EXTERNAL_ARRAY_COLLECTION__
#define __EXTERNAL_ARRAY_COLLECTION__

#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION.h>
#include <memory>
namespace PhysBAM{

class EXTERNAL_ARRAY_COLLECTION:public ARRAY_COLLECTION
{
    typedef ARRAY_COLLECTION BASE;
public:
    using BASE::Find_Attribute_Index;using BASE::buffer_size;using BASE::number;using BASE::delete_data;
    
    HASHTABLE<ATTRIBUTE_INDEX> external_elements;

    EXTERNAL_ARRAY_COLLECTION();
    ~EXTERNAL_ARRAY_COLLECTION();

    void Remove_Array_Using_Index(const ATTRIBUTE_INDEX attribute_index)
    {if(Owns_Element_From_Index(attribute_index)) ARRAY_COLLECTION::Remove_Array_Using_Index(attribute_index);}

    bool External_Element(const ATTRIBUTE_ID attribute_id)
    {return External_Element_From_Index(Find_Attribute_Index(attribute_id));}

    bool External_Element_From_Index(const ATTRIBUTE_INDEX attribute_index)
    {return external_elements.Size()>0 && external_elements.Contains(attribute_index);}

    bool Owns_Element(const ATTRIBUTE_ID attribute_id)
    {return !External_Element(attribute_id);}
    
    bool Owns_Element_From_Index(const ATTRIBUTE_INDEX attribute_index)
    {return !External_Element_From_Index(attribute_index);}

    void Set_External(const ATTRIBUTE_ID attribute_id)
    {Set_External_From_Index(Find_Attribute_Index(attribute_id));}
    
    void Set_External_From_Index(const ATTRIBUTE_INDEX attribute_index)
    {if(Owns_Element_From_Index(attribute_index)) external_elements.Insert(attribute_index);}

//#####################################################################
    void Clean_Memory();
    ATTRIBUTE_INDEX Add_Array(const ATTRIBUTE_ID attribute_id,ARRAY_COLLECTION_ELEMENT_BASE* array,bool owned_element);
protected:
    void Reallocate_Buffer(int new_size);
    void Resize(const int new_size);
//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2009, Michael Lentine, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_ARRAY_COLLECTION
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_ARRAY_COLLECTION__
#define __READ_WRITE_ARRAY_COLLECTION__

#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Arrays/ATTRIBUTE_ID.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY_COLLECTION_ELEMENT.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_HASHTABLE.h>
namespace PhysBAM{

struct READ_WRITE_ARRAY_COLLECTION_FUNCTIONS_HELPER
{
    ARRAY<ARRAY_COLLECTION_ELEMENT_BASE*> samples;
    ~READ_WRITE_ARRAY_COLLECTION_FUNCTIONS_HELPER()
    {samples.Delete_Pointers_And_Clean_Memory();}
};

struct READ_WRITE_ARRAY_COLLECTION_FUNCTIONS
{
    void (*Read)(std::istream&,ARRAY_COLLECTION_ELEMENT_BASE&);
    void (*Write)(std::ostream&,const ARRAY_COLLECTION_ELEMENT_BASE&);
    int (*Write_Size)(const ARRAY_COLLECTION_ELEMENT_BASE&);
    void (*Print)(std::ostream&,const ARRAY_COLLECTION_ELEMENT_BASE&,const int);
    std::string name;
    ARRAY_COLLECTION_ELEMENT_BASE* sample_attribute;
};

void Register_Attribute_Name(const ATTRIBUTE_ID id,const char* name);
const char* Get_Attribute_Name(const ATTRIBUTE_ID id);

template<class RW>
class Read_Write<ARRAY_COLLECTION,RW>
{
    static HASHTABLE<ATTRIBUTE_ID,READ_WRITE_ARRAY_COLLECTION_FUNCTIONS>& Read_Write_Array_Collection_Registry();
public:

    static void Read(std::istream& input,ARRAY_COLLECTION& object)
    {int size;
    Read_Binary<RW>(input,size);
    if(size<0) throw READ_ERROR(STRING_UTILITIES::string_sprintf("Invalid negative size %d",size));
    object.Clean_Memory();
    object.Resize(size);
    Read_Arrays(input,object);}

    static void Write(std::ostream& output,const ARRAY_COLLECTION& object)
    {Write_Binary<RW>(output,object.Size());Write_Arrays(output,object);}

    static ATTRIBUTE_ID Type_Only(ATTRIBUTE_ID id)
    {return ATTRIBUTE_ID(Value(id)&0xFFFF0000);}

    static ATTRIBUTE_ID Id_Only(ATTRIBUTE_ID id)
    {return ATTRIBUTE_ID(Value(id)&0x0000FFFF);}

//#####################################################################
    template<class T> static void Register_Read_Write();
    static void Read_Arrays(std::istream& input,ARRAY_COLLECTION& object);
    static void Write_Arrays(std::ostream& output,const ARRAY_COLLECTION& object);
    static void Print(std::ostream& output,const ARRAY_COLLECTION& object,const int p);
//#####################################################################
};
//#####################################################################
// Function Register_Read_Write
//#####################################################################
template<class RW> template<class T> void Read_Write<ARRAY_COLLECTION,RW>::
Register_Read_Write()
{
    static READ_WRITE_ARRAY_COLLECTION_FUNCTIONS_HELPER sample_helper;
    READ_WRITE_ARRAY_COLLECTION_FUNCTIONS functions;
    functions.Read=Read_Write<ARRAY_COLLECTION_ELEMENT<T>,RW>::Read;
    functions.Write=Read_Write<ARRAY_COLLECTION_ELEMENT<T>,RW>::Write;
    functions.Write_Size=Read_Write<ARRAY_COLLECTION_ELEMENT<T>,RW>::Write_Size;
    functions.Print=Read_Write<ARRAY_COLLECTION_ELEMENT<T>,RW>::Print;
    functions.sample_attribute=new ARRAY_COLLECTION_ELEMENT<T>();
    functions.sample_attribute->id=ATTRIBUTE_ID();
    PHYSBAM_ASSERT(Read_Write_Array_Collection_Registry().Set(Type_Only(functions.sample_attribute->Hashed_Id()),functions));
    sample_helper.samples.Append(functions.sample_attribute);
}
}
#endif
#endif

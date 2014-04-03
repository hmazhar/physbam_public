//#####################################################################
// Copyright 2008-2009, Geoffrey Irving, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_COLLECTION
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <sstream>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
ARRAY_COLLECTION::
ARRAY_COLLECTION()
    :number(0),buffer_size(0),delete_data(true)
{}
//#####################################################################
// Constructor
//#####################################################################
ARRAY_COLLECTION::
~ARRAY_COLLECTION()
{
    if(delete_data) for(ATTRIBUTE_INDEX i(1);i<=arrays.m;i++) delete arrays(i);
}
//#####################################################################
// Function Initialize
//#####################################################################
void ARRAY_COLLECTION::
Initialize(const ARRAY_COLLECTION& elements)
{
    Clean_Memory();
    Add_Arrays(elements);
    Append(elements);
}
//#####################################################################
// Function Initialize
//#####################################################################
void ARRAY_COLLECTION::
Initialize(ARRAY_VIEW<const ARRAY_COLLECTION* const> elements_per_cell)
{
    Clean_Memory();
    int total_number=0;for(int c=1;c<=elements_per_cell.Size();c++) if(elements_per_cell(c)){
        total_number+=elements_per_cell(c)->number;
        Add_Arrays(*elements_per_cell(c));} // include arrays that occur on any of the cell elements
    Preallocate(total_number);
    for(int c=1;c<=elements_per_cell.Size();c++) if(elements_per_cell(c)) Append(*elements_per_cell(c));
}
//#####################################################################
// Function Add_Arrays
//#####################################################################
void ARRAY_COLLECTION::
Add_Arrays(const ARRAY_COLLECTION& collection)
{
    ATTRIBUTE_INDEX i(1),j(1);
    for(;i<=arrays.m && j<=collection.arrays.m;i++){
        if(arrays(i)->id<collection.arrays(j)->id) continue;
        if(arrays(i)->id>collection.arrays(j)->id) Add_Array(collection.arrays(j)->id,collection.arrays(j)->Clone_Default());
        j++;}
    for(;j<=collection.arrays.m;j++)
        Add_Array(collection.arrays(j)->id,collection.arrays(j)->Clone_Default());
}
//#####################################################################
// Function Add_Elements_From_Deletion_List
//#####################################################################
void ARRAY_COLLECTION::
Add_Elements_From_Deletion_List(const int count,ARRAY<int>& added_indices)
{
    added_indices.Preallocate(added_indices.Size()+count);
    int added=min(deletion_list.m,count);
    added_indices.Append_Elements(deletion_list.Pop_Elements(added));
    added_indices.Append_Elements(Add_Elements(count-added));
}
//#####################################################################
// Function Delete_Elements_On_Deletion_List
//#####################################################################
void ARRAY_COLLECTION::
Delete_Elements_On_Deletion_List(const bool preserve_order)
{
    Sort(deletion_list);
    if(preserve_order){
        for(int k=1;k<=deletion_list.m;k++){
            int next=k<deletion_list.m?deletion_list(k+1):number+1;
            for(int i=deletion_list(k)+1;i<next;i++) Copy_Element_Helper(i,i-k);}}
    else{
        int last=number;
        for(int k=deletion_list.m;k>=1;k--)
            Copy_Element_Helper(last--,deletion_list(k));}
    Resize(number-deletion_list.m);
    deletion_list.Remove_All();
}
//#####################################################################
// Function Copy_Element_Helper
//#####################################################################
void ARRAY_COLLECTION::
Copy_Element(const ARRAY_COLLECTION& from_collection,const int from,const int to)
{
    ATTRIBUTE_INDEX i(1),j(1);
    while(i<=arrays.m && j<=from_collection.arrays.m){
        if(arrays(i)->id<from_collection.arrays(j)->id) arrays(i++)->Clear(to);
        else if(arrays(i)->id>from_collection.arrays(j)->id) j++;
        else arrays(i++)->Copy_Element(*from_collection.arrays(j++),from,to);}
    for(;i<=arrays.m;i++) arrays(i)->Clear(to);
}
//#####################################################################
// Function Copy_All_Elements_Helper
//#####################################################################
void ARRAY_COLLECTION::
Copy_All_Elements_Helper(const ARRAY_COLLECTION& from_collection,const int offset)
{
    ATTRIBUTE_INDEX i(1),j(1);
    while(i<=arrays.m && j<=from_collection.arrays.m){
        if(arrays(i)->id<from_collection.arrays(j)->id) arrays(i++)->Clear_Range(offset+1,offset+from_collection.number);
        else if(arrays(i)->id>from_collection.arrays(j)->id) j++;
        else arrays(i++)->Copy_With_Offset(*from_collection.arrays(j++),offset);}
    for(;i<=arrays.m;i++) arrays(i)->Clear_Range(offset+1,offset+from_collection.number);
}
//#####################################################################
// Function Get_Attribute_Index
//#####################################################################
ATTRIBUTE_INDEX ARRAY_COLLECTION::
Find_Attribute_Index(const ATTRIBUTE_ID attribute_id) const
{
    ATTRIBUTE_INDEX first(1),last(arrays.m+1);
    while(first<last){
        ATTRIBUTE_INDEX middle((Value(first)+Value(last))/2);
        if(arrays(middle)->id<attribute_id) first=middle+1;
        else last=middle;}
    return first;
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
void ARRAY_COLLECTION::
Clean_Memory()
{
    number=buffer_size=0;
    for(ATTRIBUTE_INDEX i(1);i<=arrays.m;i++)
        arrays(i)->Clean_Memory();
}
//#####################################################################
// operator==
//#####################################################################
bool ARRAY_COLLECTION::
operator==(const ARRAY_COLLECTION& collection) const
{
    if(this==&collection) return true;
    if(arrays.m!=collection.arrays.m) return false;
    for(ATTRIBUTE_INDEX i(1);i<=arrays.m;i++){
        if(arrays(i)->id!=collection.arrays(i)->id) return false;
        assert(typeid(arrays(i))==typeid(collection.arrays(i)));
        if(arrays(i)!=collection.arrays(i)) return false;}
    return true;
}
//#####################################################################
// Function Resize
//#####################################################################
void ARRAY_COLLECTION::
Resize(const int new_size)
{
    if(buffer_size<new_size) Reallocate_Buffer(max(4*number/3+2,new_size));
    number=new_size;
    for(ATTRIBUTE_INDEX i(1);i<=arrays.m;i++) arrays(i)->Set_Size(number);
}
//#####################################################################
// Function Add_Array
//#####################################################################
ATTRIBUTE_INDEX ARRAY_COLLECTION::
Add_Array(const ATTRIBUTE_ID attribute_id,ARRAY_COLLECTION_ELEMENT_BASE* array)
{
    ATTRIBUTE_INDEX index=Find_Attribute_Index(attribute_id);
    if(index<=arrays.m && arrays(index)->id==attribute_id){
        PHYSBAM_ASSERT(array==arrays(index));
        return index;}
    array->id=attribute_id;
    array->Reallocate(buffer_size);
    array->Set_Size(number);
    arrays.Insert(array,index);
    return index;
}
//#####################################################################
// Function Pack_Size
//#####################################################################
int ARRAY_COLLECTION::
Pack_Size() const
{
    int pack_size=0;
    for(ATTRIBUTE_INDEX i(1);i<=arrays.m;i++)
        pack_size+=arrays(i)->Pack_Size();
    return pack_size;
}
//#####################################################################
// Function Pack
//#####################################################################
void ARRAY_COLLECTION::
Pack(ARRAY_VIEW<char> buffer,int& position,const int p) const
{
    assert(1<=p && p<=number);
    for(ATTRIBUTE_INDEX i(1);i<=arrays.m;i++)
        arrays(i)->Pack(buffer,position,p);
}
//#####################################################################
// Function Unpack
//#####################################################################
void ARRAY_COLLECTION::
Unpack(ARRAY_VIEW<const char> buffer,int& position,const int p)
{
    assert(1<=p && p<=number);
    for(ATTRIBUTE_INDEX i(1);i<=arrays.m;i++)
        arrays(i)->Unpack(buffer,position,p);
}
//#####################################################################
// Function Copy_Element_Helper
//#####################################################################
void ARRAY_COLLECTION::
Copy_Element_Helper(const int from,const int to)
{
    for(ATTRIBUTE_INDEX i(1);i<=arrays.m;i++)
        arrays(i)->Copy_Element(from,to);
}
//#####################################################################
// Function Reallocate_Buffer
//#####################################################################
void ARRAY_COLLECTION::
Reallocate_Buffer(int new_size)
{
    buffer_size=new_size;
    for(ATTRIBUTE_INDEX i(1);i<=arrays.m;i++)
        arrays(i)->Reallocate(buffer_size);
}
//#####################################################################
}

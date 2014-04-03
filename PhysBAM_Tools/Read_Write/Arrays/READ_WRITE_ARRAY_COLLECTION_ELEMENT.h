//#####################################################################
// Copyright 2009, Michael Lentine, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_ARRAY_COLLECTION_ELEMENT
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_ARRAY_COLLECTION_ELEMENT__
#define __READ_WRITE_ARRAY_COLLECTION_ELEMENT__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION_ELEMENT.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY_VIEW.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_FRAME.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_TWIST.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<ARRAY_COLLECTION_ELEMENT<T>,RW>
{
    struct UNUSABLE{};
public:

    static void Read_Helper(typename IF<IS_POINTER<T>::value,std::istream&,UNUSABLE>::TYPE input,ARRAY_COLLECTION_ELEMENT<T>& object)
    {}

    static void Read_Helper(typename IF<IS_POINTER<T>::value,UNUSABLE,std::istream&>::TYPE input,ARRAY_COLLECTION_ELEMENT<T>& object)
    {Read_Binary<RW>(input,*object.array);}

    static void Write_Helper(typename IF<IS_POINTER<T>::value,std::ostream&,UNUSABLE>::TYPE output,const ARRAY_COLLECTION_ELEMENT<T>& object)
    {}

    static void Write_Helper(typename IF<IS_POINTER<T>::value,UNUSABLE,std::ostream&>::TYPE output,const ARRAY_COLLECTION_ELEMENT<T>& object)
    {Write_Binary<RW>(output,*object.array);}

    static void Read(std::istream& input,ARRAY_COLLECTION_ELEMENT_BASE& object)
    {Read_Helper(input,dynamic_cast<ARRAY_COLLECTION_ELEMENT<T>&>(object));}

    static void Write(std::ostream& output,const ARRAY_COLLECTION_ELEMENT_BASE& object)
    {Write_Helper(output,dynamic_cast<const ARRAY_COLLECTION_ELEMENT<T>&>(object));}

    static int Write_Size(typename IF<IS_POINTER<T>::value,const ARRAY_COLLECTION_ELEMENT_BASE&,UNUSABLE>::TYPE object)
    {return 0;}

    static int Write_Size(typename IF<IS_POINTER<T>::value,UNUSABLE,const ARRAY_COLLECTION_ELEMENT_BASE&>::TYPE object)
    {return sizeof(T)*dynamic_cast<const ARRAY_COLLECTION_ELEMENT<T>&>(object).array->Size()+sizeof(int);}

    static void Print_Helper(std::ostream& output,const ARRAY_COLLECTION_ELEMENT<T>& object,const int p)
    {if(const char* name=Get_Attribute_Name(object.id)) output<<name;
    else output<<"id "<<object.id;
    output<<" = "<<(*object.array)(p)<<std::endl;}

    static void Print(std::ostream& output,const ARRAY_COLLECTION_ELEMENT_BASE& object,const int p)
    {Print_Helper(output,dynamic_cast<const ARRAY_COLLECTION_ELEMENT<T>&>(object),p);}
};
//#####################################################################
}
#endif
#endif

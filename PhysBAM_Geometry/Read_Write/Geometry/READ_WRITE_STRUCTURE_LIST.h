//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_STRUCTURE_LIST
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_STRUCTURE_LIST__
#define __READ_WRITE_STRUCTURE_LIST__

#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_DYNAMIC_LIST.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_STRUCTURE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
#include <string>
namespace PhysBAM{

template<class RW,class TV,class ID>
class Read_Write<STRUCTURE_LIST<TV,ID>,RW>
{
    typedef DYNAMIC_LIST<STRUCTURE<TV>,ID> OBJECT_BASE;
public:
    static void Read(const std::string& prefix,const std::string& list_name,const int frame,STRUCTURE_LIST<TV,ID>& object)
    {ARRAY<ID> needs_init;
    Read_Write<OBJECT_BASE,RW>::Read(object,STRING_UTILITIES::string_sprintf("%s/%d/%s",prefix.c_str(),frame,list_name.c_str()),needs_init);
    FILE_UTILITIES::Read_From_File<RW>(STRING_UTILITIES::string_sprintf("%s/common/%skey",prefix.c_str(),list_name.c_str()),object.names);
    for(int i=1;i<=needs_init.m;i++){ID id=needs_init(i),index=object.Element_Index(id);
        object.Set_Active_Element(index,STRUCTURE<TV>::Create_From_Name(object.names(id)));
        FILE_UTILITIES::Read_From_File<RW>(STRING_UTILITIES::string_sprintf("%s/common/%s%d.%s",prefix.c_str(),list_name.c_str(),id,object.Active_Element(index)->Extension().c_str()),*object.Active_Element(index));}

    for(ID id=1;id<=object.Size();id++)
        if(object.Is_Active(id)){
            STRUCTURE<TV>& structure=*object.Active_Element(object.Element_Index(id));
            std::string filename=STRING_UTILITIES::string_sprintf("%s/%d/%s%d.%s",prefix.c_str(),frame,list_name.c_str(),id,structure.Extension().c_str());
            if(FILE_UTILITIES::File_Exists(filename)) {
                FILE_UTILITIES::Read_From_File<RW>(filename,structure);
                structure.update_every_frame=true;
            } else structure.update_every_frame=false;
        }
    }

    static void Write(const std::string& prefix, const std::string& list_name,const int frame,const STRUCTURE_LIST<TV,ID>& object)
    {Read_Write<OBJECT_BASE,RW>::Write(object,STRING_UTILITIES::string_sprintf("%s/%d/%s",prefix.c_str(),frame,list_name.c_str()));
    object.names.Resize(Value(object.Size()));
    ARRAY<ID>& needs_write=object.Needs_Write();
    for(int i=1;i<=needs_write.Size();i++){ID id=needs_write(i);int index=object.Element_Index(id);assert(index);
        object.names(id)=object.Active_Element(index)->Name();
        FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%s/common/%s%d.%s",prefix.c_str(),list_name.c_str(),id,object.Active_Element(index)->Extension().c_str()),*object.Active_Element(index));}
    if(frame==0 || needs_write.Size()) FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%s/common/%skey",prefix.c_str(),list_name.c_str()),object.names);
    needs_write.Remove_All();

    for(ID id=1;id<=object.Size();id++)
        if(object.Is_Active(id)){
            STRUCTURE<TV>& structure=*object.Active_Element(object.Element_Index(id));
            if(structure.update_every_frame){
                std::string filename=STRING_UTILITIES::string_sprintf("%s/%d/%s%d.%s",prefix.c_str(),frame,list_name.c_str(),id,structure.Extension().c_str());
                FILE_UTILITIES::Write_To_File<RW>(filename,structure);}}
    }
};
}

#endif
#endif

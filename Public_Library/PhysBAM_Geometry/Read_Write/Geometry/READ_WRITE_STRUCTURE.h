//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_STRUCTURE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_STRUCTURE__
#define __READ_WRITE_STRUCTURE__

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<STRUCTURE<TV>,RW>
{
    static HASHTABLE<std::string,void (*)(std::istream&,STRUCTURE<TV>&)>& Read_Registry();
    static HASHTABLE<std::string,void (*)(std::istream&,STRUCTURE<TV>&)>& Read_Structure_Registry();
    static HASHTABLE<std::string,void (*)(std::ostream&,const STRUCTURE<TV>&)>& Write_Registry();
    static HASHTABLE<std::string,void (*)(std::ostream&,const STRUCTURE<TV>&)>& Write_Structure_Registry();

public:
    static void Read(std::istream& input,STRUCTURE<TV>& object)
    {std::string name=std::string(typeid(object).name());
    (*Read_Registry().Get_Pointer(name))(input,object);}

    static void Read_Structure(std::istream& input,STRUCTURE<TV>& object)
    {std::string name;Read_Binary<RW>(input,name);
    if(name!=object.Name()){LOG::cerr<<"Trying to read in a "<<name<<" as a "<<object.Name()<<std::endl;PHYSBAM_FATAL_ERROR();}
    (*Read_Structure_Registry().Get_Pointer(std::string(typeid(object).name())))(input,object);}

    static void Read_Structure_Without_Name(std::istream& input,STRUCTURE<TV>& object)
    {(*Read_Structure_Registry().Get_Pointer(std::string(typeid(object).name())))(input,object);}

    static void Read_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
    {PHYSBAM_NOT_IMPLEMENTED();}

    static void Read_Structure_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
    {PHYSBAM_NOT_IMPLEMENTED();}

    static void Write(std::ostream& output,const STRUCTURE<TV>& object)
    {(*Write_Registry().Get_Pointer(std::string(typeid(object).name())))(output,object);}

    static void Write_Structure(std::ostream& output,const STRUCTURE<TV>& object)
    {Write_Binary<RW>(output,object.Name());
    (*Write_Structure_Registry().Get_Pointer(std::string(typeid(object).name())))(output,object);}

    static void Write_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
    {PHYSBAM_NOT_IMPLEMENTED();}

    static void Write_Structure_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
    {PHYSBAM_NOT_IMPLEMENTED();}
    
    static STRUCTURE<TV>* Create_Structure(std::istream& input,GEOMETRY_PARTICLES<TV>& particles)
    {std::string name;Read_Binary<RW>(input,name);
    STRUCTURE<TV>* structure=STRUCTURE<TV>::Create_From_Name(name,particles);
    Read_Structure_Without_Name(input,*structure);return structure;}

    static STRUCTURE<TV>* Create_From_File(const std::string& filename)
    {STRUCTURE<TV>* structure=STRUCTURE<TV>::Create_From_Extension(FILE_UTILITIES::Get_File_Extension(filename));
    FILE_UTILITIES::Read_From_File<RW>(filename,*structure);return structure;}

    template<class T_STRUCTURE> static void Register_Read_Write();
};
//#####################################################################
// Function Register_Read_Write
//#####################################################################
template<class RW,class TV> template<class T_STRUCTURE> void Read_Write<STRUCTURE<TV>,RW>::
Register_Read_Write()
{
    // TODO: string has custom comparison function, this won't work.  Need to hash it or something?
    std::string name=typeid(T_STRUCTURE).name();
    Read_Registry().Insert(name,Read_Write<T_STRUCTURE,RW>::Read_Helper);
    Read_Structure_Registry().Insert(name,Read_Write<T_STRUCTURE,RW>::Read_Structure_Helper);
    Write_Registry().Insert(name,Read_Write<T_STRUCTURE,RW>::Write_Helper);
    Write_Structure_Registry().Insert(name,Read_Write<T_STRUCTURE,RW>::Write_Structure_Helper);
}
}
#endif
#endif

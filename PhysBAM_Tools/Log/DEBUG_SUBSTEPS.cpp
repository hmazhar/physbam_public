//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace DEBUG_SUBSTEPS
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
namespace PhysBAM{
//#####################################################################
int DEBUG_SUBSTEPS::write_substeps_level=-1;
void* DEBUG_SUBSTEPS::writer_object=0;
void (*DEBUG_SUBSTEPS::writer)(void*,const std::string&,int,int)=0;
//#####################################################################
// Function Write_Substep
//#####################################################################
void DEBUG_SUBSTEPS::
Write_Substep_Helper(const std::string& title,const int substep,const int level)
{
    if(writer) (*writer)(writer_object,title,substep,level);
}
//#####################################################################
// Function Set_Substep_Writer
//#####################################################################
void DEBUG_SUBSTEPS::
Set_Substep_Writer(void* writer_object_input,void (*writer_input)(void*,const std::string&,int,int))
{
    if(!writer_object) {writer_object=writer_object_input;writer=writer_input;} // Only one driver instantiated gets to write substeps...
}
//#####################################################################
// Function Set_Substep_Writer
//#####################################################################
void DEBUG_SUBSTEPS::
Clear_Substep_Writer(void* writer_object_input)
{
    if(writer_object==writer_object_input){writer_object=0;writer=0;} // Only clear writer if caller registered it
}
//#####################################################################
// Function Write_Substeps_Level
//#####################################################################
int DEBUG_SUBSTEPS::
Write_Substeps_Level()
{
    return write_substeps_level;
}
//#####################################################################
// Function Set_Write_Substeps_Level
//#####################################################################
void DEBUG_SUBSTEPS::
Set_Write_Substeps_Level(const int level)
{
    write_substeps_level=level;
}
//#####################################################################
}

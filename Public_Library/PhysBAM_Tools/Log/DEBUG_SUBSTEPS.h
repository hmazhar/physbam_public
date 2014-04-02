//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEBUG_SUBSTEPS
//#####################################################################
#ifndef __DEBUG_SUBSTEPS__
#define __DEBUG_SUBSTEPS__

#include <string>
namespace PhysBAM{
class DEBUG_SUBSTEPS
{
    static int write_substeps_level;
    static void* writer_object;
    static void (*writer)(void *,const std::string&,int,int);
public:

//#####################################################################
    static int Write_Substeps_Level();
    static void Set_Write_Substeps_Level(const int level);
    static void Set_Substep_Writer(void* writer_object,void writer(void *,const std::string&,int,int));
    static void Clear_Substep_Writer(void* writer_object);
    static void Write_Substep_Helper(const std::string& title,const int substep,const int level=0); // don't call directly
//#####################################################################
};

#define PHYSBAM_DEBUG_WRITE_SUBSTEP(title,substep,level) \
    do \
        if(level<=PhysBAM::DEBUG_SUBSTEPS::Write_Substeps_Level()) \
            PhysBAM::DEBUG_SUBSTEPS::Write_Substep_Helper(title,substep,level); \
    while(0)

}
#endif

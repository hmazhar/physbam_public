//#####################################################################
// Copyright 2004-2007, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LOG_SCOPE
//##################################################################### 
#ifndef __LOG_SCOPE__
#define __LOG_SCOPE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Log/LOG_ENTRY.h>
#include <PhysBAM_Tools/Utilities/TIMER.h>
#include <string>
namespace PhysBAM{
namespace LOG_REAL{

class LOG_SCOPE:public LOG_ENTRY
{
public:
    HASHTABLE<std::string,int> entries;
    ARRAY<LOG_ENTRY*> children;
    std::string scope_identifier;

    LOG_SCOPE(LOG_ENTRY* parent_input,int depth_input,int timer_id_input,const std::string& scope_identifier_input,const std::string& name_input,int& verbosity_level_input)
        :LOG_ENTRY(parent_input,depth_input,timer_id_input,name_input,verbosity_level_input),scope_identifier(scope_identifier_input)
    {}

    virtual ~LOG_SCOPE()
    {children.Delete_Pointers_And_Clean_Memory();}

    virtual LOG_ENTRY* Get_Stop_Time(LOG_CLASS& instance)
    {return this;}

    LOG_ENTRY* Get_New_Scope(LOG_CLASS& instance,const std::string& new_scope_identifier,const std::string& new_name)
    {int entry;end_on_separate_line=true;log_file_end_on_separate_line=true;
    if(entries.Get(new_scope_identifier,entry)){children(entry)->name=new_name;return children(entry);}
    LOG_ENTRY* new_entry=new LOG_SCOPE(this,depth+1,timer_id,new_scope_identifier,new_name,verbosity_level);
    children.Append(new_entry);
    entries.Insert(new_scope_identifier,children.m);
    return new_entry;}

    LOG_ENTRY* Get_New_Item(LOG_CLASS& instance,const std::string& new_name)
    {int entry;end_on_separate_line=true;log_file_end_on_separate_line=true;
    if(entries.Get(new_name,entry)) return children(entry);
    LOG_ENTRY* new_entry=new LOG_ENTRY(this,depth+1,timer_id,new_name,verbosity_level);
    children.Append(new_entry);
    entries.Insert(new_name,children.m);
    return new_entry;}

    LOG_ENTRY* Get_Pop_Scope(LOG_CLASS& instance)
    {Stop(instance);return parent;}

    void Start_XML(LOG_CLASS& instance) PHYSBAM_OVERRIDE
    {fprintf(instance.log_file,"%*s<scope id=\"%s\" name=\"%s\">",2*depth,"",scope_identifier.c_str(),name.c_str());}

    void Dump_Log(FILE* output) PHYSBAM_OVERRIDE
    {fprintf(output,"%*s%-*s%8.4f\n",2*depth,"",50-2*depth,scope_identifier.c_str(),time);fflush(output);
    for(int i=1;i<=children.m;i++) children(i)->Dump_Log(output);}

    void Dump_Names(FILE* output) PHYSBAM_OVERRIDE
    {LOG_ENTRY::Dump_Names(output);for(int i=1;i<=children.m;i++)children(i)->Dump_Names(output);}

//##################################################################### 
};
}
}
#endif

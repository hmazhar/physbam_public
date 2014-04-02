//#####################################################################
// Copyright 2004. Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARAMETER_LIST_CONTAINER
//##################################################################### 
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __PARAMETER_LIST_CONTAINER__
#define __PARAMETER_LIST_CONTAINER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <string>
namespace PhysBAM{

class PARAMETER_LIST_CONTAINER
{
private:
    ARRAY<std::string> lines;
    ARRAY<std::string>::iterator line;
    std::string partial_line;
public:

    PARAMETER_LIST_CONTAINER(const std::string& filename)
    {
        Pre_Process_File(filename);
        Reset_Parser();
    }
    
    void Reset_Parser()
    {line=lines.begin();partial_line=*line;}

    bool Get_Next_Parameter_List(std::string& identifier,PARAMETER_LIST& list);

//##################################################################### 
private:
    void Pre_Process_File(const std::string& filename);
//##################################################################### 
};
}
#endif
#endif

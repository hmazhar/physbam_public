//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GENERIC_PARSER
//##################################################################### 
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __GENERIC_PARSER__
#define __GENERIC_PARSER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FORWARD.h>
#include <string>
namespace PhysBAM{

template<class T>
class GENERIC_PARSER
{
private:
    ARRAY<std::string> lines;
    ARRAY<std::string>::iterator line;
    std::string partial_line;
public:

    GENERIC_PARSER(const std::string& filename,const int frame=0)
    {
        Preprocess_File(filename,frame);
        Reset_Parser();
    }
    
//##################################################################### 
    void Reset_Parser();
    bool Get_Statement(std::string& identifier,PARAMETER_LIST& list);
private:
    void Preprocess_File(std::string raw_filename,const int frame);
//##################################################################### 
};
}
#endif
#endif

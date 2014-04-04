//#####################################################################
// Copyright 2002, Ronald Fedkiw, Duc Nguyen, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MAIN  
//##################################################################### 
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include "Generic/GENERIC_RENDER_EXAMPLE.h"
#include "RAY_TRACING_DEBUG_DRIVER.h"

using namespace PhysBAM;
int main(int argc,char *argv[]) 
{
    PARSE_ARGS parse_args;
    parse_args.Add_Option_Argument("-preview");
    parse_args.Set_Extra_Arguments(-1, "<scene file> <frame number>");
    parse_args.Parse(argc,argv);
    if(parse_args.Num_Extra_Args() != 2){parse_args.Print_Usage();exit(0);}
    std::string scene_filename=parse_args.Extra_Arg(1);
    int frame_number=atoi(parse_args.Extra_Arg(2).c_str());

    GENERIC_RENDER_EXAMPLE<double,float> example(scene_filename,frame_number);
    RAY_TRACING_DEBUG_DRIVER<double> driver(example);
    driver.Execute_Main_Program();

    return 0;
}
//#####################################################################

//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROGRESS_INDICATOR
//#####################################################################
#ifndef __PROGRESS_INDICATOR__
#define __PROGRESS_INDICATOR__

#include <PhysBAM_Tools/Log/LOG.h>
namespace PhysBAM{

class PROGRESS_INDICATOR
{
public:
    int total;
    int done;
    int percent_done;
    bool print;

    PROGRESS_INDICATOR(const int total_input=1)
        :print(true)
    {
        Initialize(total_input);
    }

    void Initialize(const int total_input)
    {total=total_input;done=percent_done=0;}

    bool Progress(const int by=1)
    {done+=by;
    int new_percent_done=100*done/total;
    if(new_percent_done>percent_done){
        percent_done=new_percent_done;
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(print){LOG::cout<<percent_done<<"% "<<std::flush;if(percent_done==100) LOG::cout<<std::endl;}
#endif
        return true;}
    return false;}

//#####################################################################
};
}
#endif

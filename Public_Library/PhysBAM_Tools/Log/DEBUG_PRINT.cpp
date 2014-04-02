//#####################################################################
// Copyright 2006, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Macro PHYSBAM_DEBUG_PRINT
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <cstdarg>
namespace PhysBAM{
//#####################################################################
// Function Debug_Print_Helper
//#####################################################################
void Debug_Print_Helper(const char* prefix,...)
{
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    LOG::cout<<prefix<<": ";
    va_list marker;va_start(marker,prefix);
    bool first=true;
    for(;;){
        char* name=va_arg(marker,char*);if(!name) break;
        char* value=va_arg(marker,char*);if(!value) break;
        if(!first) LOG::cout<<", ";first=false;
        LOG::cout<<name<<"="<<value;}
    va_end(marker);
    LOG::cout<<std::endl;
#endif
}
//#####################################################################
}

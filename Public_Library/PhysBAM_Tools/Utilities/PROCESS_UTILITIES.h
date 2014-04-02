//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROCESS_UTILITIES
//#####################################################################
#ifndef __PROCESS_UTILITIES__
#define __PROCESS_UTILITIES__

namespace PhysBAM{
namespace PROCESS_UTILITIES{

unsigned int Memory_Usage();

void Set_Floating_Point_Exception_Handling(const bool enable,const bool division_by_zero=true,const bool invalid_operation=true,
    const bool overflow=true,const bool underflow=false,const bool inexact_result=false);
void Backtrace();
void Set_Backtrace(const bool enable=true);
void Block_Interrupt_Signal(const bool block=true);
void Sleep(const double seconds);

}
}
#endif

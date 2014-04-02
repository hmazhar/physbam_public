//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// File INTERRUPTS
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/INTERRUPTS.h>
namespace PhysBAM{
//#####################################################################
namespace {
ARRAY<void(*)()> interrupt_checkers;
}
//#####################################################################
// Function Check_For_Interrupts
//#####################################################################
void Check_For_Interrupts()
{
    for(int i=1;i<=interrupt_checkers.Size();i++)
        interrupt_checkers(i)();
}
//#####################################################################
// Function Check_For_Interrupts
//#####################################################################
void Add_Interrupt_Checker(void (*checker)())
{
    interrupt_checkers.Append(checker);
}
//#####################################################################
}

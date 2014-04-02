//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// File INTERRUPTS
//#####################################################################
#ifndef __INTERRUPTS__
#define __INTERRUPTS__

namespace PhysBAM{

// check if an interrupt has been posted, and throw an exception if so
void Check_For_Interrupts();

// add a function to call from Check_Interrupts
void Add_Interrupt_Checker(void (*checker)());

}
#endif

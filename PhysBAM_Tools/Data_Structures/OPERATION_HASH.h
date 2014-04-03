//#####################################################################
// Copyright 2004-2007, Geoffrey Irving, Craig Schroeder, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPERATION_HASH
//#####################################################################
#ifndef __OPERATION_HASH__
#define __OPERATION_HASH__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
namespace PhysBAM{

template<class ID> // ID=int
class OPERATION_HASH
{
public:
    ARRAY<unsigned int,ID> operations;
    unsigned int current_operation;

    OPERATION_HASH(const ID m=ID())
        :operations(m),current_operation(1)
    {}

    void Initialize(const ID m)
    {if(m==operations.Size()){Next_Operation();return;}
    operations.Resize(m,false,false);ARRAYS_COMPUTATIONS::Fill(operations,(unsigned int)0);current_operation=1;}
    
    void Mark(const ID i)
    {operations(i)=current_operation;}

    void Unmark(const ID i)
    {operations(i)=0;}

    bool Is_Marked_Current(const ID i) const
    {return operations(i)==current_operation;}
    
    void Next_Operation()
    {current_operation++;if(current_operation==0){current_operation=1;ARRAYS_COMPUTATIONS::Fill(operations,(unsigned int)0);}} // reset everything if overflow

    void Remove_Duplicates(ARRAY<ID>& list)
    {Next_Operation();
    for(int i=list.Size();i>=1;i--) if(Is_Marked_Current(list(i))) list.Remove_Index_Lazy(i);else Mark(list(i));
    Next_Operation();}

//#####################################################################
};
}
#endif


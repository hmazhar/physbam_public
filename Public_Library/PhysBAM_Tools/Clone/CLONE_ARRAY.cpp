//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CLONEABLE 
//#####################################################################
#include <PhysBAM_Tools/Clone/CLONE_ARRAY.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
CLONE_ARRAY<CLONEABLE_BASE>::
CLONE_ARRAY(const CLONEABLE_BASE& template_object,const int count)
    :sizeof_clone(template_object.Sizeof_Clone()),count(count),data(new char[sizeof_clone*count])
{
    int i=1;
    try{
        for(;i<=count;i++)
            template_object.Placement_Clone(&(*this)(i));}
    catch(...){
        for(int j=1;j<i;j++)
            (*this)(j).~CLONEABLE_BASE();
        throw;}
}
//#####################################################################
// Destructor
//#####################################################################
CLONE_ARRAY<CLONEABLE_BASE>::
~CLONE_ARRAY()
{
    for(int i=1;i<=count;i++)
        (*this)(i).~CLONEABLE_BASE();
    delete[] data;
}
//#####################################################################
}

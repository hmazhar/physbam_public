//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMBOL
//#####################################################################
#include <PhysBAM_Tools/Utilities/SYMBOL.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
SYMBOL::
SYMBOL(const char* name)
    :key(Mapping().Name_To_Key(name))
{}
SYMBOL::
SYMBOL(const std::string& name)
    :key(Mapping().Name_To_Key(name))
{}
//#####################################################################
// Function Name
//#####################################################################
const std::string& SYMBOL::
Name() const
{
    return Mapping().key_to_name.Get(key);
}
//#####################################################################
// Constructor
//#####################################################################
SYMBOL_MAPPING::
SYMBOL_MAPPING()
{
#ifdef USE_PTHREADS
    PHYSBAM_NOT_IMPLEMENTED("Thread-safety for SYMBOL mapping");
#endif
    Name_To_Key(""); // map the empty string to 0
}
//#####################################################################
// Function Name_To_Key
//#####################################################################
int SYMBOL_MAPPING::
Name_To_Key(const std::string& name)
{
    if(int* i=name_to_key.Get_Pointer(name))
        return *i;
    int key=name_to_key.Size();
    name_to_key.Insert(name,key);
    key_to_name.Insert(key,name);
    return key;
}
//#####################################################################
// Function Mapping
//#####################################################################
SYMBOL_MAPPING&
Mapping()
{
    static SYMBOL_MAPPING mapping;
    return mapping;
}
//#####################################################################
}

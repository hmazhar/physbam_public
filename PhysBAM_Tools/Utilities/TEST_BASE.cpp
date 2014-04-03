//#####################################################################
// Copyright 2007, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TEST_BASE
//#####################################################################
#include <PhysBAM_Tools/Utilities/TEST_BASE.h>
namespace PhysBAM{

HASHTABLE<std::string,TEST_BASE*>& TEST_BASE::Test_Registry()
{
    static HASHTABLE<std::string,TEST_BASE*> registry;
    return registry;
}
}

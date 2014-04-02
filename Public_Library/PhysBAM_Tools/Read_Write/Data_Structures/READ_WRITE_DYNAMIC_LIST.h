//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_DYNAMIC_LIST
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_DYNAMIC_LIST__
#define __READ_WRITE_DYNAMIC_LIST__

#include <PhysBAM_Tools/Data_Structures/DYNAMIC_LIST.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW>
class Read_Write<DYNAMIC_LIST_CORE,RW>
{
public:
    // Reads the set of active id's and updates the index<->id maps.
    // needs_init is filled with id's of those elements which have become newly active and need to be initialized.
    static void Read(DYNAMIC_LIST_CORE& object,const std::string& prefix,ARRAY<int>& needs_init);

    // Writes the set of active id's.
    // needs_write indicates which elements have not been written since their creation, so derived classes should write those out (and reset the needs_write list).
    static void Write(const DYNAMIC_LIST_CORE& object,const std::string& prefix);
};

template<class RW,class T,class ID>
class Read_Write<DYNAMIC_LIST<T,ID>,RW>
{
public:
    static void Read(DYNAMIC_LIST<T,ID>& object,const std::string& prefix,ARRAY<ID>& needs_init)
    {Read_Write<DYNAMIC_LIST_CORE,RW>::Read(object.core,prefix,reinterpret_cast<ARRAY<int>&>(needs_init));}

    static void Write(const DYNAMIC_LIST<T,ID>& object,const std::string& prefix)
    {Read_Write<DYNAMIC_LIST_CORE,RW>::Write(object.core,prefix);}
 };
}
#endif
#endif

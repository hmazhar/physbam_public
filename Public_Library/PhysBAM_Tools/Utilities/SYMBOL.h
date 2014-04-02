//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMBOL
//#####################################################################
//
// Wraps a string so that integer comparison can be used instead of string comparison.
// Converting a string into a SYMBOL requires a hashtable lookup / modification.
// Symbols have a total order, but it is creation order dependent and unrelated to the lexicographic ordering of strings.
//
//#####################################################################
#ifndef __SYMBOL__
#define __SYMBOL__

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Math_Tools/Hash.h>
namespace PhysBAM{

class SYMBOL
{
protected:
    int key; // a key into an internal hash table
public:
    SYMBOL() // equivalent to SYMBOL("")
        :key(0)
    {}

    SYMBOL(const char* name_string);
    SYMBOL(const std::string& name_string);

    bool operator==(const SYMBOL symbol) const
    {return key==symbol.key;}

    bool operator!=(const SYMBOL symbol) const
    {return key!=symbol.key;}

    bool operator<(const SYMBOL symbol) const // not related to string comparison
    {return key<symbol.key;}

    bool operator>(const SYMBOL symbol) const
    {return key>symbol.key;}

    bool operator<=(const SYMBOL symbol) const
    {return key<=symbol.key;}

    bool operator>=(const SYMBOL symbol) const
    {return key>=symbol.key;}

private:
    struct UNUSABLE{void F(){}};
    typedef void (UNUSABLE::*SAFE_BOOL)();
public:

    operator SAFE_BOOL() const // allow conversion to bool without allowing conversion to T
    {return key?&UNUSABLE::F:0;} // empty string always maps to 0

    friend inline int Hash_Reduce(const SYMBOL symbol)
    {return symbol.key;}

//#####################################################################
    const std::string& Name() const;
//#####################################################################
};
//#####################################################################
struct SYMBOL_MAPPING
{
    HASHTABLE<int,std::string> key_to_name;
    HASHTABLE<std::string,int> name_to_key;

    SYMBOL_MAPPING();
    int Name_To_Key(const std::string& name);
};
SYMBOL_MAPPING& Mapping();
//#####################################################################
}
#endif

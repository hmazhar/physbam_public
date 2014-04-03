//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXPANSION_ARITHMETIC
//#####################################################################
#ifndef __EXPANSION_ARITHMETIC__
#define __EXPANSION_ARITHMETIC__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>

namespace PhysBAM{

template <class T>
class EXPANSION_ARITHMETIC{
public:
    static const T splitter;

//#####################################################################
    static void Check_IEEE_Compliance();
    static void Add_Expansions(const ARRAY<T>& e,const ARRAY<T>& f,ARRAY<T>& h);
    static void Scale_Expansion(const ARRAY<T>& e,const T b,ARRAY<T>& h);
    static void Compress_Expansion(ARRAY<T>& e);
//#####################################################################
};
}
#endif

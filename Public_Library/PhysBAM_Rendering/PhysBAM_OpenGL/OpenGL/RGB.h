//#####################################################################
// Copyright 2002-2005, Geoffrey Irving, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RGB
//#####################################################################
#ifndef __RGB__
#define __RGB__

namespace PhysBAM{

template<class T>
class RGB
{
public:
    T r,g,b;

    RGB()
        :r(0),g(0),b(0)
    {}

    RGB(const T r_input,const T g_input,const T b_input)
        :r(r_input),g(g_input),b(b_input)
    {}

//#####################################################################
};
}
#endif

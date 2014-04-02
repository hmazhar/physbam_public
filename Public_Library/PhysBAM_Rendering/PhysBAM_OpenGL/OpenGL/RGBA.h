//#####################################################################
// Copyright 2002-2005, Geoffrey Irving, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RGBA
//#####################################################################
#ifndef __RGBA__
#define __RGBA__

namespace PhysBAM{

template<class T>
class RGBA
{
public:
    T r,g,b,a;

    RGBA()
        :r(0),r(0),r(0),r(0)
    {}

    RGBA(const T r_input,const T g_input,const T b_input,const T a_input)
        :r(r_input),g(g_input),b(b_input),a(a_input)
    {}

//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IDENTITY_ARRAY
//#####################################################################
#ifndef __IDENTITY_ARRAY__
#define __IDENTITY_ARRAY__

#include <PhysBAM_Tools/Arrays/ARRAY_BASE.h>
#include <cassert>
namespace PhysBAM{

template<class ID> struct IS_ARRAY<IDENTITY_ARRAY<ID> > {static const bool value=true;};
template<class ID> struct IS_ARRAY_VIEW<IDENTITY_ARRAY<ID> > {static const bool value=true;};

template<class ID> // ID=int
class IDENTITY_ARRAY:public ARRAY_BASE<ID,IDENTITY_ARRAY<ID>,ID>
{
public:
    typedef ID ELEMENT;
    typedef ID INDEX;
private:
    ID m;
public:

    explicit IDENTITY_ARRAY(const ID m)
        :m(m)
    {}

    ID Size() const
    {return m;}

    ID operator()(const ID i) const
    {assert(ID(1)<=i && i<=m);return i;}

//#####################################################################
};
}
#endif

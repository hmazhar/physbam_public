//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONCOPYABLE 
//#####################################################################
#ifndef __NONCOPYABLE__
#define __NONCOPYABLE__

namespace PhysBAM{

class NONCOPYABLE
{
protected:
    NONCOPYABLE(){}
    ~NONCOPYABLE(){}
private:
    NONCOPYABLE(const NONCOPYABLE&);
    void operator=(const NONCOPYABLE&);
};
}
#endif

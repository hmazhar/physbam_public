//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class THREAD_PACKAGE
//#####################################################################
#ifndef __THREAD_PACKAGE__
#define __THREAD_PACKAGE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>

namespace PhysBAM{

class THREAD_PACKAGE
{
public:
    ARRAY<char> buffer;
    int send_tid,recv_tid;

    THREAD_PACKAGE()
        :send_tid(0),recv_tid(0)
    {}

    THREAD_PACKAGE(const int size)
        :buffer(size),send_tid(0),recv_tid(0)
    {}

    THREAD_PACKAGE(const THREAD_PACKAGE& pack)
        :buffer(pack.buffer),send_tid(pack.send_tid),recv_tid(pack.recv_tid)
    {}

//#####################################################################
};
}
#endif

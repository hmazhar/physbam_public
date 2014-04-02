//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Classes RLE_BLOCK_RUN_2D and RLE_BLOCK_RUN_3D
//##################################################################### 
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_BLOCK_RUNS__
#define __RLE_BLOCK_RUNS__

namespace PhysBAM{

struct RLE_BLOCK_RUN_2D
{
    short i,jmin,jmax;
    int block;
    int cells[2];
    int faces_x[3];
};

struct RLE_BLOCK_RUN_3D
{
    short ij,jmin,jmax;
    int block;
    int cells[4];
    int faces_x[6];
    int faces_z[6];
};
}
#endif
#endif

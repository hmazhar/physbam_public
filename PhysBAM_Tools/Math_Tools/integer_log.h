//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function integer_log and integer_log_exact
//#####################################################################
//
// base 2 logarithms of integers
//               
//#####################################################################
#ifndef __integer_log__
#define __integer_log__

#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <cassert>
namespace PhysBAM{

inline int integer_log(unsigned int v) // this works for any v, but it is slower
{int c=0;
if(v&0xffff0000){v>>=16;c|=16;}
if(v&0xff00){v>>=8;c|=8;}
if(v&0xf0){v>>=4;c|=4;}
if(v&0xc){v>>=2;c|=2;}
if(v&2)c|=1;return c;}

inline bool power_of_two(const unsigned int v)
{return v>0 && (v&(v-1))==0;}

inline unsigned int next_power_of_two(unsigned int v)
{STATIC_ASSERT(sizeof(v)==4);v--;
v|=v>>1;v|=v>>2;v|=v>>4;v|=v>>8;v|=v>>16;return v+1;}

inline int rightmost_bit(unsigned int v)
{return v&(unsigned int)-(int)v;} // TODO: VS workaround

inline int integer_log_exact(const unsigned int v) // this only works if v is a power of 2
{int log_value= (((v&0xffff0000)!=0)<<4)+(((v&0xff00ff00)!=0)<<3)+(((v&0xf0f0f0f0)!=0)<<2)+(((v&0xcccccccc)!=0)<<1)+((v&0xaaaaaaaa)!=0);
assert(v==(unsigned int)1<<log_value);return log_value;}

}
#endif

//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Geoffrey Irving, Andrew Selle, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Permutation group operations
//#####################################################################
//
// Permutations of 2, 3 or 4 values.
// Given a reference triple or quadraple, returns the ith permutation of it.
//
// Different Permuation Groups can be added here.
//
// Cyclic group is a group consisting of right cyclic shifts.
//               
//#####################################################################
#ifndef __permutation__
#define __permutation__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{


template<class T>
inline VECTOR<T,2> permute_two(const VECTOR<T,2>& r,const int ith_permutation) 
{assert(ith_permutation>=1 && ith_permutation<=2);
switch(ith_permutation){
    case 1:return r;
    default:return VECTOR<T,2>(r[2],r[1]);}} // case 2

template<class T>
inline VECTOR<T,2> permute_two_inverse(const VECTOR<T,2>& r,const int ith_permutation) 
{assert(ith_permutation>=1 && ith_permutation<=2);
    return permute_two(r,ith_permutation);} // case 2

template<class T>
inline VECTOR<T,3> permute_three(const VECTOR<T,3>& r,const int ith_permutation) 
{assert(ith_permutation>=1 && ith_permutation<=6);
switch(ith_permutation){
    case 1:return VECTOR<T,3>(r[1],r[2],r[3]);
    case 2:return VECTOR<T,3>(r[1],r[3],r[2]);
    case 3:return VECTOR<T,3>(r[2],r[1],r[3]);
    case 4:return VECTOR<T,3>(r[2],r[3],r[1]);
    case 5:return VECTOR<T,3>(r[3],r[1],r[2]);
    default:return VECTOR<T,3>(r[3],r[2],r[1]);}} // case 6

template<class T>
inline VECTOR<T,3> permute_three_inverse(const VECTOR<T,3>& r,const int ith_permutation) 
{assert(ith_permutation>=1 && ith_permutation<=6);
switch(ith_permutation){
    case 1:return VECTOR<T,3>(r[1],r[2],r[3]);
    case 2:return VECTOR<T,3>(r[1],r[3],r[2]);
    case 3:return VECTOR<T,3>(r[2],r[1],r[3]);
    case 4:return VECTOR<T,3>(r[3],r[1],r[2]);
    case 5:return VECTOR<T,3>(r[2],r[3],r[1]);
    default:return VECTOR<T,3>(r[3],r[2],r[1]);}} // case 6

template<class T>
inline VECTOR<T,4> permute_four(const VECTOR<T,4>& r,const int ith_permutation) 
{assert(ith_permutation>=1 && ith_permutation<=24);
switch(ith_permutation){
    case 1:return VECTOR<T,4>(r[1],r[2],r[3],r[4]);
    case 2:return VECTOR<T,4>(r[1],r[2],r[4],r[3]);
    case 3:return VECTOR<T,4>(r[1],r[3],r[2],r[4]);
    case 4:return VECTOR<T,4>(r[1],r[3],r[4],r[2]);
    case 5:return VECTOR<T,4>(r[1],r[4],r[2],r[3]);
    case 6:return VECTOR<T,4>(r[1],r[4],r[3],r[2]);
    case 7:return VECTOR<T,4>(r[2],r[1],r[3],r[4]);
    case 8:return VECTOR<T,4>(r[2],r[1],r[4],r[3]);
    case 9:return VECTOR<T,4>(r[2],r[3],r[1],r[4]);
    case 10:return VECTOR<T,4>(r[2],r[3],r[4],r[1]);
    case 11:return VECTOR<T,4>(r[2],r[4],r[1],r[3]);
    case 12:return VECTOR<T,4>(r[2],r[4],r[3],r[1]);
    case 13:return VECTOR<T,4>(r[3],r[1],r[2],r[4]);
    case 14:return VECTOR<T,4>(r[3],r[1],r[4],r[2]);
    case 15:return VECTOR<T,4>(r[3],r[2],r[1],r[4]);
    case 16:return VECTOR<T,4>(r[3],r[2],r[4],r[1]);
    case 17:return VECTOR<T,4>(r[3],r[4],r[1],r[2]);
    case 18:return VECTOR<T,4>(r[3],r[4],r[2],r[1]);
    case 19:return VECTOR<T,4>(r[4],r[1],r[2],r[3]);
    case 20:return VECTOR<T,4>(r[4],r[1],r[3],r[2]);
    case 21:return VECTOR<T,4>(r[4],r[2],r[1],r[3]);
    case 22:return VECTOR<T,4>(r[4],r[2],r[3],r[1]);
    case 23:return VECTOR<T,4>(r[4],r[3],r[1],r[2]);
    default:return VECTOR<T,4>(r[4],r[3],r[2],r[1]);}} // case 24

template<class T>
inline VECTOR<T,4> permute_four_inverse(const VECTOR<T,4>& r,const int ith_permutation) 
{assert(ith_permutation>=1 && ith_permutation<=24);
switch(ith_permutation){
    case 1:return VECTOR<T,4>(r[1],r[2],r[3],r[4]);
    case 2:return VECTOR<T,4>(r[1],r[2],r[4],r[3]);
    case 3:return VECTOR<T,4>(r[1],r[3],r[2],r[4]);
    case 4:return VECTOR<T,4>(r[1],r[4],r[2],r[3]);
    case 5:return VECTOR<T,4>(r[1],r[3],r[4],r[2]);
    case 6:return VECTOR<T,4>(r[1],r[4],r[3],r[2]);
    case 7:return VECTOR<T,4>(r[2],r[1],r[3],r[4]);
    case 8:return VECTOR<T,4>(r[2],r[1],r[4],r[3]);
    case 9:return VECTOR<T,4>(r[3],r[1],r[2],r[4]);
    case 10:return VECTOR<T,4>(r[4],r[1],r[2],r[3]);
    case 11:return VECTOR<T,4>(r[3],r[1],r[4],r[2]);
    case 12:return VECTOR<T,4>(r[4],r[1],r[3],r[2]);
    case 13:return VECTOR<T,4>(r[2],r[3],r[1],r[4]);
    case 14:return VECTOR<T,4>(r[2],r[4],r[1],r[3]);
    case 15:return VECTOR<T,4>(r[3],r[2],r[1],r[4]);
    case 16:return VECTOR<T,4>(r[4],r[2],r[1],r[3]);
    case 17:return VECTOR<T,4>(r[3],r[4],r[1],r[2]);
    case 18:return VECTOR<T,4>(r[4],r[3],r[1],r[2]);
    case 19:return VECTOR<T,4>(r[2],r[3],r[4],r[1]);
    case 20:return VECTOR<T,4>(r[2],r[4],r[3],r[1]);
    case 21:return VECTOR<T,4>(r[3],r[2],r[4],r[1]);
    case 22:return VECTOR<T,4>(r[4],r[2],r[3],r[1]);
    case 23:return VECTOR<T,4>(r[3],r[4],r[2],r[1]);
    default:return VECTOR<T,4>(r[4],r[3],r[2],r[1]);}} // case 24

inline void next_cyclic_permutation_of_four_index(int& i)
{switch(i){
    case 1:i=19;return;
    case 19:i=17;return;
    case 17:i=10;return;
    default:assert(i==10);i=1;return;}}

inline bool permutation_of_four_is_even(const int ith_permutation) 
{assert(ith_permutation>=1 && ith_permutation<=24);return !(ith_permutation&2);}

inline bool permutation_of_three_is_even(const int ith_permutation)
{assert(ith_permutation>=1 && ith_permutation<=6);return !(ith_permutation&2);}

inline bool permutation_of_two_is_even(const int ith_permutation)
{assert(ith_permutation>=1 && ith_permutation<=2);return !(ith_permutation&2);}
}
#endif

//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Craig Schroeder, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header DATA_STRUCTURES_FORWARD
//#####################################################################
#ifndef __DATA_STRUCTURES_FORWARD__
#define __DATA_STRUCTURES_FORWARD__

namespace PhysBAM{

template<class T1,class T2> class PAIR;
template<class T1,class T2,class T3> class TRIPLE;
template<class TK,class T=void> class HASHTABLE;
template<class TK,class T=void> class HASHTABLE_ITERATOR;
template<class ID=int> class UNION_FIND;
template<class ID=int> class SPARSE_UNION_FIND;

template<class ID,class T,int flags> class ELEMENT_ID;
template<class ID=int> class OPERATION_HASH;
template<class T> class QUEUE;
template<class T> class STACK;

template<class ID> class DIRECTED_GRAPH;
template<class ID,class EID> class UNDIRECTED_GRAPH;

struct INITIAL_SIZE;
}
#endif

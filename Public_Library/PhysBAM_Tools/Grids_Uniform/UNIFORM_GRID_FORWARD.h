//#####################################################################
// Copyright 2009, Nipun Kwatra
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __UNIFORM_GRID_FORWARD__
#define __UNIFORM_GRID_FORWARD__

namespace PhysBAM{

template<class TV> class GRID;

template<class T> class UNIFORM_GRID_ITERATOR;
template<class T> class UNIFORM_GRID_ITERATOR_NODE;
template<class T> class UNIFORM_GRID_ITERATOR_CELL;
template<class T> class UNIFORM_GRID_ITERATOR_FACE;
template<class TV> struct UNIFORM_TAG{};

template<class T_GRID> class BLOCK_UNIFORM;

}
#endif

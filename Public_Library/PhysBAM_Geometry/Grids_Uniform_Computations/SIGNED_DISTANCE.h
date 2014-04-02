//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace SIGNED_DISTANCE
//##################################################################### 
#ifndef __SIGNED_DISTANCE__
#define __SIGNED_DISTANCE__

namespace PhysBAM{
namespace SIGNED_DISTANCE{
//#####################################################################
template<class T_GRID,class T_ARRAY,class T_GEOMETRY> void Calculate(T_GEOMETRY& geometry,const T_GRID& grid,T_ARRAY& phi,bool print_progress=false){
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::VECTOR_INT TV_INT;
    for(NODE_ITERATOR iterator(grid,grid.Domain_Indices());iterator.Valid();iterator.Next()){TV_INT index=iterator.Node_Index();
        phi(index)=geometry.Signed_Distance(grid.X(index));}
}
//#####################################################################
};
};
#endif

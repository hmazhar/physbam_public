#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_QUADTREE_GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_QUADTREE_NODE_VECTOR_FIELD.h>
using namespace PhysBAM;

template<class T> void OPENGL_QUADTREE_NODE_VECTOR_FIELD<T>::
Update()
{
    ARRAY<bool> active(grid.number_of_nodes);
    int number_of_active_nodes=0;
    for(DYADIC_GRID_ITERATOR_NODE<QUADTREE_GRID<T> > iterator(grid,grid.Map_All_Nodes());iterator.Valid();iterator.Next()){
        active(iterator.Node_Index())=true;number_of_active_nodes++;}

    vector_field.Resize(number_of_active_nodes);
    vector_locations.Resize(number_of_active_nodes);

    int current_node=0;
    for(int i=1;i<=grid.number_of_nodes;i++)if(active(i)){current_node++;
        vector_locations(current_node)=grid.Node_Location(i);
        vector_field(current_node)=V(i);}
}

template<class T> RANGE<VECTOR<float,3> > OPENGL_QUADTREE_NODE_VECTOR_FIELD<T>::
Bounding_Box() const
{
    RANGE<VECTOR<float,3> > box(grid.uniform_grid.domain.min_corner.x,grid.uniform_grid.domain.max_corner.x,grid.uniform_grid.domain.min_corner.y,grid.uniform_grid.domain.max_corner.y,0,0);
    return World_Space_Box(box);
}

template<class T> void OPENGL_QUADTREE_NODE_VECTOR_FIELD<T>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::QUADTREE_NODE){
        int index=((OPENGL_SELECTION_QUADTREE_NODE<T>*)current_selection)->index;
        if(V.Valid_Index(index))
            output_stream<<V(index);}
    output_stream<<std::endl;
}


template class OPENGL_QUADTREE_NODE_VECTOR_FIELD<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_QUADTREE_NODE_VECTOR_FIELD<double>;
#endif
#endif

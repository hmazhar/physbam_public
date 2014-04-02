#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OCTREE_GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OCTREE_NODE_VECTOR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OCTREE_SLICE.h>
using namespace PhysBAM;

//#####################################################################
// Update
//#####################################################################
template<class T> void OPENGL_OCTREE_NODE_VECTOR_FIELD<T>::
Update()
{
    vector_field.Resize(0);
    vector_locations.Resize(0);

    if(V.m==0) return;

    ARRAY<bool> active(grid.number_of_nodes);
    if(this->slice && this->slice->Is_Slice_Mode()){
        int number_of_active_nodes=0;
        OPENGL_OCTREE_SLICE* slice=(OPENGL_OCTREE_SLICE*)this->slice;
        for(DYADIC_GRID_ITERATOR_NODE<OCTREE_GRID<T> > iterator(grid,grid.Map_All_Nodes());iterator.Valid();iterator.Next())
            if(abs(iterator.Location()[slice->axis]-slice->position)<.00001){number_of_active_nodes++;active(iterator.Node_Index())=true;}

        vector_field.Resize(number_of_active_nodes);
        vector_locations.Resize(number_of_active_nodes);

        int current_node=0;
        for(int i=1;i<=grid.number_of_nodes;i++)if(active(i)){current_node++;
            vector_locations(current_node)=grid.Node_Location(i);
            vector_field(current_node)=V(i);}}
    else{
        vector_field.Resize(grid.number_of_nodes);
        vector_locations.Resize(grid.number_of_nodes);
        for(int i=1;i<=grid.number_of_nodes;i++){
            vector_locations(i)=grid.Node_Location(i);
            vector_field(i)=V(i);}}
}
//#####################################################################
// Slice_Has_Changed
//#####################################################################
template<class T> void OPENGL_OCTREE_NODE_VECTOR_FIELD<T>::
Slice_Has_Changed()
{
    Update();
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_OCTREE_NODE_VECTOR_FIELD<T>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<VECTOR<float,3> >(grid.uniform_grid.domain));
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_OCTREE_NODE_VECTOR_FIELD<T>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::OCTREE_NODE){
        int index=((OPENGL_SELECTION_OCTREE_NODE<T>*)current_selection)->index;
        output_stream<<V(index);}
    output_stream<<std::endl;
}
//#####################################################################
template class OPENGL_OCTREE_NODE_VECTOR_FIELD<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_OCTREE_NODE_VECTOR_FIELD<double>;
#endif
#endif

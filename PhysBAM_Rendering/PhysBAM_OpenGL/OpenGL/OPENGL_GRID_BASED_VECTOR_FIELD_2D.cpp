//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_BASED_VECTOR_FIELD_2D.h>
using namespace PhysBAM;
//#####################################################################
// OPENGL_GRID_BASED_VECTOR_FIELD_2D
//#####################################################################
template<class T> OPENGL_GRID_BASED_VECTOR_FIELD_2D<T>::
OPENGL_GRID_BASED_VECTOR_FIELD_2D(GRID<TV>& grid, ARRAY<VECTOR<T,2>,VECTOR<int,2> >& V)
    :OPENGL_VECTOR_FIELD_2D<ARRAY<TV> >(vector_field,vector_locations),grid(grid),V(V)
{}
//#####################################################################
// ~OPENGL_GRID_BASED_VECTOR_FIELD_2D
//#####################################################################
template<class T> OPENGL_GRID_BASED_VECTOR_FIELD_2D<T>::
~OPENGL_GRID_BASED_VECTOR_FIELD_2D()
{}
//#####################################################################
// Update
//#####################################################################
template<class T> void OPENGL_GRID_BASED_VECTOR_FIELD_2D<T>::
Update()
{
    int idx=1;
    vector_field.Resize(V.counts.Product());
    vector_locations.Resize(V.counts.Product());
    for(int i=V.domain.min_corner.x;i<=V.domain.max_corner.x;i++)for(int j=V.domain.min_corner.y;j<=V.domain.max_corner.y;j++){
        vector_field(idx)=V(i,j);vector_locations(idx)=grid.X(i,j);idx++;}
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_GRID_BASED_VECTOR_FIELD_2D<T>::
Bounding_Box() const
{
    return RANGE<VECTOR<float,3> >(grid.domain.min_corner.x,grid.domain.max_corner.x,grid.domain.min_corner.y,grid.domain.max_corner.y,0,0);
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_GRID_BASED_VECTOR_FIELD_2D<T>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* current_selection) const
{
    // TODO: interpolate to particles
    if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_NODE_2D && !grid.Is_MAC_Grid()){
        VECTOR<int,2> index=((OPENGL_SELECTION_GRID_NODE_2D<T>*)current_selection)->index;
        if(V.Valid_Index(index)) stream<<V(index);}
    if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_CELL_2D && grid.Is_MAC_Grid()){
        VECTOR<int,2> index=((OPENGL_SELECTION_GRID_CELL_2D<T>*)current_selection)->index;
        if(V.Valid_Index(index)) stream<<V(index);}
    stream<<std::endl;
}
//#####################################################################
template class OPENGL_GRID_BASED_VECTOR_FIELD_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_GRID_BASED_VECTOR_FIELD_2D<double>;
#endif

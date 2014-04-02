#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OCTREE_GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OCTREE_NODE_SCALAR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OCTREE_SLICE.h>
using namespace PhysBAM;

//#####################################################################
// Display
//#####################################################################
template<class T,class T2> void OPENGL_OCTREE_NODE_SCALAR_FIELD<T,T2>::
Display(const int in_color) const
{
    if(!slice) return;

    glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT);
    glPointSize(point_size);glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);
    OPENGL_OCTREE_SLICE* slice=(OPENGL_OCTREE_SLICE*)this->slice;

    OpenGL_Begin(GL_POINTS);
    if(slice->Is_Slice_Mode()) for(int i=1;i<=grid.number_of_nodes;i++){
            if(abs(grid.Node_Location(i)[slice->axis]-(T)slice->position)<.00001){
                color_map->Lookup(value(i)).Send_To_GL_Pipeline();
                OpenGL_Vertex(grid.Node_Location(i));}}
    else for(int i=1;i<=grid.number_of_nodes;i++){
            color_map->Lookup(value(i)).Send_To_GL_Pipeline();
            OpenGL_Vertex(grid.Node_Location(i));}
    OpenGL_End();
    glPopAttrib();
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<float,3> > OPENGL_OCTREE_NODE_SCALAR_FIELD<T,T2>::
Bounding_Box() const
{
    return RANGE<VECTOR<float,3> >(grid.uniform_grid.domain);
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_OCTREE_NODE_SCALAR_FIELD<T,T2>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::OCTREE_NODE){
        int index=((OPENGL_SELECTION_OCTREE_NODE<T>*)current_selection)->index;
        output_stream<<value(index);}
    output_stream<<std::endl;
}
//#####################################################################
template class OPENGL_OCTREE_NODE_SCALAR_FIELD<float,int>;
template class OPENGL_OCTREE_NODE_SCALAR_FIELD<float,bool>;
template class OPENGL_OCTREE_NODE_SCALAR_FIELD<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_OCTREE_NODE_SCALAR_FIELD<double,int>;
template class OPENGL_OCTREE_NODE_SCALAR_FIELD<double,bool>;
template class OPENGL_OCTREE_NODE_SCALAR_FIELD<double,double>;
#endif
#endif

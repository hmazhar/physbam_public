//#####################################################################
// Copyright 2002-2005, Eilene Hao, Geoffrey Irving, Sergey Koltakov, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_HEXAHEDRALIZED_VOLUME
//##################################################################### 
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_HEXAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
using namespace PhysBAM;
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_HEXAHEDRALIZED_VOLUME<T>::
Display(const int in_color) const
{
    if(boundary_only) Draw_Boundary_Triangles(*hexahedron_mesh);
    else Draw_Wireframe_Mesh(*hexahedron_mesh);
    if(draw_subsets) {Draw_Subset_Particles();Draw_Vector_At_Hex_Center();}
}
//#####################################################################
// Function Draw_Subset_Particles
//#####################################################################
template<class T> void OPENGL_HEXAHEDRALIZED_VOLUME<T>::
Draw_Subset_Particles() const
{
    for(int p=1;p<=subset_particles.m;p++) OPENGL_SHAPES::Draw_Dot(particles->X(subset_particles(p)),OPENGL_COLOR(1,1,0));
}
//#####################################################################
// Function Draw_Vector_At_Hex_Center
//#####################################################################
template<class T> void OPENGL_HEXAHEDRALIZED_VOLUME<T>::
Draw_Vector_At_Hex_Center() const
{
    glMatrixMode(GL_MODELVIEW);glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    VECTOR<T,3> head;
    glMatrixMode(GL_MODELVIEW);

    glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_LINE_BIT | GL_CURRENT_BIT);
    OPENGL_COLOR(0,1,0).Send_To_GL_Pipeline();
    glLineWidth(1);glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);
    OpenGL_Begin(GL_LINES);
    for(int i=1;i<=vectors_at_hex_centers.m;i++){
        ARRAY<int> p(8);
        VECTOR<T,3> hex_center=VECTOR<T,3>(0,0,0);hexahedron_mesh->elements(i).Get(p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8));
        for(int j=1;j<=8;j++) hex_center+=particles->X(p(j));hex_center*=(T).125;
        head=hex_center+(T)vector_size*vectors_at_hex_centers(i);
        OpenGL_Line(hex_center,head);
        VECTOR<T,3> orth_vect=vectors_at_hex_centers(i).Orthogonal_Vector();
        orth_vect*=.15*vector_size;
        OpenGL_Line(head,head+orth_vect-(T).15*(T)vector_size*vectors_at_hex_centers(i));
        OpenGL_Line(head,head-orth_vect-(T).15*(T)vector_size*vectors_at_hex_centers(i));}
    OpenGL_End();glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Function Draw_Wireframe_Mesh
//#####################################################################
template<class T> void OPENGL_HEXAHEDRALIZED_VOLUME<T>::
Draw_Wireframe_Mesh(const HEXAHEDRON_MESH& hexahedron_mesh) const
{
    glDisable(GL_LIGHTING);
    material.diffuse.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_LINES); 
    for(int h=1;h<=hexahedron_mesh.elements.m;h++)for(int k=0;k<24;k++)
        OpenGL_Vertex(particles->X(hexahedron_mesh.elements(h)(HEXAHEDRON_MESH::edge_indices[0][k])));
    OpenGL_End();
    glEnable(GL_LIGHTING);
}   
//#####################################################################
// Function Draw_Boundary_Triangles
//#####################################################################
template<class T> void OPENGL_HEXAHEDRALIZED_VOLUME<T>::
Draw_Boundary_Triangles(const HEXAHEDRON_MESH& hexahedron_mesh) const
{   
    TRIANGLE_MESH& mesh=*hexahedron_mesh.boundary_mesh;
    material.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_TRIANGLES);
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k);
        OpenGL_Normal(TRIANGLE_3D<T>::Normal(xi,xj,xk));
        OpenGL_Triangle(xi,xj,xk);}
    OpenGL_End();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_HEXAHEDRALIZED_VOLUME<T>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<VECTOR<float,3> >(RANGE<VECTOR<T,3> >::Bounding_Box(particles->X)));
}
//#####################################################################
template class OPENGL_HEXAHEDRALIZED_VOLUME<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_HEXAHEDRALIZED_VOLUME<double>;
#endif

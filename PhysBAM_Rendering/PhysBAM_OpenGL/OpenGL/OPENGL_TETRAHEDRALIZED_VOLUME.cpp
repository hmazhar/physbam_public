//#####################################################################
// Copyright 2005-2009, Zhaosheng Bao, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Sergey Koltakov, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TETRAHEDRALIZED_VOLUME.h>
using namespace PhysBAM;
namespace{
static VECTOR<int,4> Spring_Nodes(unsigned char pair_id,const VECTOR<int,4>& n)
{
    switch(pair_id){
        case 1: return VECTOR<int,4>(n[1],n[2],n[4],n[3]); // point face
        case 2: return VECTOR<int,4>(n[2],n[1],n[3],n[4]); // point face
        case 3: return VECTOR<int,4>(n[3],n[1],n[4],n[2]); // point face
        case 4: return VECTOR<int,4>(n[4],n[1],n[2],n[3]); // point face
        case 5: return VECTOR<int,4>(n[1],n[2],n[3],n[4]); // edge edge
        case 6: return VECTOR<int,4>(n[2],n[3],n[1],n[4]); // edge edge
        case 7: return VECTOR<int,4>(n[1],n[3],n[4],n[2]); // edge edge
        default: PHYSBAM_FATAL_ERROR();}
}
}
//#####################################################################
// Function Display_In_Color
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Display(const int in_color) const
{
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    GLint mode;
    glGetIntegerv(GL_RENDER_MODE,&mode);
    glDisable(GL_CULL_FACE);
    if(mode==GL_SELECT){
        glPushName(1);
        Draw_Vertices_For_Selection();
        glLoadName(2);
        Draw_Tetrahedra_For_Selection();
        glPopName();}
    else if(cutaway_mode){
        if(boundary_only) Draw_Boundary_Triangles(cutaway_mesh);
        else Draw_Wireframe_Mesh(cutaway_mesh);}
    else{
        if(spectrum.m) Draw_In_Color_From_Spectrum();
        if(color_map) Draw_In_Color_From_Color_Map();
        if(boundary_only) Draw_Boundary_Triangles(*mesh);
        else Draw_Wireframe_Mesh(*mesh);
        if(draw_subsets){Draw_Subset();Draw_Subset_Particles();Draw_Subset_Triangles();}
        //Draw_Current_Tetrahedron();
        //Highlight_Boundary_Nodes_Of_Current_Tetrahedron();
        //Highlight_Boundary_Normal_Vectors_Of_Current_Tetrahedron();
        //Highlight_Current_Node();
        //Highlight_Nodes_Of_Minimum_Valence();
        //Highlight_Current_Boundary_Triangle();

        if(current_selection){
            if(current_selection->type==OPENGL_SELECTION::TETRAHEDRALIZED_VOLUME_VERTEX){
                int index=((OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_VERTEX<T>*)current_selection)->index;
                OPENGL_SELECTION::Draw_Highlighted_Vertex(particles->X(index),index);}
            else if(current_selection->type==OPENGL_SELECTION::TETRAHEDRALIZED_VOLUME_TETRAHEDRON){
                int index=((OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_TETRAHEDRON<T>*)current_selection)->index;
                const VECTOR<int,4>& element_nodes=mesh->elements(index);
                ARRAY_VIEW<const TV> X(particles->X);
                OPENGL_SELECTION::Draw_Highlighted_Tetrahedron_Boundary(X(element_nodes[1]),X(element_nodes[2]),X(element_nodes[3]),X(element_nodes[4]),index);
                T distance;TV min_normal,weights;
                int spring=Find_Shortest_Spring(element_nodes,distance,min_normal,weights);
                VECTOR<int,4> spring_nodes;
                spring_nodes=Spring_Nodes(spring,element_nodes);

                glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
                glDisable(GL_LIGHTING);
                glLineWidth(OPENGL_PREFERENCES::highlighted_line_width*2);
                OPENGL_COLOR colors[]={OPENGL_COLOR::Red(),OPENGL_COLOR::Blue(),OPENGL_COLOR::Green(),OPENGL_COLOR::Yellow(),OPENGL_COLOR::Cyan(),OPENGL_COLOR::Magenta(),
                                       OPENGL_COLOR::White()};
                if(spring>0){
                    colors[spring-1].Send_To_GL_Pipeline();
                    OpenGL_Begin(GL_LINES);
                    if(spring<=4) OpenGL_Line(X(spring_nodes[1]),TRIANGLE_3D<T>(X.Subset(spring_nodes.Remove_Index(1))).Point_From_Barycentric_Coordinates(weights));
                    else if(spring<=7) OpenGL_Line((1-weights.x)*X(spring_nodes[1])+weights.x*X(spring_nodes[2]),(1-weights.y)*X(spring_nodes[3])+weights.y*X(spring_nodes[4]));
                    OpenGL_End();
                    glPopAttrib();}}}}
    glPopMatrix();
}
//#####################################################################
// Function Find_Shortest_Spring
//#####################################################################
template<class T> int OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Find_Shortest_Spring(const VECTOR<int,4>& element_nodes,T& minimum_signed_distance,TV& minimum_normal,TV& weights) const
{
    ARRAY_VIEW<const TV> X(particles->X);
    // Find Shortest Altitude
    int maximum_cross_squared_index=0;T maximum_cross_squared=(T)-FLT_MAX;TV maximum_cross;
    for(unsigned char h=1;h<=4;h++){
        VECTOR<int,4> spring_nodes=Spring_Nodes(h,element_nodes);
        TV u_cross_v=TV::Cross_Product(X(spring_nodes[3])-X(spring_nodes[2]),X(spring_nodes[4])-X(spring_nodes[2]));
        T u_cross_v_squared=u_cross_v.Magnitude_Squared();
        if(u_cross_v_squared>maximum_cross_squared){maximum_cross_squared_index=h;maximum_cross_squared=u_cross_v_squared;maximum_cross=u_cross_v;}}
    for(unsigned char h=5;h<=7;h++){
        VECTOR<int,4> spring_nodes=Spring_Nodes(h,element_nodes);
        TV u_cross_v=TV::Cross_Product(X(spring_nodes[2])-X(spring_nodes[1]),X(spring_nodes[4])-X(spring_nodes[3]));
        T u_cross_v_squared=u_cross_v.Magnitude_Squared();
        if(u_cross_v_squared>maximum_cross_squared){maximum_cross_squared_index=h;maximum_cross_squared=u_cross_v_squared;maximum_cross=u_cross_v;}}
    TV u,v;VECTOR<int,4> spring_nodes=Spring_Nodes(maximum_cross_squared_index,element_nodes);
    if(maximum_cross_squared_index<5){u=X(spring_nodes[3])-X(spring_nodes[2]);v=X(spring_nodes[4])-X(spring_nodes[2]);}
    else{u=X(spring_nodes[2])-X(spring_nodes[1]);v=X(spring_nodes[4])-X(spring_nodes[3]);}
    T u_length_squared=u.Magnitude_Squared(),v_length_squared=v.Magnitude_Squared();
    if(abs(maximum_cross_squared)<sqr(sin((T)pi/(T)180))*u_length_squared*v_length_squared){
        VECTOR<int,2> edge_index;T max_distance_squared=0;
        for(int i=1;i<=3;i++) for(int j=i+1;j<=4;j++){
            T distance_squared=(X(element_nodes[i])-X(element_nodes[j])).Magnitude_Squared();
            if(distance_squared>max_distance_squared){edge_index=VECTOR<int,2>(i,j);max_distance_squared=distance_squared;}}
            VECTOR<int,2> other_nodes=element_nodes.Remove_Index(edge_index[2]).Remove_Index(edge_index[1]);
            if((X(element_nodes[edge_index[1]])-X(other_nodes[2])).Magnitude_Squared()<(X(element_nodes[edge_index[1]])-X(other_nodes[1])).Magnitude_Squared())
                other_nodes=other_nodes.Reversed();
            VECTOR<int,2> edge1(element_nodes[edge_index[1]],other_nodes[2]),edge2(other_nodes[1],element_nodes[edge_index[2]]);
            bool found=false;VECTOR<int,4> spring_nodes;
            for(maximum_cross_squared_index=5;maximum_cross_squared_index<=7;maximum_cross_squared_index++){
                spring_nodes=Spring_Nodes(maximum_cross_squared_index,element_nodes);
                if(spring_nodes.Slice<1,2>()==edge1 || spring_nodes.Slice<3,4>()==edge1 || spring_nodes.Slice<1,2>().Reversed()==edge1 || spring_nodes.Slice<3,4>().Reversed()==edge1){
                    found=true;break;}}
            PHYSBAM_ASSERT(found);
            minimum_normal=(X(edge1[1])-X(edge2[2])).Orthogonal_Vector().Normalized();
            minimum_signed_distance=0;
            TV midpoint=(T).5*(X(edge1[2])+X(edge2[1]));
            weights=VECTOR<T,3>(SEGMENT_3D<T>(X(spring_nodes[1]),X(spring_nodes[2])).Interpolation_Fraction(midpoint),
                SEGMENT_3D<T>(X(spring_nodes[3]),X(spring_nodes[4])).Interpolation_Fraction(midpoint),0);}
    else{
        minimum_normal=maximum_cross.Normalized();
        minimum_signed_distance=TV::Dot_Product(minimum_normal,X(spring_nodes[1])-X(spring_nodes[3]));
        if(maximum_cross_squared_index<5){
            weights=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(X(spring_nodes[1]),X(spring_nodes[2]),X(spring_nodes[3]),X(spring_nodes[4]));}
        else if(maximum_cross_squared_index<8){
            VECTOR<T,2> dummy_weights;
            SEGMENT_3D<T>(X(spring_nodes[1]),X(spring_nodes[2])).Shortest_Vector_Between_Segments(SEGMENT_3D<T>(X(spring_nodes[3]),X(spring_nodes[4])),dummy_weights);
            weights=dummy_weights.Append(0);}}
    return maximum_cross_squared_index;

}
//#####################################################################
// Function Highlight_Nodes_Of_Minimum_Valence
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Highlight_Nodes_Of_Minimum_Valence() const
{
    bool neighbor_nodes_defined=mesh->neighbor_nodes!=0;if(!neighbor_nodes_defined) mesh->Initialize_Neighbor_Nodes();
    for(int p=1;p<=particles->array_collection->Size();p++) if((*mesh->neighbor_nodes)(p).m==minimum_valence)
        OPENGL_SHAPES::Draw_Dot(particles->X(p),OPENGL_COLOR(1,1,0),4);
    if(!neighbor_nodes_defined){delete mesh->neighbor_nodes;mesh->neighbor_nodes=0;}
}
//#####################################################################
// Function Highlight_Current_Node
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Highlight_Current_Node() const
{
    OPENGL_SHAPES::Draw_Dot(particles->X(current_node),OPENGL_COLOR(1,0,1),7);
}
//#####################################################################
// Function Highlight_Boundary_Normal_Vectors of Current_Tetrahedron
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Highlight_Boundary_Normal_Vectors_Of_Current_Tetrahedron() const
{
    int t=current_tetrahedron,i,j,k,l;mesh->elements(t).Get(i,j,k,l);
    VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k),xl=particles->X(l);
    OPENGL_SHAPES::Draw_Dot((T)one_third*(xj+xl+xk),OPENGL_COLOR((float).7,0,0),5);
    OPENGL_SHAPES::Draw_Vector((T)one_third*(xj+xl+xk),TRIANGLE_3D<T>::Normal(xj,xl,xk),OPENGL_COLOR(1,0,1),2);
    OPENGL_SHAPES::Draw_Vector((T)one_third*(xi+xk+xl),TRIANGLE_3D<T>::Normal(xi,xk,xl),OPENGL_COLOR(1,0,1),1);
    OPENGL_SHAPES::Draw_Vector((T)one_third*(xi+xl+xj),TRIANGLE_3D<T>::Normal(xi,xl,xj),OPENGL_COLOR(1,0,1),1);
    OPENGL_SHAPES::Draw_Vector((T)one_third*(xi+xj+xk),TRIANGLE_3D<T>::Normal(xi,xj,xk),OPENGL_COLOR(1,0,1),1);
}
//#####################################################################
// Function Highlight_Boundary_Nodes_Of_Current_Tetrahedron
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Highlight_Boundary_Nodes_Of_Current_Tetrahedron() const
{
    int t=current_tetrahedron,i,j,k,l;mesh->elements(t).Get(i,j,k,l);
    VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k),xl=particles->X(l);
    if((*mesh->node_on_boundary)(i)) OPENGL_SHAPES::Draw_Dot(xi,OPENGL_COLOR(1,0,0));
    if((*mesh->node_on_boundary)(j)) OPENGL_SHAPES::Draw_Dot(xj,OPENGL_COLOR(1,0,0));
    if((*mesh->node_on_boundary)(k)) OPENGL_SHAPES::Draw_Dot(xk,OPENGL_COLOR(1,0,0));
    if((*mesh->node_on_boundary)(l)) OPENGL_SHAPES::Draw_Dot(xl,OPENGL_COLOR(1,0,0));
}
//#####################################################################
// Function Draw_Wireframe_Mesh
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Wireframe_Mesh(const TETRAHEDRON_MESH& tetrahedron_mesh) const
{
    glDisable(GL_LIGHTING);
    material.diffuse.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_LINES);
    for(int t=1;t<=tetrahedron_mesh.elements.m;t++){
        int i,j,k,l;tetrahedron_mesh.elements(t).Get(i,j,k,l);
        OpenGL_Line(particles->X(i),particles->X(j));
        OpenGL_Line(particles->X(j),particles->X(k));
        OpenGL_Line(particles->X(k),particles->X(l));
        OpenGL_Line(particles->X(l),particles->X(i));
        OpenGL_Line(particles->X(i),particles->X(k));
        OpenGL_Line(particles->X(j),particles->X(l));}
    OpenGL_End();
    glEnable(GL_LIGHTING);
}
//#####################################################################
// Function Draw_Current_Tetrahedron
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Current_Tetrahedron() const
{
    int t=current_tetrahedron,i,j,k,l;mesh->elements(t).Get(i,j,k,l);
    VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k),xl=particles->X(l);
    TETRAHEDRON<T> tet(xi,xj,xk,xl);
    OpenGL_Begin(GL_TRIANGLES);
    OpenGL_Normal(tet.triangle1.normal);OpenGL_Triangle(tet.triangle1.x1,tet.triangle1.x2,tet.triangle1.x3);
    OpenGL_Normal(tet.triangle2.normal);OpenGL_Triangle(tet.triangle2.x1,tet.triangle2.x2,tet.triangle2.x3);
    OpenGL_Normal(tet.triangle3.normal);OpenGL_Triangle(tet.triangle3.x1,tet.triangle3.x2,tet.triangle3.x3);
    OpenGL_Normal(tet.triangle4.normal);OpenGL_Triangle(tet.triangle4.x1,tet.triangle4.x2,tet.triangle4.x3);
    OpenGL_End();
}
//#####################################################################
// Function Set_Color_From_Aspect_Ratio
//#####################################################################
template<class T>
void Set_Color_From_Aspect_Ratio(const TRIANGLE_3D<T>& triangle)
{
    if(abs(triangle.Aspect_Ratio()-root_two)<1e-3) OPENGL_MATERIAL::Plastic(OPENGL_COLOR(1.f,.1f,.3f)).Send_To_GL_Pipeline();
    else OPENGL_MATERIAL::Plastic(OPENGL_COLOR(.5f,1.f,.3f)).Send_To_GL_Pipeline();
}
//#####################################################################
// Function Draw_Subset
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Subset() const
{
    OPENGL_MATERIAL::Plastic(OPENGL_COLOR(VECTOR<double,3>((T)1,(T).7,(T).7))).Send_To_GL_Pipeline();
    for(int tet_index=1;tet_index<=subset.m;tet_index++){
        int t=subset(tet_index);int i,j,k,l;mesh->elements(t).Get(i,j,k,l);
        VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k),xl=particles->X(l);
        TETRAHEDRON<T> tet(xi,xj,xk,xl);
        OpenGL_Begin(GL_TRIANGLES);
        //Set_Color_From_Aspect_Ratio(tet.triangle1);
        OpenGL_Normal(tet.triangle1.normal);OpenGL_Triangle(tet.triangle1.x1,tet.triangle1.x2,tet.triangle1.x3);
        //Set_Color_From_Aspect_Ratio(tet.triangle2);
        OpenGL_Normal(tet.triangle2.normal);OpenGL_Triangle(tet.triangle2.x1,tet.triangle2.x2,tet.triangle2.x3);
        //Set_Color_From_Aspect_Ratio(tet.triangle3);
        //OPENGL_MATERIAL::Plastic(OPENGL_COLOR(VECTOR<T,3>(T(1),T(.7),T(.7)))).Send_To_GL_Pipeline();
        OpenGL_Normal(tet.triangle3.normal);OpenGL_Triangle(tet.triangle3.x1,tet.triangle3.x2,tet.triangle3.x3);
        //Set_Color_From_Aspect_Ratio(tet.triangle4);
        OpenGL_Normal(tet.triangle4.normal);OpenGL_Triangle(tet.triangle4.x1,tet.triangle4.x2,tet.triangle4.x3);
        OpenGL_End();}
    material.Send_To_GL_Pipeline();
}
//#####################################################################
// Function Draw_Subset_Triangles
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Subset_Triangles() const
{
    OPENGL_MATERIAL::Plastic(OPENGL_COLOR(VECTOR<double,3>((T).7,(T).7,(T)1))).Send_To_GL_Pipeline();
    OpenGL_Begin(GL_TRIANGLES);
    for(int t=1;t<=subset_triangles.m;t++){
        int i,j,k,tri;tri=subset_triangles(t);
        mesh->boundary_mesh->elements(tri).Get(i,j,k);
        VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k);
        OpenGL_Normal(TRIANGLE_3D<T>::Normal(xi,xj,xk));
        OpenGL_Triangle(xi,xj,xk);}
    OpenGL_End();
    material.Send_To_GL_Pipeline();
}
//#####################################################################
// Function Draw_Subset_Particles
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Subset_Particles() const
{
    for(int p=1;p<=subset_particles.m;p++) OPENGL_SHAPES::Draw_Dot(particles->X(subset_particles(p)),OPENGL_COLOR(1,1,0));
}
//#####################################################################
// Function Set_Color_From_Spectrum
//#####################################################################
template<class T>
void Set_Color_From_Spectrum(int tetrahedron,const ARRAY<T>& spectrum,T spectrum_max,T spectrum_min)
{
    T red=spectrum(tetrahedron)-spectrum_min/(spectrum_max-spectrum_min);
    PHYSBAM_ASSERT(0<=red && red<=1);
    OPENGL_MATERIAL::Plastic(OPENGL_COLOR((float)red,.1f,.1f)).Send_To_GL_Pipeline();
}
//#####################################################################
// Function Draw_In_Color_From_Spectrum
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_In_Color_From_Spectrum() const
{
    if(spectrum.m!=mesh->elements.m) PHYSBAM_FATAL_ERROR("Spectrum has incorrect size");
    T spectrum_max=ARRAYS_COMPUTATIONS::Max(spectrum),spectrum_min=ARRAYS_COMPUTATIONS::Min(spectrum);
    for(int t=1;t<=mesh->elements.m;t++){
        int i,j,k,l;mesh->elements(t).Get(i,j,k,l);
        VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k),xl=particles->X(l);
        TETRAHEDRON<T> tet(xi,xj,xk,xl);
        OpenGL_Begin(GL_TRIANGLES);
        Set_Color_From_Spectrum(t,spectrum,spectrum_max,spectrum_min);
        OpenGL_Normal(tet.triangle1.normal);OpenGL_Triangle(tet.triangle1.x1,tet.triangle1.x2,tet.triangle1.x3);
        Set_Color_From_Spectrum(t,spectrum,spectrum_max,spectrum_min);
        OpenGL_Normal(tet.triangle2.normal);OpenGL_Triangle(tet.triangle2.x1,tet.triangle2.x2,tet.triangle2.x3);
        Set_Color_From_Spectrum(t,spectrum,spectrum_max,spectrum_min);
        OpenGL_Normal(tet.triangle3.normal);OpenGL_Triangle(tet.triangle3.x1,tet.triangle3.x2,tet.triangle3.x3);
        Set_Color_From_Spectrum(t,spectrum,spectrum_max,spectrum_min);
        OpenGL_Normal(tet.triangle4.normal);OpenGL_Triangle(tet.triangle4.x1,tet.triangle4.x2,tet.triangle4.x3);
        OpenGL_End();}
    material.Send_To_GL_Pipeline();
}
//#####################################################################
// Function Draw_In_Color_From_Spectrum
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_In_Color_From_Color_Map() const
{
    glPushAttrib(GL_ENABLE_BIT|GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    for(int t=1;t<=mesh->elements.m;t++){
        (*color_map)(t).Send_To_GL_Pipeline();
        int i,j,k,l;mesh->elements(t).Get(i,j,k,l);
        VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k),xl=particles->X(l);
        TETRAHEDRON<T> tet(xi,xj,xk,xl);
        OpenGL_Begin(GL_TRIANGLES);
        OpenGL_Normal(tet.triangle1.normal);OpenGL_Triangle(tet.triangle1.x1,tet.triangle1.x2,tet.triangle1.x3);
        OpenGL_Normal(tet.triangle2.normal);OpenGL_Triangle(tet.triangle2.x1,tet.triangle2.x2,tet.triangle2.x3);
        OpenGL_Normal(tet.triangle3.normal);OpenGL_Triangle(tet.triangle3.x1,tet.triangle3.x2,tet.triangle3.x3);
        OpenGL_Normal(tet.triangle4.normal);OpenGL_Triangle(tet.triangle4.x1,tet.triangle4.x2,tet.triangle4.x3);
        OpenGL_End();}
    glPopAttrib();
}
//#####################################################################
// Function Draw_Boundary_Triangles
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Boundary_Triangles(const TETRAHEDRON_MESH& tetrahedron_mesh) const
{
    material.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_TRIANGLES);
    if(smooth_normals)
        for(int t=1;t<=tetrahedron_mesh.boundary_mesh->elements.m;t++){
            int i,j,k;tetrahedron_mesh.boundary_mesh->elements(t).Get(i,j,k);
            OpenGL_Normal((*vertex_normals)(i));OpenGL_Vertex(particles->X(i));
            OpenGL_Normal((*vertex_normals)(j));OpenGL_Vertex(particles->X(j));
            OpenGL_Normal((*vertex_normals)(k));OpenGL_Vertex(particles->X(k));}
    else
        for(int t=1;t<=tetrahedron_mesh.boundary_mesh->elements.m;t++){
            int i,j,k;tetrahedron_mesh.boundary_mesh->elements(t).Get(i,j,k);
            VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k);
            OpenGL_Normal(TRIANGLE_3D<T>::Normal(xi,xj,xk));
            OpenGL_Triangle(xi,xj,xk);}
    OpenGL_End();
}
//#####################################################################
// Function Highlight_Current_Boundary_Triangle
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Highlight_Current_Boundary_Triangle() const
{
    glDisable(GL_LIGHTING);
    glColor3f(1,1,1);
    OpenGL_Begin(GL_TRIANGLES);
    int i,j,k;mesh->boundary_mesh->elements(current_boundary_triangle).Get(i,j,k);
    VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k);
    OpenGL_Normal(TRIANGLE_3D<T>::Normal(xi,xj,xk));
    OpenGL_Triangle(xi,xj,xk);
    OpenGL_End();
    OPENGL_SHAPES::Draw_Dot((T)1./3*(xi+xj+xk),OPENGL_COLOR(1,1,1),8);
    glEnable(GL_LIGHTING);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<VECTOR<T,3> >::Bounding_Box(particles->X));
}
//#####################################################################
// Function Turn_Smooth_Shading_On
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Turn_Smooth_Shading_On()
{
    if(!vertex_normals)Initialize_Vertex_Normals();smooth_normals=true;
}
//#####################################################################
// Function Turn_Smooth_Shading_Off
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Turn_Smooth_Shading_Off()
{
    smooth_normals=false;
}
//#####################################################################
// Function Initialize_Vertex_Normals
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Initialize_Vertex_Normals()
{
    delete vertex_normals;vertex_normals=new ARRAY<VECTOR<T,3> >(particles->array_collection->Size());
    for(int t=1;t<=mesh->boundary_mesh->elements.m;t++){
        int i,j,k;mesh->boundary_mesh->elements(t).Get(i,j,k);
        VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k);
        VECTOR<T,3> normal=TRIANGLE_3D<T>::Normal(xi,xj,xk);
        (*vertex_normals)(i)+=normal;(*vertex_normals)(j)+=normal;(*vertex_normals)(k)+=normal;}
    for(int p=1;p<=particles->array_collection->Size();p++)(*vertex_normals)(p).Normalize();
}
//#####################################################################
// Function Update_Cutaway_Plane
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Update_Cutaway_Plane()
{
    TETRAHEDRALIZED_VOLUME<T> tetrahedralized_volume(*mesh,*particles);
    tetrahedralized_volume.Update_Bounding_Box();
    BOX<VECTOR<T,3> > box=*tetrahedralized_volume.bounding_box;
    ARRAY<bool> inside(particles->array_collection->Size());T threshold;
    switch(cutaway_mode){
        case 1:threshold=box.min_corner.x+cutaway_fraction*(box.max_corner.x-box.min_corner.x);for(int p=1;p<=particles->array_collection->Size();p++)inside(p)=particles->X(p).x<threshold;break;
        case 2:threshold=box.max_corner.x+cutaway_fraction*(box.min_corner.x-box.max_corner.x);for(int p=1;p<=particles->array_collection->Size();p++)inside(p)=particles->X(p).x>threshold;break;
        case 3:threshold=box.min_corner.y+cutaway_fraction*(box.max_corner.y-box.min_corner.y);for(int p=1;p<=particles->array_collection->Size();p++)inside(p)=particles->X(p).y<threshold;break;
        case 4:threshold=box.max_corner.y+cutaway_fraction*(box.min_corner.y-box.max_corner.y);for(int p=1;p<=particles->array_collection->Size();p++)inside(p)=particles->X(p).y>threshold;break;
        case 5:threshold=box.min_corner.z+cutaway_fraction*(box.max_corner.z-box.min_corner.z);for(int p=1;p<=particles->array_collection->Size();p++)inside(p)=particles->X(p).z<threshold;break;
        case 6:threshold=box.max_corner.z+cutaway_fraction*(box.min_corner.z-box.max_corner.z);for(int p=1;p<=particles->array_collection->Size();p++)inside(p)=particles->X(p).z>threshold;break;}
    cutaway_mesh.number_nodes=mesh->number_nodes;
    cutaway_mesh.elements.Remove_All();
    for(int t=1;t<=mesh->elements.m;t++){
        int i,j,k,l;mesh->elements(t).Get(i,j,k,l);
        if(inside(i) || inside(j) || inside(k) || inside(l)) cutaway_mesh.elements.Append(VECTOR<int,4>(i,j,k,l));}
    cutaway_mesh.Initialize_Boundary_Mesh();
}
//#####################################################################
// Function Display_Subset
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Display_Subset()
{
    material.Send_To_GL_Pipeline();
    glColor3f((float).7,(float).7,(float).7);
    OpenGL_Begin(GL_TRIANGLES);
    for(int tet_index=1; tet_index<=subset.m; tet_index++){
        int t=subset(tet_index);
        int i,j,k,l;mesh->elements(t).Get(i,j,k,l);
        VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k),xl=particles->X(l);
        TETRAHEDRON<T> tet(xi,xj,xk,xl);
        OpenGL_Normal(tet.triangle1.normal);OpenGL_Triangle(tet.triangle1.x1,tet.triangle1.x2,tet.triangle1.x3);
        OpenGL_Normal(tet.triangle2.normal);OpenGL_Triangle(tet.triangle2.x1,tet.triangle2.x2,tet.triangle2.x3);
        OpenGL_Normal(tet.triangle3.normal);OpenGL_Triangle(tet.triangle3.x1,tet.triangle3.x2,tet.triangle3.x3);
        OpenGL_Normal(tet.triangle4.normal);OpenGL_Triangle(tet.triangle4.x1,tet.triangle4.x2,tet.triangle4.x3);}
    OpenGL_End();
    OpenGL_Begin(GL_LINES);
    for(int t=1;t<=mesh->boundary_mesh->elements.m;t++){
        int i,j,k;mesh->boundary_mesh->elements(t).Get(i,j,k);
        VECTOR<T,3> xi=particles->X(i),xj=particles->X(j),xk=particles->X(k);
        OpenGL_Line(xi,xj);
        OpenGL_Line(xj,xk);
        OpenGL_Line(xi,xk);}
    OpenGL_End();
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION* OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Get_Selection(GLuint* buffer,int buffer_size)
{
    OPENGL_SELECTION* selection=0;
    if(buffer_size==2){
        if(buffer[0]==1)selection=Get_Vertex_Selection(buffer[1]);
        else if(buffer[0]==2)selection=Get_Tetrahedron_Selection(buffer[1]);}
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Highlight_Selection(OPENGL_SELECTION* selection)
{
    delete current_selection;current_selection=0;
    // Make a copy of selection
    if(selection->type==OPENGL_SELECTION::TETRAHEDRALIZED_VOLUME_VERTEX)
        current_selection=new OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_VERTEX<T>(this,((OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_VERTEX<T>*)selection)->index);
    else if(selection->type==OPENGL_SELECTION::TETRAHEDRALIZED_VOLUME_TETRAHEDRON)
        current_selection=new OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_TETRAHEDRON<T>(this,((OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_TETRAHEDRON<T>*)selection)->index);
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Clear_Highlight()
{
    delete current_selection;current_selection=0;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION *selection) const
{
    Print_Selection_Info(output_stream,selection,0);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION* selection,MATRIX<T,4>* transform) const
{
    if(selection->type==OPENGL_SELECTION::TETRAHEDRALIZED_VOLUME_VERTEX){
        int index=((OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_VERTEX<T>*)selection)->index;
        output_stream<<"Vertex "<<index<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(particles->X(index))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,3> >,T>::Print(output_stream,*particles,index);}
    else if(selection->type==OPENGL_SELECTION::TETRAHEDRALIZED_VOLUME_TETRAHEDRON){
        int index=((OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_TETRAHEDRON<T>*)selection)->index;
        const VECTOR<int,4>& nodes=mesh->elements(index);
        output_stream<<"Tetrahedron "<<index<<" ("<<nodes[1]<<", "<<nodes[2]<<", "<<nodes[3]<<", "<<nodes[4]<<")"<<std::endl;
        output_stream<<"Signed Volume = "<<TETRAHEDRON<T>(particles->X.Subset(nodes)).Signed_Volume()<<std::endl;
        output_stream<<std::endl;
        output_stream<<"Vertex"<<nodes[1]<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(particles->X(nodes[1]))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,3> >,T>::Print(output_stream,*particles,nodes[1]);
        output_stream<<std::endl;
        output_stream<<"Vertex "<<nodes[2]<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(particles->X(nodes[2]))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,3> >,T>::Print(output_stream,*particles,nodes[2]);
        output_stream<<std::endl;
        output_stream<<"Vertex "<<nodes[3]<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(particles->X(nodes[3]))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,3> >,T>::Print(output_stream,*particles,nodes[3]);
        output_stream<<std::endl;
        output_stream<<"Vertex "<<nodes[4]<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(particles->X(nodes[4]))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,3> >,T>::Print(output_stream,*particles,nodes[4]);

        T distance;TV min_normal,weights;
        int spring=Find_Shortest_Spring(nodes,distance,min_normal,weights);
        output_stream<<"Shortest Altitude is ";
        if(spring<=4) output_stream<<"point-face ";
        else if(spring<=7) output_stream<<"edge-edge ";
        output_stream<<spring<<" nodes="<<Spring_Nodes(spring,nodes)<<std::endl;
        output_stream<<"   distance="<<distance<<" normal="<<min_normal<<" weights="<<weights<<std::endl;
    }

}
//#####################################################################
// Function Get_Vertex_Selection
//#####################################################################
template<class T> OPENGL_SELECTION *OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Get_Vertex_Selection(int index)
{
    return new OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_VERTEX<T>(this,index);
}
//#####################################################################
// Function Get_Tetrahedron_Selection
//#####################################################################
template<class T> OPENGL_SELECTION *OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Get_Tetrahedron_Selection(int index)
{
    return new OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_TETRAHEDRON<T>(this,index);
}
//#####################################################################
// Function Draw_Vertices_For_Selection
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Vertices_For_Selection() const
{
    OPENGL_SELECTION::Draw_Vertices_For_Selection(*mesh,*particles);
}
//#####################################################################
// Function Draw_Tetrahedra_For_Selection
//#####################################################################
template<class T> void OPENGL_TETRAHEDRALIZED_VOLUME<T>::
Draw_Tetrahedra_For_Selection() const
{
    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_CULL_FACE);
    glPushName(0);
    for(int t=1;t<=mesh->elements.m;t++){
        glLoadName(t);OpenGL_Begin(GL_TRIANGLES);
        int i,j,k,l;mesh->elements(t).Get(i,j,k,l);
        OpenGL_Triangle(particles->X(i),particles->X(j),particles->X(k));
        OpenGL_Triangle(particles->X(i),particles->X(k),particles->X(l));
        OpenGL_Triangle(particles->X(i),particles->X(l),particles->X(j));
        OpenGL_Triangle(particles->X(l),particles->X(k),particles->X(j));
        OpenGL_End();}
    glPopName();
    glPopAttrib();
}
//#####################################################################
// Selection object functions
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_VERTEX<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const OPENGL_TETRAHEDRALIZED_VOLUME<T>& volume=dynamic_cast<OPENGL_TETRAHEDRALIZED_VOLUME<T>&>(*object);
    return object->World_Space_Box(RANGE<VECTOR<float,3> >(VECTOR<float,3>(volume.particles->X(index))));
}
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_TETRAHEDRALIZED_VOLUME_TETRAHEDRON<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const OPENGL_TETRAHEDRALIZED_VOLUME<T>& volume=dynamic_cast<OPENGL_TETRAHEDRALIZED_VOLUME<T>&>(*object);
    return object->World_Space_Box(RANGE<VECTOR<float,3> >(RANGE<VECTOR<T,3> >::Bounding_Box(volume.particles->X.Subset(volume.mesh->elements(index)))));
}
//#####################################################################
template class OPENGL_TETRAHEDRALIZED_VOLUME<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_TETRAHEDRALIZED_VOLUME<double>;
#endif

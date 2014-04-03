//#####################################################################
// Copyright 2002-2005, Eran Guendelman, Eilene Hao, Neil Molino, Robert Bridson, Igor Neverov, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor from geometry, one-sided
//#####################################################################
template<class T> OPENGL_TRIANGULATED_SURFACE<T>::
OPENGL_TRIANGULATED_SURFACE(TRIANGULATED_SURFACE<T>& surface_input,bool smooth_normals_input,const OPENGL_MATERIAL& material_input)
    :surface(surface_input),two_sided(false),front_material(material_input),front_material_gray(material_input.Grayscale()),
    back_material(material_input),back_material_gray(material_input.Grayscale()),vertex_normals(0),vertex_colors(0),smooth_normals(smooth_normals_input),use_display_list(false),
    owns_display_list(false),current_selection(0),current_node(1),highlight_current_node(false),highlight_neighbors_of_current_node(true),highlight_boundary(false),
    wireframe_only(false),draw_subsets(false),draw_velocities(false),draw_particles(false),velocity_scale((T).025)
{
    if(smooth_normals) Initialize_Vertex_Normals();
}
//#####################################################################
// Constructor from geometry, two-sided
//#####################################################################
template<class T> OPENGL_TRIANGULATED_SURFACE<T>::
OPENGL_TRIANGULATED_SURFACE(TRIANGULATED_SURFACE<T>& surface_input,bool smooth_normals_input,const OPENGL_MATERIAL& front_material_input,const OPENGL_MATERIAL& back_material_input)
    :surface(surface_input),two_sided(true),front_material(front_material_input),front_material_gray(front_material_input.Grayscale()),
    back_material(back_material_input),back_material_gray(back_material_input.Grayscale()),vertex_normals(0),vertex_colors(0),smooth_normals(smooth_normals_input),use_display_list(false),
    owns_display_list(false),current_selection(0),current_node(1),highlight_current_node(false),highlight_neighbors_of_current_node(true),highlight_boundary(false),wireframe_only(false),
    draw_subsets(false),draw_velocities(false),draw_particles(false),velocity_scale((T).025)
{
    if(smooth_normals) Initialize_Vertex_Normals();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_TRIANGULATED_SURFACE<T>::
~OPENGL_TRIANGULATED_SURFACE()
{
    if(owns_display_list) glDeleteLists(display_list_id,1);
    delete vertex_normals;delete vertex_colors;delete current_selection;
}
//#####################################################################
// Function Highlight_Current_Node
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Highlight_Current_Node() const
{
    int node=(current_node-1)%surface.particles.array_collection->Size()+1;
    OPENGL_SHAPES::Draw_Dot(surface.particles.X(node),OPENGL_COLOR(1,0,1),7);
    if(highlight_neighbors_of_current_node){
        if(!surface.mesh.neighbor_nodes) surface.mesh.Initialize_Neighbor_Nodes();
        for(int t=1;t<=(*surface.mesh.neighbor_nodes)(node).m;t++){
            int nbr=(*surface.mesh.neighbor_nodes)(node)(t);
            OPENGL_SHAPES::Draw_Dot(surface.particles.X(nbr),OPENGL_COLOR(1,1,0),5);}}
}
//#####################################################################
// Function Print_Triangles_Incident_On_Current_Node
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Print_Triangles_Incident_On_Current_Node()
{
    int node=(current_node-1)%surface.particles.array_collection->Size()+1;
    if(!surface.mesh.incident_elements) surface.mesh.Initialize_Incident_Elements();
    LOG::cout<<"number of incident triangles to node "<<node<<"="<<(*surface.mesh.incident_elements)(node).m<<std::endl;
    for(int t=1;t<=(*surface.mesh.incident_elements)(node).m;t++){
        int triangle=(*surface.mesh.incident_elements)(node)(t);
        int i,j,k;surface.mesh.elements(triangle).Get(i,j,k);
        LOG::cout<<"triangle "<<triangle<<"=("<<i<<","<<j<<","<<k<<")"<<std::endl;}
}
//#####################################################################
// Function Print_Neighbor_Nodes_Of_Current_Node
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Print_Neighbor_Nodes_Of_Current_Node()
{
    int node=(current_node-1)%surface.particles.array_collection->Size()+1;
    if(!surface.mesh.neighbor_nodes) surface.mesh.Initialize_Neighbor_Nodes();
    LOG::cout<<"number of neighbors of node "<<node<<"="<<(*surface.mesh.neighbor_nodes)(node).m<<std::endl;
    LOG::cout<<"neighbors of node "<<node<<"={";
    for(int t=1;t<=(*surface.mesh.neighbor_nodes)(node).m;t++) LOG::cout<<(*surface.mesh.neighbor_nodes)(node)(t)<<",";
    LOG::cout<<"}"<<std::endl;
}
//#####################################################################
// Function Draw_Triangles_Incident_On_Current_Node
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Draw_Triangles_Incident_On_Current_Node() const
{
    int node=(current_node-1)%surface.particles.array_collection->Size()+1;
    if(!surface.mesh.incident_elements) surface.mesh.Initialize_Incident_Elements();
    glDisable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
    front_material.Send_To_GL_Pipeline();
    back_material.Send_To_GL_Pipeline();
    //glColor3f(1.0,0,0);
    for(int t=1;t<=(*surface.mesh.incident_elements)(node).m;t++){
        int triangle=(*surface.mesh.incident_elements)(node)(t);
        int i,j,k;surface.mesh.elements(triangle).Get(i,j,k);
        OpenGL_Begin(GL_TRIANGLES);
            OpenGL_Normal(PLANE<T>(surface.particles.X(i),surface.particles.X(j),surface.particles.X(k)).Normal());
            OpenGL_Triangle(surface.particles.X(i),surface.particles.X(j),surface.particles.X(k));
        OpenGL_End();}
}
//#####################################################################
// Function Turn_Smooth_Shading_On
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Turn_Smooth_Shading_On()
{
    if(!vertex_normals) Initialize_Vertex_Normals();
    smooth_normals=true;
    if(owns_display_list) Reinitialize_Display_List();
}
//#####################################################################
// Function Turn_Smooth_Shading_Off
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Turn_Smooth_Shading_Off()
{
    smooth_normals=false;
    if(owns_display_list) Reinitialize_Display_List();
}
//#####################################################################
// Function Initialize_Vertex_Normals
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Initialize_Vertex_Normals()
{
    if(!vertex_normals) vertex_normals=new ARRAY<VECTOR<T,3> >;
    vertex_normals->Resize(surface.particles.array_collection->Size());
    ARRAYS_COMPUTATIONS::Fill(*vertex_normals,VECTOR<T,3>());
    for(int t=1;t<=surface.mesh.elements.m;t++){
        int i,j,k;surface.mesh.elements(t).Get(i,j,k);
        VECTOR<T,3> normal=TRIANGLE_3D<T>::Normal(surface.particles.X(i),surface.particles.X(j),surface.particles.X(k));
        (*vertex_normals)(i)+=normal;(*vertex_normals)(j)+=normal;(*vertex_normals)(k)+=normal;}
    for(int p=1;p<=surface.particles.array_collection->Size();p++)(*vertex_normals)(p).Normalize();
}
//#####################################################################
// Function Delete_Vertex_Normals
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Delete_Vertex_Normals()
{
    delete vertex_normals;vertex_normals=0;
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Display(const int in_color) const
{
    if(two_sided){
        glDisable(GL_CULL_FACE);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
        if(in_color){front_material.Send_To_GL_Pipeline(GL_FRONT);back_material.Send_To_GL_Pipeline(GL_BACK);}
        else{front_material_gray.Send_To_GL_Pipeline(GL_FRONT);back_material_gray.Send_To_GL_Pipeline(GL_BACK);}}
    else{
        if(in_color) front_material.Send_To_GL_Pipeline();
        else front_material_gray.Send_To_GL_Pipeline();}

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    if(draw_particles) for(int i=1;i<=surface.particles.array_collection->Size();i++) OPENGL_SHAPES::Draw_Dot(surface.particles.X(i),OPENGL_COLOR(1,0,1),7);

    GLint mode;
    glGetIntegerv(GL_RENDER_MODE,&mode);

    if(wireframe_only){
        glPushAttrib(GL_POLYGON_BIT);
        glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);}

    if(mode == GL_SELECT){
        glPushName(1);
        Draw_Vertices_For_Selection();
        if(surface.mesh.segment_mesh){
            glLoadName(2);
            Draw_Segments_For_Selection();}
        glLoadName(3);
        Draw_Triangles_For_Selection();
        glPopName();}
    else if(use_display_list) glCallList(display_list_id);
    else{Draw();if(draw_subsets) Draw_Subsets();}
    if(wireframe_only) glPopAttrib();

    if(mode != GL_SELECT){
        if(current_selection){
            if(current_selection->type == OPENGL_SELECTION::TRIANGULATED_SURFACE_VERTEX){
                int index=((OPENGL_SELECTION_TRIANGULATED_SURFACE_VERTEX<T>*)current_selection)->index;
                OPENGL_SELECTION::Draw_Highlighted_Vertex(surface.particles.X(index));}
            else if(current_selection->type == OPENGL_SELECTION::TRIANGULATED_SURFACE_SEGMENT && surface.mesh.segment_mesh){
                int index=((OPENGL_SELECTION_TRIANGULATED_SURFACE_SEGMENT<T>*)current_selection)->index;
                int node1,node2;surface.mesh.segment_mesh->elements(index).Get(node1,node2);
                OPENGL_SELECTION::Draw_Highlighted_Segment(surface.particles.X(node1),surface.particles.X(node2));}
            else if(current_selection->type == OPENGL_SELECTION::TRIANGULATED_SURFACE_TRIANGLE){
                int index=((OPENGL_SELECTION_TRIANGULATED_SURFACE_TRIANGLE<T>*)current_selection)->index;
                int node1,node2,node3;surface.mesh.elements(index).Get(node1,node2,node3);
                OPENGL_SELECTION::Draw_Highlighted_Triangle_Boundary(surface.particles.X(node1),surface.particles.X(node2),surface.particles.X(node3));}}

        if(highlight_current_node){
            Highlight_Current_Node();
            Draw_Triangles_Incident_On_Current_Node();}
        else if(highlight_boundary) {
            glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
            glDisable(GL_LIGHTING);
            glLineWidth(5);
            if(!surface.mesh.boundary_mesh) surface.mesh.Initialize_Boundary_Mesh();
            if(!surface.mesh.boundary_mesh->connected_segments) surface.mesh.boundary_mesh->Initialize_Connected_Segments();
#if 0
            OpenGL_Begin(GL_LINES);
            OPENGL_COLOR::Yellow().Send_To_GL_Pipeline();
            for(int i=1;i<=surface.mesh.boundary_mesh->elements.m;i++){
                int node1,node2;surface.mesh.boundary_mesh->elements(i).Get(node1,node2);
                OpenGL_Line(surface.particles.X(node1),surface.particles.X(node2));}
            OpenGL_End();
#else
            ARRAY<ARRAY<VECTOR<int,2> > >& connected_segments=*surface.mesh.boundary_mesh->connected_segments;
            for(int i=1;i<=connected_segments.m;i++){
                OPENGL_COLOR::Random().Send_To_GL_Pipeline();
                OpenGL_Begin(GL_LINES);
                for(int j=1;j<=connected_segments(i).m;j++){
                    int node1,node2;connected_segments(i)(j).Get(node1,node2);
                    OpenGL_Line(surface.particles.X(node1),surface.particles.X(node2));}
                OpenGL_End();}
#endif
            glPopAttrib();}
        if(draw_velocities){
            glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
            glDisable(GL_LIGHTING);
            glLineWidth(5);
            for(int i=1;i<=surface.mesh.elements.m;i++){
                const VECTOR<int,3>& nodes=surface.mesh.elements(i);
                for(int j=1;j<=nodes.m;j++){
                    //OPENGL_SHAPES::Draw_Arrow(surface.particles.X(nodes[j]),surface.particles.X(nodes[j])+velocity_scale*surface.particles.V(nodes[j]));
                }
            }
            glPopAttrib();
        }
    }

    if(two_sided){
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,0);
        glEnable(GL_CULL_FACE);}

    glPopMatrix();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_TRIANGULATED_SURFACE<T>::
Bounding_Box() const
{
    RANGE<VECTOR<float,3> > box=RANGE<VECTOR<float,3> >::Empty_Box();
    for(int i=1;i<=surface.particles.array_collection->Size();i++) box.Enlarge_To_Include_Point(World_Space_Point(VECTOR<float,3>(surface.particles.X(i))));
    return box;
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION *OPENGL_TRIANGULATED_SURFACE<T>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION* selection=0;
    if(buffer_size == 2){
        if(buffer[0] == 1) selection=Get_Vertex_Selection(buffer[1]);
        else if(buffer[0] == 2) selection=Get_Segment_Selection(buffer[1]);
        else if(buffer[0] == 3) selection=Get_Triangle_Selection(buffer[1]);}
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Highlight_Selection(OPENGL_SELECTION* selection)
{
    delete current_selection;current_selection=0;
    // Make a copy of selection
    if(selection->type == OPENGL_SELECTION::TRIANGULATED_SURFACE_VERTEX)
        current_selection=new OPENGL_SELECTION_TRIANGULATED_SURFACE_VERTEX<T>(this,((OPENGL_SELECTION_TRIANGULATED_SURFACE_VERTEX<T>*)selection)->index);
    else if(selection->type == OPENGL_SELECTION::TRIANGULATED_SURFACE_SEGMENT)
        current_selection=new OPENGL_SELECTION_TRIANGULATED_SURFACE_SEGMENT<T>(this,((OPENGL_SELECTION_TRIANGULATED_SURFACE_SEGMENT<T>*)selection)->index);
    else if(selection->type == OPENGL_SELECTION::TRIANGULATED_SURFACE_TRIANGLE)
        current_selection=new OPENGL_SELECTION_TRIANGULATED_SURFACE_TRIANGLE<T>(this,((OPENGL_SELECTION_TRIANGULATED_SURFACE_TRIANGLE<T>*)selection)->index);
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Clear_Highlight()
{
    delete current_selection;current_selection=0;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const
{
    Print_Selection_Info(output_stream,selection,0);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection,MATRIX<T,4>* transform) const
{
    if(selection->type == OPENGL_SELECTION::TRIANGULATED_SURFACE_VERTEX){
        int index=((OPENGL_SELECTION_TRIANGULATED_SURFACE_VERTEX<T>*)selection)->index;
        output_stream<<"Vertex "<<index<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(surface.particles.X(index))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,3> >,T>::Print(output_stream,surface.particles,index);}
    else if(selection->type == OPENGL_SELECTION::TRIANGULATED_SURFACE_SEGMENT){
        PHYSBAM_ASSERT(surface.mesh.segment_mesh);
        int index=((OPENGL_SELECTION_TRIANGULATED_SURFACE_SEGMENT<T>*)selection)->index;
        int node1,node2;surface.mesh.segment_mesh->elements(index).Get(node1,node2);
        output_stream<<"Segment "<<index<<" ("<<node1<<", "<<node2<<")"<<std::endl;
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node1<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(surface.particles.X(node1))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,3> >,T>::Print(output_stream,surface.particles,node1);
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node2<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(surface.particles.X(node2))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,3> >,T>::Print(output_stream,surface.particles,node2);}
    else if(selection->type == OPENGL_SELECTION::TRIANGULATED_SURFACE_TRIANGLE){
        int index=((OPENGL_SELECTION_TRIANGULATED_SURFACE_TRIANGLE<T>*)selection)->index;
        int node1,node2,node3;surface.mesh.elements(index).Get(node1,node2,node3);
        output_stream<<"Triangle "<<index<<" ("<<node1<<", "<<node2<<", "<<node3<<")"<<std::endl;
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node1<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(surface.particles.X(node1))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,3> >,T>::Print(output_stream,surface.particles,node1);
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node2<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(surface.particles.X(node2))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,3> >,T>::Print(output_stream,surface.particles,node2);
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node3<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(surface.particles.X(node3))<<std::endl;}
        Read_Write<GEOMETRY_PARTICLES<VECTOR<T,3> >,T>::Print(output_stream,surface.particles,node3);}
}
//#####################################################################
// Function Get_Vertex_Selection
//#####################################################################
template<class T> OPENGL_SELECTION* OPENGL_TRIANGULATED_SURFACE<T>::
Get_Vertex_Selection(int index)
{
    return new OPENGL_SELECTION_TRIANGULATED_SURFACE_VERTEX<T>(this,index);
}
//#####################################################################
// Function Get_Segment_Selection
//#####################################################################
template<class T> OPENGL_SELECTION* OPENGL_TRIANGULATED_SURFACE<T>::
Get_Segment_Selection(int index)
{
    return new OPENGL_SELECTION_TRIANGULATED_SURFACE_SEGMENT<T>(this,index);
}
//#####################################################################
// Function Get_Triangle_Selection
//#####################################################################
template<class T> OPENGL_SELECTION* OPENGL_TRIANGULATED_SURFACE<T>::
Get_Triangle_Selection(int index)
{
    return new OPENGL_SELECTION_TRIANGULATED_SURFACE_TRIANGLE<T>(this,index);
}
//#####################################################################
// Function Set_Front_Material
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Set_Front_Material(const OPENGL_MATERIAL& material_input)
{
    front_material=material_input;
    front_material_gray=material_input.Grayscale();
}
//#####################################################################
// Function Set_Back_Material
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Set_Back_Material(const OPENGL_MATERIAL& material_input)
{
    back_material=material_input;
    back_material_gray=material_input.Grayscale();
}
//#####################################################################
// Function
//#####################################################################
template<class T> int OPENGL_TRIANGULATED_SURFACE<T>::
Create_Display_List()
{
    PHYSBAM_ASSERT(!owns_display_list);
    owns_display_list=true;
    use_display_list=true;
    display_list_id=glGenLists(1);
    Reinitialize_Display_List();
    return display_list_id;
}
//#####################################################################
// Function
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Use_Display_List(int input_display_list_id)
{
    PHYSBAM_ASSERT(!owns_display_list);
    use_display_list=true;
    display_list_id=input_display_list_id;
}
//#####################################################################
// Function Use_Vertex_Colors
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Use_Vertex_Colors()
{
    if(!vertex_colors) vertex_colors=new ARRAY<OPENGL_COLOR>(surface.particles.array_collection->Size(),false); // TODO: if using big particle pool not efficient
}
//#####################################################################
// Function Set_Vertex_Color
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Set_Vertex_Color(const int i,const OPENGL_COLOR color)
{
    PHYSBAM_ASSERT(vertex_colors);
    (*vertex_colors)(i)=color;
}
//#####################################################################
// Function Reinitialize_Display_List
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Reinitialize_Display_List()
{
    glNewList(display_list_id,GL_COMPILE);
    Draw();
    glEndList();
}
//#####################################################################
// Function Draw
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Draw() const
{
    const ARRAY<VECTOR<int,3> > &elements=surface.mesh.elements;
    if(!vertex_colors){
        //glPushAttrib(GL_ENABLE_BIT);
        //glEnable(GL_BLEND);glDisable(GL_DEPTH_TEST);glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        OpenGL_Begin(GL_TRIANGLES);
        if(smooth_normals)
            for(int t=1;t<=elements.m;t++){
                int i,j,k;elements(t).Get(i,j,k);
                OpenGL_Normal((*vertex_normals)(i));OpenGL_Vertex(surface.particles.X(i));
                OpenGL_Normal((*vertex_normals)(j));OpenGL_Vertex(surface.particles.X(j));
                OpenGL_Normal((*vertex_normals)(k));OpenGL_Vertex(surface.particles.X(k));}
        else
            for(int t=1;t<=elements.m;t++){
                int i,j,k;elements(t).Get(i,j,k);
                OpenGL_Normal(TRIANGLE_3D<T>::Normal(surface.particles.X(i),surface.particles.X(j),surface.particles.X(k)));
                OpenGL_Triangle(surface.particles.X(i),surface.particles.X(j),surface.particles.X(k));}
        OpenGL_End();
        //glPopAttrib();
    }
    else{
        glPushAttrib(GL_ENABLE_BIT);
        //glEnable(GL_BLEND);glDisable(GL_DEPTH_TEST);glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_COLOR_MATERIAL);
        OpenGL_Begin(GL_TRIANGLES);
        if(smooth_normals)
            for(int t=1;t<=elements.m;t++){
                int i,j,k;elements(t).Get(i,j,k);
                (*vertex_colors)(i).Send_To_GL_Pipeline();OpenGL_Normal((*vertex_normals)(i));OpenGL_Vertex(surface.particles.X(i));
                (*vertex_colors)(j).Send_To_GL_Pipeline();OpenGL_Normal((*vertex_normals)(j));OpenGL_Vertex(surface.particles.X(j));
                (*vertex_colors)(k).Send_To_GL_Pipeline();OpenGL_Normal((*vertex_normals)(k));OpenGL_Vertex(surface.particles.X(k));}
        else
            for(int t=1;t<=elements.m;t++){
                int i,j,k;elements(t).Get(i,j,k);
                OpenGL_Normal(TRIANGLE_3D<T>::Normal(surface.particles.X(i),surface.particles.X(j),surface.particles.X(k)));
                (*vertex_colors)(i).Send_To_GL_Pipeline();OpenGL_Vertex(surface.particles.X(i));
                (*vertex_colors)(j).Send_To_GL_Pipeline();OpenGL_Vertex(surface.particles.X(j));
                (*vertex_colors)(k).Send_To_GL_Pipeline();OpenGL_Vertex(surface.particles.X(k));}
        OpenGL_End();
        glPopAttrib();}
}
//#####################################################################
// Function Draw_Subsets
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Draw_Subsets() const
{
    OpenGL_Begin(GL_TRIANGLES);
    const ARRAY<VECTOR<int,3> > &elements=surface.mesh.elements;
    if(smooth_normals)
        for(int t=1;t<=subset.m;t++){
            int tri=subset(t),i,j,k;elements(tri).Get(i,j,k);
            OpenGL_Normal((*vertex_normals)(i));OpenGL_Vertex(surface.particles.X(i));
            OpenGL_Normal((*vertex_normals)(j));OpenGL_Vertex(surface.particles.X(j));
            OpenGL_Normal((*vertex_normals)(k));OpenGL_Vertex(surface.particles.X(k));}
    else
        for(int t=1;t<=subset.m;t++){
            int tri=subset(t),i,j,k;elements(tri).Get(i,j,k);
            OpenGL_Normal(TRIANGLE_3D<T>::Normal(surface.particles.X(i),surface.particles.X(j),surface.particles.X(k)));
            OpenGL_Triangle(surface.particles.X(i),surface.particles.X(j),surface.particles.X(k));}
    OpenGL_End();
}
//#####################################################################
// Function Draw_Vertices_For_Selection
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Draw_Vertices_For_Selection() const
{
    OPENGL_SELECTION::Draw_Vertices_For_Selection(surface.mesh,surface.particles);
}
//#####################################################################
// Function Draw_Segments_For_Selection
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Draw_Segments_For_Selection() const
{
    PHYSBAM_ASSERT(surface.mesh.segment_mesh);
    glPushAttrib(GL_LINE_BIT);
    glLineWidth(OPENGL_PREFERENCES::selection_line_width);
    glPushName(0);
    for(int s=1;s<=surface.mesh.segment_mesh->elements.m;s++){
        int i,j;surface.mesh.segment_mesh->elements(s).Get(i,j);
        glLoadName(s);
        OpenGL_Begin(GL_LINES);
        OpenGL_Line(surface.particles.X(i),surface.particles.X(j));
        OpenGL_End();}
    glPopName();
    glPopAttrib();
}
//#####################################################################
// Function Draw_Triangles_For_Selection
//#####################################################################
template<class T> void OPENGL_TRIANGULATED_SURFACE<T>::
Draw_Triangles_For_Selection() const
{
    glPushName(0);
    for(int t=1;t<=surface.mesh.elements.m;t++){
        int i,j,k;surface.mesh.elements(t).Get(i,j,k);
        glLoadName(t);
        OpenGL_Begin(GL_TRIANGLES);
        OpenGL_Triangle(surface.particles.X(i),surface.particles.X(j),surface.particles.X(k));
        OpenGL_End();}
    glPopName();
}
//#####################################################################
//
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_TRIANGULATED_SURFACE_VERTEX<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const TRIANGULATED_SURFACE<T> &surface=((OPENGL_TRIANGULATED_SURFACE<T> *)object)->surface;
    RANGE<VECTOR<float,3> > box(VECTOR<float,3>(surface.particles.X(index)));
    return object->World_Space_Box(box);
}

template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_TRIANGULATED_SURFACE_SEGMENT<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const TRIANGULATED_SURFACE<T> &surface=((OPENGL_TRIANGULATED_SURFACE<T> *)object)->surface;
    PHYSBAM_ASSERT(surface.mesh.segment_mesh);
    int node1,node2;surface.mesh.segment_mesh->elements(index).Get(node1,node2);
    RANGE<VECTOR<float,3> > box(VECTOR<float,3>(surface.particles.X(node1)));
    box.Enlarge_To_Include_Point(VECTOR<float,3>(surface.particles.X(node2)));
    return object->World_Space_Box(box);
}

template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_TRIANGULATED_SURFACE_TRIANGLE<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const TRIANGULATED_SURFACE<T> &surface=((OPENGL_TRIANGULATED_SURFACE<T> *)object)->surface;
    int node1,node2,node3;surface.mesh.elements(index).Get(node1,node2,node3);
    RANGE<VECTOR<float,3> > box(VECTOR<float,3>(surface.particles.X(node1)));
    box.Enlarge_To_Include_Point(VECTOR<float,3>(surface.particles.X(node2)));
    box.Enlarge_To_Include_Point(VECTOR<float,3>(surface.particles.X(node3)));
    return object->World_Space_Box(box);
}
//#####################################################################
template class OPENGL_TRIANGULATED_SURFACE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_TRIANGULATED_SURFACE<double>;
#endif

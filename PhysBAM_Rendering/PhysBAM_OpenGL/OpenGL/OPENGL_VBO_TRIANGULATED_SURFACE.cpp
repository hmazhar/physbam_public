//#####################################################################
// Copyright 2004-05, Jiayi Chong, Sergey Koltakov, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifdef WITH_OPENGL_EXTENSIONS

#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VBO_TRIANGULATED_SURFACE.h>
#include <glux.h>
GLUX_REQUIRE(GL_ARB_vertex_buffer_object);

using namespace PhysBAM;
using namespace std;


GLenum GL_Type(float) { return GL_FLOAT; }
GLenum GL_Type(double) { return GL_DOUBLE; }

template <class T>
void Set_GL_Vertex_Pointer(ARRAY<VECTOR<T,2> >& X)
{   if(sizeof(VECTOR<T,2>) != sizeof(T)*2) abort();
    glVertexPointer(2,GL_Type(T()),0,X.array);
}

template <class T>
void Set_GL_Vertex_Pointer(ARRAY<VECTOR<T,3> >& X)
{   if(sizeof(VECTOR<T,3>) != sizeof(T)*3) abort();
    glVertexPointer(3,GL_Type(T()),0,X.array);
}


//#####################################################################
// Constructor from geometry, one-sided
//#####################################################################
template<class T>
OPENGL_VBO_TRIANGULATED_SURFACE<T>::
OPENGL_VBO_TRIANGULATED_SURFACE(TRIANGULATED_SURFACE<T> &surface_input,bool smooth_normals_input,
                            const OPENGL_MATERIAL &material_input)
    : OPENGL_TRIANGULATED_SURFACE<T>(surface_input, smooth_normals_input, material_input)
{
    if(smooth_normals) Initialize_Vertex_Normals();
    useVBO = false;
}
//#####################################################################
// Constructor from geometry, two-sided
//#####################################################################
template<class T>
OPENGL_VBO_TRIANGULATED_SURFACE<T>::
OPENGL_VBO_TRIANGULATED_SURFACE(TRIANGULATED_SURFACE<T> &surface_input,bool smooth_normals_input,
                            const OPENGL_MATERIAL &front_material_input,const OPENGL_MATERIAL &back_material_input)
    : OPENGL_TRIANGULATED_SURFACE<T>(surface_input, smooth_normals_input, front_material_input, back_material_input)
{
    if(smooth_normals) Initialize_Vertex_Normals();
    useVBO = false;
}

template<class T>
OPENGL_VBO_TRIANGULATED_SURFACE<T>::~OPENGL_VBO_TRIANGULATED_SURFACE(void)
{

}

template<class T>
struct vecVBO {
    T x;
    T y;
    T z;
};

/* Attempt to use Vertex Buffer Objects for rendering*/
template<class T>
void OPENGL_VBO_TRIANGULATED_SURFACE<T>::Create_VBO(void)
{
    int p = 0, t = 0;
    Turn_Smooth_Shading_On();
    //first check to see if VBOs are available, if not default to Display Lists 
    if (gluxIsAvailable("GL_ARB_vertex_buffer_object") == GLUX_AVAILABLE) {
        ARRAY<int> &triangles = surface.triangle_mesh.triangles;
        // Generate And Bind The Vertex Buffer 
        useVBOnum = triangles.m * 3 / MAX_VBO_POINT_NUM;
        if( ((triangles.m * 3) % MAX_VBO_POINT_NUM) != 0) useVBOnum++;
        
        if(!useVBO){
            LOG::cout<<"VBOs to create: "<<useVBOnum<<std::endl;
            for(p = 0; p < useVBOnum; p++) {
                glGenBuffersARB(1, &nameVertexVBO[p]);      
            }
        }

        // Load Data into GPU
        vecVBO<T> *vertArray = new vecVBO<T>[triangles.m * 3];
        vecVBO<T> *curVert = vertArray;  

        for (t = 1; t <= triangles.m; t++) {
            int i=triangles(1,t),j=triangles(2,t),k=triangles(3,t);
            (curVert)->x = surface.particles.X(i).x;
            (curVert)->y = surface.particles.X(i).y;
            (curVert)->z = surface.particles.X(i).z;

            (curVert + 1)->x = surface.particles.X(j).x;
            (curVert + 1)->y = surface.particles.X(j).y;
            (curVert + 1)->z = surface.particles.X(j).z;

            (curVert + 2)->x = surface.particles.X(k).x;
            (curVert + 2)->y = surface.particles.X(k).y;
            (curVert + 2)->z = surface.particles.X(k).z;

            curVert = curVert + 3;
        }

        LOG::cout<<"Total number of primitives: "<<triangles.m<<std::endl;

        curVert = vertArray;
        int pointsLeft = triangles.m * 3;
        for(p = 0; p < useVBOnum; p++) {
            glBindBufferARB(GL_ARRAY_BUFFER_ARB, nameVertexVBO[p]);     
            if(pointsLeft > MAX_VBO_POINT_NUM) {
                if(!useVBO) 
                    glBufferDataARB(GL_ARRAY_BUFFER_ARB, MAX_VBO_POINT_NUM * sizeof(vecVBO<T>), curVert, GL_STATIC_DRAW_ARB );
                else
                    glBufferSubDataARB(GL_ARRAY_BUFFER_ARB, 0, MAX_VBO_POINT_NUM * sizeof(vecVBO<T>), curVert);
            }
            else {
                if(!useVBO) 
                    glBufferDataARB(GL_ARRAY_BUFFER_ARB, pointsLeft * sizeof(vecVBO<T>), curVert, GL_STATIC_DRAW_ARB );
                else 
                    glBufferSubDataARB(GL_ARRAY_BUFFER_ARB, 0, pointsLeft * sizeof(vecVBO<T>), curVert);
                
                lastVertexfilled = pointsLeft;
            }
            curVert = curVert + MAX_VBO_POINT_NUM;
            pointsLeft = pointsLeft - MAX_VBO_POINT_NUM;
        }


        delete [] vertArray;

        // Generate and Bind Normal Buffer 
        if(!useVBO) {
            for(p = 0; p < useVBOnum; p++) {
                glGenBuffersARB(1, &nameNormalVBO[p]);      
            }
        }

        // Load Data into GPU
        vecVBO<T> *normArray = new vecVBO<T>[triangles.m * 3];
        vecVBO<T> *curNorm = normArray;  

        for (t = 1; t <= triangles.m; t++) {
            int i=triangles(1,t),j=triangles(2,t),k=triangles(3,t);
            (curNorm)->x = (*vertex_normals)(i).x;
            (curNorm)->y = (*vertex_normals)(i).y;
            (curNorm)->z = (*vertex_normals)(i).z;

            (curNorm + 1)->x = (*vertex_normals)(j).x;
            (curNorm + 1)->y = (*vertex_normals)(j).y;
            (curNorm + 1)->z = (*vertex_normals)(j).z;

            (curNorm + 2)->x = (*vertex_normals)(k).x;
            (curNorm + 2)->y = (*vertex_normals)(k).y;
            (curNorm + 2)->z = (*vertex_normals)(k).z;

            curNorm = curNorm + 3;
        } 

        curNorm = normArray;
        pointsLeft = triangles.m * 3;
        for(p = 0; p < useVBOnum; p++) {
            glBindBufferARB(GL_ARRAY_BUFFER_ARB, nameNormalVBO[p]);     
            if(pointsLeft > MAX_VBO_POINT_NUM) {
                if(!useVBO)
                    glBufferDataARB(GL_ARRAY_BUFFER_ARB, MAX_VBO_POINT_NUM * sizeof(vecVBO<T>), curNorm, GL_STATIC_DRAW_ARB );
                else
                    glBufferSubDataARB(GL_ARRAY_BUFFER_ARB, 0, MAX_VBO_POINT_NUM * sizeof(vecVBO<T>), curNorm);
            }
            else {
                if(!useVBO)
                    glBufferDataARB(GL_ARRAY_BUFFER_ARB, pointsLeft * sizeof(vecVBO<T>), curNorm, GL_STATIC_DRAW_ARB );
                else
                    glBufferSubDataARB(GL_ARRAY_BUFFER_ARB, 0, pointsLeft * sizeof(vecVBO<T>), curNorm);
                
                lastNormalfilled = pointsLeft;
            }
            curNorm = curNorm + MAX_VBO_POINT_NUM;
            pointsLeft = pointsLeft - MAX_VBO_POINT_NUM;
        }

        delete [] normArray; 

        useVBO = true;
    }
    else {
        LOG::cout<<"Vertex Buffer Objects NOT SUPPORTED! Defaulting to Display Lists."<<std::endl;
        Create_Display_List();
        useVBO = false;
    }
}

template<class T>
void OPENGL_VBO_TRIANGULATED_SURFACE<T>::Use_VBO(void)
{
    useVBO = true;
}

template<class T>
void OPENGL_VBO_TRIANGULATED_SURFACE<T>::
Display(const int in_color) const
{
    if(useVBO) {
        if (two_sided) {
            glDisable (GL_CULL_FACE);
            glLightModeli (GL_LIGHT_MODEL_TWO_SIDE, 1);
            if (in_color) {
                front_material.Send_To_GL_Pipeline(GL_FRONT);
                back_material.Send_To_GL_Pipeline(GL_BACK);
            } else {
                front_material_gray.Send_To_GL_Pipeline(GL_FRONT);
                back_material_gray.Send_To_GL_Pipeline(GL_BACK);
            }
        } else {
            if (in_color)   front_material.Send_To_GL_Pipeline();
            else            front_material_gray.Send_To_GL_Pipeline();
        }

        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        Send_Transform_To_GL_Pipeline();
    
        GLint mode;
        glGetIntegerv(GL_RENDER_MODE, &mode);

        if(wireframe_only){
            glPushAttrib(GL_POLYGON_BIT);
            glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        }
        
        //do VBO rendering
        glEnableClientState(GL_VERTEX_ARRAY);   
        glEnableClientState(GL_NORMAL_ARRAY);

        for(int p = 0; p < useVBOnum; p++) {
            glBindBufferARB(GL_ARRAY_BUFFER_ARB, nameVertexVBO[p]);
            glVertexPointer(3,GL_Type(T()),0,(char *) 0 );
            glBindBufferARB(GL_ARRAY_BUFFER_ARB, nameNormalVBO[p]);
            glNormalPointer(GL_Type(T()),0,(char *) 0 );
            if(p == (useVBOnum - 1)) {
                glDrawArrays(GL_TRIANGLES,0,lastVertexfilled);
            }
            else {
                glDrawArrays(GL_TRIANGLES,0,MAX_VBO_POINT_NUM);
            }
        }


        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);

        if(wireframe_only) glPopAttrib();

        if(two_sided){
            glLightModeli (GL_LIGHT_MODEL_TWO_SIDE, 0);
            glEnable (GL_CULL_FACE);}

        glPopMatrix();

    }
    else {
        Org_Display(in_color);
    }
}


template<class T>
void OPENGL_VBO_TRIANGULATED_SURFACE<T>::
Org_Display(const int in_color) const
{   
    if (two_sided) {
        glDisable (GL_CULL_FACE);
        glLightModeli (GL_LIGHT_MODEL_TWO_SIDE, 1);
        if (in_color) {
            front_material.Send_To_GL_Pipeline(GL_FRONT);
            back_material.Send_To_GL_Pipeline(GL_BACK);
        } else {
            front_material_gray.Send_To_GL_Pipeline(GL_FRONT);
            back_material_gray.Send_To_GL_Pipeline(GL_BACK);
        }
    } else {
        if (in_color)   front_material.Send_To_GL_Pipeline();
        else            front_material_gray.Send_To_GL_Pipeline();
    }

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    GLint mode;
    glGetIntegerv(GL_RENDER_MODE, &mode);

    if(wireframe_only){
        glPushAttrib(GL_POLYGON_BIT);
        glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    }

    if(mode == GL_SELECT){
        glPushName(1);
        Draw_Vertices_For_Selection();
        if (surface.triangle_mesh.segment_mesh){
            glLoadName(2);
            Draw_Segments_For_Selection();
        }
        glLoadName(3);
        Draw_Triangles_For_Selection();
        glPopName();
    }
    else if(use_display_list)
        glCallList(display_list_id);
    else
        Draw();
    if(wireframe_only) glPopAttrib();

    if(mode != GL_SELECT) {
        if(current_selection){
            if(current_selection->type == OPENGL_SELECTION::TRIANGULATED_SURFACE_VERTEX){
                int index=((OPENGL_SELECTION_TRIANGULATED_SURFACE_VERTEX<T> *)current_selection)->index;
                OPENGL_SELECTION::Draw_Highlighted_Vertex(surface.particles.X(index));
            }
            else if(current_selection->type == OPENGL_SELECTION::TRIANGULATED_SURFACE_SEGMENT && surface.triangle_mesh.segment_mesh){
                int index=((OPENGL_SELECTION_TRIANGULATED_SURFACE_SEGMENT<T> *)current_selection)->index;
                int node1,node2;surface.triangle_mesh.segment_mesh->segments.Get(index,node1,node2);
                OPENGL_SELECTION::Draw_Highlighted_Segment(surface.particles.X(node1),surface.particles.X(node2));
            }
            else if(current_selection->type == OPENGL_SELECTION::TRIANGULATED_SURFACE_TRIANGLE){
                int index=((OPENGL_SELECTION_TRIANGULATED_SURFACE_TRIANGLE<T> *)current_selection)->index;
                int node1,node2,node3;surface.triangle_mesh.triangles.Get(index,node1,node2,node3);
                OPENGL_SELECTION::Draw_Highlighted_Triangle_Boundary(surface.particles.X(node1),surface.particles.X(node2),surface.particles.X(node3));
            }
        }

        if(highlight_current_node) {
            Highlight_Current_Node();
            Draw_Triangles_Incident_On_Current_Node();
        } else if(highlight_boundary) {
            glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
            glDisable(GL_LIGHTING);
            glLineWidth(5);
            if(!surface.triangle_mesh.boundary_mesh) surface.triangle_mesh.Initialize_Boundary_Mesh();
            if(!surface.triangle_mesh.boundary_mesh->connected_segments) surface.triangle_mesh.boundary_mesh->Initialize_Connected_Segments();
            glPopAttrib();
        }
    }

    if(two_sided){
        glLightModeli (GL_LIGHT_MODEL_TWO_SIDE, 0);
        glEnable (GL_CULL_FACE);}

    glPopMatrix();
}


template class OPENGL_VBO_TRIANGULATED_SURFACE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_VBO_TRIANGULATED_SURFACE<double>;
#endif

#endif

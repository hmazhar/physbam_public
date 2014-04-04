//#####################################################################
// Copyright 2004-2006, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SCENE_WINDOW
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/IRRADIANCE_SAMPLE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Lights/RENDERING_LIGHT.h>
#include "SCENE_WINDOW.h"
#include <Fl/Fl.h>
using namespace PhysBAM;
//#####################################################################
// draw
//#####################################################################
template<class T> void SCENE_WINDOW<T>::
draw()
{
    if(!valid())Init_Gl();
    // setup camera
    glMatrixMode(GL_PROJECTION);glLoadIdentity();
    gluPerspective(60,float(T(w())/T(h())),.3,1000);
    glScalef (-1,1,1); // convert to a left-handed coordinate system,consistent with PhysBAM
    glMatrixMode(GL_MODELVIEW);glLoadIdentity();
    GLfloat light_position[]={0,0,0,1};glLightfv(GL_LIGHT0,GL_POSITION,light_position);
    VECTOR<T,3> position=world.camera.position,look_vector=world.camera.look_vector;
    glTranslatef(0,0,(GLfloat)-camera_distance);glMultMatrixf(arcball_matrix.x);glMultMatrixf(rotation_matrix.x);OpenGL_Translate(-camera_target);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    
    // Draw rays
    Draw_Ray(debug_ray);
    if(selected_ray){
        glBegin(GL_LINES);
        VECTOR<T,3> hit_point=selected_ray->ray.ray.Point(selected_ray->ray.ray.t_max);
        glColor3f(1,0,0);OpenGL_Vertex(hit_point);OpenGL_Vertex(selected_ray->same_side_normal+hit_point);
        glEnd();}
    
    // Draw objects
    glEnable(GL_LIGHTING);
    glColor3f(0.7f,0.7f,0.7f);
    if(world.standard_objects.m==object_triangles.m)for(int i=1;i<=world.standard_objects.m;i++)if(object_triangles(i)){
        MATRIX<float,4> matrix;for(int k=0;k<16;k++)matrix.x[k]=world.standard_objects(i)->transform.x[k];
        glPushMatrix();glMultMatrixf(matrix.x);
        if(world.standard_objects(i)==selected_object){glColor4f(1.0f,1.0f,0.7f,0.3f);glEnable(GL_BLEND);glDepthMask(GL_FALSE);}else glColor4f(0.7f,0.7f,0.7f,0.3f);
        glBegin(GL_TRIANGLES);  
        for(int j=1;j<=object_triangles(i)->mesh.elements.m;j++){
            VECTOR<int,3> nodes=object_triangles(i)->mesh.elements(j);
            TRIANGLE_3D<T> triangle(object_triangles(i)->particles.X(nodes(1)),object_triangles(i)->particles.X(nodes(2)),object_triangles(i)->particles.X(nodes(3)));
            OpenGL_Normal(triangle.Normal());OpenGL_Vertex(triangle.x1);OpenGL_Vertex(triangle.x2);OpenGL_Vertex(triangle.x3);}
        glEnd();glPopMatrix();
        glDisable(GL_BLEND);glDepthMask(GL_TRUE);}
    glDisable(GL_LIGHTING);
    
    // Draw others
    Draw_Photons();
    Draw_Irradiance_Cache();
    // Draw Subsurface Scattering samples
    for(int i=1;i<=world.standard_objects.m;++i)if(world.standard_objects(i)->bssrdf_shader) Draw_Subsurface_Scattering_Samples(world.standard_objects(i)->bssrdf_tree);
    Draw_Subsurface_Candidate_List();
    Draw_Camera();
}
//#####################################################################
// Initialize_Objects
//#####################################################################
template<class T> void SCENE_WINDOW<T>::Initialize_Objects()
{
    object_triangles.Remove_All();
    for(int i=1;i<=world.standard_objects.m;i++){RENDERING_OBJECT<T>* object=world.standard_objects(i);object_triangles.Append(object->Generate_Triangles());}
    selected_object=0;
}
//#####################################################################
// Draw_Irradiance_Cache
//#####################################################################
template<class T> void Draw_Octree(OCTREE_GRID<T>& octree,OCTREE_CELL<T>* cell)
{
    glPushAttrib(GL_LINE_BIT);
    glLineWidth(1);
    glDisable(GL_LINE_SMOOTH);
    glBegin(GL_LINES);
    ARRAY<VECTOR<T,3> ,VECTOR<int,1> > box_points(0,7);for(int i=0;i<box_points.m;i++)box_points(i)=octree.Node_Location(i,cell);
    glBegin(GL_LINES);
    OpenGL_Vertex(box_points(0));OpenGL_Vertex(box_points(1));
    OpenGL_Vertex(box_points(0));OpenGL_Vertex(box_points(2));
    OpenGL_Vertex(box_points(0));OpenGL_Vertex(box_points(4));
    OpenGL_Vertex(box_points(1));OpenGL_Vertex(box_points(3));
    OpenGL_Vertex(box_points(1));OpenGL_Vertex(box_points(5));
    OpenGL_Vertex(box_points(2));OpenGL_Vertex(box_points(3));
    OpenGL_Vertex(box_points(2));OpenGL_Vertex(box_points(6));
    OpenGL_Vertex(box_points(3));OpenGL_Vertex(box_points(7));
    OpenGL_Vertex(box_points(4));OpenGL_Vertex(box_points(6));
    OpenGL_Vertex(box_points(4));OpenGL_Vertex(box_points(5));
    OpenGL_Vertex(box_points(5));OpenGL_Vertex(box_points(7));
    OpenGL_Vertex(box_points(6));OpenGL_Vertex(box_points(7));glEnd();
    if(cell->Has_Children())for(int i=0;i<8;i++)Draw_Octree(octree,cell->Child(i));
    glPopAttrib();
}
template<class T> void SCENE_WINDOW<T>::Draw_Irradiance_Cache()
{
    if(display_irradiance_cache)Draw_Octree(world.irradiance_cache.octree_grid,world.irradiance_cache.octree_grid.cells(1,1,1));
}
//#####################################################################
// Draw_Subsurface_Scattering_Samples
//#####################################################################
template<class T> void SCENE_WINDOW<T>::Draw_Subsurface_Scattering_Samples(SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<T>* bssrdf_tree)
{
    if(!display_subsurface_scattering_samples) return;
    glPointSize(3.0f);
    
    glBegin(GL_POINTS);
    for(int i=1;i<=bssrdf_tree->samples.m;++i){
        VECTOR<T,3> curval=bssrdf_tree->samples(i).transmitted_irradiance;
        glColor3f(curval.x, curval.y, curval.z);
        OpenGL_Vertex(bssrdf_tree->samples(i).position);
    }
    glEnd();
    
    glPointSize(1.0f);
    glBegin(GL_LINES);
    for(int i=1;i<=bssrdf_tree->samples.m;++i){
        glColor3f(0.0,0.0,1.0);
        OpenGL_Vertex(bssrdf_tree->samples(i).position);
        OpenGL_Vertex(bssrdf_tree->samples(i).position + (T)0.25 * bssrdf_tree->samples(i).normal / 5);
        glColor3f(0.0,1.0,1.0);
        OpenGL_Vertex(bssrdf_tree->samples(i).position + (T)0.25 * bssrdf_tree->samples(i).normal / 5);
        OpenGL_Vertex(bssrdf_tree->samples(i).position + bssrdf_tree->samples(i).normal);
    }
    glEnd(); 
    
//    Draw_Octree(bssrdf_tree->octree_grid,bssrdf_tree->octree_grid.cells(1,1,1));
}

//#####################################################################
// Draw_Subsurface_Candidate_List
//#####################################################################
template<class T> void SCENE_WINDOW<T>::Draw_Subsurface_Candidate_List()
{
    if(selected_ray){
        for(int i=1;i<=world.standard_objects.m;++i)
            if(world.standard_objects(i)->bssrdf_shader){
                T totalArea=0,position_magnitude=0;VECTOR<T,3>totalIrradiance(0,0,0);
                ARRAY<SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE<T> > sample_list;
                VECTOR<T,3> hit_point=selected_ray->ray.ray.Point(selected_ray->ray.ray.t_max);
                world.standard_objects(i)->bssrdf_tree->Find_Irradiance_Candidates(sample_list,hit_point);
                
                glPointSize(3.0f);
                glDisable(GL_DEPTH_TEST);
                glBegin(GL_POINTS);
                for(int i=1;i<=sample_list.m;++i){
//                VECTOR<T,3> curval=sample_list(i).transmitted_irradiance;
//                glColor3f(curval.x, curval.y, curval.z);
                    glColor3f(0.0, 1.0, 0.0);
                    OpenGL_Vertex(sample_list(i).position);
                    totalArea+=sample_list(i).area;
                    totalIrradiance+=sample_list(i).transmitted_irradiance*sample_list(i).area;
                    position_magnitude+=(sample_list(i).position-hit_point).Magnitude();
                }
                LOG::cerr<<"Area: " << totalArea << std::endl;
                LOG::cerr<<"Irradiance: " << totalIrradiance << std::endl;
                LOG::cerr<<"Samples: " << sample_list.m << std::endl;
                LOG::cerr<<"Sum of Position Magnitude: " << position_magnitude << std::endl;
                glPointSize(3.0f);
                glColor3f(0.0, 1.0, 1.0);
                OpenGL_Vertex(hit_point);
                glEnd();glEnable(GL_DEPTH_TEST);}
    } 
}

//#####################################################################
// Draw_Photons
//#####################################################################
template<class T> void SCENE_WINDOW<T>::Draw_Photons()
{
    glPointSize(1.0f);
    
    if(display_global_photon_map)Draw_Photons(world.global_photon_map);
    if(display_caustic_photon_map)Draw_Photons(world.caustic_photon_map);
    if(display_volume_photon_map)Draw_Photons(world.volume_photon_map);
    if(display_used_photons){// display photons that were used in yellow
        glDisable(GL_DEPTH_TEST);
        glBegin(GL_POINTS);
        if(selected_ray)for(int i=1;i<=selected_ray->photons_used.m;i++){
            PHOTON<T>& photon=*selected_ray->photons_used(i);
            glColor3f(1,1,0);
            OpenGL_Vertex(photon.location);}
        glEnd();
        glEnable(GL_DEPTH_TEST);}
}
//#####################################################################
// Draw_Photon_Map
//#####################################################################
template<class T> void SCENE_WINDOW<T>::Draw_Photons(const PHOTON_MAP<T>& photon_map)
{
    glPointSize(1.0f);
    glBegin(GL_POINTS);
    for(int i=1;i<=photon_map.Max_Number_Of_Photons();i++){
        VECTOR<T,3> color=(photon_map.Photon_Power(i)).Normalized();
        glColor3f(color.x,color.y,color.z);
        //printf("Drawing photon %d at %f %f %f\n",i,world.global_photon_map.Photon_Position(i).x,world.global_photon_map.Photon_Position(i).y,world.global_photon_map.Photon_Position(i).z);
        OpenGL_Vertex(photon_map.Photon_Position(i));}
    glEnd();
}
//#####################################################################
// Draw_Camera
//#####################################################################
template<class T> void SCENE_WINDOW<T>::Draw_Camera()
{
    glEnable(GL_BLEND);glColor4f(1,1,1,0.5f);
    glBegin(GL_QUADS);
    VECTOR<T,3> focal_point=world.camera.focal_point;
    T image_width=world.camera.film.width,image_height=world.camera.film.height;
    VECTOR<T,3> u=world.camera.horizontal_vector*image_width/T(2);
    VECTOR<T,3> v=world.camera.vertical_vector*image_height/T(2);
    OpenGL_Vertex(focal_point-u-v);OpenGL_Vertex(focal_point-u+v);OpenGL_Vertex(focal_point+u+v);OpenGL_Vertex(focal_point+u-v);
    glEnd();
    glLineWidth(0.1f);glBegin(GL_LINES);
    OpenGL_Vertex(world.camera.position);OpenGL_Vertex(focal_point-u-v);OpenGL_Vertex(world.camera.position);OpenGL_Vertex(focal_point-u+v);
    OpenGL_Vertex(world.camera.position);OpenGL_Vertex(focal_point+u+v);OpenGL_Vertex(world.camera.position);OpenGL_Vertex(focal_point+u-v);
    glEnd();
    glDisable(GL_BLEND);
}
//#####################################################################
// Draw_Ray
//#####################################################################
template<class T> void SCENE_WINDOW<T>::Draw_Ray(RENDERING_RAY_DEBUG<T>* debug_ray)
{
    if(!debug_ray)return;
    if(debug_ray->parent){// don't draw the top level one because it is the suck
        glBegin(GL_LINES);
        VECTOR<T,3> point1=debug_ray->ray.ray.endpoint,point2;
        if(debug_ray->ray.ray.semi_infinite)point2=debug_ray->ray.ray.Point(10000);
        else point2=debug_ray->ray.ray.Point(debug_ray->ray.ray.t_max);
        
        VECTOR<T,3> start_color,end_color;
        if(debug_ray==selected_ray){start_color=VECTOR<T,3>(1,1,0);end_color=VECTOR<T,3>(1,1,0.5);}
        else if(debug_ray->ray.ray_type==RENDERING_RAY<T>::COLOR_RAY){start_color=VECTOR<T,3>(0.3,0.3,1);end_color=VECTOR<T,3>(0.8,0.8,1.0);}
        else {start_color=VECTOR<T,3>(0.3,1,0.3);end_color=VECTOR<T,3>(0.8,1,0.8);}
        glColor3f((GLfloat)start_color.x,(GLfloat)start_color.y,(GLfloat)start_color.z);
        glVertex3f((GLfloat)point1.x,(GLfloat)point1.y,(GLfloat)point1.z);
        glColor3f((GLfloat)end_color.x,(GLfloat)end_color.y,(GLfloat)end_color.z);
        glVertex3f((GLfloat)point2.x,(GLfloat)point2.y,(GLfloat)point2.z);
        glEnd();
        glPointSize(4.0);
        glColor3f(1,0,0);
        glBegin(GL_POINTS);
        glVertex3f((GLfloat)point2.x,(GLfloat)point2.y,(GLfloat)point2.z);
        glEnd();}
    for(int i=1;i<=debug_ray->children.m;i++)Draw_Ray(debug_ray->children(i));
}
//#####################################################################
// Select_Ray
//#####################################################################
template<class T> void SCENE_WINDOW<T>::Select_Ray(RENDERING_RAY_DEBUG<T>* debug_ray)
{
    selected_ray=debug_ray;
    damage(1);
}
//#####################################################################
// Select_Object
//#####################################################################
template<class T> void SCENE_WINDOW<T>::Select_Object(RENDERING_OBJECT<T>* incoming_selection)
{
    selected_object=incoming_selection;
    damage(1);
}
//#####################################################################
// Init_Gl
//#####################################################################
template<class T> void SCENE_WINDOW<T>::
Init_Gl()
{
    glClearColor(0,0,0,1);
    glLineWidth(2.0);
    glEnable(GL_DEPTH_TEST);glEnable(GL_LINE_SMOOTH);
    glFrontFace (GL_CW); // because we're left-handed
    // lighting
    glEnable(GL_LIGHT0);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,1);
    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_NORMALIZE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE);
    // set window size
    glViewport(0,0,w(),h());
}
//#####################################################################
// Set_Debug_Ray
//#####################################################################
template<class T> void SCENE_WINDOW<T>::
Set_Debug_Ray(RENDERING_RAY_DEBUG<T>* debug_ray)
{
    if(this->debug_ray)delete this->debug_ray;
    this->debug_ray=debug_ray;
    this->selected_ray=0;
    damage(1);
}
//#####################################################################
// Convert_Mouse_Coordinates
//#####################################################################
template<class T> VECTOR<GLfloat,2> SCENE_WINDOW<T>::
Convert_Mouse_Coordinates(int x,int y)
{
    VECTOR<GLfloat,2> coord;T width=(T)w(),height=(T)h();
    if(width>=height){coord.x =-((T)x/width-T(0.5))*2*(width/height);coord.y =-((T)y/height-T(0.5))*2;}
    else{coord.x=-((T)x/width-T(0.5))*2;coord.y=-((T)y/height-T(0.5))*2*(height/width);}
    return coord;
}

//#####################################################################
// handle
//#####################################################################
template<class T> int SCENE_WINDOW<T>::
handle(int event)
{
    if(Fl::event_button()==1){
        if(event==FL_PUSH){
            arcball.Begin_Drag(Convert_Mouse_Coordinates(Fl::event_x(),Fl::event_y()));
            damage(1);camera_rotation=true;return 1;}
        else if(event==FL_DRAG){
            arcball.Update(Convert_Mouse_Coordinates(Fl::event_x(),Fl::event_y()));
            arcball_matrix=arcball.Value();
            damage(1);return 1;}
        else if(event==FL_RELEASE){
            arcball.End_Drag(Convert_Mouse_Coordinates(Fl::event_x(),Fl::event_y()));
            arcball_matrix=arcball.Value();rotation_matrix=arcball_matrix*rotation_matrix;
            arcball_matrix=MATRIX<GLfloat,4>::Identity_Matrix();
            camera_rotation=false;damage(1);return 1;}}
    else if(Fl::event_button()==2){
        if(event==FL_PUSH){
            damage(1);camera_translate=true;old_mouse_vector=Convert_Mouse_Coordinates(Fl::event_x(),Fl::event_y());return 1;}
        else if(event==FL_DRAG){
            // Find the unprojected points
            VECTOR<GLfloat,2> new_mouse_vector=Convert_Mouse_Coordinates(Fl::event_x(),Fl::event_y());
            make_current();
            camera_target=Translate_In_Camera_Plane(old_mouse_vector,new_mouse_vector);
            old_mouse_vector=new_mouse_vector;
            damage(1);return 1;}
        else if(event==FL_RELEASE){
            camera_translate=false;damage(1);return 1;}}
    else if(Fl::event_button()==3){
        if(event==FL_PUSH){
            damage(1);camera_zoom=true;old_mouse_vector=Convert_Mouse_Coordinates(Fl::event_x(),Fl::event_y());return 1;}
        else if(event==FL_DRAG){
            VECTOR<GLfloat,2> new_mouse_vector=Convert_Mouse_Coordinates(Fl::event_x(),Fl::event_y());
            T factor=pow((T)1.5,-(T)(new_mouse_vector-old_mouse_vector).y);
            camera_distance*=factor;old_mouse_vector=new_mouse_vector;
            damage(1);return 1;}
        else if(event==FL_RELEASE){
            camera_zoom=false;damage(1);return 1;}}
    return Fl_Gl_Window::handle(event);
}
//#####################################################################
// Get_Look_At
//#####################################################################
template<class T> void SCENE_WINDOW<T>::
Set_Look_At(const VECTOR<GLfloat,3> &camera,const VECTOR<GLfloat,3> &target,const VECTOR<GLfloat,3> &up)
{
    VECTOR<GLfloat,3> view_forward = target - camera;
    VECTOR<GLfloat,3> view_up = (up - up.Projected(view_forward)).Normalized();
    VECTOR<GLfloat,3> view_right = VECTOR<GLfloat,3>::Cross_Product(view_up,view_forward);
    view_right.Normalize();
    camera_target=target;
    camera_distance=view_forward.Normalize();
    arcball_matrix=MATRIX<GLfloat,4>::Identity_Matrix();
    rotation_matrix=MATRIX<GLfloat,4>(-view_right.x,view_up.x,-view_forward.x,0,-view_right.y,view_up.y,-view_forward.y,0,-view_right.z,view_up.z,-view_forward.z,0,0,0,0,1);
}
//#####################################################################
// Translate_In_Camera_Plane
//#####################################################################
template<class T> VECTOR<GLfloat,3> SCENE_WINDOW<T>::
Translate_In_Camera_Plane(VECTOR<GLfloat,2> old_mouse_vector,VECTOR<GLfloat,2> new_mouse_vector)
{
    GLint viewport[4];GLdouble mvmatrix[16],projmatrix[16];
    GLdouble wx,wy,wz;  
    GLdouble win_x,win_y,win_z; 
    glGetIntegerv (GL_VIEWPORT,viewport);
    glGetDoublev (GL_MODELVIEW_MATRIX,mvmatrix);
    glGetDoublev (GL_PROJECTION_MATRIX,projmatrix);
    float dx = (float)(100*(new_mouse_vector-old_mouse_vector).x);
    float dy = (float)(100*(old_mouse_vector-new_mouse_vector).y);
    gluProject((GLdouble) camera_target.x,(GLdouble) camera_target.y,(GLdouble) camera_target.z,mvmatrix,projmatrix,viewport,&win_x,&win_y,&win_z);
    gluUnProject((GLdouble) win_x+dx,(GLdouble) win_y+dy,win_z,mvmatrix,projmatrix,viewport,&wx,&wy,&wz);
    return VECTOR<GLfloat,3>(wx,wy,wz);
}
template class SCENE_WINDOW<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SCENE_WINDOW<double>;
#endif

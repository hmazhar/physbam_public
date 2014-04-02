//#####################################################################
// Copyright 2008, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __H_OPENGL_ARCBALL__
#define __H_OPENGL_ARCBALL__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_CYLINDER_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_PLANE_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_SPHERE_INTERSECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WINDOW.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>

namespace PhysBAM{

#define LG_NSEGS 4
#define NSEGS (1<<LG_NSEGS)

template<class T>
class OPENGL_ARCBALL: public OPENGL_OBJECT
{
    typedef VECTOR<T,3> TV;
public:
    OPENGL_WORLD& opengl_world;

    SPHERE<TV> sphere;
    ROTATION<TV> qNow,qDown,qDrag;
    VECTOR<T,2> center,vDown;
    TV vFrom,vTo,vrFrom,vrTo;
    MATRIX<T,4> mNow,mDown,mDeltaNow; //note that these only work in world space
    bool dragging,use_sphere_center,use_object_space;
    OPENGL_WORLD *world;
    int rotation_axis;
    OPENGL_COLOR x_axis,y_axis,z_axis,outer_rim,highlight;

    OPENGL_ARCBALL(OPENGL_WORLD& opengl_world_input)
        :opengl_world(opengl_world_input),dragging(false),use_sphere_center(false),use_object_space(false),rotation_axis(-1)
    {mNow=mDown=mDeltaNow=MATRIX<T,4>::Identity_Matrix();
    x_axis=OPENGL_COLOR::Red();y_axis=OPENGL_COLOR::Green();z_axis=OPENGL_COLOR::Blue();
    outer_rim=OPENGL_COLOR::Ground_Tan();highlight=OPENGL_COLOR::Yellow();}

    void Reinitialize()
    {mNow=mDown=mDeltaNow=MATRIX<T,4>::Identity_Matrix();
    sphere=SPHERE<TV>();
    qNow=qDown=qDrag=ROTATION<TV>();
    center=vDown=VECTOR<T,2>();
    vFrom=vTo=vrFrom=vrTo=TV();}

    void Update(const VECTOR<T,2> &vNow)
    {if(use_object_space) Update_Obj(vNow);
    else Update_World(vNow);}

    void Update_World(const VECTOR<T,2> &vNow)
    {if(use_sphere_center){
        GLdouble project_x,project_y,project_z;
        GLdouble model_view[16];glGetDoublev(GL_MODELVIEW_MATRIX,model_view);
        GLdouble projection[16];glGetDoublev(GL_PROJECTION_MATRIX,projection);
        GLint viewport[4];glGetIntegerv(GL_VIEWPORT,viewport);
        gluProject((GLdouble)sphere.center.x,(GLdouble)sphere.center.y,(GLdouble)sphere.center.z,model_view,projection,viewport,&project_x,&project_y,&project_z);
        center.x=(T)project_x;center.y=(T)project_y;}
    vFrom=MouseOnSphere(vDown,center,sphere.radius);vTo=MouseOnSphere(vNow,center,sphere.radius);
    if (dragging){qDrag=Qt_FromBallPoints(vFrom,vTo);qNow=qDrag*qDown;}
    Qt_ToBallPoints(qDown,vrFrom,vrTo);
    mNow=MATRIX<T,4>::From_Linear(qNow.Rotation_Matrix());
    mDeltaNow=MATRIX<T,4>::From_Linear(qDrag.Rotation_Matrix());}

    void Update_Obj(const VECTOR<T,2> &vNow)
    {if(!use_sphere_center){
        GLdouble unproject_x,unproject_y,unproject_z;
        GLdouble model_view[16];glGetDoublev(GL_MODELVIEW_MATRIX,model_view);
        GLdouble projection[16];glGetDoublev(GL_PROJECTION_MATRIX,projection);
        GLint viewport[4];glGetIntegerv(GL_VIEWPORT,viewport);
        gluProject((GLdouble)center.x,(GLdouble)center.y,(GLdouble)0,model_view,projection,viewport,&unproject_x,&unproject_y,&unproject_z);
        sphere.center.x=(T)unproject_x;sphere.center.y=(T)unproject_y;sphere.center.z=(T)unproject_z;}
    vFrom=MouseOnSphere(vDown,sphere,vFrom);vTo=MouseOnSphere(vNow,sphere,vTo);
    if (dragging){qDrag=Qt_FromBallPoints(vFrom,vTo);qNow=qDrag*qDown;}
    Qt_ToBallPoints(qDown,vrFrom,vrTo);}

    MATRIX<T,4> Value()
    {return mDeltaNow;}

    void Begin_Drag(const VECTOR<T,2> &vNow)
    {dragging=true;vDown=vNow;Update(vNow);}

    void End_Drag(const VECTOR<T,2> &vNow)
    {dragging=false;Update(vNow);qDown=qNow;mDown=mNow;rotation_axis=-1;}
    
    void Display(const int color) const
    {glDisable(GL_LIGHTING);DrawOuterRing();DrawResultArc();DrawDragArc();glEnable(GL_LIGHTING);}
    
private:
    void DrawColor(const OPENGL_COLOR &color,int roation_axis,int my_axis) const
    {if(rotation_axis==my_axis) highlight.Send_To_GL_Pipeline();
    else color.Send_To_GL_Pipeline();}

    void DrawAnyArc(const TV &vFrom,const TV &vTo) const
    {assert(vFrom!=vTo);DrawAnyArcObj(vFrom,vTo);return;}

    void DrawAnyArcWorld(const TV &vFrom,const TV &vTo) const
    {assert(vFrom!=vTo);int i;TV sphere_center;
    if(use_sphere_center){
        GLdouble model_view[16];glGetDoublev(GL_MODELVIEW_MATRIX,model_view);
        GLdouble projection[16];glGetDoublev(GL_PROJECTION_MATRIX,projection);
        GLint viewport[4];viewport[0]=viewport[1]=0;viewport[2]=viewport[3]=2;
        //glGetIntegerv(GL_VIEWPORT,viewport);
        GLdouble project_x,project_y,project_z;
        gluProject((GLdouble)sphere.center.x,(GLdouble)sphere.center.y,(GLdouble)sphere.center.z,model_view,projection,viewport,&project_x,&project_y,&project_z);
        sphere_center.x=(T)project_x-(T)1;sphere_center.y=(T)project_y-(T)1;sphere_center.z=(T)(2.*project_z-1.);}
    TV pts[NSEGS+1];
    double dot;
    pts[0]=vFrom;
    pts[1]=pts[NSEGS]=vTo;
    for (i=0;i<LG_NSEGS;i++) pts[1]=Bisect_Vectors(pts[0],pts[1]);
    dot=2.0*TV::Dot_Product(pts[0],pts[1]);
    for (i=2;i<NSEGS;i++) pts[i]=(pts[i-1]*dot)-pts[i-2];
    glMatrixMode(GL_MODELVIEW);glPushMatrix();glLoadIdentity();
    glMatrixMode(GL_PROJECTION);glPushMatrix();glLoadIdentity();
    glDepthMask(0);glDisable(GL_DEPTH_TEST);
    OpenGL_Begin(GL_LINE_STRIP);
    for (i=0;i<=NSEGS;i++){pts[i].y*=(T)world->window->Width()/(T)world->window->Height();OpenGL_Vertex(sphere_center+(T)4.8*sphere.radius*pts[i]);}
    OpenGL_End();
    glPopMatrix();glMatrixMode(GL_MODELVIEW);glPopMatrix();
    glDepthMask(1);glEnable(GL_DEPTH_TEST);}

    void DrawAnyArcObj(const TV &vFrom,const TV &vTo) const
    {assert(vFrom!=vTo);int i;TV sphere_center;
    if(use_sphere_center) sphere_center=sphere.center;
    T radius=sphere.radius*(sphere_center-world->Get_Camera_Position()).Magnitude();
    TV pts[NSEGS+1];
    double dot;
    pts[0]=vFrom;
    pts[1]=pts[NSEGS]=vTo;
    for (i=0;i<LG_NSEGS;i++) pts[1]=Bisect_Vectors(pts[0],pts[1]);
    dot=2.0*TV::Dot_Product(pts[0],pts[1]);
    for (i=2;i<NSEGS;i++) pts[i]=(pts[i-1]*dot)-pts[i-2];
    glDepthMask(0);glDisable(GL_DEPTH_TEST);
    OpenGL_Begin(GL_LINE_STRIP);
    TV camera=TV(world->Get_Camera_Position()-world->Get_Target_Position()).Normalized();
    for (i=0;i<=NSEGS;i++){
        if(TV::Dot_Product(qDown.Rotate(pts[i]),camera)>-1e-8||TV::Dot_Product(qNow.Rotate(pts[i]),camera)>1e-8) OpenGL_Vertex(sphere_center+qNow.Rotate(radius*pts[i]));
        else{OpenGL_End();OpenGL_Begin(GL_LINE_STRIP);}}
    OpenGL_End();
    glDepthMask(1);glEnable(GL_DEPTH_TEST);}

    void DrawHalfArc(const TV &n) const
    {TV p,m;p.z=0.0;
    if(n.z!=1.0){p.x=n.y;p.y=-n.x;p.Normalize();}
    else{p.x=0.0;p.y=1.0;}
    m=TV::Cross_Product(p,n);
    DrawAnyArc(p,m);DrawAnyArc(m,-p);}

    void DrawOuterRing() const
    {Circ();}

    void DrawDragArc() const
    {if(dragging&&vFrom!=vTo){
        TV sphere_center;if(use_sphere_center) sphere_center=sphere.center;
        OPENGL_COLOR::Gray().Send_To_GL_Pipeline();OpenGL_Begin(GL_LINE_STRIP);
        OpenGL_Vertex(sphere_center+vFrom);OpenGL_Vertex(sphere_center);OpenGL_Vertex(sphere_center+vTo);
        OpenGL_End();}}

    void DrawResultArc() const
    {/*RESCOLOR();if(vrFrom!=vrTo) DrawAnyArc(vrFrom,vrTo);*/}

    void Circ() const
    {TV p(0,1,0),m(1,0,0),n(0,0,1);
    DrawColor(outer_rim,rotation_axis,4);DrawAnyArcWorld(p,m);DrawAnyArcWorld(m,-p);DrawAnyArcWorld(-p,-m);DrawAnyArcWorld(-m,p);
    DrawColor(z_axis,rotation_axis,2);DrawAnyArc(p,m);DrawAnyArc(m,-p);DrawAnyArc(-p,-m);DrawAnyArc(-m,p);
    DrawColor(x_axis,rotation_axis,0);DrawAnyArc(p,n);DrawAnyArc(n,-p);DrawAnyArc(-p,-n);DrawAnyArc(-n,p);
    DrawColor(y_axis,rotation_axis,1);DrawAnyArc(m,n);DrawAnyArc(n,-m);DrawAnyArc(-m,-n);DrawAnyArc(-n,m);}

    TV MouseOnSphere(const VECTOR<T,2> &mouse,const VECTOR<T,2> &ballCenter,double ballRadius)
    {TV ballMouse;T mag;
    ballMouse.x=T((mouse.x-ballCenter.x)/ballRadius);
    ballMouse.y=T((mouse.y-ballCenter.y)/ballRadius);
    mag=ballMouse.Magnitude_Squared();
    if (mag>1.0){T scale=T(1.0/sqrt(mag));ballMouse.x*=scale;ballMouse.y*=scale;ballMouse.z=0.0;}
    else ballMouse.z=T(sqrt(1-mag));
    return ballMouse;}

    TV MouseOnSphere(const VECTOR<T,2> &mouse,const SPHERE<TV> &sphere,const TV &previous_result)
    {T visual_radius=sphere.radius*(sphere.center-world->Get_Camera_Position()).Magnitude();
    T added_diameter=2*visual_radius-2*sphere.radius;T thickness=visual_radius/4.;
    TV x_max((T)0.5*thickness,(T)0,(T)0),x_min((T)-0.5*thickness,(T)0,(T)0);
    TV y_max((T)0,(T)0.5*thickness,(T)0),y_min((T)0,(T)-0.5*thickness,(T)0);
    TV z_max((T)0,(T)0,(T)0.5*thickness),z_min((T)0,(T)0,(T)-0.5*thickness);
    CYLINDER<T> x_circ(qDown.Rotate(x_max)+sphere.center,qDown.Rotate(x_min)+sphere.center,visual_radius);
    CYLINDER<T> y_circ(qDown.Rotate(y_max)+sphere.center,qDown.Rotate(y_min)+sphere.center,visual_radius);
    CYLINDER<T> z_circ(qDown.Rotate(z_max)+sphere.center,qDown.Rotate(z_min)+sphere.center,visual_radius);
    PLANE<T> x_plane(qDown.Rotate(TV(1,0,0)),sphere.center);
    PLANE<T> y_plane(qDown.Rotate(TV(0,1,0)),sphere.center);
    PLANE<T> z_plane(qDown.Rotate(TV(0,0,1)),sphere.center);
    RAY<TV> ray=opengl_world.Ray_Through_Normalized_Image_Coordinate(mouse);
    if(rotation_axis==0){
        if(INTERSECTION::Intersects(ray,x_plane)){TV point=ray.Point(ray.t_max);return (point-sphere.center).Normalized()*visual_radius;}
        else{LOG::cout<<"Failure"<<std::endl;return previous_result;}}
    if(rotation_axis==1){
        if(INTERSECTION::Intersects(ray,y_plane)){TV point=ray.Point(ray.t_max);return (point-sphere.center).Normalized()*visual_radius;}
        else{LOG::cout<<"Failure"<<std::endl;return previous_result;}}
    if(rotation_axis==2){
        if(INTERSECTION::Intersects(ray,z_plane)){TV point=ray.Point(ray.t_max);return (point-sphere.center).Normalized()*visual_radius;}
        else{LOG::cout<<"Failure"<<std::endl;return previous_result;}}
    if(rotation_axis==3){
        if(INTERSECTION::Intersects(ray,sphere,added_diameter)){return (ray.Point(ray.t_max)-sphere.center);}
        else{qDown=qNow;vDown=mouse;rotation_axis=4;vFrom=MouseOnSphere(vDown,sphere,vFrom);}}
    if(rotation_axis==4){
        GLdouble project_x,project_y,project_z;
        GLdouble model_view[16];glGetDoublev(GL_MODELVIEW_MATRIX,model_view);
        GLdouble projection[16];glGetDoublev(GL_PROJECTION_MATRIX,projection);
        GLint viewport[4];glGetIntegerv(GL_VIEWPORT,viewport);
        gluUnProject((GLdouble)0,(GLdouble)0,(GLdouble)1,model_view,projection,viewport,&project_x,&project_y,&project_z);
        PLANE<T> view_plane(TV(project_x,project_y,project_z).Normalized(),sphere.center);
        if(INTERSECTION::Intersects(ray,view_plane)){TV point=ray.Point(ray.t_max);return (point-sphere.center).Normalized()*visual_radius*1.1;}
        else{PHYSBAM_FATAL_ERROR();}}
    if(rotation_axis==-1&&INTERSECTION::Intersects(ray,sphere,added_diameter)){
        TV sphere_point=ray.Point(ray.t_max);
        ray=opengl_world.Ray_Through_Normalized_Image_Coordinate(mouse);
        if(INTERSECTION::Intersects(ray,x_circ)){
            LOG::cout<<"X axis init"<<std::endl;rotation_axis=0;
            rotation_axis=0;
            TV point=ray.Point(ray.t_max);
            T magnitude=(point-sphere.center).Magnitude();
            if(magnitude<visual_radius+thickness&&magnitude>visual_radius-thickness&&TV::Dot_Product(point-sphere.center,world->Get_Camera_Position()-sphere.center)>-1e-8){
                TV closest=x_plane.Surface(point);return (closest-sphere.center).Normalized()*visual_radius;}}
        else if (INTERSECTION::Intersects(ray,y_circ)){
            LOG::cout<<"Y axis init"<<std::endl;rotation_axis=1;
            rotation_axis=1;
            TV point=ray.Point(ray.t_max);
            T magnitude=(point-sphere.center).Magnitude();
            if(magnitude<visual_radius+thickness&&magnitude>visual_radius-thickness&&TV::Dot_Product(point-sphere.center,world->Get_Camera_Position()-sphere.center)>-1e-8){
                TV closest=y_plane.Surface(point);return (closest-sphere.center).Normalized()*visual_radius;}}
        else if (INTERSECTION::Intersects(ray,z_circ)){
            LOG::cout<<"Z axis init"<<std::endl;rotation_axis=2;
            rotation_axis=2;
            TV point=ray.Point(ray.t_max);
            T magnitude=(point-sphere.center).Magnitude();
            if(magnitude<visual_radius+thickness&&magnitude>visual_radius-thickness&&TV::Dot_Product(point-sphere.center,world->Get_Camera_Position()-sphere.center)>-1e-8){
                TV closest=z_plane.Surface(point);return (closest-sphere.center).Normalized()*visual_radius;}}
        LOG::cout<<"Middle init"<<std::endl;rotation_axis=3;
        rotation_axis=3;
        return (sphere_point-sphere.center);
    }else if(INTERSECTION::Intersects(ray,sphere,added_diameter+visual_radius/(T)2)){
        ray=opengl_world.Ray_Through_Normalized_Image_Coordinate(mouse);
        LOG::cout<<"Outer Rim"<<std::endl;rotation_axis=4;
        rotation_axis=4;
        GLdouble project_x,project_y,project_z;
        GLdouble model_view[16];glGetDoublev(GL_MODELVIEW_MATRIX,model_view);
        GLdouble projection[16];glGetDoublev(GL_PROJECTION_MATRIX,projection);
        GLint viewport[4];glGetIntegerv(GL_VIEWPORT,viewport);
        gluUnProject((GLdouble)0,(GLdouble)0,(GLdouble)1,model_view,projection,viewport,&project_x,&project_y,&project_z);
        PLANE<T> view_plane(TV(project_x,project_y,project_z).Normalized(),sphere.center);
        if(INTERSECTION::Intersects(ray,view_plane)){TV point=ray.Point(ray.t_max);return (point-sphere.center).Normalized()*visual_radius*1.1;}
        else{PHYSBAM_FATAL_ERROR();}}
    return previous_result;}

    ROTATION<TV> Qt_FromBallPoints(const TV &from,const TV &to)
    {return ROTATION<TV>::From_Rotated_Vector(from,to);}

    void Qt_ToBallPoints(const ROTATION<TV> &q,TV &arcFrom,TV &arcTo)
    {arcFrom=q.Get_Axis().Unit_Orthogonal_Vector();arcTo=q.Rotate(arcFrom);}

    TV Bisect_Vectors(const TV &v1,const TV &v2) const
    {TV v=v1+v2;double normal=v.Magnitude_Squared();
    if(normal<1.0e-5) v=TV(0,0,1);
    return v.Normalized();}
};
}
#endif

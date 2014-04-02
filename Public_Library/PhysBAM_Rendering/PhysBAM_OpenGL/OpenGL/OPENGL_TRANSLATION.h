//#####################################################################
// Copyright 2008, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __H_OPENGL_TRANSLATION__
#define __H_OPENGL_TRANSLATION__

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_PLANE_INTERSECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_BOX_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>

#define FACTOR 0.95

namespace PhysBAM{

template<class T>
class OPENGL_TRANSLATION: public OPENGL_OBJECT
{
    typedef VECTOR<T,3> TV;
public:
    VECTOR<RAY<TV>,3> translations;
    TV qNow,qDown,qDrag;
    VECTOR<T,2> center,vDown;
    TV vFrom,vTo,vrFrom,vrTo;
    bool dragging;
    bool translate_opengl;
    OPENGL_WORLD *world;
    int translation_axis;
    OPENGL_COLOR x_axis,y_axis,z_axis,highlight;

    OPENGL_TRANSLATION()
        :dragging(false),translate_opengl(true),translation_axis(-1)
    {x_axis=OPENGL_COLOR::Red();y_axis=OPENGL_COLOR::Green();z_axis=OPENGL_COLOR::Blue();highlight=OPENGL_COLOR::Yellow();}

    void Reinitialize()
    {qNow=qDown=qDrag=TV();
    center=vDown=VECTOR<T,2>();
    vFrom=vTo=vrFrom=vrTo=TV();}

    void Initialize_Ray(TV origin,T length)
    {translations=VECTOR<RAY<TV>,3>(RAY<TV>(origin,TV(1,0,0),true),RAY<TV>(origin,TV(0,1,0),true),RAY<TV>(origin,TV(0,0,1),true));
    translations.x.semi_infinite=false;translations.x.t_max=length;
    translations.y.semi_infinite=false;translations.y.t_max=length;
    translations.z.semi_infinite=false;translations.z.t_max=length;}

    void Update(const VECTOR<T,2> &vNow)
    {vFrom=MouseOnTranslation(vDown,translations,vFrom);vTo=MouseOnTranslation(vNow,translations,vTo);
    if (dragging){qDrag=Qt_FromTransPoints(vFrom,vTo);qNow=qDrag+qDown;}}

    void Begin_Drag(const VECTOR<T,2> &vNow)
    {dragging=true;vDown=vNow;Update(vNow);}

    void End_Drag(const VECTOR<T,2> &vNow)
    {dragging=false;Update(vNow);qDown=qNow;translation_axis=-1;}
    
    void Display(const int color) const
    {DrawTranslation();}
    
private:
    void DrawColor(const OPENGL_COLOR &color,int translation_axis,int my_axis) const
    {if(translation_axis==my_axis) highlight.Send_To_GL_Pipeline();
    else color.Send_To_GL_Pipeline();}

    void DrawAnyRay(const RAY<TV> &ray,const OPENGL_COLOR &color,bool highlighted) const
    {TV current;if(highlighted) current=qNow;else current=qDown;
    T scale=(ray.endpoint-world->Get_Camera_Position()).Magnitude();
    OpenGL_Begin(GL_LINES);OpenGL_Line(qDown+ray.endpoint,current+ray.Point(ray.t_max)+ray.direction*scale*ray.t_max);
    TV center=current+ray.Point(ray.t_max)+ray.direction*scale*ray.t_max;T distance=(ray.Point(ray.t_max)-ray.Point(FACTOR*ray.t_max)).Magnitude()*scale;
    RANGE<TV> box(center-TV(distance,distance,distance),center+TV(distance,distance,distance));
    OPENGL_BOX_3D<T> opengl_box(box,highlighted?highlight:color);opengl_box.Display();}

    void DrawTranslation() const
    {glDisable(GL_LIGHTING);DrawColor(x_axis,translation_axis,0);DrawAnyRay(translations.x,x_axis,translation_axis==0);
    DrawColor(y_axis,translation_axis,1);DrawAnyRay(translations.y,y_axis,translation_axis==1);
    DrawColor(z_axis,translation_axis,2);DrawAnyRay(translations.z,z_axis,translation_axis==2);glEnable(GL_LIGHTING);}

    TV MouseOnTranslation(const VECTOR<T,2> &mouse,VECTOR<RAY<TV>,3> &translations,const TV &previous_result){
        RAY<TV> camera_ray=world->Ray_Through_Normalized_Image_Coordinate(mouse);
        TV down;if(translate_opengl) down=qDown;
        TV z(0,0,1),x(1,0,0),y(0,1,0);
        if(translation_axis==0){
            PLANE<T> z_plane(z,down+translations.x.endpoint);
            PLANE<T> y_plane(y,down+translations.x.endpoint);
            if(INTERSECTION::Intersects(camera_ray,z_plane)){TV point=camera_ray.Point(camera_ray.t_max);return TV(point.x,translations.x.endpoint.y+down.y,translations.x.endpoint.z+down.z);}
            else if(INTERSECTION::Intersects(camera_ray,y_plane)){TV point=camera_ray.Point(camera_ray.t_max);return TV(point.x,translations.x.endpoint.y+down.y,translations.x.endpoint.z+down.z);}
            else{LOG::cout<<"Failure"<<std::endl;return previous_result;}}
        if(translation_axis==1){
            PLANE<T> z_plane(z,down+translations.y.endpoint);
            PLANE<T> x_plane(x,down+translations.y.endpoint);
            if(INTERSECTION::Intersects(camera_ray,z_plane)){TV point=camera_ray.Point(camera_ray.t_max);return TV(translations.y.endpoint.x+down.x,point.y,translations.y.endpoint.z+down.z);}
            else if(INTERSECTION::Intersects(camera_ray,x_plane)){TV point=camera_ray.Point(camera_ray.t_max);return TV(translations.y.endpoint.x+down.x,point.y,translations.y.endpoint.z+down.z);}
            else{LOG::cout<<"Failure"<<std::endl;return previous_result;}}
        if(translation_axis==2){
            PLANE<T> y_plane(y,down+translations.z.endpoint);
            PLANE<T> x_plane(x,down+translations.z.endpoint);
            if(INTERSECTION::Intersects(camera_ray,y_plane)){TV point=camera_ray.Point(camera_ray.t_max);return TV(translations.z.endpoint.x+down.x,translations.z.endpoint.y+down.y,point.z);}
            else if(INTERSECTION::Intersects(camera_ray,x_plane)){TV point=camera_ray.Point(camera_ray.t_max);return TV(translations.z.endpoint.x+down.x,translations.z.endpoint.y+down.y,point.z);}
            else{LOG::cout<<"Failure"<<std::endl;return previous_result;}}
        if(translation_axis==-1){
            LOG::cout<<"No axis"<<std::endl;
            RAY<TV> &ray_x=translations.x;T scale=(ray_x.endpoint-world->Get_Camera_Position()).Magnitude();
            TV center=down+ray_x.Point(ray_x.t_max)+ray_x.direction*scale*ray_x.t_max;TV end=down+ray_x.Point((1-FACTOR)*ray_x.t_max)+ray_x.direction*scale*ray_x.t_max;
            T distance=(ray_x.Point(ray_x.t_max)-ray_x.Point(FACTOR*ray_x.t_max)).Magnitude()*scale;center+=distance*ray_x.direction;
            BOX<TV> box(end-TV(distance,distance,distance),center+TV(distance,distance,distance));
            if(INTERSECTION::Intersects(camera_ray,box)){
                LOG::cout<<"X axis"<<std::endl;
                translation_axis=0;
                TV point=camera_ray.Point(camera_ray.t_max);
                return TV(point.x,translations.x.endpoint.y+down.y,translations.x.endpoint.z+down.z);}
            RAY<TV> &ray_y=translations.y;scale=(ray_y.endpoint-world->Get_Camera_Position()).Magnitude();
            center=down+ray_y.Point(ray_y.t_max)+ray_y.direction*scale*ray_y.t_max;end=down+ray_y.Point((1-FACTOR)*ray_y.t_max)+ray_y.direction*scale*ray_y.t_max;
            distance=(ray_y.Point(ray_y.t_max)-ray_y.Point(FACTOR*ray_y.t_max)).Magnitude()*scale;center+=distance*ray_y.direction;
            box=BOX<TV>(end-TV(distance,distance,distance),center+TV(distance,distance,distance));
            if(INTERSECTION::Intersects(camera_ray,box)){
                LOG::cout<<"Y axis"<<std::endl;
                translation_axis=1;
                TV point=camera_ray.Point(camera_ray.t_max);
                return TV(translations.y.endpoint.x+down.x,point.y,translations.y.endpoint.z+down.z);}
            RAY<TV> &ray_z=translations.z;scale=(ray_z.endpoint-world->Get_Camera_Position()).Magnitude();
            center=down+ray_z.Point(ray_z.t_max)+ray_z.direction*scale*ray_z.t_max;end=down+ray_z.Point((1-FACTOR)*ray_z.t_max)+ray_z.direction*scale*ray_z.t_max;
            distance=(ray_z.Point(ray_z.t_max)-ray_z.Point(FACTOR*ray_z.t_max)).Magnitude()*scale;center+=distance*ray_z.direction;
            box=BOX<TV>(end-TV(distance,distance,distance),center+TV(distance,distance,distance));
            if(INTERSECTION::Intersects(camera_ray,box)){
                LOG::cout<<"Z axis"<<std::endl;
                translation_axis=2;
                TV point=camera_ray.Point(camera_ray.t_max);
                return TV(translations.z.endpoint.x+down.x,translations.z.endpoint.y+down.y,point.z);}}
    return previous_result;}

    TV Qt_FromTransPoints(const TV &from,const TV &to)
    {return (to-from);}
};
}
#endif

//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Duc Nguyen, Eran Guendelman, Geoffrey Irving, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RAY_TRACING_DRIVER_WITH_PREVIEW
//##################################################################### 
#ifndef __RAY_TRACING_DRIVER_WITH_PREVIEW__
#define __RAY_TRACING_DRIVER_WITH_PREVIEW__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include "RAY_TRACING_DRIVER.h"
namespace PhysBAM{

template<class T>
class RAY_TRACING_DRIVER_WITH_PREVIEW:public RAY_TRACING_DRIVER<T>
{
public:
    using RAY_TRACING_DRIVER<T>::example;using RAY_TRACING_DRIVER<T>::world;using RAY_TRACING_DRIVER<T>::progress;using RAY_TRACING_DRIVER<T>::large_pixel_size;

    static RAY_TRACING_DRIVER_WITH_PREVIEW<T>* driver;
    T one_over_gamma;
    int state_j;
    int state_pixel_size;

    RAY_TRACING_DRIVER_WITH_PREVIEW(RAY_TRACING_EXAMPLE<T>& example_input)
        :RAY_TRACING_DRIVER<T>(example_input)
    {
        driver=this;
    }

    virtual ~RAY_TRACING_DRIVER_WITH_PREVIEW()
    {}

    static void Handle_Display()
    {driver->Redraw_From_Film();}
    
    static void Handle_Idle()
    {driver->Render_More();}
    
    static void Handle_Visibility(int state)
    {if(state) Handle_Display();}

    VECTOR<T,3> Gamma_Correct(VECTOR<T,3>& color)
    {return VECTOR<T,3>(pow(max((T)0,color.x),one_over_gamma),pow(max((T)0,color.y),one_over_gamma),pow(max((T)0,color.z),one_over_gamma));}

//#####################################################################
    void Initialize_Preview();
    void Draw_Big_Pixel(int i,int j,const VECTOR<T,3> &color,int size);
    void Redraw_From_Film();
    void Render_More();
    void Execute_Main_Program() PHYSBAM_OVERRIDE;
//#####################################################################
};
//#####################################################################
template<class T> RAY_TRACING_DRIVER_WITH_PREVIEW<T>* RAY_TRACING_DRIVER_WITH_PREVIEW<T>::driver=0;
//#####################################################################
// Function Initialize_Preview
//#####################################################################
template<class T> void RAY_TRACING_DRIVER_WITH_PREVIEW<T>::
Initialize_Preview()
{
    int argc=1;const char *(argv[1]);argv[0]="render";
    glutInit(&argc,(char**)argv);
    glutInitDisplayMode(GLUT_SINGLE|GLUT_RGB);

    glutInitWindowSize(world.camera.film.grid.counts.x,world.camera.film.grid.counts.y);
    glutCreateWindow("Render Preview");
    glutDisplayFunc(Handle_Display);
    glutVisibilityFunc(Handle_Visibility);
    glutIdleFunc(Handle_Idle);

    glClearColor(0,0,0,0);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0,world.camera.film.grid.counts.x,0,world.camera.film.grid.counts.y);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(.375,.375,0);

    glutReshapeWindow(world.camera.film.grid.counts.x,world.camera.film.grid.counts.y);
}
//#####################################################################
// Function Draw_Big_Pixel
//#####################################################################
template<class T> void RAY_TRACING_DRIVER_WITH_PREVIEW<T>::
Draw_Big_Pixel(int i,int j,const VECTOR<T,3>& color,int size)
{
    glColor3d(color.x,color.y,color.z);
    glVertex2i(i-size,j-size);glVertex2i(i-size,j);glVertex2i(i,j);glVertex2i(i,j-size);
}
//#####################################################################
// Function Redraw_From_Film
//#####################################################################
template<class T> void RAY_TRACING_DRIVER_WITH_PREVIEW<T>::
Redraw_From_Film()
{
    // Regenerate image (using big pixels) from what we've already drawn to film
    glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_QUADS);
    RANGE<VECTOR<int,2> >& clip=example.clipping_region;
    for(int pixel_size=large_pixel_size;pixel_size>=state_pixel_size;pixel_size/=2)
        for(int j=clip.min_corner.y;j<=clip.max_corner.y;j++){
            if(j%pixel_size) continue;
            if(pixel_size==state_pixel_size && j>=state_j) break;
            for(int i=clip.min_corner.x;i<=clip.max_corner.x;i++){
                if(i%pixel_size || (pixel_size!=large_pixel_size && i%(pixel_size*2)==0 && j%(pixel_size*2)==0)) continue;
                Draw_Big_Pixel(i,j,Gamma_Correct(world.camera.film.colors(i,j)),pixel_size);}}
    glEnd();
    glFlush();
}
//#####################################################################
// Function Render_More
//#####################################################################
template<class T> void RAY_TRACING_DRIVER_WITH_PREVIEW<T>::
Render_More()
{
    RANGE<VECTOR<int,2> >& clip=example.clipping_region;
    while(state_j%state_pixel_size)state_j++;
    if(state_j>clip.max_corner.y){
        if(state_pixel_size==1){
            LOG::Time("Writing Images");
            world.camera.film.Print_Film_Clipped(example.Get_Output_Filename(example.frame),example.gamma_correction,clip);
            LOG::cout<<std::flush;
            exit(0);}
        state_j=clip.min_corner.y;state_pixel_size/=2; // next pass
        while(state_j%state_pixel_size)state_j++;}
    for(int i=clip.min_corner.x;i<=clip.max_corner.x;i++){
        if(i%state_pixel_size || (state_pixel_size!=large_pixel_size && state_j%(state_pixel_size*2)==0 && i%(state_pixel_size*2)==0)) continue;
        VECTOR<int,2> pixel_index(i,state_j);
        world.Render_Pixel(pixel_index);
        VECTOR<T,3> color=world.camera.film.Pixel_Color(pixel_index);
        glBegin(GL_QUADS);Draw_Big_Pixel(i,state_j,Gamma_Correct(color),state_pixel_size);glEnd();
        progress.Progress();}
    glFlush();
    state_j+=state_pixel_size;
}
//#####################################################################
// Function Execute_Main_Program
//#####################################################################
template<class T> void RAY_TRACING_DRIVER_WITH_PREVIEW<T>::
Execute_Main_Program()
{
    LOG::SCOPE scope("RENDER","Rendering frame %d",example.frame);
    RAY_TRACING_DRIVER<T>::Initialize();
    LOG::Time("Forward Ray Trace");
    one_over_gamma=(T)1/example.gamma_correction;
    large_pixel_size=1<<5;
    state_j=1;state_pixel_size=large_pixel_size;
    
    Initialize_Preview();
    glutMainLoop();
}       
//#####################################################################
}
#endif

//#####################################################################
// Copyright 2004, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMAGE_WINDOW
//#####################################################################
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Math_Tools/integer_log.h>
#include <fstream>
#include "IMAGE_WINDOW.h"
#include <Fl/Fl.h>
using namespace PhysBAM;
//#####################################################################
// IMAGE_WINDOW
//#####################################################################
template <class T> IMAGE_WINDOW<T>::
IMAGE_WINDOW(int x,int y,int width,int height,int picture_width,int picture_height)
:Fl_Gl_Window(x,y,width,height,"Image Window"),view_x_offset(0),view_y_offset(0),view_scale(1),view_scale_power(0),user_callback(0),have_last_click(false),texture_initialized(false)
{
    image.Resize(1,picture_width,1,picture_height);
    for(int i=1;i<=image.m;i++)for(int j=1;j<=image.n;j++)image(i,j)=VECTOR<T,3>((T)0.5,(T)0.3,(T)0.3);
}
//#####################################################################
// Init_Gl
//#####################################################################
template<class T> void IMAGE_WINDOW<T>::
Init_Gl()
{
    glEnable(GL_LINE_SMOOTH);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA) ;    
    glClearColor(0.5,0.5,0.5,1);
    if(!texture_initialized){
        glGenTextures(1,&image_texture);
        glBindTexture(GL_TEXTURE_2D,image_texture);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
        // compute the power of two square size for the texture
        int power_of_two_width=integer_log(image.m),power_of_two_height=integer_log(image.n);
        if(1<<power_of_two_width!=image.m)power_of_two_width++;
        if(1<<power_of_two_height!=image.n)power_of_two_height++;
        power_of_two_width=1<<power_of_two_width;power_of_two_height=1<<power_of_two_height;
        texture_dimension=max(power_of_two_width,power_of_two_height);
        int array_size=4*texture_dimension*texture_dimension;
        //printf("Image %d by %d --> texture %d by %d\n",image.m,image.n,texture_dimension,texture_dimension);
        float* dummy_image=new float[array_size];memset(dummy_image,0,array_size*sizeof(float));
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texture_dimension, texture_dimension, 0, GL_RGBA, GL_FLOAT, dummy_image);
        delete [] dummy_image;texture_initialized=true;}
    glViewport(0,0,w(),h());
    Print_Gl_Error();
}
//#####################################################################
// Set_Callback
//#####################################################################
template<class T> void IMAGE_WINDOW<T>::
Set_Callback(void (*f)(void*,int,int),void *data)
{
    user_callback=f;user_callback_data=data;
}
//#####################################################################
// hanLoad_Cameradle
//#####################################################################
template<class T> void IMAGE_WINDOW<T>::
Load_Camera()
{
    glLoadIdentity();
    glScalef((GLfloat)view_scale,(GLfloat)view_scale,(GLfloat)view_scale);
    gluOrtho2D(1,w(),1,h());
    glTranslatef((GLfloat)view_x_offset,(GLfloat)view_y_offset,0);
}
//#####################################################################
// handle
//#####################################################################
template <class T> int IMAGE_WINDOW<T>::
handle(int event)
{
    static bool translating=false;static bool scaling=false;static int origin_x;static int origin_y;
    if(event==FL_MOUSEWHEEL){
        if(Fl::event_dy()>0)view_scale*=0.75;else view_scale*=1.25;
        damage(1);
        return 1;
    }else if(event==FL_PUSH){
        if(Fl::event_button()==1){
            make_current();
            Load_Camera();
            GLdouble projection[16];glGetDoublev(GL_PROJECTION_MATRIX, projection);    
            GLdouble modelview[16];glGetDoublev(GL_MODELVIEW_MATRIX, modelview);    
            GLint viewport[4];glGetIntegerv(GL_VIEWPORT, viewport);    
            GLdouble coord_x,coord_y,coord_z;
            int xx=Fl::event_x();int yy=h()-Fl::event_y()-1;
            gluUnProject(xx,yy,0,modelview,projection,viewport,&coord_x,&coord_y,&coord_z);
            coord_x++;
            int pixel_x=(int)(coord_x-0.5),pixel_y=(int)(coord_y+0.5);
            last_click=VECTOR<int,2>(pixel_x,pixel_y);have_last_click=true;
            if(user_callback)user_callback(user_callback_data,pixel_x,pixel_y);
            damage(1);
        }
        else if(Fl::event_button()==2)translating=true;
        else if(Fl::event_button()==3)scaling=true;
        origin_x=Fl::event_x();origin_y=Fl::event_y();
        return 1;
    }else if(event==FL_DRAG){
        if(Fl::event_button()==2&&translating){view_x_offset+=(Fl::event_x()-origin_x)/view_scale;view_y_offset-=(Fl::event_y()-origin_y)/view_scale;}
        else if(Fl::event_button()==3&&scaling){
            T dy=T(origin_y-Fl::event_y());
            view_scale_power+=T(dy*0.01);
            view_scale=pow(T(2),view_scale_power);}
        origin_x=Fl::event_x();origin_y=Fl::event_y();
        damage(1);
        return 1;
    }else if(event==FL_RELEASE){
        translating=false;scaling=false;
        return 1;
    }
    return Fl_Gl_Window::handle(event);
}
//#####################################################################
// draw
//#####################################################################
template <class T> void IMAGE_WINDOW<T>::
draw()
{
    if(!valid())Init_Gl();
    if(update_image){
        float *temp_image=new float[image.m*image.n*3];memset(temp_image,0,image.m*image.n*3);
        for(int i=1;i<=image.m;i++)for(int j=1;j<=image.n;j++)for(int k=1;k<=3;k++)temp_image[3*((j-1)*image.m+(i-1))+k-1]=float(image(i,j)[k]);
        glTexSubImage2D(GL_TEXTURE_2D,0,0,0,image.m,image.n,GL_RGB,GL_FLOAT,temp_image);
        delete temp_image;
        update_image=false;
    }
    
    Load_Camera();
    
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,image_texture);
    glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
    glBegin(GL_QUADS);
    glColor3f(1,1,0);
    glTexCoord2f(0,0);glVertex2f(0.5,0.5);
    glTexCoord2f(0,1);glVertex2f(0.5,(GLfloat)texture_dimension+0.5);
    glTexCoord2f(1,1);glVertex2f((GLfloat)texture_dimension+0.5,(GLfloat)texture_dimension+0.5);
    glTexCoord2f(1,0);glVertex2f((GLfloat)texture_dimension+0.5,0.5);
    /*glTexCoord2f(0,0);glVertex2f(1,1);
      glTexCoord2f(0,1);glVertex2f(1,100);
      glTexCoord2f(1,1);glVertex2f(100,100);
      glTexCoord2f(1,0);glVertex2f(100,1);*/
    glEnd();
    glDisable(GL_TEXTURE_2D);
    if(have_last_click){
        int i=last_click.x,j=last_click.y;
        glLineWidth(4.0);
        glBegin(GL_LINES);
        glVertex2f(i-0.5,j-0.5);glVertex2f(i-0.5,j+0.5);glVertex2f(i+0.5,j-0.5);glVertex2f(i+0.5,j+0.5);
        glVertex2f(i-0.5,j-0.5);glVertex2f(i+0.5,j-0.5);glVertex2f(i-0.5,j+0.5);glVertex2f(i+0.5,j+0.5);
        //glVertex2i(i,j-10);glVertex2i(i,j+10);
        glEnd();
    }
    Print_Gl_Error();
    
}
//#####################################################################
// Update_Image
//#####################################################################
template <class T> void IMAGE_WINDOW<T>::
Update_Image()
{
    update_image=true;
    damage(1);
}
//#####################################################################
// Set_Color
//#####################################################################
template <class T> void IMAGE_WINDOW<T>::
Set_Color(int i,int j,const VECTOR<T,3>& color) 
{
    image(i,j)=color;
}
//#####################################################################
// Print_Gl_Error
//#####################################################################
template <class T> void IMAGE_WINDOW<T>::
Print_Gl_Error() 
{
    GLenum errCode;
    if((errCode=glGetError())!=GL_NO_ERROR){
        fprintf(stderr,"OpenGL Error: %s\n",gluErrorString(errCode));
    }
}
//#####################################################################
// Save_RGB
//#####################################################################
template <class T> void IMAGE_WINDOW<T>::
Save_RGB(const char *filename) 
{
    IMAGE<T>::Write(filename,image);
}
//#####################################################################
// Save_PhysBAM
//#####################################################################
template <class T> void IMAGE_WINDOW<T>::
Save_PhysBAM(const char *filename) 
{
    Save_RGB(filename); // broke
//    std::ofstream output;
//    output.open(filename,std::ios::binary);image.template Write<float>(output);output.close();
}
template class IMAGE_WINDOW<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class IMAGE_WINDOW<double>;
#endif

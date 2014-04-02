//#####################################################################
// Class OPENGL_LIGHT
//##################################################################### 
//
//#####################################################################
// Neverov - January 30, 2002
// Bridson - February 3, 2002
//#####################################################################
#ifndef __OPENGL_LIGHT__
#define __OPENGL_LIGHT__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
namespace PhysBAM{

class OPENGL_LIGHT
{
public:
    GLfloat homogeneous_position[4];
    OPENGL_COLOR color;

    OPENGL_LIGHT(const VECTOR<double,3>& position,const float value=1,const int finite_position=0)
        :color(OPENGL_COLOR(value,value,value))
    {homogeneous_position[0]=position.x;homogeneous_position[1]=position.y;homogeneous_position[2]=position.z;homogeneous_position[3]=(GLfloat)finite_position;}

    OPENGL_LIGHT(const VECTOR<double,3>& position,const OPENGL_COLOR& color_input,const int finite_position=0)
        :color(color_input)
    {homogeneous_position[0]=position.x;homogeneous_position[1]=position.y;homogeneous_position[2]=position.z;homogeneous_position[3]=(GLfloat)finite_position;}

    void Send_To_GL_Pipeline(const int id)
    {
        GLenum gl_light_id=GL_LIGHT0+id;
        glLightfv(gl_light_id,GL_DIFFUSE,color.rgba);
        glLightfv(gl_light_id,GL_SPECULAR,color.rgba);
        glLightfv(gl_light_id,GL_POSITION,homogeneous_position);
        glEnable(gl_light_id);
    }
};
}
#endif


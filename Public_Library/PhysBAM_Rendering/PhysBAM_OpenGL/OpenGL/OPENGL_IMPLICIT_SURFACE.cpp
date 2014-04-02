//#####################################################################
// Copyright 2002, Robert Bridson, Eilene Hao, Geoffrey Irving, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_IMPLICIT_SURFACE
//##################################################################### 
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_IMPLICIT_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_IMPLICIT_SURFACE<T>::
OPENGL_IMPLICIT_SURFACE(IMPLICIT_OBJECT<TV>& surface_input,const OPENGL_MATERIAL& material_input)
    :surface(surface_input),two_sided(false),front_material(material_input),front_material_gray(material_input.Grayscale()),nx(32),ny(24),nz(100),slice1(0),slice2(0),up_to_date1(0),
    up_to_date2(0)
{
    Initialize_Slices();
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_IMPLICIT_SURFACE<T>::
OPENGL_IMPLICIT_SURFACE(IMPLICIT_OBJECT<TV>& surface_input,const OPENGL_MATERIAL &front_material_input,const OPENGL_MATERIAL& back_material_input)
    :surface(surface_input),two_sided(true),front_material(front_material_input),front_material_gray(front_material_input.Grayscale()),back_material(back_material_input),
    back_material_gray(back_material_input.Grayscale()),nx(32),ny(24),nz(100),slice1(0),slice2(0),up_to_date1(0),up_to_date2(0)
{
    Initialize_Slices();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_IMPLICIT_SURFACE<T>::
~OPENGL_IMPLICIT_SURFACE()
{
    delete[] slice1;
    delete[] slice2;
    delete[] up_to_date1;
    delete[] up_to_date2;
}
//#####################################################################
// Function Set_Resolution
//#####################################################################
template<class T> void OPENGL_IMPLICIT_SURFACE<T>::
Set_Resolution(const int nx_input,const int ny_input,const int nz_input)
{
    nx=nx_input;ny=ny_input;nz=nz_input;
    delete[] slice1;delete[] slice2;delete[] up_to_date1;delete[] up_to_date2;
    Initialize_Slices();
}
//#####################################################################
// Function Initialize_Slices
//#####################################################################
template<class T> void OPENGL_IMPLICIT_SURFACE<T>::
Initialize_Slices()
{
    int size=(ny+1)*(nz+1);
    slice1=new float[size];
    slice2=new float[size];
    up_to_date1=new unsigned int[size];
    up_to_date2=new unsigned int[size];
    memset(up_to_date1,0,sizeof(float)*size);
    memset(up_to_date2,0,sizeof(float)*size);
    up_to_date_counter=1;
}
//#####################################################################
// Function Phi
//#####################################################################
template<class T> float OPENGL_IMPLICIT_SURFACE<T>::
Phi(VECTOR<T,3> x) const
{
    if (surface.box.Outside(x)) return max(x.x-surface.box.xmax,surface.box.xmin-x.x,x.y-surface.box.ymax,surface.box.ymin-x.y,x.z-surface.box.zmax,surface.box.zmin-x.z);
    else return surface(x);
}
//#####################################################################
// Function Display_Tetrahedron
//#####################################################################
template<class T>
void OPENGL_IMPLICIT_SURFACE<T>::
Display_Tetrahedron(const VECTOR<T,3>& x1,const VECTOR<T,3>& x2,const VECTOR<T,3>& x3,const VECTOR<T,3>& x4, 
                    const float phi1,const float phi2,const float phi3,const float phi4) const
{
    VECTOR<T,3> p1,p2,p3,p4;
    int positive_count=(phi1>0)+(phi2>0)+(phi3>0)+(phi4>0);
    switch(positive_count){
        case 0:
            // no triangles
            break;
        case 1:
            // one triangle
            if(phi1 > 0){
                p1=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x1,x2,phi1/(phi1-phi2));
                p2=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x1,x4,phi1/(phi1-phi4));
                p3=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x1,x3,phi1/(phi1-phi3));} 
            else if(phi2 > 0){
                p1=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x2,x1,phi2/(phi2-phi1));
                p2=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x2,x3,phi2/(phi2-phi3));
                p3=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x2,x4,phi2/(phi2-phi4));}
            else if(phi3 > 0){
                p1=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x3,x1,phi3/(phi3-phi1));
                p2=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x3,x4,phi3/(phi3-phi4));
                p3=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x3,x2,phi3/(phi3-phi2));}
            else{
                p1=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x4,x1,phi4/(phi4-phi1));
                p2=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x4,x2,phi4/(phi4-phi2));
                p3=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x4,x3, phi4/(phi4-phi3));} 
            OpenGL_Normal(surface.Normal(p1));OpenGL_Vertex(p1);OpenGL_Normal(surface.Normal(p2));OpenGL_Vertex(p2);OpenGL_Normal(surface.Normal(p3));OpenGL_Vertex(p3);
            break;
        case 2:
            // two triangles
            if(phi1 > 0 && phi2 > 0){
                p1=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x1,x3,phi1/(phi1-phi3));
                p2=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x1,x4,phi1/(phi1-phi4));
                p3=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x2,x3,phi2/(phi2-phi3));
                p4=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x2,x4,phi2/(phi2-phi4));}
            else if(phi1 > 0 && phi3 > 0){
                p1=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x1,x4,phi1/(phi1-phi4));
                p2=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x1,x2,phi1/(phi1-phi2));
                p3=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x3,x4,phi3/(phi3-phi4));
                p4=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x3,x2,phi3/(phi3-phi2));}
            else if(phi1 > 0 && phi4 > 0){
                p1=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x1,x2,phi1/(phi1-phi2));
                p2=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x1,x3,phi1/(phi1-phi3));
                p3=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x4,x2,phi4/(phi4-phi2));
                p4=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x4,x3,phi4/(phi4-phi3));}
            else if(phi2 > 0 && phi3 > 0){
                p1=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x2,x1,phi2/(phi2-phi1));
                p2=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x2,x4,phi2/(phi2-phi4));
                p3=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x3,x1,phi3/(phi3-phi1));
                p4=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x3,x4,phi3/(phi3-phi4));}
            else if(phi2 > 0 && phi4 > 0){
                p1=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x2,x3,phi2/(phi2-phi3));
                p2=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x2,x1,phi2/(phi2-phi1));
                p3=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x4,x3,phi4/(phi4-phi3));
                p4=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x4,x1,phi4/(phi4-phi1));}
            else{
                p1=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x3,x1,phi3/(phi3-phi1));
                p2=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x3,x2,phi3/(phi3-phi2));
                p3=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x4,x1,phi4/(phi4-phi1));
                p4=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x4,x2,phi4/(phi4-phi2));}
            OpenGL_Normal(surface.Normal(p1));OpenGL_Vertex(p1);OpenGL_Normal(surface.Normal(p3));OpenGL_Vertex(p3);OpenGL_Normal(surface.Normal(p2));OpenGL_Vertex(p2);
            OpenGL_Normal(surface.Normal(p2));OpenGL_Vertex(p2);OpenGL_Normal(surface.Normal(p3));OpenGL_Vertex(p3);OpenGL_Normal(surface.Normal(p4));OpenGL_Vertex(p4);
            break;
        case 3:
            // one triangle
            if(phi1 <= 0){
                p1=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x1,x2,phi1/(phi1-phi2));
                p2=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x1,x3,phi1/(phi1-phi3));
                p3=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x1,x4,phi1/(phi1-phi4));} 
            else if(phi2 <= 0){
                p1=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x2,x1,phi2/(phi2-phi1));
                p2=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x2,x4,phi2/(phi2-phi4));
                p3=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x2,x3,phi2/(phi2-phi3));}
            else if(phi3 <= 0){
                p1=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x3,x1,phi3/(phi3-phi1));
                p2=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x3,x2,phi3/(phi3-phi2));
                p3=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x3,x4,phi3/(phi3-phi4));}
            else{
                p1=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x4,x1,phi4/(phi4-phi1));
                p2=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x4,x3,phi4/(phi4-phi3));
                p3=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(x4,x2,phi4/(phi4-phi2));}  
            OpenGL_Normal(surface.Normal(p1));OpenGL_Vertex(p1);OpenGL_Normal(surface.Normal(p2));OpenGL_Vertex(p2);OpenGL_Normal(surface.Normal(p3));OpenGL_Vertex(p3);
            break;
        case 4:
            // no triangles
            break;
        default:
            LOG::cerr<<"WE HAD A PROBLEM.  THERE SHOULD NOT BE MORE THAN 4 POSITIVE PHI TETRAHEDRON VERTICES"<<std::endl;
            break;
    }

}
//#####################################################################
// Function Display_Brick
//#####################################################################
template<class T>
void OPENGL_IMPLICIT_SURFACE<T>::
Display_Brick(const VECTOR<T,3>& x1,float p1,const VECTOR<T,3>& x2,float p2,const VECTOR<T,3>& x3,float p3,const VECTOR<T,3>& x4,float p4,
              const VECTOR<T,3>& x5,float p5,const VECTOR<T,3>& x6,float p6,const VECTOR<T,3>& x7,float p7,const VECTOR<T,3>& x8,float p8,const int parity) const
{
    // is there a way to avoid computing normals redundantly (but not compute unnecessary ones as well)?
    if (!parity) { //means (i+j+k)%2==0
        Display_Tetrahedron(x1,x2,x3,x5,p1,p2,p3,p5);  //bottom left
        Display_Tetrahedron(x2,x5,x6,x8,p2,p5,p6,p8);  //bottom right
        Display_Tetrahedron(x2,x3,x8,x4,p2,p3,p8,p4);  //top towards
        Display_Tetrahedron(x3,x5,x8,x7,p3,p5,p8,p7);  //top away
        Display_Tetrahedron(x2,x3,x5,x8,p2,p3,p5,p8);  //center tetrahedron
    }
    else { // means (i+j+k)%2==1
        Display_Tetrahedron(x1,x2,x4,x6,p1,p2,p4,p6); //bottom towards
        Display_Tetrahedron(x1,x5,x6,x7,p1,p5,p6,p7); //bottom away
        Display_Tetrahedron(x1,x3,x7,x4,p1,p3,p7,p4); //top left
        Display_Tetrahedron(x4,x6,x8,x7,p4,p6,p8,p7); //top right
        Display_Tetrahedron(x1,x4,x7,x6,p1,p4,p7,p6); //center tetrahedron
    }
}
//#####################################################################
// Function Display
//#####################################################################
template<class T>
void OPENGL_IMPLICIT_SURFACE<T>::
Display(const int in_color) const
{
    if(two_sided){
        glDisable(GL_CULL_FACE);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,1);
        if(in_color) {front_material.Send_To_GL_Pipeline(GL_FRONT);back_material.Send_To_GL_Pipeline(GL_BACK);}
        else {front_material_gray.Send_To_GL_Pipeline(GL_FRONT);back_material_gray.Send_To_GL_Pipeline(GL_BACK);}}
    else{
        if(in_color) front_material.Send_To_GL_Pipeline();
        else front_material_gray.Send_To_GL_Pipeline();}

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    GLfloat p[16],m[16];
    glGetFloatv(GL_PROJECTION_MATRIX,p);
    glGetFloatv(GL_MODELVIEW_MATRIX,m);
    MATRIX<T,4> transform=MATRIX<T,4>(p[0]*m[0]+p[4]*m[1]+p[8]*m[2]+p[12]*m[3],
        p[1]*m[0]+p[5]*m[1]+p[9]*m[2]+p[13]*m[3],
        p[2]*m[0]+p[6]*m[1]+p[10]*m[2]+p[14]*m[3],
        p[3]*m[0]+p[7]*m[1]+p[11]*m[2]+p[15]*m[3],
        p[0]*m[4]+p[4]*m[5]+p[8]*m[6]+p[12]*m[7],
        p[1]*m[4]+p[5]*m[5]+p[9]*m[6]+p[13]*m[7],
        p[2]*m[4]+p[6]*m[5]+p[10]*m[6]+p[14]*m[7],
        p[3]*m[4]+p[7]*m[5]+p[11]*m[6]+p[15]*m[7],
        p[0]*m[8]+p[4]*m[9]+p[8]*m[10]+p[12]*m[11],
        p[1]*m[8]+p[5]*m[9]+p[9]*m[10]+p[13]*m[11],
        p[2]*m[8]+p[6]*m[9]+p[10]*m[10]+p[14]*m[11],
        p[3]*m[8]+p[7]*m[9]+p[11]*m[10]+p[15]*m[11],
        p[0]*m[12]+p[4]*m[13]+p[8]*m[14]+p[12]*m[15],
        p[1]*m[12]+p[5]*m[13]+p[9]*m[14]+p[13]*m[15],
        p[2]*m[12]+p[6]*m[13]+p[10]*m[14]+p[14]*m[15],
        p[3]*m[12]+p[7]*m[13]+p[11]*m[14]+p[15]*m[15]).Inverse(); // the OpenGL eye to world transform
    GRID<TV> eye_grid(nx,ny,nz,-1,1,-1,1,-1,1); // set up a 3d grid in eye space

    OpenGL_Begin(GL_TRIANGLES);
    int old_index,index;
    bool seen_outside,seen_inside,old_sign_change,sign_change;
    int i,j,k,skip;
    int bricks=0;
    VECTOR<T,3> old_x_stupid_declaration[4],x_stupid_declaration[4];
    VECTOR<T,3> *old_x=old_x_stupid_declaration,*x=x_stupid_declaration;
    double skip_length;
    for(i=1; i<nx; i++) { // loop over the vertical slices
        for(j=1; j<ny; j++){
            VECTOR<T,3> direction=(transform*eye_grid.X(i,j,nz)-transform*eye_grid.X(i,j,1)).Normalized();
            k=1;
            seen_outside=seen_inside=false;
            ++bricks;
            old_x[0]=transform*eye_grid.X(i,j,1);
            old_x[1]=transform*eye_grid.X(i+1,j,1);
            old_x[2]=transform*eye_grid.X(i,j+1,1);
            old_x[3]=transform*eye_grid.X(i+1,j+1,1);
            old_index=j-1;
            if(up_to_date1[old_index]  !=up_to_date_counter)   {slice1[old_index]=  Phi(old_x[0]);up_to_date1[old_index]=  up_to_date_counter;}
            if(up_to_date2[old_index]  !=up_to_date_counter+1) {slice2[old_index]=  Phi(old_x[1]);up_to_date2[old_index]=  up_to_date_counter+1;}
            if(up_to_date1[old_index+1]!=up_to_date_counter)   {slice1[old_index+1]=Phi(old_x[2]);up_to_date1[old_index+1]=up_to_date_counter;}
            if(up_to_date2[old_index+1]!=up_to_date_counter+1) {slice2[old_index+1]=Phi(old_x[3]);up_to_date2[old_index+1]=up_to_date_counter+1;}
            if(slice1[old_index]>0 && slice2[old_index]>0 && slice1[old_index+1]>0 && slice2[old_index+1]>0){
                old_sign_change=false; seen_outside=true;
                skip_length=min(min(slice1[old_index],slice2[old_index]),min(slice1[old_index+1],slice2[old_index+1]));
                skip=1;while(VECTOR<T,3>::Dot_Product(direction,transform*eye_grid.X(i,j,k+skip)-old_x[0])<skip_length) ++skip;}
            else if(slice1[old_index]<=0 && slice2[old_index]<=0 && slice1[old_index+1]<=0 && slice2[old_index+1]<=0){
                old_sign_change=false; seen_inside=true;
                skip_length=-max(max(slice1[old_index],slice2[old_index]),max(slice1[old_index+1],slice2[old_index+1]));
                skip=1;while(VECTOR<T,3>::Dot_Product(direction,transform*eye_grid.X(i,j,k+skip)-old_x[0])<skip_length) ++skip;}
            else {old_sign_change=true; skip=1;}
            k+=skip;
            for(;;){
                ++bricks;
                x[0]=transform*eye_grid.X(i,j,k);
                x[1]=transform*eye_grid.X(i+1,j,k);
                x[2]=transform*eye_grid.X(i,j+1,k);
                x[3]=transform*eye_grid.X(i+1,j+1,k);
                index=(ny+1)*(k-1)+(j-1);
                if(up_to_date1[index]  !=up_to_date_counter)   {slice1[index]=  Phi(x[0]);up_to_date1[index]=  up_to_date_counter;}
                if(up_to_date2[index]  !=up_to_date_counter+1) {slice2[index]=  Phi(x[1]);up_to_date2[index]=  up_to_date_counter+1;}
                if(up_to_date1[index+1]!=up_to_date_counter)   {slice1[index+1]=Phi(x[2]);up_to_date1[index+1]=up_to_date_counter;}
                if(up_to_date2[index+1]!=up_to_date_counter+1) {slice2[index+1]=Phi(x[3]);up_to_date2[index+1]=up_to_date_counter+1;}
                if(slice1[index]>0 && slice2[index]>0 && slice1[index+1]>0 && slice2[index+1]>0){
                    sign_change=false; seen_outside=true;
                    skip_length=min(min(slice1[index],slice2[index]),min(slice1[index+1],slice2[index+1]));
                    skip=1;while(VECTOR<T,3>::Dot_Product(direction,transform*eye_grid.X(i,j,k+skip)-old_x[0])<skip_length) ++skip;}
                else if(slice1[index]<=0 && slice2[index]<=0 && slice1[index+1]<=0 && slice2[index+1]<=0){
                    sign_change=false; seen_inside=true;
                    skip_length=-max(max(slice1[index],slice2[index]),max(slice1[index+1],slice2[index+1]));
                    skip=1;while(VECTOR<T,3>::Dot_Product(direction,transform*eye_grid.X(i,j,k+skip)-old_x[0])<skip_length) ++skip;}
                else {sign_change=true; skip=1;}
                if(sign_change || old_sign_change || (seen_inside && seen_outside)){
                    assert(skip==1);
                    Display_Brick(old_x[0],slice1[old_index],old_x[1],slice2[old_index],old_x[2],slice1[old_index+1],old_x[3],slice2[old_index+1],
                                  x[0],slice1[index],x[1],slice2[index],x[2],slice1[index+1],x[3],slice2[index+1],(i+j+k)%2);
                }
                if(seen_inside && seen_outside) break; // occlusion culling
                k+=skip;
                if(k>nz) break;
                old_index=index;
                old_sign_change=sign_change;
                exchange(old_x,x);
            }
        }
        exchange(slice1,slice2);exchange(up_to_date1,up_to_date2); // rename slice2 as slice1 for the next iteration
        ++up_to_date_counter;
        if(up_to_date_counter+1==0){
            memset(up_to_date1,0,sizeof(float)*(ny+1)*(nz+1));
            memset(up_to_date2,0,sizeof(float)*(ny+1)*(nz+1));
            up_to_date_counter=1;}
    }
    OpenGL_End();
    ++up_to_date_counter;
    if(up_to_date_counter+1==0){
        memset(up_to_date1,0,sizeof(float)*(ny+1)*(nz+1));
        memset(up_to_date2,0,sizeof(float)*(ny+1)*(nz+1));
        up_to_date_counter=1;}

    glPopMatrix();
    if(two_sided) {glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,0);glEnable(GL_CULL_FACE);}
}


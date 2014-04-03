//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Michael Lentine, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/DUALCONTOUR_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_LEVELSET_MULTIVIEW.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
using namespace PhysBAM;
//#####################################################################
// Function Reset
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Reset()
{
levelset_filename="";
triangulated_surface_filename="";
generate_triangulated_surface=false;
write_generated_triangulated_surface=false;
if(i_own_levelset && levelset){delete &levelset->grid;delete &levelset->phi;delete levelset;}
levelset=0;
delete levelset_implicit_surface;
levelset_implicit_surface=0;
if(i_own_triangulated_surface) delete triangulated_surface;triangulated_surface=0;
delete opengl_triangulated_surface;opengl_triangulated_surface=0;
delete opengl_scalar_field;opengl_scalar_field=0;
}
//#####################################################################
// Function Reset_Surface
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Reset_Surface()
{
triangulated_surface_filename="";
if(i_own_triangulated_surface) delete triangulated_surface;triangulated_surface=0;
delete opengl_triangulated_surface;opengl_triangulated_surface=0;
}
//#####################################################################
// Function Set_Levelset
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Set_Levelset(LEVELSET_3D<GRID<TV> >& levelset_input)
{
Reset();
levelset=&levelset_input;
i_own_levelset=false;
}
//#####################################################################
// Function Read_Levelset
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Read_Levelset(const std::string& levelset_filename_input)
{
Reset();
levelset_filename=levelset_filename_input;
}
//#####################################################################
// Function Levelset
//#####################################################################
template<class T,class RW> const LEVELSET_3D<GRID<VECTOR<T,3> > >* OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Levelset() const
{
return levelset;
}
//#####################################################################
// Function Set_Triangulated_Surface
//#####################################################################
template<class T,class RW> const TRIANGULATED_SURFACE<T>* OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Get_Triangulated_Surface() const
{
return triangulated_surface;
}
//#####################################################################
// Function Set_Triangulated_Surface
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Set_Triangulated_Surface(TRIANGULATED_SURFACE<T> &triangulated_surface_input)
{
PHYSBAM_ASSERT(!triangulated_surface);
triangulated_surface=&triangulated_surface_input;
i_own_triangulated_surface=false;
}
//#####################################################################
// Function Read_Triangulated_Surface
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Read_Triangulated_Surface(const std::string& triangulated_surface_filename_input)
{
PHYSBAM_ASSERT(!triangulated_surface);
triangulated_surface_filename=triangulated_surface_filename_input;
}
//#####################################################################
// Function Generate_Triangulated_Surface
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Generate_Triangulated_Surface(bool write_generated_triangulated_surface_input,const std::string& triangulated_surface_filename_input)
{
PHYSBAM_ASSERT(!triangulated_surface);
generate_triangulated_surface=true;
write_generated_triangulated_surface=write_generated_triangulated_surface_input;
if(write_generated_triangulated_surface) triangulated_surface_filename=triangulated_surface_filename_input;
}
//#####################################################################
// Function Initialize_Levelset
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Initialize_Levelset()
{
if(!levelset){
    if(levelset_filename.length() > 0){
        levelset=new LEVELSET_3D<GRID<TV> >(*(new GRID<TV>),*(new ARRAY<T,VECTOR<int,3> >));
        FILE_UTILITIES::Read_From_File<RW>(levelset_filename,*levelset);
            i_own_levelset=true;}}
}
//#####################################################################
// Function Initialize_Triangulated_Surface
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Initialize_Triangulated_Surface()
{
    if(!triangulated_surface){
        if(generate_triangulated_surface){
            Initialize_Levelset();
            if(levelset){
                if(!use_marching_tetrahedra)
                    triangulated_surface=DUALCONTOUR_3D<T>::Create_Triangulated_Surface_From_Levelset(*levelset);
                else{
                    triangulated_surface=TRIANGULATED_SURFACE<T>::Create();
                    levelset->Calculate_Triangulated_Surface_From_Marching_Tetrahedra(*triangulated_surface,true);levelset->Compute_Normals();
                    triangulated_surface->Use_Vertex_Normals();triangulated_surface->vertex_normals=new ARRAY<TV>(triangulated_surface->particles.array_collection->Size());
                    for(int p=1;p<=triangulated_surface->particles.array_collection->Size();p++)
                        (*triangulated_surface->vertex_normals)(p)=levelset->Normal(triangulated_surface->particles.X(p));}
                if(write_generated_triangulated_surface && triangulated_surface_filename.length() > 0){
                    const int vertex_normals_length=1; // needed for backwards compatibility, since vertex_normals used to be ARRAYS<TV,1>
                    FILE_UTILITIES::Write_To_File<RW>(triangulated_surface_filename,*triangulated_surface,vertex_normals_length,*triangulated_surface->vertex_normals);}
                i_own_triangulated_surface=true;}}
        else if(triangulated_surface_filename.length() > 0){
            triangulated_surface=TRIANGULATED_SURFACE<T>::Create();
            triangulated_surface->Use_Vertex_Normals();triangulated_surface->vertex_normals=new ARRAY<TV>;
            int vertex_normals_length; // needed for backwards compatibility, since vertex_normals used to be ARRAYS<TV,1>
            FILE_UTILITIES::Read_From_File<RW>(triangulated_surface_filename,*triangulated_surface,vertex_normals_length,*triangulated_surface->vertex_normals);
            PHYSBAM_ASSERT(vertex_normals_length==1);
            i_own_triangulated_surface=true;}}
}
//#####################################################################
// Function Initialize_OpenGL_Triangulated_Surface
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Initialize_OpenGL_Triangulated_Surface()
{
    if(!opengl_triangulated_surface){
        Initialize_Triangulated_Surface();
        if(triangulated_surface){
            opengl_triangulated_surface=new OPENGL_TRIANGULATED_SURFACE<T>(*triangulated_surface,smooth_shading,front_surface_material,back_surface_material);
            if(triangulated_surface->vertex_normals){
                opengl_triangulated_surface->Initialize_Vertex_Normals();
                *opengl_triangulated_surface->vertex_normals=*triangulated_surface->vertex_normals;}}}
}
//#####################################################################
// Function Initialize_OpenGL_Scalar_Field
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Initialize_OpenGL_Scalar_Field()
{
    Initialize_Levelset();
    if(!levelset) return;

    if(!opengl_scalar_field && levelset){
        opengl_scalar_field=new OPENGL_SCALAR_FIELD_3D<T,T>(levelset->grid,levelset->phi,OPENGL_COLOR_RAMP<T>::Levelset_Color_Constant_Ramp(inside_slice_color,outside_slice_color));
        opengl_scalar_field->color_maps.Insert(OPENGL_COLOR_RAMP<T>::Levelset_Color_Linear_Ramp(inside_slice_color,outside_slice_color,opengl_scalar_field->grid.dX.x*10),2);}
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Initialize()
{
    if(!opengl_triangulated_surface) Update();
}
//#####################################################################
// Function Set_Surface_Material
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Set_Surface_Material(const OPENGL_MATERIAL &front_surface_mat,const OPENGL_MATERIAL &back_surface_mat)
{
    front_surface_material=front_surface_mat;
    back_surface_material=back_surface_mat;
    Update();
}
//#####################################################################
// Function Set_Overlayed_Surface_Material
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Set_Overlayed_Surface_Material(const OPENGL_MATERIAL &overlayed_surface_mat)
{
    overlayed_surface_material=overlayed_surface_mat;
    Update();
}
//#####################################################################
// Function Set_Slice_Color
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Set_Slice_Color(const OPENGL_COLOR &inside_slice_color_input,const OPENGL_COLOR &outside_slice_color_input)
{
    inside_slice_color=inside_slice_color_input;
    outside_slice_color=outside_slice_color_input;
    Update();
}
//#####################################################################
// Function Set_Slice_Color_Mode
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Set_Slice_Color_Mode(COLOR_MODE color_mode_input)
{
    color_mode=color_mode_input;
    Update();
}
//#####################################################################
// Function Set_Two_Sided
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Set_Two_Sided(bool two_sided_input)
{
    two_sided=two_sided_input;
}
//#####################################################################
// Function Toggle_Slice_Color_Mode
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Toggle_Slice_Color_Mode() 
{
    COLOR_MODE new_color_mode=(COLOR_MODE)(((int)color_mode+1)%2);
    Set_Slice_Color_Mode(new_color_mode);
}
//#####################################################################
// Function Toggle_Marching_Tetrahedra
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Toggle_Marching_Tetrahedra() 
{
    use_marching_tetrahedra=!use_marching_tetrahedra;
}
//#####################################################################
// Function Update
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Update()
{
    OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)this->slice;
    if(!slice || slice->mode==OPENGL_SLICE::NO_SLICE){
        Initialize_OpenGL_Triangulated_Surface();
        if(opengl_triangulated_surface){
            opengl_triangulated_surface->Set_Front_Material(front_surface_material);
            opengl_triangulated_surface->Set_Back_Material(back_surface_material);
            opengl_triangulated_surface->Set_Two_Sided(two_sided);}}
    else{
        Initialize_OpenGL_Scalar_Field();
        if(opengl_scalar_field){
            opengl_scalar_field->current_color_map=(color_mode==COLOR_SOLID)?1:2;
            opengl_scalar_field->Set_Slice(slice);
            opengl_scalar_field->Set_Smooth_Slice_Texture(smooth_slice_texture);}
        if(display_overlay){
            Initialize_OpenGL_Triangulated_Surface();
            if(opengl_triangulated_surface){
                opengl_triangulated_surface->Set_Front_Material(overlayed_surface_material);
                opengl_triangulated_surface->Set_Two_Sided(false);}}}
}
//#####################################################################
// Function Slice_Has_Changed
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Slice_Has_Changed()
{
    if(opengl_scalar_field) {
        opengl_scalar_field->Set_Slice(slice);
        opengl_scalar_field->Slice_Has_Changed();}
    Update();
}
//#####################################################################
// Function Toggle_Display_Overlay
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Toggle_Display_Overlay()
{
    display_overlay=!display_overlay;
    Update();
}
//#####################################################################
// Function Toggle_Smooth_Slice_Texture
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Toggle_Smooth_Slice_Texture()
{
    smooth_slice_texture=!smooth_slice_texture;
    if(opengl_scalar_field) opengl_scalar_field->Set_Smooth_Slice_Texture(smooth_slice_texture);
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Display(const int in_color) const
{
    OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)this->slice;
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();
    if(implicit_object_transform) {OpenGL_Translate(implicit_object_transform->t);OpenGL_Rotate(implicit_object_transform->r);}
    if(!slice || slice->mode==OPENGL_SLICE::NO_SLICE){if(opengl_triangulated_surface) opengl_triangulated_surface->Display(in_color);}
    else{
        if(opengl_scalar_field) opengl_scalar_field->Display(in_color);
        if(display_overlay && opengl_triangulated_surface) opengl_triangulated_surface->Display(in_color);}
    glPopMatrix();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Bounding_Box() const
{
    if(levelset) return (RANGE<VECTOR<float,3> >)levelset->grid.domain;
    else if(triangulated_surface){
        triangulated_surface->Update_Bounding_Box();return RANGE<VECTOR<float,3> >(*triangulated_surface->bounding_box);}
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Turn_Smooth_Shading_On
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Turn_Smooth_Shading_On()
{
    smooth_shading=true;
    if(opengl_triangulated_surface) opengl_triangulated_surface->Turn_Smooth_Shading_On();
}
//#####################################################################
// Function Turn_Smooth_Shading_Off
//#####################################################################
template<class T,class RW> void OPENGL_LEVELSET_MULTIVIEW<T,RW>::
Turn_Smooth_Shading_Off()
{
    smooth_shading=false;
    if(opengl_triangulated_surface) opengl_triangulated_surface->Turn_Smooth_Shading_Off();
}
//#####################################################################
template class OPENGL_LEVELSET_MULTIVIEW<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_LEVELSET_MULTIVIEW<double>;
#endif

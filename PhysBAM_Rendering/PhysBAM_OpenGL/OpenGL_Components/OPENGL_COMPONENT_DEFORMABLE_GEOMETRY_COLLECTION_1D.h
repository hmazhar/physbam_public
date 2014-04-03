//#####################################################################
// Copyright 2009, Nipun Kwatra
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D
//#####################################################################
#ifndef __OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D__
#define __OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D__

#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_HEXAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POINT_SIMPLICES_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POINTS_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T> class OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_1D;

template<class T,class RW=T>
class OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D:public OPENGL_COMPONENT
{
    typedef VECTOR<T,1> TV;
protected:
    std::string prefix;
    int frame_loaded;
    bool valid,use_active_list,hide_unselected;
    int display_mode,incremented_active_object;
    bool smooth_shading;
    int selected_vertex;
    bool invalidate_deformable_objects_selection_each_frame;
    bool own_deformable_geometry;
public:
    DEFORMABLE_GEOMETRY_COLLECTION<TV> *deformable_geometry_collection;
    OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_1D<T>* real_selection;
    ARRAY<OPENGL_POINT_SIMPLICES_1D<T>*> point_simplices_1d_objects;
    ARRAY<bool> active_list;
    ARRAY<TV> positions;
    ARRAY<TV> velocity_vectors;
    //OPENGL_VECTOR_FIELD_3D<T> velocity_field;
    bool draw_velocity_vectors;
    OPENGL_INDEXED_COLOR_MAP *color_map;

    OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D(const std::string& prefix,const int start_frame,const bool initialize_geometry=true);
    ~OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D();
    
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE;
    void Draw_All_Objects() PHYSBAM_OVERRIDE;

    virtual void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION *Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Set_Selection(OPENGL_SELECTION *selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION *selection) const PHYSBAM_OVERRIDE;
    
    void Turn_Smooth_Shading_On() PHYSBAM_OVERRIDE
    {smooth_shading=true;
    for(int i=1;i<=point_simplices_1d_objects.m;i++)if(point_simplices_1d_objects(i))point_simplices_1d_objects(i)->Turn_Smooth_Shading_On();}
    
    void Turn_Smooth_Shading_Off() PHYSBAM_OVERRIDE
    {smooth_shading=false;
    for(int i=1;i<=point_simplices_1d_objects.m;i++)if(point_simplices_1d_objects(i))point_simplices_1d_objects(i)->Turn_Smooth_Shading_Off();}

    OPENGL_SELECTION* Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection) PHYSBAM_OVERRIDE;


    void Toggle_Active_Value();
    void Toggle_Use_Active_List();
    void Toggle_Selection_Mode();
    void Increment_Active_Object();
    void Cycle_Display_Mode();
    void Show_Only_First();
    void Highlight_Particle();
    void Set_Vector_Size(double size);
    void Toggle_Velocity_Vectors();
    void Decrease_Vector_Size();
    void Increase_Vector_Size();
    void Update_Velocity_Field();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D,Toggle_Active_Value,"Toggle viewing of elements");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D,Toggle_Use_Active_List,"Toggle drawing subset of the deformable objects in the list");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D,Toggle_Selection_Mode,"Toggle selecting a whole segment or just one part");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D,Increment_Active_Object,"Increment deformable object being drawn");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D,Cycle_Display_Mode,"Cycle embedded display mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D,Show_Only_First,"Show only first deformable object");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D,Highlight_Particle,"Highlight a particle");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D,Toggle_Velocity_Vectors,"Toggle particle velocity vectors");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D,Increase_Vector_Size,"Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D,Decrease_Vector_Size,"Decrease vector size");

    virtual void Reinitialize(bool force=false,bool read_geometry=true);    // Needs to be called after some state changes
    
    void Toggle_Active_Value_Response();
    void Highlight_Particle_Response();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D,Toggle_Active_Value_Response,"");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D,Highlight_Particle_Response,"");
private:
    void Initialize();    // Needs to be called after some state changes
};

template<class T>
class OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_1D : public OPENGL_SELECTION
{
public:
    int body_index;
    int subobject_type;
    OPENGL_OBJECT* subobject;
    OPENGL_SELECTION* body_selection;
    OPENGL_SELECTION* saved_selection;

    OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_1D(OPENGL_OBJECT *object) : OPENGL_SELECTION(OPENGL_SELECTION::COMPONENT_DEFORMABLE_COLLECTION_3D,object) {}
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual OPENGL_SELECTION::TYPE Actual_Type() const {return body_selection->Actual_Type();}
};

}
#endif

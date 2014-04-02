//#####################################################################
// Copyright 2004-2009, Zhaosheng Bao, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D__
#define __OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D__

#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_HEXAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POINTS_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T> class OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D;

template<class T,class RW=T>
class OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D:public OPENGL_COMPONENT
{
    typedef VECTOR<T,3> TV;
protected:
    std::string prefix;
    int frame_loaded;
    bool valid,use_active_list,hide_unselected;
    int display_mode,display_relative_velocity_mode,number_of_segmented_curve,incremented_active_object;
    bool smooth_shading;
    int selected_vertex;
    bool invalidate_deformable_objects_selection_each_frame;
    bool own_deformable_geometry;
public:
    DEFORMABLE_GEOMETRY_COLLECTION<TV> *deformable_geometry;
    OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>* real_selection;
    ARRAY<OPENGL_SEGMENTED_CURVE_3D<T>*> segmented_curve_objects;
    ARRAY<OPENGL_TRIANGULATED_SURFACE<T>*> triangulated_surface_objects;
    ARRAY<OPENGL_TETRAHEDRALIZED_VOLUME<T>*> tetrahedralized_volume_objects;
    ARRAY<OPENGL_HEXAHEDRALIZED_VOLUME<T>*> hexahedralized_volume_objects;
    ARRAY<OPENGL_POINTS_3D<T,INDIRECT_ARRAY<ARRAY_VIEW<TV> > >*> free_particles_objects;
    ARRAY<INDIRECT_ARRAY<ARRAY_VIEW<TV> >*> free_particles_indirect_arrays;
    bool has_tetrahedralized_volumes,has_hexahedralized_volumes;
    ARRAY<bool> active_list;
    OPENGL_COLOR_RAMP<T>* color_map_relative_velocity;
    ARRAY<TV> positions;
    ARRAY<TV> velocity_vectors;
    OPENGL_VECTOR_FIELD_3D<T> velocity_field;
    bool draw_velocity_vectors;
    OPENGL_INDEXED_COLOR_MAP *color_map;

    OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D(const std::string& prefix,const int start_frame,const bool initialize_geometry=true);
    ~OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D();
    
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE;
    void Draw_All_Objects() PHYSBAM_OVERRIDE;

    virtual void Set_Display_Modes_For_Geometry_Collection(bool& display_triangulated_surface_objects,bool& display_tetrahedralized_volume_objects,
            bool& display_hexahedralized_volume_objects,bool& display_free_particles_objects) const;
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
    for(int i=1;i<=segmented_curve_objects.m;i++)if(segmented_curve_objects(i))segmented_curve_objects(i)->Turn_Smooth_Shading_On();
    for(int i=1;i<=triangulated_surface_objects.m;i++)if(triangulated_surface_objects(i))triangulated_surface_objects(i)->Turn_Smooth_Shading_On();
    for(int i=1;i<=tetrahedralized_volume_objects.m;i++)if(tetrahedralized_volume_objects(i))tetrahedralized_volume_objects(i)->Turn_Smooth_Shading_On();
    for(int i=1;i<=hexahedralized_volume_objects.m;i++)if(hexahedralized_volume_objects(i))hexahedralized_volume_objects(i)->Turn_Smooth_Shading_On();}
    
    void Turn_Smooth_Shading_Off() PHYSBAM_OVERRIDE
    {smooth_shading=false;
    for(int i=1;i<=segmented_curve_objects.m;i++)if(segmented_curve_objects(i))segmented_curve_objects(i)->Turn_Smooth_Shading_Off();
    for(int i=1;i<=triangulated_surface_objects.m;i++)if(triangulated_surface_objects(i))triangulated_surface_objects(i)->Turn_Smooth_Shading_Off();
    for(int i=1;i<=tetrahedralized_volume_objects.m;i++)if(tetrahedralized_volume_objects(i))tetrahedralized_volume_objects(i)->Turn_Smooth_Shading_Off();
    for(int i=1;i<=hexahedralized_volume_objects.m;i++)if(hexahedralized_volume_objects(i))hexahedralized_volume_objects(i)->Turn_Smooth_Shading_Off();}

    void Set_Material(const int object,const OPENGL_MATERIAL& front_material,const OPENGL_MATERIAL& back_material);
    void Set_All_Materials(const OPENGL_MATERIAL& meshfront,const OPENGL_MATERIAL& front_material,const OPENGL_MATERIAL& back_material);
    OPENGL_SELECTION* Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection) PHYSBAM_OVERRIDE;


    void Toggle_Active_Value();
    void Toggle_Draw_Interior();
    void Toggle_Draw_Subsets();
    void Toggle_Use_Active_List();
    void Toggle_Selection_Mode();
    void Toggle_Hide_Unselected();
    void Increment_Active_Object();
    void Cycle_Display_Mode();
    void Cycle_Cutaway_Mode();
    void Decrease_Cutaway_Fraction();
    void Increase_Cutaway_Fraction();
    void Create_One_Big_Triangulated_Surface_And_Write_To_File();
    void Show_Only_First();
    void Highlight_Particle();
    void Set_Vector_Size(double size);
    void Toggle_Velocity_Vectors();
    void Decrease_Vector_Size();
    void Increase_Vector_Size();
    void Update_Velocity_Field();
    void Cycle_Relative_Velocity_Mode();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Toggle_Active_Value,"Toggle viewing of elements");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Toggle_Draw_Interior,"Toggle view of interior elements for tetraheralized volumes");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Toggle_Draw_Subsets,"Toggle drawing of subset tets and particles for tetraheralized volumes");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Toggle_Hide_Unselected,"Toggle drawing of the selected regions");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Toggle_Use_Active_List,"Toggle drawing subset of the deformable objects in the list");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Toggle_Selection_Mode,"Toggle selecting a whole segment or just one part");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Increment_Active_Object,"Increment deformable object being drawn");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Cycle_Display_Mode,"Cycle embedded display mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Cycle_Cutaway_Mode,"Cycle cutaway mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Decrease_Cutaway_Fraction,"Decrease_Cutaway_Fraction");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Increase_Cutaway_Fraction,"Increase_Cutaway_Fraction");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Create_One_Big_Triangulated_Surface_And_Write_To_File,"Make one big boundary tri surface");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Cycle_Relative_Velocity_Mode,"Visualize relative velocity");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Show_Only_First,"Show only first deformable object");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Highlight_Particle,"Highlight a particle");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Toggle_Velocity_Vectors,"Toggle particle velocity vectors");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Increase_Vector_Size,"Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Decrease_Vector_Size,"Decrease vector size");

    virtual void Reinitialize(bool force=false,bool read_geometry=true);    // Needs to be called after some state changes
    
    void Toggle_Active_Value_Response();
    void Highlight_Particle_Response();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Toggle_Active_Value_Response,"");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D,Highlight_Particle_Response,"");
private:
    void Initialize();    // Needs to be called after some state changes
};

template<class T>
class OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D : public OPENGL_SELECTION
{
public:
    int body_index;
    int subobject_type;
    OPENGL_OBJECT* subobject;
    OPENGL_SELECTION* body_selection;
    OPENGL_SELECTION* saved_selection;

    OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D(OPENGL_OBJECT *object) : OPENGL_SELECTION(OPENGL_SELECTION::COMPONENT_DEFORMABLE_COLLECTION_3D,object) {}
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual OPENGL_SELECTION::TYPE Actual_Type() const {return body_selection->Actual_Type();}
};

}
#endif

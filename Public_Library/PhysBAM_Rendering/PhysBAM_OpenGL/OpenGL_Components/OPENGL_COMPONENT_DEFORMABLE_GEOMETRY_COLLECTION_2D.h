//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Frank Losasso, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D__
#define __OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D__

#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POINTS_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_VECTOR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T> class OPENGL_SEGMENTED_CURVE_2D;
template<class T> class OPENGL_TRIANGULATED_AREA;
template<class T_GRID> class OPENGL_ADAPTIVE_NODE_SCALAR_FIELD;

template<class T,class RW=T>
class OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D:public OPENGL_COMPONENT
{
    typedef VECTOR<T,2> TV;
protected:
    std::string prefix;
    int frame_loaded;
    bool valid;
    int display_mode;
    bool draw_velocities;
    T velocity_scale;
    bool invalidate_deformable_objects_selection_each_frame;
public:
    DEFORMABLE_GEOMETRY_COLLECTION<TV>* deformable_geometry_collection;
    ARRAY<OPENGL_SEGMENTED_CURVE_2D<T>*> segmented_curve_objects;
    ARRAY<OPENGL_TRIANGULATED_AREA<T>*> triangulated_area_objects;
    ARRAY<OPENGL_TRIANGULATED_AREA<T>*> triangles_of_material_objects;
    ARRAY<OPENGL_POINTS_2D<T,INDIRECT_ARRAY<ARRAY_VIEW<TV> > >*> free_particles_objects;
    ARRAY<INDIRECT_ARRAY<ARRAY_VIEW<TV> >*> free_particles_indirect_arrays;
    OPENGL_VECTOR_FIELD_2D<ARRAY_VIEW<TV> > velocity_field;
    OPENGL_INDEXED_COLOR_MAP *color_map;

    OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D(const std::string& prefix,const int start_frame);
    ~OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D();
    
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE;

    virtual void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION* Get_Selection(GLuint* buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION* selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;
    OPENGL_SELECTION* Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection) PHYSBAM_OVERRIDE;

    void Set_Vector_Size(const T vector_size);

    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Draw_Velocities();
    void Cycle_Display_Mode();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D,Increase_Vector_Size,"Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D,Decrease_Vector_Size,"Decrease vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D,Toggle_Draw_Velocities,"Toggle draw velocities");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D,Cycle_Display_Mode,"Cycle embedded display mode");

    virtual void Reinitialize(bool force=false);    // Needs to be called after some state changes

private:
    void Initialize();    // Needs to be called after some state changes
};

template<class T>
class OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_2D:public OPENGL_SELECTION
{
public:
    int body_index;
    int subobject_type;
    OPENGL_OBJECT* subobject;
    OPENGL_SELECTION* body_selection;

    OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_2D(OPENGL_OBJECT *object):OPENGL_SELECTION(OPENGL_SELECTION::COMPONENT_DEFORMABLE_OBJECT_2D, object){}
    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    OPENGL_SELECTION::TYPE Actual_Type() const PHYSBAM_OVERRIDE {return body_selection->Actual_Type();}
};

}
#endif

//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_SLICE_MANAGER
//##################################################################### 
#ifndef __OPENGL_SLICE_MANAGER__
#define __OPENGL_SLICE_MANAGER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SLICE.h>
namespace PhysBAM{

class OPENGL_SLICE_MANAGER
{
public:
    OPENGL_SLICE* slice;
    ARRAY<OPENGL_OBJECT*> objects;
private:
    OPENGL_CALLBACK* slice_has_changed_callback;
   
public:
    OPENGL_SLICE_MANAGER()
        :slice(0),slice_has_changed_callback(0)
    {}

    void Set_Slice_Has_Changed_Callback(OPENGL_CALLBACK* slice_has_changed_callback_input)
    {
        slice_has_changed_callback = slice_has_changed_callback_input;
    }
    
    void Add_Object(OPENGL_OBJECT* object)
    {
        objects.Append(object);
        object->Set_Slice(slice);
    }

    void Update_Objects()
    {
        for(int i=1;i<=objects.m;i++) objects(i)->Slice_Has_Changed();
    }

    // Callbacks for convenience
    void Toggle_Slice_Mode() { slice->Toggle_Slice_Mode(); Update_Objects(); if(slice_has_changed_callback) (*slice_has_changed_callback)();}
    void Toggle_Slice_Axis() { slice->Toggle_Slice_Axis(); Update_Objects(); if(slice_has_changed_callback) (*slice_has_changed_callback)();}
    void Increment_Slice() { slice->Increment_Slice(); Update_Objects(); if(slice_has_changed_callback) (*slice_has_changed_callback)();}
    void Decrement_Slice() { slice->Decrement_Slice(); Update_Objects(); if(slice_has_changed_callback) (*slice_has_changed_callback)();}

    DEFINE_CALLBACK_CREATOR(OPENGL_SLICE_MANAGER, Toggle_Slice_Mode);
    DEFINE_CALLBACK_CREATOR(OPENGL_SLICE_MANAGER, Toggle_Slice_Axis);
    DEFINE_CALLBACK_CREATOR(OPENGL_SLICE_MANAGER, Increment_Slice);
    DEFINE_CALLBACK_CREATOR(OPENGL_SLICE_MANAGER, Decrement_Slice);
};
}
#endif

//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __OPENGL_COMPONENT_DIAGNOSTICS__
#define __OPENGL_COMPONENT_DIAGNOSTICS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

class OPENGL_COMPONENT_DIAGNOSTICS:public OPENGL_COMPONENT
{
    std::string filename;
    int frame_loaded;
    bool valid;
    ARRAY<std::string> lines;

//#####################################################################
public:
    OPENGL_COMPONENT_DIAGNOSTICS(const std::string& filename);
    virtual ~OPENGL_COMPONENT_DIAGNOSTICS();
private:
    void Reinitialize();
    void Print_Selection_Info(std::ostream& ostream,OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;    
    bool Valid_Frame(int frame) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame) PHYSBAM_OVERRIDE;
    //virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const;
//#####################################################################
};
}
#endif

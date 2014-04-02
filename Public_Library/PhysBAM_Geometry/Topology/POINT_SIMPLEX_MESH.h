//#####################################################################
// Copyright 2007-2009, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __POINT_SIMPLEX_MESH__
#define __POINT_SIMPLEX_MESH__

#include <PhysBAM_Geometry/Topology/SIMPLEX_MESH.h>
namespace PhysBAM{

class POINT_SIMPLEX_MESH:public SIMPLEX_MESH<0>
{
public:
    typedef VECTOR<int,1> ELEMENT_TYPE;
    ARRAY<bool> directions; // false for left and true for right
private:
    SEGMENT_MESH* segment_mesh; // segment mesh consisting of edges joining consecutive elements

public:
    POINT_SIMPLEX_MESH();
    ~POINT_SIMPLEX_MESH();

    SEGMENT_MESH& Get_Segment_Mesh()
    {if(!segment_mesh) Initialize_Segment_Mesh();return *segment_mesh;}

//#####################################################################
private:
    void Initialize_Segment_Mesh();
//#####################################################################
};
}
#endif

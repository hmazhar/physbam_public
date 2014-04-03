//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Andrew Selle, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HEXAHEDRALIZED_VOLUME
//#####################################################################
#ifndef __HEXAHEDRALIZED_VOLUME__
#define __HEXAHEDRALIZED_VOLUME__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Topology/HEXAHEDRON_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>
namespace PhysBAM{

template<class TV> class GRID;
template<class T> class HEXAHEDRON;
template<class T> class TETRAHEDRALIZED_VOLUME;
template<class T> class TRIANGULATED_SURFACE;

template<class T_input>
class HEXAHEDRALIZED_VOLUME:public MESH_OBJECT<VECTOR<T_input,3>,HEXAHEDRON_MESH>
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef MESH_OBJECT<TV,HEXAHEDRON_MESH> BASE;
    using BASE::mesh;using BASE::particles;

    ARRAY<HEXAHEDRON<T> >* hexahedron_list;
    TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume;
    TRIANGULATED_SURFACE<T>* triangulated_surface;
    BOX_HIERARCHY<TV>* hierarchy;

    HEXAHEDRALIZED_VOLUME(HEXAHEDRON_MESH& mesh_input,GEOMETRY_PARTICLES<TV>& particles_input);
    virtual ~HEXAHEDRALIZED_VOLUME();

    virtual std::string Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    static std::string Static_Name()
    {return "HEXAHEDRALIZED_VOLUME<T>";}

    void Initialize_Hierarchy()
    {PHYSBAM_NOT_IMPLEMENTED();}

//#####################################################################
    void Clean_Memory() PHYSBAM_OVERRIDE;
    void Update_Hexahedron_List(); // updates the hexahedrons assuming the particle positions are already updated
    void Initialize_Tetrahedralized_Volume();
    void Initialize_Triangulated_Surface();
    void Initialize_Cube_Mesh_And_Particles(const GRID<TV>& grid);
    T Total_Volume() const;
private:
    void Refresh_Auxiliary_Structures_Helper() PHYSBAM_OVERRIDE {PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
};
}
#endif

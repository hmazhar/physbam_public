//#####################################################################
// Copyright 2002-2004, Robert Bridson, Ronald Fedkiw, Eilene Hao, Geoffrey Irving, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_SUBDIVISION
//#####################################################################
#ifndef __TRIANGLE_SUBDIVISION__
#define __TRIANGLE_SUBDIVISION__

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

class TRIANGLE_MESH;

class TRIANGLE_SUBDIVISION:public NONCOPYABLE
{
public:
    TRIANGLE_MESH& triangle_mesh;
    int start_index_for_new_nodes;
    bool delete_segment_mesh; // delete if defined in Refine_Mesh() function
    bool delete_neighbor_nodes; // delete if defined in Apply_Linear_Subdivision() or Apply_Root_Three_Subdivision() functions
    bool delete_topologically_sorted_neighbor_nodes; // delete if defined in Apply_Loop_Subdivision()
    bool delete_boundary_mesh; // delete if defined in Apply_Loop_Subdivision() function
    bool delete_incident_elements; // delete if defined in Apply_Fractal_Subdivision() function

    TRIANGLE_SUBDIVISION(TRIANGLE_MESH& triangle_mesh_input)
        :triangle_mesh(triangle_mesh_input),start_index_for_new_nodes(0),delete_segment_mesh(false),delete_neighbor_nodes(false),delete_topologically_sorted_neighbor_nodes(false),
        delete_boundary_mesh(false),delete_incident_elements(false)
    {}

    virtual ~TRIANGLE_SUBDIVISION()
    {Clean_Memory();}

    virtual bool Node_Is_A_Corner(const int node_index) // customize this if you know the mesh has corners
    {return false;}

    template<class T_ARRAY1,class T_ARRAY2> void Apply_Linear_Subdivision(const T_ARRAY1& base_values,T_ARRAY2& subdivided_values)
    {Apply_Linear_Subdivision<typename T_ARRAY1::ELEMENT>(base_values,subdivided_values);}

    template<class T_ARRAY1,class T_ARRAY2> void Apply_Fractal_Subdivision(const T_ARRAY1& base_values,T_ARRAY2& subdivided_values,const float power)
    {Apply_Fractal_Subdivision<typename T_ARRAY1::ELEMENT>(base_values,subdivided_values,power);}

    template<class T_ARRAY1,class T_ARRAY2> void Apply_Loop_Subdivision(const T_ARRAY1& base_values,T_ARRAY2& subdivided_values)
    {Apply_Loop_Subdivision<typename T_ARRAY1::ELEMENT>(base_values,subdivided_values);}

    template<class T_ARRAY1,class T_ARRAY2> void Apply_Root_Three_Subdivision(const T_ARRAY1& base_values,T_ARRAY2& subdivided_values)
    {Apply_Root_Three_Subdivision<typename T_ARRAY1::ELEMENT>(base_values,subdivided_values);}

//#####################################################################
    virtual void Clean_Memory();
    void Refine_Mesh(TRIANGLE_MESH& refined_triangle_mesh,const int start_index_for_new_nodes_input=0);
    void Refine_Mesh_Dual(TRIANGLE_MESH& refined_triangle_mesh,const int start_index_for_new_nodes_input=0);
    template<class TV> void Apply_Linear_Subdivision(ARRAY_VIEW<const TV> base_values,ARRAY_VIEW<TV> subdivided_values);
    template<class TV> void Apply_Fractal_Subdivision(ARRAY_VIEW<const TV> base_values,ARRAY_VIEW<TV> subdivided_values,const float power);
    template<class TV> void Apply_Loop_Subdivision(ARRAY_VIEW<const TV> base_values,ARRAY_VIEW<TV> subdivided_values);
    template<class TV> void Apply_Root_Three_Subdivision(ARRAY_VIEW<const TV> base_values,ARRAY_VIEW<TV> subdivided_values);
//#####################################################################
};
}
#endif

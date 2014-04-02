//#####################################################################
// Copyright 2007, Avi Robinson-Mosher, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_OBJECT_INTERSECTOR
//#####################################################################
#ifndef __IMPLICIT_OBJECT_INTERSECTOR__
#define __IMPLICIT_OBJECT_INTERSECTOR__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class TV>
class IMPLICIT_OBJECT_INTERSECTOR
{
private:
    typedef typename TV::SCALAR T;typedef typename GRID<TV>::VECTOR_INT TV_INT;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::HYPERPLANE T_HYPERPLANE;typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID<TV>::NODE_ITERATOR NODE_ITERATOR;typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
    typedef VECTOR<int,TV::dimension+1> T_ELEMENT;typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension>::SIMPLEX T_SIMPLEX;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef IMPLICIT_OBJECT<TV> BASE;
    
public:
    const IMPLICIT_OBJECT<TV>& implicit_object;
    int maximum_refinement_depth;
    int minimum_refinement_depth;
    T full_cell_size;

    ARRAY<T_ELEMENT> cell_refinement_simplices;
    ARRAY<T> cell_phis;
    ARRAY<TV> cell_particle_X;
    T_ARRAYS_SCALAR grid_nodal_phis;

    IMPLICIT_OBJECT_INTERSECTOR(const IMPLICIT_OBJECT<TV>& implicit_object_input);
    virtual ~IMPLICIT_OBJECT_INTERSECTOR();

    void Set_Full_Cell_Size(const T full_cell_size_input)
    {full_cell_size=full_cell_size_input;}

    void Set_Maximum_Refinement_Depth(const int maximum_refinement_depth_input)
    {maximum_refinement_depth=maximum_refinement_depth_input;}

    void Set_Minimum_Refinement_Depth(const int minimum_refinement_depth_input)
    {minimum_refinement_depth=minimum_refinement_depth_input;}

    static bool Refinement_Condition(const ARRAY<TV>& X,const ARRAY<T>& phis,const T_ELEMENT& indices)
    {T maximum_edge_length_magnitude_squared=0; 
    for(int i=1;i<=GRID<TV>::dimension;i++)for(int j=i+1;j<=GRID<TV>::dimension+1;j++)
        maximum_edge_length_magnitude_squared=max(maximum_edge_length_magnitude_squared,(X(indices[i])-X(indices[j])).Magnitude_Squared());
    T result=sqr(phis(indices[1]));for(int i=2;i<=GRID<TV>::dimension+1;i++) result=min(result,sqr(phis(indices[i])));
    return result<=maximum_edge_length_magnitude_squared;}

    static void Refined_Object_Initialization_Helper(ARRAY<VECTOR<int,4> >& tets)
    {// corner tets
    tets.Append(VECTOR<int,4>(2,6,4,1));
    tets.Append(VECTOR<int,4>(3,7,1,4));
    tets.Append(VECTOR<int,4>(5,6,1,7));
    tets.Append(VECTOR<int,4>(8,6,7,4));
    // inside tet
    tets.Append(VECTOR<int,4>(1,7,6,4));}

    static void Refined_Object_Initialization_Helper(ARRAY<VECTOR<int,3> >& tris)
    {tris.Append(VECTOR<int,3>(1,2,3));tris.Append(VECTOR<int,3>(2,4,3));}

    static void Refined_Object_Initialization_Helper(ARRAY<VECTOR<int,2> >& segments)
    {segments.Append(VECTOR<int,2>(1,2));}

    static void Make_List(ARRAY<VECTOR<T,3> >& points,const ARRAY<VECTOR<T,3> ,VECTOR<int,3> >& corners)
    {points.Resize(8);int index=1;for(int k=0;k<=1;k++)for(int j=0;j<=1;j++)for(int i=0;i<=1;i++) points(index++)=corners(i,j,k);}

    static void Make_List(ARRAY<VECTOR<T,2> >& points,const ARRAY<VECTOR<T,2> ,VECTOR<int,2> >& corners)
    {points.Resize(4);int index=1;for(int j=0;j<=1;j++)for(int i=0;i<=1;i++) points(index++)=corners(i,j);}

    static void Make_List(ARRAY<VECTOR<T,1> >& points,const ARRAY<VECTOR<T,1> ,VECTOR<int,1> >& corners)
    {points.Resize(2);int index=1;for(int i=0;i<=1;i++) points(index++)=corners(i);}

    static void Refine_Simplex(const IMPLICIT_OBJECT_INTERSECTOR<VECTOR<T,1> >& implicit_object_intersector,ARRAY<VECTOR<int,2> >& segments,
        ARRAY<VECTOR<T,1> >& particle_X,const VECTOR<int,2>& indices)
    {const int x1=indices[1],x2=indices[2];
    const VECTOR<T,1> X1=particle_X(x1),X2=particle_X(x2);
    int midpoint=particle_X.Append(implicit_object_intersector.Iterative_Find_Interface(X1,X2));
    segments.Append(VECTOR<int,2>(x1,midpoint));
    segments.Append(VECTOR<int,2>(midpoint,x2));}

    static void Refine_Simplex(const IMPLICIT_OBJECT_INTERSECTOR<VECTOR<T,2> >& implicit_object_intersector,ARRAY<VECTOR<int,3> >& tris,
        ARRAY<VECTOR<T,2> >& particle_X,const VECTOR<int,3>& indices)
    {const int x1=indices[1],x2=indices[2],x3=indices[3];
    const VECTOR<T,2> X1=particle_X(x1),X2=particle_X(x2),X3=particle_X(x3);
    int first_midpoint=particle_X.Append(implicit_object_intersector.Iterative_Find_Interface(X1,X2));
    int second_midpoint=particle_X.Append(implicit_object_intersector.Iterative_Find_Interface(X2,X3));
    int third_midpoint=particle_X.Append(implicit_object_intersector.Iterative_Find_Interface(X3,X1));
    tris.Append(VECTOR<int,3>(x1,first_midpoint,third_midpoint));
    tris.Append(VECTOR<int,3>(x2,second_midpoint,first_midpoint));
    tris.Append(VECTOR<int,3>(x3,third_midpoint,second_midpoint));
    tris.Append(VECTOR<int,3>(first_midpoint,second_midpoint,third_midpoint));}

    static void Refine_Simplex(const IMPLICIT_OBJECT_INTERSECTOR<VECTOR<T,3> >& implicit_object_intersector,ARRAY<VECTOR<int,4> >& tets,
        ARRAY<VECTOR<T,3> >& particle_X,const VECTOR<int,4>& indices)
    {const int x1=indices[1],x2=indices[2],x3=indices[3],x4=indices[4];
    const VECTOR<T,3> X1=particle_X(x1),X2=particle_X(x2),X3=particle_X(x3),X4=particle_X(x4);
    int first_midpoint=particle_X.Append(implicit_object_intersector.Iterative_Find_Interface(X1,X2));
    int second_midpoint=particle_X.Append(implicit_object_intersector.Iterative_Find_Interface(X1,X3));
    int third_midpoint=particle_X.Append(implicit_object_intersector.Iterative_Find_Interface(X1,X4));
    int fourth_midpoint=particle_X.Append(implicit_object_intersector.Iterative_Find_Interface(X2,X3));
    int fifth_midpoint=particle_X.Append(implicit_object_intersector.Iterative_Find_Interface(X2,X4));
    int sixth_midpoint=particle_X.Append(implicit_object_intersector.Iterative_Find_Interface(X3,X4));
    tets.Append(VECTOR<int,4>(x1,first_midpoint,second_midpoint,third_midpoint));
    tets.Append(VECTOR<int,4>(x2,fourth_midpoint,first_midpoint,fifth_midpoint));
    tets.Append(VECTOR<int,4>(x3,second_midpoint,fourth_midpoint,sixth_midpoint));
    tets.Append(VECTOR<int,4>(x4,third_midpoint,sixth_midpoint,fifth_midpoint));
    // pick our octagon diagonal as first_midpoint,sixth_midpoint
    tets.Append(VECTOR<int,4>(first_midpoint,sixth_midpoint,third_midpoint,fifth_midpoint));
    tets.Append(VECTOR<int,4>(first_midpoint,sixth_midpoint,fifth_midpoint,fourth_midpoint));
    tets.Append(VECTOR<int,4>(first_midpoint,sixth_midpoint,fourth_midpoint,second_midpoint));
    tets.Append(VECTOR<int,4>(first_midpoint,sixth_midpoint,second_midpoint,third_midpoint));}

    void Precompute_Nodal_Phis(const GRID<TV>& grid)
    {grid_nodal_phis.Resize(grid.Node_Indices(3));for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) grid_nodal_phis(iterator.Node_Index())=Phi(iterator.Location());}

    T Extended_Phi(const TV& location) const
    {return implicit_object.Extended_Phi(location);}
    
//#####################################################################
    TV Iterative_Find_Interface(TV left,TV right,const int interations=7) const;
    T Negative_Material_In_Cell(const GRID<TV>& grid,const TV_INT& cell_index,const bool force_full_refinement=false);
    T Negative_Material_In_Box(const RANGE<TV>& box,const bool force_full_refinement=false);
    T Negative_Material_In_Box_Excluding_Object(const RANGE<TV>& box,const ARRAY<IMPLICIT_OBJECT<TV>*>& excluded_implicit_object,const bool force_full_refinement=false);
//#####################################################################
};
}
#endif

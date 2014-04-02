//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BINTREE_GRID
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __BINTREE_GRID__
#define __BINTREE_GRID__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

template<class T> class MAP_BINTREE_MESH;
template<class T> class BINTREE_CELL_HELPER;
template<class T_GRID> class DYADIC_GRID_ITERATOR_NODE;
template<class T_GRID> class DYADIC_GRID_ITERATOR_CELL;
template<class T_GRID> class DYADIC_GRID_ITERATOR_FACE;
template<class TV> struct DYADIC_TAG;
template<class T_GRID> class BLOCK_DYADIC;

template<class T_input>
class BINTREE_GRID:public NONCOPYABLE
{
    typedef T_input T;
    typedef VECTOR<T,1> TV;typedef VECTOR<int,1> TV_INT;
  public:
    static const int number_of_children_per_cell=2;
    static const T one_over_number_of_children_per_cell;
    static const int number_of_cells_per_block=2;
    static const int number_of_faces_per_block=3;
    static const int number_of_faces_per_cell=2;
    static const int number_of_nodes_per_cell=2;
    enum {number_of_cells_per_node=2}; // workaround for ICC bug
    static const T one_over_number_of_nodes_per_cell;
    static const int number_of_nodes_per_face=2;
    static const T one_over_number_of_nodes_per_face;
    static const int number_of_neighbors_per_node=2;
    static const int number_of_neighbors_per_cell=2;
    static const int number_of_neighbors_per_face=2;
    static const int number_of_complete_one_ring_neighbors_face_per_cell=2;
    static const int number_of_complete_one_ring_neighbors_edge_per_cell=0;
    static const int number_of_complete_one_ring_neighbors_node_per_cell=2;
    static const int number_of_complete_one_ring_neighbors_per_cell=2;
    static const int dimension=1;
    static const char* name;

    typedef T SCALAR;
    typedef VECTOR<T,1> VECTOR_T;
    typedef DYADIC_TAG<VECTOR_T> GRID_TAG;

    typedef GRID<TV> UNIFORM_GRID;
    typedef BINTREE_GRID<T> GRID_T;
    typedef VECTOR<int,1> VECTOR_INT;
    typedef int INDEX; typedef DYADIC_FACE_INDEX<TV::m> INDEX_FACE;
    typedef MAP_BINTREE_MESH<T> MAP_MESH;
    typedef BINTREE_CELL<T> CELL;
    typedef BINTREE_CELL_HELPER<T> CELL_HELPER;
    typedef BLOCK_DYADIC<BINTREE_GRID<T> > BLOCK;

    typedef DYADIC_GRID_ITERATOR_FACE<GRID_T> FACE_ITERATOR;
    typedef DYADIC_GRID_ITERATOR_NODE<GRID_T> NODE_ITERATOR;
    typedef DYADIC_GRID_ITERATOR_CELL<GRID_T> CELL_ITERATOR;

    ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> > cells;
    GRID<TV> uniform_grid;
    T minimum_cell_size;
    T one_over_minimum_cell_size;
    T minimum_cell_size_over_two;
    T one_over_minimum_cell_size_over_two;
    int number_of_ghost_cells;
    int number_of_cells;
    int number_of_nodes;
    int number_of_faces;
    int maximum_depth;

    // node iterator arrays
    mutable ARRAY<int> internal_nodes;
    mutable ARRAY<int> boundary_nodes;
    mutable ARRAY<ARRAY<int> > individual_side_boundary_nodes; // left,right,bottom,top,front,back ghost nodes
    mutable ARRAY<int> ghost_nodes;
    mutable ARRAY<ARRAY<int> > individual_side_ghost_nodes; // left,right,bottom,top,front,back ghost nodes
    mutable ARRAY<int> node_iterator_deepest_cells;
    mutable ARRAY<unsigned char> nodes;
    mutable bool node_iterator_data_up_to_date;

    // face iterator arrays
    mutable ARRAY<int> internal_faces; // the faces that are inside (not including the boundary) of the domain
    mutable ARRAY<int> boundary_faces; // the faces that lie on the domain
    mutable ARRAY<ARRAY<int> > individual_side_boundary_faces; // the faces that lie on the individual sides of the domain
    mutable ARRAY<ARRAY<int> > individual_side_domain_ghost_faces; // the faces that are are on xmin, xmax, etc, minus the boundary faces...
    mutable ARRAY<int> ghost_faces; // the faces that are in the ghost region (not including the boundary
    mutable ARRAY<ARRAY<int> > individual_side_ghost_faces; // the ghost faces on the individual sides (not including the boundary)
    mutable ARRAY<int> face_iterator_deepest_cells;
    mutable ARRAY<unsigned char> faces;
    mutable bool face_iterator_data_up_to_date;
private:
    ARRAY<BINTREE_CHILDREN<T>*> allocated_children;
    mutable bool cell_pointer_from_index_up_to_date;mutable ARRAY<BINTREE_CELL<T>*> cell_pointer_from_index;
    mutable bool neighbors_up_to_date;mutable ARRAY<VECTOR<BINTREE_CELL<T>*,2> > neighbors;
    mutable bool node_neighbors_up_to_date;mutable ARRAY<VECTOR<int,2> > node_neighbors;
public:
    mutable bool fully_refined_block_up_to_date;mutable ARRAY<bool> fully_refined_block;

    BINTREE_GRID();
    ~BINTREE_GRID();

    T Minimum_Edge_Length() const
    {return minimum_cell_size;}

    TV Minimum_Cell_DX() const
    {return TV(minimum_cell_size);}

    TV One_Over_Minimum_Cell_DX() const
    {return TV(one_over_minimum_cell_size);}

    TV Minimum_Cell_DX_Over_Two() const
    {return TV(minimum_cell_size_over_two);}

    void Compute_Minimum_Cell_Constants(int proposed_maximum_depth)
    {minimum_cell_size=uniform_grid.dX.x/(1<<(proposed_maximum_depth-1));
    one_over_minimum_cell_size=(T)1./minimum_cell_size;
    minimum_cell_size_over_two=(T).5*minimum_cell_size;
    one_over_minimum_cell_size_over_two=(T)1./minimum_cell_size_over_two;}

    void Set_Maximum_Depth(const int maximum_depth_input)
    {maximum_depth=maximum_depth_input;Compute_Minimum_Cell_Constants(maximum_depth);}

    void Update_Maximum_Depth()
    {maximum_depth=0;for(int i=1;i<=uniform_grid.numbers_of_cells.x;i++) maximum_depth=max(maximum_depth,cells(i)->Maximum_Relative_Depth());}

    TV Clamp(const TV& location) const
    {return uniform_grid.Clamp(location);}

    const RANGE<TV>& Domain() const
    {return uniform_grid.domain;}

    bool Outside(const TV& location) const
    {return uniform_grid.Outside(location);}

    // TODO(jontg): Useless functions?
    // static TV Face_Corner_To_Opposite_Corner_Normalized_Vector(const int axis,const int corner)
    // {PHYSBAM_FATAL_ERROR();}

    // BINTREE_CELL<T>* Complete_Minimal_One_Ring_Neighbor_Face(const int i,BINTREE_CELL<T>* cell) const
    // {assert(1<=i&&i<=2);BINTREE_CELL<T> *temp=neighbors(cell->Cell())(i);if(temp&&temp->Depth_Of_This_Cell()==maximum_depth)return temp;return 0;}
    // 
    // BINTREE_CELL<T>* Complete_Minimal_One_Ring_Neighbor_Edge(const int i,BINTREE_CELL<T>* cell) const
    // {PHYSBAM_FATAL_ERROR();}

    // BINTREE_CELL<T>* Complete_Minimal_One_Ring_Neighbor_Node(const int i,BINTREE_CELL<T>* cell) const
    // {assert(1<=i&&i<=4);BINTREE_CELL<T> *temp,*temp2;int lookup[4][2]={{1,3},{2,3},{1,4},{2,4}};
    // temp=neighbors(cell->Cell())(lookup[i-1][0]);
    // if(temp&&temp->Depth_Of_This_Cell()==maximum_depth){temp2=neighbors(temp->Cell())(lookup[i-1][1]);if(temp2&&temp2->Depth_Of_This_Cell()==maximum_depth)return temp2;}
    // temp=neighbors(cell->Cell())(lookup[i-1][1]);
    // if(temp&&temp->Depth_Of_This_Cell()==maximum_depth){temp2=neighbors(temp->Cell())(lookup[i-1][0]);if(temp2&&temp2->Depth_Of_This_Cell()==maximum_depth)return temp2;}
    // return 0;}

    TV Node_Location(const int global_node_index) const
    {return Node_Location(nodes(global_node_index),cell_pointer_from_index(node_iterator_deepest_cells(global_node_index)));}

    TV Node_Location(const int node_index,const BINTREE_CELL<T>* cell) const
    {assert(node_index >=0 && node_index <= 1);
    TV center((TV)cell->Center()),DX_over_two((T).5*(TV)cell->DX());
    switch(node_index){
        case 0:return TV(center.x-DX_over_two.x);
        default:return TV(center.x+DX_over_two.x);}}



    TV Face_Location(const int global_face_index) const
    {return Face_Location(faces(global_face_index),cell_pointer_from_index(face_iterator_deepest_cells(global_face_index)));}

    char Face_Axis(const int global_face_index) const
    {return 0;}

    TV Face_Location(const int face_index,const BINTREE_CELL<T>* cell) const
    {assert(face_index >=0 && face_index <= 1);
    TV center((TV)cell->Center()),DX_over_two((T).5*(TV)cell->DX());
    switch(face_index){
        case 0:return TV(center.x-DX_over_two.x);
        default:return TV(center.x+DX_over_two.x);}}



    TV Cell_Location(const int global_cell_index) const
    {return Cell_Pointer_From_Index()(global_cell_index)->Center();}

    ARRAY<VECTOR<BINTREE_CELL<T>*,2> >& Neighbors() const
    {if(!neighbors_up_to_date){Calculate_Neighbors_Array(neighbors);neighbors_up_to_date=true;}
    return neighbors;}

    ARRAY<VECTOR<int,2> >& Node_Neighbors() const
    {if(!node_neighbors_up_to_date){Calculate_Node_Neighbors_Array(node_neighbors);node_neighbors_up_to_date=true;}
    return node_neighbors;}

    ARRAY<BINTREE_CELL<T>*>& Cell_Pointer_From_Index() const
    {if(!cell_pointer_from_index_up_to_date){Calculate_Cell_Pointer_From_Index_Array(cell_pointer_from_index);cell_pointer_from_index_up_to_date=true;}
    return cell_pointer_from_index;}

    ARRAY<bool>& Fully_Refined_Block() const
    {if(!fully_refined_block_up_to_date){Calculate_Fully_Refined_Block_Array(fully_refined_block);fully_refined_block_up_to_date=true;}
    return fully_refined_block;}

    const BINTREE_CELL<T>* Base_Cell(const TV& X) const
    {BINTREE_CELL<T>* base_cell=Leaf_Cell(X-Minimum_Cell_DX_Over_Two());
    if(base_cell&&fully_refined_block(base_cell->Cell())){assert(abs(minimum_cell_size-base_cell->DX().x)<1e-6);return base_cell;}else return 0;}

    const BINTREE_CELL<T>* Base_Cell_By_Neighbor_Path(const BINTREE_CELL<T>* start_cell,const TV& location,const T tolerance=1e-3,const int max_iterations=10) const
    {const BINTREE_CELL<T>* base_cell=Leaf_Cell_By_Neighbor_Path(start_cell,location-Minimum_Cell_DX_Over_Two(),tolerance,max_iterations);
    if(base_cell&&fully_refined_block(base_cell->Cell()))return base_cell;return 0;}

    BINTREE_CELL<T>* Base_Cell_By_Neighbor_Path(BINTREE_CELL<T>* start_cell,const TV& location,const T tolerance=1e-3,const int max_iterations=10) const
    {BINTREE_CELL<T>* base_cell=Leaf_Cell_By_Neighbor_Path(start_cell,location-Minimum_Cell_DX_Over_Two(),tolerance,max_iterations);
    if(base_cell&&fully_refined_block(base_cell->Cell()))return base_cell;return 0;}

    // TODO(jontg): Don't really understand cell blocks...
    // const BINTREE_CELL<T>* Cell_From_Cell_Block(const BINTREE_CELL<T>* base_cell,const int number) const // 0 to 1 in 1D
    // {assert(fully_refined_block(base_cell->Cell()));const BINTREE_CELL<T>* cell=base_cell;
    // if(number&1) cell=neighbors(cell->Cell())(2); else cell=neighbors(cell->Cell())(4);
    // return cell;}

    void All_Cells_In_Cell_Block(const BINTREE_CELL<T>* base_cell,const BINTREE_CELL<T>* cells[number_of_cells_per_block]) const
    {assert(fully_refined_block(base_cell->Cell()));cells[0]=base_cell;int lookup[][2]={{0,0},{2,0}};
    for(int i=1;i<number_of_cells_per_block;i++) cells[i]=neighbors(cells[lookup[i][1]]->Cell())(lookup[i][0]);}

    int First_Face_Index_In_Cell(const int axis,const int cell_index) const
    {return cell_pointer_from_index(cell_index)->Face(0);}

    int Second_Face_Index_In_Cell(const int axis,const int cell_index) const
    {return cell_pointer_from_index(cell_index)->Face(1);}

    ARRAY<ARRAY<int>*> Map_All_Nodes() const
    {ARRAY<ARRAY<int>*> indirection_arrays(3);Node_Iterator_Data();
    indirection_arrays(1)=&internal_nodes;
    indirection_arrays(2)=&boundary_nodes;
    indirection_arrays(3)=&ghost_nodes;
    return indirection_arrays;}

    ARRAY<ARRAY<int>*> Map_Interior_Nodes() const
    {ARRAY<ARRAY<int>*> indirection_arrays(1);Node_Iterator_Data();
    indirection_arrays(1)=&internal_nodes;
    return indirection_arrays;}

    ARRAY<ARRAY<int>*> Map_Boundary_Nodes() const
    {ARRAY<ARRAY<int>*> indirection_arrays(1);Node_Iterator_Data();
    indirection_arrays(1)=&boundary_nodes;
    return indirection_arrays;}

    ARRAY<ARRAY<int>*> Map_Regular_Nodes() const
    {ARRAY<ARRAY<int>*> indirection_arrays(2);Node_Iterator_Data();
    indirection_arrays(1)=&internal_nodes;
    indirection_arrays(2)=&boundary_nodes;
    return indirection_arrays;}

    ARRAY<ARRAY<int>*> Map_Ghost_Nodes() const
    {ARRAY<ARRAY<int>*> indirection_arrays(1);Node_Iterator_Data();
    indirection_arrays(1)=&ghost_nodes;
    return indirection_arrays;}

    ARRAY<ARRAY<int>*> Map_Ghost_And_Boundary_Nodes() const
    {ARRAY<ARRAY<int>*> indirection_arrays(2);Node_Iterator_Data();
    indirection_arrays(1)=&ghost_nodes;
    indirection_arrays(2)=&boundary_nodes;
    return indirection_arrays;}

    ARRAY<ARRAY<int>*> Map_Individual_Side_Ghost_Nodes(const int side) const
    {ARRAY<ARRAY<int>*> indirection_arrays(1);Node_Iterator_Data();
    indirection_arrays(1)=&individual_side_ghost_nodes(side);
    return indirection_arrays;}

    ARRAY<ARRAY<int>*> Map_Individual_Side_Boundary_Nodes(const int side) const
    {ARRAY<ARRAY<int>*> indirection_arrays(1);Node_Iterator_Data();
    indirection_arrays(1)=&individual_side_boundary_nodes(side);
    return indirection_arrays;}

    ARRAY<ARRAY<int>*> Map_Individual_Side_Ghost_And_Boundary_Nodes(const int side) const
    {ARRAY<ARRAY<int>*> indirection_arrays(2);Node_Iterator_Data();
    indirection_arrays(1)=&individual_side_ghost_nodes(side);
    indirection_arrays(2)=&individual_side_boundary_nodes(side);
    return indirection_arrays;}

    ARRAY<ARRAY<int>*> Map_All_Faces() const
    {ARRAY<ARRAY<int>*> indirection_arrays(3);Face_Iterator_Data();
    indirection_arrays(1)=&internal_faces;
    indirection_arrays(2)=&boundary_faces;
    indirection_arrays(3)=&ghost_faces;
    return indirection_arrays;}

    ARRAY<ARRAY<int>*> Map_Interior_Faces() const
    {ARRAY<ARRAY<int>*> indirection_arrays(1);Face_Iterator_Data();
    indirection_arrays(1)=&internal_faces;
    return indirection_arrays;}

    ARRAY<ARRAY<int>*> Map_Boundary_Faces() const
    {ARRAY<ARRAY<int>*> indirection_arrays(1);Face_Iterator_Data();
    indirection_arrays(1)=&boundary_faces;
    return indirection_arrays;}

    ARRAY<ARRAY<int>*> Map_Regular_Faces() const
    {ARRAY<ARRAY<int>*> indirection_arrays(2);Face_Iterator_Data();
    indirection_arrays(1)=&internal_faces;
    indirection_arrays(2)=&boundary_faces;
    return indirection_arrays;}

    ARRAY<ARRAY<int>*> Map_Ghost_Faces() const
    {ARRAY<ARRAY<int>*> indirection_arrays(1);Face_Iterator_Data();
    indirection_arrays(1)=&ghost_faces;
    return indirection_arrays;}

    ARRAY<ARRAY<int>*> Map_Individual_Side_Ghost_Faces(const int side) const
    {ARRAY<ARRAY<int>*> indirection_arrays(1);Face_Iterator_Data();
    indirection_arrays(1)=&individual_side_ghost_faces(side);
    return indirection_arrays;}

    ARRAY<ARRAY<int>*> Map_Individual_Side_Boundary_Faces(const int side) const
    {ARRAY<ARRAY<int>*> indirection_arrays(1);Face_Iterator_Data();
    indirection_arrays(1)=&individual_side_boundary_faces(side);
    return indirection_arrays;}

    ARRAY<ARRAY<int>*> Map_Individual_Side_Ghost_And_Boundary_Faces(const int side) const
    {ARRAY<ARRAY<int>*> indirection_arrays(2);Face_Iterator_Data();
    indirection_arrays(1)=&individual_side_ghost_faces(side);
    indirection_arrays(2)=&individual_side_boundary_faces(side);
    return indirection_arrays;}

    ARRAY<ARRAY<int>*> Map_Individual_Domain_Faces(const int side) const
    {ARRAY<ARRAY<int>*> indirection_arrays(2);Face_Iterator_Data();
    indirection_arrays(1)=&individual_side_domain_ghost_faces(side);
    indirection_arrays(2)=&individual_side_boundary_faces(side);
    return indirection_arrays;}

    int Node_Indices(const int ghost_cells=0) const
    {return number_of_nodes;}

    int Cell_Indices(const int ghost_cells=0) const
    {return number_of_cells;}

    int Block_Indices(const int ghost_cells=0) const
    {return number_of_cells;}

    int Face_Indices(const int ghost_cells=0) const
    {return number_of_faces;}

    int Domain_Indices(const int ghost_cells=0) const
    {return number_of_cells;}

    bool Inside_Domain(const int& cell_index) const
    {return uniform_grid.Domain().Lazy_Inside(cells(cell_index)->Center());}

//#####################################################################
    void Initialize(const GRID<TV> uniform_grid_input,const int maximum_depth_input,const int number_of_ghost_cells_input,const bool use_nodes,const bool use_faces);
    void Tree_Topology_Changed();
    BINTREE_CELL<T>* Leaf_Cell(const TV& location,const T thickness=1e-3) const;
    BINTREE_CELL<T>* Clamped_Leaf_Cell(const TV& location,const T thickness=1e-3) const;
    const BINTREE_CELL<T>* Leaf_Cell_By_Neighbor_Path(const BINTREE_CELL<T>* start_cell,const TV& location,const T tolerance=1e-3,const int max_iterations=10) const;
    BINTREE_CELL<T>* Leaf_Cell_By_Neighbor_Path(BINTREE_CELL<T>* start_cell,const TV& location,const T tolerance=1e-3,const int max_iterations=10) const;
    bool Inside_Cell(const BINTREE_CELL<T>* cell,const TV& location) const;
    bool Inside_Thickened_Cell(const BINTREE_CELL<T>* cell,const TV& location,const T thickness=1e-3) const;
    BINTREE_CELL<T>* Inside_Offspring(BINTREE_CELL<T>* cell,const TV& location,const T tolerance=1e-3) const;
    const BINTREE_CELL<T>* Inside_Offspring(const BINTREE_CELL<T>* cell,const TV& location,const T tolerance=1e-3) const;
    void Refine_Cell(const int max_depth, BINTREE_CELL<T>* cell,const TV& location,ARRAY<BINTREE_CELL<T>*>* new_cells,ARRAY<PAIR<BINTREE_CELL<T>*,int> >* new_nodes, 
        ARRAY<PAIR<BINTREE_CELL<T>*,int> >* new_faces,const BINTREE_GRID<T>* grid);
    void Compact_Array_Indices(ARRAY<int>* cell_mapping_array,ARRAY<int>* node_mapping_array,ARRAY<int>* face_mapping_array);
    void Get_Cells_Intersecting_Box(const RANGE<TV>& box,ARRAY<BINTREE_CELL<T>*>& intersecting_cells) const;
    void Refine_Cells_Intersecting_Box(const RANGE<TV>& box,ARRAY<BINTREE_CELL<T>*>& refined_cells,const int refinement_depth=0);
    void Node_Iterator_Data() const;
    void Face_Iterator_Data() const;
private:
    void Calculate_Cell_Pointer_From_Index_Array(ARRAY<BINTREE_CELL<T>*>& cell_pointer_from_index) const;
    void Calculate_Neighbors_Array(ARRAY<VECTOR<BINTREE_CELL<T>*,2> >& neighbors) const;
    void Calculate_Node_Neighbors_Array(ARRAY<VECTOR<int,2> >& node_neighbors) const;
    void Calculate_Fully_Refined_Block_Array(ARRAY<bool>& fully_refined_block) const;
public:
    template<class TV> void Enslave_T_Junction_Nodes(ARRAY<TV>* nodes,int depth_to_enforce) {} // No T-Junction nodes in 1-D, so this is a no-op.
    void Check_Tree_Consistency(bool check_cells,bool check_nodes,bool check_faces,bool check_neighbor_structures);
//#####################################################################
};
template<class T> const T BINTREE_GRID<T>::one_over_number_of_children_per_cell=(T).5;
template<class T> const T BINTREE_GRID<T>::one_over_number_of_nodes_per_cell=(T).5;
template<class T> const T BINTREE_GRID<T>::one_over_number_of_nodes_per_face=(T)1;
}
#endif
#endif

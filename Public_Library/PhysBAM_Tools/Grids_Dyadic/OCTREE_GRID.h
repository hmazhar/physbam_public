#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2003-2006, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OCTREE_GRID
//#####################################################################
#ifndef __OCTREE_GRID__
#define __OCTREE_GRID__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
namespace PhysBAM{

template<class T> class MAP_OCTREE_MESH;
template<class T> class OCTREE_CELL_HELPER;
template<class TV> class GRID;
template<class T_GRID> class DYADIC_GRID_ITERATOR_NODE;
template<class T_GRID> class DYADIC_GRID_ITERATOR_CELL;
template<class T_GRID> class DYADIC_GRID_ITERATOR_FACE;
template<class TV> struct DYADIC_TAG;
template<class T_GRID> class BLOCK_DYADIC;

template<class T_input>
class OCTREE_GRID:public NONCOPYABLE
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    static const int number_of_children_per_cell=8;
    static const T one_over_number_of_children_per_cell;
    static const int number_of_cells_per_block=8;
    static const int number_of_faces_per_block=36;
    static const int number_of_faces_per_cell=6;
    static const int number_of_nodes_per_cell=8;
    enum {number_of_cells_per_node=8}; // workaround for ICC bug
    static const T one_over_number_of_nodes_per_cell;
    static const int number_of_nodes_per_face=4;
    static const T one_over_number_of_nodes_per_face;
    static const int number_of_neighbors_per_node=6;
    static const int number_of_neighbors_per_cell=6;
    static const int number_of_neighbors_per_face=6;
    static const int number_of_complete_one_ring_neighbors_face_per_cell=6;
    static const int number_of_complete_one_ring_neighbors_edge_per_cell=12;
    static const int number_of_complete_one_ring_neighbors_node_per_cell=8;
    static const int dimension=3;
    static const char* name;
    typedef T SCALAR;
    typedef VECTOR<T,3> VECTOR_T;
    typedef DYADIC_TAG<VECTOR_T> GRID_TAG;
    typedef BLOCK_DYADIC<OCTREE_GRID<T> > BLOCK;

    typedef OCTREE_GRID<T> GRID_T;
    typedef VECTOR<int,3> VECTOR_INT;
    typedef int INDEX;
    typedef MAP_OCTREE_MESH<T> MAP_MESH;
    typedef OCTREE_CELL<T> CELL;
    typedef OCTREE_CELL_HELPER<T> CELL_HELPER;
    typedef GRID<TV> UNIFORM_GRID;

    typedef DYADIC_GRID_ITERATOR_FACE<GRID_T> FACE_ITERATOR;
    typedef DYADIC_GRID_ITERATOR_NODE<GRID_T> NODE_ITERATOR;
    typedef DYADIC_GRID_ITERATOR_CELL<GRID_T> CELL_ITERATOR;

    ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> > cells;
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
    ARRAY<OCTREE_CHILDREN<T>*> allocated_children;
    mutable bool cell_pointer_from_index_up_to_date;mutable ARRAY<OCTREE_CELL<T>*> cell_pointer_from_index;
    mutable bool neighbors_up_to_date;mutable ARRAY<VECTOR<OCTREE_CELL<T>*,6> > neighbors;
    mutable bool node_neighbors_up_to_date;mutable ARRAY<VECTOR<int,6> > node_neighbors;
public:
    mutable bool fully_refined_block_up_to_date;mutable ARRAY<bool> fully_refined_block;

    OCTREE_GRID();
    ~OCTREE_GRID();

    T Minimum_Edge_Length() const
    {return minimum_cell_size;}

    TV Minimum_Cell_DX() const
    {return TV(minimum_cell_size,minimum_cell_size,minimum_cell_size);}

    TV One_Over_Minimum_Cell_DX() const
    {return TV(one_over_minimum_cell_size,one_over_minimum_cell_size,one_over_minimum_cell_size);}

    TV Minimum_Cell_DX_Over_Two() const
    {return TV(minimum_cell_size_over_two,minimum_cell_size_over_two,minimum_cell_size_over_two);}

    void Compute_Minimum_Cell_Constants(int proposed_maximum_depth)
    {minimum_cell_size=min(uniform_grid.dX.x,uniform_grid.dX.y,uniform_grid.dX.z)/(1<<(proposed_maximum_depth-1));
    one_over_minimum_cell_size=(T)1./minimum_cell_size;
    minimum_cell_size_over_two=(T).5*minimum_cell_size;
    one_over_minimum_cell_size_over_two=(T)1./minimum_cell_size_over_two;}

    void Set_Maximum_Depth(const int maximum_depth_input)
    {maximum_depth=maximum_depth_input;Compute_Minimum_Cell_Constants(maximum_depth);}

    void Update_Maximum_Depth()
    {maximum_depth=0;for(int i=1;i<=uniform_grid.numbers_of_cells.x;i++)for(int j=1;j<=uniform_grid.numbers_of_cells.y;j++)for(int ij=1;ij<=uniform_grid.numbers_of_cells.z;ij++)maximum_depth=max(maximum_depth,cells(i,j,ij)->Maximum_Relative_Depth());}

    TV Clamp(const TV& location) const
    {return uniform_grid.Clamp(location);}

    const RANGE<TV>& Domain() const
    {return uniform_grid.domain;}

    bool Outside(const TV& location) const
    {return uniform_grid.Outside(location);}

    static TV Face_Corner_To_Opposite_Corner_Normalized_Vector(const int axis,const int corner)
    {static const TV vector_table[3][4]={
        {TV(0,-(T)one_over_root_two,-(T)one_over_root_two),TV(0,(T)one_over_root_two,-(T)one_over_root_two),
         TV(0,-(T)one_over_root_two,(T)one_over_root_two),TV(0,(T)one_over_root_two,(T)one_over_root_two)},
        {TV(-(T)one_over_root_two,0,-(T)one_over_root_two),TV((T)one_over_root_two,0,-(T)one_over_root_two),
         TV(-(T)one_over_root_two,0,(T)one_over_root_two),TV((T)one_over_root_two,0,(T)one_over_root_two)},
        {TV(-(T)one_over_root_two,-(T)one_over_root_two,0),TV((T)one_over_root_two,-(T)one_over_root_two,0),
         TV(-(T)one_over_root_two,(T)one_over_root_two,0),TV((T)one_over_root_two,(T)one_over_root_two,0)}};
    assert(1<=axis&&axis<=3);assert(1<=corner&&corner<=4);return vector_table[axis-1][corner-1];}

    OCTREE_CELL<T>* Complete_Minimal_One_Ring_Neighbor_Face(const int i,OCTREE_CELL<T>* cell) const
    {assert(1<=i&&i<=6);
    OCTREE_CELL<T> *temp=neighbors(cell->Cell())(i);if(temp&&temp->Depth_Of_This_Cell()==maximum_depth)return temp;return 0;}

    OCTREE_CELL<T>* Complete_Minimal_One_Ring_Neighbor_Edge(const int i,OCTREE_CELL<T>* cell) const
    {assert(1<=i&&i<=12);OCTREE_CELL<T> *temp,*temp2;int lookup[12][2]={{1,3},{2,3},{1,4},{2,4},{1,5},{2,5},{1,6},{2,6},{3,5},{3,6},{4,5},{4,6}};
    temp=neighbors(cell->Cell())(lookup[i-1][0]);
    if(temp&&temp->Depth_Of_This_Cell()==maximum_depth){temp2=neighbors(temp->Cell())(lookup[i-1][1]);if(temp2&&temp2->Depth_Of_This_Cell()==maximum_depth)return temp2;}
    temp=neighbors(cell->Cell())(lookup[i-1][1]);
    if(temp&&temp->Depth_Of_This_Cell()==maximum_depth){temp2=neighbors(temp->Cell())(lookup[i-1][0]);if(temp2&&temp2->Depth_Of_This_Cell()==maximum_depth)return temp2;}
    return 0;}

    OCTREE_CELL<T>* Complete_Minimal_One_Ring_Neighbor_Node(const int i,OCTREE_CELL<T>* cell) const
    {assert(1<=i&&i<=8);OCTREE_CELL<T> *temp,*temp2;
    int lookup[8][3][2]={{{1,5},{5,1},{9,3}},{{2,5},{5,2},{10,3}},{{1,6},{6,1},{11,3}},{{2,6},{6,2},{12,3}},{{3,5},{7,1},{9,4}},{{4,5},{7,2},{10,4}},{{3,6},{8,1},{11,4}},{{4,6},{8,2},{12,4}}};
    for(int attempt=0;attempt<3;attempt++){
        temp=Complete_Minimal_One_Ring_Neighbor_Edge(lookup[i-1][attempt][0],cell);
        if(temp){temp2=neighbors(temp->Cell())(lookup[i-1][attempt][1]);if(temp2&&temp2->Depth_Of_This_Cell()==maximum_depth)return temp2;}}
    return 0;}

    TV Node_Location(const int global_node_index) const
    {return Node_Location(nodes(global_node_index),cell_pointer_from_index(node_iterator_deepest_cells(global_node_index)));}

    TV Node_Location(const int node_index,const OCTREE_CELL<T>* cell) const
    {assert(node_index >=0 && node_index <= 7);
    TV center(cell->Center()),DX_over_two((T).5*cell->DX());
    switch(node_index){
        case 0:return TV(center.x-DX_over_two.x,center.y-DX_over_two.y,center.z-DX_over_two.z);
        case 1:return TV(center.x+DX_over_two.x,center.y-DX_over_two.y,center.z-DX_over_two.z);
        case 2:return TV(center.x-DX_over_two.x,center.y+DX_over_two.y,center.z-DX_over_two.z);
        case 3:return TV(center.x+DX_over_two.x,center.y+DX_over_two.y,center.z-DX_over_two.z);
        case 4:return TV(center.x-DX_over_two.x,center.y-DX_over_two.y,center.z+DX_over_two.z);
        case 5:return TV(center.x+DX_over_two.x,center.y-DX_over_two.y,center.z+DX_over_two.z);
        case 6:return TV(center.x-DX_over_two.x,center.y+DX_over_two.y,center.z+DX_over_two.z);
        default:return TV(center.x+DX_over_two.x,center.y+DX_over_two.y,center.z+DX_over_two.z);}}

    TV Face_Location(const int global_face_index) const
    {return Face_Location(faces(global_face_index),cell_pointer_from_index(face_iterator_deepest_cells(global_face_index)));}

    char Face_Axis(const int global_face_index) const
    {static const int axis[6]={0,0,1,1,2,2};return axis[faces(global_face_index)];}

    TV Face_Location(const int face_index,const OCTREE_CELL<T>* cell) const
    {assert(face_index >=0 && face_index <= 5);
    TV center(cell->Center()),DX_over_two((T).5*cell->DX());
    switch(face_index){
        case 0:return TV(center.x-DX_over_two.x,center.y,center.z);
        case 1:return TV(center.x+DX_over_two.x,center.y,center.z);
        case 2:return TV(center.x,center.y-DX_over_two.y,center.z);
        case 3:return TV(center.x,center.y+DX_over_two.y,center.z);
        case 4:return TV(center.x,center.y,center.z-DX_over_two.z);
        default:return TV(center.x,center.y,center.z+DX_over_two.z);}}

    TV Cell_Location(const int global_cell_index) const
    {return Cell_Pointer_From_Index()(global_cell_index)->Center();}

    ARRAY<VECTOR<OCTREE_CELL<T>*,6> >& Neighbors() const
    {if(!neighbors_up_to_date){Calculate_Neighbors_Array(neighbors);neighbors_up_to_date=true;}
    return neighbors;}

    ARRAY<VECTOR<int,6> >& Node_Neighbors() const
    {if(!node_neighbors_up_to_date){Calculate_Node_Neighbors_Array(node_neighbors);node_neighbors_up_to_date=true;}
    return node_neighbors;}

    ARRAY<OCTREE_CELL<T>*>& Cell_Pointer_From_Index() const
    {if(!cell_pointer_from_index_up_to_date){Calculate_Cell_Pointer_From_Index_Array(cell_pointer_from_index);cell_pointer_from_index_up_to_date=true;}
    return cell_pointer_from_index;}

    ARRAY<bool>& Fully_Refined_Block() const
    {if(!fully_refined_block_up_to_date){Calculate_Fully_Refined_Block_Array(fully_refined_block);fully_refined_block_up_to_date=true;}
    return fully_refined_block;}

    const OCTREE_CELL<T>* Base_Cell(const TV& X) const
    {assert(fully_refined_block_up_to_date);
    OCTREE_CELL<T>* base_cell=Leaf_Cell(X-Minimum_Cell_DX_Over_Two());
    if(base_cell&&fully_refined_block(base_cell->Cell()))return base_cell;else return 0;}

    const OCTREE_CELL<T>* Base_Cell_By_Neighbor_Path(const OCTREE_CELL<T>* start_cell,const TV& location,const T tolerance=1e-3,const int max_iterations=10) const
    {assert(fully_refined_block_up_to_date);
    const OCTREE_CELL<T>* base_cell=Leaf_Cell_By_Neighbor_Path(start_cell,location-Minimum_Cell_DX_Over_Two(),tolerance,max_iterations);
    if(base_cell&&fully_refined_block(base_cell->Cell()))return base_cell;return 0;}

    OCTREE_CELL<T>* Base_Cell_By_Neighbor_Path(OCTREE_CELL<T>* start_cell,const TV& location,const T tolerance=1e-3,const int max_iterations=10) const
    {assert(fully_refined_block_up_to_date);
    OCTREE_CELL<T>* base_cell=Leaf_Cell_By_Neighbor_Path(start_cell,location-Minimum_Cell_DX_Over_Two(),tolerance,max_iterations);
    if(base_cell&&fully_refined_block(base_cell->Cell()))return base_cell;return 0;}

    const OCTREE_CELL<T>* Cell_From_Cell_Block(const OCTREE_CELL<T>* base_cell,const int number) // 0 to 7 in 3D
    {assert(fully_refined_block_up_to_date && fully_refined_block(base_cell->Cell()));const OCTREE_CELL<T>* cell=base_cell;
    if(number&1)cell=neighbors(cell->Cell())(2);if(number&2)cell=neighbors(cell->Cell())(4);if(number&4)cell=neighbors(cell->Cell())(6);
    return cell;}

    void All_Cells_In_Cell_Block(const OCTREE_CELL<T>* base_cell,const OCTREE_CELL<T>* cells[number_of_cells_per_block]) const
    {assert(fully_refined_block_up_to_date && fully_refined_block(base_cell->Cell()));cells[0]=base_cell;int lookup[][2]={{0,0},{2,0},{4,0},{4,1},{6,0},{6,1},{6,2},{6,3}};
    for(int i=1;i<number_of_cells_per_block;i++)cells[i]=neighbors(cells[lookup[i][1]]->Cell())(lookup[i][0]);}

    int First_Face_Index_In_Cell(const int axis,const int cell_index) const
    {return cell_pointer_from_index(cell_index)->Face(2*axis);}

    int Second_Face_Index_In_Cell(const int axis,const int cell_index) const
    {return cell_pointer_from_index(cell_index)->Face(2*axis+1);}

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

//#####################################################################
    void Initialize(const GRID<TV> uniform_grid_input,const int maximum_depth_input,const int number_of_ghost_cells_input,const bool use_nodes,const bool use_faces);
    void Tree_Topology_Changed();
    OCTREE_CELL<T>* Leaf_Cell(const TV& location,const T thickness=1e-3) const;
    OCTREE_CELL<T>* Clamped_Leaf_Cell(const TV& location,const T thickness=1e-3) const;
    const OCTREE_CELL<T>* Leaf_Cell_By_Neighbor_Path(const OCTREE_CELL<T>* start_cell,const TV& location,const T tolerance=1e-3,const int max_iterations=10) const;
    OCTREE_CELL<T>* Leaf_Cell_By_Neighbor_Path(OCTREE_CELL<T>* start_cell,const TV& location,const T tolerance=1e-3,const int max_iterations=10) const;
    bool Inside_Cell(const OCTREE_CELL<T>* cell,const TV& location) const;
    bool Inside_Thickened_Cell(const OCTREE_CELL<T>* cell,const TV& location,const T thickness=1e-3) const;
    OCTREE_CELL<T>* Inside_Offspring(OCTREE_CELL<T>* cell,const TV& location,const T thickness=1e-3) const;
    const OCTREE_CELL<T>* Inside_Offspring(const OCTREE_CELL<T>* cell,const TV& location,const T thickness=1e-3) const;
    void Refine_Cell(const int max_depth,OCTREE_CELL<T>* cell,const TV& location,ARRAY<OCTREE_CELL<T>*>* new_cells,ARRAY<PAIR<OCTREE_CELL<T>*,int> >* new_nodes,
        ARRAY<PAIR<OCTREE_CELL<T>*,int> >* new_faces,const OCTREE_GRID<T>* grid);
    void Compact_Array_Indices(ARRAY<int>* cell_mapping_array,ARRAY<int>* node_mapping_array,ARRAY<int>* face_mapping_array);
    void Get_Cells_Intersecting_Box(const RANGE<TV>& box,ARRAY<OCTREE_CELL<T>*>& intersecting_cells);
    void Refine_Cells_Intersecting_Box(const RANGE<TV>& box,ARRAY<OCTREE_CELL<T>*>& refined_cells,const int refinement_depth=0);
    void Node_Iterator_Data() const;
    void Face_Iterator_Data() const;
private:
    void Calculate_Cell_Pointer_From_Index_Array(ARRAY<OCTREE_CELL<T>*>& cell_pointer_from_index) const;
    void Calculate_Neighbors_Array(ARRAY<VECTOR<OCTREE_CELL<T>*,6> >& neighbors) const;
    void Calculate_Node_Neighbors_Array(ARRAY<VECTOR<int,6> >& node_neighbors) const;
    void Calculate_Fully_Refined_Block_Array(ARRAY<bool>& fully_refined_block) const;
public:
    template<class TV> void Enslave_T_Junction_Nodes(ARRAY<TV>* nodes,int depth_to_enforce);
    void Check_Tree_Consistency(bool check_cells,bool check_nodes,bool check_faces,bool check_neighbor_structures);
    void Print_Storage_Requirement();
    void Move_Contents_Left_Two();
    void Move_Contents_Right_Two();
    void Move_Contents_Down_Two();
    void Move_Contents_Up_Two();
    void Move_Contents_Forward_Two();
    void Move_Contents_Backward_Two();
//#####################################################################
};
template<class T> const T OCTREE_GRID<T>::one_over_number_of_children_per_cell=(T).125;
template<class T> const T OCTREE_GRID<T>::one_over_number_of_nodes_per_cell=(T).125;
template<class T> const T OCTREE_GRID<T>::one_over_number_of_nodes_per_face=(T).25;
}
#endif
#endif

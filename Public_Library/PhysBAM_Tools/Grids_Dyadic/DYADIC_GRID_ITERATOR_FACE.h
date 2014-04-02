//#####################################################################
// Copyright 2005, Ron Fedkiw, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DYADIC_GRID_ITERATOR_FACE 
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __DYADIC_GRID_ITERATOR_FACE__
#define __DYADIC_GRID_ITERATOR_FACE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_BINTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
namespace PhysBAM{

template<class T_GRID>
class DYADIC_GRID_ITERATOR_FACE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::MAP_MESH MAP_MESH;typedef typename T_GRID::CELL CELL;
    typedef typename T_GRID::INDEX_FACE INDEX_FACE;
public:
    const T_GRID& grid;
    ARRAY<CELL*>& cell_pointer_from_index;
    ARRAY<VECTOR<CELL*,T_GRID::number_of_neighbors_per_cell> >& neighbors;
    ARRAY<int>* indirection_array;
    int indirection_index;
    int current_index;
    ARRAY<ARRAY<int>*> indirection_arrays;
    int current_indirection_array;

    DYADIC_GRID_ITERATOR_FACE(const T_GRID& grid_input,const ARRAY<ARRAY<int>*>& indirection_arrays_input)
        :grid(grid_input),cell_pointer_from_index(grid.Cell_Pointer_From_Index()),neighbors(grid.Neighbors())
    {
        grid.Face_Iterator_Data();
        indirection_arrays=indirection_arrays_input;
        indirection_index=0;current_index=0;
        indirection_array=indirection_arrays(current_indirection_array=1);       
        Next();
    }

    DYADIC_GRID_ITERATOR_FACE(const T_GRID& grid_input,const int current_index_input)
        :grid(grid_input),cell_pointer_from_index(grid.Cell_Pointer_From_Index()),neighbors(grid.Neighbors()),current_index(current_index_input)
    {
        grid.Face_Iterator_Data();
    }

    void Next()
    {assert(indirection_arrays.m);if(indirection_index<indirection_array->m){current_index=(*indirection_array)(++indirection_index);return;}
    else if(++current_indirection_array>indirection_arrays.m)return;
    indirection_array=indirection_arrays(current_indirection_array);current_index=(*indirection_array)(indirection_index=1);}

    bool Valid() const
    {assert(indirection_arrays.m);return current_indirection_array<=indirection_arrays.m;}

    void Set_Index(const int index) // this function destroys the ability to call Next and Valid. Should only be used to manually set the face index
    {current_index=index;}

    int Face_Index_In_Cell() const // returns a number between 0 and number_of_faces_per_cell-1. Which of the faces in the deepest cell is being iterated over
    {return grid.faces(Face_Index());}

    int Other_Face_Index_In_Cell() const // returns a number between 0 and number_of_faces_per_cell-1. Which of the faces in the deepest cell is being iterated over
    {static const int other_face[6]={1,0,3,2,5,4};return other_face[Face_Index_In_Cell()];} // the tal

    int Face_Index() const // returns a number between 1 and grid.number_of_faces. The global face value
    {return current_index;}

    INDEX_FACE Full_Index() const
    {return INDEX_FACE(Axis(),current_index,First_Cell_Index(),Second_Cell_Index());}

    int Other_Face_Index() const // returns a number between 1 and grid.number_of_faces. The global face value
    {return Other_Cell()->Face(Other_Face_Index_In_Cell());}

    int Axis() const
    {static const int axis[6]={0,0,1,1,2,2};return axis[Face_Index_In_Cell()];}

    int Face_Node(const int node) const
    {static const bool deepest_cell_position[6]={1,0,1,0,1,0};assert(node>=0&&node<T_GRID::number_of_nodes_per_face);
    return Deepest_Cell()->Node(MAP_MESH::nodes_by_axis[Axis()][deepest_cell_position[Face_Index_In_Cell()]][node]);}

    int Other_Face_Node(const int node) const
    {static const bool other_cell_position[6]={0,1,0,1,0,1};assert(node>=0&&node<T_GRID::number_of_nodes_per_face);
    return Other_Cell()->Node(MAP_MESH::nodes_by_axis[Axis()][other_cell_position[Face_Index_In_Cell()]][node]);}

    void Face_Nodes(int* nodes) const
    {for(int i=0;i<T_GRID::number_of_nodes_per_face;i++) nodes[i]=Face_Node(i);}

    void Other_Face_Nodes(int* nodes) const
    {for(int i=0;i<T_GRID::number_of_nodes_per_face;i++) nodes[i]=Other_Face_Node(i);}

    TV Location() const
    {return grid.Face_Location(Face_Index_In_Cell(),Deepest_Cell());}

    T Face_DX() const
    {return Deepest_Cell()->DX().x;} // assumes square cells

    T Other_Face_DX() const
    {return Other_Cell()->DX().x;} // assumes square cells

    T Face_Size() const
    {return Deepest_Cell()->Face_Size();} // assumes square cells

    T Other_Face_Size() const
    {return Other_Cell()->Face_Size();} // assumes square cells

    CELL* Deepest_Cell() const
    {return cell_pointer_from_index(grid.face_iterator_deepest_cells(current_index));}

    int Deepest_Cell_Index() const
    {return grid.face_iterator_deepest_cells(current_index);}

    CELL* Other_Cell() const
    {return neighbors(Deepest_Cell_Index())(Face_Index_In_Cell()+1);}

    int Other_Cell_Index() const
    {return neighbors(Deepest_Cell_Index())(Face_Index_In_Cell()+1)->Cell();}

    int Deepest_Cell_Number() const
    {static const int deepest_cell_number[6]={1,0,1,0,1,0};return deepest_cell_number[Face_Index_In_Cell()];}

    bool Deepest_Cell_Is_First_Cell() const
    {static const bool is_first_cell[6]={0,1,0,1,0,1};return is_first_cell[Face_Index_In_Cell()];}

    CELL* Cell(const int cell) const // returns first cell if cell==0 returns second cell if cell==1
    {assert(0<=cell&&cell<=1);
    if(cell==0)return First_Cell();
    else return Second_Cell();}

    CELL* First_Cell() const
    {if(Deepest_Cell_Is_First_Cell())return Deepest_Cell();
    else return Other_Cell();}

    int First_Cell_Index() const
    {if(Deepest_Cell_Is_First_Cell())return Deepest_Cell_Index();
    else return Other_Cell_Index();}

    T First_Cell_DX() const
    {if(Deepest_Cell_Is_First_Cell())return Deepest_Cell()->DX().x;
    else return Other_Cell()->DX().x;}

    CELL* Second_Cell() const
    {if(!Deepest_Cell_Is_First_Cell())return Deepest_Cell();
    else return Other_Cell();}

    int Second_Cell_Index() const
    {if(!Deepest_Cell_Is_First_Cell())return Deepest_Cell_Index();
    else return Other_Cell_Index();}

    T Second_Cell_DX() const
    {if(!Deepest_Cell_Is_First_Cell())return Deepest_Cell()->DX().x;
    else return Other_Cell()->DX().x;}

    void Unordered_Cell_Indices_Touching_Face(int& cell1,int& cell2) const
    {cell1=Deepest_Cell_Index();cell2=Other_Cell_Index();}

    int Minimal_Neighbor(const int neighbor)
    {assert(neighbor>=0&&(neighbor<4||(neighbor<6&&T_GRID::dimension==3)));
    if(neighbor==Axis()*2){CELL* cell=First_Cell();if(cell&&cell->Depth_Of_This_Cell()==grid.maximum_depth)return cell->Face(neighbor);else return 0;}
    else if(neighbor==Axis()*2+1){CELL* cell=Second_Cell();if(cell&&cell->Depth_Of_This_Cell()==grid.maximum_depth)return cell->Face(neighbor);else return 0;}
    else{CELL* cell=neighbors(Deepest_Cell_Index())(neighbor+1);if(cell&&cell->Depth_Of_This_Cell()==grid.maximum_depth)return cell->Face(Face_Index_In_Cell());else return 0;}}

//#####################################################################
};
}
#endif
#endif

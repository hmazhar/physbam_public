//#####################################################################
// Copyright 2005, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __FACE_LOOKUP_COLLIDABLE_DYADIC__
#define __FACE_LOOKUP_COLLIDABLE_DYADIC__

#include <PhysBAM_Geometry/Grids_Dyadic_Collisions/GRID_BASED_COLLISION_GEOMETRY_DYADIC.h>
namespace PhysBAM{

template<class T_GRID,class T_NESTED_LOOKUP> // T_NESTED_LOOKUP=FACE_LOOKUP_DYADIC<T_GRID>
class FACE_LOOKUP_COLLIDABLE_DYADIC
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;typedef typename T_GRID::CELL T_CELL;
public:
    template<class T_NESTED_LOOKUP_2> struct REBIND_NESTED_LOOKUP{typedef FACE_LOOKUP_COLLIDABLE_DYADIC<T_GRID,T_NESTED_LOOKUP_2> TYPE;};
    typedef T ELEMENT;
    typedef T_NESTED_LOOKUP NESTED_LOOKUP;

    class LOOKUP;

    const T_NESTED_LOOKUP& nested_face_lookup;
    const T_GRID_BASED_COLLISION_GEOMETRY& body_list;
    const ARRAY<bool>* valid_value_mask;

    FACE_LOOKUP_COLLIDABLE_DYADIC(const T_NESTED_LOOKUP& nested_face_lookup_input,const T_GRID_BASED_COLLISION_GEOMETRY& body_list_input,const ARRAY<bool>* valid_value_mask_input)
        :nested_face_lookup(nested_face_lookup_input),body_list(body_list_input),valid_value_mask(valid_value_mask_input)
    {}

    const ARRAY<T>& Raw_Data() const
    {return nested_face_lookup.Raw_Data();}

    const T_NESTED_LOOKUP& Nested() const
    {return nested_face_lookup;}

    LOOKUP Starting_Point_Face(const int face) const
    {return LOOKUP(*this,nested_face_lookup.Starting_Point_Face(face));}

    LOOKUP Starting_Point_Cell(const int cell) const
    {return LOOKUP(*this,nested_face_lookup.Starting_Point_Cell(cell));}

    LOOKUP Region(const T phi) const
    {return LOOKUP(*this,nested_face_lookup.Region(phi));}

    class LOOKUP
    {
    private:
        const FACE_LOOKUP_COLLIDABLE_DYADIC<T_GRID,T_NESTED_LOOKUP>& face_lookup;
        typename T_NESTED_LOOKUP::LOOKUP nested_lookup;
    public:
        typedef T ELEMENT;
        mutable bool reference_point_set;
        mutable bool reference_point_inside;
        mutable bool found_valid_point;
        mutable TV reference_point;
        mutable TV object_velocity;
        mutable int cells_in_block[T_GRID::number_of_cells_per_block];

        LOOKUP(const FACE_LOOKUP_COLLIDABLE_DYADIC<T_GRID,T_NESTED_LOOKUP>& face_lookup_input,const typename T_NESTED_LOOKUP::LOOKUP& nested_lookup_input)
            :face_lookup(face_lookup_input),nested_lookup(nested_lookup_input),reference_point_set(false),reference_point_inside(false)
        {}

        void Set_Reference_Point(const TV& reference_point_input) const
        {reference_point=reference_point_input;assert(!reference_point_set && !reference_point_inside);reference_point_set=true;
        COLLISION_GEOMETRY_ID body_id;int aggregate_id;found_valid_point=false;
        const T_CELL* base_cell=face_lookup.body_list.grid.Base_Cell(reference_point);
        if(base_cell && face_lookup.body_list.Inside_Any_Simplex_Of_Any_Body(BLOCK_DYADIC<T_GRID>(face_lookup.body_list.grid,base_cell),reference_point,body_id,aggregate_id)){
            reference_point_inside=true;
            object_velocity=face_lookup.body_list.Object_Velocity(body_id,aggregate_id,reference_point);}
        else{BLOCK_DYADIC<T_GRID> block(face_lookup.body_list.grid,reference_point);block.All_Cell_Indices(cells_in_block);}}

        void Clear_Reference_Point() const
        {assert(reference_point_set);reference_point_set=false;reference_point_inside=false;}

        T operator()(const int face) const
        {T face_velocity=0;assert(reference_point_set);int axis=face_lookup.body_list.grid.Face_Axis(face);
        if(reference_point_inside){found_valid_point=true;return object_velocity[axis+1];}
        else if(face_lookup.body_list.Face_Velocity(axis+1,face,cells_in_block,T_GRID::number_of_cells_per_block,reference_point,face_velocity)){
            found_valid_point=true;return face_velocity;}
        else{
            if(face_lookup.valid_value_mask && !(*face_lookup.valid_value_mask)(face)) return 0;
            else{found_valid_point=true;return nested_lookup(face);}}}
    };

//#####################################################################
};
}
#endif
#endif

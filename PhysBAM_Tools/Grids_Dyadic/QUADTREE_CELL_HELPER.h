//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __QUADTREE_CELL_HELPER__
#define __QUADTREE_CELL_HELPER__

#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_CELL.h>
namespace PhysBAM{

template<class T>
class QUADTREE_CELL_HELPER
{
public:
    //#####################################################################
    // Function Interpolate_Face_Values_From_Direct_Children
    //#####################################################################
    template<class T_FACE_LOOKUP> static void
    Interpolate_Face_Values_From_Direct_Children(QUADTREE_CELL<T>& cell,ARRAY<typename T_FACE_LOOKUP::ELEMENT>& face_values,T_FACE_LOOKUP face_values_lookup)
    {assert(cell.Has_Children());
    {const typename T_FACE_LOOKUP::LOOKUP& lookup=face_values_lookup.Starting_Point_Face(cell.Face(0));face_values(cell.Face(0))=(T).5*(lookup(cell.children->Face(0,0))+lookup(cell.children->Face(2,0)));}
    {const typename T_FACE_LOOKUP::LOOKUP& lookup=face_values_lookup.Starting_Point_Face(cell.Face(1));face_values(cell.Face(1))=(T).5*(lookup(cell.children->Face(1,1))+lookup(cell.children->Face(3,1)));}
    {const typename T_FACE_LOOKUP::LOOKUP& lookup=face_values_lookup.Starting_Point_Face(cell.Face(2));face_values(cell.Face(2))=(T).5*(lookup(cell.children->Face(0,2))+lookup(cell.children->Face(1,2)));}
    {const typename T_FACE_LOOKUP::LOOKUP& lookup=face_values_lookup.Starting_Point_Face(cell.Face(3));face_values(cell.Face(3))=(T).5*(lookup(cell.children->Face(2,3))+lookup(cell.children->Face(3,3)));}}

//#####################################################################
};
}
#endif
#endif

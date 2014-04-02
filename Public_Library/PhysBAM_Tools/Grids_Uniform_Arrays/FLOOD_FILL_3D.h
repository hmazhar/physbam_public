//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Eran Guendelman, Frank Losasso, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLOOD_FILL_3D
//#####################################################################
#ifndef __FLOOD_FILL_3D__
#define __FLOOD_FILL_3D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/STACK.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <climits>
namespace PhysBAM{

class FLOOD_FILL_3D
{
    typedef VECTOR<int,3> TV_INT;
private:
    STACK<TV_INT> flood_fill_stack;
    bool optimize_fill_for_single_cell_regions;
    TV_INT last_uncolored_node;

public:
    FLOOD_FILL_3D();
    ~FLOOD_FILL_3D();

    void Optimize_Fill_For_Single_Cell_Regions(const bool optimize=true)
    {optimize_fill_for_single_cell_regions=optimize;}

    int Flood_Fill(ARRAYS_ND_BASE<TV_INT>& colors,const ARRAY<bool,FACE_INDEX<3> >& edge_is_blocked,ARRAY<bool>* color_touches_uncolorable_node=0)
    {return Flood_Fill(colors,edge_is_blocked.Component(1),edge_is_blocked.Component(2),edge_is_blocked.Component(3),color_touches_uncolorable_node);}

    // colors should be initialized by the user with 0's where colors will be filled and negative values for nodes which will not be colored.
    // -1 is the distinguished uncolorable node (which color_touches_uncolorable_node refers to)
    int Flood_Fill(ARRAYS_ND_BASE<TV_INT>& colors,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_x,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_y,
        const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_z,ARRAY<bool>* color_touches_uncolorable_node=0);
    void Fill_Single_Cell_Regions(ARRAYS_ND_BASE<TV_INT>& colors,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_x,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_y,
        const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_z,int& fill_color);
    bool Find_Uncolored_Node(const ARRAYS_ND_BASE<TV_INT>& colors,TV_INT& node_index);
    void Flood_Fill_From_Seed_Node(ARRAYS_ND_BASE<TV_INT>& colors,const int fill_color,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_x,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_y,
        const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_z,bool& touches_uncolorable_node,const TV_INT& seed_node);
    void Flood_Fill_Node(ARRAYS_ND_BASE<TV_INT>& colors,const int fill_color,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_x,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_y,
        const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_z,bool& touches_uncolorable_node,STACK<TV_INT>& flood_fill_stack,const TV_INT& node);
    void Identify_Colors_Touching_Boundary(const int number_of_colors,const ARRAYS_ND_BASE<TV_INT>& colors,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_x,
        const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_y,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_z,ARRAY<bool>& color_touches_boundary);
    void Identify_Colors_Touching_Color(const int color,const int number_of_colors,const ARRAYS_ND_BASE<TV_INT>& colors,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_x,
        const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_y,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_z,ARRAY<bool>& color_touches_color);
    bool Path_Between_Nodes(const RANGE<TV_INT>& domain,const TV_INT& start_node,const TV_INT& end_node,
        const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_x,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_y,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_z,ARRAY<TV_INT>* path=0);
    void Explore_Path(ARRAYS_ND_BASE<VECTOR<TV_INT,3> >& parents,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_x,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_y,
        const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_z,const TV_INT& node);
//#####################################################################
};
}
#endif

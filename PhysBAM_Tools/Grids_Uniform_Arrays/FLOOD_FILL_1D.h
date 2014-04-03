//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLOOD_FILL_1D
//#####################################################################
#ifndef __FLOOD_FILL_1D__
#define __FLOOD_FILL_1D__

#include <PhysBAM_Tools/Data_Structures/STACK.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
namespace PhysBAM{

class FLOOD_FILL_1D
{
    typedef VECTOR<int,1> TV_INT;
public:
    ARRAY<VECTOR<int,2> > region_boundaries;

    FLOOD_FILL_1D()
    {}

    int Flood_Fill(ARRAYS_ND_BASE<TV_INT>& colors,const ARRAY<bool,FACE_INDEX<1> >& edge_is_blocked,ARRAY<bool>* color_touches_uncolorable_node=0)
    {return Flood_Fill(colors,edge_is_blocked.Component(1),color_touches_uncolorable_node);}

    // colors should be initialized by the user with 0's where colors will be filled and negative values for nodes which will not be colored.
    // -1 is the distinguished uncolorable node (which color_touches_uncolorable_node refers to)
    int Flood_Fill(ARRAYS_ND_BASE<TV_INT>& colors,const ARRAYS_ND_BASE<VECTOR<bool,1> >& edge_is_blocked_x,ARRAY<bool>* color_touches_uncolorable_node=0)
    {int fill_color=1,number_of_regions=0;if(color_touches_uncolorable_node) color_touches_uncolorable_node->Remove_All();region_boundaries.Clean_Memory();
    bool touches_uncolorable_node=false;
    for(int i=colors.domain.min_corner.x;i<=colors.domain.max_corner.x;i++){
        if(colors(i)!=-1){
            colors(i)=fill_color;
            if(number_of_regions!=fill_color){
                if(color_touches_uncolorable_node) touches_uncolorable_node=i>colors.domain.min_corner.x&&colors(i-1)==-1&&!edge_is_blocked_x(i);
                region_boundaries.Append(VECTOR<int,2>(i,i));number_of_regions=fill_color;}
            if(edge_is_blocked_x(i+1)||i==colors.domain.max_corner.x||colors(i+1)==-1){
                region_boundaries(fill_color).y=i;fill_color++;
                if(color_touches_uncolorable_node)
                    color_touches_uncolorable_node->Append(touches_uncolorable_node||(i+1<=colors.domain.max_corner.x&&colors(i+1)==-1&&!edge_is_blocked_x(i+1)));}}}
    return region_boundaries.m;}

    void Identify_Colors_Touching_Boundary(const int number_of_colors,const ARRAYS_ND_BASE<TV_INT>& colors,const ARRAYS_ND_BASE<VECTOR<bool,1> >& edge_is_blocked_x,ARRAY<bool>& color_touches_boundary)
    {color_touches_boundary.Resize(number_of_colors);ARRAYS_COMPUTATIONS::Fill(color_touches_boundary,false);
    int left_color=colors(colors.domain.min_corner.x),right_color=colors(colors.domain.max_corner.x);
    if(left_color>0)color_touches_boundary(left_color)=true;
    if(right_color>0)color_touches_boundary(right_color)=true;}    
//#####################################################################
};
}
#endif

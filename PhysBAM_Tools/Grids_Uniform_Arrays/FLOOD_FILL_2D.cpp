//#####################################################################
// Copyright 2004, Ron Fedkiw, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FLOOD_FILL_2D.h>
using namespace PhysBAM;
FLOOD_FILL_2D::
FLOOD_FILL_2D()
{
    Optimize_Fill_For_Single_Cell_Regions(false);
}
FLOOD_FILL_2D::
~FLOOD_FILL_2D()
{
}
    // colors should be initialized by the user with 0's where colors will be filled and negative values for nodes which will not be colored.
    // -1 is the distinguished uncolorable node (which color_touches_uncolorable_node refers to)
int FLOOD_FILL_2D::
Flood_Fill(ARRAYS_ND_BASE<TV_INT>& colors,const ARRAYS_ND_BASE<VECTOR<bool,2> >& edge_is_blocked_x,const ARRAYS_ND_BASE<VECTOR<bool,2> >& edge_is_blocked_y,ARRAY<bool>* color_touches_uncolorable_node)
{
    TV_INT seed_node;
    int fill_color=0;
    if(color_touches_uncolorable_node) color_touches_uncolorable_node->Remove_All();
    if(optimize_fill_for_single_cell_regions){
        Fill_Single_Cell_Regions(colors,edge_is_blocked_x,edge_is_blocked_y,fill_color);
        if(color_touches_uncolorable_node) color_touches_uncolorable_node->Resize(fill_color);} // Resize sets new elements to false for us
    flood_fill_stack.Preallocate(colors.counts.Product());
    while(Find_Uncolored_Node(colors,seed_node)){
        bool touches_uncolorable_node;fill_color++;Flood_Fill_From_Seed_Node(colors,fill_color,edge_is_blocked_x,edge_is_blocked_y,touches_uncolorable_node,seed_node);
        if(color_touches_uncolorable_node)color_touches_uncolorable_node->Append(touches_uncolorable_node);}
    flood_fill_stack.Clean_Memory();
    return fill_color;
}
void FLOOD_FILL_2D::
Fill_Single_Cell_Regions(ARRAYS_ND_BASE<TV_INT>& colors,const ARRAYS_ND_BASE<VECTOR<bool,2> >& edge_is_blocked_x,const ARRAYS_ND_BASE<VECTOR<bool,2> >& edge_is_blocked_y,int& fill_color)
{
    for(int i=colors.domain.min_corner.x;i<=colors.domain.max_corner.x;i++) for(int j=colors.domain.min_corner.y;j<=colors.domain.max_corner.y;j++) if(colors(i,j)==0)
        if((i==colors.domain.min_corner.x || edge_is_blocked_x(i,j)) && (i==colors.domain.max_corner.x || edge_is_blocked_x(i+1,j)) &&
           (j==colors.domain.min_corner.y || edge_is_blocked_y(i,j)) && (j==colors.domain.max_corner.y || edge_is_blocked_y(i,j+1))) colors(i,j)=++fill_color;
}
bool  FLOOD_FILL_2D::
Find_Uncolored_Node(const ARRAYS_ND_BASE<TV_INT>& colors,TV_INT& node_index)
{
    for(int i=colors.domain.min_corner.x;i<=colors.domain.max_corner.x;i++) for(int j=colors.domain.min_corner.y;j<=colors.domain.max_corner.y;j++)
        if(colors(i,j)==0){node_index=TV_INT(i,j);return true;}
    return false;
}
void FLOOD_FILL_2D::
Flood_Fill_From_Seed_Node(ARRAYS_ND_BASE<TV_INT>& colors,const int fill_color,const ARRAYS_ND_BASE<VECTOR<bool,2> >& edge_is_blocked_x,
    const ARRAYS_ND_BASE<VECTOR<bool,2> >& edge_is_blocked_y,bool& touches_uncolorable_node,const TV_INT& seed_node)
{
    assert(colors(seed_node)==0);
    touches_uncolorable_node=false;
    flood_fill_stack.Remove_All();
    flood_fill_stack.Push(seed_node);
    while(!flood_fill_stack.Empty()){
        TV_INT node=flood_fill_stack.Pop();
        Flood_Fill_Node(colors,fill_color,edge_is_blocked_x,edge_is_blocked_y,touches_uncolorable_node,flood_fill_stack,node);}
}
void FLOOD_FILL_2D::
Flood_Fill_Node(ARRAYS_ND_BASE<TV_INT>& colors,const int fill_color,const ARRAYS_ND_BASE<VECTOR<bool,2> >& edge_is_blocked_x,const ARRAYS_ND_BASE<VECTOR<bool,2> >& edge_is_blocked_y,
    bool& touches_uncolorable_node,STACK<TV_INT >& flood_fill_stack,const TV_INT& node)
{
    if(colors(node)==-1){touches_uncolorable_node=true;return;}else if(colors(node)!=0)return;colors(node)=fill_color;
    if(node.x>colors.domain.min_corner.x&&!edge_is_blocked_x(node.x,node.y)&&colors(node.x-1,node.y)<=0) flood_fill_stack.Push(TV_INT(node.x-1,node.y));
    if(node.x<colors.domain.max_corner.x&&!edge_is_blocked_x(node.x+1,node.y)&&colors(node.x+1,node.y)<=0) flood_fill_stack.Push(TV_INT(node.x+1,node.y));
    if(node.y>colors.domain.min_corner.y&&!edge_is_blocked_y(node.x,node.y)&&colors(node.x,node.y-1)<=0) flood_fill_stack.Push(TV_INT(node.x,node.y-1));
    if(node.y<colors.domain.max_corner.y&&!edge_is_blocked_y(node.x,node.y+1)&&colors(node.x,node.y+1)<=0) flood_fill_stack.Push(TV_INT(node.x,node.y+1));
}
void FLOOD_FILL_2D::
Identify_Colors_Touching_Boundary(const int number_of_colors,const ARRAYS_ND_BASE<TV_INT>& colors,const ARRAYS_ND_BASE<VECTOR<bool,2> >& edge_is_blocked_x,
    const ARRAYS_ND_BASE<VECTOR<bool,2> >& edge_is_blocked_y,ARRAY<bool>& color_touches_boundary)
{
    color_touches_boundary.Resize(number_of_colors);
    ARRAYS_COMPUTATIONS::Fill(color_touches_boundary,false);
    for(int j=colors.domain.min_corner.y;j<=colors.domain.max_corner.y;j++){ // left and right faces
        int left_color=colors(colors.domain.min_corner.x,j),right_color=colors(colors.domain.max_corner.x,j);
        if(left_color>0) color_touches_boundary(left_color)=true;
        if(right_color>0) color_touches_boundary(right_color)=true;}
    for(int i=colors.domain.min_corner.x;i<=colors.domain.max_corner.x;i++){ // bottom and top faces
        int bottom_color=colors(i,colors.domain.min_corner.y),top_color=colors(i,colors.domain.max_corner.y);
        if(bottom_color>0) color_touches_boundary(bottom_color)=true;
        if(top_color>0) color_touches_boundary(top_color)=true;}
}
void FLOOD_FILL_2D::
Identify_Colors_Touching_Color(const int color,const int number_of_colors,const ARRAYS_ND_BASE<TV_INT>& colors,const ARRAYS_ND_BASE<VECTOR<bool,2> >& edge_is_blocked_x,
    const ARRAYS_ND_BASE<VECTOR<bool,2> >& edge_is_blocked_y,ARRAY<bool>& color_touches_color)
{
    color_touches_color.Resize(number_of_colors);
    ARRAYS_COMPUTATIONS::Fill(color_touches_color,false);
    for(int i=colors.domain.min_corner.x;i<=colors.domain.max_corner.x;i++) for(int j=colors.domain.min_corner.y;j<=colors.domain.max_corner.y;j++) if(colors(i,j)==color){
        if(i>colors.domain.min_corner.x&&!edge_is_blocked_x(i,j)&&colors(i-1,j)>0) color_touches_color(colors(i-1,j))=true;
        if(i<colors.domain.max_corner.x&&!edge_is_blocked_x(i+1,j)&&colors(i+1,j)>0) color_touches_color(colors(i+1,j))=true;
        if(j>colors.domain.min_corner.y&&!edge_is_blocked_y(i,j)&&colors(i,j-1)>0) color_touches_color(colors(i,j-1))=true;
        if(j<colors.domain.max_corner.y&&!edge_is_blocked_y(i,j+1)&&colors(i,j+1)>0) color_touches_color(colors(i,j+1))=true;}
}

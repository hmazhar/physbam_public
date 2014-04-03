//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Eran Guendelman, Frank Losasso, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FLOOD_FILL_3D.h>
using namespace PhysBAM;
FLOOD_FILL_3D::
FLOOD_FILL_3D()
{
    Optimize_Fill_For_Single_Cell_Regions(false);
}
FLOOD_FILL_3D::
~FLOOD_FILL_3D()
{
}
int FLOOD_FILL_3D::
Flood_Fill(ARRAYS_ND_BASE<TV_INT>& colors,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_x,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_y,
    const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_z,ARRAY<bool>* color_touches_uncolorable_node)
{
    TV_INT seed_node;
    int fill_color=0;
    if(color_touches_uncolorable_node) color_touches_uncolorable_node->Remove_All();
    if(optimize_fill_for_single_cell_regions){
        Fill_Single_Cell_Regions(colors,edge_is_blocked_x,edge_is_blocked_y,edge_is_blocked_z,fill_color);
        if(color_touches_uncolorable_node) color_touches_uncolorable_node->Resize(fill_color);} // Resize sets new elements to false for us
    flood_fill_stack.Preallocate(colors.counts.Product());
    last_uncolored_node=colors.domain.min_corner;
    while(Find_Uncolored_Node(colors,seed_node)){
        bool touches_uncolorable_node;
        fill_color++;
        Flood_Fill_From_Seed_Node(colors,fill_color,edge_is_blocked_x,edge_is_blocked_y,edge_is_blocked_z,touches_uncolorable_node,seed_node);
        if(color_touches_uncolorable_node)color_touches_uncolorable_node->Append(touches_uncolorable_node);}
    flood_fill_stack.Clean_Memory();
    return fill_color;
}
void FLOOD_FILL_3D::
Fill_Single_Cell_Regions(ARRAYS_ND_BASE<TV_INT>& colors,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_x,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_y,
    const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_z,int& fill_color)
{
    TV_INT i;
    const TV_INT &a(colors.domain.min_corner),&b(colors.domain.max_corner);
    for(i.x=a.x;i.x<=b.x;i.x++) for(i.y=a.y;i.y<=b.y;i.y++) for(i.z=a.z;i.z<=b.z;i.z++) if(colors(i)==0){
        if((i.x==a.x || edge_is_blocked_x(i)) && (i.x==b.x || edge_is_blocked_x(TV_INT(i.x+1,i.y,i.z))) &&
            (i.y==a.y || edge_is_blocked_y(i)) && (i.y==b.y || edge_is_blocked_y(TV_INT(i.x,i.y+1,i.z))) &&
            (i.z==a.z || edge_is_blocked_z(i)) && (i.z==b.z || edge_is_blocked_z(TV_INT(i.x,i.y,i.z+1)))) colors(i)=++fill_color;}
}
bool FLOOD_FILL_3D::
Find_Uncolored_Node(const ARRAYS_ND_BASE<TV_INT>& colors,TV_INT& node_index)
{
    bool first_time_through_loop=true;
    const TV_INT &a(colors.domain.min_corner),&b(colors.domain.max_corner);
    TV_INT i;
    for(i.x=a.x;i.x<=b.x;i.x++) for(i.y=a.y;i.y<=b.y;i.y++) for(i.z=a.z;i.z<=b.z;i.z++){
        if(first_time_through_loop){first_time_through_loop=false;i=last_uncolored_node;}
        if(colors(i)==0){last_uncolored_node=node_index=i;return true;}}
    return false;
}
void FLOOD_FILL_3D::
Flood_Fill_From_Seed_Node(ARRAYS_ND_BASE<TV_INT>& colors,const int fill_color,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_x,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_y,
    const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_z,bool& touches_uncolorable_node,const TV_INT& seed_node)
{
    assert(colors(seed_node)==0);
    touches_uncolorable_node=false;
    flood_fill_stack.Remove_All();
    flood_fill_stack.Push(seed_node);
    while(!flood_fill_stack.Empty()){
        TV_INT node=flood_fill_stack.Pop();
        Flood_Fill_Node(colors,fill_color,edge_is_blocked_x,edge_is_blocked_y,edge_is_blocked_z,touches_uncolorable_node,flood_fill_stack,node);}
}
void FLOOD_FILL_3D::
Flood_Fill_Node(ARRAYS_ND_BASE<TV_INT>& colors,const int fill_color,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_x,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_y,
    const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_z,bool& touches_uncolorable_node,STACK<TV_INT>& flood_fill_stack,const TV_INT& node)
{
    if(colors(node)==-1){touches_uncolorable_node=true;return;}
    else if(colors(node)!=0)return;colors(node)=fill_color;
    if(node.x>colors.domain.min_corner.x &&!edge_is_blocked_x(node.x,node.y,node.z)&&colors(node.x-1,node.y,node.z)<=0) flood_fill_stack.Push(TV_INT(node.x-1,node.y,node.z));
    if(node.x<colors.domain.max_corner.x&&!edge_is_blocked_x(node.x+1,node.y,node.z)&&colors(node.x+1,node.y,node.z)<=0) flood_fill_stack.Push(TV_INT(node.x+1,node.y,node.z));
    if(node.y>colors.domain.min_corner.y&&!edge_is_blocked_y(node.x,node.y,node.z)&&colors(node.x,node.y-1,node.z)<=0) flood_fill_stack.Push(TV_INT(node.x,node.y-1,node.z));
    if(node.y<colors.domain.max_corner.y&&!edge_is_blocked_y(node.x,node.y+1,node.z)&&colors(node.x,node.y+1,node.z)<=0) flood_fill_stack.Push(TV_INT(node.x,node.y+1,node.z));
    if(node.z>colors.domain.min_corner.z&&!edge_is_blocked_z(node.x,node.y,node.z)&&colors(node.x,node.y,node.z-1)<=0) flood_fill_stack.Push(TV_INT(node.x,node.y,node.z-1));
    if(node.z<colors.domain.max_corner.z&&!edge_is_blocked_z(node.x,node.y,node.z+1)&&colors(node.x,node.y,node.z+1)<=0) flood_fill_stack.Push(TV_INT(node.x,node.y,node.z+1));
}
void FLOOD_FILL_3D::
Identify_Colors_Touching_Boundary(const int number_of_colors,const ARRAYS_ND_BASE<TV_INT>& colors,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_x,
    const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_y,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_z,ARRAY<bool>& color_touches_boundary)
{
    color_touches_boundary.Resize(number_of_colors);ARRAYS_COMPUTATIONS::Fill(color_touches_boundary,false);
    for(int j=colors.domain.min_corner.y;j<=colors.domain.max_corner.y;j++) for(int k=colors.domain.min_corner.z;k<=colors.domain.max_corner.z;k++){ // left and right faces
        int left_color=colors(colors.domain.min_corner.x ,j,k),right_color=colors(colors.domain.max_corner.x,j,k);
        if(left_color>0) color_touches_boundary(left_color)=true;
        if(right_color>0) color_touches_boundary(right_color)=true;}
    for(int i=colors.domain.min_corner.x ;i<=colors.domain.max_corner.x;i++) for(int k=colors.domain.min_corner.z;k<=colors.domain.max_corner.z;k++){ // bottom and top faces
        int bottom_color=colors(i,colors.domain.min_corner.y,k),top_color=colors(i,colors.domain.max_corner.y,k);
        if(bottom_color>0) color_touches_boundary(bottom_color)=true;
        if(top_color>0) color_touches_boundary(top_color)=true;}
    for(int i=colors.domain.min_corner.x ;i<=colors.domain.max_corner.x;i++) for(int j=colors.domain.min_corner.y;j<=colors.domain.max_corner.y;j++){ // front and back faces
        int front_color=colors(i,j,colors.domain.min_corner.z),back_color=colors(i,j,colors.domain.max_corner.z);
        if(front_color>0) color_touches_boundary(front_color)=true;
        if(back_color>0) color_touches_boundary(back_color)=true;}
}
void FLOOD_FILL_3D::
Identify_Colors_Touching_Color(const int color,const int number_of_colors,const ARRAYS_ND_BASE<TV_INT>& colors,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_x,
    const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_y,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_z,ARRAY<bool>& color_touches_color)
{
    color_touches_color.Resize(number_of_colors);
    ARRAYS_COMPUTATIONS::Fill(color_touches_color,false);
    for(int i=colors.domain.min_corner.x ;i<=colors.domain.max_corner.x;i++) for(int j=colors.domain.min_corner.y;j<=colors.domain.max_corner.y;j++)
        for(int k=colors.domain.min_corner.z;k<=colors.domain.max_corner.z;k++) if(colors(i,j,k)==color){
            if(i>colors.domain.min_corner.x &&!edge_is_blocked_x(i,j,k)&&colors(i-1,j,k)>0) color_touches_color(colors(i-1,j,k))=true;
            if(i<colors.domain.max_corner.x&&!edge_is_blocked_x(i+1,j,k)&&colors(i+1,j,k)>0) color_touches_color(colors(i+1,j,k))=true;
            if(j>colors.domain.min_corner.y&&!edge_is_blocked_y(i,j,k)&&colors(i,j-1,k)>0) color_touches_color(colors(i,j-1,k))=true;
            if(j<colors.domain.max_corner.y&&!edge_is_blocked_y(i,j+1,k)&&colors(i,j+1,k)>0) color_touches_color(colors(i,j+1,k))=true;
            if(k>colors.domain.min_corner.z&&!edge_is_blocked_z(i,j,k)&&colors(i,j,k-1)>0) color_touches_color(colors(i,j,k-1))=true;
            if(k<colors.domain.max_corner.z&&!edge_is_blocked_z(i,j,k+1)&&colors(i,j,k+1)>0) color_touches_color(colors(i,j,k+1))=true;}
}
bool FLOOD_FILL_3D::
Path_Between_Nodes(const RANGE<TV_INT>& domain,const TV_INT& start_node,const TV_INT& end_node,
    const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_x,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_y,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_z,ARRAY<TV_INT>* path)
{
    ARRAY<TV_INT,VECTOR<int,3> > parents(domain,false);
    parents.Fill(TV_INT(INT_MAX,INT_MAX,INT_MAX));
    flood_fill_stack.Remove_All();
    flood_fill_stack.Preallocate(parents.counts.Product());
    flood_fill_stack.Push(end_node);
    bool success=false;
    parents(end_node)=end_node;
    while(!flood_fill_stack.Empty()){
        TV_INT node=flood_fill_stack.Pop();
        if(node==start_node){success=true;break;}
        Explore_Path(parents,edge_is_blocked_x,edge_is_blocked_y,edge_is_blocked_z,node);}
    if(success){
        if(path){path->Remove_All();path->Preallocate(20);
            TV_INT node=start_node;
            for(;;){path->Append(node);node=parents(node);if(node==end_node){path->Append(node);break;}}}
        return true;}
    else return false;
}
void FLOOD_FILL_3D::
Explore_Path(ARRAYS_ND_BASE<VECTOR<TV_INT,3> >& parents,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_x,const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_y,
    const ARRAYS_ND_BASE<VECTOR<bool,3> >& edge_is_blocked_z,const TV_INT& node)
{
    if(node.x>parents.domain.min_corner.x &&!edge_is_blocked_x(node.x,node.y,node.z)&&parents(node.x-1,node.y,node.z).x==INT_MAX){
        TV_INT new_node(node.x-1,node.y,node.z);parents(new_node)=node;flood_fill_stack.Push(new_node);}
    if(node.x<parents.domain.max_corner.x&&!edge_is_blocked_x(node.x+1,node.y,node.z)&&parents(node.x+1,node.y,node.z).x==INT_MAX){
        TV_INT new_node(node.x+1,node.y,node.z);parents(new_node)=node;flood_fill_stack.Push(new_node);}
    if(node.y>parents.domain.min_corner.y&&!edge_is_blocked_y(node.x,node.y,node.z)&&parents(node.x,node.y-1,node.z).x==INT_MAX){
        TV_INT new_node(node.x,node.y-1,node.z);parents(new_node)=node;flood_fill_stack.Push(new_node);}
    if(node.y<parents.domain.max_corner.y&&!edge_is_blocked_y(node.x,node.y+1,node.z)&&parents(node.x,node.y+1,node.z).x==INT_MAX){
        TV_INT new_node(node.x,node.y+1,node.z);parents(new_node)=node;flood_fill_stack.Push(new_node);}
    if(node.z>parents.domain.min_corner.z&&!edge_is_blocked_z(node.x,node.y,node.z)&&parents(node.x,node.y,node.z-1).x==INT_MAX){
        TV_INT new_node(node.x,node.y,node.z-1);parents(new_node)=node;flood_fill_stack.Push(new_node);}
    if(node.z<parents.domain.max_corner.z&&!edge_is_blocked_z(node.x,node.y,node.z+1)&&parents(node.x,node.y,node.z+1).x==INT_MAX){
        TV_INT new_node(node.x,node.y,node.z+1);parents(new_node)=node;flood_fill_stack.Push(new_node);}
}

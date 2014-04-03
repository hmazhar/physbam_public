//#####################################################################
// Copyright 2002-2009, Kevin Der, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Bridget Vuong.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIRECTED_GRAPH
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Data_Structures/DIRECTED_GRAPH.h>
#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
DIRECTED_GRAPH_CORE::
DIRECTED_GRAPH_CORE(const int number_nodes)
{
    Initialize(number_nodes);
}
//#####################################################################
// Function Reset
//#####################################################################
void DIRECTED_GRAPH_CORE::
Reset()
{
    for(int i(1);i<=parents.m;i++){parents(i).Remove_All();children(i).Remove_All();}
    level_of_node.Clean_Memory();nodes_in_level.Clean_Memory();
}
//#####################################################################
// Function Initialize
//#####################################################################
void DIRECTED_GRAPH_CORE::
Initialize(const int number_nodes)
{
    parents.Resize(number_nodes,false,false);children.Resize(number_nodes,false,false);Reset();
}
//#####################################################################
// Function Resize
//#####################################################################
void DIRECTED_GRAPH_CORE::
Resize(const int number_nodes)
{
    parents.Resize(number_nodes);children.Resize(number_nodes);
}
//#####################################################################
// Function Add_Edge
//#####################################################################
void DIRECTED_GRAPH_CORE::
Add_Edge(const int parent,const int child,bool ensure_distinct)
{
    if(ensure_distinct){parents(child).Append_Unique(parent);children(parent).Append_Unique(child);}
    else{parents(child).Append(parent);children(parent).Append(child);}
}
//#####################################################################
// Function Remove_Edge
//#####################################################################
void DIRECTED_GRAPH_CORE::
Remove_Edge(const int parent,const int child)
{
    children(parent).Remove_Index_Lazy(children(parent).Find(child));parents(child).Remove_Index_Lazy(parents(child).Find(parent));
}
//#####################################################################
// Function Topological_Sort
//#####################################################################
bool DIRECTED_GRAPH_CORE::
Topological_Sort(ARRAY<int>& finish_time,ARRAY<int>& node_index) // returns false if cycles are found
{
    finish_time.Resize(parents.m);node_index.Resize(parents.m);
    ARRAY<short,int> visit_tag(parents.m);int time(1);bool cycle_free=true;
    for(int i(1);i<=parents.m;i++) if(!visit_tag(i)) Visit(i,finish_time,node_index,visit_tag,time,cycle_free);
    return cycle_free;
}
//#####################################################################
// Function Topological_Sort_Assuming_Cycle_Free
//#####################################################################
void DIRECTED_GRAPH_CORE::
Topological_Sort_Assuming_Cycle_Free(ARRAY<int>& finish_time,ARRAY<int>& node_index)
{
    finish_time=CONSTANT_ARRAY<int>(parents.m,0);node_index.Resize(parents.m,false,false);
    int time(1);for(int i(1);i<=parents.m;i++)if(!finish_time(i)) Visit_Assuming_Cycle_Free(i,finish_time,node_index,time);
}
//#####################################################################
// Function Strongly_Connected_Components
//#####################################################################
int DIRECTED_GRAPH_CORE::
Strongly_Connected_Components(ARRAY<int>& component_id) // returns total number of components
{
    component_id.Resize(parents.m);
    ARRAYS_COMPUTATIONS::Fill(component_id,0);ARRAY<int> finish_time;ARRAY<int> node_index;Topological_Sort(finish_time,node_index);
    int total_components=0;ARRAY<int> component;component.Preallocate(Value(parents.Size()));ARRAY<short,int> visit_tag(parents.m); // use parents.m for DIRECTED_GRAPH_CORE<int>
    for(int time=parents.m;time>=1;time--) if(!visit_tag(node_index(time))){
        total_components++;component.Remove_All();Visit_Transpose(node_index(time),visit_tag,component);
        for(int j=1;j<=component.m;j++) component_id(component(j))=total_components;}
    return total_components;
}
//#####################################################################
// Function Generate_Levels
//#####################################################################
void DIRECTED_GRAPH_CORE::
Generate_Levels() // collapses all cycles into a single level
{
    level_of_node.Resize(parents.m);
    ARRAY<int> component_id;int number_components=Strongly_Connected_Components(component_id);
    DIRECTED_GRAPH_CORE component_graph(number_components);
    for(int i(1);i<=parents.m;i++){int parent=component_id(i);for(int j=1;j<=children(i).m;j++){
        int child=component_id(children(i)(j));if(parent != child) component_graph.Add_Edge(parent,child,true);}}
    ARRAY<int> finish_time,node_index;component_graph.Topological_Sort_Assuming_Cycle_Free(finish_time,node_index);
    nodes_in_level.Resize(number_components);for(int i=1;i<=number_components;i++) nodes_in_level(i).Remove_All();
    for(int i(1);i<=parents.m;i++){int level=number_components+1-finish_time(component_id(i));level_of_node(i)=level;nodes_in_level(level).Append(i);} // level 1 has the highest finish time
}
//#####################################################################
// Function Maximal_Depth_On_Acyclic_Graph
//#####################################################################
namespace {
template<class T_DEPTHS> int
Get_Depth_In_Graph(const DIRECTED_GRAPH_CORE& self,const int index,T_DEPTHS& depths PHYSBAM_DEBUG_ONLY(,const int recursion_depth=0))
{
    if(depths(index)) return depths(index);
    assert(int(recursion_depth)<self.parents.m); // detect cycles
    int max_depth=0;for(int k=1;k<=self.parents(index).m;k++) max_depth=max(max_depth,Get_Depth_In_Graph(self,self.parents(index)(k),depths PHYSBAM_DEBUG_ONLY(,recursion_depth+1)));
    return depths(index)=max_depth+1;
}
}
void DIRECTED_GRAPH_CORE::
Maximal_Depth_On_Acyclic_Graph(ARRAY<int>& depths) const // depth is 1-based
{
    depths=CONSTANT_ARRAY<int>(parents.m,0);
    for(int i(1);i<=parents.m;i++) Get_Depth_In_Graph(*this,i,depths);
}
//#####################################################################
// Function Visit
//#####################################################################
void DIRECTED_GRAPH_CORE::
Visit(const int index,ARRAY<int>& finish_time,ARRAY<int>& node_index,ARRAY<short,int>& visit_tag,int& time,bool& cycle_free)
{
    visit_tag(index)=1;
    for(int i=1;i<=children(index).m;i++){
        if(!visit_tag(children(index)(i))) Visit(children(index)(i),finish_time,node_index,visit_tag,time,cycle_free);
        else if(visit_tag(children(index)(i)) == 1) cycle_free=false;}
    visit_tag(index)=2;finish_time(index)=time;node_index(time)=index;time++;
}
//#####################################################################
// Function Visit_Transpose
//#####################################################################
void DIRECTED_GRAPH_CORE::
Visit_Transpose(const int index,ARRAY<short,int>& visit_tag,ARRAY<int>& component)
{
    visit_tag(index)=1;component.Append(index);
    for(int i=1;i<=parents(index).m;i++) if(!visit_tag(parents(index)(i))) Visit_Transpose(parents(index)(i),visit_tag,component);
    visit_tag(index)=2;
}
//#####################################################################
// Function Visit_Assuming_Cycle_Free
//#####################################################################
void DIRECTED_GRAPH_CORE::
Visit_Assuming_Cycle_Free(const int index,ARRAY<int>& finish_time,ARRAY<int>& node_index,int& time)
{
    for(int i=1;i<=children(index).m;i++)if(!finish_time(children(index)(i)))
        Visit_Assuming_Cycle_Free(children(index)(i),finish_time,node_index,time);
    finish_time(index)=time;node_index(time)=index;time++;
}
//#####################################################################
}

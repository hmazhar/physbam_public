//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Grids_Dyadic/READ_WRITE_BINTREE_CELL.h>
#include <PhysBAM_Tools/Read_Write/Grids_Dyadic/READ_WRITE_BINTREE_CHILDREN.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function Read
//#####################################################################
template<class RW,class T> void Read_Write<BINTREE_CELL<T>,RW>::
Read(std::istream& input,BINTREE_CELL<T>& object)
{
    assert(object.Parent() == 0);if(object.children) delete object.children;object.children=0;
    Read_Binary<RW>(input,object.owner->parents_center,object.owner->childrens_DX,object.owner->childrens_depth);
    bool has_nodes;Read_Binary<RW>(input,has_nodes);
    bool has_faces;Read_Binary<RW>(input,has_faces);
    if(has_nodes){
        if(!object.owner->nodes) object.owner->nodes=new int[3];
        for(int i=0;i<2;i++) Read_Binary<RW>(input,object.Node(i));}
    else delete[] object.owner->nodes;
    if(has_faces){
        if(!object.owner->faces) object.owner->faces=new int[3];
        for(int i=0;i<2;i++) Read_Binary<RW>(input,object.Face(i));}
    else delete[] object.owner->faces;
    Read_Helper(input,object);
}
//#####################################################################
// Function Read_Helper
//#####################################################################
template<class RW,class T> void Read_Write<BINTREE_CELL<T>,RW>::
Read_Helper(std::istream& input,BINTREE_CELL<T>& object)
{
    bool has_children;
    Read_Binary<RW>(input,has_children);
    Read_Binary<RW>(input,object.child_index);
    Read_Binary<RW>(input,object.cell);assert(object.children==0);
    if(has_children){
        object.children=new BINTREE_CHILDREN<T>(object.owner->nodes!=0,object.owner->faces!=0);
        object.children->parent=&object;Read_Write<BINTREE_CHILDREN<T>,RW>::Read(input,*object.children);}
    else object.children=0;
}
//#####################################################################
// Function Write
//#####################################################################
template<class RW,class T> void Read_Write<BINTREE_CELL<T>,RW>::
Write(std::ostream& output,const BINTREE_CELL<T>& object)
{
    Write_Binary<RW>(output,object.owner->parents_center,object.owner->childrens_DX,object.owner->childrens_depth);
    bool has_nodes=object.owner->nodes!=0;Write_Binary<RW>(output,has_nodes);
    bool has_faces=object.owner->faces!=0;Write_Binary<RW>(output,has_faces);
    if(has_nodes) for(int i=0;i<2;i++) Write_Binary<RW>(output,object.Node(i));
    if(has_faces) for(int i=0;i<2;i++) Write_Binary<RW>(output,object.Face(i));
    Write_Helper(output,object);
}
//#####################################################################
// Function Write_Helper
//#####################################################################
template<class RW,class T> void Read_Write<BINTREE_CELL<T>,RW>::
Write_Helper(std::ostream& output,const BINTREE_CELL<T>& object)
{
    bool has_children=object.Has_Children();
    Write_Binary<RW>(output,has_children);
    Write_Binary<RW>(output,object.child_index);
    Write_Binary<RW>(output,object.cell);
    if(has_children) Read_Write<BINTREE_CHILDREN<T>,RW>::Write(output,*object.children);
}
template class Read_Write<BINTREE_CELL<float>,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class Read_Write<BINTREE_CELL<float>,double>;
template class Read_Write<BINTREE_CELL<double>,double>;
template class Read_Write<BINTREE_CELL<double>,float>;
#endif
#endif
#endif

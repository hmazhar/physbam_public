#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Grids_Dyadic/READ_WRITE_OCTREE_CELL.h>
#include <PhysBAM_Tools/Read_Write/Grids_Dyadic/READ_WRITE_OCTREE_CHILDREN.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_3D.h>
using namespace PhysBAM;
//#####################################################################
// Function Read
//#####################################################################
template<class RW,class T> void Read_Write<OCTREE_CELL<T>,RW>::
Read(std::istream& input,OCTREE_CELL<T>& object)
{
    assert(object.Parent() == 0);if(object.children) delete object.children;object.children=0;
    Read_Binary<RW>(input,object.owner->parents_center,object.owner->childrens_DX,object.owner->childrens_depth);
    bool has_nodes;Read_Binary<RW>(input,has_nodes);if(has_nodes){
        if(!object.owner->nodes)object.owner->nodes=new int[27];for(int i=0;i<8;i++) Read_Binary<RW>(input,object.Node(i));} else delete[] object.owner->nodes;
    bool has_faces;Read_Binary<RW>(input,has_faces);if(has_faces){
        if(!object.owner->faces)object.owner->faces=new int[36];for(int i=0;i<6;i++) Read_Binary<RW>(input,object.Face(i));} else delete[] object.owner->faces;
    Read_Helper(input,object);
}
//#####################################################################
// Function Read_Helper
//#####################################################################
template<class RW,class T> void Read_Write<OCTREE_CELL<T>,RW>::
Read_Helper(std::istream& input,OCTREE_CELL<T>& object)
{
    bool has_children;Read_Binary<RW>(input,has_children);
    Read_Binary<RW>(input,object.child_index);Read_Binary<RW>(input,object.cell);assert(object.children==0);
    if(has_children){object.children=new OCTREE_CHILDREN<T>(object.owner->nodes!=0,object.owner->faces!=0);object.children->parent=&object;Read_Binary<RW>(input,*object.children);}
    else object.children=0;
}
//#####################################################################
// Function Write
//#####################################################################
template<class RW,class T> void Read_Write<OCTREE_CELL<T>,RW>::
Write(std::ostream& output,const OCTREE_CELL<T>& object)
{
    Write_Binary<RW>(output,object.owner->parents_center,object.owner->childrens_DX,object.owner->childrens_depth);
    bool has_nodes=object.owner->nodes!=0;Write_Binary<RW>(output,has_nodes);if(has_nodes) for(int i=0;i<8;i++) Write_Binary<RW>(output,object.Node(i));
    bool has_faces=object.owner->faces!=0;Write_Binary<RW>(output,has_faces);if(has_faces) for(int i=0;i<6;i++) Write_Binary<RW>(output,object.Face(i));
    Write_Helper(output,object);
}
//#####################################################################
// Function Write_Helper
//#####################################################################
template<class RW,class T> void Read_Write<OCTREE_CELL<T>,RW>::
Write_Helper(std::ostream& output,const OCTREE_CELL<T>& object)
{
    bool has_children=object.Has_Children();Write_Binary<RW>(output,has_children);
    Write_Binary<RW>(output,object.child_index);Write_Binary<RW>(output,object.cell);
    if(has_children) Write_Binary<RW>(output,*object.children);
}
template class Read_Write<OCTREE_CELL<float>,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class Read_Write<OCTREE_CELL<float>,double>;
template class Read_Write<OCTREE_CELL<double>,double>;
template class Read_Write<OCTREE_CELL<double>,float>;
#endif
#endif
#endif

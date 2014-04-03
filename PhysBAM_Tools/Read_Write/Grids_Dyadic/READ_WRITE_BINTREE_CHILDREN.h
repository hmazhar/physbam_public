//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_BINTREE_CHILDREN
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_BINTREE_CHILDREN__
#define __READ_WRITE_BINTREE_CHILDREN__

#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_CHILDREN.h>
#include <PhysBAM_Tools/Read_Write/Grids_Dyadic/READ_WRITE_BINTREE_CELL.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<BINTREE_CHILDREN<T>,RW>
{
public:
    static void Read(std::istream& input,BINTREE_CHILDREN<T>& object)
    {Read_Binary<RW>(input,object.parents_center,object.childrens_DX,object.childrens_depth);
    bool has_nodes=false;Read_Binary<RW>(input,has_nodes); if(has_nodes) for(int i=0;i<3;i++) Read_Binary<RW>(input,object.nodes[i]);
    bool has_faces=false;Read_Binary<RW>(input,has_faces); if(has_faces) for(int i=0;i<3;i++) Read_Binary<RW>(input,object.faces[i]);
    for(int i=0;i<2;i++){
        int dummy=-1;object.children[i].Initialize(&object,dummy,i);Read_Write<BINTREE_CELL<T>,RW>::Read_Helper(input,object.children[i]);}}

    static void Write(std::ostream& output,const BINTREE_CHILDREN<T>& object)
    {Write_Binary<RW>(output,object.parents_center,object.childrens_DX,object.childrens_depth);
    bool has_nodes=object.nodes!=0;Write_Binary<RW>(output,has_nodes); if(has_nodes) for(int i=0;i<3;i++) Write_Binary<RW>(output,object.nodes[i]);
    bool has_faces=object.faces!=0;Write_Binary<RW>(output,has_faces); if(has_faces) for(int i=0;i<3;i++) Write_Binary<RW>(output,object.faces[i]);
    for(int i=0;i<2;i++) Read_Write<BINTREE_CELL<T>,RW>::Write_Helper(output,object.children[i]);}
};
}
#endif
#endif
#endif

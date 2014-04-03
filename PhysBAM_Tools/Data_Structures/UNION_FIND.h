//#####################################################################
// Copyright 2003-2007, Don Hatch, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNION_FIND
//#####################################################################
#ifndef __UNION_FIND__
#define __UNION_FIND__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/CONSTANT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class ID> // ID=int
class UNION_FIND
{
    typedef unsigned char T_RANK;
public:
    mutable ARRAY<ID,ID> parents;
    ARRAY<T_RANK,ID> ranks; 

    explicit UNION_FIND(const ID entries=ID())
        :parents(entries),ranks(entries)
    {}

    void Initialize(const ID entries)
    {parents=CONSTANT_ARRAY<ID,ID>(entries,ID());ranks=CONSTANT_ARRAY<T_RANK,ID>(entries,0);}

    ID Size() const
    {return parents.Size();}

    void Clear_Connectivity()
    {ARRAYS_COMPUTATIONS::Fill(parents,ID());ARRAYS_COMPUTATIONS::Fill(ranks,0);}

    ID Add_Entry()
    {parents.Append(ID());ranks.Append(0);return Size();}

    bool Is_Root(const ID i) const
    {return !parents(i);}

    ID Find(const ID i) const
    {ID root=Find_Without_Path_Compression(i);Path_Compress(i,root);return root;}

    ID Union(const ID i,const ID j)
    {ID root_i=Find_Without_Path_Compression(i),root_j=Find_Without_Path_Compression(j);
    ID root=ranks(root_i)>=ranks(root_j)?root_i:root_j;
    Path_Compress(i,root);Path_Compress(j,root);
    if(ranks(root_i)==ranks(root_j) && root_i!=root_j) ranks(root)++;
    return root;}

    template<int d>
    ID Union(const VECTOR<ID,d>& indices)
    {ID root=Find_Without_Path_Compression(indices[1]);T_RANK max_rank=ranks(root);bool max_tie=false;
    for(int i=2;i<=d;i++){
        ID root_i=Find_Without_Path_Compression(indices[i]);
        if(max_rank<ranks(root_i)){max_rank=ranks(root_i);root=root_i;max_tie=false;}
        else if(max_rank==ranks(root_i) && root_i!=root) max_tie=true;}
    for(int i=1;i<=d;i++) Path_Compress(indices[i],root);
    if(max_tie) ranks(root)++;return root;}

    void Merge(const UNION_FIND<ID>& union_find)
    {assert(Size()==union_find.Size());
    for(ID i(1);i<=Size();i++){ID j=union_find.parents(i);if(j) Union(i,j);}}

    // Okay for map to yield invalid indices for isolated elements
    template<class ID2,class T_ARRAY>
    void Mapped_Merge(const UNION_FIND<ID2>& union_find,const T_ARRAY& map)
    {for(ID2 i(1);i<=union_find.Size();i++){ID2 root=union_find.Find(i);if(i!=root) Union(map(i),map(root));}}

    void Forest_Edges(ARRAY<PAIR<ID,ID> >& pairs) const
    {pairs.Remove_All();for(ID i(1);i<=Size();i++){ID j=Find(i);if(i!=j) pairs.Append(PAIR<ID,ID>(i,j));}}

    void Merge_Forest_Edges(const ARRAY<PAIR<ID,ID> >& pairs)
    {for(int i=1;i<=pairs.m;i++) Union(pairs(i).x,pairs(i).y);}

private:
    ID Find_Without_Path_Compression(const ID i) const
    {ID j=i;while(parents(j)) j=parents(j);return j;}

    void Path_Compress(const ID i,const ID root) const
    {ID j=i;while(j!=root){ID parent=parents(j);parents(j)=root;j=parent;if(!j) break;}}

//#####################################################################
};
}
#endif

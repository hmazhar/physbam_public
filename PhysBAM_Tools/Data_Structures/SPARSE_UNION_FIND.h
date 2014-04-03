//#####################################################################
// Copyright 2008, Don Hatch, Geoffrey Irving, Michael Lentine, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPARSE_UNION_FIND
//#####################################################################
#ifndef __SPARSE_UNION_FIND__
#define __SPARSE_UNION_FIND__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/CONSTANT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class ID> // ID=int
class SPARSE_UNION_FIND
{
    typedef unsigned char T_RANK;
public:
    mutable HASHTABLE<ID,ID> parents;
    HASHTABLE<ID,T_RANK> ranks;
    ID size;

    explicit SPARSE_UNION_FIND(const ID entries=ID()) 
        :size(entries)
    {}

    void Initialize(const ID entries)
    {Clear_Connectivity();size=entries;}

    ID Size() const
    {return size;}

    void Clear_Connectivity()
    {parents.Remove_All();ranks.Remove_All();}

    ID Add_Entry()
    {return ++size;}

    bool Is_Root(const ID i) const
    {return !parents.Contains(i);}

    ID Find(const ID i) const
    {ID root=Find_Without_Path_Compression(i);Path_Compress(i,root);return root;}

    ID Union(const ID i,const ID j)
    {ID root_i=Find_Without_Path_Compression(i),root_j=Find_Without_Path_Compression(j);
    T_RANK rank_i(0),rank_j(0);ranks.Get(root_i,rank_i);ranks.Get(root_j,rank_j);
    ID root=rank_i>=rank_j?root_i:root_j;
    Path_Compress(i,root);Path_Compress(j,root);
    if(rank_i==rank_j && root_i!=root_j) ranks.Get_Or_Insert(root)++;return root;}

    template<class T_ARRAY>
    int Union(const T_ARRAY& array)
    {int root=0;typename T_ARRAY::ELEMENT i(1);for(;i<=array.Size();i++){root=Find(array(i));break;}if(!root) return 0;
    for(;i<=array.Size();i++) Union(root,array(i));return Find(root);}

    template<int d>
    ID Union(const VECTOR<ID,d>& indices)
    {T_RANK max_rank(0);ID root=Find_Without_Path_Compression(indices[1]);bool max_tie=false;ranks.Get(root,max_rank);
    for(int i=2;i<=d;i++) {
        ID root_i=Find_Without_Path_Compression(indices[i]);
        T_RANK tmp_rank(0);ranks.Get(root_i,tmp_rank);
        if(max_rank<tmp_rank){max_rank=tmp_rank;root=root_i;max_tie=false;}
        else if(max_rank==tmp_rank && root!=root_i) max_tie=true;}
    for(int i=1;i<=d;i++) Path_Compress(indices[i],root);
    if(max_tie) ranks.Get_Or_Insert(root)++;return root;}

    void Merge(const SPARSE_UNION_FIND<ID>& union_find)
    {assert(Size()==union_find.Size());
    for(HASHTABLE_ITERATOR<ID,ID> iterator(union_find.parents);iterator.Valid();iterator.Next()) Union(iterator.Key(),iterator.Data());}

    // Okay for map to yield invalid indices for isolated elements
    template<class ID2,class T_ARRAY>
    void Mapped_Merge(const SPARSE_UNION_FIND<ID2>& union_find,const T_ARRAY& map)
    {for(HASHTABLE_ITERATOR<ID2,ID2> iterator(union_find.parents);iterator.Valid();iterator.Next()) Union(map(iterator.Key()),map(iterator.Data()));}

    void Forest_Edges(ARRAY<PAIR<ID,ID> >& pairs) const
    {pairs.Remove_All();for(HASHTABLE_ITERATOR<ID,ID> iterator(parents);iterator.Valid();iterator.Next()) pairs.Append(PAIR<ID,ID>(iterator.Key(),iterator.Data()));}

    void Merge_Forest_Edges(const ARRAY<PAIR<ID,ID> >& pairs)
    {for(int i=1;i<=pairs.m;i++) Union(pairs(i).x,pairs(i).y);}

private:
    ID Find_Without_Path_Compression(const ID i) const
    {ID k,j=i;while(parents.Get(j,k)) j=k;return j;}

    void Path_Compress(const ID i,const ID root) const
    {ID j=i;while(j && j!=root){ID &ref_parent=parents.Get_Or_Insert(j),parent=ref_parent;ref_parent=root;j=parent;}}

//#####################################################################
};
}
#endif

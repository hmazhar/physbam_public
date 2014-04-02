//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Fracture/CUTTING_SIMPLEX.h>
using namespace PhysBAM;
//#####################################################################
// Function Hash_Reduce
//#####################################################################
template<class T,int d> int
Hash_Reduce(const CUTTING_SIMPLEX<T,d>& s)
{
    return HASH(HASH(s.nodes,s.type,s.weights,s.parent,s.element_owner),HASH(s.element_original_coordinates,s.simplex_original_coordinates)).value;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> CUTTING_SIMPLEX<T,d>::
CUTTING_SIMPLEX()
    :parent(0),element_owner(0)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> CUTTING_SIMPLEX<T,d>::
CUTTING_SIMPLEX(const VECTOR<int,d>& nodes,SIMPLEX_TYPE type)
    : type(type),parent(0),element_owner(0),nodes(nodes)
{
    assert(type==GLOBAL_EMBEDDING_FACE||type==GLOBAL_CUT_FACE);
}
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> CUTTING_SIMPLEX<T,d>::
CUTTING_SIMPLEX(const VECTOR<int,d>& nodes,const SIMPLEX_TYPE type,const VECTOR<VECTOR<T,d>,d>& weights,const int parent,const int element_owner)
    :type(type),parent(parent),element_owner(element_owner),nodes(nodes),weights(weights),abs_tol(0)
{
    assert(type==LOCAL_EMBEDDING_FACE);
    for(int i=1;i<=d;++i){assert(nodes[i]!=0);for(int j=1;j<=d;++j) assert(weights[i][j]==0||weights[i][j]==1);}
}
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> CUTTING_SIMPLEX<T,d>::
CUTTING_SIMPLEX(const VECTOR<int,d>& nodes,SIMPLEX_TYPE type,T abs_tol,int parent,int element_owner,const VECTOR<VECTOR<T,d>,d+1>& element_original_coordinates,
    const VECTOR<VECTOR<T,d>,d>& simplex_original_coordinates)
    :type(type),parent(parent),element_owner(element_owner),nodes(nodes),abs_tol(abs_tol),element_original_coordinates(element_original_coordinates),
    simplex_original_coordinates(simplex_original_coordinates)
{
    assert(type==LOCAL_CUT_FACE);
    for(int i=1;i<d;++i) assert(nodes[i]!=0);
    if(nodes[d]==0){
        this->simplex_original_coordinates[d].Fill(0);
        for(int i=1;i<d;++i) this->simplex_original_coordinates[d]+=simplex_original_coordinates[i];
        this->simplex_original_coordinates[d]/=(d-1);}
    // compute weights for each (real) node
    if(nodes[d]!=0){
        VECTOR<VECTOR<T,d+1>,d> temp_weights;
        Barycentric_Coordinates<EXACT_FLOAT<T> >(element_original_coordinates,simplex_original_coordinates,temp_weights,abs_tol);
        for(int i=1;i<=d;++i) for(int j=1;j<=d;++j) weights[i][j]=temp_weights[i][j];}
    else{
        VECTOR<VECTOR<T,d>,d-1> simplex2;
        for(int i=1;i<d;++i) simplex2[i]=simplex_original_coordinates[i];
        VECTOR<VECTOR<T,d+1>,d-1> temp_weights;
        Barycentric_Coordinates<EXACT_FLOAT<T> >(element_original_coordinates,simplex2,temp_weights,abs_tol);
        // set the fake node's weights to be at the center of the real nodes
        weights[d].Fill(0);
        for(int i=1;i<d;++i){
            for(int j=1;j<=d;++j) weights[i][j]=temp_weights[i][j];
            weights[d]+=weights[i];}
        weights[d]/=(T)(d-1);}
    // determine whether each node is inside the embedding element
    VECTOR<VECTOR<GET_ADAPTIVE_WEIGHTS_RESULT_TYPE,d>,d> adaptive_weights;
    Get_Adaptive_Weights(adaptive_weights);
    for(int i=1;i<=d;++i){
        bool& inside=node_in_embedded_simplex[i];
        inside=true;
        for(int j=1;j<=d&&inside;++j) {
            int sign=adaptive_weights[i][j].Sign();
            if(sign==0) PHYSBAM_FATAL_ERROR("degeneracy");
            inside=(sign>0);}
        if(inside){
            int sign=((T)1-Adaptive_Vector_Sum<void>(adaptive_weights[i])).Sign();
            if(sign==0) PHYSBAM_FATAL_ERROR("degeneracy");
            inside=(sign>0);}}
}
//#####################################################################
// Function Get_Adaptive_Weights
//#####################################################################
template<class T,int d> void CUTTING_SIMPLEX<T,d>::
Get_Adaptive_Weights(VECTOR<VECTOR<GET_ADAPTIVE_WEIGHTS_RESULT_TYPE,d>,d>& adaptive_weights) const
{
    VECTOR<int,d> node_indices;
    for(int i=1;i<=d;++i) node_indices[i]=i;
    Get_Adaptive_Weights(adaptive_weights,node_indices);
}
//#####################################################################
// Function Get_Adaptive_Weights
//#####################################################################
template<class T,int d> template<int N> void CUTTING_SIMPLEX<T,d>::
Get_Adaptive_Weights(VECTOR<VECTOR<GET_ADAPTIVE_WEIGHTS_RESULT_TYPE,d>,N>& adaptive_weights,const VECTOR<int,N>& node_indices) const
{
    STATIC_ASSERT((N<=d));
    assert(type==LOCAL_EMBEDDING_FACE||type==LOCAL_CUT_FACE);
    // by non-const reference to allow updating of weights
    CUTTING_SIMPLEX_COORDINATE_SHARED<T,d>* shared=new CUTTING_SIMPLEX_COORDINATE_SHARED<T,d>(*const_cast<CUTTING_SIMPLEX<T,d>*>(this));shared->Add_Ref();
    for(int i=1;i<=N;++i) Get_Adaptive_Weights(adaptive_weights[i],node_indices[i],shared);
    assert((shared->Ref_Count()-1)%d==0);
    if(shared->Release()==0) delete shared;
}
//#####################################################################
// Function Get_Adaptive_Weights
//#####################################################################
template<class T,int d> void CUTTING_SIMPLEX<T,d>::
Get_Adaptive_Weights(VECTOR<GET_ADAPTIVE_WEIGHTS_RESULT_TYPE,d>& adaptive_weights,int node_index) const
{
    assert(type==LOCAL_EMBEDDING_FACE||type==LOCAL_CUT_FACE);
    // by non-const reference to allow updating of weights
    CUTTING_SIMPLEX_COORDINATE_SHARED<T,d>* shared=new CUTTING_SIMPLEX_COORDINATE_SHARED<T,d>(*const_cast<CUTTING_SIMPLEX<T,d>*>(this));shared->Add_Ref();
    Get_Adaptive_Weights(adaptive_weights,node_index,shared);
    assert(shared->Ref_Count()==1||shared->Ref_Count()==1+d);
    if(shared->Release()==0) delete shared;
}
//#####################################################################
// Function Get_Adaptive_Weights
//#####################################################################
template<class T,int d> void CUTTING_SIMPLEX<T,d>::
Get_Adaptive_Weights(VECTOR<GET_ADAPTIVE_WEIGHTS_RESULT_TYPE,d>& adaptive_weights,int node_index,CUTTING_SIMPLEX_COORDINATE_SHARED<T,d>* shared) const
{
    assert(1<=node_index&&node_index<=d);
    // for a fake node, we shouldn't be using the coordinate, so no need to resort to exact arithmetic
    T error=(nodes[node_index]!=0?abs_tol:0);
    for(int j=1;j<=d;++j){
        T local_weight=weights[node_index][j],local_error=(error>0?max(error,ulp(local_weight)/2):0);
        adaptive_weights[j].Set(node_index,j,local_weight,local_error,shared);}
}
//#####################################################################
template class CUTTING_SIMPLEX<float,2>;
template class CUTTING_SIMPLEX<float,3>;
template void CUTTING_SIMPLEX<float,3>::Get_Adaptive_Weights<2>(VECTOR<VECTOR<CUTTING_SIMPLEX_COORDINATE<float,3>,3>,2>&,VECTOR<int,2> const&) const;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CUTTING_SIMPLEX<double,2>;
template class CUTTING_SIMPLEX<double,3>;
template void CUTTING_SIMPLEX<double,3>::Get_Adaptive_Weights<2>(VECTOR<VECTOR<CUTTING_SIMPLEX_COORDINATE<double,3>,3>,2>&,VECTOR<int,2> const&) const;
#endif

//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __CUTTING_SIMPLEX_COORDINATE__
#define __CUTTING_SIMPLEX_COORDINATE__

#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_BASE.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_OBJECT.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/EXACT_FLOAT.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/EXACT_RATIONAL.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Adaptive_Geometry/Barycentric_Coordinates.h>

#include <cassert>

namespace PhysBAM{
template<class T,int d>
class CUTTING_SIMPLEX_COORDINATE;

namespace ADAPTIVE_DETAIL
{
template<class T,int d> struct IS_ADAPTIVE<PhysBAM::CUTTING_SIMPLEX_COORDINATE<T,d> > {static const bool value=true;};
}

template<class T,int d> struct CUTTING_SIMPLEX;
template<class T,int d> class CUTTING_SIMPLEX_COORDINATE_SHARED;

template<class T,int d>
class CUTTING_SIMPLEX_COORDINATE : public ADAPTIVE_BASE< CUTTING_SIMPLEX_COORDINATE<T,d>, EXACT_RATIONAL<T> >
{
    int node_index;
    int coordinate_index;
    T value;
    T abs_tol;
    CUTTING_SIMPLEX_COORDINATE_SHARED<T,d>* shared;

public:
    CUTTING_SIMPLEX_COORDINATE()
        :abs_tol(0),shared(0)
    {}

    CUTTING_SIMPLEX_COORDINATE(T value_input)
        :value(value_input),abs_tol(0)
    {}

    CUTTING_SIMPLEX_COORDINATE(int node_index_input,int coordinate_index_input,T value_input,T abs_tol_input,const CUTTING_SIMPLEX_COORDINATE_SHARED<T,d>* shared_input)
        :node_index(node_index_input),coordinate_index(coordinate_index_input),value(value_input),abs_tol(abs_tol_input),shared(shared_input)
    {
        assert(abs_tol==0||shared!=0);
        if(shared) shared->Add_Ref();
    }

    ~CUTTING_SIMPLEX_COORDINATE()
    {if(shared && shared->Release()==0) delete shared;}

    void Set(T value_input)
    {value=value_input;abs_tol=0;shared=0;}

    void Set(int node_index_input,int coordinate_index_input,T value_input,T abs_tol_input,CUTTING_SIMPLEX_COORDINATE_SHARED<T,d>* shared_input)
    {if(abs_tol_input==0){Set(value_input);return;}
    assert(abs_tol_input>0);
    assert(shared_input!=0);
    node_index=node_index_input;
    coordinate_index=coordinate_index_input;
    value=value_input;
    abs_tol=abs_tol_input;
    if(shared && shared->Release()==0) delete shared;
    shared=shared_input;shared->Add_Ref();}

    PAIR<T,T> Estimate_And_Error_Implementation() const
    {return PAIR<T,T>(value,abs_tol);}

    EXACT_RATIONAL<T> Exact_Implementation() const
    {return (abs_tol==0?EXACT_RATIONAL<T>(value):shared->Exact(node_index,coordinate_index));}
//#####################################################################
};



template<class T,int d>
class CUTTING_SIMPLEX_COORDINATE_SHARED : NONCOPYABLE
{
    typedef VECTOR<VECTOR<ADAPTIVE_OBJECT< EXACT_RATIONAL<T> >,d+1>,d> COORDINATES_COLLECTION_TYPE;

    int reference_count;
    CUTTING_SIMPLEX<T,d>& cutting_simplex;
    COORDINATES_COLLECTION_TYPE* adaptive_coordinates;
public:

    CUTTING_SIMPLEX_COORDINATE_SHARED(CUTTING_SIMPLEX<T,d>& cutting_simplex_input)
        : reference_count(0),cutting_simplex(cutting_simplex_input),adaptive_coordinates(0) {}
    ~CUTTING_SIMPLEX_COORDINATE_SHARED() {assert(reference_count==0);delete adaptive_coordinates;}

    int Ref_Count()
    {return reference_count;}

    int Add_Ref()
    {assert(reference_count>=0);return ++reference_count;}

    int Release()
    {assert(reference_count>0);return --reference_count;}

    EXACT_RATIONAL<T> Exact(int node_index, int coordinate_index)
    {assert(1<=node_index&&node_index<=d);
    assert(1<=coordinate_index&&coordinate_index<=d+1);
    assert(cutting_simplex.nodes[node_index]!=0);
    if(adaptive_coordinates==0){
        adaptive_coordinates=new COORDINATES_COLLECTION_TYPE;
        Barycentric_Coordinates< EXACT_FLOAT<T> >(cutting_simplex.element_original_coordinates,cutting_simplex.simplex_original_coordinates,*adaptive_coordinates);}
        EXACT_RATIONAL<T> exact=(*adaptive_coordinates)[node_index][coordinate_index].Exact();
        if(coordinate_index<=d) cutting_simplex.weights[node_index][coordinate_index]=exact.Compress_And_Estimate();
        return exact;}
//#####################################################################
};
}
#endif

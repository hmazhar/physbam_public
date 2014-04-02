//#####################################################################
// Copyright 2002-2006, Ron Fedkiw, Sergey Koltakov, Robert Bridson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUADRILATERAL_MESH
//#####################################################################
#ifndef __QUADRILATERAL_MESH__
#define __QUADRILATERAL_MESH__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{ 

class QUADRILATERAL_MESH
{
public:
    int number_nodes;
    ARRAY<VECTOR<int,4> > quadrilaterals; // array of 4 indices for each quadrilateral - quadrilaterals(i)(j) is j'th index in quadrilateral i

    QUADRILATERAL_MESH()
        :number_nodes(0)
    {}

    QUADRILATERAL_MESH(const int number_nodes_input,const ARRAY<VECTOR<int,4> >& quadrilateral_list)
    {
        Initialize_Quadrilateral_Mesh(number_nodes_input,quadrilateral_list);
    }

    ~QUADRILATERAL_MESH()
    {}

private:
    void operator=(const QUADRILATERAL_MESH&); // use Initialize_Mesh instead
public:

    void Clean_Memory()
    {quadrilaterals.Clean_Memory();}

    void Initialize_Quadrilateral_Mesh(const int number_nodes_input,const ARRAY<VECTOR<int,4> >& quadrilateral_list)
    {Clean_Memory();number_nodes=number_nodes_input;quadrilaterals=quadrilateral_list;}

    void Initialize_Square_Mesh(const int m,const int n)
    {Clean_Memory();
    number_nodes=m*n;quadrilaterals.Resize((m-1)*(n-1));
    int t=0;for(int j=1;j<=n-1;j++) for(int i=1;i<=m-1;i++){t++;quadrilaterals(t)(1)=i+m*(j-1);quadrilaterals(t)(2)=(i+1)+m*(j-1);quadrilaterals(t)(3)=i+m*j;quadrilaterals(t)(4)=(i+1)+m*j;}}

//#####################################################################
};   
}
#endif

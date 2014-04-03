//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_GRID_TRANSFER_ITERATOR
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_GRID_TRANSFER_ITERATOR__
#define __RLE_GRID_TRANSFER_ITERATOR__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_HORIZONTAL.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_Y.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_LINEAR_PROFILE.h>
namespace PhysBAM{

template<class T_GRID,class T_ITERATOR>
class RLE_GRID_TRANSFER_ITERATOR
{
public:
    typedef typename T_GRID::BOX_HORIZONTAL_INT T_BOX_HORIZONTAL_INT;

    T_ITERATOR source,destination;

    RLE_GRID_TRANSFER_ITERATOR(const T_GRID& source_grid,const T_GRID& destination_grid,const int number_of_ghost_cells)
        :source(source_grid,number_of_ghost_cells),destination(destination_grid,number_of_ghost_cells)
    {}

    RLE_GRID_TRANSFER_ITERATOR(const T_GRID& source_grid,const T_BOX_HORIZONTAL_INT& source_region,const T_GRID& destination_grid,const T_BOX_HORIZONTAL_INT& destination_region)
        :source(source_grid,source_region),destination(destination_grid,destination_region)
    {}

    operator bool() const
    {assert((bool)source==(bool)destination);return source;}

    void operator++(int)
    {if(source.jmax()==destination.jmax()){source++;destination++;}
    else if(source.jmax()<destination.jmax()) source++;else destination++;}
//#####################################################################
};

//#####################################################################
// Function Transfer_Cells
//##################5###################################################
template<class T_GRID,class T2> void Transfer_Cells(const RLE_GRID_TRANSFER_ITERATOR<T_GRID,typename T_GRID::CELL_ITERATOR>& cells,const ARRAY<T2>& source,ARRAY<T2>& destination)
{
    assert(cells.source.grid.long_run_cells==2);
    int c1=cells.source.Cell(),c2=cells.destination.Cell();
    destination(c2)=source(c1);if(cells.destination.Long()) destination(c2+1)=source(c1);
}
//#####################################################################
// Function Transfer_Constant_Horizontal_Faces
//#####################################################################
template<class T_GRID,class T_FACE,class T2> void Transfer_Constant_Horizontal_Faces(const RLE_GRID_TRANSFER_ITERATOR<T_GRID,T_FACE>& faces,const ARRAY<T2>& source,ARRAY<T2>& destination)
{
    assert(faces.source.cell1.grid.long_run_faces_horizontal==1);
    int f1=faces.source.Face(),f2=faces.destination.Face();
    destination(f2)+=source(f1)*(min(faces.source.jmax(),faces.destination.jmax())-max(faces.source.j(),faces.destination.j()))/faces.destination.Length();
}
//#####################################################################
// Function Transfer_Linear_Horizontal_Faces
//#####################################################################
template<class T,class T_GRID,class T_FACE> void Transfer_Linear_Horizontal_Faces(const RLE_GRID_TRANSFER_ITERATOR<T_GRID,T_FACE>& faces,const ARRAY<T>& source,ARRAY<T>& destination)
{
    assert(T_FACE::Axis()!=2 && faces.source.cell1.grid.long_run_faces_horizontal==2);
    int f1=faces.source.Face(),f2=faces.destination.Face();
    int jmin=max(faces.source.j(),faces.destination.j()),jmax=min(faces.source.jmax(),faces.destination.jmax()),len_f=jmax-jmin-1;
    RLE_LINEAR_PROFILE<T> velocity_profile;
    if(faces.source.Long()){T V_slope=(source(f1+1)-source(f1))/(faces.source.Length()-1);
        velocity_profile.Add_Linear_Profile(jmin-faces.destination.j(),len_f,source(f1)+(jmin-faces.source.j())*V_slope,V_slope);}
    else velocity_profile.Add_Single_Value(jmin-faces.destination.j(),source(f1));
    velocity_profile.Update_Linear_Profile(0,faces.destination.Length()-1,1,destination(f2),destination(f2+1));
}
//#####################################################################
// Function Transfer_Vertical_Faces
//#####################################################################
template<class T_GRID,class T2> void Transfer_Vertical_Faces(const RLE_GRID_TRANSFER_ITERATOR<T_GRID,RLE_GRID_ITERATOR_FACE_Y<T_GRID> >& faces,const ARRAY<T2>& source,ARRAY<T2>& destination)
{
    typedef typename T_GRID::SCALAR T;
    assert(faces.source.cell1.grid.long_run_cells==2);
    int f1=faces.source.Face(),f2=faces.destination.Face(),j1=faces.source.j(),j2=faces.destination.j();
    if(j1==j2) destination(f2)=source(f1);
    else if(j1<j2) destination(f2)=source(f1+1);
    if(faces.destination.Long()){
        T one_over_len2=(T)1/(faces.destination.cell2.length-1);
        if(j1>j2) destination(f2+1)+=one_over_len2*source(f1);
        destination(f2+1)+=one_over_len2*max(0,min(faces.source.cell2.jmax(),faces.destination.cell2.jmax())-max(j1,j2)-1)*source(f1+1);}
}
//#####################################################################
// Function Transfer
//#####################################################################
template<class T_GRID,class T2> void Transfer(const RLE_GRID_TRANSFER_ITERATOR<T_GRID,typename T_GRID::CELL_ITERATOR>& cells,const ARRAY<T2>& source,ARRAY<T2>& destination)
{Transfer_Cells(cells,source,destination);}
template<class T,class T_GRID,int axis> void Transfer(const RLE_GRID_TRANSFER_ITERATOR<T_GRID,RLE_GRID_ITERATOR_FACE_HORIZONTAL<T_GRID,axis> >& faces,const ARRAY<T>& source,ARRAY<T>& destination)
{Transfer_Linear_Horizontal_Faces(faces,source,destination);}
template<class T_GRID,class T2> void Transfer(const RLE_GRID_TRANSFER_ITERATOR<T_GRID,RLE_GRID_ITERATOR_FACE_Y<T_GRID> >& faces,const ARRAY<T2>& source,ARRAY<T2>& destination)
{Transfer_Vertical_Faces(faces,source,destination);}
//#####################################################################
}
#endif
#endif

//#####################################################################
// Copyright 2005, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FAST_MARCHING
//#####################################################################
#ifndef __FAST_MARCHING__
#define __FAST_MARCHING__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
namespace PhysBAM{

template<class T>
class FAST_MARCHING
{
public:

    template<class T2,class T3,class T4>
    static void Up_Heap(const T2& phi,T3& close_k,ARRAY<T4>& heap,int index)
    {for(;;){int parent=index/2;
        if(parent >= 1 && abs(phi(heap(index))) < abs(phi(heap(parent)))){
            exchange(heap(index),heap(parent));
            close_k(heap(index))=index; // update child k
            close_k(heap(parent))=parent; // update parent k
            index=parent;} // move up one
        else return;}}

    template<class T2,class T3,class T4>
    static void Down_Heap(const T2& phi,T3& close_k,ARRAY<T4>& heap,const int heap_length)
    {int index=1;
    for(;;){int left=2*index,right=2*index+1;
        if(right > heap_length) break;
        else if(abs(phi(heap(left))) <= abs(phi(heap(right)))){
            heap(index)=heap(left); // update i, j, ij
            close_k(heap(index))=index; // update k
            index=left;}  // move down one
        else{
            heap(index)=heap(right); // update i, j, ij
            close_k(heap(index))=index; // update k
            index=right;}} // move down one
    // fill the hole with the last element
    if(index != heap_length){
        heap(index)=heap(heap_length); // update i, j, ij
        close_k(heap(index))=index; // update k
        Up_Heap(phi,close_k,heap,index);}}

    static T Solve_Quadratic(const T phi,const T value_x,const T value_y,const T dx,const T dy)
    {assert(LEVELSET_UTILITIES<T>::Sign(value_x)==LEVELSET_UTILITIES<T>::Sign(value_y));
    if(abs(value_x) >= abs(value_y)+dy) return value_y+LEVELSET_UTILITIES<T>::Sign(phi)*dy;
    else if(abs(value_y) >= abs(value_x)+dx) return value_x+LEVELSET_UTILITIES<T>::Sign(phi)*dx;
    else{T dx2=sqr(dx),dy2=sqr(dy);return (dy2*value_x+dx2*value_y+LEVELSET_UTILITIES<T>::Sign(phi)*dx*dy*sqrt(dx2+dy2-sqr(value_x-value_y)))/(dx2+dy2);}}

    static T Solve_Close_Point(const T phi,const bool no_x,const T value_x,const bool no_y,const T value_y,const T dx,const T dy)
    {T sign=LEVELSET_UTILITIES<T>::Sign(phi);assert(!no_x || !no_y);
    if(no_x) return value_y+sign*dy;
    if(no_y) return value_x+sign*dx;
    return Solve_Quadratic(phi,value_x,value_y,dx,dy);} // candidates exist in both directions

    static T Solve_Close_Point(const T phi,const bool no_x,const T value_x,const bool no_y,const T value_y,const bool no_z,const T value_z,const T dx,const T dy,const T dz)
    {T sign=LEVELSET_UTILITIES<T>::Sign(phi);assert(!no_x || !no_y || !no_z);
    if(no_x){
        if(no_y) return value_z+sign*dz;
        if(no_z) return value_y+sign*dy;
        return Solve_Quadratic(phi,value_y,value_z,dy,dz);}
    if(no_y){
        if(no_z) return value_x+sign*dx;
        return Solve_Quadratic(phi,value_x,value_z,dx,dz);}
    if(no_z) return Solve_Quadratic(phi,value_x,value_y,dx,dy);
    // candidates exist in all three directions
    T value_yz=Solve_Quadratic(phi,value_y,value_z,dy,dz);
    if(abs(value_x) >= abs(value_yz)) return value_yz;
    T value_xz=Solve_Quadratic(phi,value_x,value_z,dx,dz);
    if(abs(value_y) >= abs(value_xz)) return value_xz;
    T value_xy=Solve_Quadratic(phi,value_x,value_y,dx,dy);
    if(abs(value_z) >= abs(value_xy)) return value_xy;
    // use the candidates in all three directions
    T dx2=sqr(dx),dy2=sqr(dy),dz2=sqr(dz),dx2dy2=dx2*dy2,dx2dz2=dx2*dz2,dy2dz2=dy2*dz2;
    return (dy2dz2*value_x+dx2dz2*value_y+dx2dy2*value_z+sign*dx*dy*dz*
        sqrt(dx2dy2+dx2dz2+dy2dz2-dx2*sqr(value_y-value_z)-dy2*sqr(value_x-value_z)-dz2*sqr(value_x-value_y)))/(dx2dy2+dx2dz2+dy2dz2);}

    template<int dimension>
    static T Solve_Close_Point(const T phi,const int number_of_axis,const T value[dimension],const T dx[dimension])
    {assert(number_of_axis);
    if(dimension==1 || number_of_axis==1) return value[0]+LEVELSET_UTILITIES<T>::Sign(phi)*dx[0];
    if(dimension==2 || number_of_axis==2) return Solve_Quadratic(phi,value[0],value[1],dx[0],dx[1]);
    assert(dimension==3); // candidates exist in all three directions (must be in 3d)
    T value_yz=Solve_Quadratic(phi,value[1],value[2],dx[1],dx[2]);
    if(abs(value[0]) >= abs(value_yz)) return value_yz;
    T value_xz=Solve_Quadratic(phi,value[0],value[2],dx[0],dx[2]);
    if(abs(value[1]) >= abs(value_xz)) return value_xz;
    T value_xy=Solve_Quadratic(phi,value[0],value[1],dx[0],dx[1]);
    if(abs(value[2]) >= abs(value_xy)) return value_xy;
    // use the candidates in all three directions
    T dx2=sqr(dx[0]),dy2=sqr(dx[1]),dz2=sqr(dx[2]),dx2dy2=dx2*dy2,dx2dz2=dx2*dz2,dy2dz2=dy2*dz2;
    return (dy2dz2*value[0]+dx2dz2*value[1]+dx2dy2*value[2]+LEVELSET_UTILITIES<T>::Sign(phi)*dx[0]*dx[1]*dx[2]*
        sqrt(dx2dy2+dx2dz2+dy2dz2-dx2*sqr(value[1]-value[2])-dy2*sqr(value[0]-value[2])-dz2*sqr(value[0]-value[1])))/(dx2dy2+dx2dz2+dy2dz2);}

//#####################################################################
};
}
#endif

//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_UTILITIES
//#####################################################################
#include <PhysBAM_Tools/Polynomials/CUBIC.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/POLYGON.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Function Theta_Quadratic
//#####################################################################
template<class T> T LEVELSET_UTILITIES<T>::
Theta_Quadratic(const T phi_left_left,const T phi_left,const T phi_right,const T dx)
{
    T one_over_dx=1/dx;
    QUADRATIC<T> quadratic((T)-.5*(2*phi_left-phi_right-phi_left_left)*sqr(one_over_dx),(T)-.5*(-phi_right+phi_left_left)*one_over_dx,phi_left);
    quadratic.Compute_Roots_In_Interval(0,dx);assert(quadratic.roots==1);return quadratic.root1*one_over_dx;
}
//#####################################################################
// Function Theta_Cubic
//#####################################################################
template<class T> T LEVELSET_UTILITIES<T>::
Theta_Cubic(const T phi_left_left,const T phi_left,const T phi_right,const T phi_right_right,const T dx)
{
    T one_over_dx=1/dx,one_over_dx_squared=sqr(one_over_dx);
    CUBIC<T> cubic((T)one_sixth*(3*(phi_left-phi_right)-phi_left_left+phi_right_right)*one_over_dx_squared*one_over_dx,(T).5*(-2*phi_left+phi_right+phi_left_left)*one_over_dx_squared,
        ((T)-.5*phi_left+phi_right-(T)one_third*phi_left_left-(T)one_sixth*phi_right_right)*one_over_dx,phi_left);
    cubic.Compute_Roots_In_Interval(0,dx);assert(cubic.roots==1);return cubic.root1*one_over_dx;
}
//#####################################################################
// Function Negative_Cell_Fraction
//#####################################################################
// finds the fraction of the 2D cell that has phi <= 0 - aspect ratio is height/width, i.e. y/
template<class T> T LEVELSET_UTILITIES<T>::
Negative_Cell_Fraction(const T phi_lower_left,const T phi_lower_right,const T phi_upper_left,const T phi_upper_right,const T aspect_ratio)
{
    // all corners same sign
    if(phi_lower_left <= 0 && phi_lower_right <= 0 && phi_upper_left <= 0 && phi_upper_right <= 0) return 1; // all negative
    else if(phi_lower_left > 0 && phi_lower_right > 0 && phi_upper_left > 0 && phi_upper_right > 0) return 0; // all positive
    else{ // more work required
        T dy=aspect_ratio,dxy=sqrt(1+sqr(dy)),area=dy; // note that dx=1
        // one corner negative
        if(phi_lower_left <= 0 && phi_lower_right > 0 && phi_upper_left > 0 && phi_upper_right > 0){
            T x_bottom=Theta(phi_lower_left,phi_lower_right);
            T diagonal=Theta(phi_lower_left,phi_upper_right)*dxy,x_diagonal=diagonal/dxy,y_diagonal=diagonal*dy/dxy;
            T y_left=Theta(phi_lower_left,phi_upper_left)*dy;
            POLYGON<VECTOR<T,2> > quadrilateral(4);
            quadrilateral.X(1)=VECTOR<T,2>(0,0);quadrilateral.X(2)=VECTOR<T,2>(x_bottom,0);
            quadrilateral.X(3)=VECTOR<T,2>(x_diagonal,y_diagonal);quadrilateral.X(4)=VECTOR<T,2>(0,y_left);
            return quadrilateral.Area()/area;}
        else if(phi_lower_left > 0 && phi_lower_right <= 0 && phi_upper_left > 0 && phi_upper_right > 0){
            T y_right=Theta(phi_lower_right,phi_upper_right)*dy;
            T diagonal=Theta(phi_upper_left,phi_lower_right)*dxy,x_diagonal=diagonal/dxy,y_diagonal=dy-diagonal*dy/dxy;
            T x_bottom=Theta(phi_lower_left,phi_lower_right);
            POLYGON<VECTOR<T,2> > quadrilateral(4);
            quadrilateral.X(1)=VECTOR<T,2>(1,0);quadrilateral.X(2)=VECTOR<T,2>(1,y_right);
            quadrilateral.X(3)=VECTOR<T,2>(x_diagonal,y_diagonal);quadrilateral.X(4)=VECTOR<T,2>(x_bottom,0);
            return quadrilateral.Area()/area;}
        else if(phi_lower_left > 0 && phi_lower_right > 0 && phi_upper_left <= 0 && phi_upper_right > 0){
            T y_left=Theta(phi_lower_left,phi_upper_left)*dy;
            T diagonal=Theta(phi_upper_left,phi_lower_right)*dxy,x_diagonal=diagonal/dxy,y_diagonal=dy-diagonal*dy/dxy;
            T x_top=Theta(phi_upper_left,phi_upper_right);
            POLYGON<VECTOR<T,2> > quadrilateral(4);
            quadrilateral.X(1)=VECTOR<T,2>(0,dy);quadrilateral.X(2)=VECTOR<T,2>(0,y_left);
            quadrilateral.X(3)=VECTOR<T,2>(x_diagonal,y_diagonal);quadrilateral.X(4)=VECTOR<T,2>(x_top,dy);
            return quadrilateral.Area()/area;}
        else if(phi_lower_left > 0 && phi_lower_right > 0 && phi_upper_left > 0 && phi_upper_right <= 0){
            T x_top=Theta(phi_upper_left,phi_upper_right);
            T diagonal=Theta(phi_lower_left,phi_upper_right)*dxy,x_diagonal=diagonal/dxy,y_diagonal=diagonal*dy/dxy;
            T y_right=Theta(phi_lower_right,phi_upper_right)*dy;
            POLYGON<VECTOR<T,2> > quadrilateral(4);
            quadrilateral.X(1)=VECTOR<T,2>(1,dy);quadrilateral.X(2)=VECTOR<T,2>(x_top,dy);
            quadrilateral.X(3)=VECTOR<T,2>(x_diagonal,y_diagonal);quadrilateral.X(4)=VECTOR<T,2>(1,y_right);
            return quadrilateral.Area()/area;}
        // two corners negative - opposite corners
        else if(phi_lower_left <= 0 && phi_lower_right > 0 && phi_upper_left > 0 && phi_upper_right <= 0){
            T x_bottom=Theta(phi_lower_left,phi_lower_right);
            T x_top=Theta(phi_upper_left,phi_upper_right);
            T y_left=Theta(phi_lower_left,phi_upper_left)*dy;
            T y_right=Theta(phi_lower_right,phi_upper_right)*dy;
            T x_center=(x_bottom+x_top+1)/4,y_center=(y_left+y_right+dy)/4; // average of the four points
            POLYGON<VECTOR<T,2> > quadrilateral(4);
            quadrilateral.X(1)=VECTOR<T,2>(0,0);quadrilateral.X(2)=VECTOR<T,2>(x_bottom,0);
            quadrilateral.X(3)=VECTOR<T,2>(x_center,y_center);quadrilateral.X(4)=VECTOR<T,2>(0,y_left);
            T area1=quadrilateral.Area();
            quadrilateral.X(1)=VECTOR<T,2>(1,dy);quadrilateral.X(2)=VECTOR<T,2>(x_top,dy);
            quadrilateral.X(3)=VECTOR<T,2>(x_center,y_center);quadrilateral.X(4)=VECTOR<T,2>(1,y_right);
            return (area1+quadrilateral.Area())/area;}
        else if(phi_lower_left > 0 && phi_lower_right <= 0 && phi_upper_left <= 0 && phi_upper_right > 0){
            T x_bottom=Theta(phi_lower_left,phi_lower_right);
            T x_top=Theta(phi_upper_left,phi_upper_right);
            T y_left=Theta(phi_lower_left,phi_upper_left)*dy;
            T y_right=Theta(phi_lower_right,phi_upper_right)*dy;
            T x_center=(x_bottom+x_top+1)/4,y_center=(y_left+y_right+dy)/4; // average of the four points
            POLYGON<VECTOR<T,2> > quadrilateral(4);
            quadrilateral.X(1)=VECTOR<T,2>(1,0);quadrilateral.X(2)=VECTOR<T,2>(1,y_right);
            quadrilateral.X(3)=VECTOR<T,2>(x_center,y_center);quadrilateral.X(4)=VECTOR<T,2>(x_bottom,0);
            T area1=quadrilateral.Area();
            quadrilateral.X(1)=VECTOR<T,2>(0,dy);quadrilateral.X(2)=VECTOR<T,2>(0,y_left);
            quadrilateral.X(3)=VECTOR<T,2>(x_center,y_center);quadrilateral.X(4)=VECTOR<T,2>(x_top,dy);
            return (area1+quadrilateral.Area())/area;}
        // two corners negative - adjacent corners
        else if(phi_lower_left <= 0 && phi_lower_right <= 0 && phi_upper_left > 0 && phi_upper_right > 0){
            T y_left=Theta(phi_lower_left,phi_upper_left)*dy;
            T y_right=Theta(phi_lower_right,phi_upper_right)*dy;
            T diagonal_1=Theta(phi_lower_left,phi_upper_right)*dxy,x_diagonal_1=diagonal_1/dxy,y_diagonal_1=diagonal_1*dy/dxy;
            T diagonal_2=Theta(phi_upper_left,phi_lower_right)*dxy,x_diagonal_2=diagonal_2/dxy,y_diagonal_2=dy-diagonal_2*dy/dxy;
            T x_diagonal=(x_diagonal_1+x_diagonal_2)/2,y_diagonal=(y_diagonal_1+y_diagonal_2)/2;
            POLYGON<VECTOR<T,2> > pentagon(5);
            pentagon.X(1)=VECTOR<T,2>(0,0);pentagon.X(2)=VECTOR<T,2>(1,0);
            pentagon.X(3)=VECTOR<T,2>(1,y_right);pentagon.X(4)=VECTOR<T,2>(x_diagonal,y_diagonal);
            pentagon.X(5)=VECTOR<T,2>(0,y_left);
            return pentagon.Area()/area;}
        else if(phi_lower_left <= 0 && phi_lower_right > 0 && phi_upper_left <= 0 && phi_upper_right > 0){
            T x_bottom=Theta(phi_lower_left,phi_lower_right);
            T x_top=Theta(phi_upper_left,phi_upper_right);
            T diagonal_1=Theta(phi_lower_left,phi_upper_right)*dxy,x_diagonal_1=diagonal_1/dxy,y_diagonal_1=diagonal_1*dy/dxy;
            T diagonal_2=Theta(phi_upper_left,phi_lower_right)*dxy,x_diagonal_2=diagonal_2/dxy,y_diagonal_2=dy-diagonal_2*dy/dxy;
            T x_diagonal=(x_diagonal_1+x_diagonal_2)/2,y_diagonal=(y_diagonal_1+y_diagonal_2)/2;
            POLYGON<VECTOR<T,2> > pentagon(5);
            pentagon.X(1)=VECTOR<T,2>(0,0);pentagon.X(2)=VECTOR<T,2>(x_bottom,0);
            pentagon.X(3)=VECTOR<T,2>(x_diagonal,y_diagonal);pentagon.X(4)=VECTOR<T,2>(x_top,dy);
            pentagon.X(5)=VECTOR<T,2>(0,dy);
            return pentagon.Area()/area;}
        else if(phi_lower_left > 0 && phi_lower_right <= 0 && phi_upper_left > 0 && phi_upper_right <= 0){
            T x_bottom=Theta(phi_lower_left,phi_lower_right);
            T x_top=Theta(phi_upper_left,phi_upper_right);
            T diagonal_1=Theta(phi_lower_left,phi_upper_right)*dxy,x_diagonal_1=diagonal_1/dxy,y_diagonal_1=diagonal_1*dy/dxy;
            T diagonal_2=Theta(phi_upper_left,phi_lower_right)*dxy,x_diagonal_2=diagonal_2/dxy,y_diagonal_2=dy-diagonal_2*dy/dxy;
            T x_diagonal=(x_diagonal_1+x_diagonal_2)/2,y_diagonal=(y_diagonal_1+y_diagonal_2)/2;
            POLYGON<VECTOR<T,2> > pentagon(5);
            pentagon.X(1)=VECTOR<T,2>(1,0);pentagon.X(2)=VECTOR<T,2>(1,dy);
            pentagon.X(3)=VECTOR<T,2>(x_top,dy);pentagon.X(4)=VECTOR<T,2>(x_diagonal,y_diagonal);
            pentagon.X(5)=VECTOR<T,2>(x_bottom,0);
            return pentagon.Area()/area;}
        else if(phi_lower_left > 0 && phi_lower_right > 0 && phi_upper_left <= 0 && phi_upper_right <= 0){
            T y_left=Theta(phi_lower_left,phi_upper_left)*dy;
            T y_right=Theta(phi_lower_right,phi_upper_right)*dy;
            T diagonal_1=Theta(phi_lower_left,phi_upper_right)*dxy,x_diagonal_1=diagonal_1/dxy,y_diagonal_1=diagonal_1*dy/dxy;
            T diagonal_2=Theta(phi_upper_left,phi_lower_right)*dxy,x_diagonal_2=diagonal_2/dxy,y_diagonal_2=dy-diagonal_2*dy/dxy;
            T x_diagonal=(x_diagonal_1+x_diagonal_2)/2,y_diagonal=(y_diagonal_1+y_diagonal_2)/2;
            POLYGON<VECTOR<T,2> > pentagon(5);
            pentagon.X(1)=VECTOR<T,2>(0,dy);pentagon.X(2)=VECTOR<T,2>(0,y_left);
            pentagon.X(3)=VECTOR<T,2>(x_diagonal,y_diagonal);pentagon.X(4)=VECTOR<T,2>(1,y_right);
            pentagon.X(5)=VECTOR<T,2>(1,dy);
            return pentagon.Area()/area;}
        // three corners negative
        else if(phi_lower_left > 0 && phi_lower_right <= 0 && phi_upper_left <= 0 && phi_upper_right <= 0){
            T x_bottom=Theta(phi_lower_left,phi_lower_right);
            T diagonal=Theta(phi_lower_left,phi_upper_right)*dxy,x_diagonal=diagonal/dxy,y_diagonal=diagonal*dy/dxy;
            T y_left=Theta(phi_lower_left,phi_upper_left)*dy;
            POLYGON<VECTOR<T,2> > quadrilateral(4);
            quadrilateral.X(1)=VECTOR<T,2>(0,0);quadrilateral.X(2)=VECTOR<T,2>(x_bottom,0);
            quadrilateral.X(3)=VECTOR<T,2>(x_diagonal,y_diagonal);quadrilateral.X(4)=VECTOR<T,2>(0,y_left);
            return (dy-quadrilateral.Area())/area;}
        else if(phi_lower_left <= 0 && phi_lower_right > 0 && phi_upper_left <= 0 && phi_upper_right <= 0){
            T y_right=Theta(phi_lower_right,phi_upper_right)*dy;
            T diagonal=Theta(phi_upper_left,phi_lower_right)*dxy,x_diagonal=diagonal/dxy,y_diagonal=dy-diagonal*dy/dxy;
            T x_bottom=Theta(phi_lower_left,phi_lower_right);
            POLYGON<VECTOR<T,2> > quadrilateral(4);
            quadrilateral.X(1)=VECTOR<T,2>(1,0);quadrilateral.X(2)=VECTOR<T,2>(1,y_right);
            quadrilateral.X(3)=VECTOR<T,2>(x_diagonal,y_diagonal);quadrilateral.X(4)=VECTOR<T,2>(x_bottom,0);
            return (dy-quadrilateral.Area())/area;}
        else if(phi_lower_left <= 0 && phi_lower_right <= 0 && phi_upper_left > 0 && phi_upper_right <= 0){
            T y_left=Theta(phi_lower_left,phi_upper_left)*dy;
            T diagonal=Theta(phi_upper_left,phi_lower_right)*dxy,x_diagonal=diagonal/dxy,y_diagonal=dy-diagonal*dy/dxy;
            T x_top=Theta(phi_upper_left,phi_upper_right);
            POLYGON<VECTOR<T,2> > quadrilateral(4);
            quadrilateral.X(1)=VECTOR<T,2>(0,dy);quadrilateral.X(2)=VECTOR<T,2>(0,y_left);
            quadrilateral.X(3)=VECTOR<T,2>(x_diagonal,y_diagonal);quadrilateral.X(4)=VECTOR<T,2>(x_top,dy);
            return (dy-quadrilateral.Area())/area;}
        else if(phi_lower_left <= 0 && phi_lower_right <= 0 && phi_upper_left <= 0 && phi_upper_right > 0){
            T x_top=Theta(phi_upper_left,phi_upper_right);
            T diagonal=Theta(phi_lower_left,phi_upper_right)*dxy,x_diagonal=diagonal/dxy,y_diagonal=diagonal*dy/dxy;
            T y_right=Theta(phi_lower_right,phi_upper_right)*dy;
            POLYGON<VECTOR<T,2> > quadrilateral(4);
            quadrilateral.X(1)=VECTOR<T,2>(1,dy);quadrilateral.X(2)=VECTOR<T,2>(x_top,dy);
            quadrilateral.X(3)=VECTOR<T,2>(x_diagonal,y_diagonal);quadrilateral.X(4)=VECTOR<T,2>(1,y_right);
            return (dy-quadrilateral.Area())/area;}}

    return 0; // should never be called
}
//#####################################################################
template class LEVELSET_UTILITIES<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_UTILITIES<double>;
#endif

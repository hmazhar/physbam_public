//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_FACE_VELOCITY_DONOR_CELL_RLE
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_HORIZONTAL.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_Y.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_LINEAR_PROFILE.h>
#include <PhysBAM_Tools/Grids_RLE_Advection/ADVECTION_FACE_VELOCITY_DONOR_CELL_RLE.h>
#include <PhysBAM_Tools/Log/LOG.h>
using namespace PhysBAM;
//#####################################################################
// Function Apply_Horizontal_Velocity_Flux
//#####################################################################
template<class T_FACE_H,class T,class T_GRID> static void
Apply_Horizontal_Velocity_Flux(const T_GRID& grid,ARRAY<T>& Z,const ARRAY<T>& Z_ghost,const ARRAY<T>& V,const T dt,ARRAY<T>& divergence)
{
    // TODO: this uses twice as many cell iterators as necessary, speed up
    assert(grid.long_run_cells==2);
    T dt_over_length=dt/grid.uniform_grid.dX.x;
    T_FACE_H face1(grid,0,-T_FACE_H::Sentinels()),face2(grid,0,T_FACE_H::Sentinels());
    if(grid.long_run_faces_horizontal==1) PHYSBAM_FATAL_ERROR();
    else{
        while(face1){
            int jmin=max(face1.j(),face2.j()),jmax=min(face1.jmax(),face2.jmax()),len_f=jmax-jmin-1,dj1=face1.j()-jmin,dj2=face2.j()-jmin;
            int len_f1=face1.Length()-1,len_f2=face2.Length()-1;
            int f1=face1.Face(),f2=face2.Face();

            if(!(len_f1 || len_f2)){ // handle all short case specially
                T Vh_mid=(T).5*(V(f1)+V(f2));
                int donor=Vh_mid>0?f1:f2;
                T flux_volume=dt_over_length*V(donor);
                T flux=flux_volume*Z_ghost(donor);
                Z(f1)-=flux;Z(f2)+=flux;
                divergence(f1)-=flux_volume;divergence(f2)+=flux_volume;}
            else{
                T Vh1_lo=V(f1),Vh2_lo=V(f2),Vh1_slope=0,Vh2_slope=0;
                if(len_f1){Vh1_slope=(V(f1+1)-V(f1))/len_f1;Vh1_lo+=(jmin-face1.j())*Vh1_slope;}
                if(len_f2){Vh2_slope=(V(f2+1)-V(f2))/len_f2;Vh2_lo+=(jmin-face2.j())*Vh2_slope;}
                T Vh_mid_lo=(T).5*(Vh1_lo+Vh2_lo),Vh_mid_slope=(T).5*(Vh1_slope+Vh2_slope); // find velocity profile at flux

                RLE_LINEAR_PROFILE<T> velocity_profile;
                RLE_LINEAR_PROFILE<T> divergence_profile;

                if(Vh_mid_lo*(Vh_mid_lo+len_f*Vh_mid_slope)<0){ // have zero crossing
                    assert(Vh_mid_slope);
                    // jmin to jmin+k goes one direction, and jmin+k+1 to jmax-1 goes the other direction
                    int k=min((int)(-Vh_mid_lo/Vh_mid_slope),len_f-1);assert(0<=k&&k<len_f);
                    T Zh1_slope=len_f1?(Z_ghost(f1+1)-Z_ghost(f1))/len_f1:0,Zh2_slope=len_f2?(Z_ghost(f2+1)-Z_ghost(f2))/len_f2:0;
                    if(Vh_mid_lo>0){
                        velocity_profile.Add_Linear_Product_Profile(0,k,Vh1_lo,Vh1_slope,Z_ghost(f1)+(jmin-face1.j())*Zh1_slope,Zh1_slope);
                        divergence_profile.Add_Linear_Profile(0,k,Vh1_lo,Vh1_slope);
                        velocity_profile.Add_Linear_Product_Profile(k+1,len_f-k-1,Vh2_lo+(k+1)*Vh2_slope,Vh2_slope,Z_ghost(f2)+(jmin+k+1-face2.j())*Zh2_slope,Zh2_slope);
                        divergence_profile.Add_Linear_Profile(k+1,len_f-k-1,Vh2_lo+(k+1)*Vh2_slope,Vh2_slope);}
                    else{
                        velocity_profile.Add_Linear_Product_Profile(0,k,Vh2_lo,Vh2_slope,Z_ghost(f2)+(jmin-face2.j())*Zh2_slope,Zh2_slope);
                        divergence_profile.Add_Linear_Profile(0,k,Vh2_lo,Vh2_slope);
                        velocity_profile.Add_Linear_Product_Profile(k+1,len_f-k-1,Vh1_lo+(k+1)*Vh1_slope,Vh1_slope,Z_ghost(f1)+(jmin+k+1-face1.j())*Zh1_slope,Zh1_slope);
                        divergence_profile.Add_Linear_Profile(k+1,len_f-k-1,Vh1_lo+(k+1)*Vh1_slope,Vh1_slope);}}
                else{ // no zero crossing
                    if(Vh_mid_lo>0){
                        T Zh1_slope=len_f1?(Z_ghost(f1+1)-Z_ghost(f1))/len_f1:0;
                        velocity_profile.Add_Linear_Product_Profile(0,len_f,Vh1_lo,Vh1_slope,Z_ghost(f1)+(jmin-face1.j())*Zh1_slope,Zh1_slope);
                        divergence_profile.Add_Linear_Profile(0,len_f,Vh1_lo,Vh1_slope);}
                    else{
                        T Zh2_slope=len_f2?(Z_ghost(f2+1)-Z_ghost(f2))/len_f2:0;
                        velocity_profile.Add_Linear_Product_Profile(0,len_f,Vh2_lo,Vh2_slope,Z_ghost(f2)+(jmin-face2.j())*Zh2_slope,Zh2_slope);
                        divergence_profile.Add_Linear_Profile(0,len_f,Vh2_lo,Vh2_slope);}}

                velocity_profile.Update_Linear_Profile(dj1,len_f1,-dt_over_length,Z(f1),Z(f1+1));
                divergence_profile.Update_Linear_Profile(dj1,len_f1,-dt_over_length,divergence(f1),divergence(f1+1));
                velocity_profile.Update_Linear_Profile(dj2,len_f2,+dt_over_length,Z(f2),Z(f2+1));
                divergence_profile.Update_Linear_Profile(dj2,len_f2,+dt_over_length,divergence(f2),divergence(f2+1));}

            if(face1.jmax()==face2.jmax()){face1++;face2++;}else if(face1.jmax()<face2.jmax()) face1++;else face2++;}}
}
//#####################################################################
// Function Apply_Horizontal_Horizontal_Velocity_Flux
//#####################################################################
template<class T_FACE_H1,class T_FACE_H2,class T,class T_GRID> static void
Apply_Horizontal_Horizontal_Velocity_Flux(const T_GRID& grid,ARRAY<T>& Z,const ARRAY<T>& Z_ghost,const ARRAY<T>& V,const T dt,ARRAY<T>& divergence)
{
    // TODO: this uses twice as many cell iterators as necessary, speed up
    assert(grid.long_run_cells==2);
    T dt_over_length=dt/grid.uniform_grid.dX.x;
    T_FACE_H1 face11(grid,0,-T_FACE_H2::Sentinels()),face12(grid,0,T_FACE_H2::Sentinels());
    T_FACE_H2 face21(grid,0,-T_FACE_H1::Sentinels()),face22(grid,0,T_FACE_H1::Sentinels());
    if(grid.long_run_faces_horizontal==1){
        while(face11){
            int f11=face11.Face(),f12=face12.Face(),f21=face21.Face(),f22=face22.Face();
            T Vh1_mid=(T).5*(V(f11)+V(f12));
            T Zh2_donor=Vh1_mid>0?Z_ghost(f21):Z_ghost(f22);
            T flux_volume=dt_over_length*Vh1_mid*(min(face11.jmax(),face12.jmax())-max(face11.j(),face12.j()));
            T flux=flux_volume*Zh2_donor;
            Z(f21)-=flux/face21.Length();Z(f22)+=flux/face22.Length();
            divergence(f21)-=flux_volume/face21.Length();divergence(f22)+=flux_volume/face22.Length();
            if(face11.jmax()==face12.jmax()){face11++;face12++;}else if(face11.jmax()<face12.jmax()) face11++;else face12++;
            if(face21.jmax()==face22.jmax()){face21++;face22++;}else if(face21.jmax()<face22.jmax()) face21++;else face22++;}}
    else{
        while(face11){
            int jmin=max(face11.j(),face12.j()),jmax=min(face11.jmax(),face12.jmax()),len_f=jmax-jmin-1,dj1=face21.j()-jmin,dj2=face22.j()-jmin;
            int len_f11=face11.Length()-1,len_f12=face12.Length()-1,len_f21=face21.Length()-1,len_f22=face22.Length()-1;
            int f11=face11.Face(),f12=face12.Face(),f21=face21.Face(),f22=face22.Face();

            if(!(len_f11 || len_f12 || len_f21 || len_f22)){ // handle all short case specially
                T Vh1_mid=(T).5*(V(f11)+V(f12));
                T Zh2_donor=Vh1_mid>0?Z_ghost(f21):Z_ghost(f22);
                T flux_volume=dt_over_length*Vh1_mid;
                T flux=flux_volume*Zh2_donor;
                Z(f21)-=flux;Z(f22)+=flux;
                divergence(f21)-=flux_volume;divergence(f22)+=flux_volume;}
            else{
                T Vh1_lo=V(f11)+V(f22),Vh1_slope=0; // find velocity profile at flux
                if(len_f11){T slope=(V(f11+1)-V(f11))/len_f11;Vh1_lo+=(jmin-face11.j())*slope;Vh1_slope+=slope;}
                if(len_f12){T slope=(V(f12+1)-V(f12))/len_f12;Vh1_lo+=(jmin-face12.j())*slope;Vh1_slope+=slope;}
                Vh1_lo*=(T).5;Vh1_slope*=(T).5;

                RLE_LINEAR_PROFILE<T> velocity_profile;

                if(Vh1_lo*(Vh1_lo+len_f*Vh1_slope)<0){ // have zero crossing
                    assert(Vh1_slope);
                    // jmin to jmin+k goes one direction, and jmin+k+1 to jmax-1 goes the other direction
                    int k=min((int)(-Vh1_lo/Vh1_slope),len_f-1);assert(0<=k&&k<len_f);
                    T Zh21_slope=len_f21?(Z_ghost(f21+1)-Z_ghost(f21))/len_f21:0,Zh22_slope=len_f22?(Z_ghost(f22+1)-Z_ghost(f22))/len_f22:0,Vh1_kp1=Vh1_lo+(k+1)*Vh1_slope;
                    if(Vh1_lo>0){
                        velocity_profile.Add_Linear_Product_Profile(0,k,Vh1_lo,Vh1_slope,Z_ghost(f21)+(jmin-face21.j())*Zh21_slope,Zh21_slope);
                        velocity_profile.Add_Linear_Product_Profile(k+1,len_f-k-1,Vh1_kp1,Vh1_slope,Z_ghost(f22)+(jmin+k+1-face22.j())*Zh22_slope,Zh22_slope);}
                    else{
                        velocity_profile.Add_Linear_Product_Profile(0,k,Vh1_lo,Vh1_slope,Z_ghost(f22)+(jmin-face22.j())*Zh22_slope,Zh22_slope);
                        velocity_profile.Add_Linear_Product_Profile(k+1,len_f-k-1,Vh1_kp1,Vh1_slope,Z_ghost(f21)+(jmin+k+1-face21.j())*Zh21_slope,Zh21_slope);}}
                else{ // no zero crossing
                    if(Vh1_lo>0){
                        T Zh21_slope=len_f21?(Z_ghost(f21+1)-Z_ghost(f21))/len_f21:0;
                        velocity_profile.Add_Linear_Product_Profile(0,len_f,Vh1_lo,Vh1_slope,Z_ghost(f21)+(jmin-face21.j())*Zh21_slope,Zh21_slope);}
                    else{
                        T Zh22_slope=len_f22?(Z_ghost(f22+1)-Z_ghost(f22))/len_f22:0;
                        velocity_profile.Add_Linear_Product_Profile(0,len_f,Vh1_lo,Vh1_slope,Z_ghost(f22)+(jmin-face22.j())*Zh22_slope,Zh22_slope);}}

                RLE_LINEAR_PROFILE<T> divergence_profile;
                divergence_profile.Add_Linear_Profile(0,len_f,Vh1_lo,Vh1_slope);
                velocity_profile.Update_Linear_Profile(dj1,len_f21,-dt_over_length,Z(f21),Z(f21+1));
                divergence_profile.Update_Linear_Profile(dj1,len_f21,-dt_over_length,divergence(f21),divergence(f21+1));
                velocity_profile.Update_Linear_Profile(dj2,len_f22,+dt_over_length,Z(f22),Z(f22+1));
                divergence_profile.Update_Linear_Profile(dj2,len_f22,+dt_over_length,divergence(f22),divergence(f22+1));}

            if(face11.jmax()==face12.jmax()){face11++;face12++;}else if(face11.jmax()<face12.jmax()) face11++;else face12++;
            if(face21.jmax()==face22.jmax()){face21++;face22++;}else if(face21.jmax()<face22.jmax()) face21++;else face22++;}}
}
//#####################################################################
// Function Apply_Vertical_Horizontal_Velocity_Flux
//#####################################################################
template<class T_FACE_H,class T,class T_GRID> static void
Apply_Vertical_Horizontal_Velocity_Flux(const T_GRID& grid,ARRAY<T>& Z,const ARRAY<T>& Z_ghost,const ARRAY<T>& V,const T dt,ARRAY<T>& divergence)
{
    assert(grid.long_run_cells==2);
    T dt_over_dy=dt/grid.uniform_grid.dX.y;
    int previous_face_length=0;
    if(grid.long_run_faces_horizontal==1){
        for(T_FACE_H face(grid,0);face;face++){
            if(!face.First_In_Column()){
                int fh1=face.Face(),fh0=fh1-1,fy1=face.cell1.Face_Y(),fy2=face.cell2.Face_Y(),j1=face.cell1.j,j2=face.cell2.j;
                T v_mid=(T).5*(V(fy1+(j1<j2))+V(fy2+(j2<j1)));
                T Zh_donor=v_mid>0?Z_ghost(fh0):Z_ghost(fh1);
                T flux_volume=dt_over_dy*v_mid;
                T flux=flux_volume*Zh_donor;
                Z(fh0)-=flux/previous_face_length;Z(fh1)+=flux/face.Length();
                divergence(fh0)-=flux_volume/previous_face_length;divergence(fh1)+=flux_volume/face.Length();}
            previous_face_length=face.Length();}}
    else{
        for(T_FACE_H face(grid,0);face;face++){
            if(!face.First_In_Column()){
                int fh1=face.Face(),fh0=fh1-1,fy1=face.cell1.Face_Y(),fy2=face.cell2.Face_Y(),j1=face.cell1.j,j2=face.cell2.j;
                T v_mid_bottom=(T).5*(V(fy1+(j1<j2))+V(fy2+(j2<j1)));
                T flux_volume_bottom=dt_over_dy*v_mid_bottom,flux_bottom=flux_volume_bottom*(v_mid_bottom>0?Z_ghost(fh0):Z_ghost(fh1));

                // deal with previous face
                if(previous_face_length==1){Z(fh0)-=flux_bottom;divergence(fh0)-=flux_volume_bottom;}
                else{
                    RLE_LINEAR_PROFILE<T> previous_velocity_profile,previous_divergence_profile;
                    previous_velocity_profile.Add_Single_Value(previous_face_length-1,flux_bottom);
                    previous_velocity_profile.Update_Long_Linear_Profile(0,previous_face_length-1,-1,Z(fh0-1),Z(fh0));
                    previous_divergence_profile.Add_Single_Value(previous_face_length-1,flux_volume_bottom);
                    previous_divergence_profile.Update_Long_Linear_Profile(0,previous_face_length-1,-1,divergence(fh0-1),divergence(fh0));}

                if(face.Short()){Z(fh1)+=flux_bottom;divergence(fh1)+=flux_volume_bottom;}
                else{
                    int len_f=face.Length()-1;
                    T v_mid_interior=(T).5*(V(fy1+1)+V(fy2+1));
                    T Zh_slope=(Z(fh1+1)-Z(fh1))/len_f;
                    T flux_volume_interior=dt_over_dy*v_mid_interior;
                    T flux_differential=-flux_volume_interior*Zh_slope;

                    RLE_LINEAR_PROFILE<T> velocity_profile,divergence_profile;
                    if(face.Length()>2) velocity_profile.Add_Constant_Profile(1,face.Length()-3,flux_differential);
                    // fix ends
                    velocity_profile.Add_Single_Value(0,flux_bottom);
                    velocity_profile.Add_Single_Value(0,-flux_volume_interior*(v_mid_interior>0?Z(fh1):Z(fh1)+Zh_slope));
                    velocity_profile.Add_Single_Value(len_f,flux_volume_interior*(v_mid_interior>0?Z(fh1)+(len_f-1)*Zh_slope:Z(fh1+1)));
                    divergence_profile.Add_Single_Value(0,flux_volume_bottom);
                    divergence_profile.Add_Single_Value(0,-flux_volume_interior);
                    divergence_profile.Add_Single_Value(len_f,flux_volume_interior);
                    velocity_profile.Update_Long_Linear_Profile(0,len_f,+1,Z(fh1),Z(fh1+1));
                    divergence_profile.Update_Long_Linear_Profile(0,len_f,+1,divergence(fh1),divergence(fh1+1));}}
            previous_face_length=face.Length();}}
}
//#####################################################################
// Function Apply_Horizontal_Vertical_Velocity_Flux
//#####################################################################
template<class T_FACE_H,class T,class T_GRID> static void
Apply_Horizontal_Vertical_Velocity_Flux(const T_GRID& grid,ARRAY<T>& Z,const ARRAY<T>& Z_ghost,const ARRAY<T>& V,const T dt,ARRAY<T>& divergence)
{
    assert(grid.long_run_cells==2);
    T dt_over_length=dt/grid.uniform_grid.dX.x;
    if(grid.long_run_faces_horizontal==1){
        for(T_FACE_H face(grid,0);face;face++)if(!face.First_In_Column()){
            int fh1=face.Face(),fh0=fh1-1,fy1=face.cell1.Face_Y(),fy2=face.cell2.Face_Y(),j1=face.cell1.j,j2=face.cell2.j;
            T Vh_mid1=(T).5*(V(fh0)+V(fh1));
            T Zv_donor1=Vh_mid1>0?Z_ghost(fy1+(j1<j2)):Z_ghost(fy2+(j2<j1));
            T flux_volume1=dt_over_length*Vh_mid1;
            T flux1=flux_volume1*Zv_donor1;
            if(j1<j2) Z(fy1+1)-=flux1/(face.cell1.length-1);else Z(fy1)-=flux1;
            if(j2<j1) Z(fy2+1)+=flux1/(face.cell2.length-1);else Z(fy2)+=flux1;
            if(j1<j2) divergence(fy1+1)-=flux_volume1/(face.cell1.length-1);else divergence(fy1)-=flux_volume1;
            if(j2<j1) divergence(fy2+1)+=flux_volume1/(face.cell2.length-1);else divergence(fy2)+=flux_volume1;
            if(face.Long()){
                T Vh_mid2=V(fh1);
                T Zv_donor2=Vh_mid2>0?Z_ghost(fy1+1):Z_ghost(fy2+1);
                T flux_volume2=dt_over_length*Vh_mid2*(face.Length()-1);
                T flux2=flux_volume2*Zv_donor2;
                Z(fy1+1)-=flux2/(face.cell1.length-1);Z(fy2+1)+=flux2/(face.cell2.length-1);
                divergence(fy1+1)-=flux_volume2/(face.cell1.length-1);divergence(fy2+1)+=flux_volume2/(face.cell2.length-1);}}}
    else{
        for(T_FACE_H face(grid,0);face;face++)if(!face.First_In_Column()){
            int fh1=face.Face(),fh0=fh1-1,fy1=face.cell1.Face_Y(),fy2=face.cell2.Face_Y(),j1=face.cell1.j,j2=face.cell2.j;
            T Vh_mid1=(T).5*(V(fh0)+V(fh1));
            T Zv_donor1=Vh_mid1>0?Z_ghost(fy1+(j1<j2)):Z_ghost(fy2+(j2<j1));
            T flux_volume1=dt_over_length*Vh_mid1;
            T flux1=flux_volume1*Zv_donor1;
            if(j1<j2) Z(fy1+1)-=flux1/(face.cell1.length-1);else Z(fy1)-=flux1;
            if(j2<j1) Z(fy2+1)+=flux1/(face.cell2.length-1);else Z(fy2)+=flux1;
            if(j1<j2) divergence(fy1+1)-=flux_volume1/(face.cell1.length-1);else divergence(fy1)-=flux_volume1;
            if(j2<j1) divergence(fy2+1)+=flux_volume1/(face.cell2.length-1);else divergence(fy2)+=flux_volume1;
            if(face.Long()){
                int len_f=face.Length()-1;
                T Vh_mid_slope=(V(fh1+1)-V(fh1))/len_f,Vh_mid_half_slope=(T).5*Vh_mid_slope,Vh_mid_lo=V(fh1)+Vh_mid_half_slope;
                // sum up 0 to len_f-1
                T Vh_sum=len_f*(Vh_mid_lo+(len_f-1)*Vh_mid_half_slope);
                T flux_volume2=dt_over_length*Vh_sum,flux2;
                if(Vh_mid_lo*(V(fh1+1)-Vh_mid_half_slope)<0){ // have zero crossing
                    int k=(int)(-Vh_mid_lo/Vh_mid_slope);assert(0<=k&&k<len_f-1);
                    // 0 to k goes one direction, and k+1 to len_f-1 goes the other direction
                    T Vh_sum_sign1=(k+1)*(Vh_mid_lo+k*Vh_mid_half_slope),Vh_sum_sign2=Vh_sum-Vh_sum_sign1;
                    flux2=dt_over_length*(Vh_sum_sign1>0?(Vh_sum_sign1*Z_ghost(fy1+1)+Vh_sum_sign2*Z_ghost(fy2+1)):(Vh_sum_sign1*Z_ghost(fy2+1)+Vh_sum_sign2*Z_ghost(fy1+1)));}
                else{flux2=flux_volume2*(Vh_sum>0?Z_ghost(fy1+1):Z_ghost(fy2+1));}
                Z(fy1+1)-=flux2/(face.cell1.length-1);Z(fy2+1)+=flux2/(face.cell2.length-1);
                divergence(fy1+1)-=flux_volume2/(face.cell1.length-1);divergence(fy2+1)+=flux_volume2/(face.cell2.length-1);}}}
}
//#####################################################################
// Function Apply_Vertical_Vertical_Velocity_Flux
//#####################################################################
template<class T,class T_GRID> static void
Apply_Vertical_Vertical_Velocity_Flux(const T_GRID& grid,ARRAY<T>& Z,const ARRAY<T>& Z_ghost,const ARRAY<T>& V,const T dt,ARRAY<T>& divergence)
{
    T dt_over_dy=dt/grid.uniform_grid.dX.y;
    for(typename T_GRID::CELL_ITERATOR cell(grid,0);cell;cell++){
        int fy1=cell.Face_Y(),fy2=fy1+1,fy3=fy1+2;
        T v_mid1=(T).5*(V(fy1)+V(fy2));
        int donor1=v_mid1>0?fy1:fy2;
        T flux_volume1=dt_over_dy*V(donor1);
        T flux1=flux_volume1*Z_ghost(donor1);
        Z(fy1)-=flux1;
        divergence(fy1)-=flux_volume1;
        if(cell.Short()){
            Z(fy2)+=flux1;
            divergence(fy2)+=flux_volume1;}
        else{
            T v_mid2=(T).5*(V(fy2)+V(fy3));
            int donor2=v_mid2>0?fy2:fy3;
            T flux_volume2=dt_over_dy*V(donor2);
            T flux2=flux_volume2*Z_ghost(donor2);
            Z(fy2)+=(flux1-flux2)/(cell.length-1);Z(fy3)+=flux2;
            divergence(fy2)+=(flux_volume1-flux_volume2)/(cell.length-1);divergence(fy3)+=flux_volume2;}}
}
//#####################################################################
// Function Euler_Step
//#####################################################################
static const bool upwind=true; // TODO: remove

template<class T> void ADVECTION_FACE_VELOCITY_DONOR_CELL_RLE<T>::
Euler_Step(const RLE_GRID_2D<T>& grid,ARRAY<T>& Z,const ARRAY<T>& Z_ghost,const ARRAY<T>& V,const T dt)
{
    ARRAY<T> divergence(grid.number_of_faces);

    if(!upwind) Apply_Horizontal_Horizontal_Velocity_Flux<typename RLE_GRID_2D<T>::FACE_X_ITERATOR,typename RLE_GRID_2D<T>::FACE_X_ITERATOR>(grid,Z,Z_ghost,V,dt,divergence);
    else Apply_Horizontal_Velocity_Flux<typename RLE_GRID_2D<T>::FACE_X_ITERATOR>(grid,Z,Z_ghost,V,dt,divergence);

    Apply_Vertical_Horizontal_Velocity_Flux<typename RLE_GRID_2D<T>::FACE_X_ITERATOR>(grid,Z,Z_ghost,V,dt,divergence);
    Apply_Horizontal_Vertical_Velocity_Flux<typename RLE_GRID_2D<T>::FACE_X_ITERATOR>(grid,Z,Z_ghost,V,dt,divergence);
    Apply_Vertical_Vertical_Velocity_Flux(grid,Z,Z_ghost,V,dt,divergence);

    if(clamp_divergence_fix) for(int f=1;f<=Z.m;f++)Z(f)/=max((T).9,1+divergence(f));
    else for(int f=1;f<=Z.m;f++)Z(f)/=1+divergence(f);

    LOG::cout<<"maximum advection divergence = "<<divergence.Max()<<std::endl;
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void ADVECTION_FACE_VELOCITY_DONOR_CELL_RLE<T>::
Euler_Step(const RLE_GRID_3D<T>& grid,ARRAY<T>& Z,const ARRAY<T>& Z_ghost,const ARRAY<T>& V,const T dt)
{
    ARRAY<T> divergence(grid.number_of_faces);

    if(!upwind){
        Apply_Horizontal_Horizontal_Velocity_Flux<typename RLE_GRID_3D<T>::FACE_X_ITERATOR,typename RLE_GRID_3D<T>::FACE_X_ITERATOR>(grid,Z,Z_ghost,V,dt,divergence);
        Apply_Horizontal_Horizontal_Velocity_Flux<typename RLE_GRID_3D<T>::FACE_Z_ITERATOR,typename RLE_GRID_3D<T>::FACE_Z_ITERATOR>(grid,Z,Z_ghost,V,dt,divergence);}
    else{
        Apply_Horizontal_Velocity_Flux<typename RLE_GRID_3D<T>::FACE_X_ITERATOR>(grid,Z,Z_ghost,V,dt,divergence);
        Apply_Horizontal_Velocity_Flux<typename RLE_GRID_3D<T>::FACE_Z_ITERATOR>(grid,Z,Z_ghost,V,dt,divergence);}

    Apply_Horizontal_Horizontal_Velocity_Flux<typename RLE_GRID_3D<T>::FACE_Z_ITERATOR,typename RLE_GRID_3D<T>::FACE_X_ITERATOR>(grid,Z,Z_ghost,V,dt,divergence);
    Apply_Horizontal_Horizontal_Velocity_Flux<typename RLE_GRID_3D<T>::FACE_X_ITERATOR,typename RLE_GRID_3D<T>::FACE_Z_ITERATOR>(grid,Z,Z_ghost,V,dt,divergence);
    Apply_Vertical_Horizontal_Velocity_Flux<typename RLE_GRID_3D<T>::FACE_X_ITERATOR>(grid,Z,Z_ghost,V,dt,divergence);
    Apply_Vertical_Horizontal_Velocity_Flux<typename RLE_GRID_3D<T>::FACE_Z_ITERATOR>(grid,Z,Z_ghost,V,dt,divergence);
    Apply_Horizontal_Vertical_Velocity_Flux<typename RLE_GRID_3D<T>::FACE_X_ITERATOR>(grid,Z,Z_ghost,V,dt,divergence);
    Apply_Horizontal_Vertical_Velocity_Flux<typename RLE_GRID_3D<T>::FACE_Z_ITERATOR>(grid,Z,Z_ghost,V,dt,divergence);
    Apply_Vertical_Vertical_Velocity_Flux(grid,Z,Z_ghost,V,dt,divergence);

    if(clamp_divergence_fix) for(int f=1;f<=Z.m;f++)Z(f)/=max((T).9,1+divergence(f));
    else for(int f=1;f<=Z.m;f++)Z(f)/=1+divergence(f);

    LOG::cout<<"maximum advection divergence = "<<divergence.Max()<<std::endl;
}
//#####################################################################
template class ADVECTION_FACE_VELOCITY_DONOR_CELL_RLE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ADVECTION_FACE_VELOCITY_DONOR_CELL_RLE<double>;
#endif
#endif

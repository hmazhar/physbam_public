//#####################################################################
// Copyright 2002-2006, Ron Fedkiw, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Fourier_Transforms/FFT.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/exchange.h>
using namespace PhysBAM;
using ::std::sin;
//#####################################################################
// Function NR_fourn
//#####################################################################
// from Numerical Recipes
void PhysBAM::NR_fourn(const int isign,const ARRAY<int>& dim,ARRAY<float>& data)
{
    int idim;
    unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
    unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
    float tempi,tempr;
    double theta,wi,wpi,wpr,wr,wtemp;

    for(int i=1;i<=dim.m;i++)if(dim(i)&(dim(i)-1)){
        LOG::cerr<<"Input dimensions must be powers of two: "<<dim;PHYSBAM_FATAL_ERROR();}

    for(ntot=1,idim=1;idim<=dim.m;idim++) ntot*=dim(idim);
    nprev=1;
    for(idim=dim.m;idim>=1;idim--){
        n=dim(idim);nrem=ntot/(n*nprev);ip1=nprev<<1;ip2=ip1*n;ip3=ip2*nrem;i2rev=1;
        for(i2=1;i2<=ip2;i2+=ip1){
            if(i2 < i2rev) for(i1=i2;i1<=i2+ip1-2;i1+=2) for(i3=i1;i3<=ip3;i3+=ip2){i3rev=i2rev+i3-i2;exchange(data(i3),data(i3rev));exchange(data(i3+1),data(i3rev+1));}
            ibit=ip2 >> 1;
            while(ibit >= ip1 && i2rev > ibit){i2rev-=ibit;ibit>>=1;}
            i2rev += ibit;}
        ifp1=ip1;
        while(ifp1 < ip2){
            ifp2=ifp1<<1;theta=isign*6.28318530717959/(ifp2/ip1);wtemp=sin(.5*theta);wpr=-2*wtemp*wtemp;wpi=sin(theta);wr=1;wi=0;
            for(i3=1;i3<=ifp1;i3+=ip1){
                for(i1=i3;i1<=i3+ip1-2;i1+=2) for(i2=i1;i2<=ip3;i2+=ifp2){
                    k1=i2;k2=k1+ifp1;tempr=(float)wr*data(k2)-(float)wi*data(k2+1);tempi=(float)wr*data(k2+1)+(float)wi*data(k2);
                    data(k2)=data(k1)-tempr;data(k2+1)=data(k1+1)-tempi;data(k1)+=tempr;data(k1+1)+=tempi;}
                wr=(wtemp=wr)*wpr-wi*wpi+wr;wi=wi*wpr+wtemp*wpi+wi;}
            ifp1=ifp2;}
        nprev*=n;}
}

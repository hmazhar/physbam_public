//#####################################################################
// Copyright 2002-2008, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Duc Nguyen, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RUNGEKUTTA
//#####################################################################
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
namespace PhysBAM{
//#####################################################################
// Function Start
//#####################################################################
template<class T> void RUNGEKUTTA_CORE<T>::
Start(const T dt_input)
{
    dt=dt_input;substep=0;
    if(order>1)
        u_copy=u;
}
//#####################################################################
// Function Main
//#####################################################################
template<class T> T RUNGEKUTTA_CORE<T>::
Main()
{
    substep+=1;
    if(substep==1){if(!pseudo_time) time+=dt;}
    else{
        if(substep==2 && order==2)
            u=(T).5*(u_copy+u);
        else if(substep==2 && order==3){
            if(!pseudo_time) time-=dt/2;
            u=(T).75*u_copy+(T).25*u;}
        else if(substep==3){
            if(!pseudo_time) time+=dt/2;
            u=(T)1./3*u_copy+(T)2./3*u;}
        if(enforce_boundary_condition) (*enforce_boundary_condition)(u,time);}
    return time;
}
//#####################################################################
template class RUNGEKUTTA_CORE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RUNGEKUTTA_CORE<double>;
#endif
}

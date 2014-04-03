//#####################################################################
// Copyright 2002-2008, Ronald Fedkiw, Jon Gretarsson, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POISSON
//#####################################################################
#ifndef __POISSON__
#define __POISSON__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

template<class T>
class POISSON:public NONCOPYABLE
{
protected:
    T beta_minus,beta_plus; // constant beta for phi <= 0 and phi > 0
    ARRAY<T> beta_multiphase; // constant beta for multiphase flows
    bool GFM; // (true) for GFM, (false) for delta function smearing
    T number_of_interface_cells; // e.g. 3 interface cells gives half_width = 1.5*dx
    bool smear_beta; // (true) for smearing beta, (false) for smearing 1/beta
public:
    bool u_jumps,beta_un_jumps; // true when [u] != 0 or [beta un] != 0
    bool use_variable_beta; // otherwise using piecewise constant beta
    bool use_weighted_divergence;
    bool beta_given_on_faces;
    bool multiphase;

public:
    POISSON(bool multiphase_input=false)
        :u_jumps(false),beta_un_jumps(false),use_variable_beta(false),use_weighted_divergence(false),beta_given_on_faces(false),multiphase(multiphase_input)
    {
        Use_GFM();
        Set_Constant_beta();
    }

    void Set_Constant_beta(const T beta_minus_input=1,const T beta_plus_input=1)
    {beta_minus=beta_minus_input;beta_plus=beta_plus_input;use_variable_beta=false;}

    void Set_Constant_beta(const ARRAY<T>& beta_multiphase_input)
    {beta_multiphase=beta_multiphase_input;}

    void Use_GFM()
    {GFM=true;number_of_interface_cells=0;}

    void Use_Delta_Function_Method(const T number_of_interface_cells_input=3)
    {number_of_interface_cells=number_of_interface_cells_input;GFM=false;Smear_beta();}

    void Smear_beta()
    {smear_beta=true;}

    void Smear_One_Over_beta()
    {smear_beta=false;}

//#####################################################################
};
}
#endif


//#####################################################################
// Copyright 2011, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ROTATED_STRESS_DERIVATIVE
//#####################################################################
#ifndef __ROTATED_STRESS_DERIVATIVE__
#define __ROTATED_STRESS_DERIVATIVE__

#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>

namespace PhysBAM{

template<class T,int d> class ROTATED_STRESS_DERIVATIVE;

template<class T>
class ROTATED_STRESS_DERIVATIVE<T,2>
{
public:
    T a1111,a1122,a2222;
    T a1212,a1221;

    MATRIX<T,2> dP_hat(const MATRIX<T,2>& dF_hat) const // Contraction
    {
        MATRIX<T,2> dP_hat;
        dP_hat(1,1)=a1111*dF_hat(1,1)+a1122*dF_hat(2,2);
        dP_hat(2,2)=a1122*dF_hat(1,1)+a2222*dF_hat(2,2);
        dP_hat(1,2)=a1212*dF_hat(1,2)+a1221*dF_hat(2,1);
        dP_hat(2,1)=a1221*dF_hat(1,2)+a1212*dF_hat(2,1);

        return dP_hat;
    }

    void Make_Positive_Definite()
    {
        SYMMETRIC_MATRIX<T,2> A(a1111,a1122,a2222);
        A=A.Positive_Definite_Part();
        a1111=A(1,1);a1122=A(1,2);a2222=A(2,2);

        T lambda1=std::max(.5*(a1212+a1221),0.),lambda2=std::max(.5*(a1212-a1221),0.);
        a1212=lambda1+lambda2;a1221=lambda1-lambda2;
    }

    ROTATED_STRESS_DERIVATIVE operator+=(const ROTATED_STRESS_DERIVATIVE& input)
    {
        a1111+=input.a1111;a1122+=input.a1122;a2222+=input.a2222;
        a1212+=input.a1212;a1221+=input.a1221;
        return *this;
    }

    ROTATED_STRESS_DERIVATIVE operator-=(const ROTATED_STRESS_DERIVATIVE& input)
    {
        a1111-=input.a1111;a1122-=input.a1122;a2222-=input.a2222;
        a1212-=input.a1212;a1221-=input.a1221;
        return *this;
    }

};

template<class T>
class ROTATED_STRESS_DERIVATIVE<T,3>
{
public:
    T a1111,a1122,a1133,a2222,a2233,a3333;
    T a1212,a1221,a1313,a1331,a2323,a2332;

    MATRIX<T,3> dP_hat(const MATRIX<T,3>& dF_hat) // Contraction
    {
        MATRIX<T,3> dP_hat;
        dP_hat(1,1)=a1111*dF_hat(1,1)+a1122*dF_hat(2,2)+a1133*dF_hat(3,3);
        dP_hat(2,2)=a1122*dF_hat(1,1)+a2222*dF_hat(2,2)+a2233*dF_hat(3,3);
        dP_hat(3,3)=a1133*dF_hat(1,1)+a2233*dF_hat(2,2)+a3333*dF_hat(3,3);
        dP_hat(1,2)=a1212*dF_hat(1,2)+a1221*dF_hat(2,1);
        dP_hat(2,1)=a1221*dF_hat(1,2)+a1212*dF_hat(2,1);
        dP_hat(1,3)=a1313*dF_hat(1,3)+a1331*dF_hat(3,1);
        dP_hat(3,1)=a1331*dF_hat(1,3)+a1313*dF_hat(3,1);
        dP_hat(2,3)=a2323*dF_hat(2,3)+a2332*dF_hat(3,2);
        dP_hat(3,2)=a2332*dF_hat(2,3)+a2323*dF_hat(3,2);

        return dP_hat;
    }

    void Make_Positive_Definite()
    {
        SYMMETRIC_MATRIX<T,3> A(a1111,a1122,a1133,a2222,a2233,a3333);
        A=A.Positive_Definite_Part();
	a1111=A(1,1);a1122=A(1,2);a1133=A(1,3);a2222=A(2,2);a2233=A(2,3);a3333=A(3,3);

        T lambda1=std::max(.5*(a1212+a1221),0.),lambda2=std::max(.5*(a1212-a1221),0.);
        a1212=lambda1+lambda2;a1221=lambda1-lambda2;

        lambda1=std::max(.5*(a1313+a1331),0.);lambda2=std::max(.5*(a1313-a1331),0.);
        a1313=lambda1+lambda2;a1331=lambda1-lambda2;

        lambda1=std::max(.5*(a2323+a2332),0.);lambda2=std::max(.5*(a2323-a2332),0.);
        a2323=lambda1+lambda2;a2332=lambda1-lambda2;
    }

    ROTATED_STRESS_DERIVATIVE operator+=(const ROTATED_STRESS_DERIVATIVE& input)
    {
        a1111+=input.a1111;a1122+=input.a1122;a1133+=input.a1133;
        a2222+=input.a2222;a2233+=input.a2233;a3333+=input.a3333;
        a1212+=input.a1212;a1221+=input.a1221;
        a1313+=input.a1313;a1331+=input.a1331;
        a2323+=input.a2323;a2332+=input.a2332;
        return *this;
    }

    ROTATED_STRESS_DERIVATIVE operator-=(const ROTATED_STRESS_DERIVATIVE& input)
    {
        a1111-=input.a1111;a1122-=input.a1122;a1133-=input.a1133;
        a2222-=input.a2222;a2233-=input.a2233;a3333-=input.a3333;
        a1212-=input.a1212;a1221-=input.a1221;
        a1313-=input.a1313;a1331-=input.a1331;
        a2323-=input.a2323;a2332-=input.a2332;
        return *this;
    }

};

}

#endif

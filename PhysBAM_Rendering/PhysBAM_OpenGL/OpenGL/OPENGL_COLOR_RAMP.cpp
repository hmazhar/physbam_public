//#####################################################################
// Copyright 2004, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
using namespace PhysBAM;
template<class T> OPENGL_COLOR_RAMP<T>::
OPENGL_COLOR_RAMP() 
{
}
template<class T> OPENGL_COLOR_RAMP<T>::
~OPENGL_COLOR_RAMP() 
{
}
template<class T> OPENGL_COLOR OPENGL_COLOR_RAMP<T>::
Lookup(T x) const
{
    int left_index=0,right_index=0;
    for(int i=1;i<=color_x.m;i++){
        if(x>color_x(i))left_index=i;
        else if(x<color_x(i)){right_index=i;break;}
        else return equal_colors(i);}
    if(left_index&&right_index){T alpha=(x-color_x(left_index))/(color_x(right_index)-color_x(left_index));
        return T(alpha)*less_colors(right_index)+T(1.0-alpha)*greater_colors(left_index);}
    else if(left_index)return greater_colors(left_index);
    else if(right_index)return less_colors(right_index);
    return OPENGL_COLOR(1,0,0);
}
template<class T> void OPENGL_COLOR_RAMP<T>::
Add_Color(T x,const OPENGL_COLOR& less_color,const OPENGL_COLOR& exact_color,const OPENGL_COLOR& greater_color)
{
    assert(color_x.m==0||x>color_x(color_x.m));
    color_x.Append(x);
    less_colors.Append(less_color);
    equal_colors.Append(exact_color);
    greater_colors.Append(greater_color);
}
template<class T> OPENGL_COLOR_RAMP<T>* OPENGL_COLOR_RAMP<T>::
Matlab_Jet(T value_min,T value_max)
{
    OPENGL_COLOR_RAMP<T> *jet=new OPENGL_COLOR_RAMP<T>;
    T interval_width=value_max-value_min;
    jet->Add_Color(interval_width*0+value_min,OPENGL_COLOR(0,0,0.5608));
    jet->Add_Color(interval_width*0.1406+value_min,OPENGL_COLOR(0,0,1));
    jet->Add_Color(interval_width*0.3594+value_min,OPENGL_COLOR(0,1,1));
    jet->Add_Color(interval_width*0.6094+value_min,OPENGL_COLOR(1,1,0));
    jet->Add_Color(interval_width*0.8594+value_min,OPENGL_COLOR(1,0,0));
    jet->Add_Color(interval_width*1+value_min,OPENGL_COLOR(0.5,0,0));
    return jet;
}
template<class T> OPENGL_COLOR_RAMP<T>* OPENGL_COLOR_RAMP<T>::
Matlab_Hot(T value_min,T value_max)
{
    OPENGL_COLOR_RAMP<T> *hot=new OPENGL_COLOR_RAMP<T>;
    T interval_width=value_max-value_min;
    hot->Add_Color(interval_width*0+value_min,OPENGL_COLOR(0,0,0,0));
    hot->Add_Color(interval_width*0.3750+value_min,OPENGL_COLOR(1,0,0));
    hot->Add_Color(interval_width*0.7656+value_min,OPENGL_COLOR(1,1,0));
    hot->Add_Color(interval_width*1+value_min,OPENGL_COLOR(1,1,1));
    return hot;
}
template<class T> OPENGL_COLOR_RAMP<T>* OPENGL_COLOR_RAMP<T>::
Two_Color_Ramp(T value_min,T value_max,const OPENGL_COLOR& color_min,const OPENGL_COLOR& color_min_exact,const OPENGL_COLOR& color_max,const OPENGL_COLOR& color_max_exact)
{
    OPENGL_COLOR_RAMP<T> *ramp=new OPENGL_COLOR_RAMP<T>;
    ramp->Add_Color(value_min,color_min,color_min_exact);ramp->Add_Color(value_max,color_max,color_max_exact);
    return ramp;
}
template<class T> OPENGL_COLOR_RAMP<T>* OPENGL_COLOR_RAMP<T>::
Two_Color_Ramp(T value_min,T value_max,const OPENGL_COLOR& color_min,const OPENGL_COLOR& color_max)
{
    OPENGL_COLOR_RAMP<T> *ramp=new OPENGL_COLOR_RAMP<T>;
    ramp->Add_Color(value_min,color_min,color_min);ramp->Add_Color(value_max,color_max,color_max);
    return ramp;
}
template<class T> OPENGL_COLOR_RAMP<T>* OPENGL_COLOR_RAMP<T>::
Levelset_Color_Constant_Ramp(const OPENGL_COLOR& negative_color,const OPENGL_COLOR& positive_color)
    {OPENGL_COLOR_RAMP<T> *ramp=new OPENGL_COLOR_RAMP<T>;
    ramp->Add_Color(0,negative_color,negative_color,positive_color);
    return ramp;
}
template<class T> OPENGL_COLOR_RAMP<T>* OPENGL_COLOR_RAMP<T>::
Levelset_Color_Linear_Ramp(const OPENGL_COLOR& negative_color,const OPENGL_COLOR& positive_color,T abs_value_max)
{
    OPENGL_COLOR_RAMP<T> *ramp=new OPENGL_COLOR_RAMP<T>;
    if(abs_value_max>0) ramp->Add_Color(-abs_value_max,OPENGL_COLOR::Gray(0,0));
    ramp->Add_Color(0,negative_color,negative_color,positive_color);
    if(abs_value_max>0) ramp->Add_Color(abs_value_max,OPENGL_COLOR::Gray(0,0));
    return ramp;
}
template class OPENGL_COLOR_RAMP<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COLOR_RAMP<double>;
#endif

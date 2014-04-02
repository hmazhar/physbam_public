#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
using namespace PhysBAM;
template<class T,class T2> T2 LINEAR_INTERPOLATION<T,T2>::
Bilinear(const T2& u1,const T2& u2,const T2& u3,const T2& u4,const VECTOR<T,2>& minimum_corner,const VECTOR<T,2>& maximum_corner,const VECTOR<T,2>& X)
{
    T one_over_x_right_minus_x_left=1/(maximum_corner.x-minimum_corner.x);
    T2 u_bottom=Linear_Predivided(minimum_corner.x,one_over_x_right_minus_x_left,u1,u2,X.x),
        u_top=Linear_Predivided(minimum_corner.x,one_over_x_right_minus_x_left,u3,u4,X.x);
    return Linear(minimum_corner.y,maximum_corner.y,u_bottom,u_top,X.y);
}
template<class T,class T2> T2 LINEAR_INTERPOLATION<T,T2>::
Bilinear(const T2& u1,const T2& u2,const T2& u3,const T2& u4,const VECTOR<T,2>& X)
{
    T2 u_bottom=Linear_Normalized(u1,u2-u1,X.x),u_top=Linear_Normalized(u3,u4-u3,X.x);
    return Linear_Normalized(u_bottom,u_top-u_bottom,X.y);
}
template<class T,class T2> T2 LINEAR_INTERPOLATION<T,T2>::
Bilinear(const T2& u1,const T2& u3,T one_over_y_top_minus_y_bottom,const T x_left,const T y_bottom,const T2& slope12,const T2& slope34,const VECTOR<T,2>& X)
{
    T2 u_bottom=Linear(x_left,u1,slope12,X.x),u_top=Linear(x_left,u3,slope34,X.x);
    return Linear_Predivided(y_bottom,one_over_y_top_minus_y_bottom,u_bottom,u_top,X.y);
}
template<class T,class T2> T2 LINEAR_INTERPOLATION<T,T2>::
Trilinear(const T2& u1,const T2& u2,const T2& u3,const T2& u4,const T2& u5,const T2& u6,const T2& u7,const T2& u8,const VECTOR<T,3>& minimum_corner,const VECTOR<T,3>& maximum_corner,
    const VECTOR<T,3>& X)
{
    T one_over_x_right_minus_x_left=1/(maximum_corner.x-minimum_corner.x),one_over_y_right_minus_y_left=1/(maximum_corner.y-minimum_corner.y);
    T2 u_bottom=Linear_Predivided(minimum_corner.x,one_over_x_right_minus_x_left,u1,u2,X.x),
        u_top=Linear_Predivided(minimum_corner.x,one_over_x_right_minus_x_left,u3,u4,X.x);
    T2 u_front=Linear_Predivided(minimum_corner.y,one_over_y_right_minus_y_left,u_bottom,u_top,X.y);
    u_bottom=Linear_Predivided(minimum_corner.x,one_over_x_right_minus_x_left,u5,u6,X.x);
    u_top=Linear_Predivided(minimum_corner.x,one_over_x_right_minus_x_left,u7,u8,X.x);    
    T2 u_back=Linear_Predivided(minimum_corner.y,one_over_y_right_minus_y_left,u_bottom,u_top,X.y);
    return Linear(minimum_corner.z,maximum_corner.z,u_front,u_back,X.z);
}
template<class T,class T2> T2 LINEAR_INTERPOLATION<T,T2>::
Trilinear(const T2& u1,const T2& u2,const T2& u3,const T2& u4,const T2& u5,const T2& u6,const T2& u7,const T2& u8,const VECTOR<T,3>& X)
{
    T2 u_bottom=Linear_Normalized(u1,u2-u1,X.x),u_top=Linear_Normalized(u3,u4-u3,X.x),u_front=Linear_Normalized(u_bottom,u_top-u_bottom,X.y);
    u_bottom=Linear_Normalized(u5,u6-u5,X.x);u_top=Linear_Normalized(u7,u8-u7,X.x);T2 u_back=Linear_Normalized(u_bottom,u_top-u_bottom,X.y);
    return Linear_Normalized(u_front,u_back-u_front,X.z);
}
template<class T,class T2> T2 LINEAR_INTERPOLATION<T,T2>::
Trilinear(const T2& u1,const T2& u3,const T2& u5,const T2& u7,T one_over_y_top_minus_y_bottom,T one_over_z_back_minus_z_front,const T x_left,const T y_bottom,const T z_front,
    const T2& slope12,const T2& slope34,const T2& slope56,const T2& slope78,const VECTOR<T,3>& X)
{
    T2 u_bottom=Linear(x_left,u1,slope12,X.x),u_top=Linear(x_left,u3,slope34,X.x);
    T2 u_front=Linear_Predivided(y_bottom,one_over_y_top_minus_y_bottom,u_bottom,u_top,X.y);
    u_bottom=Linear(x_left,u5,slope56,X.x);u_top=Linear(x_left,u7,slope78,X.x);    
    T2 u_back=Linear_Predivided(y_bottom,one_over_y_top_minus_y_bottom,u_bottom,u_top,X.y);
    return Linear_Predivided(z_front,one_over_z_back_minus_z_front,u_front,u_back,X.z);
}

template class LINEAR_INTERPOLATION<float,MATRIX<float,3,3> >;
template class LINEAR_INTERPOLATION<float,SYMMETRIC_MATRIX<float,2> >;
template class LINEAR_INTERPOLATION<float,SYMMETRIC_MATRIX<float,3> >;
template class LINEAR_INTERPOLATION<float,VECTOR<float,2> >;
template class LINEAR_INTERPOLATION<float,VECTOR<float,3> >;
template class LINEAR_INTERPOLATION<float,VECTOR<float,4> >;
template class LINEAR_INTERPOLATION<float,VECTOR<float,5> >;
template class LINEAR_INTERPOLATION<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LINEAR_INTERPOLATION<double,MATRIX<double,3,3> >;
template class LINEAR_INTERPOLATION<double,SYMMETRIC_MATRIX<double,2> >;
template class LINEAR_INTERPOLATION<double,SYMMETRIC_MATRIX<double,3> >;
template class LINEAR_INTERPOLATION<double,VECTOR<double,2> >;
template class LINEAR_INTERPOLATION<double,VECTOR<double,3> >;
template class LINEAR_INTERPOLATION<double,VECTOR<double,4> >;
template class LINEAR_INTERPOLATION<double,VECTOR<double,5> >;
template class LINEAR_INTERPOLATION<double,double>;
#endif

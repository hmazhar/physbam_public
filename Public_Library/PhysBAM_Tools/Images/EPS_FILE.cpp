#include <PhysBAM_Tools/Images/EPS_FILE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
using namespace PhysBAM;
template<class T> EPS_FILE<T>::
//#####################################################################
// Constructor
//#####################################################################
EPS_FILE(const std::string& filename,const RANGE<TV>& box)
    :stream(FILE_UTILITIES::Safe_Open_Output(filename,false,false)),bounding_box(RANGE<TV>::Empty_Box()),output_box(box)
{
    Write_Head();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> EPS_FILE<T>::
~EPS_FILE()
{
    Finish();
    delete stream;
}
//#####################################################################
// Function Finish
//#####################################################################
template<class T> void EPS_FILE<T>::
Finish()
{
    Write_Tail();
}
//#####################################################################
// Function Emit
//#####################################################################
template<class T> void EPS_FILE<T>::
Emit(const std::string& str)
{
    (*stream)<<str<<" ";
}
//#####################################################################
// Function Emit
//#####################################################################
template<class T> void EPS_FILE<T>::
Emit(const TV &pt)
{
    (*stream)<<pt.x<<" "<<pt.y<<" ";
}
//#####################################################################
// Function Bound
//#####################################################################
template<class T> void EPS_FILE<T>::
Bound(const TV& pt)
{
    bounding_box.Enlarge_To_Include_Point(pt);
}
//#####################################################################
// Function Bound
//#####################################################################
template<class T> template<int d> void EPS_FILE<T>::
Bound(const VECTOR<TV,d>& pts)
{
    for(int i=1;i<=d;i++) Bound(pts(i));
}
//#####################################################################
// Function Draw_Object_Colored
//#####################################################################
template<class T> template<class T_OBJECT> void EPS_FILE<T>::
Draw_Object_Colored(const T_OBJECT& object,const VECTOR<T,3>& color)
{
    (*stream)<<"gsave ";
    Line_Color(color);
    Draw_Object(object);
    (*stream)<<"grestore"<<std::endl;
}
//#####################################################################
// Function Line_Color
//#####################################################################
template<class T> void EPS_FILE<T>::
Line_Color(const VECTOR<T,3>& color)
{
    (*stream)<<color.x<<" "<<color.y<<" "<<color.z<<" setrgbcolor"<<std::endl;
}
//#####################################################################
// Function Write_Head
//#####################################################################
template<class T> void EPS_FILE<T>::
Write_Head()
{
    (*stream)<<"%!PS-Adobe-3.0 EPSF-3.0"<<std::endl;
    (*stream)<<"%%BoundingBox: ";
    Emit(output_box.min_corner);
    Emit(output_box.max_corner);
    (*stream)<<output_box.max_corner.y<<std::endl;
    (*stream)<<"/untransformed-image {"<<std::endl;
}
//#####################################################################
// Function Write_Tail
//#####################################################################
template<class T> void EPS_FILE<T>::
Write_Tail()
{
    T scale;
    TV shift;
    Compute_Transform(scale,shift);
    (*stream)<<"} def"<<std::endl;
    Emit(shift);
    (*stream)<<"translate"<<std::endl;
    (*stream)<<scale<<" dup dup scale /transform-scale exch def"<<std::endl;
    (*stream)<<"1 transform-scale div dup setlinewidth /line-width exch def"<<std::endl;
    (*stream)<<"untransformed-image"<<std::endl;
}
//#####################################################################
// Function Compute_Transform
//#####################################################################
template<class T> void EPS_FILE<T>::
Compute_Transform(T& scale,TV& shift)
{
    if(output_box.Empty() || bounding_box.Empty()){
        scale=1;
        shift=TV();}
    else{
        scale=(output_box.Edge_Lengths()/bounding_box.Edge_Lengths()).Min();
        shift=output_box.min_corner-bounding_box.min_corner*scale;}
}
//#####################################################################
// Function Draw_Point
//#####################################################################
template<class T> void EPS_FILE<T>::
Draw_Point(const TV &pt)
{
    Emit(pt);
    (*stream)<<"pointradius 0 360 arc stroke"<<std::endl;
    Bound(pt);
}
//#####################################################################
// Function Draw_Line
//#####################################################################
template<class T> void EPS_FILE<T>::
Draw_Line(const TV &a,const TV &b)
{
    Emit(a);
    (*stream)<<"moveto ";
    Emit(b);
    (*stream)<<"lineto stroke"<<std::endl;
    Bound(a);
    Bound(b);
}
//#####################################################################
// Function Draw_Object
//#####################################################################
template<class T> void EPS_FILE<T>::
Draw_Object(const RANGE<TV>& box)
{
    TV a(box.min_corner.x,box.max_corner.y),b(box.max_corner.x,box.min_corner.y);
    Draw_Line(box.min_corner,a);
    Draw_Line(box.min_corner,b);
    Draw_Line(box.max_corner,a);
    Draw_Line(box.max_corner,b);
}
template class EPS_FILE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EPS_FILE<double>;
#endif

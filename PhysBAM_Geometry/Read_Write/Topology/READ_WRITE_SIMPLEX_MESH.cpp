#include <PhysBAM_Geometry/Read_Write/Topology/READ_WRITE_SIMPLEX_MESH.h>
#include <PhysBAM_Geometry/Topology/POINT_SIMPLEX_MESH.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TETRAHEDRON_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>
using namespace PhysBAM;
template<class RW,class T_MESH> void Read_Write<T_MESH,RW,typename ENABLE_IF<IS_BASE_OF<SIMPLEX_MESH<T_MESH::dimension>,T_MESH>::value>::TYPE>::
Read(std::istream& input,SIMPLEX_MESH<d>& object)
{
    object.Clean_Memory();
    int backward_compatible;Read_Binary<RW>(input,object.number_nodes,backward_compatible,object.elements);
    if(object.elements.m){
        int min_index=ARRAYS_COMPUTATIONS::Min(object.elements.Flattened()),max_index=ARRAYS_COMPUTATIONS::Min(object.elements.Flattened());
        if(object.number_nodes<0) throw READ_ERROR(STRING_UTILITIES::string_sprintf("Invalid negative number_nodes = %d",object.number_nodes));
        if(min_index<1) throw READ_ERROR(STRING_UTILITIES::string_sprintf("Invalid vertex index %d",min_index));
        if(max_index>object.number_nodes) throw READ_ERROR(STRING_UTILITIES::string_sprintf("Invalid vertex index %d (number_nodes = %d)",min_index,object.number_nodes));}
}
template<class RW,class T_MESH> void Read_Write<T_MESH,RW,typename ENABLE_IF<IS_BASE_OF<SIMPLEX_MESH<T_MESH::dimension>,T_MESH>::value>::TYPE>::
Write(std::ostream& output,const SIMPLEX_MESH<d>& object)
{
    Write_Binary<RW>(output,object.number_nodes,d+1,object.elements);
}
template void Read_Write<POINT_SIMPLEX_MESH,float,void>::Read(std::basic_istream<char,std::char_traits<char> >&,SIMPLEX_MESH<0>&);
template void Read_Write<POINT_SIMPLEX_MESH,float,void>::Write(std::basic_ostream<char,std::char_traits<char> >&,SIMPLEX_MESH<0> const&);
template void Read_Write<SEGMENT_MESH,float,void>::Read(std::basic_istream<char,std::char_traits<char> >&,SIMPLEX_MESH<1>&);
template void Read_Write<SEGMENT_MESH,float,void>::Write(std::basic_ostream<char,std::char_traits<char> >&,SIMPLEX_MESH<1> const&);
template void Read_Write<TETRAHEDRON_MESH,float,void>::Read(std::basic_istream<char,std::char_traits<char> >&,SIMPLEX_MESH<3>&);
template void Read_Write<TETRAHEDRON_MESH,float,void>::Write(std::basic_ostream<char,std::char_traits<char> >&,SIMPLEX_MESH<3> const&);
template void Read_Write<TRIANGLE_MESH,float,void>::Read(std::basic_istream<char,std::char_traits<char> >&,SIMPLEX_MESH<2>&);
template void Read_Write<TRIANGLE_MESH,float,void>::Write(std::basic_ostream<char,std::char_traits<char> >&,SIMPLEX_MESH<2> const&);
template void Read_Write<POINT_SIMPLEX_MESH,double,void>::Read(std::basic_istream<char,std::char_traits<char> >&,SIMPLEX_MESH<0>&);
template void Read_Write<POINT_SIMPLEX_MESH,double,void>::Write(std::basic_ostream<char,std::char_traits<char> >&,SIMPLEX_MESH<0> const&);
template void Read_Write<SEGMENT_MESH,double,void>::Read(std::basic_istream<char,std::char_traits<char> >&,SIMPLEX_MESH<1>&);
template void Read_Write<SEGMENT_MESH,double,void>::Write(std::basic_ostream<char,std::char_traits<char> >&,SIMPLEX_MESH<1> const&);
template void Read_Write<TETRAHEDRON_MESH,double,void>::Read(std::basic_istream<char,std::char_traits<char> >&,SIMPLEX_MESH<3>&);
template void Read_Write<TETRAHEDRON_MESH,double,void>::Write(std::basic_ostream<char,std::char_traits<char> >&,SIMPLEX_MESH<3> const&);
template void Read_Write<TRIANGLE_MESH,double,void>::Read(std::basic_istream<char,std::char_traits<char> >&,SIMPLEX_MESH<2>&);
template void Read_Write<TRIANGLE_MESH,double,void>::Write(std::basic_ostream<char,std::char_traits<char> >&,SIMPLEX_MESH<2> const&);

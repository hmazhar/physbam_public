//#####################################################################
// Copyright 2003-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPARSE_VECTOR_ND
//#####################################################################
#ifndef __SPARSE_VECTOR_ND__
#define __SPARSE_VECTOR_ND__

#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
namespace PhysBAM{

template<class T>
class SPARSE_VECTOR_ND
{
public:
    int n;
    int number_of_active_indices;
    int* indices; // stores the locations of the entries
    T* x; // stores the data

     SPARSE_VECTOR_ND(const int n_input)
        :n(n_input)
    {
        number_of_active_indices=0; // no elements originally
        indices=new int[number_of_active_indices+1];
        x=new T[number_of_active_indices+1];
    }

    SPARSE_VECTOR_ND(const SPARSE_VECTOR_ND& vector)
        :n(vector.n)
    {
        number_of_active_indices=vector.number_of_active_indices;
        indices=new int[number_of_active_indices+1];
        x=new T[number_of_active_indices+1];
        for(int i=1;i<=number_of_active_indices;i++){indices[i]=vector.indices[i];x[i]=vector.x[i];}
    }

    ~SPARSE_VECTOR_ND()
    {delete[] x;delete[] indices;}

private:
    void operator=(const SPARSE_VECTOR_ND&);
public:

    const T operator()(const int i) const
    {assert(i>=1 && i<=n);for(int j=1;j<=number_of_active_indices;j++) if(indices[j] == i) return x[j];else if(indices[j]>i) break;
    return T();}

    void Set_Element(const int i,const T& element)
    {assert(i>=1 && i<=n);
    int j=1;for(;j<=number_of_active_indices;j++) if(indices[j] == i){x[j]=element;return;}else if(indices[j]>i) break;
    Insert_New_Element(i,j,element);}

    void Insert_New_Element(const int index,const int array_position,const T element=T())
    {number_of_active_indices++;
    T* new_x=new T[number_of_active_indices+1];int* new_indices=new int[number_of_active_indices+1];
    for(int i=1;i<array_position;i++){new_indices[i]=indices[i];new_x[i]=x[i];}
    new_indices[array_position]=index;new_x[array_position]=element;
    for(int i=array_position+1;i<=number_of_active_indices;i++){new_indices[i]=indices[i-1];new_x[i]=x[i-1];}
    delete[] indices;delete[] x;indices=new_indices;x=new_x;}

    void Add_Element(const int i,const T& element)
    {assert(i>=1 && i<=n);
    int j=1;for(;j<=number_of_active_indices;j++) if(indices[j] == i){x[j]+=element;return;}else if(indices[j]>i) break;
    Insert_New_Element(i,j,element);}

    bool Element_Present(const int i)
    {assert(i>=1 && i<=n);for(int j=1;j<=number_of_active_indices;j++) if(indices[j] == i) return true;else if(indices[j]>i) return false;
    return false;}

    bool Element_Present_And_Location(const int i,int& location)
    {assert(i>=1 && i<=n);for(int j=1;j<=number_of_active_indices;j++) if(indices[j] == i){location=j;return true;}else if(indices[j]>i) return false;
    return false;}

    void Clear()
    {delete[] indices;delete[] x;number_of_active_indices=0;indices=new int[number_of_active_indices+1];x=new T[number_of_active_indices+1];}

    T Dot_Product(const VECTOR_ND<T>& vector)
    {T sum=T();for(int i=1;i<=number_of_active_indices;i++) sum+=x[i]*vector(indices[i]);return sum;}

    void Negate()
    {for(int i=1;i<=number_of_active_indices;i++) x[i]=-x[i];}

    SPARSE_VECTOR_ND<T>& operator*=(const T a)
    {for(int i=1;i<=number_of_active_indices;i++)x[i]*=a;return *this;}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    void Write_Internal_Arrays(std::ostream& output_stream)
    {for(int i=0;i<=number_of_active_indices;i++)output_stream<<indices[i]<<", "<<x[i]<<std::endl;}
#endif
//#####################################################################
};
}
#endif

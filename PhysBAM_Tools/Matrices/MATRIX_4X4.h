//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_4X4
//#####################################################################
#ifndef __MATRIX_4X4__
#define __MATRIX_4X4__

#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T>
class MATRIX<T,4>:public MATRIX_BASE<T,MATRIX<T,4> >
{
public:
    typedef T SCALAR;typedef MATRIX_BASE<T,MATRIX<T,4> > BASE;
    using BASE::operator*;using BASE::Transpose_Times;using BASE::Times_Transpose;

    T x[16];

    explicit MATRIX(INITIAL_SIZE m=INITIAL_SIZE(4),INITIAL_SIZE n=INITIAL_SIZE(4))
    {
        assert(m==INITIAL_SIZE(4) && n==INITIAL_SIZE(4));for(int i=0;i<16;i++) x[i]=T();
    }

    MATRIX(const MATRIX& matrix_input)
    {
        for(int i=0;i<16;i++) x[i]=matrix_input.x[i];
    }

    MATRIX(const T x11,const T x21,const T x31,const T x41,const T x12,const T x22,const T x32,const T x42,const T x13,const T x23,const T x33,const T x43,const T x14,const T x24,const T x34,
        const T x44)
    {
        x[0]=x11;x[1]=x21;x[2]=x31;x[3]=x41;x[4]=x12;x[5]=x22;x[6]=x32;x[7]=x42;x[8]=x13;x[9]=x23;x[10]=x33;x[11]=x43;x[12]=x14;x[13]=x24;x[14]=x34;x[15]=x44;
    }

    template<class T_MATRIX>
    explicit MATRIX(const MATRIX_BASE<T,T_MATRIX>& A)
    {
        assert(A.Rows()==4 && A.Columns()==4);for(int j=1;j<=4;j++) for(int i=1;i<=4;i++) (*this)(i,j)=A(i,j);
    }

    template<class T2> explicit
    MATRIX(const MATRIX<T2,4>& matrix_input)
    {
        for(int i=0;i<16;i++) x[i]=(T)matrix_input.x[i];
    }

    MATRIX& operator=(const MATRIX& matrix_input)
    {
        for(int i=0;i<16;i++) x[i]=matrix_input.x[i];return *this;
    }

    int Rows() const
    {return 4;}

    int Columns() const
    {return 4;}

    T& operator()(const int i,const int j)
    {assert(i>=1 && i<=4);assert(j>=1 && j<=4);return x[i-1+4*(j-1)];}

    const T& operator()(const int i,const int j) const
    {assert(i>=1 && i<=4);assert(j>=1 && j<=4);return x[i-1+4*(j-1)];}

    bool Valid_Index(const int i,const int j) const
    {return 1<=i && i<=4 && 1<=j && j<=4;}

    VECTOR<T,4>& Column(const int j)
    {assert(1<=j && j<=4);return *(VECTOR<T,4>*)(x+4*(j-1));}

    const VECTOR<T,4>& Column(const int j) const
    {assert(1<=j && j<=4);return *(const VECTOR<T,4>*)(x+4*(j-1));}

    bool operator==(const MATRIX& A) const
    {for(int i=0;i<16;i++) if(x[i]!=A.x[i]) return false;return true;}

    bool operator!=(const MATRIX& A) const
    {return !(*this==A);}

    MATRIX operator-() const
    {return MATRIX(-x[0],-x[1],-x[2],-x[3],-x[4],-x[5],-x[6],-x[7],-x[8],-x[9],-x[10],-x[11],-x[12],-x[13],-x[14],-x[15]);}

    MATRIX& operator+=(const MATRIX& A)
    {for(int i=0;i<16;i++) x[i]+=A.x[i];return *this;}

    MATRIX& operator-=(const MATRIX& A)
    {for(int i=0;i<16;i++) x[i]-=A.x[i];return *this;}

    MATRIX& operator*=(const MATRIX& A)
    {return *this=*this*A;}

    MATRIX& operator*=(const T a)
    {for(int i=0;i<16;i++) x[i]*=a;return *this;}

    MATRIX& operator/=(const T a)
    {assert(a!=0);T s=1/a;for(int i=0;i<16;i++) x[i]*=s;return *this;}

    MATRIX operator+(const MATRIX& A) const
    {return MATRIX(x[0]+A.x[0],x[1]+A.x[1],x[2]+A.x[2],x[3]+A.x[3],x[4]+A.x[4],x[5]+A.x[5],x[6]+A.x[6],x[7]+A.x[7],
                   x[8]+A.x[8],x[9]+A.x[9],x[10]+A.x[10],x[11]+A.x[11],x[12]+A.x[12],x[13]+A.x[13],x[14]+A.x[14],x[15]+A.x[15]);}

    MATRIX operator-(const MATRIX& A) const
    {return MATRIX(x[0]-A.x[0],x[1]-A.x[1],x[2]-A.x[2],x[3]-A.x[3],x[4]-A.x[4],x[5]-A.x[5],x[6]-A.x[6],x[7]-A.x[7],
                   x[8]-A.x[8],x[9]-A.x[9],x[10]-A.x[10],x[11]-A.x[11],x[12]-A.x[12],x[13]-A.x[13],x[14]-A.x[14],x[15]-A.x[15]);}

    VECTOR<T,4> operator*(const VECTOR<T,4>& v) const
    {return VECTOR<T,4>(x[0]*v(1)+x[4]*v(2)+x[8]*v(3)+x[12]*v(4),x[1]*v(1)+x[5]*v(2)+x[9]*v(3)+x[13]*v(4),
                   x[2]*v(1)+x[6]*v(2)+x[10]*v(3)+x[14]*v(4),x[3]*v(1)+x[7]*v(2)+x[11]*v(3)+x[15]*v(4));}

    MATRIX operator*(const T a) const
    {return MATRIX(a*x[0],a*x[1],a*x[2],a*x[3],a*x[4],a*x[5],a*x[6],a*x[7],a*x[8],a*x[9],a*x[10],a*x[11],a*x[12],a*x[13],a*x[14],a*x[15]);}

    MATRIX operator/(const T a) const
    {assert(a!=0);T s=1/a;
    return MATRIX(s*x[0],s*x[1],s*x[2],s*x[3],s*x[4],s*x[5],s*x[6],s*x[7],s*x[8],s*x[9],s*x[10],s*x[11],s*x[12],s*x[13],s*x[14],s*x[15]);}

    VECTOR<T,3> Homogeneous_Times(const VECTOR<T,3>& v) const // assumes w=1 is the 4th coordinate of v
    {T w=x[3]*v.x+x[7]*v.y+x[11]*v.z+x[15];assert(w!=0);
    if(w==1) return VECTOR<T,3>(x[0]*v.x+x[4]*v.y+x[8]*v.z+x[12],x[1]*v.x+x[5]*v.y+x[9]*v.z+x[13],x[2]*v.x+x[6]*v.y+x[10]*v.z+x[14]);
    T s=1/w;// rescale so w=1
    return VECTOR<T,3>(s*(x[0]*v.x+x[4]*v.y+x[8]*v.z+x[12]),s*(x[1]*v.x+x[5]*v.y+x[9]*v.z+x[13]),s*(x[2]*v.x+x[6]*v.y+x[10]*v.z+x[14]));}

    VECTOR<T,3> Transform_3X3(const VECTOR<T,3>& v) const // multiplies vector by upper 3x3 of matrix only
    {return VECTOR<T,3>(x[0]*v.x+x[4]*v.y+x[8]*v.z,x[1]*v.x+x[5]*v.y+x[9]*v.z,x[2]*v.x+x[6]*v.y+x[10]*v.z);}

    void Invert()
    {*this=Inverse();}

    void Invert_Rotation_And_Translation()
    {*this=Inverse_Rotation_And_Translation();}

    MATRIX Inverse_Rotation_And_Translation()
    {return MATRIX(x[0],x[4],x[8],0,x[1],x[5],x[9],0,x[2],x[6],x[10],0,-x[0]*x[12]-x[1]*x[13]-x[2]*x[14],-x[4]*x[12]-x[5]*x[13]-x[6]*x[14],-x[8]*x[12]-x[9]*x[13]-x[10]*x[14],1);}

    MATRIX Rotation_Only() const
    {return MATRIX(x[0],x[1],x[2],0,x[4],x[5],x[6],0,x[8],x[9],x[10],0,0,0,0,1);}

    const VECTOR<T,3>& Translation() const
    {return *(const VECTOR<T,3>*)(x+12);} // x[12],x[13],x[14]

    VECTOR<T,3>& Translation()
    {return *(VECTOR<T,3>*)(x+12);} // x[12],x[13],x[14]

    MATRIX<T,3> Extract_Rotation() const
    {return MATRIX<T,3>(x[0],x[1],x[2],x[4],x[5],x[6],x[8],x[9],x[10]);}

    static MATRIX From_Linear(const MATRIX<T,3>& M) // Create a homogeneous 4x4 matrix corresponding to a 3x3 transform
    {return MATRIX(M.x[0],M.x[1],M.x[2],0,M.x[3],M.x[4],M.x[5],0,M.x[6],M.x[7],M.x[8],0,0,0,0,1);}

    void Transpose()
    {exchange(x[1],x[4]);exchange(x[2],x[8]);exchange(x[3],x[12]);exchange(x[6],x[9]);exchange(x[7],x[13]);exchange(x[11],x[14]);}

    MATRIX Transposed() const
    {return MATRIX(x[0],x[4],x[8],x[12],x[1],x[5],x[9],x[13],x[2],x[6],x[10],x[14],x[3],x[7],x[11],x[15]);}

    static MATRIX Translation_Matrix(const VECTOR<T,3>& translation)
    {return MATRIX(1,0,0,0,0,1,0,0,0,0,1,0,translation.x,translation.y,translation.z,1);}

    static MATRIX Identity_Matrix()
    {return MATRIX(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);}

    static MATRIX Rotation_Matrix_X_Axis(const T radians)
    {T c=cos(radians),s=sin(radians);return MATRIX(1,0,0,0,0,c,s,0,0,-s,c,0,0,0,0,1);}

    static MATRIX Rotation_Matrix_Y_Axis(const T radians)
    {T c=cos(radians),s=sin(radians);return MATRIX(c,0,-s,0,0,1,0,0,s,0,c,0,0,0,0,1);}

    static MATRIX Rotation_Matrix_Z_Axis(const T radians)
    {T c=cos(radians),s=sin(radians);return MATRIX(c,s,0,0,-s,c,0,0,0,0,1,0,0,0,0,1);}

    static MATRIX Rotation_Matrix(const VECTOR<T,3>& axis,const T radians)
    {return From_Linear(MATRIX<T,3>::Rotation_Matrix(axis,radians));}

    static MATRIX Rotation_Matrix(const VECTOR<T,3>& x_final,const VECTOR<T,3>& y_final,const VECTOR<T,3>& z_final)
    {return MATRIX(x_final.x,x_final.y,x_final.z,0,y_final.x,y_final.y,y_final.z,0,z_final.x,z_final.y,z_final.z,0,0,0,0,1);}

    static MATRIX Scale_Matrix(const VECTOR<T,3>& scale_vector)
    {return MATRIX(scale_vector.x,0,0,0,0,scale_vector.y,0,0,0,0,scale_vector.z,0,0,0,0,1);}

    static MATRIX Scale_Matrix(const T scale)
    {return MATRIX(scale,0,0,0,0,scale,0,0,0,0,scale,0,0,0,0,1);}

    template<class T_MATRIX>
    void Get_Submatrix(const int istart,const int jstart,T_MATRIX& a) const
    {for(int i=1;i<=a.Rows();i++) for(int j=1;j<=a.Columns();j++) a(i,j)=(*this)(istart+i-1,jstart+j-1);}

    template<class T_MATRIX>
    void Set_Submatrix(const int istart,const int jstart,const T_MATRIX& a)
    {for(int i=1;i<=a.Rows();i++) for(int j=1;j<=a.Columns();j++) (*this)(istart+i-1,jstart+j-1)=a(i,j);}

//#####################################################################
    MATRIX operator*(const MATRIX& A) const;
    MATRIX Inverse() const;
    MATRIX Cofactor_Matrix() const;
    T Frobenius_Norm_Squared() const;
};
template<class T>
inline MATRIX<T,4> operator*(const T a,const MATRIX<T,4>& A)
{return MATRIX<T,4>(a*A.x[0],a*A.x[1],a*A.x[2],a*A.x[3],a*A.x[4],a*A.x[5],a*A.x[6],a*A.x[7],a*A.x[8],a*A.x[9],a*A.x[10],a*A.x[11],a*A.x[12],a*A.x[13],a*A.x[14],a*A.x[15]);}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
template<class T>
inline std::istream& operator>>(std::istream& input,MATRIX<T,4>& A)
{FILE_UTILITIES::Ignore(input,'[');for(int i=0;i<4;i++){for(int j=0;j<4;j++) input>>A.x[i+j*4];FILE_UTILITIES::Ignore(input,';');}FILE_UTILITIES::Ignore(input,']');return input;}
#endif
}
#endif

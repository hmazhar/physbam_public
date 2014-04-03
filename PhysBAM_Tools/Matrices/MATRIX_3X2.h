//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_3X2
//#####################################################################
#ifndef __MATRIX_3X2__
#define __MATRIX_3X2__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/exchange.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Tools/Matrices/MATRIX_BASE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T>
class MATRIX<T,3,2>:public MATRIX_BASE<T,MATRIX<T,3,2> >
{
public:
    typedef T SCALAR;typedef MATRIX_BASE<T,MATRIX<T,3,2> > BASE;
    enum WORKAROUND1 {m=3,n=2};
    using BASE::operator*;using BASE::Times_Transpose;using BASE::Transpose_Times;

    T x[6];

    MATRIX()
    {
        for(int i=0;i<6;i++) x[i]=0;
    }

    explicit MATRIX(INITIAL_SIZE mm,INITIAL_SIZE nn)
    {
        assert(mm==INITIAL_SIZE(3) && nn==INITIAL_SIZE(2));for(int i=0;i<6;i++) x[i]=0;
    }

    MATRIX(const MATRIX& matrix_input)
    {
        for(int i=0;i<6;i++) x[i]=matrix_input.x[i];
    }

    template<class T2> explicit
    MATRIX(const MATRIX<T2,3,2>& matrix_input)
    {
        for(int i=0;i<6;i++) x[i]=(T)matrix_input.x[i];
    }

    MATRIX(const T x11,const T x21,const T x31,const T x12,const T x22,const T x32)
    {
        x[0]=x11;x[1]=x21;x[2]=x31;x[3]=x12;x[4]=x22;x[5]=x32;
    }

    MATRIX(const VECTOR<T,3>& column1,const VECTOR<T,3>& column2)
    {
        x[0]=column1.x;x[1]=column1.y;x[2]=column1.z;x[3]=column2.x;x[4]=column2.y;x[5]=column2.z;
    }

    template<class T_MATRIX>
    explicit MATRIX(const MATRIX_BASE<T,T_MATRIX>& A)
    {
        assert(A.Rows()==3 && A.Columns()==2);for(int j=1;j<=2;j++) for(int i=1;i<=3;i++) (*this)(i,j)=A(i,j);
    }

    MATRIX& operator=(const MATRIX& matrix_input)
    {
        for(int i=0;i<6;i++) x[i]=matrix_input.x[i];return *this;
    }

    int Rows() const
    {return 3;}

    int Columns() const
    {return 2;}

    T& operator()(const int i,const int j)
    {assert(i>=1 && i<=3);assert(j>=1 && j<=2);return x[i-1+3*(j-1)];}

    const T& operator()(const int i,const int j) const
    {assert(i>=1 && i<=3);assert(j>=1 && j<=2);return x[i-1+3*(j-1)];}

    bool Valid_Index(const int i,const int j) const
    {return 1<=i && i<=3 && 1<=j && j<=2;}

    VECTOR<T,3>& Column(const int j)
    {assert(1<=j && j<=2);return *(VECTOR<T,3>*)(x+3*(j-1));}

    const VECTOR<T,3>& Column(const int j) const
    {assert(1<=j && j<=2);return *(const VECTOR<T,3>*)(x+3*(j-1));}

    bool operator==(const MATRIX& A) const
    {for(int i=0;i<6;i++) if(x[i]!=A.x[i]) return false;return true;}

    bool operator!=(const MATRIX& A) const
    {return !(*this==A);}

    VECTOR<T,3> Column_Sum() const
    {return VECTOR<T,3>(x[0]+x[3],x[1]+x[4],x[2]+x[5]);}

    MATRIX operator-() const
    {return MATRIX(-x[0],-x[1],-x[2],-x[3],-x[4],-x[5]);}

    MATRIX& operator+=(const MATRIX& A)
    {for(int i=0;i<6;i++) x[i]+=A.x[i];return *this;}

    MATRIX& operator-=(const MATRIX& A)
    {for(int i=0;i<6;i++) x[i]-=A.x[i];return *this;}

    MATRIX& operator*=(const MATRIX<T,2>& A)
    {return *this=*this*A;}

    MATRIX& operator*=(const T a)
    {for(int i=0;i<6;i++) x[i]*=a;return *this;}

    MATRIX& operator/=(const T a)
    {assert(a!=0);T s=1/a;for(int i=0;i<6;i++) x[i]*=s;return *this;}

    MATRIX operator+(const MATRIX& A) const // 6 adds
    {return MATRIX(x[0]+A.x[0],x[1]+A.x[1],x[2]+A.x[2],x[3]+A.x[3],x[4]+A.x[4],x[5]+A.x[5]);}

    MATRIX operator-(const MATRIX& A) const // 6 adds
    {return MATRIX(x[0]-A.x[0],x[1]-A.x[1],x[2]-A.x[2],x[3]-A.x[3],x[4]-A.x[4],x[5]-A.x[5]);}

    MATRIX operator*(const MATRIX<T,2>& A) const // 12 mults, 6 adds
    {return MATRIX(x[0]*A.x[0]+x[3]*A.x[1],x[1]*A.x[0]+x[4]*A.x[1],x[2]*A.x[0]+x[5]*A.x[1],x[0]*A.x[2]+x[3]*A.x[3],x[1]*A.x[2]+x[4]*A.x[3],x[2]*A.x[2]+x[5]*A.x[3]);}

    MATRIX_MXN<T> operator*(const MATRIX_MXN<T>& A) const
    {assert(2==A.m);MATRIX_MXN<T> matrix(3,A.n);for(int j=1;j<=A.n;j++) for(int i=1;i<=3;i++) for(int k=1;k<=2;k++) matrix(i,j)+=(*this)(i,k)*A(k,j);return matrix;}

    MATRIX Times_Transpose(const UPPER_TRIANGULAR_MATRIX<T,2>& A) const // 9 mults, 3 adds
    {return MATRIX(x[0]*A.x11+x[3]*A.x12,x[1]*A.x11+x[4]*A.x12,x[2]*A.x11+x[5]*A.x12,x[3]*A.x22,x[4]*A.x22,x[5]*A.x22);}

    MATRIX Times_Transpose(const MATRIX<T,2>& A) const // 12 mults, 6 adds
    {return MATRIX(x[0]*A.x[0]+x[3]*A.x[2],x[1]*A.x[0]+x[4]*A.x[2],x[2]*A.x[0]+x[5]*A.x[2],x[0]*A.x[1]+x[3]*A.x[3],x[1]*A.x[1]+x[4]*A.x[3],x[2]*A.x[1]+x[5]*A.x[3]);}

    MATRIX<T,3> Times_Transpose(const MATRIX& A) const // 18 mults, 9 adds
    {return MATRIX<T,3>(x[0]*A.x[0]+x[3]*A.x[3],x[1]*A.x[0]+x[4]*A.x[3],x[2]*A.x[0]+x[5]*A.x[3],
                        x[0]*A.x[1]+x[3]*A.x[4],x[1]*A.x[1]+x[4]*A.x[4],x[2]*A.x[1]+x[5]*A.x[4],
                        x[0]*A.x[2]+x[3]*A.x[5],x[1]*A.x[2]+x[4]*A.x[5],x[2]*A.x[2]+x[5]*A.x[5]);}

    MATRIX operator*(const T a) const // 6 mults
    {return MATRIX(a*x[0],a*x[1],a*x[2],a*x[3],a*x[4],a*x[5]);}

    MATRIX operator/(const T a) const // 6 mults, 1 div
    {assert(a!=0);T s=1/a;return MATRIX(s*x[0],s*x[1],s*x[2],s*x[3],s*x[4],s*x[5]);}

    VECTOR<T,3> operator*(const VECTOR<T,2>& v) const // 6 mults, 3 adds
    {return VECTOR<T,3>(x[0]*v.x+x[3]*v.y,x[1]*v.x+x[4]*v.y,x[2]*v.x+x[5]*v.y);}

    UPPER_TRIANGULAR_MATRIX<T,2> R_From_QR_Factorization() const // Gram Schmidt
    {T x_dot_x=Column(1).Magnitude_Squared(),x_dot_y=VECTOR<T,3>::Dot_Product(Column(1),Column(2)),y_dot_y=Column(2).Magnitude_Squared();
    T r11=sqrt(x_dot_x),r12=r11?x_dot_y/r11:0,r22=sqrt(max((T)0,y_dot_y-r12*r12));
    return UPPER_TRIANGULAR_MATRIX<T,2>(r11,r12,r22);}

    SYMMETRIC_MATRIX<T,2> Normal_Equations_Matrix() const // 9 mults, 6 adds
    {return SYMMETRIC_MATRIX<T,2>(x[0]*x[0]+x[1]*x[1]+x[2]*x[2],x[3]*x[0]+x[4]*x[1]+x[5]*x[2],x[3]*x[3]+x[4]*x[4]+x[5]*x[5]);}

    MATRIX<T,2> Transpose_Times(const MATRIX& A) const // 12 mults, 8 adds
    {return MATRIX<T,2>(x[0]*A.x[0]+x[1]*A.x[1]+x[2]*A.x[2],x[3]*A.x[0]+x[4]*A.x[1]+x[5]*A.x[2],x[0]*A.x[3]+x[1]*A.x[4]+x[2]*A.x[5],x[3]*A.x[3]+x[4]*A.x[4]+x[5]*A.x[5]);}

    MATRIX<T,2,3> Transpose_Times(const MATRIX<T,3>& A) const
    {return MATRIX<T,2,3>(x[0]*A.x[0]+x[1]*A.x[1]+x[2]*A.x[2],x[3]*A.x[0]+x[4]*A.x[1]+x[5]*A.x[2],
                          x[0]*A.x[3]+x[1]*A.x[4]+x[2]*A.x[5],x[3]*A.x[3]+x[4]*A.x[4]+x[5]*A.x[5],
                          x[0]*A.x[6]+x[1]*A.x[7]+x[2]*A.x[8],x[3]*A.x[6]+x[4]*A.x[7]+x[5]*A.x[8]);}

    MATRIX<T,2,3> Transpose_Times(const SYMMETRIC_MATRIX<T,3>& A) const
    {return MATRIX<T,2,3>(x[0]*A.x11+x[1]*A.x21+x[2]*A.x31,x[3]*A.x11+x[4]*A.x21+x[5]*A.x31,
                          x[0]*A.x21+x[1]*A.x22+x[2]*A.x32,x[3]*A.x21+x[4]*A.x22+x[5]*A.x32,
                          x[0]*A.x31+x[1]*A.x32+x[2]*A.x33,x[3]*A.x31+x[4]*A.x32+x[5]*A.x33);}

    MATRIX<T,2,3> Transpose_Times(const UPPER_TRIANGULAR_MATRIX<T,3>& A) const
    {return MATRIX<T,2,3>(x[0]*A.x11,x[3]*A.x11,x[0]*A.x12+x[1]*A.x22,x[3]*A.x12+x[4]*A.x22,
                          x[0]*A.x13+x[1]*A.x23+x[2]*A.x33,x[3]*A.x13+x[4]*A.x23+x[5]*A.x33);}

    MATRIX<T,2,3> Transpose_Times(const DIAGONAL_MATRIX<T,3>& A) const
    {return MATRIX<T,2,3>(x[0]*A.x11,x[3]*A.x11,x[1]*A.x22,x[4]*A.x22,x[2]*A.x33,x[5]*A.x33);}

    VECTOR<T,2> Transpose_Times(const VECTOR<T,3>& v) const // 6 mults, 4 adds
    {return VECTOR<T,2>(x[0]*v.x+x[1]*v.y+x[2]*v.z,x[3]*v.x+x[4]*v.y+x[5]*v.z);}

    T Max_Abs() const
    {return maxabs(x[0],x[1],x[2],x[3],x[4],x[5]);}

    T Frobenius_Norm_Squared() const
    {return sqr(x[0])+sqr(x[1])+sqr(x[2])+sqr(x[3])+sqr(x[4])+sqr(x[5]);}

    MATRIX operator*(const UPPER_TRIANGULAR_MATRIX<T,2>& A) const // 9 mults, 3 adds
    {return MATRIX(x[0]*A.x11,x[1]*A.x11,x[2]*A.x11,x[0]*A.x12+x[3]*A.x22,x[1]*A.x12+x[4]*A.x22,x[2]*A.x12+x[5]*A.x22);}

    MATRIX operator*(const SYMMETRIC_MATRIX<T,2>& A) const // 12 mults, 6 adds
    {return MATRIX(x[0]*A.x11+x[3]*A.x21,x[1]*A.x11+x[4]*A.x21,x[2]*A.x11+x[5]*A.x21,x[0]*A.x21+x[3]*A.x22,x[1]*A.x21+x[4]*A.x22,x[2]*A.x21+x[5]*A.x22);}

    MATRIX operator*(const DIAGONAL_MATRIX<T,2>& A) const // 6 mults
    {return MATRIX(x[0]*A.x11,x[1]*A.x11,x[2]*A.x11,x[3]*A.x22,x[4]*A.x22,x[5]*A.x22);}

    static T Inner_Product(const MATRIX& A,const MATRIX& B) // 6 mults, 5 adds
    {return A.x[0]*B.x[0]+A.x[1]*B.x[1]+A.x[2]*B.x[2]+A.x[3]*B.x[3]+A.x[4]*B.x[4]+A.x[5]*B.x[5];}

    VECTOR<T,3> Weighted_Normal() const
    {return VECTOR<T,3>::Cross_Product(Column(1),Column(2));}

    MATRIX Cofactor_Matrix() const
    {VECTOR<T,3> normal=Weighted_Normal().Normalized();
    return MATRIX(VECTOR<T,3>::Cross_Product(Column(2),normal),VECTOR<T,3>::Cross_Product(normal,Column(1)));}

    T Parallelepiped_Measure() const
    {return Weighted_Normal().Magnitude();}

    MATRIX<T,2,3> Transposed() const
    {return MATRIX<T,2,3>::Transposed(*this);}

    static MATRIX Outer_Product(const VECTOR<T,3>& u,const VECTOR<T,2>& v)
    {MATRIX result;for(int i=1;i<=3;i++) for(int j=1;j<=2;j++) result(i,j)=u(i)*v(j);return result;}

//#####################################################################
    void Fast_Singular_Value_Decomposition(MATRIX& U,DIAGONAL_MATRIX<T,2>& singular_values,MATRIX<T,2>& V) const;
    void Fast_Indefinite_Polar_Decomposition(MATRIX<T,3,2>& Q,SYMMETRIC_MATRIX<T,2>& S) const;
//#####################################################################
};
// global functions
template<class T>
inline MATRIX<T,3,2> operator*(const MATRIX<T,3>& B,const MATRIX<T,3,2>& A) // 18 mults, 12 adds
{return MATRIX<T,3,2>(B.x[0]*A.x[0]+B.x[3]*A.x[1]+B.x[6]*A.x[2],B.x[1]*A.x[0]+B.x[4]*A.x[1]+B.x[7]*A.x[2],B.x[2]*A.x[0]+B.x[5]*A.x[1]+B.x[8]*A.x[2],
    B.x[0]*A.x[3]+B.x[3]*A.x[4]+B.x[6]*A.x[5],B.x[1]*A.x[3]+B.x[4]*A.x[4]+B.x[7]*A.x[5],B.x[2]*A.x[3]+B.x[5]*A.x[4]+B.x[8]*A.x[5]);}

template<class T>
inline MATRIX<T,3,2> operator*(const T a,const MATRIX<T,3,2>& A) // 6 mults
{return MATRIX<T,3,2>(a*A.x[0],a*A.x[1],a*A.x[2],a*A.x[3],a*A.x[4],a*A.x[5]);}

template<class T>
inline VECTOR<T,3> operator*(const VECTOR<T,2>& v,const MATRIX<T,3,2>& A) // 6 mults, 3 adds
{return VECTOR<T,3>(v.x*A.x[0]+v.y*A.x[3],v.x*A.x[1]+v.y*A.x[4],v.x*A.x[2]+v.y*A.x[5]);}

template<class T>
inline MATRIX<T,3,2> operator*(const SYMMETRIC_MATRIX<T,3>& A,const MATRIX<T,3,2>& B) // 18 mults, 12 adds
{return MATRIX<T,3,2>(A.x11*B.x[0]+A.x21*B.x[1]+A.x31*B.x[2],A.x21*B.x[0]+A.x22*B.x[1]+A.x32*B.x[2],A.x31*B.x[0]+A.x32*B.x[1]+A.x33*B.x[2],
    A.x11*B.x[3]+A.x21*B.x[4]+A.x31*B.x[5],A.x21*B.x[3]+A.x22*B.x[4]+A.x32*B.x[5],A.x31*B.x[3]+A.x32*B.x[4]+A.x33*B.x[5]);}

template<class T>
inline MATRIX<T,3,2> operator*(const UPPER_TRIANGULAR_MATRIX<T,3>& A,const MATRIX<T,3,2>& B)
{return MATRIX<T,3,2>(A.x11*B.x[0]+A.x12*B.x[1]+A.x13*B.x[2],A.x22*B.x[1]+A.x23*B.x[2],A.x33*B.x[2],
    A.x11*B.x[3]+A.x12*B.x[4]+A.x13*B.x[5],A.x22*B.x[4]+A.x23*B.x[5],A.x33*B.x[5]);}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
template<class T>
inline std::istream& operator>>(std::istream& input,MATRIX<T,3,2>& A)
{FILE_UTILITIES::Ignore(input,'[');for(int i=0;i<3;i++){for(int j=0;j<2;j++) input>>A.x[i+j*3];FILE_UTILITIES::Ignore(input,';');}FILE_UTILITIES::Ignore(input,']');return input;}
#endif
}
#endif

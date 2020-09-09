#include "field_m0.hh"
#include<cmath>
void field_m0_P2M_1(double x, double y, double z, double q, double * M) {
M[0] += q;
M[1] += -q*x;
M[2] += -q*y;
M[3] += -q*z;

}
void field_m0_M2M_1(double x, double y, double z, double * M, double * Ms) {
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += x*M[0] + M[1];
#pragma omp atomic
Ms[2] += y*M[0] + M[2];
#pragma omp atomic
Ms[3] += z*M[0] + M[3];

}

void field_m0_M2L_1(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double D[4];
double Dtmp0 = 1.0*pow(R, -3.0);
D[0] = 1.0/R;
D[1] = -Dtmp0*x;
D[2] = -Dtmp0*y;
D[3] = -Dtmp0*z;
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3];
#pragma omp atomic
L[1] += D[1]*M[0];
#pragma omp atomic
L[2] += D[2]*M[0];
#pragma omp atomic
L[3] += D[3]*M[0];

}

void field_m0_L2L_1(double x, double y, double z, double * L, double * Ls) {
#pragma omp atomic
Ls[0] += x*L[1] + y*L[2] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += L[1];
#pragma omp atomic
Ls[2] += L[2];
#pragma omp atomic
Ls[3] += L[3];

}

void field_m0_L2P_1(double x, double y, double z, double * L, double * F) {
#pragma omp atomic
F[0] += -L[1];
#pragma omp atomic
F[1] += -L[2];
#pragma omp atomic
F[2] += -L[3];

}

void field_m0_M2P_1(double x, double y, double z, double * M, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = 1.0*pow(R, -3.0);
double Ftmp1 = Ftmp0*M[0];
double Ftmp2 = 3.0*pow(R, -5.0);
double Ftmp3 = Ftmp2*M[2];
double Ftmp4 = x*y;
double Ftmp5 = Ftmp2*M[3];
double Ftmp6 = Ftmp5*z;
double Ftmp7 = Ftmp2*M[1];
#pragma omp atomic
F[0] += Ftmp0*M[1] + Ftmp1*x - Ftmp3*Ftmp4 - Ftmp6*x - Ftmp7*(x*x);
#pragma omp atomic
F[1] += Ftmp0*M[2] + Ftmp1*y - Ftmp3*(y*y) - Ftmp4*Ftmp7 - Ftmp6*y;
#pragma omp atomic
F[2] += Ftmp0*M[3] + Ftmp1*z - Ftmp3*y*z - Ftmp5*(z*z) - Ftmp7*x*z;

}

template<>
void P2P<0, 3>(double x, double y, double z, double * S, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = 1.0*pow(R, -3.0)*S[0];
#pragma omp atomic
F[0] += Ftmp0*x;
#pragma omp atomic
F[1] += Ftmp0*y;
#pragma omp atomic
F[2] += Ftmp0*z;

}

void field_m0_P2M_2(double x, double y, double z, double q, double * M) {
double Mtmp0 = q*x;
double Mtmp1 = q*y;
double Mtmp2 = (1.0/2.0)*q;
M[0] += q;
M[1] += -Mtmp0;
M[2] += -Mtmp1;
M[3] += -q*z;
M[4] += Mtmp2*(x*x);
M[5] += Mtmp0*y;
M[6] += Mtmp0*z;
M[7] += Mtmp2*(y*y);
M[8] += Mtmp1*z;
M[9] += Mtmp2*(z*z);

}
void field_m0_M2M_2(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = y*M[0];
double Mstmp2 = (1.0/2.0)*M[0];
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += Mstmp0 + M[1];
#pragma omp atomic
Ms[2] += Mstmp1 + M[2];
#pragma omp atomic
Ms[3] += z*M[0] + M[3];
#pragma omp atomic
Ms[4] += Mstmp2*(x*x) + x*M[1] + M[4];
#pragma omp atomic
Ms[5] += Mstmp0*y + x*M[2] + y*M[1] + M[5];
#pragma omp atomic
Ms[6] += Mstmp0*z + x*M[3] + z*M[1] + M[6];
#pragma omp atomic
Ms[7] += Mstmp2*(y*y) + y*M[2] + M[7];
#pragma omp atomic
Ms[8] += Mstmp1*z + y*M[3] + z*M[2] + M[8];
#pragma omp atomic
Ms[9] += Mstmp2*(z*z) + z*M[3] + M[9];

}

void field_m0_M2L_2(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double D[10];
double Dtmp0 = 1.0*pow(R, -3.0);
double Dtmp1 = -Dtmp0;
double Dtmp2 = 3.0*pow(R, -5.0);
double Dtmp3 = Dtmp2*x;
D[0] = 1.0/R;
D[1] = -Dtmp0*x;
D[2] = -Dtmp0*y;
D[3] = -Dtmp0*z;
D[4] = Dtmp1 + Dtmp2*(x*x);
D[5] = Dtmp3*y;
D[6] = Dtmp3*z;
D[7] = Dtmp1 + Dtmp2*(y*y);
D[8] = Dtmp2*y*z;
D[9] = -D[4] - D[7];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[8]*M[8] + D[9]*M[9];
#pragma omp atomic
L[1] += D[1]*M[0] + D[4]*M[1] + D[5]*M[2] + D[6]*M[3];
#pragma omp atomic
L[2] += D[2]*M[0] + D[5]*M[1] + D[7]*M[2] + D[8]*M[3];
#pragma omp atomic
L[3] += D[3]*M[0] + D[6]*M[1] + D[8]*M[2] + D[9]*M[3];
#pragma omp atomic
L[4] += D[4]*M[0];
#pragma omp atomic
L[5] += D[5]*M[0];
#pragma omp atomic
L[6] += D[6]*M[0];
#pragma omp atomic
L[7] += D[7]*M[0];
#pragma omp atomic
L[8] += D[8]*M[0];
#pragma omp atomic
L[9] += D[9]*M[0];

}

void field_m0_L2L_2(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x + Lstmp2*y + (1.0/2.0)*(x*x)*L[4] + x*L[1] + (1.0/2.0)*(y*y)*L[7] + y*L[2] + (1.0/2.0)*(z*z)*L[9] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp2 + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += x*L[6] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += L[4];
#pragma omp atomic
Ls[5] += L[5];
#pragma omp atomic
Ls[6] += L[6];
#pragma omp atomic
Ls[7] += L[7];
#pragma omp atomic
Ls[8] += L[8];
#pragma omp atomic
Ls[9] += L[9];

}

void field_m0_L2P_2(double x, double y, double z, double * L, double * F) {
#pragma omp atomic
F[0] += -x*L[4] - y*L[5] - z*L[6] - L[1];
#pragma omp atomic
F[1] += -x*L[5] - y*L[7] - z*L[8] - L[2];
#pragma omp atomic
F[2] += -x*L[6] - y*L[8] - z*L[9] - L[3];

}

void field_m0_M2P_2(double x, double y, double z, double * M, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = 1.0*pow(R, -3.0);
double Ftmp1 = pow(R, -5.0);
double Ftmp2 = 3.0*Ftmp1;
double Ftmp3 = Ftmp2*M[5];
double Ftmp4 = Ftmp2*z;
double Ftmp5 = Ftmp0*M[0];
double Ftmp6 = Ftmp2*M[2];
double Ftmp7 = x*y;
double Ftmp8 = Ftmp4*M[3];
double Ftmp9 = (x*x);
double Ftmp10 = Ftmp2*M[1];
double Ftmp11 = y*M[8];
double Ftmp12 = 15.0*pow(R, -7.0);
double Ftmp13 = Ftmp12*z;
double Ftmp14 = Ftmp13*x;
double Ftmp15 = Ftmp12*Ftmp9;
double Ftmp16 = y*M[5];
double Ftmp17 = -9.0*Ftmp1;
double Ftmp18 = -Ftmp2;
double Ftmp19 = (y*y);
double Ftmp20 = Ftmp12*Ftmp19;
double Ftmp21 = (Ftmp18 + Ftmp20)*M[7];
double Ftmp22 = (z*z);
double Ftmp23 = Ftmp12*Ftmp22;
double Ftmp24 = (Ftmp18 + Ftmp23)*M[9];
double Ftmp25 = x*M[6];
double Ftmp26 = (Ftmp15 + Ftmp18)*M[4];
#pragma omp atomic
F[0] += Ftmp0*M[1] - Ftmp10*Ftmp9 + Ftmp11*Ftmp14 + Ftmp15*Ftmp16 + Ftmp15*z*M[6] + Ftmp21*x + Ftmp24*x - Ftmp3*y - Ftmp4*M[6] + Ftmp5*x - Ftmp6*Ftmp7 - Ftmp8*x + x*(Ftmp15 + Ftmp17)*M[4];
#pragma omp atomic
F[1] += Ftmp0*M[2] - Ftmp10*Ftmp7 + Ftmp13*Ftmp25*y - Ftmp19*Ftmp6 + Ftmp20*x*M[5] + Ftmp20*z*M[8] + Ftmp24*y + Ftmp26*y - Ftmp3*x - Ftmp4*M[8] + Ftmp5*y - Ftmp8*y + y*(Ftmp17 + Ftmp20)*M[7];
#pragma omp atomic
F[2] += Ftmp0*M[3] - Ftmp11*Ftmp2 + Ftmp11*Ftmp23 + Ftmp14*Ftmp16 - Ftmp2*Ftmp22*M[3] - Ftmp2*Ftmp25 + Ftmp21*z + Ftmp23*Ftmp25 + Ftmp26*z - Ftmp4*x*M[1] - Ftmp4*y*M[2] + Ftmp5*z + z*(Ftmp17 + Ftmp23)*M[9];

}

void field_m0_P2M_3(double x, double y, double z, double q, double * M) {
double Mtmp0 = q*x;
double Mtmp1 = q*y;
double Mtmp2 = q*z;
double Mtmp3 = (x*x);
double Mtmp4 = (1.0/2.0)*q;
double Mtmp5 = Mtmp0*y;
double Mtmp6 = (y*y);
double Mtmp7 = (z*z);
double Mtmp8 = (1.0/6.0)*q;
double Mtmp9 = (1.0/2.0)*Mtmp3;
double Mtmp10 = (1.0/2.0)*Mtmp0;
M[0] += q;
M[1] += -Mtmp0;
M[2] += -Mtmp1;
M[3] += -Mtmp2;
M[4] += Mtmp3*Mtmp4;
M[5] += Mtmp5;
M[6] += Mtmp0*z;
M[7] += Mtmp4*Mtmp6;
M[8] += Mtmp1*z;
M[9] += Mtmp4*Mtmp7;
M[10] += -Mtmp8*(x*x*x);
M[11] += -Mtmp1*Mtmp9;
M[12] += -Mtmp2*Mtmp9;
M[13] += -Mtmp10*Mtmp6;
M[14] += -Mtmp5*z;
M[15] += -Mtmp10*Mtmp7;
M[16] += -Mtmp8*(y*y*y);
M[17] += -1.0/2.0*Mtmp2*Mtmp6;
M[18] += -1.0/2.0*Mtmp1*Mtmp7;
M[19] += -Mtmp8*(z*z*z);

}
void field_m0_M2M_3(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = y*M[0];
double Mstmp2 = z*M[0];
double Mstmp3 = x*M[1];
double Mstmp4 = (x*x);
double Mstmp5 = (1.0/2.0)*M[0];
double Mstmp6 = x*M[2];
double Mstmp7 = y*M[1];
double Mstmp8 = Mstmp0*y;
double Mstmp9 = x*M[3];
double Mstmp10 = y*M[2];
double Mstmp11 = (y*y);
double Mstmp12 = y*M[3];
double Mstmp13 = (z*z);
double Mstmp14 = (1.0/2.0)*Mstmp4;
double Mstmp15 = (1.0/6.0)*M[0];
double Mstmp16 = (1.0/2.0)*M[1];
double Mstmp17 = (1.0/2.0)*Mstmp11;
double Mstmp18 = (1.0/2.0)*Mstmp13;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += Mstmp0 + M[1];
#pragma omp atomic
Ms[2] += Mstmp1 + M[2];
#pragma omp atomic
Ms[3] += Mstmp2 + M[3];
#pragma omp atomic
Ms[4] += Mstmp3 + Mstmp4*Mstmp5 + M[4];
#pragma omp atomic
Ms[5] += Mstmp6 + Mstmp7 + Mstmp8 + M[5];
#pragma omp atomic
Ms[6] += Mstmp0*z + Mstmp9 + z*M[1] + M[6];
#pragma omp atomic
Ms[7] += Mstmp10 + Mstmp11*Mstmp5 + M[7];
#pragma omp atomic
Ms[8] += Mstmp1*z + Mstmp12 + z*M[2] + M[8];
#pragma omp atomic
Ms[9] += Mstmp13*Mstmp5 + z*M[3] + M[9];
#pragma omp atomic
Ms[10] += Mstmp14*M[1] + Mstmp15*(x*x*x) + x*M[4] + M[10];
#pragma omp atomic
Ms[11] += Mstmp1*Mstmp14 + Mstmp14*M[2] + Mstmp3*y + x*M[5] + y*M[4] + M[11];
#pragma omp atomic
Ms[12] += Mstmp14*Mstmp2 + Mstmp14*M[3] + Mstmp3*z + x*M[6] + z*M[4] + M[12];
#pragma omp atomic
Ms[13] += Mstmp0*Mstmp17 + Mstmp11*Mstmp16 + Mstmp6*y + x*M[7] + y*M[5] + M[13];
#pragma omp atomic
Ms[14] += Mstmp6*z + Mstmp7*z + Mstmp8*z + Mstmp9*y + x*M[8] + y*M[6] + z*M[5] + M[14];
#pragma omp atomic
Ms[15] += Mstmp0*Mstmp18 + Mstmp13*Mstmp16 + Mstmp9*z + x*M[9] + z*M[6] + M[15];
#pragma omp atomic
Ms[16] += Mstmp15*(y*y*y) + Mstmp17*M[2] + y*M[7] + M[16];
#pragma omp atomic
Ms[17] += Mstmp10*z + Mstmp17*Mstmp2 + Mstmp17*M[3] + y*M[8] + z*M[7] + M[17];
#pragma omp atomic
Ms[18] += Mstmp1*Mstmp18 + Mstmp12*z + Mstmp18*M[2] + y*M[9] + z*M[8] + M[18];
#pragma omp atomic
Ms[19] += Mstmp15*(z*z*z) + Mstmp18*M[3] + z*M[9] + M[19];

}

void field_m0_M2L_3(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double D[20];
double Dtmp0 = 1.0*pow(R, -3.0);
double Dtmp1 = -Dtmp0;
double Dtmp2 = (x*x);
double Dtmp3 = pow(R, -5.0);
double Dtmp4 = 3.0*Dtmp3;
double Dtmp5 = Dtmp4*x;
double Dtmp6 = (y*y);
double Dtmp7 = y*z;
double Dtmp8 = -9.0*Dtmp3;
double Dtmp9 = 15.0*pow(R, -7.0);
double Dtmp10 = Dtmp2*Dtmp9;
double Dtmp11 = -Dtmp4;
double Dtmp12 = Dtmp10 + Dtmp11;
double Dtmp13 = Dtmp6*Dtmp9;
double Dtmp14 = Dtmp11 + Dtmp13;
D[0] = 1.0/R;
D[1] = -Dtmp0*x;
D[2] = -Dtmp0*y;
D[3] = -Dtmp0*z;
D[4] = Dtmp1 + Dtmp2*Dtmp4;
D[5] = Dtmp5*y;
D[6] = Dtmp5*z;
D[7] = Dtmp1 + Dtmp4*Dtmp6;
D[8] = Dtmp4*Dtmp7;
D[9] = -D[4] - D[7];
D[10] = -x*(Dtmp10 + Dtmp8);
D[11] = -Dtmp12*y;
D[12] = -Dtmp12*z;
D[13] = -1.0*Dtmp14*x;
D[14] = -Dtmp7*Dtmp9*x;
D[15] = -D[10] - D[13];
D[16] = -y*(Dtmp13 + Dtmp8);
D[17] = -Dtmp14*z;
D[18] = -D[11] - D[16];
D[19] = -D[12] - D[17];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[8]*M[8] + D[9]*M[9] + D[10]*M[10] + D[11]*M[11] + D[12]*M[12] + D[13]*M[13] + D[14]*M[14] + D[15]*M[15] + D[16]*M[16] + D[17]*M[17] + D[18]*M[18] + D[19]*M[19];
#pragma omp atomic
L[1] += D[1]*M[0] + D[4]*M[1] + D[5]*M[2] + D[6]*M[3] + D[10]*M[4] + D[11]*M[5] + D[12]*M[6] + D[13]*M[7] + D[14]*M[8] + D[15]*M[9];
#pragma omp atomic
L[2] += D[2]*M[0] + D[5]*M[1] + D[7]*M[2] + D[8]*M[3] + D[11]*M[4] + D[13]*M[5] + D[14]*M[6] + D[16]*M[7] + D[17]*M[8] + D[18]*M[9];
#pragma omp atomic
L[3] += D[3]*M[0] + D[6]*M[1] + D[8]*M[2] + D[9]*M[3] + D[12]*M[4] + D[14]*M[5] + D[15]*M[6] + D[17]*M[7] + D[18]*M[8] + D[19]*M[9];
#pragma omp atomic
L[4] += D[4]*M[0] + D[10]*M[1] + D[11]*M[2] + D[12]*M[3];
#pragma omp atomic
L[5] += D[5]*M[0] + D[11]*M[1] + D[13]*M[2] + D[14]*M[3];
#pragma omp atomic
L[6] += D[6]*M[0] + D[12]*M[1] + D[14]*M[2] + D[15]*M[3];
#pragma omp atomic
L[7] += D[7]*M[0] + D[13]*M[1] + D[16]*M[2] + D[17]*M[3];
#pragma omp atomic
L[8] += D[8]*M[0] + D[14]*M[1] + D[17]*M[2] + D[18]*M[3];
#pragma omp atomic
L[9] += D[9]*M[0] + D[15]*M[1] + D[18]*M[2] + D[19]*M[3];
#pragma omp atomic
L[10] += D[10]*M[0];
#pragma omp atomic
L[11] += D[11]*M[0];
#pragma omp atomic
L[12] += D[12]*M[0];
#pragma omp atomic
L[13] += D[13]*M[0];
#pragma omp atomic
L[14] += D[14]*M[0];
#pragma omp atomic
L[15] += D[15]*M[0];
#pragma omp atomic
L[16] += D[16]*M[0];
#pragma omp atomic
L[17] += D[17]*M[0];
#pragma omp atomic
L[18] += D[18]*M[0];
#pragma omp atomic
L[19] += D[19]*M[0];

}

void field_m0_L2L_3(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
double Lstmp3 = z*L[14];
double Lstmp4 = Lstmp3*y;
double Lstmp5 = (1.0/2.0)*(x*x);
double Lstmp6 = (1.0/2.0)*(y*y);
double Lstmp7 = (1.0/2.0)*(z*z);
double Lstmp8 = x*L[13];
double Lstmp9 = x*L[15];
double Lstmp10 = y*L[11];
double Lstmp11 = z*L[12];
double Lstmp12 = y*L[18];
double Lstmp13 = z*L[17];
double Lstmp14 = y*L[13];
double Lstmp15 = y*L[14];
double Lstmp16 = z*L[15];
double Lstmp17 = z*L[18];
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x + Lstmp10*Lstmp5 + Lstmp11*Lstmp5 + Lstmp12*Lstmp7 + Lstmp13*Lstmp6 + Lstmp2*y + Lstmp4*x + Lstmp5*L[4] + Lstmp6*Lstmp8 + Lstmp6*L[7] + Lstmp7*Lstmp9 + Lstmp7*L[9] + (1.0/6.0)*(x*x*x)*L[10] + x*L[1] + (1.0/6.0)*(y*y*y)*L[16] + y*L[2] + (1.0/6.0)*(z*z*z)*L[19] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 + Lstmp10*x + Lstmp11*x + Lstmp4 + Lstmp5*L[10] + Lstmp6*L[13] + Lstmp7*L[15] + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp13*y + Lstmp14*x + Lstmp2 + Lstmp3*x + Lstmp5*L[11] + Lstmp6*L[16] + Lstmp7*L[18] + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += Lstmp15*x + Lstmp16*x + Lstmp17*y + Lstmp5*L[12] + Lstmp6*L[17] + Lstmp7*L[19] + x*L[6] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += Lstmp10 + Lstmp11 + x*L[10] + L[4];
#pragma omp atomic
Ls[5] += Lstmp14 + Lstmp3 + x*L[11] + L[5];
#pragma omp atomic
Ls[6] += Lstmp15 + Lstmp16 + x*L[12] + L[6];
#pragma omp atomic
Ls[7] += Lstmp13 + Lstmp8 + y*L[16] + L[7];
#pragma omp atomic
Ls[8] += Lstmp17 + x*L[14] + y*L[17] + L[8];
#pragma omp atomic
Ls[9] += Lstmp12 + Lstmp9 + z*L[19] + L[9];
#pragma omp atomic
Ls[10] += L[10];
#pragma omp atomic
Ls[11] += L[11];
#pragma omp atomic
Ls[12] += L[12];
#pragma omp atomic
Ls[13] += L[13];
#pragma omp atomic
Ls[14] += L[14];
#pragma omp atomic
Ls[15] += L[15];
#pragma omp atomic
Ls[16] += L[16];
#pragma omp atomic
Ls[17] += L[17];
#pragma omp atomic
Ls[18] += L[18];
#pragma omp atomic
Ls[19] += L[19];

}

void field_m0_L2P_3(double x, double y, double z, double * L, double * F) {
double Ftmp0 = x*y;
double Ftmp1 = x*z;
double Ftmp2 = y*z;
double Ftmp3 = (1.0/2.0)*(x*x);
double Ftmp4 = (1.0/2.0)*(y*y);
double Ftmp5 = (1.0/2.0)*(z*z);
#pragma omp atomic
F[0] += -Ftmp0*L[11] - Ftmp1*L[12] - Ftmp2*L[14] - Ftmp3*L[10] - Ftmp4*L[13] - Ftmp5*L[15] - x*L[4] - y*L[5] - z*L[6] - L[1];
#pragma omp atomic
F[1] += -Ftmp0*L[13] - Ftmp1*L[14] - Ftmp2*L[17] - Ftmp3*L[11] - Ftmp4*L[16] - Ftmp5*L[18] - x*L[5] - y*L[7] - z*L[8] - L[2];
#pragma omp atomic
F[2] += -Ftmp0*L[14] - Ftmp1*L[15] - Ftmp2*L[18] - Ftmp3*L[12] - Ftmp4*L[17] - Ftmp5*L[19] - x*L[6] - y*L[8] - z*L[9] - L[3];

}

void field_m0_M2P_3(double x, double y, double z, double * M, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = 1.0*pow(R, -3.0);
double Ftmp1 = pow(R, -5.0);
double Ftmp2 = 3.0*Ftmp1;
double Ftmp3 = Ftmp2*M[5];
double Ftmp4 = Ftmp2*z;
double Ftmp5 = Ftmp0*M[0];
double Ftmp6 = y*z;
double Ftmp7 = pow(R, -7.0);
double Ftmp8 = 15.0*Ftmp7;
double Ftmp9 = Ftmp8*M[14];
double Ftmp10 = x*y;
double Ftmp11 = Ftmp2*M[2];
double Ftmp12 = Ftmp4*M[3];
double Ftmp13 = (x*x);
double Ftmp14 = Ftmp2*M[1];
double Ftmp15 = z*M[8];
double Ftmp16 = Ftmp10*Ftmp8;
double Ftmp17 = Ftmp13*Ftmp8;
double Ftmp18 = z*M[6];
double Ftmp19 = 105.0*pow(R, -9.0);
double Ftmp20 = Ftmp13*Ftmp19;
double Ftmp21 = -9.0*Ftmp1;
double Ftmp22 = Ftmp17 + Ftmp21;
double Ftmp23 = -Ftmp2;
double Ftmp24 = (y*y);
double Ftmp25 = Ftmp24*Ftmp8;
double Ftmp26 = Ftmp23 + Ftmp25;
double Ftmp27 = (z*z);
double Ftmp28 = Ftmp27*Ftmp8;
double Ftmp29 = Ftmp23 + Ftmp28;
double Ftmp30 = Ftmp26*M[7];
double Ftmp31 = Ftmp29*M[9];
double Ftmp32 = -45.0*Ftmp7;
double Ftmp33 = Ftmp20 + Ftmp32;
double Ftmp34 = Ftmp10*Ftmp33;
double Ftmp35 = Ftmp19*Ftmp24;
double Ftmp36 = Ftmp32 + Ftmp35;
double Ftmp37 = Ftmp36*M[16];
double Ftmp38 = -Ftmp8;
double Ftmp39 = Ftmp19*Ftmp27;
double Ftmp40 = Ftmp38 + Ftmp39;
double Ftmp41 = 1.0*M[18];
double Ftmp42 = Ftmp40*Ftmp41;
double Ftmp43 = x*z;
double Ftmp44 = Ftmp33*Ftmp43;
double Ftmp45 = Ftmp35 + Ftmp38;
double Ftmp46 = Ftmp45*M[17];
double Ftmp47 = Ftmp32 + Ftmp39;
double Ftmp48 = Ftmp47*M[19];
double Ftmp49 = -75.0*Ftmp7;
double Ftmp50 = 1.0*Ftmp13;
double Ftmp51 = Ftmp45*M[13];
double Ftmp52 = Ftmp40*M[15];
double Ftmp53 = Ftmp17 + Ftmp23;
double Ftmp54 = Ftmp21 + Ftmp25;
double Ftmp55 = Ftmp53*M[4];
double Ftmp56 = 1.0*Ftmp10;
double Ftmp57 = Ftmp20 + Ftmp38;
double Ftmp58 = Ftmp57*M[12];
double Ftmp59 = Ftmp57*M[11];
double Ftmp60 = x*M[6];
double Ftmp61 = y*M[8];
double Ftmp62 = Ftmp21 + Ftmp28;
double Ftmp63 = 1.0*Ftmp43;
#pragma omp atomic
F[0] += Ftmp0*M[1] - Ftmp10*Ftmp11 - Ftmp10*Ftmp37 - Ftmp10*Ftmp42 - Ftmp12*x - Ftmp13*Ftmp14 - Ftmp13*(Ftmp20 + Ftmp49)*M[10] + Ftmp15*Ftmp16 + Ftmp17*Ftmp18 + Ftmp17*y*M[5] - Ftmp20*Ftmp6*M[14] + Ftmp22*x*M[4] + Ftmp22*M[10] + Ftmp26*M[13] + Ftmp29*M[15] - Ftmp3*y + Ftmp30*x + Ftmp31*x - Ftmp34*M[11] - Ftmp4*M[6] - Ftmp43*Ftmp46 - Ftmp43*Ftmp48 - Ftmp44*M[12] + Ftmp5*x - Ftmp50*Ftmp51 - Ftmp50*Ftmp52 + Ftmp6*Ftmp9;
#pragma omp atomic
F[1] += Ftmp0*M[2] - Ftmp10*Ftmp14 - Ftmp11*Ftmp24 - Ftmp12*y + Ftmp15*Ftmp25 + Ftmp16*Ftmp18 - Ftmp24*Ftmp42 - Ftmp24*Ftmp59 - Ftmp24*(Ftmp35 + Ftmp49)*M[16] + Ftmp25*x*M[5] + Ftmp29*M[18] - Ftmp3*x + Ftmp31*y - Ftmp34*M[10] - Ftmp35*Ftmp43*M[14] - Ftmp36*Ftmp56*M[13] - Ftmp36*Ftmp6*M[17] - Ftmp4*M[8] + Ftmp43*Ftmp9 - Ftmp48*Ftmp6 + Ftmp5*y - Ftmp52*Ftmp56 + Ftmp53*M[11] + Ftmp54*y*M[7] + Ftmp54*M[16] + Ftmp55*y - Ftmp58*Ftmp6;
#pragma omp atomic
F[2] += Ftmp0*M[3] - Ftmp10*Ftmp39*M[14] + Ftmp10*Ftmp9 + Ftmp16*z*M[5] - Ftmp2*Ftmp27*M[3] - Ftmp2*Ftmp60 - Ftmp2*Ftmp61 + Ftmp26*M[17] - Ftmp27*Ftmp46 - Ftmp27*Ftmp58 - Ftmp27*(Ftmp39 + Ftmp49)*M[19] + Ftmp28*Ftmp60 + Ftmp28*Ftmp61 + Ftmp30*z - Ftmp37*Ftmp6 - Ftmp4*x*M[1] - Ftmp4*y*M[2] - Ftmp41*Ftmp47*Ftmp6 - Ftmp44*M[10] - Ftmp47*Ftmp63*M[15] + Ftmp5*z - Ftmp51*Ftmp63 + Ftmp53*M[12] + Ftmp55*z - Ftmp59*Ftmp6 + Ftmp62*z*M[9] + Ftmp62*M[19];

}

void field_m0_P2M_4(double x, double y, double z, double q, double * M) {
double Mtmp0 = q*x;
double Mtmp1 = q*y;
double Mtmp2 = q*z;
double Mtmp3 = (x*x);
double Mtmp4 = (1.0/2.0)*q;
double Mtmp5 = Mtmp0*y;
double Mtmp6 = Mtmp0*z;
double Mtmp7 = (y*y);
double Mtmp8 = Mtmp1*z;
double Mtmp9 = (z*z);
double Mtmp10 = (x*x*x);
double Mtmp11 = (1.0/6.0)*q;
double Mtmp12 = (1.0/2.0)*Mtmp3;
double Mtmp13 = (1.0/2.0)*Mtmp0;
double Mtmp14 = (y*y*y);
double Mtmp15 = (1.0/2.0)*Mtmp7;
double Mtmp16 = (1.0/2.0)*Mtmp9;
double Mtmp17 = (z*z*z);
double Mtmp18 = (1.0/24.0)*q;
double Mtmp19 = (1.0/6.0)*Mtmp10;
double Mtmp20 = (1.0/4.0)*Mtmp3*q;
double Mtmp21 = (1.0/6.0)*Mtmp0;
M[0] += q;
M[1] += -Mtmp0;
M[2] += -Mtmp1;
M[3] += -Mtmp2;
M[4] += Mtmp3*Mtmp4;
M[5] += Mtmp5;
M[6] += Mtmp6;
M[7] += Mtmp4*Mtmp7;
M[8] += Mtmp8;
M[9] += Mtmp4*Mtmp9;
M[10] += -Mtmp10*Mtmp11;
M[11] += -Mtmp1*Mtmp12;
M[12] += -Mtmp12*Mtmp2;
M[13] += -Mtmp13*Mtmp7;
M[14] += -Mtmp5*z;
M[15] += -Mtmp13*Mtmp9;
M[16] += -Mtmp11*Mtmp14;
M[17] += -Mtmp15*Mtmp2;
M[18] += -Mtmp1*Mtmp16;
M[19] += -Mtmp11*Mtmp17;
M[20] += Mtmp18*(x*x*x*x);
M[21] += Mtmp1*Mtmp19;
M[22] += Mtmp19*Mtmp2;
M[23] += Mtmp20*Mtmp7;
M[24] += Mtmp12*Mtmp8;
M[25] += Mtmp20*Mtmp9;
M[26] += Mtmp14*Mtmp21;
M[27] += Mtmp15*Mtmp6;
M[28] += Mtmp16*Mtmp5;
M[29] += Mtmp17*Mtmp21;
M[30] += Mtmp18*(y*y*y*y);
M[31] += (1.0/6.0)*Mtmp14*Mtmp2;
M[32] += (1.0/4.0)*Mtmp7*Mtmp9*q;
M[33] += (1.0/6.0)*Mtmp1*Mtmp17;
M[34] += Mtmp18*(z*z*z*z);

}
void field_m0_M2M_4(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = y*M[0];
double Mstmp2 = z*M[0];
double Mstmp3 = x*M[1];
double Mstmp4 = (x*x);
double Mstmp5 = (1.0/2.0)*M[0];
double Mstmp6 = x*M[2];
double Mstmp7 = y*M[1];
double Mstmp8 = Mstmp0*y;
double Mstmp9 = x*M[3];
double Mstmp10 = z*M[1];
double Mstmp11 = Mstmp0*z;
double Mstmp12 = y*M[2];
double Mstmp13 = (y*y);
double Mstmp14 = y*M[3];
double Mstmp15 = z*M[2];
double Mstmp16 = Mstmp1*z;
double Mstmp17 = z*M[3];
double Mstmp18 = (z*z);
double Mstmp19 = x*M[4];
double Mstmp20 = (1.0/2.0)*Mstmp4;
double Mstmp21 = (x*x*x);
double Mstmp22 = (1.0/6.0)*M[0];
double Mstmp23 = x*M[5];
double Mstmp24 = y*M[4];
double Mstmp25 = Mstmp3*y;
double Mstmp26 = x*M[6];
double Mstmp27 = x*M[7];
double Mstmp28 = y*M[5];
double Mstmp29 = Mstmp6*y;
double Mstmp30 = (1.0/2.0)*M[1];
double Mstmp31 = (1.0/2.0)*Mstmp13;
double Mstmp32 = x*M[8];
double Mstmp33 = y*M[6];
double Mstmp34 = Mstmp9*y;
double Mstmp35 = x*M[9];
double Mstmp36 = (1.0/2.0)*Mstmp18;
double Mstmp37 = y*M[7];
double Mstmp38 = (y*y*y);
double Mstmp39 = y*M[8];
double Mstmp40 = y*M[9];
double Mstmp41 = (z*z*z);
double Mstmp42 = (1.0/6.0)*Mstmp21;
double Mstmp43 = (1.0/24.0)*M[0];
double Mstmp44 = (1.0/4.0)*Mstmp4*M[0];
double Mstmp45 = (1.0/6.0)*M[1];
double Mstmp46 = (1.0/6.0)*Mstmp38;
double Mstmp47 = (1.0/6.0)*Mstmp41;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += Mstmp0 + M[1];
#pragma omp atomic
Ms[2] += Mstmp1 + M[2];
#pragma omp atomic
Ms[3] += Mstmp2 + M[3];
#pragma omp atomic
Ms[4] += Mstmp3 + Mstmp4*Mstmp5 + M[4];
#pragma omp atomic
Ms[5] += Mstmp6 + Mstmp7 + Mstmp8 + M[5];
#pragma omp atomic
Ms[6] += Mstmp10 + Mstmp11 + Mstmp9 + M[6];
#pragma omp atomic
Ms[7] += Mstmp12 + Mstmp13*Mstmp5 + M[7];
#pragma omp atomic
Ms[8] += Mstmp14 + Mstmp15 + Mstmp16 + M[8];
#pragma omp atomic
Ms[9] += Mstmp17 + Mstmp18*Mstmp5 + M[9];
#pragma omp atomic
Ms[10] += Mstmp19 + Mstmp20*M[1] + Mstmp21*Mstmp22 + M[10];
#pragma omp atomic
Ms[11] += Mstmp1*Mstmp20 + Mstmp20*M[2] + Mstmp23 + Mstmp24 + Mstmp25 + M[11];
#pragma omp atomic
Ms[12] += Mstmp2*Mstmp20 + Mstmp20*M[3] + Mstmp26 + Mstmp3*z + z*M[4] + M[12];
#pragma omp atomic
Ms[13] += Mstmp0*Mstmp31 + Mstmp13*Mstmp30 + Mstmp27 + Mstmp28 + Mstmp29 + M[13];
#pragma omp atomic
Ms[14] += Mstmp32 + Mstmp33 + Mstmp34 + Mstmp6*z + Mstmp7*z + Mstmp8*z + z*M[5] + M[14];
#pragma omp atomic
Ms[15] += Mstmp0*Mstmp36 + Mstmp18*Mstmp30 + Mstmp35 + Mstmp9*z + z*M[6] + M[15];
#pragma omp atomic
Ms[16] += Mstmp22*Mstmp38 + Mstmp31*M[2] + Mstmp37 + M[16];
#pragma omp atomic
Ms[17] += Mstmp12*z + Mstmp2*Mstmp31 + Mstmp31*M[3] + Mstmp39 + z*M[7] + M[17];
#pragma omp atomic
Ms[18] += Mstmp1*Mstmp36 + Mstmp14*z + Mstmp36*M[2] + Mstmp40 + z*M[8] + M[18];
#pragma omp atomic
Ms[19] += Mstmp22*Mstmp41 + Mstmp36*M[3] + z*M[9] + M[19];
#pragma omp atomic
Ms[20] += Mstmp20*M[4] + Mstmp42*M[1] + Mstmp43*(x*x*x*x) + x*M[10] + M[20];
#pragma omp atomic
Ms[21] += Mstmp1*Mstmp42 + Mstmp19*y + Mstmp20*Mstmp7 + Mstmp20*M[5] + Mstmp42*M[2] + x*M[11] + y*M[10] + M[21];
#pragma omp atomic
Ms[22] += Mstmp10*Mstmp20 + Mstmp19*z + Mstmp2*Mstmp42 + Mstmp20*M[6] + Mstmp42*M[3] + x*M[12] + z*M[10] + M[22];
#pragma omp atomic
Ms[23] += Mstmp12*Mstmp20 + Mstmp13*Mstmp44 + Mstmp20*M[7] + Mstmp23*y + Mstmp3*Mstmp31 + Mstmp31*M[4] + x*M[13] + y*M[11] + M[23];
#pragma omp atomic
Ms[24] += Mstmp14*Mstmp20 + Mstmp15*Mstmp20 + Mstmp16*Mstmp20 + Mstmp20*M[8] + Mstmp23*z + Mstmp24*z + Mstmp25*z + Mstmp26*y + x*M[14] + y*M[12] + z*M[11] + M[24];
#pragma omp atomic
Ms[25] += Mstmp17*Mstmp20 + Mstmp18*Mstmp44 + Mstmp20*M[9] + Mstmp26*z + Mstmp3*Mstmp36 + Mstmp36*M[4] + x*M[15] + z*M[12] + M[25];
#pragma omp atomic
Ms[26] += Mstmp0*Mstmp46 + Mstmp27*y + Mstmp31*Mstmp6 + Mstmp31*M[5] + Mstmp38*Mstmp45 + x*M[16] + y*M[13] + M[26];
#pragma omp atomic
Ms[27] += Mstmp10*Mstmp31 + Mstmp11*Mstmp31 + Mstmp27*z + Mstmp28*z + Mstmp29*z + Mstmp31*Mstmp9 + Mstmp31*M[6] + Mstmp32*y + x*M[17] + y*M[14] + z*M[13] + M[27];
#pragma omp atomic
Ms[28] += Mstmp32*z + Mstmp33*z + Mstmp34*z + Mstmp35*y + Mstmp36*Mstmp6 + Mstmp36*Mstmp7 + Mstmp36*Mstmp8 + Mstmp36*M[5] + x*M[18] + y*M[15] + z*M[14] + M[28];
#pragma omp atomic
Ms[29] += Mstmp0*Mstmp47 + Mstmp35*z + Mstmp36*Mstmp9 + Mstmp36*M[6] + Mstmp41*Mstmp45 + x*M[19] + z*M[15] + M[29];
#pragma omp atomic
Ms[30] += Mstmp31*M[7] + Mstmp43*(y*y*y*y) + Mstmp46*M[2] + y*M[16] + M[30];
#pragma omp atomic
Ms[31] += Mstmp15*Mstmp31 + Mstmp2*Mstmp46 + Mstmp31*M[8] + Mstmp37*z + Mstmp46*M[3] + y*M[17] + z*M[16] + M[31];
#pragma omp atomic
Ms[32] += Mstmp12*Mstmp36 + (1.0/4.0)*Mstmp13*Mstmp18*M[0] + Mstmp17*Mstmp31 + Mstmp31*M[9] + Mstmp36*M[7] + Mstmp39*z + y*M[18] + z*M[17] + M[32];
#pragma omp atomic
Ms[33] += Mstmp1*Mstmp47 + Mstmp14*Mstmp36 + Mstmp36*M[8] + Mstmp40*z + Mstmp47*M[2] + y*M[19] + z*M[18] + M[33];
#pragma omp atomic
Ms[34] += Mstmp36*M[9] + Mstmp43*(z*z*z*z) + Mstmp47*M[3] + z*M[19] + M[34];

}

void field_m0_M2L_4(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double D[35];
double Dtmp0 = 1.0*pow(R, -3.0);
double Dtmp1 = -Dtmp0;
double Dtmp2 = (x*x);
double Dtmp3 = pow(R, -5.0);
double Dtmp4 = 3.0*Dtmp3;
double Dtmp5 = Dtmp4*x;
double Dtmp6 = (y*y);
double Dtmp7 = y*z;
double Dtmp8 = 9.0*Dtmp3;
double Dtmp9 = -Dtmp8;
double Dtmp10 = pow(R, -7.0);
double Dtmp11 = 15.0*Dtmp10;
double Dtmp12 = Dtmp11*Dtmp2;
double Dtmp13 = -Dtmp4;
double Dtmp14 = Dtmp12 + Dtmp13;
double Dtmp15 = Dtmp11*Dtmp6;
double Dtmp16 = Dtmp13 + Dtmp15;
double Dtmp17 = 1.0*x;
double Dtmp18 = 105.0*pow(R, -9.0);
double Dtmp19 = 90.0*Dtmp10;
double Dtmp20 = -45.0*Dtmp10;
double Dtmp21 = Dtmp18*Dtmp2;
double Dtmp22 = x*(Dtmp20 + Dtmp21);
double Dtmp23 = -Dtmp11;
double Dtmp24 = Dtmp18*Dtmp6;
double Dtmp25 = Dtmp20 + Dtmp24;
D[0] = 1.0/R;
D[1] = -Dtmp0*x;
D[2] = -Dtmp0*y;
D[3] = -Dtmp0*z;
D[4] = Dtmp1 + Dtmp2*Dtmp4;
D[5] = Dtmp5*y;
D[6] = Dtmp5*z;
D[7] = Dtmp1 + Dtmp4*Dtmp6;
D[8] = Dtmp4*Dtmp7;
D[9] = -D[4] - D[7];
D[10] = -x*(Dtmp12 + Dtmp9);
D[11] = -Dtmp14*y;
D[12] = -Dtmp14*z;
D[13] = -Dtmp16*Dtmp17;
D[14] = -Dtmp11*Dtmp7*x;
D[15] = -D[10] - D[13];
D[16] = -y*(Dtmp15 + Dtmp9);
D[17] = -Dtmp16*z;
D[18] = -D[11] - D[16];
D[19] = -D[12] - D[17];
D[20] = Dtmp18*(x*x*x*x) - Dtmp19*Dtmp2 + Dtmp8;
D[21] = Dtmp22*y;
D[22] = Dtmp22*z;
D[23] = -Dtmp12 - Dtmp15 + Dtmp21*Dtmp6 + Dtmp4;
D[24] = Dtmp7*(Dtmp21 + Dtmp23);
D[25] = -D[20] - D[23];
D[26] = Dtmp17*Dtmp25*y;
D[27] = Dtmp17*z*(Dtmp23 + Dtmp24);
D[28] = -D[21] - D[26];
D[29] = -D[22] - D[27];
D[30] = Dtmp18*(y*y*y*y) - Dtmp19*Dtmp6 + Dtmp8;
D[31] = Dtmp25*Dtmp7;
D[32] = -D[23] - D[30];
D[33] = -D[24] - D[31];
D[34] = -D[25] - D[32];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[8]*M[8] + D[9]*M[9] + D[10]*M[10] + D[11]*M[11] + D[12]*M[12] + D[13]*M[13] + D[14]*M[14] + D[15]*M[15] + D[16]*M[16] + D[17]*M[17] + D[18]*M[18] + D[19]*M[19] + D[20]*M[20] + D[21]*M[21] + D[22]*M[22] + D[23]*M[23] + D[24]*M[24] + D[25]*M[25] + D[26]*M[26] + D[27]*M[27] + D[28]*M[28] + D[29]*M[29] + D[30]*M[30] + D[31]*M[31] + D[32]*M[32] + D[33]*M[33] + D[34]*M[34];
#pragma omp atomic
L[1] += D[1]*M[0] + D[4]*M[1] + D[5]*M[2] + D[6]*M[3] + D[10]*M[4] + D[11]*M[5] + D[12]*M[6] + D[13]*M[7] + D[14]*M[8] + D[15]*M[9] + D[20]*M[10] + D[21]*M[11] + D[22]*M[12] + D[23]*M[13] + D[24]*M[14] + D[25]*M[15] + D[26]*M[16] + D[27]*M[17] + D[28]*M[18] + D[29]*M[19];
#pragma omp atomic
L[2] += D[2]*M[0] + D[5]*M[1] + D[7]*M[2] + D[8]*M[3] + D[11]*M[4] + D[13]*M[5] + D[14]*M[6] + D[16]*M[7] + D[17]*M[8] + D[18]*M[9] + D[21]*M[10] + D[23]*M[11] + D[24]*M[12] + D[26]*M[13] + D[27]*M[14] + D[28]*M[15] + D[30]*M[16] + D[31]*M[17] + D[32]*M[18] + D[33]*M[19];
#pragma omp atomic
L[3] += D[3]*M[0] + D[6]*M[1] + D[8]*M[2] + D[9]*M[3] + D[12]*M[4] + D[14]*M[5] + D[15]*M[6] + D[17]*M[7] + D[18]*M[8] + D[19]*M[9] + D[22]*M[10] + D[24]*M[11] + D[25]*M[12] + D[27]*M[13] + D[28]*M[14] + D[29]*M[15] + D[31]*M[16] + D[32]*M[17] + D[33]*M[18] + D[34]*M[19];
#pragma omp atomic
L[4] += D[4]*M[0] + D[10]*M[1] + D[11]*M[2] + D[12]*M[3] + D[20]*M[4] + D[21]*M[5] + D[22]*M[6] + D[23]*M[7] + D[24]*M[8] + D[25]*M[9];
#pragma omp atomic
L[5] += D[5]*M[0] + D[11]*M[1] + D[13]*M[2] + D[14]*M[3] + D[21]*M[4] + D[23]*M[5] + D[24]*M[6] + D[26]*M[7] + D[27]*M[8] + D[28]*M[9];
#pragma omp atomic
L[6] += D[6]*M[0] + D[12]*M[1] + D[14]*M[2] + D[15]*M[3] + D[22]*M[4] + D[24]*M[5] + D[25]*M[6] + D[27]*M[7] + D[28]*M[8] + D[29]*M[9];
#pragma omp atomic
L[7] += D[7]*M[0] + D[13]*M[1] + D[16]*M[2] + D[17]*M[3] + D[23]*M[4] + D[26]*M[5] + D[27]*M[6] + D[30]*M[7] + D[31]*M[8] + D[32]*M[9];
#pragma omp atomic
L[8] += D[8]*M[0] + D[14]*M[1] + D[17]*M[2] + D[18]*M[3] + D[24]*M[4] + D[27]*M[5] + D[28]*M[6] + D[31]*M[7] + D[32]*M[8] + D[33]*M[9];
#pragma omp atomic
L[9] += D[9]*M[0] + D[15]*M[1] + D[18]*M[2] + D[19]*M[3] + D[25]*M[4] + D[28]*M[5] + D[29]*M[6] + D[32]*M[7] + D[33]*M[8] + D[34]*M[9];
#pragma omp atomic
L[10] += D[10]*M[0] + D[20]*M[1] + D[21]*M[2] + D[22]*M[3];
#pragma omp atomic
L[11] += D[11]*M[0] + D[21]*M[1] + D[23]*M[2] + D[24]*M[3];
#pragma omp atomic
L[12] += D[12]*M[0] + D[22]*M[1] + D[24]*M[2] + D[25]*M[3];
#pragma omp atomic
L[13] += D[13]*M[0] + D[23]*M[1] + D[26]*M[2] + D[27]*M[3];
#pragma omp atomic
L[14] += D[14]*M[0] + D[24]*M[1] + D[27]*M[2] + D[28]*M[3];
#pragma omp atomic
L[15] += D[15]*M[0] + D[25]*M[1] + D[28]*M[2] + D[29]*M[3];
#pragma omp atomic
L[16] += D[16]*M[0] + D[26]*M[1] + D[30]*M[2] + D[31]*M[3];
#pragma omp atomic
L[17] += D[17]*M[0] + D[27]*M[1] + D[31]*M[2] + D[32]*M[3];
#pragma omp atomic
L[18] += D[18]*M[0] + D[28]*M[1] + D[32]*M[2] + D[33]*M[3];
#pragma omp atomic
L[19] += D[19]*M[0] + D[29]*M[1] + D[33]*M[2] + D[34]*M[3];
#pragma omp atomic
L[20] += D[20]*M[0];
#pragma omp atomic
L[21] += D[21]*M[0];
#pragma omp atomic
L[22] += D[22]*M[0];
#pragma omp atomic
L[23] += D[23]*M[0];
#pragma omp atomic
L[24] += D[24]*M[0];
#pragma omp atomic
L[25] += D[25]*M[0];
#pragma omp atomic
L[26] += D[26]*M[0];
#pragma omp atomic
L[27] += D[27]*M[0];
#pragma omp atomic
L[28] += D[28]*M[0];
#pragma omp atomic
L[29] += D[29]*M[0];
#pragma omp atomic
L[30] += D[30]*M[0];
#pragma omp atomic
L[31] += D[31]*M[0];
#pragma omp atomic
L[32] += D[32]*M[0];
#pragma omp atomic
L[33] += D[33]*M[0];
#pragma omp atomic
L[34] += D[34]*M[0];

}

void field_m0_L2L_4(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
double Lstmp3 = z*L[14];
double Lstmp4 = Lstmp3*y;
double Lstmp5 = (x*x);
double Lstmp6 = (1.0/2.0)*Lstmp5;
double Lstmp7 = (1.0/6.0)*(x*x*x);
double Lstmp8 = (y*y);
double Lstmp9 = (1.0/2.0)*Lstmp8;
double Lstmp10 = (1.0/6.0)*(y*y*y);
double Lstmp11 = (z*z);
double Lstmp12 = (1.0/2.0)*Lstmp11;
double Lstmp13 = (1.0/6.0)*(z*z*z);
double Lstmp14 = x*L[13];
double Lstmp15 = x*L[26];
double Lstmp16 = x*L[15];
double Lstmp17 = x*L[29];
double Lstmp18 = y*L[11];
double Lstmp19 = z*L[12];
double Lstmp20 = y*L[21];
double Lstmp21 = z*L[22];
double Lstmp22 = y*L[18];
double Lstmp23 = y*L[33];
double Lstmp24 = z*L[17];
double Lstmp25 = z*L[31];
double Lstmp26 = y*L[28];
double Lstmp27 = Lstmp26*x;
double Lstmp28 = z*L[27];
double Lstmp29 = Lstmp28*x;
double Lstmp30 = z*L[24];
double Lstmp31 = Lstmp30*y;
double Lstmp32 = (1.0/4.0)*Lstmp5;
double Lstmp33 = x*L[23];
double Lstmp34 = x*L[25];
double Lstmp35 = y*L[13];
double Lstmp36 = Lstmp28*y;
double Lstmp37 = x*L[28];
double Lstmp38 = y*L[23];
double Lstmp39 = y*L[32];
double Lstmp40 = y*L[14];
double Lstmp41 = z*L[15];
double Lstmp42 = z*L[18];
double Lstmp43 = z*L[28];
double Lstmp44 = Lstmp43*y;
double Lstmp45 = x*L[27];
double Lstmp46 = y*L[24];
double Lstmp47 = z*L[25];
double Lstmp48 = z*L[32];
double Lstmp49 = y*L[26];
double Lstmp50 = y*L[27];
double Lstmp51 = z*L[29];
double Lstmp52 = z*L[33];
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x + Lstmp10*Lstmp15 + Lstmp10*Lstmp25 + Lstmp10*L[16] + Lstmp11*Lstmp32*L[25] + (1.0/4.0)*Lstmp11*Lstmp8*L[32] + Lstmp12*Lstmp16 + Lstmp12*Lstmp22 + Lstmp12*Lstmp27 + Lstmp12*L[9] + Lstmp13*Lstmp17 + Lstmp13*Lstmp23 + Lstmp13*L[19] + Lstmp14*Lstmp9 + Lstmp18*Lstmp6 + Lstmp19*Lstmp6 + Lstmp2*y + Lstmp20*Lstmp7 + Lstmp21*Lstmp7 + Lstmp24*Lstmp9 + Lstmp29*Lstmp9 + Lstmp31*Lstmp6 + Lstmp32*Lstmp8*L[23] + Lstmp4*x + Lstmp6*L[4] + Lstmp7*L[10] + Lstmp9*L[7] + (1.0/24.0)*(x*x*x*x)*L[20] + x*L[1] + (1.0/24.0)*(y*y*y*y)*L[30] + y*L[2] + (1.0/24.0)*(z*z*z*z)*L[34] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 + Lstmp10*L[26] + Lstmp12*Lstmp26 + Lstmp12*Lstmp34 + Lstmp12*L[15] + Lstmp13*L[29] + Lstmp18*x + Lstmp19*x + Lstmp20*Lstmp6 + Lstmp21*Lstmp6 + Lstmp28*Lstmp9 + Lstmp31*x + Lstmp33*Lstmp9 + Lstmp4 + Lstmp6*L[10] + Lstmp7*L[20] + Lstmp9*L[13] + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp10*L[30] + Lstmp12*Lstmp37 + Lstmp12*Lstmp39 + Lstmp12*L[18] + Lstmp13*L[33] + Lstmp15*Lstmp9 + Lstmp2 + Lstmp24*y + Lstmp25*Lstmp9 + Lstmp3*x + Lstmp30*Lstmp6 + Lstmp35*x + Lstmp36*x + Lstmp38*Lstmp6 + Lstmp6*L[11] + Lstmp7*L[21] + Lstmp9*L[16] + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += Lstmp10*L[31] + Lstmp12*Lstmp17 + Lstmp12*Lstmp23 + Lstmp12*L[19] + Lstmp13*L[34] + Lstmp40*x + Lstmp41*x + Lstmp42*y + Lstmp44*x + Lstmp45*Lstmp9 + Lstmp46*Lstmp6 + Lstmp47*Lstmp6 + Lstmp48*Lstmp9 + Lstmp6*L[12] + Lstmp7*L[22] + Lstmp9*L[17] + x*L[6] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += Lstmp12*L[25] + Lstmp18 + Lstmp19 + Lstmp20*x + Lstmp21*x + Lstmp31 + Lstmp6*L[20] + Lstmp9*L[23] + x*L[10] + L[4];
#pragma omp atomic
Ls[5] += Lstmp12*L[28] + Lstmp3 + Lstmp30*x + Lstmp35 + Lstmp36 + Lstmp38*x + Lstmp6*L[21] + Lstmp9*L[26] + x*L[11] + L[5];
#pragma omp atomic
Ls[6] += Lstmp12*L[29] + Lstmp40 + Lstmp41 + Lstmp44 + Lstmp46*x + Lstmp47*x + Lstmp6*L[22] + Lstmp9*L[27] + x*L[12] + L[6];
#pragma omp atomic
Ls[7] += Lstmp12*L[32] + Lstmp14 + Lstmp24 + Lstmp25*y + Lstmp29 + Lstmp49*x + Lstmp6*L[23] + Lstmp9*L[30] + y*L[16] + L[7];
#pragma omp atomic
Ls[8] += Lstmp12*L[33] + Lstmp42 + Lstmp43*x + Lstmp48*y + Lstmp50*x + Lstmp6*L[24] + Lstmp9*L[31] + x*L[14] + y*L[17] + L[8];
#pragma omp atomic
Ls[9] += Lstmp12*L[34] + Lstmp16 + Lstmp22 + Lstmp27 + Lstmp51*x + Lstmp52*y + Lstmp6*L[25] + Lstmp9*L[32] + z*L[19] + L[9];
#pragma omp atomic
Ls[10] += Lstmp20 + Lstmp21 + x*L[20] + L[10];
#pragma omp atomic
Ls[11] += Lstmp30 + Lstmp38 + x*L[21] + L[11];
#pragma omp atomic
Ls[12] += Lstmp46 + Lstmp47 + x*L[22] + L[12];
#pragma omp atomic
Ls[13] += Lstmp28 + Lstmp33 + Lstmp49 + L[13];
#pragma omp atomic
Ls[14] += Lstmp43 + Lstmp50 + x*L[24] + L[14];
#pragma omp atomic
Ls[15] += Lstmp26 + Lstmp34 + Lstmp51 + L[15];
#pragma omp atomic
Ls[16] += Lstmp15 + Lstmp25 + y*L[30] + L[16];
#pragma omp atomic
Ls[17] += Lstmp45 + Lstmp48 + y*L[31] + L[17];
#pragma omp atomic
Ls[18] += Lstmp37 + Lstmp39 + Lstmp52 + L[18];
#pragma omp atomic
Ls[19] += Lstmp17 + Lstmp23 + z*L[34] + L[19];
#pragma omp atomic
Ls[20] += L[20];
#pragma omp atomic
Ls[21] += L[21];
#pragma omp atomic
Ls[22] += L[22];
#pragma omp atomic
Ls[23] += L[23];
#pragma omp atomic
Ls[24] += L[24];
#pragma omp atomic
Ls[25] += L[25];
#pragma omp atomic
Ls[26] += L[26];
#pragma omp atomic
Ls[27] += L[27];
#pragma omp atomic
Ls[28] += L[28];
#pragma omp atomic
Ls[29] += L[29];
#pragma omp atomic
Ls[30] += L[30];
#pragma omp atomic
Ls[31] += L[31];
#pragma omp atomic
Ls[32] += L[32];
#pragma omp atomic
Ls[33] += L[33];
#pragma omp atomic
Ls[34] += L[34];

}

void field_m0_L2P_4(double x, double y, double z, double * L, double * F) {
double Ftmp0 = x*y;
double Ftmp1 = x*z;
double Ftmp2 = y*z;
double Ftmp3 = Ftmp0*z;
double Ftmp4 = (1.0/2.0)*(x*x);
double Ftmp5 = (1.0/6.0)*(x*x*x);
double Ftmp6 = (1.0/2.0)*(y*y);
double Ftmp7 = (1.0/6.0)*(y*y*y);
double Ftmp8 = (1.0/2.0)*(z*z);
double Ftmp9 = (1.0/6.0)*(z*z*z);
double Ftmp10 = Ftmp6*x;
double Ftmp11 = Ftmp8*x;
double Ftmp12 = Ftmp4*y;
double Ftmp13 = Ftmp4*z;
double Ftmp14 = Ftmp8*y;
double Ftmp15 = Ftmp6*z;
#pragma omp atomic
F[0] += -Ftmp0*L[11] - Ftmp1*L[12] - Ftmp10*L[23] - Ftmp11*L[25] - Ftmp12*L[21] - Ftmp13*L[22] - Ftmp14*L[28] - Ftmp15*L[27] - Ftmp2*L[14] - Ftmp3*L[24] - Ftmp4*L[10] - Ftmp5*L[20] - Ftmp6*L[13] - Ftmp7*L[26] - Ftmp8*L[15] - Ftmp9*L[29] - x*L[4] - y*L[5] - z*L[6] - L[1];
#pragma omp atomic
F[1] += -Ftmp0*L[13] - Ftmp1*L[14] - Ftmp10*L[26] - Ftmp11*L[28] - Ftmp12*L[23] - Ftmp13*L[24] - Ftmp14*L[32] - Ftmp15*L[31] - Ftmp2*L[17] - Ftmp3*L[27] - Ftmp4*L[11] - Ftmp5*L[21] - Ftmp6*L[16] - Ftmp7*L[30] - Ftmp8*L[18] - Ftmp9*L[33] - x*L[5] - y*L[7] - z*L[8] - L[2];
#pragma omp atomic
F[2] += -Ftmp0*L[14] - Ftmp1*L[15] - Ftmp10*L[27] - Ftmp11*L[29] - Ftmp12*L[24] - Ftmp13*L[25] - Ftmp14*L[33] - Ftmp15*L[32] - Ftmp2*L[18] - Ftmp3*L[28] - Ftmp4*L[12] - Ftmp5*L[22] - Ftmp6*L[17] - Ftmp7*L[31] - Ftmp8*L[19] - Ftmp9*L[34] - x*L[6] - y*L[8] - z*L[9] - L[3];

}

void field_m0_M2P_4(double x, double y, double z, double * M, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = 1.0*pow(R, -3.0);
double Ftmp1 = pow(R, -5.0);
double Ftmp2 = 3.0*Ftmp1;
double Ftmp3 = Ftmp2*M[5];
double Ftmp4 = Ftmp2*z;
double Ftmp5 = Ftmp0*M[0];
double Ftmp6 = y*z;
double Ftmp7 = pow(R, -7.0);
double Ftmp8 = 15.0*Ftmp7;
double Ftmp9 = Ftmp8*M[14];
double Ftmp10 = Ftmp2*M[2];
double Ftmp11 = x*y;
double Ftmp12 = Ftmp4*M[3];
double Ftmp13 = (x*x);
double Ftmp14 = Ftmp2*M[1];
double Ftmp15 = y*M[8];
double Ftmp16 = x*z;
double Ftmp17 = Ftmp16*Ftmp8;
double Ftmp18 = Ftmp13*Ftmp8;
double Ftmp19 = y*M[5];
double Ftmp20 = pow(R, -9.0);
double Ftmp21 = 105.0*Ftmp20;
double Ftmp22 = Ftmp13*Ftmp21;
double Ftmp23 = -9.0*Ftmp1;
double Ftmp24 = Ftmp18 + Ftmp23;
double Ftmp25 = -Ftmp2;
double Ftmp26 = (y*y);
double Ftmp27 = Ftmp26*Ftmp8;
double Ftmp28 = Ftmp25 + Ftmp27;
double Ftmp29 = (z*z);
double Ftmp30 = Ftmp29*Ftmp8;
double Ftmp31 = Ftmp25 + Ftmp30;
double Ftmp32 = Ftmp28*M[7];
double Ftmp33 = Ftmp31*M[9];
double Ftmp34 = 45.0*Ftmp7;
double Ftmp35 = -Ftmp34;
double Ftmp36 = Ftmp22 + Ftmp35;
double Ftmp37 = Ftmp36*M[21];
double Ftmp38 = Ftmp21*Ftmp26;
double Ftmp39 = Ftmp35 + Ftmp38;
double Ftmp40 = Ftmp39*y;
double Ftmp41 = 1.0*M[26];
double Ftmp42 = 35.0*Ftmp20;
double Ftmp43 = 3.0*M[28];
double Ftmp44 = Ftmp43*(Ftmp29*Ftmp42 - 5.0*Ftmp7);
double Ftmp45 = Ftmp36*M[22];
double Ftmp46 = 1.0*z;
double Ftmp47 = -Ftmp8;
double Ftmp48 = Ftmp38 + Ftmp47;
double Ftmp49 = Ftmp48*M[27];
double Ftmp50 = Ftmp21*Ftmp29;
double Ftmp51 = Ftmp35 + Ftmp50;
double Ftmp52 = Ftmp46*Ftmp51;
double Ftmp53 = Ftmp36*x;
double Ftmp54 = Ftmp53*y;
double Ftmp55 = Ftmp40*M[16];
double Ftmp56 = Ftmp47 + Ftmp50;
double Ftmp57 = 1.0*Ftmp56*M[18];
double Ftmp58 = Ftmp53*z;
double Ftmp59 = Ftmp48*M[17];
double Ftmp60 = Ftmp51*M[19];
double Ftmp61 = 315.0*Ftmp20;
double Ftmp62 = -Ftmp61;
double Ftmp63 = pow(R, -11.0);
double Ftmp64 = 945.0*Ftmp63;
double Ftmp65 = Ftmp13*Ftmp64;
double Ftmp66 = Ftmp62 + Ftmp65;
double Ftmp67 = Ftmp16*y;
double Ftmp68 = Ftmp26*Ftmp64;
double Ftmp69 = Ftmp62 + Ftmp68;
double Ftmp70 = Ftmp69*M[31];
double Ftmp71 = -75.0*Ftmp7;
double Ftmp72 = 1.0*Ftmp13;
double Ftmp73 = Ftmp48*M[13];
double Ftmp74 = Ftmp56*M[15];
double Ftmp75 = -525.0*Ftmp20;
double Ftmp76 = Ftmp13*(Ftmp65 + Ftmp75);
double Ftmp77 = y*M[21];
double Ftmp78 = z*M[22];
double Ftmp79 = y*M[33];
double Ftmp80 = Ftmp29*Ftmp64;
double Ftmp81 = Ftmp62 + Ftmp80;
double Ftmp82 = Ftmp46*Ftmp81;
double Ftmp83 = Ftmp13*y;
double Ftmp84 = Ftmp41*Ftmp69;
double Ftmp85 = 315.0*Ftmp29*Ftmp63;
double Ftmp86 = Ftmp43*(-Ftmp42 + Ftmp85);
double Ftmp87 = Ftmp13*Ftmp46;
double Ftmp88 = -Ftmp21;
double Ftmp89 = (Ftmp68 + Ftmp88)*M[27];
double Ftmp90 = Ftmp81*M[29];
double Ftmp91 = 225.0*Ftmp7;
double Ftmp92 = Ftmp64*(x*x*x*x);
double Ftmp93 = 1050.0*Ftmp20;
double Ftmp94 = Ftmp64*(y*y*y*y);
double Ftmp95 = 630.0*Ftmp20;
double Ftmp96 = (-Ftmp26*Ftmp95 + Ftmp34 + Ftmp94)*M[30];
double Ftmp97 = Ftmp64*(z*z*z*z);
double Ftmp98 = (-Ftmp29*Ftmp95 + Ftmp34 + Ftmp97)*M[34];
double Ftmp99 = -Ftmp26*Ftmp61;
double Ftmp100 = Ftmp26*Ftmp65;
double Ftmp101 = -Ftmp22;
double Ftmp102 = Ftmp101 + Ftmp34;
double Ftmp103 = -Ftmp29*Ftmp61;
double Ftmp104 = Ftmp29*Ftmp65;
double Ftmp105 = -Ftmp50;
double Ftmp106 = Ftmp105 + Ftmp8;
double Ftmp107 = -Ftmp38;
double Ftmp108 = Ftmp29*Ftmp68;
double Ftmp109 = Ftmp107 + Ftmp108;
double Ftmp110 = Ftmp18 + Ftmp25;
double Ftmp111 = Ftmp23 + Ftmp27;
double Ftmp112 = Ftmp110*M[4];
double Ftmp113 = Ftmp41*x;
double Ftmp114 = Ftmp22 + Ftmp47;
double Ftmp115 = Ftmp114*M[24];
double Ftmp116 = Ftmp39*M[31];
double Ftmp117 = 1.0*x;
double Ftmp118 = Ftmp114*M[12];
double Ftmp119 = Ftmp66*x;
double Ftmp120 = Ftmp114*M[11];
double Ftmp121 = Ftmp26*z;
double Ftmp122 = (Ftmp65 + Ftmp88)*M[24];
double Ftmp123 = Ftmp68 + Ftmp75;
double Ftmp124 = Ftmp11*Ftmp46;
double Ftmp125 = (-Ftmp13*Ftmp95 + Ftmp34 + Ftmp92)*M[20];
double Ftmp126 = Ftmp100 + Ftmp107;
double Ftmp127 = -Ftmp13*Ftmp61 + Ftmp34;
double Ftmp128 = x*M[6];
double Ftmp129 = Ftmp23 + Ftmp30;
double Ftmp130 = Ftmp117*M[29];
double Ftmp131 = 1.0*Ftmp79;
double Ftmp132 = Ftmp29*y;
double Ftmp133 = Ftmp29*(Ftmp75 + Ftmp80);
#pragma omp atomic
F[0] += Ftmp0*M[1] - Ftmp10*Ftmp11 - Ftmp11*Ftmp57 - Ftmp12*x - Ftmp13*Ftmp14 - Ftmp13*(Ftmp22 + Ftmp71)*M[10] + Ftmp15*Ftmp17 - Ftmp16*Ftmp59 - Ftmp16*Ftmp60 + Ftmp18*Ftmp19 + Ftmp18*z*M[6] - Ftmp22*Ftmp6*M[14] + Ftmp24*x*M[4] + Ftmp24*M[10] + Ftmp28*M[13] - Ftmp3*y + Ftmp31*M[15] + Ftmp32*x + Ftmp33*x - Ftmp37*y - Ftmp4*M[6] - Ftmp40*Ftmp41 - Ftmp44*y - Ftmp45*z - Ftmp46*Ftmp49 + Ftmp5*x - Ftmp52*M[29] - Ftmp54*M[11] - Ftmp55*x - Ftmp58*M[12] + Ftmp6*Ftmp9 + Ftmp66*Ftmp67*M[24] + Ftmp67*Ftmp70 - Ftmp72*Ftmp73 - Ftmp72*Ftmp74 + Ftmp76*Ftmp77 + Ftmp76*Ftmp78 + Ftmp79*Ftmp82*x + Ftmp83*Ftmp84 + Ftmp83*Ftmp86 + Ftmp87*Ftmp89 + Ftmp87*Ftmp90 + Ftmp96*x + Ftmp98*x + x*(Ftmp106 + Ftmp109)*M[32] + x*(Ftmp100 + Ftmp102 + Ftmp99)*M[23] + x*(Ftmp102 + Ftmp103 + Ftmp104)*M[25] + x*(-Ftmp13*Ftmp93 + Ftmp91 + Ftmp92)*M[20];
#pragma omp atomic
F[1] += Ftmp0*M[2] - Ftmp10*Ftmp26 - Ftmp11*Ftmp14 + Ftmp110*M[11] + Ftmp111*y*M[7] + Ftmp111*M[16] + Ftmp112*y + Ftmp113*Ftmp123*Ftmp26 - Ftmp113*Ftmp39 - Ftmp115*z - Ftmp116*z - Ftmp117*Ftmp40*M[13] - Ftmp117*Ftmp74*y - Ftmp118*Ftmp6 + Ftmp119*Ftmp26*M[21] + Ftmp119*Ftmp78*y - Ftmp12*y - Ftmp120*Ftmp26 + Ftmp121*Ftmp122 + Ftmp121*Ftmp123*M[31] + Ftmp124*Ftmp69*M[27] + Ftmp124*Ftmp90 + Ftmp125*y - Ftmp16*Ftmp38*M[14] + Ftmp16*Ftmp9 - Ftmp26*Ftmp57 + Ftmp26*Ftmp82*M[33] + Ftmp26*Ftmp86*x - Ftmp26*(Ftmp38 + Ftmp71)*M[16] + Ftmp27*x*M[5] + Ftmp27*z*M[8] - Ftmp3*x + Ftmp31*M[18] + Ftmp33*y - Ftmp37*x - Ftmp4*M[8] - Ftmp40*z*M[17] - Ftmp44*x + Ftmp5*y - Ftmp52*M[33] - Ftmp54*M[10] - Ftmp6*Ftmp60 + Ftmp67*Ftmp8*M[6] + Ftmp98*y + y*(Ftmp126 + Ftmp127)*M[23] + y*(Ftmp101 + Ftmp104 + Ftmp106)*M[25] + y*(Ftmp103 + Ftmp109 + Ftmp34)*M[32] + y*(-Ftmp26*Ftmp93 + Ftmp91 + Ftmp94)*M[30];
#pragma omp atomic
F[2] += Ftmp0*M[3] - Ftmp11*Ftmp50*M[14] + Ftmp11*Ftmp9 + Ftmp110*M[12] + Ftmp112*z - Ftmp115*y - Ftmp116*y + Ftmp117*Ftmp29*Ftmp89 - Ftmp117*Ftmp49 - Ftmp118*Ftmp29 + Ftmp119*Ftmp29*M[22] - Ftmp120*Ftmp6 + Ftmp122*Ftmp132 + Ftmp125*z - Ftmp128*Ftmp2 + Ftmp128*Ftmp30 + Ftmp129*z*M[9] + Ftmp129*M[19] + Ftmp130*Ftmp133 - Ftmp130*Ftmp51 + Ftmp131*Ftmp133 - Ftmp131*Ftmp51 + Ftmp132*Ftmp70 - Ftmp15*Ftmp2 + Ftmp15*Ftmp30 + Ftmp16*Ftmp66*Ftmp77 + Ftmp17*Ftmp19 - Ftmp2*Ftmp29*M[3] + Ftmp28*M[17] - Ftmp29*Ftmp59 - Ftmp29*(Ftmp50 + Ftmp71)*M[19] + Ftmp32*z - Ftmp4*x*M[1] - Ftmp4*y*M[2] + Ftmp43*Ftmp67*(Ftmp85 + Ftmp88) - Ftmp45*x - Ftmp46*Ftmp73*x + Ftmp5*z - Ftmp52*x*M[15] - Ftmp52*y*M[18] - Ftmp55*z - Ftmp58*M[10] + Ftmp67*Ftmp84 + Ftmp96*z + z*(Ftmp101 + Ftmp126 + Ftmp8)*M[23] + z*(Ftmp104 + Ftmp105 + Ftmp127)*M[25] + z*(-Ftmp29*Ftmp93 + Ftmp91 + Ftmp97)*M[34] + z*(Ftmp105 + Ftmp108 + Ftmp34 + Ftmp99)*M[32];

}

void field_m0_P2M_5(double x, double y, double z, double q, double * M) {
double Mtmp0 = q*x;
double Mtmp1 = q*y;
double Mtmp2 = q*z;
double Mtmp3 = (x*x);
double Mtmp4 = (1.0/2.0)*q;
double Mtmp5 = Mtmp0*y;
double Mtmp6 = Mtmp0*z;
double Mtmp7 = (y*y);
double Mtmp8 = Mtmp1*z;
double Mtmp9 = (z*z);
double Mtmp10 = (x*x*x);
double Mtmp11 = (1.0/6.0)*q;
double Mtmp12 = (1.0/2.0)*Mtmp3;
double Mtmp13 = (1.0/2.0)*Mtmp0;
double Mtmp14 = (y*y*y);
double Mtmp15 = (1.0/2.0)*Mtmp7;
double Mtmp16 = (1.0/2.0)*Mtmp9;
double Mtmp17 = (z*z*z);
double Mtmp18 = (x*x*x*x);
double Mtmp19 = (1.0/24.0)*q;
double Mtmp20 = (1.0/6.0)*Mtmp10;
double Mtmp21 = Mtmp7*q;
double Mtmp22 = (1.0/4.0)*Mtmp3;
double Mtmp23 = Mtmp9*q;
double Mtmp24 = (1.0/6.0)*Mtmp0;
double Mtmp25 = (y*y*y*y);
double Mtmp26 = (1.0/6.0)*Mtmp14;
double Mtmp27 = (1.0/4.0)*Mtmp9;
double Mtmp28 = (1.0/6.0)*Mtmp17;
double Mtmp29 = (z*z*z*z);
double Mtmp30 = (1.0/120.0)*q;
double Mtmp31 = (1.0/24.0)*Mtmp18;
double Mtmp32 = (1.0/12.0)*Mtmp10;
double Mtmp33 = (1.0/12.0)*Mtmp14;
double Mtmp34 = Mtmp3*q;
double Mtmp35 = (1.0/12.0)*Mtmp17;
double Mtmp36 = (1.0/24.0)*Mtmp0;
M[0] += q;
M[1] += -Mtmp0;
M[2] += -Mtmp1;
M[3] += -Mtmp2;
M[4] += Mtmp3*Mtmp4;
M[5] += Mtmp5;
M[6] += Mtmp6;
M[7] += Mtmp4*Mtmp7;
M[8] += Mtmp8;
M[9] += Mtmp4*Mtmp9;
M[10] += -Mtmp10*Mtmp11;
M[11] += -Mtmp1*Mtmp12;
M[12] += -Mtmp12*Mtmp2;
M[13] += -Mtmp13*Mtmp7;
M[14] += -Mtmp5*z;
M[15] += -Mtmp13*Mtmp9;
M[16] += -Mtmp11*Mtmp14;
M[17] += -Mtmp15*Mtmp2;
M[18] += -Mtmp1*Mtmp16;
M[19] += -Mtmp11*Mtmp17;
M[20] += Mtmp18*Mtmp19;
M[21] += Mtmp1*Mtmp20;
M[22] += Mtmp2*Mtmp20;
M[23] += Mtmp21*Mtmp22;
M[24] += Mtmp12*Mtmp8;
M[25] += Mtmp22*Mtmp23;
M[26] += Mtmp14*Mtmp24;
M[27] += Mtmp15*Mtmp6;
M[28] += Mtmp16*Mtmp5;
M[29] += Mtmp17*Mtmp24;
M[30] += Mtmp19*Mtmp25;
M[31] += Mtmp2*Mtmp26;
M[32] += Mtmp21*Mtmp27;
M[33] += Mtmp1*Mtmp28;
M[34] += Mtmp19*Mtmp29;
M[35] += -Mtmp30*(x*x*x*x*x);
M[36] += -Mtmp1*Mtmp31;
M[37] += -Mtmp2*Mtmp31;
M[38] += -Mtmp21*Mtmp32;
M[39] += -Mtmp20*Mtmp8;
M[40] += -Mtmp23*Mtmp32;
M[41] += -Mtmp33*Mtmp34;
M[42] += -Mtmp2*Mtmp22*Mtmp7;
M[43] += -Mtmp1*Mtmp22*Mtmp9;
M[44] += -Mtmp34*Mtmp35;
M[45] += -Mtmp25*Mtmp36;
M[46] += -Mtmp26*Mtmp6;
M[47] += -Mtmp0*Mtmp27*Mtmp7;
M[48] += -Mtmp28*Mtmp5;
M[49] += -Mtmp29*Mtmp36;
M[50] += -Mtmp30*(y*y*y*y*y);
M[51] += -1.0/24.0*Mtmp2*Mtmp25;
M[52] += -Mtmp23*Mtmp33;
M[53] += -Mtmp21*Mtmp35;
M[54] += -1.0/24.0*Mtmp1*Mtmp29;
M[55] += -Mtmp30*(z*z*z*z*z);

}
void field_m0_M2M_5(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = y*M[0];
double Mstmp2 = z*M[0];
double Mstmp3 = x*M[1];
double Mstmp4 = (x*x);
double Mstmp5 = (1.0/2.0)*M[0];
double Mstmp6 = x*M[2];
double Mstmp7 = y*M[1];
double Mstmp8 = Mstmp0*y;
double Mstmp9 = x*M[3];
double Mstmp10 = z*M[1];
double Mstmp11 = Mstmp0*z;
double Mstmp12 = y*M[2];
double Mstmp13 = (y*y);
double Mstmp14 = y*M[3];
double Mstmp15 = z*M[2];
double Mstmp16 = Mstmp1*z;
double Mstmp17 = z*M[3];
double Mstmp18 = (z*z);
double Mstmp19 = x*M[4];
double Mstmp20 = (1.0/2.0)*Mstmp4;
double Mstmp21 = (x*x*x);
double Mstmp22 = (1.0/6.0)*M[0];
double Mstmp23 = x*M[5];
double Mstmp24 = y*M[4];
double Mstmp25 = Mstmp3*y;
double Mstmp26 = x*M[6];
double Mstmp27 = z*M[4];
double Mstmp28 = Mstmp3*z;
double Mstmp29 = x*M[7];
double Mstmp30 = y*M[5];
double Mstmp31 = Mstmp6*y;
double Mstmp32 = (1.0/2.0)*M[1];
double Mstmp33 = (1.0/2.0)*Mstmp13;
double Mstmp34 = x*M[8];
double Mstmp35 = y*M[6];
double Mstmp36 = z*M[5];
double Mstmp37 = Mstmp9*y;
double Mstmp38 = Mstmp6*z;
double Mstmp39 = Mstmp7*z;
double Mstmp40 = x*M[9];
double Mstmp41 = z*M[6];
double Mstmp42 = Mstmp9*z;
double Mstmp43 = (1.0/2.0)*Mstmp18;
double Mstmp44 = y*M[7];
double Mstmp45 = (y*y*y);
double Mstmp46 = y*M[8];
double Mstmp47 = z*M[7];
double Mstmp48 = Mstmp12*z;
double Mstmp49 = y*M[9];
double Mstmp50 = z*M[8];
double Mstmp51 = Mstmp14*z;
double Mstmp52 = z*M[9];
double Mstmp53 = (z*z*z);
double Mstmp54 = x*M[10];
double Mstmp55 = (1.0/6.0)*Mstmp21;
double Mstmp56 = (x*x*x*x);
double Mstmp57 = (1.0/24.0)*M[0];
double Mstmp58 = x*M[11];
double Mstmp59 = y*M[10];
double Mstmp60 = Mstmp19*y;
double Mstmp61 = x*M[12];
double Mstmp62 = x*M[13];
double Mstmp63 = y*M[11];
double Mstmp64 = Mstmp23*y;
double Mstmp65 = Mstmp13*M[0];
double Mstmp66 = (1.0/4.0)*Mstmp4;
double Mstmp67 = x*M[14];
double Mstmp68 = y*M[12];
double Mstmp69 = Mstmp26*y;
double Mstmp70 = x*M[15];
double Mstmp71 = Mstmp18*Mstmp66;
double Mstmp72 = x*M[16];
double Mstmp73 = y*M[13];
double Mstmp74 = Mstmp29*y;
double Mstmp75 = (1.0/6.0)*M[1];
double Mstmp76 = (1.0/6.0)*Mstmp45;
double Mstmp77 = x*M[17];
double Mstmp78 = y*M[14];
double Mstmp79 = Mstmp34*y;
double Mstmp80 = x*M[18];
double Mstmp81 = y*M[15];
double Mstmp82 = Mstmp40*y;
double Mstmp83 = x*M[19];
double Mstmp84 = (1.0/6.0)*Mstmp53;
double Mstmp85 = y*M[16];
double Mstmp86 = (y*y*y*y);
double Mstmp87 = y*M[17];
double Mstmp88 = y*M[18];
double Mstmp89 = (1.0/4.0)*Mstmp18;
double Mstmp90 = y*M[19];
double Mstmp91 = (z*z*z*z);
double Mstmp92 = (1.0/24.0)*Mstmp56;
double Mstmp93 = (1.0/120.0)*M[0];
double Mstmp94 = Mstmp13*Mstmp66;
double Mstmp95 = (1.0/12.0)*Mstmp21;
double Mstmp96 = Mstmp18*M[0];
double Mstmp97 = (1.0/12.0)*Mstmp45;
double Mstmp98 = Mstmp4*M[0];
double Mstmp99 = (1.0/12.0)*Mstmp53;
double Mstmp100 = (1.0/24.0)*M[1];
double Mstmp101 = (1.0/24.0)*Mstmp86;
double Mstmp102 = Mstmp13*Mstmp89;
double Mstmp103 = (1.0/24.0)*Mstmp91;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += Mstmp0 + M[1];
#pragma omp atomic
Ms[2] += Mstmp1 + M[2];
#pragma omp atomic
Ms[3] += Mstmp2 + M[3];
#pragma omp atomic
Ms[4] += Mstmp3 + Mstmp4*Mstmp5 + M[4];
#pragma omp atomic
Ms[5] += Mstmp6 + Mstmp7 + Mstmp8 + M[5];
#pragma omp atomic
Ms[6] += Mstmp10 + Mstmp11 + Mstmp9 + M[6];
#pragma omp atomic
Ms[7] += Mstmp12 + Mstmp13*Mstmp5 + M[7];
#pragma omp atomic
Ms[8] += Mstmp14 + Mstmp15 + Mstmp16 + M[8];
#pragma omp atomic
Ms[9] += Mstmp17 + Mstmp18*Mstmp5 + M[9];
#pragma omp atomic
Ms[10] += Mstmp19 + Mstmp20*M[1] + Mstmp21*Mstmp22 + M[10];
#pragma omp atomic
Ms[11] += Mstmp1*Mstmp20 + Mstmp20*M[2] + Mstmp23 + Mstmp24 + Mstmp25 + M[11];
#pragma omp atomic
Ms[12] += Mstmp2*Mstmp20 + Mstmp20*M[3] + Mstmp26 + Mstmp27 + Mstmp28 + M[12];
#pragma omp atomic
Ms[13] += Mstmp0*Mstmp33 + Mstmp13*Mstmp32 + Mstmp29 + Mstmp30 + Mstmp31 + M[13];
#pragma omp atomic
Ms[14] += Mstmp34 + Mstmp35 + Mstmp36 + Mstmp37 + Mstmp38 + Mstmp39 + Mstmp8*z + M[14];
#pragma omp atomic
Ms[15] += Mstmp0*Mstmp43 + Mstmp18*Mstmp32 + Mstmp40 + Mstmp41 + Mstmp42 + M[15];
#pragma omp atomic
Ms[16] += Mstmp22*Mstmp45 + Mstmp33*M[2] + Mstmp44 + M[16];
#pragma omp atomic
Ms[17] += Mstmp2*Mstmp33 + Mstmp33*M[3] + Mstmp46 + Mstmp47 + Mstmp48 + M[17];
#pragma omp atomic
Ms[18] += Mstmp1*Mstmp43 + Mstmp43*M[2] + Mstmp49 + Mstmp50 + Mstmp51 + M[18];
#pragma omp atomic
Ms[19] += Mstmp22*Mstmp53 + Mstmp43*M[3] + Mstmp52 + M[19];
#pragma omp atomic
Ms[20] += Mstmp20*M[4] + Mstmp54 + Mstmp55*M[1] + Mstmp56*Mstmp57 + M[20];
#pragma omp atomic
Ms[21] += Mstmp1*Mstmp55 + Mstmp20*Mstmp7 + Mstmp20*M[5] + Mstmp55*M[2] + Mstmp58 + Mstmp59 + Mstmp60 + M[21];
#pragma omp atomic
Ms[22] += Mstmp10*Mstmp20 + Mstmp19*z + Mstmp2*Mstmp55 + Mstmp20*M[6] + Mstmp55*M[3] + Mstmp61 + z*M[10] + M[22];
#pragma omp atomic
Ms[23] += Mstmp12*Mstmp20 + Mstmp20*M[7] + Mstmp3*Mstmp33 + Mstmp33*M[4] + Mstmp62 + Mstmp63 + Mstmp64 + Mstmp65*Mstmp66 + M[23];
#pragma omp atomic
Ms[24] += Mstmp14*Mstmp20 + Mstmp15*Mstmp20 + Mstmp16*Mstmp20 + Mstmp20*M[8] + Mstmp23*z + Mstmp24*z + Mstmp25*z + Mstmp67 + Mstmp68 + Mstmp69 + z*M[11] + M[24];
#pragma omp atomic
Ms[25] += Mstmp17*Mstmp20 + Mstmp20*M[9] + Mstmp26*z + Mstmp3*Mstmp43 + Mstmp43*M[4] + Mstmp70 + Mstmp71*M[0] + z*M[12] + M[25];
#pragma omp atomic
Ms[26] += Mstmp0*Mstmp76 + Mstmp33*Mstmp6 + Mstmp33*M[5] + Mstmp45*Mstmp75 + Mstmp72 + Mstmp73 + Mstmp74 + M[26];
#pragma omp atomic
Ms[27] += Mstmp10*Mstmp33 + Mstmp11*Mstmp33 + Mstmp29*z + Mstmp30*z + Mstmp31*z + Mstmp33*Mstmp9 + Mstmp33*M[6] + Mstmp77 + Mstmp78 + Mstmp79 + z*M[13] + M[27];
#pragma omp atomic
Ms[28] += Mstmp34*z + Mstmp35*z + Mstmp37*z + Mstmp43*Mstmp6 + Mstmp43*Mstmp7 + Mstmp43*Mstmp8 + Mstmp43*M[5] + Mstmp80 + Mstmp81 + Mstmp82 + z*M[14] + M[28];
#pragma omp atomic
Ms[29] += Mstmp0*Mstmp84 + Mstmp40*z + Mstmp43*Mstmp9 + Mstmp43*M[6] + Mstmp53*Mstmp75 + Mstmp83 + z*M[15] + M[29];
#pragma omp atomic
Ms[30] += Mstmp33*M[7] + Mstmp57*Mstmp86 + Mstmp76*M[2] + Mstmp85 + M[30];
#pragma omp atomic
Ms[31] += Mstmp15*Mstmp33 + Mstmp2*Mstmp76 + Mstmp33*M[8] + Mstmp44*z + Mstmp76*M[3] + Mstmp87 + z*M[16] + M[31];
#pragma omp atomic
Ms[32] += Mstmp12*Mstmp43 + Mstmp17*Mstmp33 + Mstmp33*M[9] + Mstmp43*M[7] + Mstmp46*z + Mstmp65*Mstmp89 + Mstmp88 + z*M[17] + M[32];
#pragma omp atomic
Ms[33] += Mstmp1*Mstmp84 + Mstmp14*Mstmp43 + Mstmp43*M[8] + Mstmp49*z + Mstmp84*M[2] + Mstmp90 + z*M[18] + M[33];
#pragma omp atomic
Ms[34] += Mstmp43*M[9] + Mstmp57*Mstmp91 + Mstmp84*M[3] + z*M[19] + M[34];
#pragma omp atomic
Ms[35] += Mstmp20*M[10] + Mstmp55*M[4] + Mstmp92*M[1] + Mstmp93*(x*x*x*x*x) + x*M[20] + M[35];
#pragma omp atomic
Ms[36] += Mstmp1*Mstmp92 + Mstmp20*Mstmp24 + Mstmp20*M[11] + Mstmp54*y + Mstmp55*Mstmp7 + Mstmp55*M[5] + Mstmp92*M[2] + x*M[21] + y*M[20] + M[36];
#pragma omp atomic
Ms[37] += Mstmp10*Mstmp55 + Mstmp2*Mstmp92 + Mstmp20*Mstmp27 + Mstmp20*M[12] + Mstmp54*z + Mstmp55*M[6] + Mstmp92*M[3] + x*M[22] + z*M[20] + M[37];
#pragma omp atomic
Ms[38] += Mstmp12*Mstmp55 + Mstmp19*Mstmp33 + Mstmp20*Mstmp30 + Mstmp20*M[13] + Mstmp33*M[10] + Mstmp55*M[7] + Mstmp58*y + Mstmp65*Mstmp95 + Mstmp94*M[1] + x*M[23] + y*M[21] + M[38];
#pragma omp atomic
Ms[39] += Mstmp14*Mstmp55 + Mstmp15*Mstmp55 + Mstmp16*Mstmp55 + Mstmp20*Mstmp35 + Mstmp20*Mstmp36 + Mstmp20*Mstmp39 + Mstmp20*M[14] + Mstmp55*M[8] + Mstmp58*z + Mstmp59*z + Mstmp60*z + Mstmp61*y + x*M[24] + y*M[22] + z*M[21] + M[39];
#pragma omp atomic
Ms[40] += Mstmp17*Mstmp55 + Mstmp19*Mstmp43 + Mstmp20*Mstmp41 + Mstmp20*M[15] + Mstmp43*M[10] + Mstmp55*M[9] + Mstmp61*z + Mstmp71*M[1] + Mstmp95*Mstmp96 + x*M[25] + z*M[22] + M[40];
#pragma omp atomic
Ms[41] += Mstmp20*Mstmp44 + Mstmp20*M[16] + Mstmp23*Mstmp33 + Mstmp3*Mstmp76 + Mstmp33*M[11] + Mstmp62*y + Mstmp76*M[4] + Mstmp94*M[2] + Mstmp97*Mstmp98 + x*M[26] + y*M[23] + M[41];
#pragma omp atomic
Ms[42] += Mstmp2*Mstmp94 + Mstmp20*Mstmp46 + Mstmp20*Mstmp47 + Mstmp20*Mstmp48 + Mstmp20*M[17] + Mstmp26*Mstmp33 + Mstmp27*Mstmp33 + Mstmp28*Mstmp33 + Mstmp33*M[12] + Mstmp62*z + Mstmp63*z + Mstmp64*z + Mstmp67*y + Mstmp94*M[3] + x*M[27] + y*M[24] + z*M[23] + M[42];
#pragma omp atomic
Ms[43] += Mstmp1*Mstmp71 + Mstmp20*Mstmp49 + Mstmp20*Mstmp50 + Mstmp20*Mstmp51 + Mstmp20*M[18] + Mstmp23*Mstmp43 + Mstmp24*Mstmp43 + Mstmp25*Mstmp43 + Mstmp43*M[11] + Mstmp67*z + Mstmp68*z + Mstmp69*z + Mstmp70*y + Mstmp71*M[2] + x*M[28] + y*M[25] + z*M[24] + M[43];
#pragma omp atomic
Ms[44] += Mstmp20*Mstmp52 + Mstmp20*M[19] + Mstmp26*Mstmp43 + Mstmp3*Mstmp84 + Mstmp43*M[12] + Mstmp70*z + Mstmp71*M[3] + Mstmp84*M[4] + Mstmp98*Mstmp99 + x*M[29] + z*M[25] + M[44];
#pragma omp atomic
Ms[45] += Mstmp0*Mstmp101 + Mstmp100*Mstmp86 + Mstmp29*Mstmp33 + Mstmp33*M[13] + Mstmp6*Mstmp76 + Mstmp72*y + Mstmp76*M[5] + x*M[30] + y*M[26] + M[45];
#pragma omp atomic
Ms[46] += Mstmp10*Mstmp76 + Mstmp11*Mstmp76 + Mstmp33*Mstmp34 + Mstmp33*Mstmp36 + Mstmp33*Mstmp38 + Mstmp33*M[14] + Mstmp72*z + Mstmp73*z + Mstmp74*z + Mstmp76*Mstmp9 + Mstmp76*M[6] + Mstmp77*y + x*M[31] + y*M[27] + z*M[26] + M[46];
#pragma omp atomic
Ms[47] += Mstmp0*Mstmp102 + Mstmp102*M[1] + Mstmp29*Mstmp43 + Mstmp30*Mstmp43 + Mstmp31*Mstmp43 + Mstmp33*Mstmp40 + Mstmp33*Mstmp41 + Mstmp33*Mstmp42 + Mstmp33*M[15] + Mstmp43*M[13] + Mstmp77*z + Mstmp78*z + Mstmp79*z + Mstmp80*y + x*M[32] + y*M[28] + z*M[27] + M[47];
#pragma omp atomic
Ms[48] += Mstmp34*Mstmp43 + Mstmp35*Mstmp43 + Mstmp37*Mstmp43 + Mstmp43*M[14] + Mstmp6*Mstmp84 + Mstmp7*Mstmp84 + Mstmp8*Mstmp84 + Mstmp80*z + Mstmp81*z + Mstmp82*z + Mstmp83*y + Mstmp84*M[5] + x*M[33] + y*M[29] + z*M[28] + M[48];
#pragma omp atomic
Ms[49] += Mstmp0*Mstmp103 + Mstmp100*Mstmp91 + Mstmp40*Mstmp43 + Mstmp43*M[15] + Mstmp83*z + Mstmp84*Mstmp9 + Mstmp84*M[6] + x*M[34] + z*M[29] + M[49];
#pragma omp atomic
Ms[50] += Mstmp101*M[2] + Mstmp33*M[16] + Mstmp76*M[7] + Mstmp93*(y*y*y*y*y) + y*M[30] + M[50];
#pragma omp atomic
Ms[51] += Mstmp101*Mstmp2 + Mstmp101*M[3] + Mstmp15*Mstmp76 + Mstmp33*Mstmp47 + Mstmp33*M[17] + Mstmp76*M[8] + Mstmp85*z + y*M[31] + z*M[30] + M[51];
#pragma omp atomic
Ms[52] += Mstmp102*M[2] + Mstmp17*Mstmp76 + Mstmp33*Mstmp50 + Mstmp33*M[18] + Mstmp43*Mstmp44 + Mstmp43*M[16] + Mstmp76*M[9] + Mstmp87*z + Mstmp96*Mstmp97 + y*M[32] + z*M[31] + M[52];
#pragma omp atomic
Ms[53] += Mstmp102*M[3] + Mstmp12*Mstmp84 + Mstmp33*Mstmp52 + Mstmp33*M[19] + Mstmp43*Mstmp46 + Mstmp43*M[17] + Mstmp65*Mstmp99 + Mstmp84*M[7] + Mstmp88*z + y*M[33] + z*M[32] + M[53];
#pragma omp atomic
Ms[54] += Mstmp1*Mstmp103 + Mstmp103*M[2] + Mstmp14*Mstmp84 + Mstmp43*Mstmp49 + Mstmp43*M[18] + Mstmp84*M[8] + Mstmp90*z + y*M[34] + z*M[33] + M[54];
#pragma omp atomic
Ms[55] += Mstmp103*M[3] + Mstmp43*M[19] + Mstmp84*M[9] + Mstmp93*(z*z*z*z*z) + z*M[34] + M[55];

}

void field_m0_M2L_5(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double D[56];
double Dtmp0 = 1.0*pow(R, -3.0);
double Dtmp1 = -Dtmp0;
double Dtmp2 = (x*x);
double Dtmp3 = pow(R, -5.0);
double Dtmp4 = 3.0*Dtmp3;
double Dtmp5 = Dtmp4*x;
double Dtmp6 = (y*y);
double Dtmp7 = y*z;
double Dtmp8 = 9.0*Dtmp3;
double Dtmp9 = -Dtmp8;
double Dtmp10 = pow(R, -7.0);
double Dtmp11 = 15.0*Dtmp10;
double Dtmp12 = Dtmp11*Dtmp2;
double Dtmp13 = -Dtmp4;
double Dtmp14 = Dtmp12 + Dtmp13;
double Dtmp15 = Dtmp11*Dtmp6;
double Dtmp16 = Dtmp13 + Dtmp15;
double Dtmp17 = 1.0*x;
double Dtmp18 = Dtmp7*x;
double Dtmp19 = (x*x*x*x);
double Dtmp20 = pow(R, -9.0);
double Dtmp21 = 105.0*Dtmp20;
double Dtmp22 = 90.0*Dtmp10;
double Dtmp23 = 45.0*Dtmp10;
double Dtmp24 = -Dtmp23;
double Dtmp25 = Dtmp2*Dtmp21;
double Dtmp26 = x*(Dtmp24 + Dtmp25);
double Dtmp27 = -Dtmp11;
double Dtmp28 = Dtmp21*Dtmp6;
double Dtmp29 = Dtmp24 + Dtmp28;
double Dtmp30 = (y*y*y*y);
double Dtmp31 = 225.0*Dtmp10;
double Dtmp32 = 945.0*pow(R, -11.0);
double Dtmp33 = Dtmp19*Dtmp32;
double Dtmp34 = Dtmp2*Dtmp20;
double Dtmp35 = Dtmp23 + Dtmp33 - 630.0*Dtmp34;
double Dtmp36 = -Dtmp25;
double Dtmp37 = 315.0*Dtmp20;
double Dtmp38 = Dtmp2*Dtmp32;
double Dtmp39 = Dtmp38*Dtmp6;
double Dtmp40 = Dtmp23 + Dtmp39;
double Dtmp41 = -Dtmp37;
double Dtmp42 = -Dtmp28;
double Dtmp43 = Dtmp30*Dtmp32;
double Dtmp44 = Dtmp20*Dtmp6;
double Dtmp45 = Dtmp23 + Dtmp43 - 630.0*Dtmp44;
D[0] = 1.0/R;
D[1] = -Dtmp0*x;
D[2] = -Dtmp0*y;
D[3] = -Dtmp0*z;
D[4] = Dtmp1 + Dtmp2*Dtmp4;
D[5] = Dtmp5*y;
D[6] = Dtmp5*z;
D[7] = Dtmp1 + Dtmp4*Dtmp6;
D[8] = Dtmp4*Dtmp7;
D[9] = -D[4] - D[7];
D[10] = -x*(Dtmp12 + Dtmp9);
D[11] = -Dtmp14*y;
D[12] = -Dtmp14*z;
D[13] = -Dtmp16*Dtmp17;
D[14] = -Dtmp11*Dtmp18;
D[15] = -D[10] - D[13];
D[16] = -y*(Dtmp15 + Dtmp9);
D[17] = -Dtmp16*z;
D[18] = -D[11] - D[16];
D[19] = -D[12] - D[17];
D[20] = Dtmp19*Dtmp21 - Dtmp2*Dtmp22 + Dtmp8;
D[21] = Dtmp26*y;
D[22] = Dtmp26*z;
D[23] = -Dtmp12 - Dtmp15 + Dtmp25*Dtmp6 + Dtmp4;
D[24] = Dtmp7*(Dtmp25 + Dtmp27);
D[25] = -D[20] - D[23];
D[26] = Dtmp17*Dtmp29*y;
D[27] = Dtmp17*z*(Dtmp27 + Dtmp28);
D[28] = -D[21] - D[26];
D[29] = -D[22] - D[27];
D[30] = Dtmp21*Dtmp30 - Dtmp22*Dtmp6 + Dtmp8;
D[31] = Dtmp29*Dtmp7;
D[32] = -D[23] - D[30];
D[33] = -D[24] - D[31];
D[34] = -D[25] - D[32];
D[35] = -x*(Dtmp31 + Dtmp33 - 1050.0*Dtmp34);
D[36] = -Dtmp35*y;
D[37] = -Dtmp35*z;
D[38] = -x*(Dtmp36 - Dtmp37*Dtmp6 + Dtmp40);
D[39] = -Dtmp18*(Dtmp38 + Dtmp41);
D[40] = -D[35] - D[38];
D[41] = -y*(-Dtmp2*Dtmp37 + Dtmp40 + Dtmp42);
D[42] = -z*(Dtmp11 + Dtmp36 + Dtmp39 + Dtmp42);
D[43] = -D[36] - D[41];
D[44] = -D[37] - D[42];
D[45] = -Dtmp17*Dtmp45;
D[46] = -Dtmp17*Dtmp7*(Dtmp32*Dtmp6 + Dtmp41);
D[47] = -D[38] - D[45];
D[48] = -D[39] - D[46];
D[49] = -D[40] - D[47];
D[50] = -y*(Dtmp31 + Dtmp43 - 1050.0*Dtmp44);
D[51] = -Dtmp45*z;
D[52] = -D[41] - D[50];
D[53] = -D[42] - D[51];
D[54] = -D[43] - D[52];
D[55] = -D[44] - D[53];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[8]*M[8] + D[9]*M[9] + D[10]*M[10] + D[11]*M[11] + D[12]*M[12] + D[13]*M[13] + D[14]*M[14] + D[15]*M[15] + D[16]*M[16] + D[17]*M[17] + D[18]*M[18] + D[19]*M[19] + D[20]*M[20] + D[21]*M[21] + D[22]*M[22] + D[23]*M[23] + D[24]*M[24] + D[25]*M[25] + D[26]*M[26] + D[27]*M[27] + D[28]*M[28] + D[29]*M[29] + D[30]*M[30] + D[31]*M[31] + D[32]*M[32] + D[33]*M[33] + D[34]*M[34] + D[35]*M[35] + D[36]*M[36] + D[37]*M[37] + D[38]*M[38] + D[39]*M[39] + D[40]*M[40] + D[41]*M[41] + D[42]*M[42] + D[43]*M[43] + D[44]*M[44] + D[45]*M[45] + D[46]*M[46] + D[47]*M[47] + D[48]*M[48] + D[49]*M[49] + D[50]*M[50] + D[51]*M[51] + D[52]*M[52] + D[53]*M[53] + D[54]*M[54] + D[55]*M[55];
#pragma omp atomic
L[1] += D[1]*M[0] + D[4]*M[1] + D[5]*M[2] + D[6]*M[3] + D[10]*M[4] + D[11]*M[5] + D[12]*M[6] + D[13]*M[7] + D[14]*M[8] + D[15]*M[9] + D[20]*M[10] + D[21]*M[11] + D[22]*M[12] + D[23]*M[13] + D[24]*M[14] + D[25]*M[15] + D[26]*M[16] + D[27]*M[17] + D[28]*M[18] + D[29]*M[19] + D[35]*M[20] + D[36]*M[21] + D[37]*M[22] + D[38]*M[23] + D[39]*M[24] + D[40]*M[25] + D[41]*M[26] + D[42]*M[27] + D[43]*M[28] + D[44]*M[29] + D[45]*M[30] + D[46]*M[31] + D[47]*M[32] + D[48]*M[33] + D[49]*M[34];
#pragma omp atomic
L[2] += D[2]*M[0] + D[5]*M[1] + D[7]*M[2] + D[8]*M[3] + D[11]*M[4] + D[13]*M[5] + D[14]*M[6] + D[16]*M[7] + D[17]*M[8] + D[18]*M[9] + D[21]*M[10] + D[23]*M[11] + D[24]*M[12] + D[26]*M[13] + D[27]*M[14] + D[28]*M[15] + D[30]*M[16] + D[31]*M[17] + D[32]*M[18] + D[33]*M[19] + D[36]*M[20] + D[38]*M[21] + D[39]*M[22] + D[41]*M[23] + D[42]*M[24] + D[43]*M[25] + D[45]*M[26] + D[46]*M[27] + D[47]*M[28] + D[48]*M[29] + D[50]*M[30] + D[51]*M[31] + D[52]*M[32] + D[53]*M[33] + D[54]*M[34];
#pragma omp atomic
L[3] += D[3]*M[0] + D[6]*M[1] + D[8]*M[2] + D[9]*M[3] + D[12]*M[4] + D[14]*M[5] + D[15]*M[6] + D[17]*M[7] + D[18]*M[8] + D[19]*M[9] + D[22]*M[10] + D[24]*M[11] + D[25]*M[12] + D[27]*M[13] + D[28]*M[14] + D[29]*M[15] + D[31]*M[16] + D[32]*M[17] + D[33]*M[18] + D[34]*M[19] + D[37]*M[20] + D[39]*M[21] + D[40]*M[22] + D[42]*M[23] + D[43]*M[24] + D[44]*M[25] + D[46]*M[26] + D[47]*M[27] + D[48]*M[28] + D[49]*M[29] + D[51]*M[30] + D[52]*M[31] + D[53]*M[32] + D[54]*M[33] + D[55]*M[34];
#pragma omp atomic
L[4] += D[4]*M[0] + D[10]*M[1] + D[11]*M[2] + D[12]*M[3] + D[20]*M[4] + D[21]*M[5] + D[22]*M[6] + D[23]*M[7] + D[24]*M[8] + D[25]*M[9] + D[35]*M[10] + D[36]*M[11] + D[37]*M[12] + D[38]*M[13] + D[39]*M[14] + D[40]*M[15] + D[41]*M[16] + D[42]*M[17] + D[43]*M[18] + D[44]*M[19];
#pragma omp atomic
L[5] += D[5]*M[0] + D[11]*M[1] + D[13]*M[2] + D[14]*M[3] + D[21]*M[4] + D[23]*M[5] + D[24]*M[6] + D[26]*M[7] + D[27]*M[8] + D[28]*M[9] + D[36]*M[10] + D[38]*M[11] + D[39]*M[12] + D[41]*M[13] + D[42]*M[14] + D[43]*M[15] + D[45]*M[16] + D[46]*M[17] + D[47]*M[18] + D[48]*M[19];
#pragma omp atomic
L[6] += D[6]*M[0] + D[12]*M[1] + D[14]*M[2] + D[15]*M[3] + D[22]*M[4] + D[24]*M[5] + D[25]*M[6] + D[27]*M[7] + D[28]*M[8] + D[29]*M[9] + D[37]*M[10] + D[39]*M[11] + D[40]*M[12] + D[42]*M[13] + D[43]*M[14] + D[44]*M[15] + D[46]*M[16] + D[47]*M[17] + D[48]*M[18] + D[49]*M[19];
#pragma omp atomic
L[7] += D[7]*M[0] + D[13]*M[1] + D[16]*M[2] + D[17]*M[3] + D[23]*M[4] + D[26]*M[5] + D[27]*M[6] + D[30]*M[7] + D[31]*M[8] + D[32]*M[9] + D[38]*M[10] + D[41]*M[11] + D[42]*M[12] + D[45]*M[13] + D[46]*M[14] + D[47]*M[15] + D[50]*M[16] + D[51]*M[17] + D[52]*M[18] + D[53]*M[19];
#pragma omp atomic
L[8] += D[8]*M[0] + D[14]*M[1] + D[17]*M[2] + D[18]*M[3] + D[24]*M[4] + D[27]*M[5] + D[28]*M[6] + D[31]*M[7] + D[32]*M[8] + D[33]*M[9] + D[39]*M[10] + D[42]*M[11] + D[43]*M[12] + D[46]*M[13] + D[47]*M[14] + D[48]*M[15] + D[51]*M[16] + D[52]*M[17] + D[53]*M[18] + D[54]*M[19];
#pragma omp atomic
L[9] += D[9]*M[0] + D[15]*M[1] + D[18]*M[2] + D[19]*M[3] + D[25]*M[4] + D[28]*M[5] + D[29]*M[6] + D[32]*M[7] + D[33]*M[8] + D[34]*M[9] + D[40]*M[10] + D[43]*M[11] + D[44]*M[12] + D[47]*M[13] + D[48]*M[14] + D[49]*M[15] + D[52]*M[16] + D[53]*M[17] + D[54]*M[18] + D[55]*M[19];
#pragma omp atomic
L[10] += D[10]*M[0] + D[20]*M[1] + D[21]*M[2] + D[22]*M[3] + D[35]*M[4] + D[36]*M[5] + D[37]*M[6] + D[38]*M[7] + D[39]*M[8] + D[40]*M[9];
#pragma omp atomic
L[11] += D[11]*M[0] + D[21]*M[1] + D[23]*M[2] + D[24]*M[3] + D[36]*M[4] + D[38]*M[5] + D[39]*M[6] + D[41]*M[7] + D[42]*M[8] + D[43]*M[9];
#pragma omp atomic
L[12] += D[12]*M[0] + D[22]*M[1] + D[24]*M[2] + D[25]*M[3] + D[37]*M[4] + D[39]*M[5] + D[40]*M[6] + D[42]*M[7] + D[43]*M[8] + D[44]*M[9];
#pragma omp atomic
L[13] += D[13]*M[0] + D[23]*M[1] + D[26]*M[2] + D[27]*M[3] + D[38]*M[4] + D[41]*M[5] + D[42]*M[6] + D[45]*M[7] + D[46]*M[8] + D[47]*M[9];
#pragma omp atomic
L[14] += D[14]*M[0] + D[24]*M[1] + D[27]*M[2] + D[28]*M[3] + D[39]*M[4] + D[42]*M[5] + D[43]*M[6] + D[46]*M[7] + D[47]*M[8] + D[48]*M[9];
#pragma omp atomic
L[15] += D[15]*M[0] + D[25]*M[1] + D[28]*M[2] + D[29]*M[3] + D[40]*M[4] + D[43]*M[5] + D[44]*M[6] + D[47]*M[7] + D[48]*M[8] + D[49]*M[9];
#pragma omp atomic
L[16] += D[16]*M[0] + D[26]*M[1] + D[30]*M[2] + D[31]*M[3] + D[41]*M[4] + D[45]*M[5] + D[46]*M[6] + D[50]*M[7] + D[51]*M[8] + D[52]*M[9];
#pragma omp atomic
L[17] += D[17]*M[0] + D[27]*M[1] + D[31]*M[2] + D[32]*M[3] + D[42]*M[4] + D[46]*M[5] + D[47]*M[6] + D[51]*M[7] + D[52]*M[8] + D[53]*M[9];
#pragma omp atomic
L[18] += D[18]*M[0] + D[28]*M[1] + D[32]*M[2] + D[33]*M[3] + D[43]*M[4] + D[47]*M[5] + D[48]*M[6] + D[52]*M[7] + D[53]*M[8] + D[54]*M[9];
#pragma omp atomic
L[19] += D[19]*M[0] + D[29]*M[1] + D[33]*M[2] + D[34]*M[3] + D[44]*M[4] + D[48]*M[5] + D[49]*M[6] + D[53]*M[7] + D[54]*M[8] + D[55]*M[9];
#pragma omp atomic
L[20] += D[20]*M[0] + D[35]*M[1] + D[36]*M[2] + D[37]*M[3];
#pragma omp atomic
L[21] += D[21]*M[0] + D[36]*M[1] + D[38]*M[2] + D[39]*M[3];
#pragma omp atomic
L[22] += D[22]*M[0] + D[37]*M[1] + D[39]*M[2] + D[40]*M[3];
#pragma omp atomic
L[23] += D[23]*M[0] + D[38]*M[1] + D[41]*M[2] + D[42]*M[3];
#pragma omp atomic
L[24] += D[24]*M[0] + D[39]*M[1] + D[42]*M[2] + D[43]*M[3];
#pragma omp atomic
L[25] += D[25]*M[0] + D[40]*M[1] + D[43]*M[2] + D[44]*M[3];
#pragma omp atomic
L[26] += D[26]*M[0] + D[41]*M[1] + D[45]*M[2] + D[46]*M[3];
#pragma omp atomic
L[27] += D[27]*M[0] + D[42]*M[1] + D[46]*M[2] + D[47]*M[3];
#pragma omp atomic
L[28] += D[28]*M[0] + D[43]*M[1] + D[47]*M[2] + D[48]*M[3];
#pragma omp atomic
L[29] += D[29]*M[0] + D[44]*M[1] + D[48]*M[2] + D[49]*M[3];
#pragma omp atomic
L[30] += D[30]*M[0] + D[45]*M[1] + D[50]*M[2] + D[51]*M[3];
#pragma omp atomic
L[31] += D[31]*M[0] + D[46]*M[1] + D[51]*M[2] + D[52]*M[3];
#pragma omp atomic
L[32] += D[32]*M[0] + D[47]*M[1] + D[52]*M[2] + D[53]*M[3];
#pragma omp atomic
L[33] += D[33]*M[0] + D[48]*M[1] + D[53]*M[2] + D[54]*M[3];
#pragma omp atomic
L[34] += D[34]*M[0] + D[49]*M[1] + D[54]*M[2] + D[55]*M[3];
#pragma omp atomic
L[35] += D[35]*M[0];
#pragma omp atomic
L[36] += D[36]*M[0];
#pragma omp atomic
L[37] += D[37]*M[0];
#pragma omp atomic
L[38] += D[38]*M[0];
#pragma omp atomic
L[39] += D[39]*M[0];
#pragma omp atomic
L[40] += D[40]*M[0];
#pragma omp atomic
L[41] += D[41]*M[0];
#pragma omp atomic
L[42] += D[42]*M[0];
#pragma omp atomic
L[43] += D[43]*M[0];
#pragma omp atomic
L[44] += D[44]*M[0];
#pragma omp atomic
L[45] += D[45]*M[0];
#pragma omp atomic
L[46] += D[46]*M[0];
#pragma omp atomic
L[47] += D[47]*M[0];
#pragma omp atomic
L[48] += D[48]*M[0];
#pragma omp atomic
L[49] += D[49]*M[0];
#pragma omp atomic
L[50] += D[50]*M[0];
#pragma omp atomic
L[51] += D[51]*M[0];
#pragma omp atomic
L[52] += D[52]*M[0];
#pragma omp atomic
L[53] += D[53]*M[0];
#pragma omp atomic
L[54] += D[54]*M[0];
#pragma omp atomic
L[55] += D[55]*M[0];

}

void field_m0_L2L_5(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
double Lstmp3 = z*L[14];
double Lstmp4 = Lstmp3*y;
double Lstmp5 = (x*x);
double Lstmp6 = (1.0/2.0)*Lstmp5;
double Lstmp7 = (x*x*x);
double Lstmp8 = (1.0/6.0)*Lstmp7;
double Lstmp9 = (1.0/24.0)*(x*x*x*x);
double Lstmp10 = (y*y);
double Lstmp11 = (1.0/2.0)*Lstmp10;
double Lstmp12 = (y*y*y);
double Lstmp13 = (1.0/6.0)*Lstmp12;
double Lstmp14 = (1.0/24.0)*(y*y*y*y);
double Lstmp15 = (z*z);
double Lstmp16 = (1.0/2.0)*Lstmp15;
double Lstmp17 = (z*z*z);
double Lstmp18 = (1.0/6.0)*Lstmp17;
double Lstmp19 = (1.0/24.0)*(z*z*z*z);
double Lstmp20 = x*L[13];
double Lstmp21 = x*L[26];
double Lstmp22 = x*L[45];
double Lstmp23 = x*L[15];
double Lstmp24 = x*L[29];
double Lstmp25 = x*L[49];
double Lstmp26 = y*L[11];
double Lstmp27 = z*L[12];
double Lstmp28 = y*L[21];
double Lstmp29 = z*L[22];
double Lstmp30 = y*L[36];
double Lstmp31 = z*L[37];
double Lstmp32 = y*L[18];
double Lstmp33 = y*L[33];
double Lstmp34 = y*L[54];
double Lstmp35 = z*L[17];
double Lstmp36 = z*L[31];
double Lstmp37 = z*L[51];
double Lstmp38 = y*L[28];
double Lstmp39 = Lstmp38*x;
double Lstmp40 = y*L[48];
double Lstmp41 = Lstmp40*x;
double Lstmp42 = z*L[27];
double Lstmp43 = Lstmp42*x;
double Lstmp44 = z*L[46];
double Lstmp45 = Lstmp44*x;
double Lstmp46 = z*L[24];
double Lstmp47 = Lstmp46*y;
double Lstmp48 = z*L[39];
double Lstmp49 = Lstmp48*y;
double Lstmp50 = (1.0/4.0)*Lstmp5;
double Lstmp51 = Lstmp10*Lstmp50;
double Lstmp52 = (1.0/12.0)*Lstmp5;
double Lstmp53 = Lstmp15*Lstmp50;
double Lstmp54 = (1.0/12.0)*Lstmp7;
double Lstmp55 = (1.0/4.0)*Lstmp10*Lstmp15;
double Lstmp56 = x*L[47];
double Lstmp57 = y*L[43];
double Lstmp58 = z*L[42];
double Lstmp59 = x*L[23];
double Lstmp60 = x*L[41];
double Lstmp61 = x*L[25];
double Lstmp62 = x*L[44];
double Lstmp63 = Lstmp57*x;
double Lstmp64 = Lstmp58*x;
double Lstmp65 = y*L[13];
double Lstmp66 = Lstmp42*y;
double Lstmp67 = x*L[28];
double Lstmp68 = x*L[48];
double Lstmp69 = y*L[23];
double Lstmp70 = y*L[38];
double Lstmp71 = y*L[32];
double Lstmp72 = y*L[53];
double Lstmp73 = y*L[47];
double Lstmp74 = Lstmp73*x;
double Lstmp75 = Lstmp58*y;
double Lstmp76 = y*L[14];
double Lstmp77 = z*L[15];
double Lstmp78 = z*L[18];
double Lstmp79 = z*L[28];
double Lstmp80 = Lstmp79*y;
double Lstmp81 = x*L[27];
double Lstmp82 = x*L[46];
double Lstmp83 = y*L[24];
double Lstmp84 = z*L[25];
double Lstmp85 = y*L[39];
double Lstmp86 = z*L[40];
double Lstmp87 = z*L[32];
double Lstmp88 = z*L[52];
double Lstmp89 = z*L[47];
double Lstmp90 = Lstmp89*x;
double Lstmp91 = z*L[43];
double Lstmp92 = Lstmp91*y;
double Lstmp93 = x*L[38];
double Lstmp94 = x*L[40];
double Lstmp95 = x*L[43];
double Lstmp96 = x*L[42];
double Lstmp97 = y*L[26];
double Lstmp98 = Lstmp44*y;
double Lstmp99 = y*L[41];
double Lstmp100 = y*L[52];
double Lstmp101 = y*L[27];
double Lstmp102 = Lstmp89*y;
double Lstmp103 = y*L[42];
double Lstmp104 = z*L[29];
double Lstmp105 = z*L[33];
double Lstmp106 = z*L[48];
double Lstmp107 = Lstmp106*y;
double Lstmp108 = z*L[44];
double Lstmp109 = z*L[53];
double Lstmp110 = y*L[45];
double Lstmp111 = y*L[46];
double Lstmp112 = z*L[49];
double Lstmp113 = z*L[54];
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x + (1.0/12.0)*Lstmp10*Lstmp17*L[53] + Lstmp10*Lstmp54*L[38] + Lstmp11*Lstmp20 + Lstmp11*Lstmp35 + Lstmp11*Lstmp43 + Lstmp11*L[7] + (1.0/12.0)*Lstmp12*Lstmp15*L[52] + Lstmp12*Lstmp52*L[41] + Lstmp13*Lstmp21 + Lstmp13*Lstmp36 + Lstmp13*Lstmp45 + Lstmp13*L[16] + Lstmp14*Lstmp22 + Lstmp14*Lstmp37 + Lstmp14*L[30] + Lstmp15*Lstmp54*L[40] + Lstmp16*Lstmp23 + Lstmp16*Lstmp32 + Lstmp16*Lstmp39 + Lstmp16*L[9] + Lstmp17*Lstmp52*L[44] + Lstmp18*Lstmp24 + Lstmp18*Lstmp33 + Lstmp18*Lstmp41 + Lstmp18*L[19] + Lstmp19*Lstmp25 + Lstmp19*Lstmp34 + Lstmp19*L[34] + Lstmp2*y + Lstmp26*Lstmp6 + Lstmp27*Lstmp6 + Lstmp28*Lstmp8 + Lstmp29*Lstmp8 + Lstmp30*Lstmp9 + Lstmp31*Lstmp9 + Lstmp4*x + Lstmp47*Lstmp6 + Lstmp49*Lstmp8 + Lstmp51*Lstmp58 + Lstmp51*L[23] + Lstmp53*Lstmp57 + Lstmp53*L[25] + Lstmp55*Lstmp56 + Lstmp55*L[32] + Lstmp6*L[4] + Lstmp8*L[10] + Lstmp9*L[20] + (1.0/120.0)*(x*x*x*x*x)*L[35] + x*L[1] + (1.0/120.0)*(y*y*y*y*y)*L[50] + y*L[2] + (1.0/120.0)*(z*z*z*z*z)*L[55] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 + Lstmp11*Lstmp42 + Lstmp11*Lstmp59 + Lstmp11*Lstmp64 + Lstmp11*L[13] + Lstmp13*Lstmp44 + Lstmp13*Lstmp60 + Lstmp13*L[26] + Lstmp14*L[45] + Lstmp16*Lstmp38 + Lstmp16*Lstmp61 + Lstmp16*Lstmp63 + Lstmp16*L[15] + Lstmp18*Lstmp40 + Lstmp18*Lstmp62 + Lstmp18*L[29] + Lstmp19*L[49] + Lstmp26*x + Lstmp27*x + Lstmp28*Lstmp6 + Lstmp29*Lstmp6 + Lstmp30*Lstmp8 + Lstmp31*Lstmp8 + Lstmp4 + Lstmp47*x + Lstmp49*Lstmp6 + Lstmp51*L[38] + Lstmp53*L[40] + Lstmp55*L[47] + Lstmp6*L[10] + Lstmp8*L[20] + Lstmp9*L[35] + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp11*Lstmp21 + Lstmp11*Lstmp36 + Lstmp11*Lstmp45 + Lstmp11*L[16] + Lstmp13*Lstmp22 + Lstmp13*Lstmp37 + Lstmp13*L[30] + Lstmp14*L[50] + Lstmp16*Lstmp67 + Lstmp16*Lstmp71 + Lstmp16*Lstmp74 + Lstmp16*L[18] + Lstmp18*Lstmp68 + Lstmp18*Lstmp72 + Lstmp18*L[33] + Lstmp19*L[54] + Lstmp2 + Lstmp3*x + Lstmp35*y + Lstmp46*Lstmp6 + Lstmp48*Lstmp8 + Lstmp51*L[41] + Lstmp53*L[43] + Lstmp55*L[52] + Lstmp6*Lstmp69 + Lstmp6*Lstmp75 + Lstmp6*L[11] + Lstmp65*x + Lstmp66*x + Lstmp70*Lstmp8 + Lstmp8*L[21] + Lstmp9*L[36] + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += Lstmp11*Lstmp81 + Lstmp11*Lstmp87 + Lstmp11*Lstmp90 + Lstmp11*L[17] + Lstmp13*Lstmp82 + Lstmp13*Lstmp88 + Lstmp13*L[31] + Lstmp14*L[51] + Lstmp16*Lstmp24 + Lstmp16*Lstmp33 + Lstmp16*Lstmp41 + Lstmp16*L[19] + Lstmp18*Lstmp25 + Lstmp18*Lstmp34 + Lstmp18*L[34] + Lstmp19*L[55] + Lstmp51*L[42] + Lstmp53*L[44] + Lstmp55*L[53] + Lstmp6*Lstmp83 + Lstmp6*Lstmp84 + Lstmp6*Lstmp92 + Lstmp6*L[12] + Lstmp76*x + Lstmp77*x + Lstmp78*y + Lstmp8*Lstmp85 + Lstmp8*Lstmp86 + Lstmp8*L[22] + Lstmp80*x + Lstmp9*L[37] + x*L[6] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += Lstmp11*Lstmp58 + Lstmp11*Lstmp93 + Lstmp11*L[23] + Lstmp13*L[41] + Lstmp16*Lstmp57 + Lstmp16*Lstmp94 + Lstmp16*L[25] + Lstmp18*L[44] + Lstmp26 + Lstmp27 + Lstmp28*x + Lstmp29*x + Lstmp30*Lstmp6 + Lstmp31*Lstmp6 + Lstmp47 + Lstmp49*x + Lstmp6*L[20] + Lstmp8*L[35] + x*L[10] + L[4];
#pragma omp atomic
Ls[5] += Lstmp11*Lstmp44 + Lstmp11*Lstmp60 + Lstmp11*L[26] + Lstmp13*L[45] + Lstmp16*Lstmp73 + Lstmp16*Lstmp95 + Lstmp16*L[28] + Lstmp18*L[48] + Lstmp3 + Lstmp46*x + Lstmp48*Lstmp6 + Lstmp6*Lstmp70 + Lstmp6*L[21] + Lstmp65 + Lstmp66 + Lstmp69*x + Lstmp75*x + Lstmp8*L[36] + x*L[11] + L[5];
#pragma omp atomic
Ls[6] += Lstmp11*Lstmp89 + Lstmp11*Lstmp96 + Lstmp11*L[27] + Lstmp13*L[46] + Lstmp16*Lstmp40 + Lstmp16*Lstmp62 + Lstmp16*L[29] + Lstmp18*L[49] + Lstmp6*Lstmp85 + Lstmp6*Lstmp86 + Lstmp6*L[22] + Lstmp76 + Lstmp77 + Lstmp8*L[37] + Lstmp80 + Lstmp83*x + Lstmp84*x + Lstmp92*x + x*L[12] + L[6];
#pragma omp atomic
Ls[7] += Lstmp100*Lstmp16 + Lstmp11*Lstmp22 + Lstmp11*Lstmp37 + Lstmp11*L[30] + Lstmp13*L[50] + Lstmp16*Lstmp56 + Lstmp16*L[32] + Lstmp18*L[53] + Lstmp20 + Lstmp35 + Lstmp36*y + Lstmp43 + Lstmp58*Lstmp6 + Lstmp6*Lstmp99 + Lstmp6*L[23] + Lstmp8*L[38] + Lstmp97*x + Lstmp98*x + y*L[16] + L[7];
#pragma omp atomic
Ls[8] += Lstmp101*x + Lstmp102*x + Lstmp103*Lstmp6 + Lstmp11*Lstmp82 + Lstmp11*Lstmp88 + Lstmp11*L[31] + Lstmp13*L[51] + Lstmp16*Lstmp68 + Lstmp16*Lstmp72 + Lstmp16*L[33] + Lstmp18*L[54] + Lstmp6*Lstmp91 + Lstmp6*L[24] + Lstmp78 + Lstmp79*x + Lstmp8*L[39] + Lstmp87*y + x*L[14] + y*L[17] + L[8];
#pragma omp atomic
Ls[9] += Lstmp104*x + Lstmp105*y + Lstmp107*x + Lstmp108*Lstmp6 + Lstmp109*Lstmp11 + Lstmp11*Lstmp56 + Lstmp11*L[32] + Lstmp13*L[52] + Lstmp16*Lstmp25 + Lstmp16*Lstmp34 + Lstmp16*L[34] + Lstmp18*L[55] + Lstmp23 + Lstmp32 + Lstmp39 + Lstmp57*Lstmp6 + Lstmp6*L[25] + Lstmp8*L[40] + z*L[19] + L[9];
#pragma omp atomic
Ls[10] += Lstmp11*L[38] + Lstmp16*L[40] + Lstmp28 + Lstmp29 + Lstmp30*x + Lstmp31*x + Lstmp49 + Lstmp6*L[35] + x*L[20] + L[10];
#pragma omp atomic
Ls[11] += Lstmp11*L[41] + Lstmp16*L[43] + Lstmp46 + Lstmp48*x + Lstmp6*L[36] + Lstmp69 + Lstmp70*x + Lstmp75 + x*L[21] + L[11];
#pragma omp atomic
Ls[12] += Lstmp11*L[42] + Lstmp16*L[44] + Lstmp6*L[37] + Lstmp83 + Lstmp84 + Lstmp85*x + Lstmp86*x + Lstmp92 + x*L[22] + L[12];
#pragma omp atomic
Ls[13] += Lstmp11*L[45] + Lstmp16*L[47] + Lstmp42 + Lstmp59 + Lstmp6*L[38] + Lstmp64 + Lstmp97 + Lstmp98 + Lstmp99*x + L[13];
#pragma omp atomic
Ls[14] += Lstmp101 + Lstmp102 + Lstmp103*x + Lstmp11*L[46] + Lstmp16*L[48] + Lstmp6*L[39] + Lstmp79 + Lstmp91*x + x*L[24] + L[14];
#pragma omp atomic
Ls[15] += Lstmp104 + Lstmp107 + Lstmp108*x + Lstmp11*L[47] + Lstmp16*L[49] + Lstmp38 + Lstmp6*L[40] + Lstmp61 + Lstmp63 + L[15];
#pragma omp atomic
Ls[16] += Lstmp11*L[50] + Lstmp110*x + Lstmp16*L[52] + Lstmp21 + Lstmp36 + Lstmp37*y + Lstmp45 + Lstmp6*L[41] + y*L[30] + L[16];
#pragma omp atomic
Ls[17] += Lstmp11*L[51] + Lstmp111*x + Lstmp16*L[53] + Lstmp6*L[42] + Lstmp81 + Lstmp87 + Lstmp88*y + Lstmp90 + y*L[31] + L[17];
#pragma omp atomic
Ls[18] += Lstmp105 + Lstmp106*x + Lstmp109*y + Lstmp11*L[52] + Lstmp16*L[54] + Lstmp6*L[43] + Lstmp67 + Lstmp71 + Lstmp74 + L[18];
#pragma omp atomic
Ls[19] += Lstmp11*L[53] + Lstmp112*x + Lstmp113*y + Lstmp16*L[55] + Lstmp24 + Lstmp33 + Lstmp41 + Lstmp6*L[44] + z*L[34] + L[19];
#pragma omp atomic
Ls[20] += Lstmp30 + Lstmp31 + x*L[35] + L[20];
#pragma omp atomic
Ls[21] += Lstmp48 + Lstmp70 + x*L[36] + L[21];
#pragma omp atomic
Ls[22] += Lstmp85 + Lstmp86 + x*L[37] + L[22];
#pragma omp atomic
Ls[23] += Lstmp58 + Lstmp93 + Lstmp99 + L[23];
#pragma omp atomic
Ls[24] += Lstmp103 + Lstmp91 + x*L[39] + L[24];
#pragma omp atomic
Ls[25] += Lstmp108 + Lstmp57 + Lstmp94 + L[25];
#pragma omp atomic
Ls[26] += Lstmp110 + Lstmp44 + Lstmp60 + L[26];
#pragma omp atomic
Ls[27] += Lstmp111 + Lstmp89 + Lstmp96 + L[27];
#pragma omp atomic
Ls[28] += Lstmp106 + Lstmp73 + Lstmp95 + L[28];
#pragma omp atomic
Ls[29] += Lstmp112 + Lstmp40 + Lstmp62 + L[29];
#pragma omp atomic
Ls[30] += Lstmp22 + Lstmp37 + y*L[50] + L[30];
#pragma omp atomic
Ls[31] += Lstmp82 + Lstmp88 + y*L[51] + L[31];
#pragma omp atomic
Ls[32] += Lstmp100 + Lstmp109 + Lstmp56 + L[32];
#pragma omp atomic
Ls[33] += Lstmp113 + Lstmp68 + Lstmp72 + L[33];
#pragma omp atomic
Ls[34] += Lstmp25 + Lstmp34 + z*L[55] + L[34];
#pragma omp atomic
Ls[35] += L[35];
#pragma omp atomic
Ls[36] += L[36];
#pragma omp atomic
Ls[37] += L[37];
#pragma omp atomic
Ls[38] += L[38];
#pragma omp atomic
Ls[39] += L[39];
#pragma omp atomic
Ls[40] += L[40];
#pragma omp atomic
Ls[41] += L[41];
#pragma omp atomic
Ls[42] += L[42];
#pragma omp atomic
Ls[43] += L[43];
#pragma omp atomic
Ls[44] += L[44];
#pragma omp atomic
Ls[45] += L[45];
#pragma omp atomic
Ls[46] += L[46];
#pragma omp atomic
Ls[47] += L[47];
#pragma omp atomic
Ls[48] += L[48];
#pragma omp atomic
Ls[49] += L[49];
#pragma omp atomic
Ls[50] += L[50];
#pragma omp atomic
Ls[51] += L[51];
#pragma omp atomic
Ls[52] += L[52];
#pragma omp atomic
Ls[53] += L[53];
#pragma omp atomic
Ls[54] += L[54];
#pragma omp atomic
Ls[55] += L[55];

}

void field_m0_L2P_5(double x, double y, double z, double * L, double * F) {
double Ftmp0 = x*y;
double Ftmp1 = x*z;
double Ftmp2 = y*z;
double Ftmp3 = Ftmp0*z;
double Ftmp4 = (x*x);
double Ftmp5 = (1.0/2.0)*Ftmp4;
double Ftmp6 = (1.0/6.0)*(x*x*x);
double Ftmp7 = (1.0/24.0)*(x*x*x*x);
double Ftmp8 = (y*y);
double Ftmp9 = (1.0/2.0)*Ftmp8;
double Ftmp10 = (1.0/6.0)*(y*y*y);
double Ftmp11 = (1.0/24.0)*(y*y*y*y);
double Ftmp12 = (z*z);
double Ftmp13 = (1.0/2.0)*Ftmp12;
double Ftmp14 = (1.0/6.0)*(z*z*z);
double Ftmp15 = (1.0/24.0)*(z*z*z*z);
double Ftmp16 = Ftmp9*x;
double Ftmp17 = Ftmp10*x;
double Ftmp18 = Ftmp13*x;
double Ftmp19 = Ftmp14*x;
double Ftmp20 = Ftmp5*y;
double Ftmp21 = Ftmp5*z;
double Ftmp22 = Ftmp6*y;
double Ftmp23 = Ftmp6*z;
double Ftmp24 = Ftmp13*y;
double Ftmp25 = Ftmp14*y;
double Ftmp26 = Ftmp9*z;
double Ftmp27 = Ftmp10*z;
double Ftmp28 = Ftmp0*Ftmp13;
double Ftmp29 = Ftmp1*Ftmp9;
double Ftmp30 = Ftmp2*Ftmp5;
double Ftmp31 = (1.0/4.0)*Ftmp4;
double Ftmp32 = Ftmp31*Ftmp8;
double Ftmp33 = Ftmp12*Ftmp31;
double Ftmp34 = (1.0/4.0)*Ftmp12*Ftmp8;
#pragma omp atomic
F[0] += -Ftmp0*L[11] - Ftmp1*L[12] - Ftmp10*L[26] - Ftmp11*L[45] - Ftmp13*L[15] - Ftmp14*L[29] - Ftmp15*L[49] - Ftmp16*L[23] - Ftmp17*L[41] - Ftmp18*L[25] - Ftmp19*L[44] - Ftmp2*L[14] - Ftmp20*L[21] - Ftmp21*L[22] - Ftmp22*L[36] - Ftmp23*L[37] - Ftmp24*L[28] - Ftmp25*L[48] - Ftmp26*L[27] - Ftmp27*L[46] - Ftmp28*L[43] - Ftmp29*L[42] - Ftmp3*L[24] - Ftmp30*L[39] - Ftmp32*L[38] - Ftmp33*L[40] - Ftmp34*L[47] - Ftmp5*L[10] - Ftmp6*L[20] - Ftmp7*L[35] - Ftmp9*L[13] - x*L[4] - y*L[5] - z*L[6] - L[1];
#pragma omp atomic
F[1] += -Ftmp0*L[13] - Ftmp1*L[14] - Ftmp10*L[30] - Ftmp11*L[50] - Ftmp13*L[18] - Ftmp14*L[33] - Ftmp15*L[54] - Ftmp16*L[26] - Ftmp17*L[45] - Ftmp18*L[28] - Ftmp19*L[48] - Ftmp2*L[17] - Ftmp20*L[23] - Ftmp21*L[24] - Ftmp22*L[38] - Ftmp23*L[39] - Ftmp24*L[32] - Ftmp25*L[53] - Ftmp26*L[31] - Ftmp27*L[51] - Ftmp28*L[47] - Ftmp29*L[46] - Ftmp3*L[27] - Ftmp30*L[42] - Ftmp32*L[41] - Ftmp33*L[43] - Ftmp34*L[52] - Ftmp5*L[11] - Ftmp6*L[21] - Ftmp7*L[36] - Ftmp9*L[16] - x*L[5] - y*L[7] - z*L[8] - L[2];
#pragma omp atomic
F[2] += -Ftmp0*L[14] - Ftmp1*L[15] - Ftmp10*L[31] - Ftmp11*L[51] - Ftmp13*L[19] - Ftmp14*L[34] - Ftmp15*L[55] - Ftmp16*L[27] - Ftmp17*L[46] - Ftmp18*L[29] - Ftmp19*L[49] - Ftmp2*L[18] - Ftmp20*L[24] - Ftmp21*L[25] - Ftmp22*L[39] - Ftmp23*L[40] - Ftmp24*L[33] - Ftmp25*L[54] - Ftmp26*L[32] - Ftmp27*L[52] - Ftmp28*L[48] - Ftmp29*L[47] - Ftmp3*L[28] - Ftmp30*L[43] - Ftmp32*L[42] - Ftmp33*L[44] - Ftmp34*L[53] - Ftmp5*L[12] - Ftmp6*L[22] - Ftmp7*L[37] - Ftmp9*L[17] - x*L[6] - y*L[8] - z*L[9] - L[3];

}

void field_m0_M2P_5(double x, double y, double z, double * M, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = 1.0*pow(R, -3.0);
double Ftmp1 = pow(R, -5.0);
double Ftmp2 = 3.0*Ftmp1;
double Ftmp3 = Ftmp2*M[5];
double Ftmp4 = Ftmp2*z;
double Ftmp5 = Ftmp0*M[0];
double Ftmp6 = y*z;
double Ftmp7 = pow(R, -7.0);
double Ftmp8 = 15.0*Ftmp7;
double Ftmp9 = Ftmp8*M[14];
double Ftmp10 = x*y;
double Ftmp11 = Ftmp2*M[2];
double Ftmp12 = Ftmp4*M[3];
double Ftmp13 = (x*x);
double Ftmp14 = Ftmp2*M[1];
double Ftmp15 = Ftmp6*x;
double Ftmp16 = Ftmp15*Ftmp8;
double Ftmp17 = Ftmp13*Ftmp8;
double Ftmp18 = pow(R, -9.0);
double Ftmp19 = 105.0*Ftmp18;
double Ftmp20 = Ftmp13*Ftmp19;
double Ftmp21 = -9.0*Ftmp1;
double Ftmp22 = Ftmp17 + Ftmp21;
double Ftmp23 = -Ftmp2;
double Ftmp24 = (y*y);
double Ftmp25 = Ftmp24*Ftmp8;
double Ftmp26 = Ftmp23 + Ftmp25;
double Ftmp27 = (z*z);
double Ftmp28 = Ftmp27*Ftmp8;
double Ftmp29 = Ftmp23 + Ftmp28;
double Ftmp30 = Ftmp26*M[7];
double Ftmp31 = Ftmp29*M[9];
double Ftmp32 = 45.0*Ftmp7;
double Ftmp33 = -Ftmp32;
double Ftmp34 = Ftmp20 + Ftmp33;
double Ftmp35 = Ftmp34*M[21];
double Ftmp36 = Ftmp19*Ftmp24;
double Ftmp37 = Ftmp33 + Ftmp36;
double Ftmp38 = Ftmp37*y;
double Ftmp39 = 1.0*M[26];
double Ftmp40 = 3.0*y;
double Ftmp41 = 35.0*Ftmp18;
double Ftmp42 = (Ftmp27*Ftmp41 - 5.0*Ftmp7)*M[28];
double Ftmp43 = Ftmp34*M[22];
double Ftmp44 = 1.0*z;
double Ftmp45 = -Ftmp8;
double Ftmp46 = Ftmp36 + Ftmp45;
double Ftmp47 = Ftmp46*M[27];
double Ftmp48 = Ftmp19*Ftmp27;
double Ftmp49 = Ftmp33 + Ftmp48;
double Ftmp50 = Ftmp44*Ftmp49;
double Ftmp51 = 315.0*Ftmp18;
double Ftmp52 = -Ftmp51;
double Ftmp53 = pow(R, -11.0);
double Ftmp54 = 945.0*Ftmp53;
double Ftmp55 = Ftmp13*Ftmp54;
double Ftmp56 = Ftmp52 + Ftmp55;
double Ftmp57 = Ftmp56*M[39];
double Ftmp58 = Ftmp10*Ftmp34;
double Ftmp59 = Ftmp38*M[16];
double Ftmp60 = Ftmp45 + Ftmp48;
double Ftmp61 = Ftmp60*M[18];
double Ftmp62 = 1.0*Ftmp10;
double Ftmp63 = x*z;
double Ftmp64 = Ftmp34*Ftmp63;
double Ftmp65 = Ftmp46*M[17];
double Ftmp66 = Ftmp49*M[19];
double Ftmp67 = Ftmp24*Ftmp54;
double Ftmp68 = Ftmp52 + Ftmp67;
double Ftmp69 = Ftmp68*y;
double Ftmp70 = Ftmp44*M[46];
double Ftmp71 = -Ftmp19;
double Ftmp72 = Ftmp27*Ftmp53;
double Ftmp73 = 315.0*Ftmp72;
double Ftmp74 = Ftmp71 + Ftmp73;
double Ftmp75 = Ftmp74*M[48];
double Ftmp76 = Ftmp40*Ftmp75;
double Ftmp77 = Ftmp56*x;
double Ftmp78 = Ftmp68*M[31];
double Ftmp79 = -75.0*Ftmp7;
double Ftmp80 = 1.0*Ftmp13;
double Ftmp81 = Ftmp46*M[13];
double Ftmp82 = Ftmp60*M[15];
double Ftmp83 = 525.0*Ftmp18;
double Ftmp84 = -Ftmp83;
double Ftmp85 = Ftmp13*(Ftmp55 + Ftmp84);
double Ftmp86 = Ftmp27*Ftmp54;
double Ftmp87 = Ftmp52 + Ftmp86;
double Ftmp88 = Ftmp44*Ftmp87*M[33];
double Ftmp89 = (-Ftmp41 + Ftmp73)*M[28];
double Ftmp90 = Ftmp13*Ftmp40;
double Ftmp91 = Ftmp13*Ftmp44;
double Ftmp92 = (Ftmp67 + Ftmp71)*M[27];
double Ftmp93 = Ftmp87*M[29];
double Ftmp94 = 4725.0*Ftmp53;
double Ftmp95 = -Ftmp94;
double Ftmp96 = pow(R, -13.0);
double Ftmp97 = 10395.0*Ftmp96;
double Ftmp98 = Ftmp13*Ftmp97;
double Ftmp99 = 2835.0*Ftmp53;
double Ftmp100 = -Ftmp99;
double Ftmp101 = Ftmp24*Ftmp97;
double Ftmp102 = Ftmp100 + Ftmp101;
double Ftmp103 = 3465.0*Ftmp27*Ftmp96;
double Ftmp104 = (Ftmp103 - Ftmp54)*M[48];
double Ftmp105 = 225.0*Ftmp7;
double Ftmp106 = (x*x*x*x);
double Ftmp107 = Ftmp106*Ftmp54;
double Ftmp108 = 1050.0*Ftmp18;
double Ftmp109 = Ftmp105 + Ftmp107 - Ftmp108*Ftmp13;
double Ftmp110 = (y*y*y*y);
double Ftmp111 = Ftmp110*Ftmp54;
double Ftmp112 = 630.0*Ftmp18;
double Ftmp113 = Ftmp111 - Ftmp112*Ftmp24 + Ftmp32;
double Ftmp114 = (z*z*z*z);
double Ftmp115 = Ftmp114*Ftmp54;
double Ftmp116 = -Ftmp112*Ftmp27 + Ftmp115 + Ftmp32;
double Ftmp117 = Ftmp113*M[30];
double Ftmp118 = Ftmp116*M[34];
double Ftmp119 = 1575.0*Ftmp18;
double Ftmp120 = Ftmp106*Ftmp97;
double Ftmp121 = Ftmp13*Ftmp53;
double Ftmp122 = Ftmp119 + Ftmp120 - 9450.0*Ftmp121;
double Ftmp123 = Ftmp10*Ftmp122;
double Ftmp124 = Ftmp110*Ftmp97;
double Ftmp125 = 9450.0*Ftmp53;
double Ftmp126 = Ftmp119 + Ftmp124 - Ftmp125*Ftmp24;
double Ftmp127 = Ftmp126*M[50];
double Ftmp128 = Ftmp114*Ftmp97;
double Ftmp129 = 5670.0*Ftmp53;
double Ftmp130 = Ftmp128 - Ftmp129*Ftmp27 + Ftmp51;
double Ftmp131 = Ftmp130*M[54];
double Ftmp132 = Ftmp122*Ftmp63;
double Ftmp133 = Ftmp124 - Ftmp129*Ftmp24 + Ftmp51;
double Ftmp134 = Ftmp133*M[51];
double Ftmp135 = Ftmp119 - Ftmp125*Ftmp27 + Ftmp128;
double Ftmp136 = Ftmp135*M[55];
double Ftmp137 = 3675.0*Ftmp18;
double Ftmp138 = Ftmp133*M[45];
double Ftmp139 = Ftmp130*M[49];
double Ftmp140 = -Ftmp24*Ftmp51;
double Ftmp141 = Ftmp24*Ftmp55;
double Ftmp142 = -Ftmp20;
double Ftmp143 = Ftmp142 + Ftmp32;
double Ftmp144 = Ftmp140 + Ftmp141 + Ftmp143;
double Ftmp145 = -Ftmp27*Ftmp51;
double Ftmp146 = Ftmp27*Ftmp55;
double Ftmp147 = Ftmp143 + Ftmp145 + Ftmp146;
double Ftmp148 = -Ftmp48;
double Ftmp149 = Ftmp148 + Ftmp8;
double Ftmp150 = -Ftmp36;
double Ftmp151 = Ftmp27*Ftmp67;
double Ftmp152 = Ftmp150 + Ftmp151;
double Ftmp153 = Ftmp149 + Ftmp152;
double Ftmp154 = Ftmp24*Ftmp98;
double Ftmp155 = -Ftmp24*Ftmp99;
double Ftmp156 = Ftmp154 + Ftmp155;
double Ftmp157 = 945.0*Ftmp18;
double Ftmp158 = -Ftmp13*Ftmp99;
double Ftmp159 = Ftmp157 + Ftmp158;
double Ftmp160 = Ftmp10*(Ftmp156 + Ftmp159);
double Ftmp161 = -Ftmp27*Ftmp99;
double Ftmp162 = Ftmp161 + Ftmp51;
double Ftmp163 = -Ftmp55;
double Ftmp164 = Ftmp27*Ftmp98;
double Ftmp165 = Ftmp163 + Ftmp164;
double Ftmp166 = Ftmp10*(Ftmp162 + Ftmp165);
double Ftmp167 = -Ftmp67;
double Ftmp168 = Ftmp101*Ftmp27;
double Ftmp169 = Ftmp167 + Ftmp168;
double Ftmp170 = Ftmp10*(Ftmp162 + Ftmp169);
double Ftmp171 = Ftmp63*(Ftmp156 + Ftmp163 + Ftmp51);
double Ftmp172 = Ftmp63*(Ftmp159 + Ftmp161 + Ftmp164);
double Ftmp173 = -Ftmp86;
double Ftmp174 = Ftmp173 + Ftmp51;
double Ftmp175 = Ftmp155 + Ftmp168;
double Ftmp176 = Ftmp174 + Ftmp175;
double Ftmp177 = -Ftmp24*Ftmp94;
double Ftmp178 = Ftmp163 + Ftmp83;
double Ftmp179 = -Ftmp27*Ftmp94;
double Ftmp180 = Ftmp173 + Ftmp19;
double Ftmp181 = x*M[6];
double Ftmp182 = Ftmp17 + Ftmp23;
double Ftmp183 = Ftmp21 + Ftmp25;
double Ftmp184 = Ftmp182*M[4];
double Ftmp185 = Ftmp39*x;
double Ftmp186 = 3.0*x;
double Ftmp187 = Ftmp20 + Ftmp45;
double Ftmp188 = Ftmp187*M[24];
double Ftmp189 = Ftmp37*M[31];
double Ftmp190 = 1.0*x;
double Ftmp191 = Ftmp70*x;
double Ftmp192 = 3.0*Ftmp63;
double Ftmp193 = Ftmp187*M[12];
double Ftmp194 = Ftmp77*M[22];
double Ftmp195 = Ftmp187*M[11];
double Ftmp196 = 1.0*Ftmp24;
double Ftmp197 = Ftmp77*M[21];
double Ftmp198 = Ftmp24*z;
double Ftmp199 = (Ftmp55 + Ftmp71)*M[24];
double Ftmp200 = Ftmp67 + Ftmp84;
double Ftmp201 = Ftmp10*Ftmp44;
double Ftmp202 = (Ftmp100 + Ftmp98)*M[39];
double Ftmp203 = Ftmp107 - Ftmp112*Ftmp13 + Ftmp32;
double Ftmp204 = Ftmp105 - Ftmp108*Ftmp24 + Ftmp111;
double Ftmp205 = Ftmp203*M[20];
double Ftmp206 = Ftmp120 - 5670.0*Ftmp121 + Ftmp51;
double Ftmp207 = Ftmp206*M[37];
double Ftmp208 = Ftmp206*M[36];
double Ftmp209 = Ftmp141 + Ftmp150;
double Ftmp210 = -Ftmp13*Ftmp51 + Ftmp32;
double Ftmp211 = Ftmp209 + Ftmp210;
double Ftmp212 = Ftmp142 + Ftmp146 + Ftmp149;
double Ftmp213 = Ftmp145 + Ftmp152 + Ftmp32;
double Ftmp214 = Ftmp154 + Ftmp167;
double Ftmp215 = Ftmp6*(Ftmp158 + Ftmp214 + Ftmp51);
double Ftmp216 = Ftmp6*(Ftmp158 + Ftmp164 + Ftmp174);
double Ftmp217 = Ftmp6*(Ftmp157 + Ftmp161 + Ftmp175);
double Ftmp218 = -Ftmp13*Ftmp94 + Ftmp83;
double Ftmp219 = y*M[8];
double Ftmp220 = Ftmp21 + Ftmp28;
double Ftmp221 = Ftmp190*M[29];
double Ftmp222 = 1.0*y*M[33];
double Ftmp223 = Ftmp62*M[46];
double Ftmp224 = Ftmp44*x;
double Ftmp225 = Ftmp27*y;
double Ftmp226 = Ftmp27*(Ftmp84 + Ftmp86);
double Ftmp227 = Ftmp105 - Ftmp108*Ftmp27 + Ftmp115;
double Ftmp228 = Ftmp142 + Ftmp209 + Ftmp8;
double Ftmp229 = Ftmp146 + Ftmp148 + Ftmp210;
double Ftmp230 = Ftmp140 + Ftmp148 + Ftmp151 + Ftmp32;
#pragma omp atomic
F[0] += Ftmp0*M[1] - Ftmp10*Ftmp11 - Ftmp10*Ftmp127 + Ftmp10*Ftmp88 - Ftmp102*Ftmp13*Ftmp70*y - Ftmp104*Ftmp90*z + Ftmp109*x*M[20] + Ftmp109*M[35] + Ftmp113*M[45] + Ftmp116*M[49] + Ftmp117*x + Ftmp118*x - Ftmp12*x - Ftmp123*M[36] - Ftmp13*Ftmp14 + Ftmp13*Ftmp39*Ftmp69 - Ftmp13*Ftmp6*(Ftmp95 + Ftmp98)*M[39] - Ftmp13*(Ftmp20 + Ftmp79)*M[10] - Ftmp13*(Ftmp120 - 13230.0*Ftmp121 + Ftmp137)*M[35] - Ftmp13*(Ftmp154 + Ftmp177 + Ftmp178)*M[38] - Ftmp13*(Ftmp164 + Ftmp178 + Ftmp179)*M[40] - Ftmp131*Ftmp62 - Ftmp132*M[37] - Ftmp134*Ftmp63 - Ftmp136*Ftmp63 - Ftmp138*Ftmp80 - Ftmp139*Ftmp80 + Ftmp144*x*M[23] + Ftmp144*M[38] + Ftmp147*x*M[25] + Ftmp147*M[40] + Ftmp15*Ftmp78 + Ftmp153*x*M[32] + Ftmp153*M[47] + Ftmp16*M[8] - Ftmp160*M[41] - Ftmp166*M[43] + Ftmp17*y*M[5] + Ftmp17*z*M[6] - Ftmp170*M[52] - Ftmp171*M[42] - Ftmp172*M[44] - Ftmp176*Ftmp63*M[53] - Ftmp20*Ftmp6*M[14] + Ftmp22*x*M[4] + Ftmp22*M[10] + Ftmp26*M[13] + Ftmp29*M[15] - Ftmp3*y + Ftmp30*x + Ftmp31*x - Ftmp35*y - Ftmp38*Ftmp39 - Ftmp4*M[6] - Ftmp40*Ftmp42 - Ftmp43*z - Ftmp44*Ftmp47 + Ftmp5*x - Ftmp50*M[29] + Ftmp57*Ftmp6 - Ftmp58*M[11] - Ftmp59*x + Ftmp6*Ftmp77*M[24] + Ftmp6*Ftmp9 - Ftmp61*Ftmp62 - Ftmp63*Ftmp65 - Ftmp63*Ftmp66 - Ftmp64*M[12] + Ftmp69*Ftmp70 + Ftmp76*z - Ftmp80*Ftmp81 - Ftmp80*Ftmp82 - Ftmp80*(Ftmp169 + Ftmp180)*M[47] + Ftmp85*y*M[21] + Ftmp85*z*M[22] + Ftmp89*Ftmp90 + Ftmp91*Ftmp92 + Ftmp91*Ftmp93;
#pragma omp atomic
F[1] += Ftmp0*M[2] - Ftmp10*Ftmp14 - Ftmp104*Ftmp192*Ftmp24 - Ftmp11*Ftmp24 + Ftmp116*M[54] + Ftmp118*y - Ftmp12*y - Ftmp123*M[35] - Ftmp126*Ftmp6*M[51] - Ftmp126*Ftmp62*M[45] - Ftmp131*Ftmp196 - Ftmp136*Ftmp6 - Ftmp139*Ftmp62 - Ftmp160*M[38] - Ftmp166*M[40] - 1.0*Ftmp170*M[47] + Ftmp181*Ftmp6*Ftmp8 + Ftmp182*M[11] + Ftmp183*y*M[7] + Ftmp183*M[16] + Ftmp184*y + Ftmp185*Ftmp200*Ftmp24 - Ftmp185*Ftmp37 + Ftmp186*Ftmp24*Ftmp89 - Ftmp186*Ftmp42 - Ftmp188*z - Ftmp189*z - Ftmp190*Ftmp38*M[13] - Ftmp191*Ftmp24*(Ftmp101 + Ftmp95) + Ftmp191*Ftmp68 + Ftmp192*Ftmp75 - Ftmp193*Ftmp6 + Ftmp194*Ftmp6 - Ftmp195*Ftmp24 - Ftmp196*Ftmp61 + Ftmp197*Ftmp24 + Ftmp198*Ftmp199 + Ftmp198*Ftmp200*M[31] + Ftmp201*Ftmp68*M[27] + Ftmp201*Ftmp93 - Ftmp202*Ftmp24*Ftmp63 + Ftmp203*M[36] + Ftmp204*y*M[30] + Ftmp204*M[50] + Ftmp205*y - Ftmp207*Ftmp6 - Ftmp208*Ftmp24 + Ftmp211*y*M[23] + Ftmp211*M[41] + Ftmp212*y*M[25] + Ftmp212*M[43] + Ftmp213*y*M[32] + Ftmp213*M[52] - Ftmp215*M[42] - Ftmp216*M[44] - Ftmp217*M[53] + Ftmp24*Ftmp88 - Ftmp24*(Ftmp165 + Ftmp180)*M[43] - Ftmp24*(Ftmp214 + Ftmp218)*M[41] - Ftmp24*(Ftmp36 + Ftmp79)*M[16] - Ftmp24*(Ftmp124 + Ftmp137 - 13230.0*Ftmp24*Ftmp53)*M[50] - Ftmp24*(Ftmp169 + Ftmp179 + Ftmp83)*M[52] + Ftmp25*x*M[5] + Ftmp25*z*M[8] + Ftmp29*M[18] - Ftmp3*x + Ftmp31*y - Ftmp35*x - Ftmp36*Ftmp63*M[14] - Ftmp38*z*M[17] - Ftmp4*M[8] + Ftmp5*y - Ftmp50*M[33] + Ftmp57*Ftmp63 - Ftmp58*M[10] - Ftmp6*Ftmp66 - Ftmp62*Ftmp82 + Ftmp63*Ftmp9;
#pragma omp atomic
F[2] += Ftmp0*M[3] - Ftmp10*Ftmp202*Ftmp27 - Ftmp10*Ftmp48*M[14] + Ftmp10*Ftmp57 + Ftmp10*Ftmp9 - Ftmp102*Ftmp223*Ftmp27 + Ftmp113*M[51] + Ftmp117*z - Ftmp127*Ftmp6 - Ftmp132*M[35] - Ftmp134*Ftmp27 - Ftmp135*Ftmp224*M[49] - Ftmp135*Ftmp44*y*M[54] - Ftmp138*Ftmp224 + Ftmp16*M[5] - Ftmp171*M[38] - Ftmp172*M[40] - Ftmp176*Ftmp224*M[47] - Ftmp181*Ftmp2 + Ftmp181*Ftmp28 + Ftmp182*M[12] + Ftmp184*z + Ftmp185*Ftmp6*Ftmp68 - Ftmp188*y - Ftmp189*y + Ftmp190*Ftmp27*Ftmp92 - Ftmp190*Ftmp47 - Ftmp193*Ftmp27 + Ftmp194*Ftmp27 - Ftmp195*Ftmp6 + Ftmp197*Ftmp6 + Ftmp199*Ftmp225 - Ftmp2*Ftmp219 - Ftmp2*Ftmp27*M[3] + Ftmp203*M[37] + Ftmp205*z - Ftmp207*Ftmp27 - Ftmp208*Ftmp6 - Ftmp215*M[41] - Ftmp216*M[43] - Ftmp217*M[52] + Ftmp219*Ftmp28 + Ftmp220*z*M[9] + Ftmp220*M[19] + Ftmp221*Ftmp226 - Ftmp221*Ftmp49 + Ftmp222*Ftmp226 - Ftmp222*Ftmp49 + Ftmp223*Ftmp68 - Ftmp224*Ftmp81 + Ftmp225*Ftmp78 + Ftmp227*z*M[34] + Ftmp227*M[55] + Ftmp228*z*M[23] + Ftmp228*M[42] + Ftmp229*z*M[25] + Ftmp229*M[44] + Ftmp230*z*M[32] + Ftmp230*M[53] + Ftmp26*M[17] - Ftmp27*Ftmp40*x*(Ftmp103 - 1575.0*Ftmp53)*M[48] - Ftmp27*Ftmp65 - Ftmp27*(Ftmp48 + Ftmp79)*M[19] - Ftmp27*(Ftmp128 + Ftmp137 - 13230.0*Ftmp72)*M[55] - Ftmp27*(Ftmp163 + Ftmp19 + Ftmp214)*M[42] - Ftmp27*(Ftmp164 + Ftmp173 + Ftmp218)*M[44] - Ftmp27*(Ftmp168 + Ftmp173 + Ftmp177 + Ftmp83)*M[53] + Ftmp30*z - Ftmp4*x*M[1] - Ftmp4*y*M[2] + Ftmp40*Ftmp63*Ftmp74*M[28] - Ftmp43*x + Ftmp5*z - Ftmp50*x*M[15] - Ftmp50*y*M[18] - Ftmp59*z - Ftmp64*M[10] + Ftmp76*x;

}

void field_m0_P2M_6(double x, double y, double z, double q, double * M) {
double Mtmp0 = q*x;
double Mtmp1 = q*y;
double Mtmp2 = q*z;
double Mtmp3 = (x*x);
double Mtmp4 = (1.0/2.0)*q;
double Mtmp5 = Mtmp0*y;
double Mtmp6 = Mtmp0*z;
double Mtmp7 = (y*y);
double Mtmp8 = Mtmp1*z;
double Mtmp9 = (z*z);
double Mtmp10 = (x*x*x);
double Mtmp11 = (1.0/6.0)*q;
double Mtmp12 = (1.0/2.0)*Mtmp3;
double Mtmp13 = (1.0/2.0)*Mtmp0;
double Mtmp14 = (y*y*y);
double Mtmp15 = (1.0/2.0)*Mtmp7;
double Mtmp16 = (1.0/2.0)*Mtmp9;
double Mtmp17 = (z*z*z);
double Mtmp18 = (x*x*x*x);
double Mtmp19 = (1.0/24.0)*q;
double Mtmp20 = (1.0/6.0)*Mtmp10;
double Mtmp21 = Mtmp7*q;
double Mtmp22 = (1.0/4.0)*Mtmp3;
double Mtmp23 = Mtmp9*q;
double Mtmp24 = (1.0/6.0)*Mtmp0;
double Mtmp25 = (y*y*y*y);
double Mtmp26 = (1.0/6.0)*Mtmp14;
double Mtmp27 = (1.0/4.0)*Mtmp9;
double Mtmp28 = (1.0/6.0)*Mtmp17;
double Mtmp29 = (z*z*z*z);
double Mtmp30 = (x*x*x*x*x);
double Mtmp31 = (1.0/120.0)*q;
double Mtmp32 = (1.0/24.0)*Mtmp18;
double Mtmp33 = (1.0/12.0)*Mtmp10;
double Mtmp34 = (1.0/12.0)*Mtmp14;
double Mtmp35 = Mtmp3*q;
double Mtmp36 = Mtmp2*Mtmp7;
double Mtmp37 = Mtmp1*Mtmp9;
double Mtmp38 = (1.0/12.0)*Mtmp17;
double Mtmp39 = (1.0/24.0)*Mtmp0;
double Mtmp40 = Mtmp0*Mtmp7;
double Mtmp41 = (y*y*y*y*y);
double Mtmp42 = (1.0/24.0)*Mtmp25;
double Mtmp43 = (1.0/24.0)*Mtmp29;
double Mtmp44 = (z*z*z*z*z);
double Mtmp45 = (1.0/720.0)*q;
double Mtmp46 = (1.0/120.0)*Mtmp30;
double Mtmp47 = (1.0/48.0)*Mtmp18;
double Mtmp48 = (1.0/36.0)*Mtmp10*q;
double Mtmp49 = (1.0/48.0)*Mtmp35;
double Mtmp50 = (1.0/120.0)*Mtmp0;
M[0] += q;
M[1] += -Mtmp0;
M[2] += -Mtmp1;
M[3] += -Mtmp2;
M[4] += Mtmp3*Mtmp4;
M[5] += Mtmp5;
M[6] += Mtmp6;
M[7] += Mtmp4*Mtmp7;
M[8] += Mtmp8;
M[9] += Mtmp4*Mtmp9;
M[10] += -Mtmp10*Mtmp11;
M[11] += -Mtmp1*Mtmp12;
M[12] += -Mtmp12*Mtmp2;
M[13] += -Mtmp13*Mtmp7;
M[14] += -Mtmp5*z;
M[15] += -Mtmp13*Mtmp9;
M[16] += -Mtmp11*Mtmp14;
M[17] += -Mtmp15*Mtmp2;
M[18] += -Mtmp1*Mtmp16;
M[19] += -Mtmp11*Mtmp17;
M[20] += Mtmp18*Mtmp19;
M[21] += Mtmp1*Mtmp20;
M[22] += Mtmp2*Mtmp20;
M[23] += Mtmp21*Mtmp22;
M[24] += Mtmp12*Mtmp8;
M[25] += Mtmp22*Mtmp23;
M[26] += Mtmp14*Mtmp24;
M[27] += Mtmp15*Mtmp6;
M[28] += Mtmp16*Mtmp5;
M[29] += Mtmp17*Mtmp24;
M[30] += Mtmp19*Mtmp25;
M[31] += Mtmp2*Mtmp26;
M[32] += Mtmp21*Mtmp27;
M[33] += Mtmp1*Mtmp28;
M[34] += Mtmp19*Mtmp29;
M[35] += -Mtmp30*Mtmp31;
M[36] += -Mtmp1*Mtmp32;
M[37] += -Mtmp2*Mtmp32;
M[38] += -Mtmp21*Mtmp33;
M[39] += -Mtmp20*Mtmp8;
M[40] += -Mtmp23*Mtmp33;
M[41] += -Mtmp34*Mtmp35;
M[42] += -Mtmp22*Mtmp36;
M[43] += -Mtmp22*Mtmp37;
M[44] += -Mtmp35*Mtmp38;
M[45] += -Mtmp25*Mtmp39;
M[46] += -Mtmp26*Mtmp6;
M[47] += -Mtmp27*Mtmp40;
M[48] += -Mtmp28*Mtmp5;
M[49] += -Mtmp29*Mtmp39;
M[50] += -Mtmp31*Mtmp41;
M[51] += -Mtmp2*Mtmp42;
M[52] += -Mtmp23*Mtmp34;
M[53] += -Mtmp21*Mtmp38;
M[54] += -Mtmp1*Mtmp43;
M[55] += -Mtmp31*Mtmp44;
M[56] += Mtmp45*(x*x*x*x*x*x);
M[57] += Mtmp1*Mtmp46;
M[58] += Mtmp2*Mtmp46;
M[59] += Mtmp21*Mtmp47;
M[60] += Mtmp32*Mtmp8;
M[61] += Mtmp23*Mtmp47;
M[62] += Mtmp14*Mtmp48;
M[63] += Mtmp33*Mtmp36;
M[64] += Mtmp33*Mtmp37;
M[65] += Mtmp17*Mtmp48;
M[66] += Mtmp25*Mtmp49;
M[67] += Mtmp2*Mtmp3*Mtmp34;
M[68] += (1.0/8.0)*Mtmp21*Mtmp3*Mtmp9;
M[69] += Mtmp1*Mtmp3*Mtmp38;
M[70] += Mtmp29*Mtmp49;
M[71] += Mtmp41*Mtmp50;
M[72] += Mtmp42*Mtmp6;
M[73] += Mtmp0*Mtmp34*Mtmp9;
M[74] += Mtmp38*Mtmp40;
M[75] += Mtmp43*Mtmp5;
M[76] += Mtmp44*Mtmp50;
M[77] += Mtmp45*(y*y*y*y*y*y);
M[78] += (1.0/120.0)*Mtmp2*Mtmp41;
M[79] += (1.0/48.0)*Mtmp23*Mtmp25;
M[80] += (1.0/36.0)*Mtmp14*Mtmp17*q;
M[81] += (1.0/48.0)*Mtmp21*Mtmp29;
M[82] += (1.0/120.0)*Mtmp1*Mtmp44;
M[83] += Mtmp45*(z*z*z*z*z*z);

}
void field_m0_M2M_6(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = y*M[0];
double Mstmp2 = z*M[0];
double Mstmp3 = x*M[1];
double Mstmp4 = (x*x);
double Mstmp5 = (1.0/2.0)*M[0];
double Mstmp6 = x*M[2];
double Mstmp7 = y*M[1];
double Mstmp8 = Mstmp0*y;
double Mstmp9 = x*M[3];
double Mstmp10 = z*M[1];
double Mstmp11 = Mstmp0*z;
double Mstmp12 = y*M[2];
double Mstmp13 = (y*y);
double Mstmp14 = y*M[3];
double Mstmp15 = z*M[2];
double Mstmp16 = Mstmp1*z;
double Mstmp17 = z*M[3];
double Mstmp18 = (z*z);
double Mstmp19 = x*M[4];
double Mstmp20 = (1.0/2.0)*Mstmp4;
double Mstmp21 = (x*x*x);
double Mstmp22 = (1.0/6.0)*M[0];
double Mstmp23 = x*M[5];
double Mstmp24 = y*M[4];
double Mstmp25 = Mstmp3*y;
double Mstmp26 = x*M[6];
double Mstmp27 = z*M[4];
double Mstmp28 = Mstmp3*z;
double Mstmp29 = x*M[7];
double Mstmp30 = y*M[5];
double Mstmp31 = Mstmp6*y;
double Mstmp32 = (1.0/2.0)*M[1];
double Mstmp33 = (1.0/2.0)*Mstmp13;
double Mstmp34 = x*M[8];
double Mstmp35 = y*M[6];
double Mstmp36 = z*M[5];
double Mstmp37 = Mstmp9*y;
double Mstmp38 = Mstmp6*z;
double Mstmp39 = Mstmp7*z;
double Mstmp40 = x*M[9];
double Mstmp41 = z*M[6];
double Mstmp42 = Mstmp9*z;
double Mstmp43 = (1.0/2.0)*Mstmp18;
double Mstmp44 = y*M[7];
double Mstmp45 = (y*y*y);
double Mstmp46 = y*M[8];
double Mstmp47 = z*M[7];
double Mstmp48 = Mstmp12*z;
double Mstmp49 = y*M[9];
double Mstmp50 = z*M[8];
double Mstmp51 = Mstmp14*z;
double Mstmp52 = z*M[9];
double Mstmp53 = (z*z*z);
double Mstmp54 = x*M[10];
double Mstmp55 = (1.0/6.0)*Mstmp21;
double Mstmp56 = (x*x*x*x);
double Mstmp57 = (1.0/24.0)*M[0];
double Mstmp58 = x*M[11];
double Mstmp59 = y*M[10];
double Mstmp60 = Mstmp19*y;
double Mstmp61 = x*M[12];
double Mstmp62 = z*M[10];
double Mstmp63 = Mstmp19*z;
double Mstmp64 = x*M[13];
double Mstmp65 = y*M[11];
double Mstmp66 = Mstmp23*y;
double Mstmp67 = Mstmp13*M[0];
double Mstmp68 = (1.0/4.0)*Mstmp4;
double Mstmp69 = x*M[14];
double Mstmp70 = y*M[12];
double Mstmp71 = z*M[11];
double Mstmp72 = Mstmp26*y;
double Mstmp73 = Mstmp23*z;
double Mstmp74 = Mstmp24*z;
double Mstmp75 = x*M[15];
double Mstmp76 = z*M[12];
double Mstmp77 = Mstmp26*z;
double Mstmp78 = Mstmp18*Mstmp68;
double Mstmp79 = x*M[16];
double Mstmp80 = y*M[13];
double Mstmp81 = Mstmp29*y;
double Mstmp82 = (1.0/6.0)*M[1];
double Mstmp83 = (1.0/6.0)*Mstmp45;
double Mstmp84 = x*M[17];
double Mstmp85 = y*M[14];
double Mstmp86 = z*M[13];
double Mstmp87 = Mstmp34*y;
double Mstmp88 = Mstmp29*z;
double Mstmp89 = Mstmp30*z;
double Mstmp90 = x*M[18];
double Mstmp91 = y*M[15];
double Mstmp92 = z*M[14];
double Mstmp93 = Mstmp40*y;
double Mstmp94 = Mstmp34*z;
double Mstmp95 = Mstmp35*z;
double Mstmp96 = x*M[19];
double Mstmp97 = z*M[15];
double Mstmp98 = Mstmp40*z;
double Mstmp99 = (1.0/6.0)*Mstmp53;
double Mstmp100 = y*M[16];
double Mstmp101 = (y*y*y*y);
double Mstmp102 = y*M[17];
double Mstmp103 = z*M[16];
double Mstmp104 = Mstmp44*z;
double Mstmp105 = y*M[18];
double Mstmp106 = z*M[17];
double Mstmp107 = Mstmp46*z;
double Mstmp108 = (1.0/4.0)*Mstmp18;
double Mstmp109 = y*M[19];
double Mstmp110 = z*M[18];
double Mstmp111 = Mstmp49*z;
double Mstmp112 = z*M[19];
double Mstmp113 = (z*z*z*z);
double Mstmp114 = x*M[20];
double Mstmp115 = (1.0/24.0)*Mstmp56;
double Mstmp116 = (x*x*x*x*x);
double Mstmp117 = (1.0/120.0)*M[0];
double Mstmp118 = x*M[21];
double Mstmp119 = y*M[20];
double Mstmp120 = Mstmp54*y;
double Mstmp121 = x*M[22];
double Mstmp122 = x*M[23];
double Mstmp123 = y*M[21];
double Mstmp124 = Mstmp58*y;
double Mstmp125 = Mstmp13*Mstmp68;
double Mstmp126 = (1.0/12.0)*Mstmp21;
double Mstmp127 = x*M[24];
double Mstmp128 = y*M[22];
double Mstmp129 = Mstmp61*y;
double Mstmp130 = x*M[25];
double Mstmp131 = Mstmp18*M[0];
double Mstmp132 = x*M[26];
double Mstmp133 = y*M[23];
double Mstmp134 = Mstmp64*y;
double Mstmp135 = (1.0/12.0)*Mstmp45;
double Mstmp136 = Mstmp4*M[0];
double Mstmp137 = x*M[27];
double Mstmp138 = y*M[24];
double Mstmp139 = Mstmp69*y;
double Mstmp140 = x*M[28];
double Mstmp141 = y*M[25];
double Mstmp142 = Mstmp75*y;
double Mstmp143 = x*M[29];
double Mstmp144 = (1.0/12.0)*Mstmp53;
double Mstmp145 = x*M[30];
double Mstmp146 = y*M[26];
double Mstmp147 = Mstmp79*y;
double Mstmp148 = (1.0/24.0)*M[1];
double Mstmp149 = (1.0/24.0)*Mstmp101;
double Mstmp150 = x*M[31];
double Mstmp151 = y*M[27];
double Mstmp152 = Mstmp84*y;
double Mstmp153 = x*M[32];
double Mstmp154 = y*M[28];
double Mstmp155 = Mstmp90*y;
double Mstmp156 = Mstmp108*Mstmp13;
double Mstmp157 = x*M[33];
double Mstmp158 = y*M[29];
double Mstmp159 = Mstmp96*y;
double Mstmp160 = x*M[34];
double Mstmp161 = (1.0/24.0)*Mstmp113;
double Mstmp162 = y*M[30];
double Mstmp163 = (y*y*y*y*y);
double Mstmp164 = y*M[31];
double Mstmp165 = y*M[32];
double Mstmp166 = y*M[33];
double Mstmp167 = y*M[34];
double Mstmp168 = (z*z*z*z*z);
double Mstmp169 = (1.0/120.0)*Mstmp116;
double Mstmp170 = (1.0/720.0)*M[0];
double Mstmp171 = Mstmp126*M[1];
double Mstmp172 = (1.0/48.0)*Mstmp56;
double Mstmp173 = Mstmp4*M[1];
double Mstmp174 = Mstmp126*Mstmp13;
double Mstmp175 = (1.0/36.0)*Mstmp21*M[0];
double Mstmp176 = Mstmp126*Mstmp18;
double Mstmp177 = Mstmp135*Mstmp4;
double Mstmp178 = (1.0/48.0)*Mstmp136;
double Mstmp179 = Mstmp144*Mstmp4;
double Mstmp180 = (1.0/120.0)*M[1];
double Mstmp181 = (1.0/120.0)*Mstmp163;
double Mstmp182 = Mstmp135*Mstmp18;
double Mstmp183 = Mstmp13*Mstmp144;
double Mstmp184 = (1.0/120.0)*Mstmp168;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += Mstmp0 + M[1];
#pragma omp atomic
Ms[2] += Mstmp1 + M[2];
#pragma omp atomic
Ms[3] += Mstmp2 + M[3];
#pragma omp atomic
Ms[4] += Mstmp3 + Mstmp4*Mstmp5 + M[4];
#pragma omp atomic
Ms[5] += Mstmp6 + Mstmp7 + Mstmp8 + M[5];
#pragma omp atomic
Ms[6] += Mstmp10 + Mstmp11 + Mstmp9 + M[6];
#pragma omp atomic
Ms[7] += Mstmp12 + Mstmp13*Mstmp5 + M[7];
#pragma omp atomic
Ms[8] += Mstmp14 + Mstmp15 + Mstmp16 + M[8];
#pragma omp atomic
Ms[9] += Mstmp17 + Mstmp18*Mstmp5 + M[9];
#pragma omp atomic
Ms[10] += Mstmp19 + Mstmp20*M[1] + Mstmp21*Mstmp22 + M[10];
#pragma omp atomic
Ms[11] += Mstmp1*Mstmp20 + Mstmp20*M[2] + Mstmp23 + Mstmp24 + Mstmp25 + M[11];
#pragma omp atomic
Ms[12] += Mstmp2*Mstmp20 + Mstmp20*M[3] + Mstmp26 + Mstmp27 + Mstmp28 + M[12];
#pragma omp atomic
Ms[13] += Mstmp0*Mstmp33 + Mstmp13*Mstmp32 + Mstmp29 + Mstmp30 + Mstmp31 + M[13];
#pragma omp atomic
Ms[14] += Mstmp34 + Mstmp35 + Mstmp36 + Mstmp37 + Mstmp38 + Mstmp39 + Mstmp8*z + M[14];
#pragma omp atomic
Ms[15] += Mstmp0*Mstmp43 + Mstmp18*Mstmp32 + Mstmp40 + Mstmp41 + Mstmp42 + M[15];
#pragma omp atomic
Ms[16] += Mstmp22*Mstmp45 + Mstmp33*M[2] + Mstmp44 + M[16];
#pragma omp atomic
Ms[17] += Mstmp2*Mstmp33 + Mstmp33*M[3] + Mstmp46 + Mstmp47 + Mstmp48 + M[17];
#pragma omp atomic
Ms[18] += Mstmp1*Mstmp43 + Mstmp43*M[2] + Mstmp49 + Mstmp50 + Mstmp51 + M[18];
#pragma omp atomic
Ms[19] += Mstmp22*Mstmp53 + Mstmp43*M[3] + Mstmp52 + M[19];
#pragma omp atomic
Ms[20] += Mstmp20*M[4] + Mstmp54 + Mstmp55*M[1] + Mstmp56*Mstmp57 + M[20];
#pragma omp atomic
Ms[21] += Mstmp1*Mstmp55 + Mstmp20*Mstmp7 + Mstmp20*M[5] + Mstmp55*M[2] + Mstmp58 + Mstmp59 + Mstmp60 + M[21];
#pragma omp atomic
Ms[22] += Mstmp10*Mstmp20 + Mstmp2*Mstmp55 + Mstmp20*M[6] + Mstmp55*M[3] + Mstmp61 + Mstmp62 + Mstmp63 + M[22];
#pragma omp atomic
Ms[23] += Mstmp12*Mstmp20 + Mstmp20*M[7] + Mstmp3*Mstmp33 + Mstmp33*M[4] + Mstmp64 + Mstmp65 + Mstmp66 + Mstmp67*Mstmp68 + M[23];
#pragma omp atomic
Ms[24] += Mstmp14*Mstmp20 + Mstmp15*Mstmp20 + Mstmp16*Mstmp20 + Mstmp20*M[8] + Mstmp25*z + Mstmp69 + Mstmp70 + Mstmp71 + Mstmp72 + Mstmp73 + Mstmp74 + M[24];
#pragma omp atomic
Ms[25] += Mstmp17*Mstmp20 + Mstmp20*M[9] + Mstmp3*Mstmp43 + Mstmp43*M[4] + Mstmp75 + Mstmp76 + Mstmp77 + Mstmp78*M[0] + M[25];
#pragma omp atomic
Ms[26] += Mstmp0*Mstmp83 + Mstmp33*Mstmp6 + Mstmp33*M[5] + Mstmp45*Mstmp82 + Mstmp79 + Mstmp80 + Mstmp81 + M[26];
#pragma omp atomic
Ms[27] += Mstmp10*Mstmp33 + Mstmp11*Mstmp33 + Mstmp31*z + Mstmp33*Mstmp9 + Mstmp33*M[6] + Mstmp84 + Mstmp85 + Mstmp86 + Mstmp87 + Mstmp88 + Mstmp89 + M[27];
#pragma omp atomic
Ms[28] += Mstmp37*z + Mstmp43*Mstmp6 + Mstmp43*Mstmp7 + Mstmp43*Mstmp8 + Mstmp43*M[5] + Mstmp90 + Mstmp91 + Mstmp92 + Mstmp93 + Mstmp94 + Mstmp95 + M[28];
#pragma omp atomic
Ms[29] += Mstmp0*Mstmp99 + Mstmp43*Mstmp9 + Mstmp43*M[6] + Mstmp53*Mstmp82 + Mstmp96 + Mstmp97 + Mstmp98 + M[29];
#pragma omp atomic
Ms[30] += Mstmp100 + Mstmp101*Mstmp57 + Mstmp33*M[7] + Mstmp83*M[2] + M[30];
#pragma omp atomic
Ms[31] += Mstmp102 + Mstmp103 + Mstmp104 + Mstmp15*Mstmp33 + Mstmp2*Mstmp83 + Mstmp33*M[8] + Mstmp83*M[3] + M[31];
#pragma omp atomic
Ms[32] += Mstmp105 + Mstmp106 + Mstmp107 + Mstmp108*Mstmp67 + Mstmp12*Mstmp43 + Mstmp17*Mstmp33 + Mstmp33*M[9] + Mstmp43*M[7] + M[32];
#pragma omp atomic
Ms[33] += Mstmp1*Mstmp99 + Mstmp109 + Mstmp110 + Mstmp111 + Mstmp14*Mstmp43 + Mstmp43*M[8] + Mstmp99*M[2] + M[33];
#pragma omp atomic
Ms[34] += Mstmp112 + Mstmp113*Mstmp57 + Mstmp43*M[9] + Mstmp99*M[3] + M[34];
#pragma omp atomic
Ms[35] += Mstmp114 + Mstmp115*M[1] + Mstmp116*Mstmp117 + Mstmp20*M[10] + Mstmp55*M[4] + M[35];
#pragma omp atomic
Ms[36] += Mstmp1*Mstmp115 + Mstmp115*M[2] + Mstmp118 + Mstmp119 + Mstmp120 + Mstmp20*Mstmp24 + Mstmp20*M[11] + Mstmp55*Mstmp7 + Mstmp55*M[5] + M[36];
#pragma omp atomic
Ms[37] += Mstmp10*Mstmp55 + Mstmp115*Mstmp2 + Mstmp115*M[3] + Mstmp121 + Mstmp20*Mstmp27 + Mstmp20*M[12] + Mstmp54*z + Mstmp55*M[6] + z*M[20] + M[37];
#pragma omp atomic
Ms[38] += Mstmp12*Mstmp55 + Mstmp122 + Mstmp123 + Mstmp124 + Mstmp125*M[1] + Mstmp126*Mstmp67 + Mstmp19*Mstmp33 + Mstmp20*Mstmp30 + Mstmp20*M[13] + Mstmp33*M[10] + Mstmp55*M[7] + M[38];
#pragma omp atomic
Ms[39] += Mstmp127 + Mstmp128 + Mstmp129 + Mstmp14*Mstmp55 + Mstmp15*Mstmp55 + Mstmp16*Mstmp55 + Mstmp20*Mstmp35 + Mstmp20*Mstmp36 + Mstmp20*Mstmp39 + Mstmp20*M[14] + Mstmp55*M[8] + Mstmp58*z + Mstmp59*z + Mstmp60*z + z*M[21] + M[39];
#pragma omp atomic
Ms[40] += Mstmp126*Mstmp131 + Mstmp130 + Mstmp17*Mstmp55 + Mstmp19*Mstmp43 + Mstmp20*Mstmp41 + Mstmp20*M[15] + Mstmp43*M[10] + Mstmp55*M[9] + Mstmp61*z + Mstmp78*M[1] + z*M[22] + M[40];
#pragma omp atomic
Ms[41] += Mstmp125*M[2] + Mstmp132 + Mstmp133 + Mstmp134 + Mstmp135*Mstmp136 + Mstmp20*Mstmp44 + Mstmp20*M[16] + Mstmp23*Mstmp33 + Mstmp3*Mstmp83 + Mstmp33*M[11] + Mstmp83*M[4] + M[41];
#pragma omp atomic
Ms[42] += Mstmp125*Mstmp2 + Mstmp125*M[3] + Mstmp137 + Mstmp138 + Mstmp139 + Mstmp20*Mstmp46 + Mstmp20*Mstmp47 + Mstmp20*Mstmp48 + Mstmp20*M[17] + Mstmp26*Mstmp33 + Mstmp27*Mstmp33 + Mstmp28*Mstmp33 + Mstmp33*M[12] + Mstmp64*z + Mstmp65*z + Mstmp66*z + z*M[23] + M[42];
#pragma omp atomic
Ms[43] += Mstmp1*Mstmp78 + Mstmp140 + Mstmp141 + Mstmp142 + Mstmp20*Mstmp49 + Mstmp20*Mstmp50 + Mstmp20*Mstmp51 + Mstmp20*M[18] + Mstmp23*Mstmp43 + Mstmp24*Mstmp43 + Mstmp25*Mstmp43 + Mstmp43*M[11] + Mstmp69*z + Mstmp70*z + Mstmp72*z + Mstmp78*M[2] + z*M[24] + M[43];
#pragma omp atomic
Ms[44] += Mstmp136*Mstmp144 + Mstmp143 + Mstmp20*Mstmp52 + Mstmp20*M[19] + Mstmp26*Mstmp43 + Mstmp3*Mstmp99 + Mstmp43*M[12] + Mstmp75*z + Mstmp78*M[3] + Mstmp99*M[4] + z*M[25] + M[44];
#pragma omp atomic
Ms[45] += Mstmp0*Mstmp149 + Mstmp101*Mstmp148 + Mstmp145 + Mstmp146 + Mstmp147 + Mstmp29*Mstmp33 + Mstmp33*M[13] + Mstmp6*Mstmp83 + Mstmp83*M[5] + M[45];
#pragma omp atomic
Ms[46] += Mstmp10*Mstmp83 + Mstmp11*Mstmp83 + Mstmp150 + Mstmp151 + Mstmp152 + Mstmp33*Mstmp34 + Mstmp33*Mstmp36 + Mstmp33*Mstmp38 + Mstmp33*M[14] + Mstmp79*z + Mstmp80*z + Mstmp81*z + Mstmp83*Mstmp9 + Mstmp83*M[6] + z*M[26] + M[46];
#pragma omp atomic
Ms[47] += Mstmp0*Mstmp156 + Mstmp153 + Mstmp154 + Mstmp155 + Mstmp156*M[1] + Mstmp29*Mstmp43 + Mstmp30*Mstmp43 + Mstmp31*Mstmp43 + Mstmp33*Mstmp40 + Mstmp33*Mstmp41 + Mstmp33*Mstmp42 + Mstmp33*M[15] + Mstmp43*M[13] + Mstmp84*z + Mstmp85*z + Mstmp87*z + z*M[27] + M[47];
#pragma omp atomic
Ms[48] += Mstmp157 + Mstmp158 + Mstmp159 + Mstmp34*Mstmp43 + Mstmp35*Mstmp43 + Mstmp37*Mstmp43 + Mstmp43*M[14] + Mstmp6*Mstmp99 + Mstmp7*Mstmp99 + Mstmp8*Mstmp99 + Mstmp90*z + Mstmp91*z + Mstmp93*z + Mstmp99*M[5] + z*M[28] + M[48];
#pragma omp atomic
Ms[49] += Mstmp0*Mstmp161 + Mstmp113*Mstmp148 + Mstmp160 + Mstmp40*Mstmp43 + Mstmp43*M[15] + Mstmp9*Mstmp99 + Mstmp96*z + Mstmp99*M[6] + z*M[29] + M[49];
#pragma omp atomic
Ms[50] += Mstmp117*Mstmp163 + Mstmp149*M[2] + Mstmp162 + Mstmp33*M[16] + Mstmp83*M[7] + M[50];
#pragma omp atomic
Ms[51] += Mstmp100*z + Mstmp149*Mstmp2 + Mstmp149*M[3] + Mstmp15*Mstmp83 + Mstmp164 + Mstmp33*Mstmp47 + Mstmp33*M[17] + Mstmp83*M[8] + z*M[30] + M[51];
#pragma omp atomic
Ms[52] += Mstmp102*z + Mstmp131*Mstmp135 + Mstmp156*M[2] + Mstmp165 + Mstmp17*Mstmp83 + Mstmp33*Mstmp50 + Mstmp33*M[18] + Mstmp43*Mstmp44 + Mstmp43*M[16] + Mstmp83*M[9] + z*M[31] + M[52];
#pragma omp atomic
Ms[53] += Mstmp105*z + Mstmp12*Mstmp99 + Mstmp144*Mstmp67 + Mstmp156*M[3] + Mstmp166 + Mstmp33*Mstmp52 + Mstmp33*M[19] + Mstmp43*Mstmp46 + Mstmp43*M[17] + Mstmp99*M[7] + z*M[32] + M[53];
#pragma omp atomic
Ms[54] += Mstmp1*Mstmp161 + Mstmp109*z + Mstmp14*Mstmp99 + Mstmp161*M[2] + Mstmp167 + Mstmp43*Mstmp49 + Mstmp43*M[18] + Mstmp99*M[8] + z*M[33] + M[54];
#pragma omp atomic
Ms[55] += Mstmp117*Mstmp168 + Mstmp161*M[3] + Mstmp43*M[19] + Mstmp99*M[9] + z*M[34] + M[55];
#pragma omp atomic
Ms[56] += Mstmp115*M[4] + Mstmp169*M[1] + Mstmp170*(x*x*x*x*x*x) + Mstmp20*M[20] + Mstmp55*M[10] + x*M[35] + M[56];
#pragma omp atomic
Ms[57] += Mstmp1*Mstmp169 + Mstmp114*y + Mstmp115*Mstmp7 + Mstmp115*M[5] + Mstmp169*M[2] + Mstmp20*Mstmp59 + Mstmp20*M[21] + Mstmp24*Mstmp55 + Mstmp55*M[11] + x*M[36] + y*M[35] + M[57];
#pragma omp atomic
Ms[58] += Mstmp10*Mstmp115 + Mstmp114*z + Mstmp115*M[6] + Mstmp169*Mstmp2 + Mstmp169*M[3] + Mstmp20*Mstmp62 + Mstmp20*M[22] + Mstmp27*Mstmp55 + Mstmp55*M[12] + x*M[37] + z*M[35] + M[58];
#pragma omp atomic
Ms[59] += Mstmp115*Mstmp12 + Mstmp115*M[7] + Mstmp118*y + Mstmp125*M[4] + Mstmp13*Mstmp171 + Mstmp172*Mstmp67 + Mstmp20*Mstmp65 + Mstmp20*M[23] + Mstmp30*Mstmp55 + Mstmp33*Mstmp54 + Mstmp33*M[20] + Mstmp55*M[13] + x*M[38] + y*M[36] + M[59];
#pragma omp atomic
Ms[60] += Mstmp115*Mstmp14 + Mstmp115*Mstmp15 + Mstmp115*Mstmp16 + Mstmp115*M[8] + Mstmp118*z + Mstmp119*z + Mstmp120*z + Mstmp121*y + Mstmp20*Mstmp70 + Mstmp20*Mstmp71 + Mstmp20*Mstmp74 + Mstmp20*M[24] + Mstmp35*Mstmp55 + Mstmp36*Mstmp55 + Mstmp39*Mstmp55 + Mstmp55*M[14] + x*M[39] + y*M[37] + z*M[36] + M[60];
#pragma omp atomic
Ms[61] += Mstmp115*Mstmp17 + Mstmp115*M[9] + Mstmp121*z + Mstmp131*Mstmp172 + Mstmp171*Mstmp18 + Mstmp20*Mstmp76 + Mstmp20*M[25] + Mstmp41*Mstmp55 + Mstmp43*Mstmp54 + Mstmp43*M[20] + Mstmp55*M[15] + Mstmp78*M[4] + x*M[40] + z*M[37] + M[61];
#pragma omp atomic
Ms[62] += Mstmp122*y + Mstmp125*M[5] + Mstmp135*Mstmp173 + Mstmp174*M[2] + Mstmp175*Mstmp45 + Mstmp19*Mstmp83 + Mstmp20*Mstmp80 + Mstmp20*M[26] + Mstmp33*Mstmp58 + Mstmp33*M[21] + Mstmp44*Mstmp55 + Mstmp55*M[16] + Mstmp83*M[10] + x*M[41] + y*M[38] + M[62];
#pragma omp atomic
Ms[63] += Mstmp10*Mstmp125 + Mstmp122*z + Mstmp123*z + Mstmp124*z + Mstmp125*M[6] + Mstmp127*y + Mstmp174*Mstmp2 + Mstmp174*M[3] + Mstmp20*Mstmp85 + Mstmp20*Mstmp86 + Mstmp20*Mstmp89 + Mstmp20*M[27] + Mstmp33*Mstmp61 + Mstmp33*Mstmp62 + Mstmp33*Mstmp63 + Mstmp33*M[22] + Mstmp46*Mstmp55 + Mstmp47*Mstmp55 + Mstmp48*Mstmp55 + Mstmp55*M[17] + x*M[42] + y*M[39] + z*M[38] + M[63];
#pragma omp atomic
Ms[64] += Mstmp1*Mstmp176 + Mstmp127*z + Mstmp128*z + Mstmp129*z + Mstmp130*y + Mstmp176*M[2] + Mstmp20*Mstmp91 + Mstmp20*Mstmp92 + Mstmp20*Mstmp95 + Mstmp20*M[28] + Mstmp43*Mstmp58 + Mstmp43*Mstmp59 + Mstmp43*Mstmp60 + Mstmp43*M[21] + Mstmp49*Mstmp55 + Mstmp50*Mstmp55 + Mstmp51*Mstmp55 + Mstmp55*M[18] + Mstmp7*Mstmp78 + Mstmp78*M[5] + x*M[43] + y*M[40] + z*M[39] + M[64];
#pragma omp atomic
Ms[65] += Mstmp130*z + Mstmp144*Mstmp173 + Mstmp175*Mstmp53 + Mstmp176*M[3] + Mstmp19*Mstmp99 + Mstmp20*Mstmp97 + Mstmp20*M[29] + Mstmp43*Mstmp61 + Mstmp43*M[22] + Mstmp52*Mstmp55 + Mstmp55*M[19] + Mstmp78*M[6] + Mstmp99*M[10] + x*M[44] + z*M[40] + M[65];
#pragma omp atomic
Ms[66] += Mstmp100*Mstmp20 + Mstmp101*Mstmp178 + Mstmp125*M[7] + Mstmp132*y + Mstmp149*Mstmp3 + Mstmp149*M[4] + Mstmp177*M[2] + Mstmp20*M[30] + Mstmp23*Mstmp83 + Mstmp33*Mstmp64 + Mstmp33*M[23] + Mstmp83*M[11] + x*M[45] + y*M[41] + M[66];
#pragma omp atomic
Ms[67] += Mstmp102*Mstmp20 + Mstmp103*Mstmp20 + Mstmp104*Mstmp20 + Mstmp125*Mstmp15 + Mstmp125*M[8] + Mstmp132*z + Mstmp133*z + Mstmp134*z + Mstmp137*y + Mstmp177*Mstmp2 + Mstmp177*M[3] + Mstmp20*M[31] + Mstmp26*Mstmp83 + Mstmp27*Mstmp83 + Mstmp28*Mstmp83 + Mstmp33*Mstmp69 + Mstmp33*Mstmp71 + Mstmp33*Mstmp73 + Mstmp33*M[24] + Mstmp83*M[12] + x*M[46] + y*M[42] + z*M[41] + M[67];
#pragma omp atomic
Ms[68] += Mstmp105*Mstmp20 + Mstmp106*Mstmp20 + Mstmp107*Mstmp20 + Mstmp12*Mstmp78 + Mstmp125*Mstmp17 + Mstmp125*M[9] + Mstmp137*z + Mstmp138*z + Mstmp139*z + Mstmp140*y + Mstmp156*Mstmp3 + Mstmp156*M[4] + (1.0/8.0)*Mstmp18*Mstmp4*Mstmp67 + Mstmp20*M[32] + Mstmp33*Mstmp75 + Mstmp33*Mstmp76 + Mstmp33*Mstmp77 + Mstmp33*M[25] + Mstmp43*Mstmp64 + Mstmp43*Mstmp65 + Mstmp43*Mstmp66 + Mstmp43*M[23] + Mstmp78*M[7] + x*M[47] + y*M[43] + z*M[42] + M[68];
#pragma omp atomic
Ms[69] += Mstmp1*Mstmp179 + Mstmp109*Mstmp20 + Mstmp110*Mstmp20 + Mstmp111*Mstmp20 + Mstmp14*Mstmp78 + Mstmp140*z + Mstmp141*z + Mstmp142*z + Mstmp143*y + Mstmp179*M[2] + Mstmp20*M[33] + Mstmp23*Mstmp99 + Mstmp24*Mstmp99 + Mstmp25*Mstmp99 + Mstmp43*Mstmp69 + Mstmp43*Mstmp70 + Mstmp43*Mstmp72 + Mstmp43*M[24] + Mstmp78*M[8] + Mstmp99*M[11] + x*M[48] + y*M[44] + z*M[43] + M[69];
#pragma omp atomic
Ms[70] += Mstmp112*Mstmp20 + Mstmp113*Mstmp178 + Mstmp143*z + Mstmp161*Mstmp3 + Mstmp161*M[4] + Mstmp179*M[3] + Mstmp20*M[34] + Mstmp26*Mstmp99 + Mstmp43*Mstmp75 + Mstmp43*M[25] + Mstmp78*M[9] + Mstmp99*M[12] + x*M[49] + z*M[44] + M[70];
#pragma omp atomic
Ms[71] += Mstmp0*Mstmp181 + Mstmp145*y + Mstmp149*Mstmp6 + Mstmp149*M[5] + Mstmp163*Mstmp180 + Mstmp29*Mstmp83 + Mstmp33*Mstmp79 + Mstmp33*M[26] + Mstmp83*M[13] + x*M[50] + y*M[45] + M[71];
#pragma omp atomic
Ms[72] += Mstmp10*Mstmp149 + Mstmp11*Mstmp149 + Mstmp145*z + Mstmp146*z + Mstmp147*z + Mstmp149*Mstmp9 + Mstmp149*M[6] + Mstmp150*y + Mstmp33*Mstmp84 + Mstmp33*Mstmp86 + Mstmp33*Mstmp88 + Mstmp33*M[27] + Mstmp34*Mstmp83 + Mstmp36*Mstmp83 + Mstmp38*Mstmp83 + Mstmp83*M[14] + x*M[51] + y*M[46] + z*M[45] + M[72];
#pragma omp atomic
Ms[73] += Mstmp0*Mstmp182 + Mstmp150*z + Mstmp151*z + Mstmp152*z + Mstmp153*y + Mstmp156*Mstmp6 + Mstmp156*M[5] + Mstmp182*M[1] + Mstmp33*Mstmp90 + Mstmp33*Mstmp92 + Mstmp33*Mstmp94 + Mstmp33*M[28] + Mstmp40*Mstmp83 + Mstmp41*Mstmp83 + Mstmp42*Mstmp83 + Mstmp43*Mstmp79 + Mstmp43*Mstmp80 + Mstmp43*Mstmp81 + Mstmp43*M[26] + Mstmp83*M[15] + x*M[52] + y*M[47] + z*M[46] + M[73];
#pragma omp atomic
Ms[74] += Mstmp0*Mstmp183 + Mstmp153*z + Mstmp154*z + Mstmp155*z + Mstmp156*Mstmp9 + Mstmp156*M[6] + Mstmp157*y + Mstmp183*M[1] + Mstmp29*Mstmp99 + Mstmp30*Mstmp99 + Mstmp31*Mstmp99 + Mstmp33*Mstmp96 + Mstmp33*Mstmp97 + Mstmp33*Mstmp98 + Mstmp33*M[29] + Mstmp43*Mstmp84 + Mstmp43*Mstmp85 + Mstmp43*Mstmp87 + Mstmp43*M[27] + Mstmp99*M[13] + x*M[53] + y*M[48] + z*M[47] + M[74];
#pragma omp atomic
Ms[75] += Mstmp157*z + Mstmp158*z + Mstmp159*z + Mstmp160*y + Mstmp161*Mstmp6 + Mstmp161*Mstmp7 + Mstmp161*Mstmp8 + Mstmp161*M[5] + Mstmp34*Mstmp99 + Mstmp35*Mstmp99 + Mstmp37*Mstmp99 + Mstmp43*Mstmp90 + Mstmp43*Mstmp91 + Mstmp43*Mstmp93 + Mstmp43*M[28] + Mstmp99*M[14] + x*M[54] + y*M[49] + z*M[48] + M[75];
#pragma omp atomic
Ms[76] += Mstmp0*Mstmp184 + Mstmp160*z + Mstmp161*Mstmp9 + Mstmp161*M[6] + Mstmp168*Mstmp180 + Mstmp40*Mstmp99 + Mstmp43*Mstmp96 + Mstmp43*M[29] + Mstmp99*M[15] + x*M[55] + z*M[49] + M[76];
#pragma omp atomic
Ms[77] += Mstmp149*M[7] + Mstmp170*(y*y*y*y*y*y) + Mstmp181*M[2] + Mstmp33*M[30] + Mstmp83*M[16] + y*M[50] + M[77];
#pragma omp atomic
Ms[78] += Mstmp103*Mstmp33 + Mstmp149*Mstmp15 + Mstmp149*M[8] + Mstmp162*z + Mstmp181*Mstmp2 + Mstmp181*M[3] + Mstmp33*M[31] + Mstmp47*Mstmp83 + Mstmp83*M[17] + y*M[51] + z*M[50] + M[78];
#pragma omp atomic
Ms[79] += Mstmp100*Mstmp43 + (1.0/48.0)*Mstmp101*Mstmp131 + Mstmp106*Mstmp33 + Mstmp149*Mstmp17 + Mstmp149*M[9] + Mstmp156*M[7] + Mstmp164*z + Mstmp182*M[2] + Mstmp33*M[32] + Mstmp43*M[30] + Mstmp50*Mstmp83 + Mstmp83*M[18] + y*M[52] + z*M[51] + M[79];
#pragma omp atomic
Ms[80] += Mstmp102*Mstmp43 + Mstmp110*Mstmp33 + Mstmp156*M[8] + Mstmp165*z + Mstmp182*M[3] + Mstmp183*M[2] + Mstmp33*M[33] + Mstmp43*M[31] + Mstmp44*Mstmp99 + (1.0/36.0)*Mstmp45*Mstmp53*M[0] + Mstmp52*Mstmp83 + Mstmp83*M[19] + Mstmp99*M[16] + y*M[53] + z*M[52] + M[80];
#pragma omp atomic
Ms[81] += Mstmp105*Mstmp43 + Mstmp112*Mstmp33 + (1.0/48.0)*Mstmp113*Mstmp67 + Mstmp12*Mstmp161 + Mstmp156*M[9] + Mstmp161*M[7] + Mstmp166*z + Mstmp183*M[3] + Mstmp33*M[34] + Mstmp43*M[32] + Mstmp46*Mstmp99 + Mstmp99*M[17] + y*M[54] + z*M[53] + M[81];
#pragma omp atomic
Ms[82] += Mstmp1*Mstmp184 + Mstmp109*Mstmp43 + Mstmp14*Mstmp161 + Mstmp161*M[8] + Mstmp167*z + Mstmp184*M[2] + Mstmp43*M[33] + Mstmp49*Mstmp99 + Mstmp99*M[18] + y*M[55] + z*M[54] + M[82];
#pragma omp atomic
Ms[83] += Mstmp161*M[9] + Mstmp170*(z*z*z*z*z*z) + Mstmp184*M[3] + Mstmp43*M[34] + Mstmp99*M[19] + z*M[55] + M[83];

}

void field_m0_M2L_6(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double D[84];
double Dtmp0 = 1.0*pow(R, -3.0);
double Dtmp1 = -Dtmp0;
double Dtmp2 = (x*x);
double Dtmp3 = pow(R, -5.0);
double Dtmp4 = 3.0*Dtmp3;
double Dtmp5 = x*y;
double Dtmp6 = x*z;
double Dtmp7 = (y*y);
double Dtmp8 = y*z;
double Dtmp9 = 9.0*Dtmp3;
double Dtmp10 = -Dtmp9;
double Dtmp11 = pow(R, -7.0);
double Dtmp12 = 15.0*Dtmp11;
double Dtmp13 = Dtmp12*Dtmp2;
double Dtmp14 = -Dtmp4;
double Dtmp15 = Dtmp13 + Dtmp14;
double Dtmp16 = Dtmp12*Dtmp7;
double Dtmp17 = Dtmp14 + Dtmp16;
double Dtmp18 = 1.0*x;
double Dtmp19 = Dtmp8*x;
double Dtmp20 = (x*x*x*x);
double Dtmp21 = pow(R, -9.0);
double Dtmp22 = 105.0*Dtmp21;
double Dtmp23 = 90.0*Dtmp11;
double Dtmp24 = 45.0*Dtmp11;
double Dtmp25 = -Dtmp24;
double Dtmp26 = Dtmp2*Dtmp22;
double Dtmp27 = x*(Dtmp25 + Dtmp26);
double Dtmp28 = -Dtmp12;
double Dtmp29 = Dtmp22*Dtmp7;
double Dtmp30 = Dtmp25 + Dtmp29;
double Dtmp31 = Dtmp18*y;
double Dtmp32 = Dtmp18*z;
double Dtmp33 = (y*y*y*y);
double Dtmp34 = 225.0*Dtmp11;
double Dtmp35 = pow(R, -11.0);
double Dtmp36 = 945.0*Dtmp35;
double Dtmp37 = Dtmp20*Dtmp36;
double Dtmp38 = Dtmp2*Dtmp21;
double Dtmp39 = 630.0*Dtmp38;
double Dtmp40 = Dtmp24 + Dtmp37 - Dtmp39;
double Dtmp41 = -Dtmp26;
double Dtmp42 = 315.0*Dtmp21;
double Dtmp43 = Dtmp42*Dtmp7;
double Dtmp44 = Dtmp2*Dtmp36;
double Dtmp45 = Dtmp44*Dtmp7;
double Dtmp46 = Dtmp24 + Dtmp45;
double Dtmp47 = -Dtmp42;
double Dtmp48 = -Dtmp29;
double Dtmp49 = Dtmp2*Dtmp42;
double Dtmp50 = Dtmp33*Dtmp36;
double Dtmp51 = Dtmp21*Dtmp7;
double Dtmp52 = 630.0*Dtmp51;
double Dtmp53 = Dtmp24 + Dtmp50 - Dtmp52;
double Dtmp54 = Dtmp36*Dtmp7;
double Dtmp55 = -Dtmp34;
double Dtmp56 = 10395.0*pow(R, -13.0);
double Dtmp57 = 14175.0*Dtmp35;
double Dtmp58 = 1575.0*Dtmp21;
double Dtmp59 = Dtmp20*Dtmp56;
double Dtmp60 = Dtmp2*Dtmp35;
double Dtmp61 = x*(Dtmp58 + Dtmp59 - 9450.0*Dtmp60);
double Dtmp62 = 5670.0*Dtmp60;
double Dtmp63 = Dtmp25 - Dtmp62*Dtmp7;
double Dtmp64 = -2835.0*Dtmp60;
double Dtmp65 = Dtmp35*Dtmp7;
double Dtmp66 = Dtmp2*Dtmp56*Dtmp7;
double Dtmp67 = -2835.0*Dtmp65 + Dtmp66;
double Dtmp68 = Dtmp33*Dtmp56;
double Dtmp69 = Dtmp58 - 9450.0*Dtmp65 + Dtmp68;
D[0] = 1.0/R;
D[1] = -Dtmp0*x;
D[2] = -Dtmp0*y;
D[3] = -Dtmp0*z;
D[4] = Dtmp1 + Dtmp2*Dtmp4;
D[5] = Dtmp4*Dtmp5;
D[6] = Dtmp4*Dtmp6;
D[7] = Dtmp1 + Dtmp4*Dtmp7;
D[8] = Dtmp4*Dtmp8;
D[9] = -D[4] - D[7];
D[10] = -x*(Dtmp10 + Dtmp13);
D[11] = -Dtmp15*y;
D[12] = -Dtmp15*z;
D[13] = -Dtmp17*Dtmp18;
D[14] = -Dtmp12*Dtmp19;
D[15] = -D[10] - D[13];
D[16] = -y*(Dtmp10 + Dtmp16);
D[17] = -Dtmp17*z;
D[18] = -D[11] - D[16];
D[19] = -D[12] - D[17];
D[20] = -Dtmp2*Dtmp23 + Dtmp20*Dtmp22 + Dtmp9;
D[21] = Dtmp27*y;
D[22] = Dtmp27*z;
D[23] = -Dtmp13 - Dtmp16 + Dtmp26*Dtmp7 + Dtmp4;
D[24] = Dtmp8*(Dtmp26 + Dtmp28);
D[25] = -D[20] - D[23];
D[26] = Dtmp30*Dtmp31;
D[27] = Dtmp32*(Dtmp28 + Dtmp29);
D[28] = -D[21] - D[26];
D[29] = -D[22] - D[27];
D[30] = Dtmp22*Dtmp33 - Dtmp23*Dtmp7 + Dtmp9;
D[31] = Dtmp30*Dtmp8;
D[32] = -D[23] - D[30];
D[33] = -D[24] - D[31];
D[34] = -D[25] - D[32];
D[35] = -x*(Dtmp34 + Dtmp37 - 1050.0*Dtmp38);
D[36] = -Dtmp40*y;
D[37] = -Dtmp40*z;
D[38] = -x*(Dtmp41 - Dtmp43 + Dtmp46);
D[39] = -Dtmp19*(Dtmp44 + Dtmp47);
D[40] = -D[35] - D[38];
D[41] = -y*(Dtmp46 + Dtmp48 - Dtmp49);
D[42] = -z*(Dtmp12 + Dtmp41 + Dtmp45 + Dtmp48);
D[43] = -D[36] - D[41];
D[44] = -D[37] - D[42];
D[45] = -Dtmp18*Dtmp53;
D[46] = -Dtmp18*Dtmp8*(Dtmp47 + Dtmp54);
D[47] = -D[38] - D[45];
D[48] = -D[39] - D[46];
D[49] = -D[40] - D[47];
D[50] = -y*(Dtmp34 + Dtmp50 - 1050.0*Dtmp51);
D[51] = -Dtmp53*z;
D[52] = -D[41] - D[50];
D[53] = -D[42] - D[51];
D[54] = -D[43] - D[52];
D[55] = -D[44] - D[53];
D[56] = -Dtmp20*Dtmp57 + 4725.0*Dtmp38 + Dtmp55 + Dtmp56*(x*x*x*x*x*x);
D[57] = Dtmp61*y;
D[58] = Dtmp61*z;
D[59] = -Dtmp37 + Dtmp39 + Dtmp43 + Dtmp59*Dtmp7 + Dtmp63;
D[60] = Dtmp8*(Dtmp42 + Dtmp59 - Dtmp62);
D[61] = -D[56] - D[59];
D[62] = Dtmp5*(945.0*Dtmp21 + Dtmp64 + Dtmp67);
D[63] = Dtmp6*(Dtmp42 - Dtmp44 + Dtmp67);
D[64] = -D[57] - D[62];
D[65] = -D[58] - D[63];
D[66] = Dtmp2*Dtmp68 + Dtmp49 - Dtmp50 + Dtmp52 + Dtmp63;
D[67] = Dtmp8*(Dtmp42 - Dtmp54 + Dtmp64 + Dtmp66);
D[68] = -D[59] - D[66];
D[69] = -D[60] - D[67];
D[70] = -D[61] - D[68];
D[71] = Dtmp31*Dtmp69;
D[72] = Dtmp32*(Dtmp42 - 5670.0*Dtmp65 + Dtmp68);
D[73] = -D[62] - D[71];
D[74] = -D[63] - D[72];
D[75] = -D[64] - D[73];
D[76] = -D[65] - D[74];
D[77] = -Dtmp33*Dtmp57 + 4725.0*Dtmp51 + Dtmp55 + Dtmp56*(y*y*y*y*y*y);
D[78] = Dtmp69*Dtmp8;
D[79] = -D[66] - D[77];
D[80] = -D[67] - D[78];
D[81] = -D[68] - D[79];
D[82] = -D[69] - D[80];
D[83] = -D[70] - D[81];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[8]*M[8] + D[9]*M[9] + D[10]*M[10] + D[11]*M[11] + D[12]*M[12] + D[13]*M[13] + D[14]*M[14] + D[15]*M[15] + D[16]*M[16] + D[17]*M[17] + D[18]*M[18] + D[19]*M[19] + D[20]*M[20] + D[21]*M[21] + D[22]*M[22] + D[23]*M[23] + D[24]*M[24] + D[25]*M[25] + D[26]*M[26] + D[27]*M[27] + D[28]*M[28] + D[29]*M[29] + D[30]*M[30] + D[31]*M[31] + D[32]*M[32] + D[33]*M[33] + D[34]*M[34] + D[35]*M[35] + D[36]*M[36] + D[37]*M[37] + D[38]*M[38] + D[39]*M[39] + D[40]*M[40] + D[41]*M[41] + D[42]*M[42] + D[43]*M[43] + D[44]*M[44] + D[45]*M[45] + D[46]*M[46] + D[47]*M[47] + D[48]*M[48] + D[49]*M[49] + D[50]*M[50] + D[51]*M[51] + D[52]*M[52] + D[53]*M[53] + D[54]*M[54] + D[55]*M[55] + D[56]*M[56] + D[57]*M[57] + D[58]*M[58] + D[59]*M[59] + D[60]*M[60] + D[61]*M[61] + D[62]*M[62] + D[63]*M[63] + D[64]*M[64] + D[65]*M[65] + D[66]*M[66] + D[67]*M[67] + D[68]*M[68] + D[69]*M[69] + D[70]*M[70] + D[71]*M[71] + D[72]*M[72] + D[73]*M[73] + D[74]*M[74] + D[75]*M[75] + D[76]*M[76] + D[77]*M[77] + D[78]*M[78] + D[79]*M[79] + D[80]*M[80] + D[81]*M[81] + D[82]*M[82] + D[83]*M[83];
#pragma omp atomic
L[1] += D[1]*M[0] + D[4]*M[1] + D[5]*M[2] + D[6]*M[3] + D[10]*M[4] + D[11]*M[5] + D[12]*M[6] + D[13]*M[7] + D[14]*M[8] + D[15]*M[9] + D[20]*M[10] + D[21]*M[11] + D[22]*M[12] + D[23]*M[13] + D[24]*M[14] + D[25]*M[15] + D[26]*M[16] + D[27]*M[17] + D[28]*M[18] + D[29]*M[19] + D[35]*M[20] + D[36]*M[21] + D[37]*M[22] + D[38]*M[23] + D[39]*M[24] + D[40]*M[25] + D[41]*M[26] + D[42]*M[27] + D[43]*M[28] + D[44]*M[29] + D[45]*M[30] + D[46]*M[31] + D[47]*M[32] + D[48]*M[33] + D[49]*M[34] + D[56]*M[35] + D[57]*M[36] + D[58]*M[37] + D[59]*M[38] + D[60]*M[39] + D[61]*M[40] + D[62]*M[41] + D[63]*M[42] + D[64]*M[43] + D[65]*M[44] + D[66]*M[45] + D[67]*M[46] + D[68]*M[47] + D[69]*M[48] + D[70]*M[49] + D[71]*M[50] + D[72]*M[51] + D[73]*M[52] + D[74]*M[53] + D[75]*M[54] + D[76]*M[55];
#pragma omp atomic
L[2] += D[2]*M[0] + D[5]*M[1] + D[7]*M[2] + D[8]*M[3] + D[11]*M[4] + D[13]*M[5] + D[14]*M[6] + D[16]*M[7] + D[17]*M[8] + D[18]*M[9] + D[21]*M[10] + D[23]*M[11] + D[24]*M[12] + D[26]*M[13] + D[27]*M[14] + D[28]*M[15] + D[30]*M[16] + D[31]*M[17] + D[32]*M[18] + D[33]*M[19] + D[36]*M[20] + D[38]*M[21] + D[39]*M[22] + D[41]*M[23] + D[42]*M[24] + D[43]*M[25] + D[45]*M[26] + D[46]*M[27] + D[47]*M[28] + D[48]*M[29] + D[50]*M[30] + D[51]*M[31] + D[52]*M[32] + D[53]*M[33] + D[54]*M[34] + D[57]*M[35] + D[59]*M[36] + D[60]*M[37] + D[62]*M[38] + D[63]*M[39] + D[64]*M[40] + D[66]*M[41] + D[67]*M[42] + D[68]*M[43] + D[69]*M[44] + D[71]*M[45] + D[72]*M[46] + D[73]*M[47] + D[74]*M[48] + D[75]*M[49] + D[77]*M[50] + D[78]*M[51] + D[79]*M[52] + D[80]*M[53] + D[81]*M[54] + D[82]*M[55];
#pragma omp atomic
L[3] += D[3]*M[0] + D[6]*M[1] + D[8]*M[2] + D[9]*M[3] + D[12]*M[4] + D[14]*M[5] + D[15]*M[6] + D[17]*M[7] + D[18]*M[8] + D[19]*M[9] + D[22]*M[10] + D[24]*M[11] + D[25]*M[12] + D[27]*M[13] + D[28]*M[14] + D[29]*M[15] + D[31]*M[16] + D[32]*M[17] + D[33]*M[18] + D[34]*M[19] + D[37]*M[20] + D[39]*M[21] + D[40]*M[22] + D[42]*M[23] + D[43]*M[24] + D[44]*M[25] + D[46]*M[26] + D[47]*M[27] + D[48]*M[28] + D[49]*M[29] + D[51]*M[30] + D[52]*M[31] + D[53]*M[32] + D[54]*M[33] + D[55]*M[34] + D[58]*M[35] + D[60]*M[36] + D[61]*M[37] + D[63]*M[38] + D[64]*M[39] + D[65]*M[40] + D[67]*M[41] + D[68]*M[42] + D[69]*M[43] + D[70]*M[44] + D[72]*M[45] + D[73]*M[46] + D[74]*M[47] + D[75]*M[48] + D[76]*M[49] + D[78]*M[50] + D[79]*M[51] + D[80]*M[52] + D[81]*M[53] + D[82]*M[54] + D[83]*M[55];
#pragma omp atomic
L[4] += D[4]*M[0] + D[10]*M[1] + D[11]*M[2] + D[12]*M[3] + D[20]*M[4] + D[21]*M[5] + D[22]*M[6] + D[23]*M[7] + D[24]*M[8] + D[25]*M[9] + D[35]*M[10] + D[36]*M[11] + D[37]*M[12] + D[38]*M[13] + D[39]*M[14] + D[40]*M[15] + D[41]*M[16] + D[42]*M[17] + D[43]*M[18] + D[44]*M[19] + D[56]*M[20] + D[57]*M[21] + D[58]*M[22] + D[59]*M[23] + D[60]*M[24] + D[61]*M[25] + D[62]*M[26] + D[63]*M[27] + D[64]*M[28] + D[65]*M[29] + D[66]*M[30] + D[67]*M[31] + D[68]*M[32] + D[69]*M[33] + D[70]*M[34];
#pragma omp atomic
L[5] += D[5]*M[0] + D[11]*M[1] + D[13]*M[2] + D[14]*M[3] + D[21]*M[4] + D[23]*M[5] + D[24]*M[6] + D[26]*M[7] + D[27]*M[8] + D[28]*M[9] + D[36]*M[10] + D[38]*M[11] + D[39]*M[12] + D[41]*M[13] + D[42]*M[14] + D[43]*M[15] + D[45]*M[16] + D[46]*M[17] + D[47]*M[18] + D[48]*M[19] + D[57]*M[20] + D[59]*M[21] + D[60]*M[22] + D[62]*M[23] + D[63]*M[24] + D[64]*M[25] + D[66]*M[26] + D[67]*M[27] + D[68]*M[28] + D[69]*M[29] + D[71]*M[30] + D[72]*M[31] + D[73]*M[32] + D[74]*M[33] + D[75]*M[34];
#pragma omp atomic
L[6] += D[6]*M[0] + D[12]*M[1] + D[14]*M[2] + D[15]*M[3] + D[22]*M[4] + D[24]*M[5] + D[25]*M[6] + D[27]*M[7] + D[28]*M[8] + D[29]*M[9] + D[37]*M[10] + D[39]*M[11] + D[40]*M[12] + D[42]*M[13] + D[43]*M[14] + D[44]*M[15] + D[46]*M[16] + D[47]*M[17] + D[48]*M[18] + D[49]*M[19] + D[58]*M[20] + D[60]*M[21] + D[61]*M[22] + D[63]*M[23] + D[64]*M[24] + D[65]*M[25] + D[67]*M[26] + D[68]*M[27] + D[69]*M[28] + D[70]*M[29] + D[72]*M[30] + D[73]*M[31] + D[74]*M[32] + D[75]*M[33] + D[76]*M[34];
#pragma omp atomic
L[7] += D[7]*M[0] + D[13]*M[1] + D[16]*M[2] + D[17]*M[3] + D[23]*M[4] + D[26]*M[5] + D[27]*M[6] + D[30]*M[7] + D[31]*M[8] + D[32]*M[9] + D[38]*M[10] + D[41]*M[11] + D[42]*M[12] + D[45]*M[13] + D[46]*M[14] + D[47]*M[15] + D[50]*M[16] + D[51]*M[17] + D[52]*M[18] + D[53]*M[19] + D[59]*M[20] + D[62]*M[21] + D[63]*M[22] + D[66]*M[23] + D[67]*M[24] + D[68]*M[25] + D[71]*M[26] + D[72]*M[27] + D[73]*M[28] + D[74]*M[29] + D[77]*M[30] + D[78]*M[31] + D[79]*M[32] + D[80]*M[33] + D[81]*M[34];
#pragma omp atomic
L[8] += D[8]*M[0] + D[14]*M[1] + D[17]*M[2] + D[18]*M[3] + D[24]*M[4] + D[27]*M[5] + D[28]*M[6] + D[31]*M[7] + D[32]*M[8] + D[33]*M[9] + D[39]*M[10] + D[42]*M[11] + D[43]*M[12] + D[46]*M[13] + D[47]*M[14] + D[48]*M[15] + D[51]*M[16] + D[52]*M[17] + D[53]*M[18] + D[54]*M[19] + D[60]*M[20] + D[63]*M[21] + D[64]*M[22] + D[67]*M[23] + D[68]*M[24] + D[69]*M[25] + D[72]*M[26] + D[73]*M[27] + D[74]*M[28] + D[75]*M[29] + D[78]*M[30] + D[79]*M[31] + D[80]*M[32] + D[81]*M[33] + D[82]*M[34];
#pragma omp atomic
L[9] += D[9]*M[0] + D[15]*M[1] + D[18]*M[2] + D[19]*M[3] + D[25]*M[4] + D[28]*M[5] + D[29]*M[6] + D[32]*M[7] + D[33]*M[8] + D[34]*M[9] + D[40]*M[10] + D[43]*M[11] + D[44]*M[12] + D[47]*M[13] + D[48]*M[14] + D[49]*M[15] + D[52]*M[16] + D[53]*M[17] + D[54]*M[18] + D[55]*M[19] + D[61]*M[20] + D[64]*M[21] + D[65]*M[22] + D[68]*M[23] + D[69]*M[24] + D[70]*M[25] + D[73]*M[26] + D[74]*M[27] + D[75]*M[28] + D[76]*M[29] + D[79]*M[30] + D[80]*M[31] + D[81]*M[32] + D[82]*M[33] + D[83]*M[34];
#pragma omp atomic
L[10] += D[10]*M[0] + D[20]*M[1] + D[21]*M[2] + D[22]*M[3] + D[35]*M[4] + D[36]*M[5] + D[37]*M[6] + D[38]*M[7] + D[39]*M[8] + D[40]*M[9] + D[56]*M[10] + D[57]*M[11] + D[58]*M[12] + D[59]*M[13] + D[60]*M[14] + D[61]*M[15] + D[62]*M[16] + D[63]*M[17] + D[64]*M[18] + D[65]*M[19];
#pragma omp atomic
L[11] += D[11]*M[0] + D[21]*M[1] + D[23]*M[2] + D[24]*M[3] + D[36]*M[4] + D[38]*M[5] + D[39]*M[6] + D[41]*M[7] + D[42]*M[8] + D[43]*M[9] + D[57]*M[10] + D[59]*M[11] + D[60]*M[12] + D[62]*M[13] + D[63]*M[14] + D[64]*M[15] + D[66]*M[16] + D[67]*M[17] + D[68]*M[18] + D[69]*M[19];
#pragma omp atomic
L[12] += D[12]*M[0] + D[22]*M[1] + D[24]*M[2] + D[25]*M[3] + D[37]*M[4] + D[39]*M[5] + D[40]*M[6] + D[42]*M[7] + D[43]*M[8] + D[44]*M[9] + D[58]*M[10] + D[60]*M[11] + D[61]*M[12] + D[63]*M[13] + D[64]*M[14] + D[65]*M[15] + D[67]*M[16] + D[68]*M[17] + D[69]*M[18] + D[70]*M[19];
#pragma omp atomic
L[13] += D[13]*M[0] + D[23]*M[1] + D[26]*M[2] + D[27]*M[3] + D[38]*M[4] + D[41]*M[5] + D[42]*M[6] + D[45]*M[7] + D[46]*M[8] + D[47]*M[9] + D[59]*M[10] + D[62]*M[11] + D[63]*M[12] + D[66]*M[13] + D[67]*M[14] + D[68]*M[15] + D[71]*M[16] + D[72]*M[17] + D[73]*M[18] + D[74]*M[19];
#pragma omp atomic
L[14] += D[14]*M[0] + D[24]*M[1] + D[27]*M[2] + D[28]*M[3] + D[39]*M[4] + D[42]*M[5] + D[43]*M[6] + D[46]*M[7] + D[47]*M[8] + D[48]*M[9] + D[60]*M[10] + D[63]*M[11] + D[64]*M[12] + D[67]*M[13] + D[68]*M[14] + D[69]*M[15] + D[72]*M[16] + D[73]*M[17] + D[74]*M[18] + D[75]*M[19];
#pragma omp atomic
L[15] += D[15]*M[0] + D[25]*M[1] + D[28]*M[2] + D[29]*M[3] + D[40]*M[4] + D[43]*M[5] + D[44]*M[6] + D[47]*M[7] + D[48]*M[8] + D[49]*M[9] + D[61]*M[10] + D[64]*M[11] + D[65]*M[12] + D[68]*M[13] + D[69]*M[14] + D[70]*M[15] + D[73]*M[16] + D[74]*M[17] + D[75]*M[18] + D[76]*M[19];
#pragma omp atomic
L[16] += D[16]*M[0] + D[26]*M[1] + D[30]*M[2] + D[31]*M[3] + D[41]*M[4] + D[45]*M[5] + D[46]*M[6] + D[50]*M[7] + D[51]*M[8] + D[52]*M[9] + D[62]*M[10] + D[66]*M[11] + D[67]*M[12] + D[71]*M[13] + D[72]*M[14] + D[73]*M[15] + D[77]*M[16] + D[78]*M[17] + D[79]*M[18] + D[80]*M[19];
#pragma omp atomic
L[17] += D[17]*M[0] + D[27]*M[1] + D[31]*M[2] + D[32]*M[3] + D[42]*M[4] + D[46]*M[5] + D[47]*M[6] + D[51]*M[7] + D[52]*M[8] + D[53]*M[9] + D[63]*M[10] + D[67]*M[11] + D[68]*M[12] + D[72]*M[13] + D[73]*M[14] + D[74]*M[15] + D[78]*M[16] + D[79]*M[17] + D[80]*M[18] + D[81]*M[19];
#pragma omp atomic
L[18] += D[18]*M[0] + D[28]*M[1] + D[32]*M[2] + D[33]*M[3] + D[43]*M[4] + D[47]*M[5] + D[48]*M[6] + D[52]*M[7] + D[53]*M[8] + D[54]*M[9] + D[64]*M[10] + D[68]*M[11] + D[69]*M[12] + D[73]*M[13] + D[74]*M[14] + D[75]*M[15] + D[79]*M[16] + D[80]*M[17] + D[81]*M[18] + D[82]*M[19];
#pragma omp atomic
L[19] += D[19]*M[0] + D[29]*M[1] + D[33]*M[2] + D[34]*M[3] + D[44]*M[4] + D[48]*M[5] + D[49]*M[6] + D[53]*M[7] + D[54]*M[8] + D[55]*M[9] + D[65]*M[10] + D[69]*M[11] + D[70]*M[12] + D[74]*M[13] + D[75]*M[14] + D[76]*M[15] + D[80]*M[16] + D[81]*M[17] + D[82]*M[18] + D[83]*M[19];
#pragma omp atomic
L[20] += D[20]*M[0] + D[35]*M[1] + D[36]*M[2] + D[37]*M[3] + D[56]*M[4] + D[57]*M[5] + D[58]*M[6] + D[59]*M[7] + D[60]*M[8] + D[61]*M[9];
#pragma omp atomic
L[21] += D[21]*M[0] + D[36]*M[1] + D[38]*M[2] + D[39]*M[3] + D[57]*M[4] + D[59]*M[5] + D[60]*M[6] + D[62]*M[7] + D[63]*M[8] + D[64]*M[9];
#pragma omp atomic
L[22] += D[22]*M[0] + D[37]*M[1] + D[39]*M[2] + D[40]*M[3] + D[58]*M[4] + D[60]*M[5] + D[61]*M[6] + D[63]*M[7] + D[64]*M[8] + D[65]*M[9];
#pragma omp atomic
L[23] += D[23]*M[0] + D[38]*M[1] + D[41]*M[2] + D[42]*M[3] + D[59]*M[4] + D[62]*M[5] + D[63]*M[6] + D[66]*M[7] + D[67]*M[8] + D[68]*M[9];
#pragma omp atomic
L[24] += D[24]*M[0] + D[39]*M[1] + D[42]*M[2] + D[43]*M[3] + D[60]*M[4] + D[63]*M[5] + D[64]*M[6] + D[67]*M[7] + D[68]*M[8] + D[69]*M[9];
#pragma omp atomic
L[25] += D[25]*M[0] + D[40]*M[1] + D[43]*M[2] + D[44]*M[3] + D[61]*M[4] + D[64]*M[5] + D[65]*M[6] + D[68]*M[7] + D[69]*M[8] + D[70]*M[9];
#pragma omp atomic
L[26] += D[26]*M[0] + D[41]*M[1] + D[45]*M[2] + D[46]*M[3] + D[62]*M[4] + D[66]*M[5] + D[67]*M[6] + D[71]*M[7] + D[72]*M[8] + D[73]*M[9];
#pragma omp atomic
L[27] += D[27]*M[0] + D[42]*M[1] + D[46]*M[2] + D[47]*M[3] + D[63]*M[4] + D[67]*M[5] + D[68]*M[6] + D[72]*M[7] + D[73]*M[8] + D[74]*M[9];
#pragma omp atomic
L[28] += D[28]*M[0] + D[43]*M[1] + D[47]*M[2] + D[48]*M[3] + D[64]*M[4] + D[68]*M[5] + D[69]*M[6] + D[73]*M[7] + D[74]*M[8] + D[75]*M[9];
#pragma omp atomic
L[29] += D[29]*M[0] + D[44]*M[1] + D[48]*M[2] + D[49]*M[3] + D[65]*M[4] + D[69]*M[5] + D[70]*M[6] + D[74]*M[7] + D[75]*M[8] + D[76]*M[9];
#pragma omp atomic
L[30] += D[30]*M[0] + D[45]*M[1] + D[50]*M[2] + D[51]*M[3] + D[66]*M[4] + D[71]*M[5] + D[72]*M[6] + D[77]*M[7] + D[78]*M[8] + D[79]*M[9];
#pragma omp atomic
L[31] += D[31]*M[0] + D[46]*M[1] + D[51]*M[2] + D[52]*M[3] + D[67]*M[4] + D[72]*M[5] + D[73]*M[6] + D[78]*M[7] + D[79]*M[8] + D[80]*M[9];
#pragma omp atomic
L[32] += D[32]*M[0] + D[47]*M[1] + D[52]*M[2] + D[53]*M[3] + D[68]*M[4] + D[73]*M[5] + D[74]*M[6] + D[79]*M[7] + D[80]*M[8] + D[81]*M[9];
#pragma omp atomic
L[33] += D[33]*M[0] + D[48]*M[1] + D[53]*M[2] + D[54]*M[3] + D[69]*M[4] + D[74]*M[5] + D[75]*M[6] + D[80]*M[7] + D[81]*M[8] + D[82]*M[9];
#pragma omp atomic
L[34] += D[34]*M[0] + D[49]*M[1] + D[54]*M[2] + D[55]*M[3] + D[70]*M[4] + D[75]*M[5] + D[76]*M[6] + D[81]*M[7] + D[82]*M[8] + D[83]*M[9];
#pragma omp atomic
L[35] += D[35]*M[0] + D[56]*M[1] + D[57]*M[2] + D[58]*M[3];
#pragma omp atomic
L[36] += D[36]*M[0] + D[57]*M[1] + D[59]*M[2] + D[60]*M[3];
#pragma omp atomic
L[37] += D[37]*M[0] + D[58]*M[1] + D[60]*M[2] + D[61]*M[3];
#pragma omp atomic
L[38] += D[38]*M[0] + D[59]*M[1] + D[62]*M[2] + D[63]*M[3];
#pragma omp atomic
L[39] += D[39]*M[0] + D[60]*M[1] + D[63]*M[2] + D[64]*M[3];
#pragma omp atomic
L[40] += D[40]*M[0] + D[61]*M[1] + D[64]*M[2] + D[65]*M[3];
#pragma omp atomic
L[41] += D[41]*M[0] + D[62]*M[1] + D[66]*M[2] + D[67]*M[3];
#pragma omp atomic
L[42] += D[42]*M[0] + D[63]*M[1] + D[67]*M[2] + D[68]*M[3];
#pragma omp atomic
L[43] += D[43]*M[0] + D[64]*M[1] + D[68]*M[2] + D[69]*M[3];
#pragma omp atomic
L[44] += D[44]*M[0] + D[65]*M[1] + D[69]*M[2] + D[70]*M[3];
#pragma omp atomic
L[45] += D[45]*M[0] + D[66]*M[1] + D[71]*M[2] + D[72]*M[3];
#pragma omp atomic
L[46] += D[46]*M[0] + D[67]*M[1] + D[72]*M[2] + D[73]*M[3];
#pragma omp atomic
L[47] += D[47]*M[0] + D[68]*M[1] + D[73]*M[2] + D[74]*M[3];
#pragma omp atomic
L[48] += D[48]*M[0] + D[69]*M[1] + D[74]*M[2] + D[75]*M[3];
#pragma omp atomic
L[49] += D[49]*M[0] + D[70]*M[1] + D[75]*M[2] + D[76]*M[3];
#pragma omp atomic
L[50] += D[50]*M[0] + D[71]*M[1] + D[77]*M[2] + D[78]*M[3];
#pragma omp atomic
L[51] += D[51]*M[0] + D[72]*M[1] + D[78]*M[2] + D[79]*M[3];
#pragma omp atomic
L[52] += D[52]*M[0] + D[73]*M[1] + D[79]*M[2] + D[80]*M[3];
#pragma omp atomic
L[53] += D[53]*M[0] + D[74]*M[1] + D[80]*M[2] + D[81]*M[3];
#pragma omp atomic
L[54] += D[54]*M[0] + D[75]*M[1] + D[81]*M[2] + D[82]*M[3];
#pragma omp atomic
L[55] += D[55]*M[0] + D[76]*M[1] + D[82]*M[2] + D[83]*M[3];
#pragma omp atomic
L[56] += D[56]*M[0];
#pragma omp atomic
L[57] += D[57]*M[0];
#pragma omp atomic
L[58] += D[58]*M[0];
#pragma omp atomic
L[59] += D[59]*M[0];
#pragma omp atomic
L[60] += D[60]*M[0];
#pragma omp atomic
L[61] += D[61]*M[0];
#pragma omp atomic
L[62] += D[62]*M[0];
#pragma omp atomic
L[63] += D[63]*M[0];
#pragma omp atomic
L[64] += D[64]*M[0];
#pragma omp atomic
L[65] += D[65]*M[0];
#pragma omp atomic
L[66] += D[66]*M[0];
#pragma omp atomic
L[67] += D[67]*M[0];
#pragma omp atomic
L[68] += D[68]*M[0];
#pragma omp atomic
L[69] += D[69]*M[0];
#pragma omp atomic
L[70] += D[70]*M[0];
#pragma omp atomic
L[71] += D[71]*M[0];
#pragma omp atomic
L[72] += D[72]*M[0];
#pragma omp atomic
L[73] += D[73]*M[0];
#pragma omp atomic
L[74] += D[74]*M[0];
#pragma omp atomic
L[75] += D[75]*M[0];
#pragma omp atomic
L[76] += D[76]*M[0];
#pragma omp atomic
L[77] += D[77]*M[0];
#pragma omp atomic
L[78] += D[78]*M[0];
#pragma omp atomic
L[79] += D[79]*M[0];
#pragma omp atomic
L[80] += D[80]*M[0];
#pragma omp atomic
L[81] += D[81]*M[0];
#pragma omp atomic
L[82] += D[82]*M[0];
#pragma omp atomic
L[83] += D[83]*M[0];

}

void field_m0_L2L_6(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
double Lstmp3 = z*L[14];
double Lstmp4 = Lstmp3*y;
double Lstmp5 = (x*x);
double Lstmp6 = (1.0/2.0)*Lstmp5;
double Lstmp7 = (x*x*x);
double Lstmp8 = (1.0/6.0)*Lstmp7;
double Lstmp9 = (x*x*x*x);
double Lstmp10 = (1.0/24.0)*Lstmp9;
double Lstmp11 = (1.0/120.0)*(x*x*x*x*x);
double Lstmp12 = (y*y);
double Lstmp13 = (1.0/2.0)*Lstmp12;
double Lstmp14 = (y*y*y);
double Lstmp15 = (1.0/6.0)*Lstmp14;
double Lstmp16 = (y*y*y*y);
double Lstmp17 = (1.0/24.0)*Lstmp16;
double Lstmp18 = (1.0/120.0)*(y*y*y*y*y);
double Lstmp19 = (z*z);
double Lstmp20 = (1.0/2.0)*Lstmp19;
double Lstmp21 = (z*z*z);
double Lstmp22 = (1.0/6.0)*Lstmp21;
double Lstmp23 = (z*z*z*z);
double Lstmp24 = (1.0/24.0)*Lstmp23;
double Lstmp25 = (1.0/120.0)*(z*z*z*z*z);
double Lstmp26 = x*L[13];
double Lstmp27 = x*L[26];
double Lstmp28 = x*L[45];
double Lstmp29 = x*L[71];
double Lstmp30 = x*L[15];
double Lstmp31 = x*L[29];
double Lstmp32 = x*L[49];
double Lstmp33 = x*L[76];
double Lstmp34 = y*L[11];
double Lstmp35 = z*L[12];
double Lstmp36 = y*L[21];
double Lstmp37 = z*L[22];
double Lstmp38 = y*L[36];
double Lstmp39 = z*L[37];
double Lstmp40 = y*L[57];
double Lstmp41 = z*L[58];
double Lstmp42 = y*L[18];
double Lstmp43 = y*L[33];
double Lstmp44 = y*L[54];
double Lstmp45 = y*L[82];
double Lstmp46 = z*L[17];
double Lstmp47 = z*L[31];
double Lstmp48 = z*L[51];
double Lstmp49 = z*L[78];
double Lstmp50 = y*L[28];
double Lstmp51 = Lstmp50*x;
double Lstmp52 = y*L[48];
double Lstmp53 = Lstmp52*x;
double Lstmp54 = y*L[75];
double Lstmp55 = Lstmp54*x;
double Lstmp56 = z*L[27];
double Lstmp57 = Lstmp56*x;
double Lstmp58 = z*L[46];
double Lstmp59 = Lstmp58*x;
double Lstmp60 = z*L[72];
double Lstmp61 = Lstmp60*x;
double Lstmp62 = z*L[24];
double Lstmp63 = Lstmp62*y;
double Lstmp64 = z*L[39];
double Lstmp65 = Lstmp64*y;
double Lstmp66 = z*L[60];
double Lstmp67 = Lstmp66*y;
double Lstmp68 = (1.0/4.0)*Lstmp5;
double Lstmp69 = Lstmp12*Lstmp68;
double Lstmp70 = (1.0/12.0)*Lstmp5;
double Lstmp71 = Lstmp14*Lstmp70;
double Lstmp72 = (1.0/48.0)*Lstmp5;
double Lstmp73 = Lstmp19*Lstmp68;
double Lstmp74 = Lstmp21*Lstmp70;
double Lstmp75 = (1.0/12.0)*Lstmp7;
double Lstmp76 = Lstmp12*Lstmp75;
double Lstmp77 = (1.0/36.0)*Lstmp7;
double Lstmp78 = Lstmp19*Lstmp75;
double Lstmp79 = (1.0/48.0)*Lstmp9;
double Lstmp80 = Lstmp12*Lstmp19;
double Lstmp81 = (1.0/4.0)*Lstmp80;
double Lstmp82 = (1.0/12.0)*Lstmp12*Lstmp21;
double Lstmp83 = (1.0/12.0)*Lstmp14*Lstmp19;
double Lstmp84 = x*L[47];
double Lstmp85 = x*L[74];
double Lstmp86 = x*L[73];
double Lstmp87 = y*L[43];
double Lstmp88 = y*L[69];
double Lstmp89 = z*L[42];
double Lstmp90 = z*L[67];
double Lstmp91 = y*L[64];
double Lstmp92 = z*L[63];
double Lstmp93 = x*L[23];
double Lstmp94 = x*L[41];
double Lstmp95 = x*L[66];
double Lstmp96 = x*L[25];
double Lstmp97 = x*L[44];
double Lstmp98 = x*L[70];
double Lstmp99 = Lstmp87*x;
double Lstmp100 = Lstmp88*x;
double Lstmp101 = Lstmp89*x;
double Lstmp102 = Lstmp90*x;
double Lstmp103 = x*L[68];
double Lstmp104 = y*L[13];
double Lstmp105 = Lstmp56*y;
double Lstmp106 = x*L[28];
double Lstmp107 = x*L[48];
double Lstmp108 = x*L[75];
double Lstmp109 = y*L[23];
double Lstmp110 = y*L[38];
double Lstmp111 = y*L[59];
double Lstmp112 = y*L[32];
double Lstmp113 = y*L[53];
double Lstmp114 = y*L[81];
double Lstmp115 = y*L[47];
double Lstmp116 = Lstmp115*x;
double Lstmp117 = y*L[74];
double Lstmp118 = Lstmp117*x;
double Lstmp119 = Lstmp89*y;
double Lstmp120 = Lstmp92*y;
double Lstmp121 = y*L[68];
double Lstmp122 = y*L[14];
double Lstmp123 = z*L[15];
double Lstmp124 = z*L[18];
double Lstmp125 = z*L[28];
double Lstmp126 = Lstmp125*y;
double Lstmp127 = x*L[27];
double Lstmp128 = x*L[46];
double Lstmp129 = x*L[72];
double Lstmp130 = y*L[24];
double Lstmp131 = z*L[25];
double Lstmp132 = y*L[39];
double Lstmp133 = z*L[40];
double Lstmp134 = y*L[60];
double Lstmp135 = z*L[61];
double Lstmp136 = z*L[32];
double Lstmp137 = z*L[52];
double Lstmp138 = z*L[79];
double Lstmp139 = z*L[47];
double Lstmp140 = Lstmp139*x;
double Lstmp141 = z*L[73];
double Lstmp142 = Lstmp141*x;
double Lstmp143 = z*L[43];
double Lstmp144 = Lstmp143*y;
double Lstmp145 = z*L[64];
double Lstmp146 = Lstmp145*y;
double Lstmp147 = z*L[68];
double Lstmp148 = x*L[38];
double Lstmp149 = x*L[62];
double Lstmp150 = x*L[40];
double Lstmp151 = x*L[65];
double Lstmp152 = Lstmp91*x;
double Lstmp153 = Lstmp92*x;
double Lstmp154 = x*L[43];
double Lstmp155 = x*L[69];
double Lstmp156 = Lstmp121*x;
double Lstmp157 = x*L[42];
double Lstmp158 = x*L[67];
double Lstmp159 = Lstmp147*x;
double Lstmp160 = y*L[26];
double Lstmp161 = Lstmp58*y;
double Lstmp162 = y*L[41];
double Lstmp163 = y*L[62];
double Lstmp164 = y*L[52];
double Lstmp165 = y*L[80];
double Lstmp166 = y*L[73];
double Lstmp167 = Lstmp166*x;
double Lstmp168 = Lstmp90*y;
double Lstmp169 = y*L[27];
double Lstmp170 = Lstmp139*y;
double Lstmp171 = y*L[42];
double Lstmp172 = y*L[63];
double Lstmp173 = Lstmp147*y;
double Lstmp174 = z*L[29];
double Lstmp175 = z*L[33];
double Lstmp176 = z*L[48];
double Lstmp177 = Lstmp176*y;
double Lstmp178 = z*L[44];
double Lstmp179 = z*L[65];
double Lstmp180 = z*L[53];
double Lstmp181 = z*L[80];
double Lstmp182 = z*L[74];
double Lstmp183 = Lstmp182*x;
double Lstmp184 = z*L[69];
double Lstmp185 = Lstmp184*y;
double Lstmp186 = x*L[59];
double Lstmp187 = x*L[61];
double Lstmp188 = x*L[64];
double Lstmp189 = x*L[63];
double Lstmp190 = y*L[45];
double Lstmp191 = Lstmp60*y;
double Lstmp192 = y*L[66];
double Lstmp193 = y*L[79];
double Lstmp194 = y*L[46];
double Lstmp195 = Lstmp141*y;
double Lstmp196 = y*L[67];
double Lstmp197 = Lstmp182*y;
double Lstmp198 = z*L[49];
double Lstmp199 = z*L[54];
double Lstmp200 = z*L[75];
double Lstmp201 = Lstmp200*y;
double Lstmp202 = z*L[70];
double Lstmp203 = z*L[81];
double Lstmp204 = y*L[71];
double Lstmp205 = y*L[72];
double Lstmp206 = z*L[76];
double Lstmp207 = z*L[82];
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x + Lstmp10*Lstmp38 + Lstmp10*Lstmp39 + Lstmp10*Lstmp67 + Lstmp10*L[20] + Lstmp11*Lstmp40 + Lstmp11*Lstmp41 + Lstmp11*L[35] + (1.0/48.0)*Lstmp12*Lstmp23*L[81] + Lstmp12*Lstmp79*L[59] + Lstmp13*Lstmp26 + Lstmp13*Lstmp46 + Lstmp13*Lstmp57 + Lstmp13*L[7] + (1.0/36.0)*Lstmp14*Lstmp21*L[80] + Lstmp14*Lstmp77*L[62] + Lstmp15*Lstmp27 + Lstmp15*Lstmp47 + Lstmp15*Lstmp59 + Lstmp15*L[16] + (1.0/48.0)*Lstmp16*Lstmp19*L[79] + Lstmp16*Lstmp72*L[66] + Lstmp17*Lstmp28 + Lstmp17*Lstmp48 + Lstmp17*Lstmp61 + Lstmp17*L[30] + Lstmp18*Lstmp29 + Lstmp18*Lstmp49 + Lstmp18*L[50] + Lstmp19*Lstmp79*L[61] + Lstmp2*y + Lstmp20*Lstmp30 + Lstmp20*Lstmp42 + Lstmp20*Lstmp51 + Lstmp20*L[9] + Lstmp21*Lstmp77*L[65] + Lstmp22*Lstmp31 + Lstmp22*Lstmp43 + Lstmp22*Lstmp53 + Lstmp22*L[19] + Lstmp23*Lstmp72*L[70] + Lstmp24*Lstmp32 + Lstmp24*Lstmp44 + Lstmp24*Lstmp55 + Lstmp24*L[34] + Lstmp25*Lstmp33 + Lstmp25*Lstmp45 + Lstmp25*L[55] + Lstmp34*Lstmp6 + Lstmp35*Lstmp6 + Lstmp36*Lstmp8 + Lstmp37*Lstmp8 + Lstmp4*x + (1.0/8.0)*Lstmp5*Lstmp80*L[68] + Lstmp6*Lstmp63 + Lstmp6*L[4] + Lstmp65*Lstmp8 + Lstmp69*Lstmp89 + Lstmp69*L[23] + Lstmp71*Lstmp90 + Lstmp71*L[41] + Lstmp73*Lstmp87 + Lstmp73*L[25] + Lstmp74*Lstmp88 + Lstmp74*L[44] + Lstmp76*Lstmp92 + Lstmp76*L[38] + Lstmp78*Lstmp91 + Lstmp78*L[40] + Lstmp8*L[10] + Lstmp81*Lstmp84 + Lstmp81*L[32] + Lstmp82*Lstmp85 + Lstmp82*L[53] + Lstmp83*Lstmp86 + Lstmp83*L[52] + (1.0/720.0)*(x*x*x*x*x*x)*L[56] + x*L[1] + (1.0/720.0)*(y*y*y*y*y*y)*L[77] + y*L[2] + (1.0/720.0)*(z*z*z*z*z*z)*L[83] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 + Lstmp10*Lstmp40 + Lstmp10*Lstmp41 + Lstmp10*L[35] + Lstmp100*Lstmp22 + Lstmp101*Lstmp13 + Lstmp102*Lstmp15 + Lstmp103*Lstmp81 + Lstmp11*L[56] + Lstmp13*Lstmp56 + Lstmp13*Lstmp93 + Lstmp13*L[13] + Lstmp15*Lstmp58 + Lstmp15*Lstmp94 + Lstmp15*L[26] + Lstmp17*Lstmp60 + Lstmp17*Lstmp95 + Lstmp17*L[45] + Lstmp18*L[71] + Lstmp20*Lstmp50 + Lstmp20*Lstmp96 + Lstmp20*Lstmp99 + Lstmp20*L[15] + Lstmp22*Lstmp52 + Lstmp22*Lstmp97 + Lstmp22*L[29] + Lstmp24*Lstmp54 + Lstmp24*Lstmp98 + Lstmp24*L[49] + Lstmp25*L[76] + Lstmp34*x + Lstmp35*x + Lstmp36*Lstmp6 + Lstmp37*Lstmp6 + Lstmp38*Lstmp8 + Lstmp39*Lstmp8 + Lstmp4 + Lstmp6*Lstmp65 + Lstmp6*L[10] + Lstmp63*x + Lstmp67*Lstmp8 + Lstmp69*Lstmp92 + Lstmp69*L[38] + Lstmp71*L[62] + Lstmp73*Lstmp91 + Lstmp73*L[40] + Lstmp74*L[65] + Lstmp76*L[59] + Lstmp78*L[61] + Lstmp8*L[20] + Lstmp81*L[47] + Lstmp82*L[74] + Lstmp83*L[73] + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp10*Lstmp111 + Lstmp10*Lstmp66 + Lstmp10*L[36] + Lstmp104*x + Lstmp105*x + Lstmp106*Lstmp20 + Lstmp107*Lstmp22 + Lstmp108*Lstmp24 + Lstmp109*Lstmp6 + Lstmp11*L[57] + Lstmp110*Lstmp8 + Lstmp112*Lstmp20 + Lstmp113*Lstmp22 + Lstmp114*Lstmp24 + Lstmp116*Lstmp20 + Lstmp118*Lstmp22 + Lstmp119*Lstmp6 + Lstmp120*Lstmp8 + Lstmp121*Lstmp73 + Lstmp13*Lstmp27 + Lstmp13*Lstmp47 + Lstmp13*Lstmp59 + Lstmp13*L[16] + Lstmp15*Lstmp28 + Lstmp15*Lstmp48 + Lstmp15*Lstmp61 + Lstmp15*L[30] + Lstmp17*Lstmp29 + Lstmp17*Lstmp49 + Lstmp17*L[50] + Lstmp18*L[77] + Lstmp2 + Lstmp20*L[18] + Lstmp22*L[33] + Lstmp24*L[54] + Lstmp25*L[82] + Lstmp3*x + Lstmp46*y + Lstmp6*Lstmp62 + Lstmp6*L[11] + Lstmp64*Lstmp8 + Lstmp69*Lstmp90 + Lstmp69*L[41] + Lstmp71*L[66] + Lstmp73*L[43] + Lstmp74*L[69] + Lstmp76*L[62] + Lstmp78*L[64] + Lstmp8*L[21] + Lstmp81*Lstmp86 + Lstmp81*L[52] + Lstmp82*L[80] + Lstmp83*L[79] + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += Lstmp10*Lstmp134 + Lstmp10*Lstmp135 + Lstmp10*L[37] + Lstmp11*L[58] + Lstmp122*x + Lstmp123*x + Lstmp124*y + Lstmp126*x + Lstmp127*Lstmp13 + Lstmp128*Lstmp15 + Lstmp129*Lstmp17 + Lstmp13*Lstmp136 + Lstmp13*Lstmp140 + Lstmp13*L[17] + Lstmp130*Lstmp6 + Lstmp131*Lstmp6 + Lstmp132*Lstmp8 + Lstmp133*Lstmp8 + Lstmp137*Lstmp15 + Lstmp138*Lstmp17 + Lstmp142*Lstmp15 + Lstmp144*Lstmp6 + Lstmp146*Lstmp8 + Lstmp147*Lstmp69 + Lstmp15*L[31] + Lstmp17*L[51] + Lstmp18*L[78] + Lstmp20*Lstmp31 + Lstmp20*Lstmp43 + Lstmp20*Lstmp53 + Lstmp20*L[19] + Lstmp22*Lstmp32 + Lstmp22*Lstmp44 + Lstmp22*Lstmp55 + Lstmp22*L[34] + Lstmp24*Lstmp33 + Lstmp24*Lstmp45 + Lstmp24*L[55] + Lstmp25*L[83] + Lstmp6*L[12] + Lstmp69*L[42] + Lstmp71*L[67] + Lstmp73*Lstmp88 + Lstmp73*L[44] + Lstmp74*L[70] + Lstmp76*L[63] + Lstmp78*L[65] + Lstmp8*L[22] + Lstmp81*Lstmp85 + Lstmp81*L[53] + Lstmp82*L[81] + Lstmp83*L[80] + x*L[6] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += Lstmp10*L[56] + Lstmp13*Lstmp148 + Lstmp13*Lstmp153 + Lstmp13*Lstmp89 + Lstmp13*L[23] + Lstmp149*Lstmp15 + Lstmp15*Lstmp90 + Lstmp15*L[41] + Lstmp150*Lstmp20 + Lstmp151*Lstmp22 + Lstmp152*Lstmp20 + Lstmp17*L[66] + Lstmp20*Lstmp87 + Lstmp20*L[25] + Lstmp22*Lstmp88 + Lstmp22*L[44] + Lstmp24*L[70] + Lstmp34 + Lstmp35 + Lstmp36*x + Lstmp37*x + Lstmp38*Lstmp6 + Lstmp39*Lstmp6 + Lstmp40*Lstmp8 + Lstmp41*Lstmp8 + Lstmp6*Lstmp67 + Lstmp6*L[20] + Lstmp63 + Lstmp65*x + Lstmp69*L[59] + Lstmp73*L[61] + Lstmp8*L[35] + Lstmp81*L[68] + x*L[10] + L[4];
#pragma omp atomic
Ls[5] += Lstmp10*L[57] + Lstmp102*Lstmp13 + Lstmp104 + Lstmp105 + Lstmp109*x + Lstmp110*Lstmp6 + Lstmp111*Lstmp8 + Lstmp115*Lstmp20 + Lstmp117*Lstmp22 + Lstmp119*x + Lstmp120*Lstmp6 + Lstmp13*Lstmp58 + Lstmp13*Lstmp94 + Lstmp13*L[26] + Lstmp15*Lstmp60 + Lstmp15*Lstmp95 + Lstmp15*L[45] + Lstmp154*Lstmp20 + Lstmp155*Lstmp22 + Lstmp156*Lstmp20 + Lstmp17*L[71] + Lstmp20*L[28] + Lstmp22*L[48] + Lstmp24*L[75] + Lstmp3 + Lstmp6*Lstmp64 + Lstmp6*L[21] + Lstmp62*x + Lstmp66*Lstmp8 + Lstmp69*L[62] + Lstmp73*L[64] + Lstmp8*L[36] + Lstmp81*L[73] + x*L[11] + L[5];
#pragma omp atomic
Ls[6] += Lstmp10*L[58] + Lstmp100*Lstmp20 + Lstmp122 + Lstmp123 + Lstmp126 + Lstmp13*Lstmp139 + Lstmp13*Lstmp157 + Lstmp13*Lstmp159 + Lstmp13*L[27] + Lstmp130*x + Lstmp131*x + Lstmp132*Lstmp6 + Lstmp133*Lstmp6 + Lstmp134*Lstmp8 + Lstmp135*Lstmp8 + Lstmp141*Lstmp15 + Lstmp144*x + Lstmp146*Lstmp6 + Lstmp15*Lstmp158 + Lstmp15*L[46] + Lstmp17*L[72] + Lstmp20*Lstmp52 + Lstmp20*Lstmp97 + Lstmp20*L[29] + Lstmp22*Lstmp54 + Lstmp22*Lstmp98 + Lstmp22*L[49] + Lstmp24*L[76] + Lstmp6*L[22] + Lstmp69*L[63] + Lstmp73*L[65] + Lstmp8*L[37] + Lstmp81*L[74] + x*L[12] + L[6];
#pragma omp atomic
Ls[7] += Lstmp10*L[59] + Lstmp13*Lstmp28 + Lstmp13*Lstmp48 + Lstmp13*Lstmp61 + Lstmp13*L[30] + Lstmp15*Lstmp29 + Lstmp15*Lstmp49 + Lstmp15*L[50] + Lstmp160*x + Lstmp161*x + Lstmp162*Lstmp6 + Lstmp163*Lstmp8 + Lstmp164*Lstmp20 + Lstmp165*Lstmp22 + Lstmp167*Lstmp20 + Lstmp168*Lstmp6 + Lstmp17*L[77] + Lstmp20*Lstmp84 + Lstmp20*L[32] + Lstmp22*Lstmp85 + Lstmp22*L[53] + Lstmp24*L[81] + Lstmp26 + Lstmp46 + Lstmp47*y + Lstmp57 + Lstmp6*Lstmp89 + Lstmp6*L[23] + Lstmp69*L[66] + Lstmp73*L[68] + Lstmp8*Lstmp92 + Lstmp8*L[38] + Lstmp81*L[79] + y*L[16] + L[7];
#pragma omp atomic
Ls[8] += Lstmp10*L[60] + Lstmp107*Lstmp20 + Lstmp108*Lstmp22 + Lstmp113*Lstmp20 + Lstmp114*Lstmp22 + Lstmp118*Lstmp20 + Lstmp124 + Lstmp125*x + Lstmp128*Lstmp13 + Lstmp129*Lstmp15 + Lstmp13*Lstmp137 + Lstmp13*Lstmp142 + Lstmp13*L[31] + Lstmp136*y + Lstmp138*Lstmp15 + Lstmp143*Lstmp6 + Lstmp145*Lstmp8 + Lstmp15*L[51] + Lstmp169*x + Lstmp17*L[78] + Lstmp170*x + Lstmp171*Lstmp6 + Lstmp172*Lstmp8 + Lstmp173*Lstmp6 + Lstmp20*L[33] + Lstmp22*L[54] + Lstmp24*L[82] + Lstmp6*L[24] + Lstmp69*L[67] + Lstmp73*L[69] + Lstmp8*L[39] + Lstmp81*L[80] + x*L[14] + y*L[17] + L[8];
#pragma omp atomic
Ls[9] += Lstmp10*L[61] + Lstmp13*Lstmp180 + Lstmp13*Lstmp183 + Lstmp13*Lstmp84 + Lstmp13*L[32] + Lstmp15*Lstmp181 + Lstmp15*Lstmp86 + Lstmp15*L[52] + Lstmp17*L[79] + Lstmp174*x + Lstmp175*y + Lstmp177*x + Lstmp178*Lstmp6 + Lstmp179*Lstmp8 + Lstmp185*Lstmp6 + Lstmp20*Lstmp32 + Lstmp20*Lstmp44 + Lstmp20*Lstmp55 + Lstmp20*L[34] + Lstmp22*Lstmp33 + Lstmp22*Lstmp45 + Lstmp22*L[55] + Lstmp24*L[83] + Lstmp30 + Lstmp42 + Lstmp51 + Lstmp6*Lstmp87 + Lstmp6*L[25] + Lstmp69*L[68] + Lstmp73*L[70] + Lstmp8*Lstmp91 + Lstmp8*L[40] + Lstmp81*L[81] + z*L[19] + L[9];
#pragma omp atomic
Ls[10] += Lstmp13*Lstmp186 + Lstmp13*Lstmp92 + Lstmp13*L[38] + Lstmp15*L[62] + Lstmp187*Lstmp20 + Lstmp20*Lstmp91 + Lstmp20*L[40] + Lstmp22*L[65] + Lstmp36 + Lstmp37 + Lstmp38*x + Lstmp39*x + Lstmp40*Lstmp6 + Lstmp41*Lstmp6 + Lstmp6*L[35] + Lstmp65 + Lstmp67*x + Lstmp8*L[56] + x*L[20] + L[10];
#pragma omp atomic
Ls[11] += Lstmp109 + Lstmp110*x + Lstmp111*Lstmp6 + Lstmp119 + Lstmp120*x + Lstmp121*Lstmp20 + Lstmp13*Lstmp149 + Lstmp13*Lstmp90 + Lstmp13*L[41] + Lstmp15*L[66] + Lstmp188*Lstmp20 + Lstmp20*L[43] + Lstmp22*L[69] + Lstmp6*Lstmp66 + Lstmp6*L[36] + Lstmp62 + Lstmp64*x + Lstmp8*L[57] + x*L[21] + L[11];
#pragma omp atomic
Ls[12] += Lstmp13*Lstmp147 + Lstmp13*Lstmp189 + Lstmp13*L[42] + Lstmp130 + Lstmp131 + Lstmp132*x + Lstmp133*x + Lstmp134*Lstmp6 + Lstmp135*Lstmp6 + Lstmp144 + Lstmp146*x + Lstmp15*L[67] + Lstmp151*Lstmp20 + Lstmp20*Lstmp88 + Lstmp20*L[44] + Lstmp22*L[70] + Lstmp6*L[37] + Lstmp8*L[58] + x*L[22] + L[12];
#pragma omp atomic
Ls[13] += Lstmp101 + Lstmp103*Lstmp20 + Lstmp13*Lstmp60 + Lstmp13*Lstmp95 + Lstmp13*L[45] + Lstmp15*L[71] + Lstmp160 + Lstmp161 + Lstmp162*x + Lstmp163*Lstmp6 + Lstmp166*Lstmp20 + Lstmp168*x + Lstmp20*L[47] + Lstmp22*L[74] + Lstmp56 + Lstmp6*Lstmp92 + Lstmp6*L[38] + Lstmp8*L[59] + Lstmp93 + L[13];
#pragma omp atomic
Ls[14] += Lstmp117*Lstmp20 + Lstmp125 + Lstmp13*Lstmp141 + Lstmp13*Lstmp158 + Lstmp13*L[46] + Lstmp143*x + Lstmp145*Lstmp6 + Lstmp15*L[72] + Lstmp155*Lstmp20 + Lstmp169 + Lstmp170 + Lstmp171*x + Lstmp172*Lstmp6 + Lstmp173*x + Lstmp20*L[48] + Lstmp22*L[75] + Lstmp6*L[39] + Lstmp8*L[60] + x*L[24] + L[14];
#pragma omp atomic
Ls[15] += Lstmp103*Lstmp13 + Lstmp13*Lstmp182 + Lstmp13*L[47] + Lstmp15*L[73] + Lstmp174 + Lstmp177 + Lstmp178*x + Lstmp179*Lstmp6 + Lstmp185*x + Lstmp20*Lstmp54 + Lstmp20*Lstmp98 + Lstmp20*L[49] + Lstmp22*L[76] + Lstmp50 + Lstmp6*Lstmp91 + Lstmp6*L[40] + Lstmp8*L[61] + Lstmp96 + Lstmp99 + L[15];
#pragma omp atomic
Ls[16] += Lstmp13*Lstmp29 + Lstmp13*Lstmp49 + Lstmp13*L[50] + Lstmp15*L[77] + Lstmp190*x + Lstmp191*x + Lstmp192*Lstmp6 + Lstmp193*Lstmp20 + Lstmp20*Lstmp86 + Lstmp20*L[52] + Lstmp22*L[80] + Lstmp27 + Lstmp47 + Lstmp48*y + Lstmp59 + Lstmp6*Lstmp90 + Lstmp6*L[41] + Lstmp8*L[62] + y*L[30] + L[16];
#pragma omp atomic
Ls[17] += Lstmp127 + Lstmp129*Lstmp13 + Lstmp13*Lstmp138 + Lstmp13*L[51] + Lstmp136 + Lstmp137*y + Lstmp140 + Lstmp147*Lstmp6 + Lstmp15*L[78] + Lstmp165*Lstmp20 + Lstmp194*x + Lstmp195*x + Lstmp196*Lstmp6 + Lstmp20*Lstmp85 + Lstmp20*L[53] + Lstmp22*L[81] + Lstmp6*L[42] + Lstmp8*L[63] + y*L[31] + L[17];
#pragma omp atomic
Ls[18] += Lstmp106 + Lstmp108*Lstmp20 + Lstmp112 + Lstmp114*Lstmp20 + Lstmp116 + Lstmp121*Lstmp6 + Lstmp13*Lstmp181 + Lstmp13*Lstmp86 + Lstmp13*L[52] + Lstmp15*L[79] + Lstmp175 + Lstmp176*x + Lstmp180*y + Lstmp184*Lstmp6 + Lstmp197*x + Lstmp20*L[54] + Lstmp22*L[82] + Lstmp6*L[43] + Lstmp8*L[64] + L[18];
#pragma omp atomic
Ls[19] += Lstmp13*Lstmp203 + Lstmp13*Lstmp85 + Lstmp13*L[53] + Lstmp15*L[80] + Lstmp198*x + Lstmp199*y + Lstmp20*Lstmp33 + Lstmp20*Lstmp45 + Lstmp20*L[55] + Lstmp201*x + Lstmp202*Lstmp6 + Lstmp22*L[83] + Lstmp31 + Lstmp43 + Lstmp53 + Lstmp6*Lstmp88 + Lstmp6*L[44] + Lstmp8*L[65] + z*L[34] + L[19];
#pragma omp atomic
Ls[20] += Lstmp13*L[59] + Lstmp20*L[61] + Lstmp38 + Lstmp39 + Lstmp40*x + Lstmp41*x + Lstmp6*L[56] + Lstmp67 + x*L[35] + L[20];
#pragma omp atomic
Ls[21] += Lstmp110 + Lstmp111*x + Lstmp120 + Lstmp13*L[62] + Lstmp20*L[64] + Lstmp6*L[57] + Lstmp64 + Lstmp66*x + x*L[36] + L[21];
#pragma omp atomic
Ls[22] += Lstmp13*L[63] + Lstmp132 + Lstmp133 + Lstmp134*x + Lstmp135*x + Lstmp146 + Lstmp20*L[65] + Lstmp6*L[58] + x*L[37] + L[22];
#pragma omp atomic
Ls[23] += Lstmp13*L[66] + Lstmp148 + Lstmp153 + Lstmp162 + Lstmp163*x + Lstmp168 + Lstmp20*L[68] + Lstmp6*L[59] + Lstmp89 + L[23];
#pragma omp atomic
Ls[24] += Lstmp13*L[67] + Lstmp143 + Lstmp145*x + Lstmp171 + Lstmp172*x + Lstmp173 + Lstmp20*L[69] + Lstmp6*L[60] + x*L[39] + L[24];
#pragma omp atomic
Ls[25] += Lstmp13*L[68] + Lstmp150 + Lstmp152 + Lstmp178 + Lstmp179*x + Lstmp185 + Lstmp20*L[70] + Lstmp6*L[61] + Lstmp87 + L[25];
#pragma omp atomic
Ls[26] += Lstmp102 + Lstmp13*L[71] + Lstmp190 + Lstmp191 + Lstmp192*x + Lstmp20*L[73] + Lstmp58 + Lstmp6*L[62] + Lstmp94 + L[26];
#pragma omp atomic
Ls[27] += Lstmp13*L[72] + Lstmp139 + Lstmp157 + Lstmp159 + Lstmp194 + Lstmp195 + Lstmp196*x + Lstmp20*L[74] + Lstmp6*L[63] + L[27];
#pragma omp atomic
Ls[28] += Lstmp115 + Lstmp13*L[73] + Lstmp154 + Lstmp156 + Lstmp176 + Lstmp184*x + Lstmp197 + Lstmp20*L[75] + Lstmp6*L[64] + L[28];
#pragma omp atomic
Ls[29] += Lstmp100 + Lstmp13*L[74] + Lstmp198 + Lstmp20*L[76] + Lstmp201 + Lstmp202*x + Lstmp52 + Lstmp6*L[65] + Lstmp97 + L[29];
#pragma omp atomic
Ls[30] += Lstmp13*L[77] + Lstmp20*L[79] + Lstmp204*x + Lstmp28 + Lstmp48 + Lstmp49*y + Lstmp6*L[66] + Lstmp61 + y*L[50] + L[30];
#pragma omp atomic
Ls[31] += Lstmp128 + Lstmp13*L[78] + Lstmp137 + Lstmp138*y + Lstmp142 + Lstmp20*L[80] + Lstmp205*x + Lstmp6*L[67] + y*L[51] + L[31];
#pragma omp atomic
Ls[32] += Lstmp13*L[79] + Lstmp164 + Lstmp167 + Lstmp180 + Lstmp181*y + Lstmp183 + Lstmp20*L[81] + Lstmp6*L[68] + Lstmp84 + L[32];
#pragma omp atomic
Ls[33] += Lstmp107 + Lstmp113 + Lstmp118 + Lstmp13*L[80] + Lstmp199 + Lstmp20*L[82] + Lstmp200*x + Lstmp203*y + Lstmp6*L[69] + L[33];
#pragma omp atomic
Ls[34] += Lstmp13*L[81] + Lstmp20*L[83] + Lstmp206*x + Lstmp207*y + Lstmp32 + Lstmp44 + Lstmp55 + Lstmp6*L[70] + z*L[55] + L[34];
#pragma omp atomic
Ls[35] += Lstmp40 + Lstmp41 + x*L[56] + L[35];
#pragma omp atomic
Ls[36] += Lstmp111 + Lstmp66 + x*L[57] + L[36];
#pragma omp atomic
Ls[37] += Lstmp134 + Lstmp135 + x*L[58] + L[37];
#pragma omp atomic
Ls[38] += Lstmp163 + Lstmp186 + Lstmp92 + L[38];
#pragma omp atomic
Ls[39] += Lstmp145 + Lstmp172 + x*L[60] + L[39];
#pragma omp atomic
Ls[40] += Lstmp179 + Lstmp187 + Lstmp91 + L[40];
#pragma omp atomic
Ls[41] += Lstmp149 + Lstmp192 + Lstmp90 + L[41];
#pragma omp atomic
Ls[42] += Lstmp147 + Lstmp189 + Lstmp196 + L[42];
#pragma omp atomic
Ls[43] += Lstmp121 + Lstmp184 + Lstmp188 + L[43];
#pragma omp atomic
Ls[44] += Lstmp151 + Lstmp202 + Lstmp88 + L[44];
#pragma omp atomic
Ls[45] += Lstmp204 + Lstmp60 + Lstmp95 + L[45];
#pragma omp atomic
Ls[46] += Lstmp141 + Lstmp158 + Lstmp205 + L[46];
#pragma omp atomic
Ls[47] += Lstmp103 + Lstmp166 + Lstmp182 + L[47];
#pragma omp atomic
Ls[48] += Lstmp117 + Lstmp155 + Lstmp200 + L[48];
#pragma omp atomic
Ls[49] += Lstmp206 + Lstmp54 + Lstmp98 + L[49];
#pragma omp atomic
Ls[50] += Lstmp29 + Lstmp49 + y*L[77] + L[50];
#pragma omp atomic
Ls[51] += Lstmp129 + Lstmp138 + y*L[78] + L[51];
#pragma omp atomic
Ls[52] += Lstmp181 + Lstmp193 + Lstmp86 + L[52];
#pragma omp atomic
Ls[53] += Lstmp165 + Lstmp203 + Lstmp85 + L[53];
#pragma omp atomic
Ls[54] += Lstmp108 + Lstmp114 + Lstmp207 + L[54];
#pragma omp atomic
Ls[55] += Lstmp33 + Lstmp45 + z*L[83] + L[55];
#pragma omp atomic
Ls[56] += L[56];
#pragma omp atomic
Ls[57] += L[57];
#pragma omp atomic
Ls[58] += L[58];
#pragma omp atomic
Ls[59] += L[59];
#pragma omp atomic
Ls[60] += L[60];
#pragma omp atomic
Ls[61] += L[61];
#pragma omp atomic
Ls[62] += L[62];
#pragma omp atomic
Ls[63] += L[63];
#pragma omp atomic
Ls[64] += L[64];
#pragma omp atomic
Ls[65] += L[65];
#pragma omp atomic
Ls[66] += L[66];
#pragma omp atomic
Ls[67] += L[67];
#pragma omp atomic
Ls[68] += L[68];
#pragma omp atomic
Ls[69] += L[69];
#pragma omp atomic
Ls[70] += L[70];
#pragma omp atomic
Ls[71] += L[71];
#pragma omp atomic
Ls[72] += L[72];
#pragma omp atomic
Ls[73] += L[73];
#pragma omp atomic
Ls[74] += L[74];
#pragma omp atomic
Ls[75] += L[75];
#pragma omp atomic
Ls[76] += L[76];
#pragma omp atomic
Ls[77] += L[77];
#pragma omp atomic
Ls[78] += L[78];
#pragma omp atomic
Ls[79] += L[79];
#pragma omp atomic
Ls[80] += L[80];
#pragma omp atomic
Ls[81] += L[81];
#pragma omp atomic
Ls[82] += L[82];
#pragma omp atomic
Ls[83] += L[83];

}

void field_m0_L2P_6(double x, double y, double z, double * L, double * F) {
double Ftmp0 = x*y;
double Ftmp1 = x*z;
double Ftmp2 = y*z;
double Ftmp3 = Ftmp0*z;
double Ftmp4 = (x*x);
double Ftmp5 = (1.0/2.0)*Ftmp4;
double Ftmp6 = (x*x*x);
double Ftmp7 = (1.0/6.0)*Ftmp6;
double Ftmp8 = (1.0/24.0)*(x*x*x*x);
double Ftmp9 = (1.0/120.0)*(x*x*x*x*x);
double Ftmp10 = (y*y);
double Ftmp11 = (1.0/2.0)*Ftmp10;
double Ftmp12 = (y*y*y);
double Ftmp13 = (1.0/6.0)*Ftmp12;
double Ftmp14 = (1.0/24.0)*(y*y*y*y);
double Ftmp15 = (1.0/120.0)*(y*y*y*y*y);
double Ftmp16 = (z*z);
double Ftmp17 = (1.0/2.0)*Ftmp16;
double Ftmp18 = (z*z*z);
double Ftmp19 = (1.0/6.0)*Ftmp18;
double Ftmp20 = (1.0/24.0)*(z*z*z*z);
double Ftmp21 = (1.0/120.0)*(z*z*z*z*z);
double Ftmp22 = Ftmp11*x;
double Ftmp23 = Ftmp13*x;
double Ftmp24 = Ftmp14*x;
double Ftmp25 = Ftmp17*x;
double Ftmp26 = Ftmp19*x;
double Ftmp27 = Ftmp20*x;
double Ftmp28 = Ftmp5*y;
double Ftmp29 = Ftmp5*z;
double Ftmp30 = Ftmp7*y;
double Ftmp31 = Ftmp7*z;
double Ftmp32 = Ftmp8*y;
double Ftmp33 = Ftmp8*z;
double Ftmp34 = Ftmp17*y;
double Ftmp35 = Ftmp19*y;
double Ftmp36 = Ftmp20*y;
double Ftmp37 = Ftmp11*z;
double Ftmp38 = Ftmp13*z;
double Ftmp39 = Ftmp14*z;
double Ftmp40 = Ftmp0*Ftmp17;
double Ftmp41 = Ftmp0*Ftmp19;
double Ftmp42 = Ftmp1*Ftmp11;
double Ftmp43 = Ftmp1*Ftmp13;
double Ftmp44 = Ftmp2*Ftmp5;
double Ftmp45 = Ftmp2*Ftmp7;
double Ftmp46 = (1.0/4.0)*Ftmp4;
double Ftmp47 = Ftmp10*Ftmp46;
double Ftmp48 = (1.0/12.0)*Ftmp4;
double Ftmp49 = Ftmp12*Ftmp48;
double Ftmp50 = Ftmp16*Ftmp46;
double Ftmp51 = Ftmp18*Ftmp48;
double Ftmp52 = (1.0/12.0)*Ftmp6;
double Ftmp53 = Ftmp10*Ftmp52;
double Ftmp54 = Ftmp16*Ftmp52;
double Ftmp55 = (1.0/4.0)*Ftmp10*Ftmp16;
double Ftmp56 = (1.0/12.0)*Ftmp10*Ftmp18;
double Ftmp57 = (1.0/12.0)*Ftmp12*Ftmp16;
double Ftmp58 = Ftmp55*x;
double Ftmp59 = Ftmp50*y;
double Ftmp60 = Ftmp47*z;
#pragma omp atomic
F[0] += -Ftmp0*L[11] - Ftmp1*L[12] - Ftmp11*L[13] - Ftmp13*L[26] - Ftmp14*L[45] - Ftmp15*L[71] - Ftmp17*L[15] - Ftmp19*L[29] - Ftmp2*L[14] - Ftmp20*L[49] - Ftmp21*L[76] - Ftmp22*L[23] - Ftmp23*L[41] - Ftmp24*L[66] - Ftmp25*L[25] - Ftmp26*L[44] - Ftmp27*L[70] - Ftmp28*L[21] - Ftmp29*L[22] - Ftmp3*L[24] - Ftmp30*L[36] - Ftmp31*L[37] - Ftmp32*L[57] - Ftmp33*L[58] - Ftmp34*L[28] - Ftmp35*L[48] - Ftmp36*L[75] - Ftmp37*L[27] - Ftmp38*L[46] - Ftmp39*L[72] - Ftmp40*L[43] - Ftmp41*L[69] - Ftmp42*L[42] - Ftmp43*L[67] - Ftmp44*L[39] - Ftmp45*L[60] - Ftmp47*L[38] - Ftmp49*L[62] - Ftmp5*L[10] - Ftmp50*L[40] - Ftmp51*L[65] - Ftmp53*L[59] - Ftmp54*L[61] - Ftmp55*L[47] - Ftmp56*L[74] - Ftmp57*L[73] - Ftmp58*L[68] - Ftmp59*L[64] - Ftmp60*L[63] - Ftmp7*L[20] - Ftmp8*L[35] - Ftmp9*L[56] - x*L[4] - y*L[5] - z*L[6] - L[1];
#pragma omp atomic
F[1] += -Ftmp0*L[13] - Ftmp1*L[14] - Ftmp11*L[16] - Ftmp13*L[30] - Ftmp14*L[50] - Ftmp15*L[77] - Ftmp17*L[18] - Ftmp19*L[33] - Ftmp2*L[17] - Ftmp20*L[54] - Ftmp21*L[82] - Ftmp22*L[26] - Ftmp23*L[45] - Ftmp24*L[71] - Ftmp25*L[28] - Ftmp26*L[48] - Ftmp27*L[75] - Ftmp28*L[23] - Ftmp29*L[24] - Ftmp3*L[27] - Ftmp30*L[38] - Ftmp31*L[39] - Ftmp32*L[59] - Ftmp33*L[60] - Ftmp34*L[32] - Ftmp35*L[53] - Ftmp36*L[81] - Ftmp37*L[31] - Ftmp38*L[51] - Ftmp39*L[78] - Ftmp40*L[47] - Ftmp41*L[74] - Ftmp42*L[46] - Ftmp43*L[72] - Ftmp44*L[42] - Ftmp45*L[63] - Ftmp47*L[41] - Ftmp49*L[66] - Ftmp5*L[11] - Ftmp50*L[43] - Ftmp51*L[69] - Ftmp53*L[62] - Ftmp54*L[64] - Ftmp55*L[52] - Ftmp56*L[80] - Ftmp57*L[79] - Ftmp58*L[73] - Ftmp59*L[68] - Ftmp60*L[67] - Ftmp7*L[21] - Ftmp8*L[36] - Ftmp9*L[57] - x*L[5] - y*L[7] - z*L[8] - L[2];
#pragma omp atomic
F[2] += -Ftmp0*L[14] - Ftmp1*L[15] - Ftmp11*L[17] - Ftmp13*L[31] - Ftmp14*L[51] - Ftmp15*L[78] - Ftmp17*L[19] - Ftmp19*L[34] - Ftmp2*L[18] - Ftmp20*L[55] - Ftmp21*L[83] - Ftmp22*L[27] - Ftmp23*L[46] - Ftmp24*L[72] - Ftmp25*L[29] - Ftmp26*L[49] - Ftmp27*L[76] - Ftmp28*L[24] - Ftmp29*L[25] - Ftmp3*L[28] - Ftmp30*L[39] - Ftmp31*L[40] - Ftmp32*L[60] - Ftmp33*L[61] - Ftmp34*L[33] - Ftmp35*L[54] - Ftmp36*L[82] - Ftmp37*L[32] - Ftmp38*L[52] - Ftmp39*L[79] - Ftmp40*L[48] - Ftmp41*L[75] - Ftmp42*L[47] - Ftmp43*L[73] - Ftmp44*L[43] - Ftmp45*L[64] - Ftmp47*L[42] - Ftmp49*L[67] - Ftmp5*L[12] - Ftmp50*L[44] - Ftmp51*L[70] - Ftmp53*L[63] - Ftmp54*L[65] - Ftmp55*L[53] - Ftmp56*L[81] - Ftmp57*L[80] - Ftmp58*L[74] - Ftmp59*L[69] - Ftmp60*L[68] - Ftmp7*L[22] - Ftmp8*L[37] - Ftmp9*L[58] - x*L[6] - y*L[8] - z*L[9] - L[3];

}

void field_m0_M2P_6(double x, double y, double z, double * M, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = 1.0*pow(R, -3.0);
double Ftmp1 = pow(R, -5.0);
double Ftmp2 = 3.0*Ftmp1;
double Ftmp3 = Ftmp2*M[5];
double Ftmp4 = Ftmp2*z;
double Ftmp5 = Ftmp0*M[0];
double Ftmp6 = y*z;
double Ftmp7 = pow(R, -7.0);
double Ftmp8 = 15.0*Ftmp7;
double Ftmp9 = Ftmp8*M[14];
double Ftmp10 = x*y;
double Ftmp11 = Ftmp2*M[2];
double Ftmp12 = Ftmp4*M[3];
double Ftmp13 = (x*x);
double Ftmp14 = Ftmp2*M[1];
double Ftmp15 = Ftmp6*x;
double Ftmp16 = Ftmp15*Ftmp8;
double Ftmp17 = Ftmp13*Ftmp8;
double Ftmp18 = pow(R, -9.0);
double Ftmp19 = 105.0*Ftmp18;
double Ftmp20 = Ftmp13*Ftmp19;
double Ftmp21 = -9.0*Ftmp1;
double Ftmp22 = Ftmp17 + Ftmp21;
double Ftmp23 = -Ftmp2;
double Ftmp24 = (y*y);
double Ftmp25 = Ftmp24*Ftmp8;
double Ftmp26 = Ftmp23 + Ftmp25;
double Ftmp27 = (z*z);
double Ftmp28 = Ftmp27*Ftmp8;
double Ftmp29 = Ftmp23 + Ftmp28;
double Ftmp30 = Ftmp26*M[7];
double Ftmp31 = Ftmp29*M[9];
double Ftmp32 = 45.0*Ftmp7;
double Ftmp33 = -Ftmp32;
double Ftmp34 = Ftmp20 + Ftmp33;
double Ftmp35 = Ftmp34*M[21];
double Ftmp36 = 1.0*y;
double Ftmp37 = Ftmp19*Ftmp24;
double Ftmp38 = Ftmp33 + Ftmp37;
double Ftmp39 = Ftmp38*M[26];
double Ftmp40 = 3.0*y;
double Ftmp41 = 35.0*Ftmp18;
double Ftmp42 = (Ftmp27*Ftmp41 - 5.0*Ftmp7)*M[28];
double Ftmp43 = Ftmp34*M[22];
double Ftmp44 = 1.0*z;
double Ftmp45 = -Ftmp8;
double Ftmp46 = Ftmp37 + Ftmp45;
double Ftmp47 = Ftmp46*M[27];
double Ftmp48 = Ftmp19*Ftmp27;
double Ftmp49 = Ftmp33 + Ftmp48;
double Ftmp50 = Ftmp44*Ftmp49;
double Ftmp51 = 315.0*Ftmp18;
double Ftmp52 = -Ftmp51;
double Ftmp53 = pow(R, -11.0);
double Ftmp54 = 945.0*Ftmp53;
double Ftmp55 = Ftmp13*Ftmp54;
double Ftmp56 = Ftmp52 + Ftmp55;
double Ftmp57 = Ftmp56*M[39];
double Ftmp58 = Ftmp10*Ftmp34;
double Ftmp59 = Ftmp38*M[16];
double Ftmp60 = Ftmp45 + Ftmp48;
double Ftmp61 = Ftmp60*M[18];
double Ftmp62 = Ftmp36*x;
double Ftmp63 = x*z;
double Ftmp64 = Ftmp34*Ftmp63;
double Ftmp65 = Ftmp46*M[17];
double Ftmp66 = Ftmp49*M[19];
double Ftmp67 = Ftmp24*Ftmp54;
double Ftmp68 = Ftmp52 + Ftmp67;
double Ftmp69 = Ftmp36*Ftmp68;
double Ftmp70 = Ftmp69*M[46];
double Ftmp71 = -Ftmp19;
double Ftmp72 = Ftmp27*Ftmp53;
double Ftmp73 = 315.0*Ftmp72;
double Ftmp74 = Ftmp71 + Ftmp73;
double Ftmp75 = Ftmp74*M[48];
double Ftmp76 = Ftmp40*Ftmp75;
double Ftmp77 = Ftmp68*M[31];
double Ftmp78 = -75.0*Ftmp7;
double Ftmp79 = 1.0*Ftmp13;
double Ftmp80 = Ftmp46*M[13];
double Ftmp81 = Ftmp60*M[15];
double Ftmp82 = 525.0*Ftmp18;
double Ftmp83 = -Ftmp82;
double Ftmp84 = Ftmp55 + Ftmp83;
double Ftmp85 = Ftmp13*y;
double Ftmp86 = Ftmp13*z;
double Ftmp87 = Ftmp27*Ftmp54;
double Ftmp88 = Ftmp52 + Ftmp87;
double Ftmp89 = Ftmp36*M[33];
double Ftmp90 = Ftmp69*M[26];
double Ftmp91 = Ftmp13*Ftmp40;
double Ftmp92 = (-Ftmp41 + Ftmp73)*M[28];
double Ftmp93 = Ftmp13*Ftmp44;
double Ftmp94 = (Ftmp67 + Ftmp71)*M[27];
double Ftmp95 = Ftmp88*M[29];
double Ftmp96 = 4725.0*Ftmp53;
double Ftmp97 = -Ftmp96;
double Ftmp98 = pow(R, -13.0);
double Ftmp99 = 10395.0*Ftmp98;
double Ftmp100 = Ftmp13*Ftmp99;
double Ftmp101 = 2835.0*Ftmp53;
double Ftmp102 = -Ftmp101;
double Ftmp103 = Ftmp24*Ftmp99;
double Ftmp104 = Ftmp36*(Ftmp102 + Ftmp103)*M[46];
double Ftmp105 = 3465.0*Ftmp98;
double Ftmp106 = Ftmp105*Ftmp27;
double Ftmp107 = (Ftmp106 - Ftmp54)*M[48];
double Ftmp108 = 225.0*Ftmp7;
double Ftmp109 = (x*x*x*x);
double Ftmp110 = Ftmp109*Ftmp54;
double Ftmp111 = 1050.0*Ftmp18;
double Ftmp112 = Ftmp108 + Ftmp110 - Ftmp111*Ftmp13;
double Ftmp113 = (y*y*y*y);
double Ftmp114 = Ftmp113*Ftmp54;
double Ftmp115 = 630.0*Ftmp18;
double Ftmp116 = Ftmp114 - Ftmp115*Ftmp24 + Ftmp32;
double Ftmp117 = (z*z*z*z);
double Ftmp118 = Ftmp117*Ftmp54;
double Ftmp119 = -Ftmp115*Ftmp27 + Ftmp118 + Ftmp32;
double Ftmp120 = Ftmp116*M[30];
double Ftmp121 = Ftmp119*M[34];
double Ftmp122 = 1575.0*Ftmp18;
double Ftmp123 = Ftmp109*Ftmp98;
double Ftmp124 = 10395.0*Ftmp123;
double Ftmp125 = Ftmp13*Ftmp53;
double Ftmp126 = 9450.0*Ftmp125;
double Ftmp127 = Ftmp122 + Ftmp124 - Ftmp126;
double Ftmp128 = Ftmp127*M[57];
double Ftmp129 = Ftmp113*Ftmp99;
double Ftmp130 = 9450.0*Ftmp53;
double Ftmp131 = Ftmp130*Ftmp24;
double Ftmp132 = Ftmp122 + Ftmp129 - Ftmp131;
double Ftmp133 = Ftmp132*M[71];
double Ftmp134 = (Ftmp105*Ftmp117 + Ftmp19 - 1890.0*Ftmp72)*M[75];
double Ftmp135 = Ftmp127*M[58];
double Ftmp136 = 5670.0*Ftmp53;
double Ftmp137 = Ftmp136*Ftmp24;
double Ftmp138 = Ftmp129 - Ftmp137 + Ftmp51;
double Ftmp139 = Ftmp138*M[72];
double Ftmp140 = Ftmp117*Ftmp99;
double Ftmp141 = Ftmp130*Ftmp27;
double Ftmp142 = Ftmp122 + Ftmp140 - Ftmp141;
double Ftmp143 = Ftmp142*Ftmp44;
double Ftmp144 = Ftmp10*Ftmp127;
double Ftmp145 = Ftmp132*M[50];
double Ftmp146 = Ftmp136*Ftmp27;
double Ftmp147 = Ftmp140 - Ftmp146 + Ftmp51;
double Ftmp148 = Ftmp147*M[54];
double Ftmp149 = Ftmp127*Ftmp63;
double Ftmp150 = Ftmp138*M[51];
double Ftmp151 = Ftmp142*M[55];
double Ftmp152 = 14175.0*Ftmp53;
double Ftmp153 = pow(R, -15.0);
double Ftmp154 = 135135.0*Ftmp153;
double Ftmp155 = Ftmp109*Ftmp154;
double Ftmp156 = Ftmp13*Ftmp98;
double Ftmp157 = 103950.0*Ftmp156;
double Ftmp158 = Ftmp152 + Ftmp155 - Ftmp157;
double Ftmp159 = Ftmp113*Ftmp154;
double Ftmp160 = Ftmp24*Ftmp98;
double Ftmp161 = 103950.0*Ftmp160;
double Ftmp162 = Ftmp152 + Ftmp159 - Ftmp161;
double Ftmp163 = Ftmp162*M[78];
double Ftmp164 = 3675.0*Ftmp18;
double Ftmp165 = Ftmp138*M[45];
double Ftmp166 = Ftmp147*M[49];
double Ftmp167 = 33075.0*Ftmp53;
double Ftmp168 = Ftmp155 - 145530.0*Ftmp156 + Ftmp167;
double Ftmp169 = Ftmp117*Ftmp154;
double Ftmp170 = Ftmp27*Ftmp98;
double Ftmp171 = Ftmp152 + Ftmp169 - 103950.0*Ftmp170;
double Ftmp172 = Ftmp36*M[82];
double Ftmp173 = Ftmp13*Ftmp36;
double Ftmp174 = Ftmp162*M[71];
double Ftmp175 = 45045.0*Ftmp117*Ftmp153;
double Ftmp176 = (-20790.0*Ftmp170 + Ftmp175 + Ftmp54)*M[75];
double Ftmp177 = 62370.0*Ftmp160;
double Ftmp178 = (Ftmp101 + Ftmp159 - Ftmp177)*M[72];
double Ftmp179 = Ftmp171*M[76];
double Ftmp180 = -11025.0*Ftmp18;
double Ftmp181 = Ftmp154*(x*x*x*x*x*x);
double Ftmp182 = -Ftmp122;
double Ftmp183 = Ftmp154*(y*y*y*y*y*y);
double Ftmp184 = 155925.0*Ftmp98;
double Ftmp185 = 42525.0*Ftmp53;
double Ftmp186 = (-Ftmp113*Ftmp184 + Ftmp182 + Ftmp183 + Ftmp185*Ftmp24)*M[77];
double Ftmp187 = Ftmp154*(z*z*z*z*z*z);
double Ftmp188 = (-Ftmp117*Ftmp184 + Ftmp182 + Ftmp185*Ftmp27 + Ftmp187)*M[83];
double Ftmp189 = -Ftmp24*Ftmp51;
double Ftmp190 = Ftmp24*Ftmp55;
double Ftmp191 = -Ftmp20;
double Ftmp192 = Ftmp191 + Ftmp32;
double Ftmp193 = Ftmp189 + Ftmp190 + Ftmp192;
double Ftmp194 = -Ftmp27*Ftmp51;
double Ftmp195 = Ftmp27*Ftmp55;
double Ftmp196 = Ftmp192 + Ftmp194 + Ftmp195;
double Ftmp197 = -Ftmp48;
double Ftmp198 = Ftmp197 + Ftmp8;
double Ftmp199 = -Ftmp37;
double Ftmp200 = Ftmp27*Ftmp67;
double Ftmp201 = Ftmp199 + Ftmp200;
double Ftmp202 = Ftmp198 + Ftmp201;
double Ftmp203 = Ftmp101*Ftmp24;
double Ftmp204 = -Ftmp203;
double Ftmp205 = Ftmp100*Ftmp24;
double Ftmp206 = Ftmp204 + Ftmp205;
double Ftmp207 = 945.0*Ftmp18;
double Ftmp208 = Ftmp101*Ftmp13;
double Ftmp209 = -Ftmp208;
double Ftmp210 = Ftmp207 + Ftmp209;
double Ftmp211 = Ftmp206 + Ftmp210;
double Ftmp212 = Ftmp211*M[62];
double Ftmp213 = -Ftmp55;
double Ftmp214 = Ftmp213 + Ftmp51;
double Ftmp215 = Ftmp101*Ftmp27;
double Ftmp216 = -Ftmp215;
double Ftmp217 = Ftmp100*Ftmp27;
double Ftmp218 = Ftmp216 + Ftmp217;
double Ftmp219 = Ftmp214 + Ftmp218;
double Ftmp220 = Ftmp219*M[64];
double Ftmp221 = -Ftmp67;
double Ftmp222 = Ftmp103*Ftmp27;
double Ftmp223 = Ftmp222 + Ftmp51;
double Ftmp224 = Ftmp216 + Ftmp221 + Ftmp223;
double Ftmp225 = Ftmp224*M[73];
double Ftmp226 = Ftmp206 + Ftmp214;
double Ftmp227 = Ftmp226*M[63];
double Ftmp228 = Ftmp210 + Ftmp218;
double Ftmp229 = Ftmp228*M[65];
double Ftmp230 = -Ftmp87;
double Ftmp231 = Ftmp204 + Ftmp223 + Ftmp230;
double Ftmp232 = Ftmp231*M[74];
double Ftmp233 = Ftmp10*Ftmp211;
double Ftmp234 = Ftmp10*Ftmp219;
double Ftmp235 = Ftmp226*Ftmp63;
double Ftmp236 = Ftmp228*Ftmp63;
double Ftmp237 = 31185.0*Ftmp98;
double Ftmp238 = Ftmp13*Ftmp237;
double Ftmp239 = -Ftmp238;
double Ftmp240 = 8505.0*Ftmp53;
double Ftmp241 = Ftmp239 + Ftmp240;
double Ftmp242 = Ftmp13*Ftmp154;
double Ftmp243 = Ftmp24*Ftmp242;
double Ftmp244 = -Ftmp237*Ftmp24;
double Ftmp245 = Ftmp243 + Ftmp244;
double Ftmp246 = Ftmp15*(Ftmp241 + Ftmp245);
double Ftmp247 = Ftmp242*Ftmp27;
double Ftmp248 = Ftmp237*Ftmp27;
double Ftmp249 = -Ftmp248;
double Ftmp250 = Ftmp247 + Ftmp249;
double Ftmp251 = Ftmp15*(Ftmp241 + Ftmp250);
double Ftmp252 = Ftmp154*Ftmp24*Ftmp27;
double Ftmp253 = Ftmp244 + Ftmp252;
double Ftmp254 = Ftmp240 + Ftmp249 + Ftmp253;
double Ftmp255 = -Ftmp24*Ftmp96;
double Ftmp256 = Ftmp213 + Ftmp82;
double Ftmp257 = -Ftmp27*Ftmp96;
double Ftmp258 = Ftmp19 + Ftmp230;
double Ftmp259 = Ftmp221 + Ftmp222;
double Ftmp260 = 51975.0*Ftmp98;
double Ftmp261 = -Ftmp24*Ftmp260;
double Ftmp262 = Ftmp243 + Ftmp261;
double Ftmp263 = Ftmp152 + Ftmp239;
double Ftmp264 = -Ftmp100;
double Ftmp265 = Ftmp264 + Ftmp96;
double Ftmp266 = -Ftmp260*Ftmp27;
double Ftmp267 = Ftmp247 + Ftmp266;
double Ftmp268 = -Ftmp103;
double Ftmp269 = Ftmp101 + Ftmp252;
double Ftmp270 = -Ftmp27*Ftmp99;
double Ftmp271 = -Ftmp207;
double Ftmp272 = Ftmp208 + Ftmp271;
double Ftmp273 = Ftmp13*Ftmp159;
double Ftmp274 = 62370.0*Ftmp156;
double Ftmp275 = -Ftmp24*Ftmp274;
double Ftmp276 = Ftmp273 + Ftmp275;
double Ftmp277 = 17010.0*Ftmp53;
double Ftmp278 = -Ftmp113*Ftmp237 + Ftmp24*Ftmp277;
double Ftmp279 = Ftmp13*Ftmp169;
double Ftmp280 = -Ftmp27*Ftmp274;
double Ftmp281 = Ftmp279 + Ftmp280;
double Ftmp282 = -Ftmp117*Ftmp237 + Ftmp27*Ftmp277;
double Ftmp283 = Ftmp152*Ftmp24;
double Ftmp284 = -Ftmp157*Ftmp24 + Ftmp182;
double Ftmp285 = -Ftmp124;
double Ftmp286 = Ftmp155*Ftmp24;
double Ftmp287 = Ftmp285 + Ftmp286;
double Ftmp288 = -Ftmp157*Ftmp27;
double Ftmp289 = Ftmp155*Ftmp27;
double Ftmp290 = Ftmp285 + Ftmp289;
double Ftmp291 = Ftmp152*Ftmp27 + Ftmp182;
double Ftmp292 = -Ftmp177*Ftmp27;
double Ftmp293 = Ftmp292 + Ftmp52;
double Ftmp294 = -Ftmp140;
double Ftmp295 = Ftmp169*Ftmp24;
double Ftmp296 = Ftmp294 + Ftmp295;
double Ftmp297 = -Ftmp129;
double Ftmp298 = Ftmp159*Ftmp27;
double Ftmp299 = Ftmp297 + Ftmp298;
double Ftmp300 = Ftmp24*Ftmp247;
double Ftmp301 = -Ftmp205 + Ftmp215 + Ftmp300;
double Ftmp302 = Ftmp203 - Ftmp217;
double Ftmp303 = x*M[6];
double Ftmp304 = Ftmp17 + Ftmp23;
double Ftmp305 = Ftmp21 + Ftmp25;
double Ftmp306 = Ftmp304*M[4];
double Ftmp307 = 1.0*x;
double Ftmp308 = 3.0*x;
double Ftmp309 = Ftmp20 + Ftmp45;
double Ftmp310 = Ftmp309*M[24];
double Ftmp311 = Ftmp38*M[31];
double Ftmp312 = Ftmp44*x;
double Ftmp313 = 3.0*Ftmp63;
double Ftmp314 = Ftmp309*M[12];
double Ftmp315 = Ftmp56*M[22];
double Ftmp316 = Ftmp309*M[11];
double Ftmp317 = 1.0*Ftmp24;
double Ftmp318 = Ftmp24*x;
double Ftmp319 = Ftmp56*M[21];
double Ftmp320 = Ftmp24*z;
double Ftmp321 = (Ftmp55 + Ftmp71)*M[24];
double Ftmp322 = Ftmp67 + Ftmp83;
double Ftmp323 = Ftmp36*Ftmp63;
double Ftmp324 = Ftmp24*Ftmp307;
double Ftmp325 = Ftmp24*Ftmp308;
double Ftmp326 = Ftmp24*Ftmp44;
double Ftmp327 = (Ftmp100 + Ftmp102)*M[39];
double Ftmp328 = Ftmp110 - Ftmp115*Ftmp13 + Ftmp32;
double Ftmp329 = Ftmp108 - Ftmp111*Ftmp24 + Ftmp114;
double Ftmp330 = Ftmp328*M[20];
double Ftmp331 = 5670.0*Ftmp125;
double Ftmp332 = Ftmp124 - Ftmp331 + Ftmp51;
double Ftmp333 = Ftmp332*M[60];
double Ftmp334 = Ftmp132*M[78];
double Ftmp335 = Ftmp332*M[37];
double Ftmp336 = Ftmp158*M[58];
double Ftmp337 = Ftmp332*M[36];
double Ftmp338 = Ftmp24*Ftmp53;
double Ftmp339 = Ftmp158*M[57];
double Ftmp340 = (Ftmp101 + Ftmp155 - Ftmp274)*M[60];
double Ftmp341 = Ftmp159 - 145530.0*Ftmp160 + Ftmp167;
double Ftmp342 = (-Ftmp109*Ftmp184 + Ftmp13*Ftmp185 + Ftmp181 + Ftmp182)*M[56];
double Ftmp343 = 218295.0*Ftmp98;
double Ftmp344 = Ftmp190 + Ftmp199;
double Ftmp345 = -Ftmp13*Ftmp51 + Ftmp32;
double Ftmp346 = Ftmp344 + Ftmp345;
double Ftmp347 = Ftmp191 + Ftmp195 + Ftmp198;
double Ftmp348 = Ftmp194 + Ftmp201 + Ftmp32;
double Ftmp349 = Ftmp209 + Ftmp51;
double Ftmp350 = Ftmp205 + Ftmp221;
double Ftmp351 = Ftmp349 + Ftmp350;
double Ftmp352 = Ftmp351*M[67];
double Ftmp353 = Ftmp217 + Ftmp230;
double Ftmp354 = Ftmp349 + Ftmp353;
double Ftmp355 = Ftmp354*M[69];
double Ftmp356 = Ftmp204 + Ftmp207 + Ftmp216 + Ftmp222;
double Ftmp357 = Ftmp356*M[80];
double Ftmp358 = Ftmp351*Ftmp6;
double Ftmp359 = Ftmp354*Ftmp6;
double Ftmp360 = Ftmp356*Ftmp6;
double Ftmp361 = -Ftmp13*Ftmp96 + Ftmp82;
double Ftmp362 = -51975.0*Ftmp156;
double Ftmp363 = Ftmp152 + Ftmp362;
double Ftmp364 = Ftmp101 + Ftmp264;
double Ftmp365 = Ftmp268 + Ftmp96;
double Ftmp366 = Ftmp101 + Ftmp239;
double Ftmp367 = Ftmp247 + Ftmp270;
double Ftmp368 = Ftmp254*Ftmp323;
double Ftmp369 = Ftmp13*Ftmp152;
double Ftmp370 = Ftmp208 + Ftmp52;
double Ftmp371 = Ftmp203 + Ftmp271;
double Ftmp372 = -31185.0*Ftmp123 + 17010.0*Ftmp125;
double Ftmp373 = Ftmp331 + Ftmp52;
double Ftmp374 = Ftmp215 + Ftmp280;
double Ftmp375 = -Ftmp161*Ftmp27;
double Ftmp376 = Ftmp208 - Ftmp222;
double Ftmp377 = y*M[8];
double Ftmp378 = Ftmp21 + Ftmp28;
double Ftmp379 = Ftmp307*M[29];
double Ftmp380 = Ftmp36*z;
double Ftmp381 = Ftmp27*x;
double Ftmp382 = Ftmp27*y;
double Ftmp383 = Ftmp40*Ftmp63;
double Ftmp384 = Ftmp27*Ftmp307;
double Ftmp385 = Ftmp27*(Ftmp83 + Ftmp87);
double Ftmp386 = Ftmp108 - Ftmp111*Ftmp27 + Ftmp118;
double Ftmp387 = Ftmp307*M[76];
double Ftmp388 = Ftmp27*(Ftmp167 + Ftmp169 - 145530.0*Ftmp170);
double Ftmp389 = Ftmp191 + Ftmp344 + Ftmp8;
double Ftmp390 = Ftmp195 + Ftmp197 + Ftmp345;
double Ftmp391 = Ftmp189 + Ftmp197 + Ftmp200 + Ftmp32;
double Ftmp392 = Ftmp252 + Ftmp261;
double Ftmp393 = Ftmp141 + Ftmp182;
#pragma omp atomic
F[0] += Ftmp0*M[1] - Ftmp10*Ftmp11 - Ftmp10*Ftmp145 - Ftmp10*Ftmp224*M[52] - Ftmp10*Ftmp59 - Ftmp104*Ftmp86 - Ftmp107*Ftmp40*Ftmp86 + Ftmp112*x*M[20] + Ftmp112*M[35] + Ftmp116*M[45] + Ftmp119*M[49] - Ftmp12*x + Ftmp120*x + Ftmp121*x - Ftmp128*y - Ftmp13*Ftmp14 - Ftmp13*Ftmp6*(Ftmp100 + Ftmp97)*M[39] + Ftmp13*Ftmp90 - Ftmp13*(Ftmp20 + Ftmp78)*M[10] - Ftmp13*(Ftmp124 - 13230.0*Ftmp125 + Ftmp164)*M[35] - Ftmp13*(Ftmp205 + Ftmp255 + Ftmp256)*M[38] - Ftmp13*(Ftmp217 + Ftmp256 + Ftmp257)*M[40] - Ftmp133*Ftmp36 - Ftmp134*Ftmp40 - Ftmp135*z - Ftmp139*Ftmp44 - Ftmp143*M[76] - Ftmp144*M[36] - Ftmp148*Ftmp62 - Ftmp149*M[37] + Ftmp15*Ftmp158*M[60] + Ftmp15*Ftmp163 + Ftmp15*Ftmp254*M[80] + Ftmp15*Ftmp56*M[24] + Ftmp15*Ftmp77 - Ftmp150*Ftmp63 - Ftmp151*Ftmp63 + Ftmp16*M[8] - Ftmp165*Ftmp79 - Ftmp166*Ftmp79 + Ftmp168*Ftmp85*M[57] + Ftmp168*Ftmp86*M[58] + Ftmp17*y*M[5] + Ftmp17*z*M[6] + Ftmp171*Ftmp172*Ftmp63 + Ftmp173*Ftmp174 + Ftmp173*(Ftmp249 + Ftmp268 + Ftmp269)*M[73] + Ftmp176*Ftmp91 + Ftmp178*Ftmp93 + Ftmp179*Ftmp93 + Ftmp186*x + Ftmp188*x + Ftmp193*x*M[23] + Ftmp193*M[38] + Ftmp196*x*M[25] + Ftmp196*M[40] - Ftmp20*Ftmp6*M[14] + Ftmp202*x*M[32] + Ftmp202*M[47] - Ftmp212*y + Ftmp22*x*M[4] + Ftmp22*M[10] - Ftmp220*y - Ftmp225*Ftmp36 - Ftmp227*z - Ftmp229*z - Ftmp231*Ftmp63*M[53] - Ftmp232*Ftmp44 - Ftmp233*M[41] - Ftmp234*M[43] - Ftmp235*M[42] - Ftmp236*M[44] + Ftmp246*M[67] + Ftmp251*M[69] + Ftmp26*M[13] + Ftmp29*M[15] - Ftmp3*y + Ftmp30*x + Ftmp31*x - Ftmp35*y - Ftmp36*Ftmp39 - Ftmp4*M[6] - Ftmp40*Ftmp42 - Ftmp43*z - Ftmp44*Ftmp47 + Ftmp5*x - Ftmp50*M[29] + Ftmp57*Ftmp6 - Ftmp58*M[11] + Ftmp6*Ftmp9 - Ftmp61*Ftmp62 - Ftmp63*Ftmp65 - Ftmp63*Ftmp66 + Ftmp63*Ftmp88*Ftmp89 - Ftmp64*M[12] + Ftmp70*z + Ftmp76*z - Ftmp79*Ftmp80 - Ftmp79*Ftmp81 - Ftmp79*(Ftmp258 + Ftmp259)*M[47] + Ftmp84*Ftmp85*M[21] + Ftmp84*Ftmp86*M[22] + Ftmp85*(Ftmp262 + Ftmp263)*M[62] + Ftmp85*(Ftmp265 + Ftmp267)*M[64] + Ftmp86*(Ftmp262 + Ftmp265)*M[63] + Ftmp86*(Ftmp263 + Ftmp267)*M[65] + Ftmp91*Ftmp92 + Ftmp93*Ftmp94 + Ftmp93*Ftmp95 + Ftmp93*(Ftmp244 + Ftmp269 + Ftmp270)*M[74] + x*(Ftmp272 + Ftmp276 + Ftmp278)*M[66] + x*(Ftmp272 + Ftmp281 + Ftmp282)*M[70] + x*(-218295.0*Ftmp123 + 99225.0*Ftmp125 + Ftmp180 + Ftmp181)*M[56] + x*(Ftmp126 + Ftmp283 + Ftmp284 + Ftmp287)*M[59] + x*(Ftmp126 + Ftmp288 + Ftmp290 + Ftmp291)*M[61] + x*(Ftmp137 + Ftmp215 + Ftmp293 + Ftmp299)*M[79] + x*(Ftmp146 + Ftmp203 + Ftmp293 + Ftmp296)*M[81] + x*(-Ftmp24*Ftmp248 + Ftmp301 + Ftmp302 + Ftmp56)*M[68];
#pragma omp atomic
F[1] += Ftmp0*M[2] - Ftmp10*Ftmp14 - Ftmp107*Ftmp24*Ftmp313 - Ftmp11*Ftmp24 + Ftmp119*M[54] - Ftmp12*y + Ftmp121*y - Ftmp128*x - Ftmp132*Ftmp6*M[51] - Ftmp132*Ftmp62*M[45] - Ftmp133*Ftmp307 - Ftmp134*Ftmp308 - Ftmp143*M[82] - Ftmp144*M[35] - Ftmp148*Ftmp317 + Ftmp15*Ftmp315 + Ftmp15*Ftmp336 - Ftmp151*Ftmp6 + Ftmp162*Ftmp323*M[72] - Ftmp166*Ftmp62 + Ftmp171*Ftmp326*M[82] + Ftmp176*Ftmp325 + Ftmp179*Ftmp323 + Ftmp188*y - Ftmp212*x - Ftmp220*x - Ftmp224*Ftmp62*M[47] - Ftmp225*Ftmp307 - Ftmp233*M[38] - Ftmp234*M[40] - Ftmp24*Ftmp316 - Ftmp24*Ftmp327*Ftmp63 - Ftmp24*Ftmp337 - Ftmp24*(Ftmp350 + Ftmp361)*M[41] - Ftmp24*(Ftmp37 + Ftmp78)*M[16] - Ftmp24*(Ftmp129 + Ftmp164 - 13230.0*Ftmp338)*M[50] - Ftmp24*(Ftmp213 + Ftmp217 + Ftmp258)*M[43] - Ftmp24*(Ftmp257 + Ftmp259 + Ftmp82)*M[52] + Ftmp246*M[63] + Ftmp25*x*M[5] + Ftmp25*z*M[8] + Ftmp251*M[65] + Ftmp29*M[18] - Ftmp3*x + Ftmp303*Ftmp6*Ftmp8 + Ftmp304*M[11] + Ftmp305*y*M[7] + Ftmp305*M[16] + Ftmp306*y - Ftmp307*Ftmp39 - Ftmp308*Ftmp42 + Ftmp31*y - Ftmp310*z - Ftmp311*z + Ftmp312*Ftmp68*M[46] + Ftmp313*Ftmp75 - Ftmp314*Ftmp6 - Ftmp317*Ftmp61 + Ftmp318*Ftmp319 + Ftmp318*Ftmp339 - Ftmp318*Ftmp44*(Ftmp103 + Ftmp97)*M[46] + Ftmp318*(Ftmp245 + Ftmp363)*M[62] + Ftmp318*(Ftmp250 + Ftmp364)*M[64] + Ftmp320*Ftmp321 + Ftmp320*Ftmp322*M[31] + Ftmp320*Ftmp340 + Ftmp320*Ftmp341*M[78] + Ftmp320*(Ftmp366 + Ftmp367)*M[69] + Ftmp320*(Ftmp152 + Ftmp253 + Ftmp266)*M[80] + Ftmp320*(Ftmp243 + Ftmp362 + Ftmp365)*M[67] + Ftmp322*Ftmp324*M[26] + Ftmp323*Ftmp95 + Ftmp324*Ftmp341*M[71] + Ftmp324*(Ftmp252 + Ftmp266 + Ftmp365)*M[73] + Ftmp325*Ftmp92 + Ftmp326*Ftmp88*M[33] + Ftmp328*M[36] + Ftmp329*y*M[30] + Ftmp329*M[50] + Ftmp330*y - Ftmp333*z - Ftmp334*z - Ftmp335*Ftmp6 + Ftmp342*y + Ftmp346*y*M[23] + Ftmp346*M[41] + Ftmp347*y*M[25] + Ftmp347*M[43] + Ftmp348*y*M[32] + Ftmp348*M[52] - Ftmp35*x - Ftmp352*z - Ftmp355*z - Ftmp357*z - Ftmp358*M[42] - Ftmp359*M[44] - Ftmp360*M[53] + Ftmp368*M[74] - Ftmp37*Ftmp63*M[14] - Ftmp38*Ftmp6*M[17] - Ftmp38*Ftmp62*M[13] - Ftmp4*M[8] + Ftmp5*y - Ftmp50*M[33] + Ftmp57*Ftmp63 - Ftmp58*M[10] - Ftmp6*Ftmp66 - Ftmp62*Ftmp81 + Ftmp63*Ftmp69*M[27] + Ftmp63*Ftmp9 + y*(Ftmp290 + Ftmp373 + Ftmp374)*M[61] + y*(Ftmp131 + Ftmp291 + Ftmp299 + Ftmp375)*M[79] + y*(Ftmp146 + Ftmp281 + Ftmp294 + Ftmp370)*M[70] + y*(Ftmp275 + Ftmp286 + Ftmp371 + Ftmp372)*M[59] + y*(Ftmp282 + Ftmp292 + Ftmp295 + Ftmp371)*M[81] + y*(-Ftmp113*Ftmp343 + Ftmp180 + Ftmp183 + 99225.0*Ftmp338)*M[77] + y*(-Ftmp238*Ftmp27 + Ftmp301 + Ftmp376 + Ftmp68)*M[68] + y*(Ftmp131 + Ftmp273 + Ftmp284 + Ftmp297 + Ftmp369)*M[66];
#pragma omp atomic
F[2] += Ftmp0*M[3] - Ftmp10*Ftmp27*Ftmp327 - Ftmp10*Ftmp48*M[14] + Ftmp10*Ftmp57 + Ftmp10*Ftmp9 - Ftmp104*Ftmp381 + Ftmp116*M[51] + Ftmp120*z - Ftmp135*x - Ftmp139*Ftmp307 - Ftmp142*Ftmp172 - Ftmp142*Ftmp380*M[54] - Ftmp142*Ftmp387 - Ftmp143*x*M[49] - Ftmp145*Ftmp6 - Ftmp149*M[35] + Ftmp15*Ftmp319 + Ftmp15*Ftmp339 - Ftmp150*Ftmp27 + Ftmp16*M[5] + Ftmp163*Ftmp382 - Ftmp165*Ftmp312 + Ftmp172*Ftmp388 + Ftmp174*Ftmp323 + Ftmp178*Ftmp384 + Ftmp186*z - Ftmp2*Ftmp27*M[3] - Ftmp2*Ftmp303 - Ftmp2*Ftmp377 - Ftmp227*x - Ftmp229*x - Ftmp231*Ftmp312*M[47] - Ftmp232*Ftmp307 - Ftmp235*M[38] - Ftmp236*M[40] + Ftmp246*M[62] + Ftmp251*M[64] + Ftmp26*M[17] - Ftmp27*Ftmp314 - Ftmp27*Ftmp335 - Ftmp27*Ftmp65 - Ftmp27*(Ftmp353 + Ftmp361)*M[44] - Ftmp27*(Ftmp48 + Ftmp78)*M[19] - Ftmp27*(Ftmp140 + Ftmp164 - 13230.0*Ftmp72)*M[55] - Ftmp27*(Ftmp19 + Ftmp213 + Ftmp350)*M[42] - Ftmp27*(Ftmp222 + Ftmp230 + Ftmp255 + Ftmp82)*M[53] + Ftmp28*Ftmp303 + Ftmp28*Ftmp377 + Ftmp30*z + Ftmp304*M[12] + Ftmp306*z - Ftmp307*Ftmp47 - Ftmp310*y - Ftmp311*y - Ftmp312*Ftmp80 + Ftmp315*Ftmp381 - Ftmp316*Ftmp6 + Ftmp321*Ftmp382 + Ftmp328*M[37] + Ftmp330*z - Ftmp333*y - Ftmp334*y + Ftmp336*Ftmp381 - Ftmp337*Ftmp6 + Ftmp340*Ftmp382 + Ftmp342*z - Ftmp352*y - Ftmp355*y - Ftmp357*y - Ftmp358*M[41] - Ftmp359*M[43] - Ftmp360*M[52] + Ftmp368*M[73] + Ftmp378*z*M[9] + Ftmp378*M[19] + Ftmp379*Ftmp385 - Ftmp379*Ftmp49 - Ftmp380*Ftmp49*M[18] - Ftmp381*Ftmp40*(Ftmp106 - 1575.0*Ftmp53)*M[48] + Ftmp381*(Ftmp245 + Ftmp364)*M[63] + Ftmp381*(Ftmp250 + Ftmp363)*M[65] + Ftmp382*Ftmp77 + Ftmp382*(Ftmp152 + Ftmp249 + Ftmp392)*M[80] + Ftmp382*(Ftmp243 + Ftmp268 + Ftmp366)*M[67] + Ftmp382*(Ftmp362 + Ftmp367 + Ftmp96)*M[69] + Ftmp383*Ftmp74*M[28] + Ftmp383*(-34650.0*Ftmp170 + Ftmp175 + Ftmp96)*M[75] + Ftmp384*Ftmp94 + Ftmp384*(Ftmp270 + Ftmp392 + Ftmp96)*M[74] + Ftmp385*Ftmp89 + Ftmp386*z*M[34] + Ftmp386*M[55] + Ftmp387*Ftmp388 + Ftmp389*z*M[23] + Ftmp389*M[42] + Ftmp390*z*M[25] + Ftmp390*M[44] + Ftmp391*z*M[32] + Ftmp391*M[53] - Ftmp4*x*M[1] - Ftmp4*y*M[2] - Ftmp43*x - Ftmp49*Ftmp89 + Ftmp5*z - Ftmp50*x*M[15] - Ftmp59*Ftmp6 + Ftmp63*Ftmp90 - Ftmp64*M[10] + Ftmp70*x + Ftmp76*x + z*(Ftmp137 + Ftmp276 + Ftmp297 + Ftmp370)*M[66] + z*(Ftmp203 + Ftmp275 + Ftmp287 + Ftmp373)*M[59] + z*(Ftmp271 + Ftmp289 + Ftmp372 + Ftmp374)*M[61] + z*(Ftmp283 + Ftmp296 + Ftmp375 + Ftmp393)*M[81] + z*(-Ftmp117*Ftmp343 + Ftmp180 + Ftmp187 + 99225.0*Ftmp72)*M[83] + z*(Ftmp215 + Ftmp271 + Ftmp278 + Ftmp292 + Ftmp298)*M[79] + z*(Ftmp279 + Ftmp288 + Ftmp294 + Ftmp369 + Ftmp393)*M[70] + z*(-Ftmp238*Ftmp24 + Ftmp300 + Ftmp302 + Ftmp376 + Ftmp88)*M[68];

}

void field_m0_P2M_7(double x, double y, double z, double q, double * M) {
double Mtmp0 = q*x;
double Mtmp1 = q*y;
double Mtmp2 = q*z;
double Mtmp3 = (x*x);
double Mtmp4 = (1.0/2.0)*q;
double Mtmp5 = Mtmp0*y;
double Mtmp6 = Mtmp0*z;
double Mtmp7 = (y*y);
double Mtmp8 = Mtmp1*z;
double Mtmp9 = (z*z);
double Mtmp10 = (x*x*x);
double Mtmp11 = (1.0/6.0)*q;
double Mtmp12 = (1.0/2.0)*Mtmp3;
double Mtmp13 = (1.0/2.0)*Mtmp0;
double Mtmp14 = (y*y*y);
double Mtmp15 = (1.0/2.0)*Mtmp7;
double Mtmp16 = (1.0/2.0)*Mtmp9;
double Mtmp17 = (z*z*z);
double Mtmp18 = (x*x*x*x);
double Mtmp19 = (1.0/24.0)*q;
double Mtmp20 = (1.0/6.0)*Mtmp10;
double Mtmp21 = Mtmp7*q;
double Mtmp22 = (1.0/4.0)*Mtmp3;
double Mtmp23 = Mtmp9*q;
double Mtmp24 = (1.0/6.0)*Mtmp0;
double Mtmp25 = (y*y*y*y);
double Mtmp26 = (1.0/6.0)*Mtmp14;
double Mtmp27 = (1.0/4.0)*Mtmp9;
double Mtmp28 = (1.0/6.0)*Mtmp17;
double Mtmp29 = (z*z*z*z);
double Mtmp30 = (x*x*x*x*x);
double Mtmp31 = (1.0/120.0)*q;
double Mtmp32 = (1.0/24.0)*Mtmp18;
double Mtmp33 = (1.0/12.0)*Mtmp10;
double Mtmp34 = (1.0/12.0)*Mtmp14;
double Mtmp35 = Mtmp3*q;
double Mtmp36 = Mtmp2*Mtmp7;
double Mtmp37 = Mtmp1*Mtmp9;
double Mtmp38 = (1.0/12.0)*Mtmp17;
double Mtmp39 = (1.0/24.0)*Mtmp0;
double Mtmp40 = Mtmp0*Mtmp7;
double Mtmp41 = (y*y*y*y*y);
double Mtmp42 = (1.0/24.0)*Mtmp25;
double Mtmp43 = (1.0/24.0)*Mtmp29;
double Mtmp44 = (z*z*z*z*z);
double Mtmp45 = (x*x*x*x*x*x);
double Mtmp46 = (1.0/720.0)*q;
double Mtmp47 = (1.0/120.0)*Mtmp30;
double Mtmp48 = (1.0/48.0)*Mtmp18;
double Mtmp49 = Mtmp14*q;
double Mtmp50 = (1.0/36.0)*Mtmp10;
double Mtmp51 = Mtmp17*q;
double Mtmp52 = (1.0/48.0)*Mtmp35;
double Mtmp53 = Mtmp2*Mtmp3;
double Mtmp54 = Mtmp3*Mtmp9;
double Mtmp55 = Mtmp1*Mtmp3;
double Mtmp56 = (1.0/120.0)*Mtmp0;
double Mtmp57 = Mtmp0*Mtmp9;
double Mtmp58 = (y*y*y*y*y*y);
double Mtmp59 = (1.0/120.0)*Mtmp41;
double Mtmp60 = (1.0/48.0)*Mtmp25;
double Mtmp61 = (1.0/36.0)*Mtmp17;
double Mtmp62 = (1.0/48.0)*Mtmp29;
double Mtmp63 = (1.0/120.0)*Mtmp44;
double Mtmp64 = (z*z*z*z*z*z);
double Mtmp65 = (1.0/5040.0)*q;
double Mtmp66 = (1.0/720.0)*Mtmp45;
double Mtmp67 = (1.0/240.0)*Mtmp30;
double Mtmp68 = (1.0/144.0)*Mtmp18;
double Mtmp69 = (1.0/144.0)*Mtmp25;
double Mtmp70 = Mtmp10*q;
double Mtmp71 = Mtmp19*Mtmp7;
double Mtmp72 = (1.0/144.0)*Mtmp29;
double Mtmp73 = (1.0/240.0)*Mtmp35;
double Mtmp74 = (1.0/720.0)*Mtmp0;
M[0] += q;
M[1] += -Mtmp0;
M[2] += -Mtmp1;
M[3] += -Mtmp2;
M[4] += Mtmp3*Mtmp4;
M[5] += Mtmp5;
M[6] += Mtmp6;
M[7] += Mtmp4*Mtmp7;
M[8] += Mtmp8;
M[9] += Mtmp4*Mtmp9;
M[10] += -Mtmp10*Mtmp11;
M[11] += -Mtmp1*Mtmp12;
M[12] += -Mtmp12*Mtmp2;
M[13] += -Mtmp13*Mtmp7;
M[14] += -Mtmp5*z;
M[15] += -Mtmp13*Mtmp9;
M[16] += -Mtmp11*Mtmp14;
M[17] += -Mtmp15*Mtmp2;
M[18] += -Mtmp1*Mtmp16;
M[19] += -Mtmp11*Mtmp17;
M[20] += Mtmp18*Mtmp19;
M[21] += Mtmp1*Mtmp20;
M[22] += Mtmp2*Mtmp20;
M[23] += Mtmp21*Mtmp22;
M[24] += Mtmp12*Mtmp8;
M[25] += Mtmp22*Mtmp23;
M[26] += Mtmp14*Mtmp24;
M[27] += Mtmp15*Mtmp6;
M[28] += Mtmp16*Mtmp5;
M[29] += Mtmp17*Mtmp24;
M[30] += Mtmp19*Mtmp25;
M[31] += Mtmp2*Mtmp26;
M[32] += Mtmp21*Mtmp27;
M[33] += Mtmp1*Mtmp28;
M[34] += Mtmp19*Mtmp29;
M[35] += -Mtmp30*Mtmp31;
M[36] += -Mtmp1*Mtmp32;
M[37] += -Mtmp2*Mtmp32;
M[38] += -Mtmp21*Mtmp33;
M[39] += -Mtmp20*Mtmp8;
M[40] += -Mtmp23*Mtmp33;
M[41] += -Mtmp34*Mtmp35;
M[42] += -Mtmp22*Mtmp36;
M[43] += -Mtmp22*Mtmp37;
M[44] += -Mtmp35*Mtmp38;
M[45] += -Mtmp25*Mtmp39;
M[46] += -Mtmp26*Mtmp6;
M[47] += -Mtmp27*Mtmp40;
M[48] += -Mtmp28*Mtmp5;
M[49] += -Mtmp29*Mtmp39;
M[50] += -Mtmp31*Mtmp41;
M[51] += -Mtmp2*Mtmp42;
M[52] += -Mtmp23*Mtmp34;
M[53] += -Mtmp21*Mtmp38;
M[54] += -Mtmp1*Mtmp43;
M[55] += -Mtmp31*Mtmp44;
M[56] += Mtmp45*Mtmp46;
M[57] += Mtmp1*Mtmp47;
M[58] += Mtmp2*Mtmp47;
M[59] += Mtmp21*Mtmp48;
M[60] += Mtmp32*Mtmp8;
M[61] += Mtmp23*Mtmp48;
M[62] += Mtmp49*Mtmp50;
M[63] += Mtmp33*Mtmp36;
M[64] += Mtmp33*Mtmp37;
M[65] += Mtmp50*Mtmp51;
M[66] += Mtmp25*Mtmp52;
M[67] += Mtmp34*Mtmp53;
M[68] += (1.0/8.0)*Mtmp21*Mtmp54;
M[69] += Mtmp38*Mtmp55;
M[70] += Mtmp29*Mtmp52;
M[71] += Mtmp41*Mtmp56;
M[72] += Mtmp42*Mtmp6;
M[73] += Mtmp34*Mtmp57;
M[74] += Mtmp38*Mtmp40;
M[75] += Mtmp43*Mtmp5;
M[76] += Mtmp44*Mtmp56;
M[77] += Mtmp46*Mtmp58;
M[78] += Mtmp2*Mtmp59;
M[79] += Mtmp23*Mtmp60;
M[80] += Mtmp49*Mtmp61;
M[81] += Mtmp21*Mtmp62;
M[82] += Mtmp1*Mtmp63;
M[83] += Mtmp46*Mtmp64;
M[84] += -Mtmp65*(x*x*x*x*x*x*x);
M[85] += -Mtmp1*Mtmp66;
M[86] += -Mtmp2*Mtmp66;
M[87] += -Mtmp21*Mtmp67;
M[88] += -Mtmp47*Mtmp8;
M[89] += -Mtmp23*Mtmp67;
M[90] += -Mtmp49*Mtmp68;
M[91] += -Mtmp36*Mtmp48;
M[92] += -Mtmp37*Mtmp48;
M[93] += -Mtmp51*Mtmp68;
M[94] += -Mtmp69*Mtmp70;
M[95] += -Mtmp14*Mtmp2*Mtmp50;
M[96] += -Mtmp10*Mtmp71*Mtmp9;
M[97] += -Mtmp1*Mtmp17*Mtmp50;
M[98] += -Mtmp70*Mtmp72;
M[99] += -Mtmp41*Mtmp73;
M[100] += -Mtmp53*Mtmp60;
M[101] += -Mtmp14*Mtmp19*Mtmp54;
M[102] += -Mtmp17*Mtmp3*Mtmp71;
M[103] += -Mtmp55*Mtmp62;
M[104] += -Mtmp44*Mtmp73;
M[105] += -Mtmp58*Mtmp74;
M[106] += -Mtmp59*Mtmp6;
M[107] += -Mtmp57*Mtmp60;
M[108] += -Mtmp0*Mtmp14*Mtmp61;
M[109] += -Mtmp40*Mtmp62;
M[110] += -Mtmp5*Mtmp63;
M[111] += -Mtmp64*Mtmp74;
M[112] += -Mtmp65*(y*y*y*y*y*y*y);
M[113] += -1.0/720.0*Mtmp2*Mtmp58;
M[114] += -1.0/240.0*Mtmp23*Mtmp41;
M[115] += -Mtmp51*Mtmp69;
M[116] += -Mtmp49*Mtmp72;
M[117] += -1.0/240.0*Mtmp21*Mtmp44;
M[118] += -1.0/720.0*Mtmp1*Mtmp64;
M[119] += -Mtmp65*(z*z*z*z*z*z*z);

}
void field_m0_M2M_7(double x, double y, double z, double * M, double * Ms) {
double Mstmp0 = x*M[0];
double Mstmp1 = y*M[0];
double Mstmp2 = z*M[0];
double Mstmp3 = x*M[1];
double Mstmp4 = (x*x);
double Mstmp5 = (1.0/2.0)*M[0];
double Mstmp6 = x*M[2];
double Mstmp7 = y*M[1];
double Mstmp8 = Mstmp0*y;
double Mstmp9 = x*M[3];
double Mstmp10 = z*M[1];
double Mstmp11 = Mstmp0*z;
double Mstmp12 = y*M[2];
double Mstmp13 = (y*y);
double Mstmp14 = y*M[3];
double Mstmp15 = z*M[2];
double Mstmp16 = Mstmp1*z;
double Mstmp17 = z*M[3];
double Mstmp18 = (z*z);
double Mstmp19 = x*M[4];
double Mstmp20 = (1.0/2.0)*Mstmp4;
double Mstmp21 = (x*x*x);
double Mstmp22 = (1.0/6.0)*M[0];
double Mstmp23 = x*M[5];
double Mstmp24 = y*M[4];
double Mstmp25 = Mstmp3*y;
double Mstmp26 = x*M[6];
double Mstmp27 = z*M[4];
double Mstmp28 = Mstmp3*z;
double Mstmp29 = x*M[7];
double Mstmp30 = y*M[5];
double Mstmp31 = Mstmp6*y;
double Mstmp32 = (1.0/2.0)*M[1];
double Mstmp33 = (1.0/2.0)*Mstmp13;
double Mstmp34 = x*M[8];
double Mstmp35 = y*M[6];
double Mstmp36 = z*M[5];
double Mstmp37 = Mstmp9*y;
double Mstmp38 = Mstmp6*z;
double Mstmp39 = Mstmp7*z;
double Mstmp40 = x*M[9];
double Mstmp41 = z*M[6];
double Mstmp42 = Mstmp9*z;
double Mstmp43 = (1.0/2.0)*Mstmp18;
double Mstmp44 = y*M[7];
double Mstmp45 = (y*y*y);
double Mstmp46 = y*M[8];
double Mstmp47 = z*M[7];
double Mstmp48 = Mstmp12*z;
double Mstmp49 = y*M[9];
double Mstmp50 = z*M[8];
double Mstmp51 = Mstmp14*z;
double Mstmp52 = z*M[9];
double Mstmp53 = (z*z*z);
double Mstmp54 = x*M[10];
double Mstmp55 = (1.0/6.0)*Mstmp21;
double Mstmp56 = (x*x*x*x);
double Mstmp57 = (1.0/24.0)*M[0];
double Mstmp58 = x*M[11];
double Mstmp59 = y*M[10];
double Mstmp60 = Mstmp19*y;
double Mstmp61 = x*M[12];
double Mstmp62 = z*M[10];
double Mstmp63 = Mstmp19*z;
double Mstmp64 = x*M[13];
double Mstmp65 = y*M[11];
double Mstmp66 = Mstmp23*y;
double Mstmp67 = Mstmp13*M[0];
double Mstmp68 = (1.0/4.0)*Mstmp4;
double Mstmp69 = x*M[14];
double Mstmp70 = y*M[12];
double Mstmp71 = z*M[11];
double Mstmp72 = Mstmp26*y;
double Mstmp73 = Mstmp23*z;
double Mstmp74 = Mstmp24*z;
double Mstmp75 = x*M[15];
double Mstmp76 = z*M[12];
double Mstmp77 = Mstmp26*z;
double Mstmp78 = Mstmp18*Mstmp68;
double Mstmp79 = x*M[16];
double Mstmp80 = y*M[13];
double Mstmp81 = Mstmp29*y;
double Mstmp82 = (1.0/6.0)*M[1];
double Mstmp83 = (1.0/6.0)*Mstmp45;
double Mstmp84 = x*M[17];
double Mstmp85 = y*M[14];
double Mstmp86 = z*M[13];
double Mstmp87 = Mstmp34*y;
double Mstmp88 = Mstmp29*z;
double Mstmp89 = Mstmp30*z;
double Mstmp90 = x*M[18];
double Mstmp91 = y*M[15];
double Mstmp92 = z*M[14];
double Mstmp93 = Mstmp40*y;
double Mstmp94 = Mstmp34*z;
double Mstmp95 = Mstmp35*z;
double Mstmp96 = x*M[19];
double Mstmp97 = z*M[15];
double Mstmp98 = Mstmp40*z;
double Mstmp99 = (1.0/6.0)*Mstmp53;
double Mstmp100 = y*M[16];
double Mstmp101 = (y*y*y*y);
double Mstmp102 = y*M[17];
double Mstmp103 = z*M[16];
double Mstmp104 = Mstmp44*z;
double Mstmp105 = y*M[18];
double Mstmp106 = z*M[17];
double Mstmp107 = Mstmp46*z;
double Mstmp108 = (1.0/4.0)*Mstmp18;
double Mstmp109 = y*M[19];
double Mstmp110 = z*M[18];
double Mstmp111 = Mstmp49*z;
double Mstmp112 = z*M[19];
double Mstmp113 = (z*z*z*z);
double Mstmp114 = x*M[20];
double Mstmp115 = (1.0/24.0)*Mstmp56;
double Mstmp116 = (x*x*x*x*x);
double Mstmp117 = (1.0/120.0)*M[0];
double Mstmp118 = x*M[21];
double Mstmp119 = y*M[20];
double Mstmp120 = Mstmp54*y;
double Mstmp121 = x*M[22];
double Mstmp122 = z*M[20];
double Mstmp123 = Mstmp54*z;
double Mstmp124 = x*M[23];
double Mstmp125 = y*M[21];
double Mstmp126 = Mstmp58*y;
double Mstmp127 = Mstmp13*Mstmp68;
double Mstmp128 = (1.0/12.0)*Mstmp21;
double Mstmp129 = x*M[24];
double Mstmp130 = y*M[22];
double Mstmp131 = z*M[21];
double Mstmp132 = Mstmp61*y;
double Mstmp133 = Mstmp58*z;
double Mstmp134 = Mstmp59*z;
double Mstmp135 = x*M[25];
double Mstmp136 = z*M[22];
double Mstmp137 = Mstmp61*z;
double Mstmp138 = Mstmp18*M[0];
double Mstmp139 = x*M[26];
double Mstmp140 = y*M[23];
double Mstmp141 = Mstmp64*y;
double Mstmp142 = (1.0/12.0)*Mstmp45;
double Mstmp143 = Mstmp4*M[0];
double Mstmp144 = x*M[27];
double Mstmp145 = y*M[24];
double Mstmp146 = z*M[23];
double Mstmp147 = Mstmp69*y;
double Mstmp148 = Mstmp64*z;
double Mstmp149 = Mstmp65*z;
double Mstmp150 = x*M[28];
double Mstmp151 = y*M[25];
double Mstmp152 = z*M[24];
double Mstmp153 = Mstmp75*y;
double Mstmp154 = Mstmp69*z;
double Mstmp155 = Mstmp70*z;
double Mstmp156 = x*M[29];
double Mstmp157 = z*M[25];
double Mstmp158 = Mstmp75*z;
double Mstmp159 = (1.0/12.0)*Mstmp53;
double Mstmp160 = x*M[30];
double Mstmp161 = y*M[26];
double Mstmp162 = Mstmp79*y;
double Mstmp163 = (1.0/24.0)*M[1];
double Mstmp164 = (1.0/24.0)*Mstmp101;
double Mstmp165 = x*M[31];
double Mstmp166 = y*M[27];
double Mstmp167 = z*M[26];
double Mstmp168 = Mstmp84*y;
double Mstmp169 = Mstmp79*z;
double Mstmp170 = Mstmp80*z;
double Mstmp171 = x*M[32];
double Mstmp172 = y*M[28];
double Mstmp173 = z*M[27];
double Mstmp174 = Mstmp90*y;
double Mstmp175 = Mstmp84*z;
double Mstmp176 = Mstmp85*z;
double Mstmp177 = Mstmp108*Mstmp13;
double Mstmp178 = x*M[33];
double Mstmp179 = y*M[29];
double Mstmp180 = z*M[28];
double Mstmp181 = Mstmp96*y;
double Mstmp182 = Mstmp90*z;
double Mstmp183 = Mstmp91*z;
double Mstmp184 = x*M[34];
double Mstmp185 = z*M[29];
double Mstmp186 = Mstmp96*z;
double Mstmp187 = (1.0/24.0)*Mstmp113;
double Mstmp188 = y*M[30];
double Mstmp189 = (y*y*y*y*y);
double Mstmp190 = y*M[31];
double Mstmp191 = z*M[30];
double Mstmp192 = Mstmp100*z;
double Mstmp193 = y*M[32];
double Mstmp194 = z*M[31];
double Mstmp195 = Mstmp102*z;
double Mstmp196 = y*M[33];
double Mstmp197 = z*M[32];
double Mstmp198 = Mstmp105*z;
double Mstmp199 = y*M[34];
double Mstmp200 = z*M[33];
double Mstmp201 = Mstmp109*z;
double Mstmp202 = z*M[34];
double Mstmp203 = (z*z*z*z*z);
double Mstmp204 = x*M[35];
double Mstmp205 = (1.0/120.0)*Mstmp116;
double Mstmp206 = (x*x*x*x*x*x);
double Mstmp207 = (1.0/720.0)*M[0];
double Mstmp208 = x*M[36];
double Mstmp209 = y*M[35];
double Mstmp210 = Mstmp114*y;
double Mstmp211 = x*M[37];
double Mstmp212 = x*M[38];
double Mstmp213 = y*M[36];
double Mstmp214 = Mstmp118*y;
double Mstmp215 = Mstmp128*M[1];
double Mstmp216 = (1.0/48.0)*Mstmp56;
double Mstmp217 = x*M[39];
double Mstmp218 = y*M[37];
double Mstmp219 = Mstmp121*y;
double Mstmp220 = x*M[40];
double Mstmp221 = x*M[41];
double Mstmp222 = y*M[38];
double Mstmp223 = Mstmp124*y;
double Mstmp224 = Mstmp4*M[1];
double Mstmp225 = Mstmp128*Mstmp13;
double Mstmp226 = Mstmp45*M[0];
double Mstmp227 = (1.0/36.0)*Mstmp21;
double Mstmp228 = x*M[42];
double Mstmp229 = y*M[39];
double Mstmp230 = Mstmp129*y;
double Mstmp231 = x*M[43];
double Mstmp232 = y*M[40];
double Mstmp233 = Mstmp135*y;
double Mstmp234 = Mstmp128*Mstmp18;
double Mstmp235 = x*M[44];
double Mstmp236 = Mstmp227*Mstmp53;
double Mstmp237 = x*M[45];
double Mstmp238 = y*M[41];
double Mstmp239 = Mstmp139*y;
double Mstmp240 = Mstmp142*Mstmp4;
double Mstmp241 = (1.0/48.0)*Mstmp143;
double Mstmp242 = x*M[46];
double Mstmp243 = y*M[42];
double Mstmp244 = Mstmp144*y;
double Mstmp245 = x*M[47];
double Mstmp246 = y*M[43];
double Mstmp247 = Mstmp150*y;
double Mstmp248 = (1.0/8.0)*Mstmp18;
double Mstmp249 = Mstmp248*Mstmp4;
double Mstmp250 = x*M[48];
double Mstmp251 = y*M[44];
double Mstmp252 = Mstmp156*y;
double Mstmp253 = Mstmp159*Mstmp4;
double Mstmp254 = x*M[49];
double Mstmp255 = x*M[50];
double Mstmp256 = y*M[45];
double Mstmp257 = Mstmp160*y;
double Mstmp258 = (1.0/120.0)*M[1];
double Mstmp259 = (1.0/120.0)*Mstmp189;
double Mstmp260 = x*M[51];
double Mstmp261 = y*M[46];
double Mstmp262 = Mstmp165*y;
double Mstmp263 = x*M[52];
double Mstmp264 = y*M[47];
double Mstmp265 = Mstmp171*y;
double Mstmp266 = Mstmp142*Mstmp18;
double Mstmp267 = x*M[53];
double Mstmp268 = y*M[48];
double Mstmp269 = Mstmp178*y;
double Mstmp270 = Mstmp13*Mstmp159;
double Mstmp271 = x*M[54];
double Mstmp272 = y*M[49];
double Mstmp273 = Mstmp184*y;
double Mstmp274 = x*M[55];
double Mstmp275 = (1.0/120.0)*Mstmp203;
double Mstmp276 = y*M[50];
double Mstmp277 = (y*y*y*y*y*y);
double Mstmp278 = y*M[51];
double Mstmp279 = y*M[52];
double Mstmp280 = (1.0/48.0)*Mstmp101;
double Mstmp281 = y*M[53];
double Mstmp282 = (1.0/36.0)*Mstmp53;
double Mstmp283 = y*M[54];
double Mstmp284 = (1.0/48.0)*Mstmp113;
double Mstmp285 = y*M[55];
double Mstmp286 = (z*z*z*z*z*z);
double Mstmp287 = (1.0/720.0)*Mstmp206;
double Mstmp288 = (1.0/5040.0)*M[0];
double Mstmp289 = Mstmp216*M[1];
double Mstmp290 = (1.0/240.0)*Mstmp116;
double Mstmp291 = Mstmp227*Mstmp45;
double Mstmp292 = Mstmp13*Mstmp216;
double Mstmp293 = (1.0/144.0)*Mstmp56;
double Mstmp294 = Mstmp18*Mstmp216;
double Mstmp295 = Mstmp53*M[0];
double Mstmp296 = (1.0/144.0)*Mstmp101;
double Mstmp297 = Mstmp21*M[0];
double Mstmp298 = Mstmp18*Mstmp57;
double Mstmp299 = (1.0/144.0)*Mstmp113;
double Mstmp300 = Mstmp280*Mstmp4;
double Mstmp301 = (1.0/240.0)*Mstmp143;
double Mstmp302 = Mstmp13*Mstmp249;
double Mstmp303 = Mstmp284*Mstmp4;
double Mstmp304 = (1.0/720.0)*M[1];
double Mstmp305 = (1.0/720.0)*Mstmp277;
double Mstmp306 = Mstmp18*Mstmp280;
double Mstmp307 = Mstmp282*Mstmp45;
double Mstmp308 = Mstmp13*Mstmp284;
double Mstmp309 = (1.0/720.0)*Mstmp286;
#pragma omp atomic
Ms[0] += M[0];
#pragma omp atomic
Ms[1] += Mstmp0 + M[1];
#pragma omp atomic
Ms[2] += Mstmp1 + M[2];
#pragma omp atomic
Ms[3] += Mstmp2 + M[3];
#pragma omp atomic
Ms[4] += Mstmp3 + Mstmp4*Mstmp5 + M[4];
#pragma omp atomic
Ms[5] += Mstmp6 + Mstmp7 + Mstmp8 + M[5];
#pragma omp atomic
Ms[6] += Mstmp10 + Mstmp11 + Mstmp9 + M[6];
#pragma omp atomic
Ms[7] += Mstmp12 + Mstmp13*Mstmp5 + M[7];
#pragma omp atomic
Ms[8] += Mstmp14 + Mstmp15 + Mstmp16 + M[8];
#pragma omp atomic
Ms[9] += Mstmp17 + Mstmp18*Mstmp5 + M[9];
#pragma omp atomic
Ms[10] += Mstmp19 + Mstmp20*M[1] + Mstmp21*Mstmp22 + M[10];
#pragma omp atomic
Ms[11] += Mstmp1*Mstmp20 + Mstmp20*M[2] + Mstmp23 + Mstmp24 + Mstmp25 + M[11];
#pragma omp atomic
Ms[12] += Mstmp2*Mstmp20 + Mstmp20*M[3] + Mstmp26 + Mstmp27 + Mstmp28 + M[12];
#pragma omp atomic
Ms[13] += Mstmp0*Mstmp33 + Mstmp13*Mstmp32 + Mstmp29 + Mstmp30 + Mstmp31 + M[13];
#pragma omp atomic
Ms[14] += Mstmp34 + Mstmp35 + Mstmp36 + Mstmp37 + Mstmp38 + Mstmp39 + Mstmp8*z + M[14];
#pragma omp atomic
Ms[15] += Mstmp0*Mstmp43 + Mstmp18*Mstmp32 + Mstmp40 + Mstmp41 + Mstmp42 + M[15];
#pragma omp atomic
Ms[16] += Mstmp22*Mstmp45 + Mstmp33*M[2] + Mstmp44 + M[16];
#pragma omp atomic
Ms[17] += Mstmp2*Mstmp33 + Mstmp33*M[3] + Mstmp46 + Mstmp47 + Mstmp48 + M[17];
#pragma omp atomic
Ms[18] += Mstmp1*Mstmp43 + Mstmp43*M[2] + Mstmp49 + Mstmp50 + Mstmp51 + M[18];
#pragma omp atomic
Ms[19] += Mstmp22*Mstmp53 + Mstmp43*M[3] + Mstmp52 + M[19];
#pragma omp atomic
Ms[20] += Mstmp20*M[4] + Mstmp54 + Mstmp55*M[1] + Mstmp56*Mstmp57 + M[20];
#pragma omp atomic
Ms[21] += Mstmp1*Mstmp55 + Mstmp20*Mstmp7 + Mstmp20*M[5] + Mstmp55*M[2] + Mstmp58 + Mstmp59 + Mstmp60 + M[21];
#pragma omp atomic
Ms[22] += Mstmp10*Mstmp20 + Mstmp2*Mstmp55 + Mstmp20*M[6] + Mstmp55*M[3] + Mstmp61 + Mstmp62 + Mstmp63 + M[22];
#pragma omp atomic
Ms[23] += Mstmp12*Mstmp20 + Mstmp20*M[7] + Mstmp3*Mstmp33 + Mstmp33*M[4] + Mstmp64 + Mstmp65 + Mstmp66 + Mstmp67*Mstmp68 + M[23];
#pragma omp atomic
Ms[24] += Mstmp14*Mstmp20 + Mstmp15*Mstmp20 + Mstmp16*Mstmp20 + Mstmp20*M[8] + Mstmp25*z + Mstmp69 + Mstmp70 + Mstmp71 + Mstmp72 + Mstmp73 + Mstmp74 + M[24];
#pragma omp atomic
Ms[25] += Mstmp17*Mstmp20 + Mstmp20*M[9] + Mstmp3*Mstmp43 + Mstmp43*M[4] + Mstmp75 + Mstmp76 + Mstmp77 + Mstmp78*M[0] + M[25];
#pragma omp atomic
Ms[26] += Mstmp0*Mstmp83 + Mstmp33*Mstmp6 + Mstmp33*M[5] + Mstmp45*Mstmp82 + Mstmp79 + Mstmp80 + Mstmp81 + M[26];
#pragma omp atomic
Ms[27] += Mstmp10*Mstmp33 + Mstmp11*Mstmp33 + Mstmp31*z + Mstmp33*Mstmp9 + Mstmp33*M[6] + Mstmp84 + Mstmp85 + Mstmp86 + Mstmp87 + Mstmp88 + Mstmp89 + M[27];
#pragma omp atomic
Ms[28] += Mstmp37*z + Mstmp43*Mstmp6 + Mstmp43*Mstmp7 + Mstmp43*Mstmp8 + Mstmp43*M[5] + Mstmp90 + Mstmp91 + Mstmp92 + Mstmp93 + Mstmp94 + Mstmp95 + M[28];
#pragma omp atomic
Ms[29] += Mstmp0*Mstmp99 + Mstmp43*Mstmp9 + Mstmp43*M[6] + Mstmp53*Mstmp82 + Mstmp96 + Mstmp97 + Mstmp98 + M[29];
#pragma omp atomic
Ms[30] += Mstmp100 + Mstmp101*Mstmp57 + Mstmp33*M[7] + Mstmp83*M[2] + M[30];
#pragma omp atomic
Ms[31] += Mstmp102 + Mstmp103 + Mstmp104 + Mstmp15*Mstmp33 + Mstmp2*Mstmp83 + Mstmp33*M[8] + Mstmp83*M[3] + M[31];
#pragma omp atomic
Ms[32] += Mstmp105 + Mstmp106 + Mstmp107 + Mstmp108*Mstmp67 + Mstmp12*Mstmp43 + Mstmp17*Mstmp33 + Mstmp33*M[9] + Mstmp43*M[7] + M[32];
#pragma omp atomic
Ms[33] += Mstmp1*Mstmp99 + Mstmp109 + Mstmp110 + Mstmp111 + Mstmp14*Mstmp43 + Mstmp43*M[8] + Mstmp99*M[2] + M[33];
#pragma omp atomic
Ms[34] += Mstmp112 + Mstmp113*Mstmp57 + Mstmp43*M[9] + Mstmp99*M[3] + M[34];
#pragma omp atomic
Ms[35] += Mstmp114 + Mstmp115*M[1] + Mstmp116*Mstmp117 + Mstmp20*M[10] + Mstmp55*M[4] + M[35];
#pragma omp atomic
Ms[36] += Mstmp1*Mstmp115 + Mstmp115*M[2] + Mstmp118 + Mstmp119 + Mstmp120 + Mstmp20*Mstmp24 + Mstmp20*M[11] + Mstmp55*Mstmp7 + Mstmp55*M[5] + M[36];
#pragma omp atomic
Ms[37] += Mstmp10*Mstmp55 + Mstmp115*Mstmp2 + Mstmp115*M[3] + Mstmp121 + Mstmp122 + Mstmp123 + Mstmp20*Mstmp27 + Mstmp20*M[12] + Mstmp55*M[6] + M[37];
#pragma omp atomic
Ms[38] += Mstmp12*Mstmp55 + Mstmp124 + Mstmp125 + Mstmp126 + Mstmp127*M[1] + Mstmp128*Mstmp67 + Mstmp19*Mstmp33 + Mstmp20*Mstmp30 + Mstmp20*M[13] + Mstmp33*M[10] + Mstmp55*M[7] + M[38];
#pragma omp atomic
Ms[39] += Mstmp129 + Mstmp130 + Mstmp131 + Mstmp132 + Mstmp133 + Mstmp134 + Mstmp14*Mstmp55 + Mstmp15*Mstmp55 + Mstmp16*Mstmp55 + Mstmp20*Mstmp35 + Mstmp20*Mstmp36 + Mstmp20*Mstmp39 + Mstmp20*M[14] + Mstmp55*M[8] + Mstmp60*z + M[39];
#pragma omp atomic
Ms[40] += Mstmp128*Mstmp138 + Mstmp135 + Mstmp136 + Mstmp137 + Mstmp17*Mstmp55 + Mstmp19*Mstmp43 + Mstmp20*Mstmp41 + Mstmp20*M[15] + Mstmp43*M[10] + Mstmp55*M[9] + Mstmp78*M[1] + M[40];
#pragma omp atomic
Ms[41] += Mstmp127*M[2] + Mstmp139 + Mstmp140 + Mstmp141 + Mstmp142*Mstmp143 + Mstmp20*Mstmp44 + Mstmp20*M[16] + Mstmp23*Mstmp33 + Mstmp3*Mstmp83 + Mstmp33*M[11] + Mstmp83*M[4] + M[41];
#pragma omp atomic
Ms[42] += Mstmp127*Mstmp2 + Mstmp127*M[3] + Mstmp144 + Mstmp145 + Mstmp146 + Mstmp147 + Mstmp148 + Mstmp149 + Mstmp20*Mstmp46 + Mstmp20*Mstmp47 + Mstmp20*Mstmp48 + Mstmp20*M[17] + Mstmp26*Mstmp33 + Mstmp27*Mstmp33 + Mstmp28*Mstmp33 + Mstmp33*M[12] + Mstmp66*z + M[42];
#pragma omp atomic
Ms[43] += Mstmp1*Mstmp78 + Mstmp150 + Mstmp151 + Mstmp152 + Mstmp153 + Mstmp154 + Mstmp155 + Mstmp20*Mstmp49 + Mstmp20*Mstmp50 + Mstmp20*Mstmp51 + Mstmp20*M[18] + Mstmp23*Mstmp43 + Mstmp24*Mstmp43 + Mstmp25*Mstmp43 + Mstmp43*M[11] + Mstmp72*z + Mstmp78*M[2] + M[43];
#pragma omp atomic
Ms[44] += Mstmp143*Mstmp159 + Mstmp156 + Mstmp157 + Mstmp158 + Mstmp20*Mstmp52 + Mstmp20*M[19] + Mstmp26*Mstmp43 + Mstmp3*Mstmp99 + Mstmp43*M[12] + Mstmp78*M[3] + Mstmp99*M[4] + M[44];
#pragma omp atomic
Ms[45] += Mstmp0*Mstmp164 + Mstmp101*Mstmp163 + Mstmp160 + Mstmp161 + Mstmp162 + Mstmp29*Mstmp33 + Mstmp33*M[13] + Mstmp6*Mstmp83 + Mstmp83*M[5] + M[45];
#pragma omp atomic
Ms[46] += Mstmp10*Mstmp83 + Mstmp11*Mstmp83 + Mstmp165 + Mstmp166 + Mstmp167 + Mstmp168 + Mstmp169 + Mstmp170 + Mstmp33*Mstmp34 + Mstmp33*Mstmp36 + Mstmp33*Mstmp38 + Mstmp33*M[14] + Mstmp81*z + Mstmp83*Mstmp9 + Mstmp83*M[6] + M[46];
#pragma omp atomic
Ms[47] += Mstmp0*Mstmp177 + Mstmp171 + Mstmp172 + Mstmp173 + Mstmp174 + Mstmp175 + Mstmp176 + Mstmp177*M[1] + Mstmp29*Mstmp43 + Mstmp30*Mstmp43 + Mstmp31*Mstmp43 + Mstmp33*Mstmp40 + Mstmp33*Mstmp41 + Mstmp33*Mstmp42 + Mstmp33*M[15] + Mstmp43*M[13] + Mstmp87*z + M[47];
#pragma omp atomic
Ms[48] += Mstmp178 + Mstmp179 + Mstmp180 + Mstmp181 + Mstmp182 + Mstmp183 + Mstmp34*Mstmp43 + Mstmp35*Mstmp43 + Mstmp37*Mstmp43 + Mstmp43*M[14] + Mstmp6*Mstmp99 + Mstmp7*Mstmp99 + Mstmp8*Mstmp99 + Mstmp93*z + Mstmp99*M[5] + M[48];
#pragma omp atomic
Ms[49] += Mstmp0*Mstmp187 + Mstmp113*Mstmp163 + Mstmp184 + Mstmp185 + Mstmp186 + Mstmp40*Mstmp43 + Mstmp43*M[15] + Mstmp9*Mstmp99 + Mstmp99*M[6] + M[49];
#pragma omp atomic
Ms[50] += Mstmp117*Mstmp189 + Mstmp164*M[2] + Mstmp188 + Mstmp33*M[16] + Mstmp83*M[7] + M[50];
#pragma omp atomic
Ms[51] += Mstmp15*Mstmp83 + Mstmp164*Mstmp2 + Mstmp164*M[3] + Mstmp190 + Mstmp191 + Mstmp192 + Mstmp33*Mstmp47 + Mstmp33*M[17] + Mstmp83*M[8] + M[51];
#pragma omp atomic
Ms[52] += Mstmp138*Mstmp142 + Mstmp17*Mstmp83 + Mstmp177*M[2] + Mstmp193 + Mstmp194 + Mstmp195 + Mstmp33*Mstmp50 + Mstmp33*M[18] + Mstmp43*Mstmp44 + Mstmp43*M[16] + Mstmp83*M[9] + M[52];
#pragma omp atomic
Ms[53] += Mstmp12*Mstmp99 + Mstmp159*Mstmp67 + Mstmp177*M[3] + Mstmp196 + Mstmp197 + Mstmp198 + Mstmp33*Mstmp52 + Mstmp33*M[19] + Mstmp43*Mstmp46 + Mstmp43*M[17] + Mstmp99*M[7] + M[53];
#pragma omp atomic
Ms[54] += Mstmp1*Mstmp187 + Mstmp14*Mstmp99 + Mstmp187*M[2] + Mstmp199 + Mstmp200 + Mstmp201 + Mstmp43*Mstmp49 + Mstmp43*M[18] + Mstmp99*M[8] + M[54];
#pragma omp atomic
Ms[55] += Mstmp117*Mstmp203 + Mstmp187*M[3] + Mstmp202 + Mstmp43*M[19] + Mstmp99*M[9] + M[55];
#pragma omp atomic
Ms[56] += Mstmp115*M[4] + Mstmp20*M[20] + Mstmp204 + Mstmp205*M[1] + Mstmp206*Mstmp207 + Mstmp55*M[10] + M[56];
#pragma omp atomic
Ms[57] += Mstmp1*Mstmp205 + Mstmp115*Mstmp7 + Mstmp115*M[5] + Mstmp20*Mstmp59 + Mstmp20*M[21] + Mstmp205*M[2] + Mstmp208 + Mstmp209 + Mstmp210 + Mstmp24*Mstmp55 + Mstmp55*M[11] + M[57];
#pragma omp atomic
Ms[58] += Mstmp10*Mstmp115 + Mstmp114*z + Mstmp115*M[6] + Mstmp2*Mstmp205 + Mstmp20*Mstmp62 + Mstmp20*M[22] + Mstmp205*M[3] + Mstmp211 + Mstmp27*Mstmp55 + Mstmp55*M[12] + z*M[35] + M[58];
#pragma omp atomic
Ms[59] += Mstmp115*Mstmp12 + Mstmp115*M[7] + Mstmp127*M[4] + Mstmp13*Mstmp215 + Mstmp20*Mstmp65 + Mstmp20*M[23] + Mstmp212 + Mstmp213 + Mstmp214 + Mstmp216*Mstmp67 + Mstmp30*Mstmp55 + Mstmp33*Mstmp54 + Mstmp33*M[20] + Mstmp55*M[13] + M[59];
#pragma omp atomic
Ms[60] += Mstmp115*Mstmp14 + Mstmp115*Mstmp15 + Mstmp115*Mstmp16 + Mstmp115*M[8] + Mstmp118*z + Mstmp119*z + Mstmp120*z + Mstmp20*Mstmp70 + Mstmp20*Mstmp71 + Mstmp20*Mstmp74 + Mstmp20*M[24] + Mstmp217 + Mstmp218 + Mstmp219 + Mstmp35*Mstmp55 + Mstmp36*Mstmp55 + Mstmp39*Mstmp55 + Mstmp55*M[14] + z*M[36] + M[60];
#pragma omp atomic
Ms[61] += Mstmp115*Mstmp17 + Mstmp115*M[9] + Mstmp121*z + Mstmp138*Mstmp216 + Mstmp18*Mstmp215 + Mstmp20*Mstmp76 + Mstmp20*M[25] + Mstmp220 + Mstmp41*Mstmp55 + Mstmp43*Mstmp54 + Mstmp43*M[20] + Mstmp55*M[15] + Mstmp78*M[4] + z*M[37] + M[61];
#pragma omp atomic
Ms[62] += Mstmp127*M[5] + Mstmp142*Mstmp224 + Mstmp19*Mstmp83 + Mstmp20*Mstmp80 + Mstmp20*M[26] + Mstmp221 + Mstmp222 + Mstmp223 + Mstmp225*M[2] + Mstmp226*Mstmp227 + Mstmp33*Mstmp58 + Mstmp33*M[21] + Mstmp44*Mstmp55 + Mstmp55*M[16] + Mstmp83*M[10] + M[62];
#pragma omp atomic
Ms[63] += Mstmp10*Mstmp127 + Mstmp124*z + Mstmp125*z + Mstmp126*z + Mstmp127*M[6] + Mstmp2*Mstmp225 + Mstmp20*Mstmp85 + Mstmp20*Mstmp86 + Mstmp20*Mstmp89 + Mstmp20*M[27] + Mstmp225*M[3] + Mstmp228 + Mstmp229 + Mstmp230 + Mstmp33*Mstmp61 + Mstmp33*Mstmp62 + Mstmp33*Mstmp63 + Mstmp33*M[22] + Mstmp46*Mstmp55 + Mstmp47*Mstmp55 + Mstmp48*Mstmp55 + Mstmp55*M[17] + z*M[38] + M[63];
#pragma omp atomic
Ms[64] += Mstmp1*Mstmp234 + Mstmp129*z + Mstmp130*z + Mstmp132*z + Mstmp20*Mstmp91 + Mstmp20*Mstmp92 + Mstmp20*Mstmp95 + Mstmp20*M[28] + Mstmp231 + Mstmp232 + Mstmp233 + Mstmp234*M[2] + Mstmp43*Mstmp58 + Mstmp43*Mstmp59 + Mstmp43*Mstmp60 + Mstmp43*M[21] + Mstmp49*Mstmp55 + Mstmp50*Mstmp55 + Mstmp51*Mstmp55 + Mstmp55*M[18] + Mstmp7*Mstmp78 + Mstmp78*M[5] + z*M[39] + M[64];
#pragma omp atomic
Ms[65] += Mstmp135*z + Mstmp159*Mstmp224 + Mstmp19*Mstmp99 + Mstmp20*Mstmp97 + Mstmp20*M[29] + Mstmp234*M[3] + Mstmp235 + Mstmp236*M[0] + Mstmp43*Mstmp61 + Mstmp43*M[22] + Mstmp52*Mstmp55 + Mstmp55*M[19] + Mstmp78*M[6] + Mstmp99*M[10] + z*M[40] + M[65];
#pragma omp atomic
Ms[66] += Mstmp100*Mstmp20 + Mstmp101*Mstmp241 + Mstmp127*M[7] + Mstmp164*Mstmp3 + Mstmp164*M[4] + Mstmp20*M[30] + Mstmp23*Mstmp83 + Mstmp237 + Mstmp238 + Mstmp239 + Mstmp240*M[2] + Mstmp33*Mstmp64 + Mstmp33*M[23] + Mstmp83*M[11] + M[66];
#pragma omp atomic
Ms[67] += Mstmp102*Mstmp20 + Mstmp103*Mstmp20 + Mstmp104*Mstmp20 + Mstmp127*Mstmp15 + Mstmp127*M[8] + Mstmp139*z + Mstmp140*z + Mstmp141*z + Mstmp2*Mstmp240 + Mstmp20*M[31] + Mstmp240*M[3] + Mstmp242 + Mstmp243 + Mstmp244 + Mstmp26*Mstmp83 + Mstmp27*Mstmp83 + Mstmp28*Mstmp83 + Mstmp33*Mstmp69 + Mstmp33*Mstmp71 + Mstmp33*Mstmp73 + Mstmp33*M[24] + Mstmp83*M[12] + z*M[41] + M[67];
#pragma omp atomic
Ms[68] += Mstmp105*Mstmp20 + Mstmp106*Mstmp20 + Mstmp107*Mstmp20 + Mstmp12*Mstmp78 + Mstmp127*Mstmp17 + Mstmp127*M[9] + Mstmp144*z + Mstmp145*z + Mstmp147*z + Mstmp177*Mstmp3 + Mstmp177*M[4] + Mstmp20*M[32] + Mstmp245 + Mstmp246 + Mstmp247 + Mstmp249*Mstmp67 + Mstmp33*Mstmp75 + Mstmp33*Mstmp76 + Mstmp33*Mstmp77 + Mstmp33*M[25] + Mstmp43*Mstmp64 + Mstmp43*Mstmp65 + Mstmp43*Mstmp66 + Mstmp43*M[23] + Mstmp78*M[7] + z*M[42] + M[68];
#pragma omp atomic
Ms[69] += Mstmp1*Mstmp253 + Mstmp109*Mstmp20 + Mstmp110*Mstmp20 + Mstmp111*Mstmp20 + Mstmp14*Mstmp78 + Mstmp150*z + Mstmp151*z + Mstmp153*z + Mstmp20*M[33] + Mstmp23*Mstmp99 + Mstmp24*Mstmp99 + Mstmp25*Mstmp99 + Mstmp250 + Mstmp251 + Mstmp252 + Mstmp253*M[2] + Mstmp43*Mstmp69 + Mstmp43*Mstmp70 + Mstmp43*Mstmp72 + Mstmp43*M[24] + Mstmp78*M[8] + Mstmp99*M[11] + z*M[43] + M[69];
#pragma omp atomic
Ms[70] += Mstmp112*Mstmp20 + Mstmp113*Mstmp241 + Mstmp156*z + Mstmp187*Mstmp3 + Mstmp187*M[4] + Mstmp20*M[34] + Mstmp253*M[3] + Mstmp254 + Mstmp26*Mstmp99 + Mstmp43*Mstmp75 + Mstmp43*M[25] + Mstmp78*M[9] + Mstmp99*M[12] + z*M[44] + M[70];
#pragma omp atomic
Ms[71] += Mstmp0*Mstmp259 + Mstmp164*Mstmp6 + Mstmp164*M[5] + Mstmp189*Mstmp258 + Mstmp255 + Mstmp256 + Mstmp257 + Mstmp29*Mstmp83 + Mstmp33*Mstmp79 + Mstmp33*M[26] + Mstmp83*M[13] + M[71];
#pragma omp atomic
Ms[72] += Mstmp10*Mstmp164 + Mstmp11*Mstmp164 + Mstmp160*z + Mstmp161*z + Mstmp162*z + Mstmp164*Mstmp9 + Mstmp164*M[6] + Mstmp260 + Mstmp261 + Mstmp262 + Mstmp33*Mstmp84 + Mstmp33*Mstmp86 + Mstmp33*Mstmp88 + Mstmp33*M[27] + Mstmp34*Mstmp83 + Mstmp36*Mstmp83 + Mstmp38*Mstmp83 + Mstmp83*M[14] + z*M[45] + M[72];
#pragma omp atomic
Ms[73] += Mstmp0*Mstmp266 + Mstmp165*z + Mstmp166*z + Mstmp168*z + Mstmp177*Mstmp6 + Mstmp177*M[5] + Mstmp263 + Mstmp264 + Mstmp265 + Mstmp266*M[1] + Mstmp33*Mstmp90 + Mstmp33*Mstmp92 + Mstmp33*Mstmp94 + Mstmp33*M[28] + Mstmp40*Mstmp83 + Mstmp41*Mstmp83 + Mstmp42*Mstmp83 + Mstmp43*Mstmp79 + Mstmp43*Mstmp80 + Mstmp43*Mstmp81 + Mstmp43*M[26] + Mstmp83*M[15] + z*M[46] + M[73];
#pragma omp atomic
Ms[74] += Mstmp0*Mstmp270 + Mstmp171*z + Mstmp172*z + Mstmp174*z + Mstmp177*Mstmp9 + Mstmp177*M[6] + Mstmp267 + Mstmp268 + Mstmp269 + Mstmp270*M[1] + Mstmp29*Mstmp99 + Mstmp30*Mstmp99 + Mstmp31*Mstmp99 + Mstmp33*Mstmp96 + Mstmp33*Mstmp97 + Mstmp33*Mstmp98 + Mstmp33*M[29] + Mstmp43*Mstmp84 + Mstmp43*Mstmp85 + Mstmp43*Mstmp87 + Mstmp43*M[27] + Mstmp99*M[13] + z*M[47] + M[74];
#pragma omp atomic
Ms[75] += Mstmp178*z + Mstmp179*z + Mstmp181*z + Mstmp187*Mstmp6 + Mstmp187*Mstmp7 + Mstmp187*Mstmp8 + Mstmp187*M[5] + Mstmp271 + Mstmp272 + Mstmp273 + Mstmp34*Mstmp99 + Mstmp35*Mstmp99 + Mstmp37*Mstmp99 + Mstmp43*Mstmp90 + Mstmp43*Mstmp91 + Mstmp43*Mstmp93 + Mstmp43*M[28] + Mstmp99*M[14] + z*M[48] + M[75];
#pragma omp atomic
Ms[76] += Mstmp0*Mstmp275 + Mstmp184*z + Mstmp187*Mstmp9 + Mstmp187*M[6] + Mstmp203*Mstmp258 + Mstmp274 + Mstmp40*Mstmp99 + Mstmp43*Mstmp96 + Mstmp43*M[29] + Mstmp99*M[15] + z*M[49] + M[76];
#pragma omp atomic
Ms[77] += Mstmp164*M[7] + Mstmp207*Mstmp277 + Mstmp259*M[2] + Mstmp276 + Mstmp33*M[30] + Mstmp83*M[16] + M[77];
#pragma omp atomic
Ms[78] += Mstmp103*Mstmp33 + Mstmp15*Mstmp164 + Mstmp164*M[8] + Mstmp188*z + Mstmp2*Mstmp259 + Mstmp259*M[3] + Mstmp278 + Mstmp33*M[31] + Mstmp47*Mstmp83 + Mstmp83*M[17] + z*M[50] + M[78];
#pragma omp atomic
Ms[79] += Mstmp100*Mstmp43 + Mstmp106*Mstmp33 + Mstmp138*Mstmp280 + Mstmp164*Mstmp17 + Mstmp164*M[9] + Mstmp177*M[7] + Mstmp190*z + Mstmp266*M[2] + Mstmp279 + Mstmp33*M[32] + Mstmp43*M[30] + Mstmp50*Mstmp83 + Mstmp83*M[18] + z*M[51] + M[79];
#pragma omp atomic
Ms[80] += Mstmp102*Mstmp43 + Mstmp110*Mstmp33 + Mstmp177*M[8] + Mstmp193*z + Mstmp226*Mstmp282 + Mstmp266*M[3] + Mstmp270*M[2] + Mstmp281 + Mstmp33*M[33] + Mstmp43*M[31] + Mstmp44*Mstmp99 + Mstmp52*Mstmp83 + Mstmp83*M[19] + Mstmp99*M[16] + z*M[52] + M[80];
#pragma omp atomic
Ms[81] += Mstmp105*Mstmp43 + Mstmp112*Mstmp33 + Mstmp12*Mstmp187 + Mstmp177*M[9] + Mstmp187*M[7] + Mstmp196*z + Mstmp270*M[3] + Mstmp283 + Mstmp284*Mstmp67 + Mstmp33*M[34] + Mstmp43*M[32] + Mstmp46*Mstmp99 + Mstmp99*M[17] + z*M[53] + M[81];
#pragma omp atomic
Ms[82] += Mstmp1*Mstmp275 + Mstmp109*Mstmp43 + Mstmp14*Mstmp187 + Mstmp187*M[8] + Mstmp199*z + Mstmp275*M[2] + Mstmp285 + Mstmp43*M[33] + Mstmp49*Mstmp99 + Mstmp99*M[18] + z*M[54] + M[82];
#pragma omp atomic
Ms[83] += Mstmp187*M[9] + Mstmp207*Mstmp286 + Mstmp275*M[3] + Mstmp43*M[34] + Mstmp99*M[19] + z*M[55] + M[83];
#pragma omp atomic
Ms[84] += Mstmp115*M[10] + Mstmp20*M[35] + Mstmp205*M[4] + Mstmp287*M[1] + Mstmp288*(x*x*x*x*x*x*x) + Mstmp55*M[20] + x*M[56] + M[84];
#pragma omp atomic
Ms[85] += Mstmp1*Mstmp287 + Mstmp115*Mstmp24 + Mstmp115*M[11] + Mstmp119*Mstmp20 + Mstmp20*M[36] + Mstmp204*y + Mstmp205*Mstmp7 + Mstmp205*M[5] + Mstmp287*M[2] + Mstmp55*Mstmp59 + Mstmp55*M[21] + x*M[57] + y*M[56] + M[85];
#pragma omp atomic
Ms[86] += Mstmp10*Mstmp205 + Mstmp115*Mstmp27 + Mstmp115*M[12] + Mstmp122*Mstmp20 + Mstmp2*Mstmp287 + Mstmp20*M[37] + Mstmp204*z + Mstmp205*M[6] + Mstmp287*M[3] + Mstmp55*Mstmp62 + Mstmp55*M[22] + x*M[58] + z*M[56] + M[86];
#pragma omp atomic
Ms[87] += Mstmp114*Mstmp33 + Mstmp115*Mstmp30 + Mstmp115*M[13] + Mstmp12*Mstmp205 + Mstmp125*Mstmp20 + Mstmp127*M[10] + Mstmp13*Mstmp289 + Mstmp20*M[38] + Mstmp205*M[7] + Mstmp208*y + Mstmp225*M[4] + Mstmp290*Mstmp67 + Mstmp33*M[35] + Mstmp55*Mstmp65 + Mstmp55*M[23] + x*M[59] + y*M[57] + M[87];
#pragma omp atomic
Ms[88] += Mstmp115*Mstmp35 + Mstmp115*Mstmp36 + Mstmp115*Mstmp39 + Mstmp115*M[14] + Mstmp130*Mstmp20 + Mstmp131*Mstmp20 + Mstmp134*Mstmp20 + Mstmp14*Mstmp205 + Mstmp15*Mstmp205 + Mstmp16*Mstmp205 + Mstmp20*M[39] + Mstmp205*M[8] + Mstmp208*z + Mstmp209*z + Mstmp210*z + Mstmp211*y + Mstmp55*Mstmp70 + Mstmp55*Mstmp71 + Mstmp55*Mstmp74 + Mstmp55*M[24] + x*M[60] + y*M[58] + z*M[57] + M[88];
#pragma omp atomic
Ms[89] += Mstmp114*Mstmp43 + Mstmp115*Mstmp41 + Mstmp115*M[15] + Mstmp136*Mstmp20 + Mstmp138*Mstmp290 + Mstmp17*Mstmp205 + Mstmp18*Mstmp289 + Mstmp20*M[40] + Mstmp205*M[9] + Mstmp211*z + Mstmp234*M[4] + Mstmp43*M[35] + Mstmp55*Mstmp76 + Mstmp55*M[25] + Mstmp78*M[10] + x*M[61] + z*M[58] + M[89];
#pragma omp atomic
Ms[90] += Mstmp115*Mstmp44 + Mstmp115*M[16] + Mstmp118*Mstmp33 + Mstmp127*M[11] + Mstmp140*Mstmp20 + Mstmp20*M[41] + Mstmp212*y + Mstmp225*M[5] + Mstmp226*Mstmp293 + Mstmp240*M[4] + Mstmp291*M[1] + Mstmp292*M[2] + Mstmp33*M[36] + Mstmp54*Mstmp83 + Mstmp55*Mstmp80 + Mstmp55*M[26] + Mstmp83*M[20] + x*M[62] + y*M[59] + M[90];
#pragma omp atomic
Ms[91] += Mstmp10*Mstmp225 + Mstmp115*Mstmp46 + Mstmp115*Mstmp47 + Mstmp115*Mstmp48 + Mstmp115*M[17] + Mstmp121*Mstmp33 + Mstmp122*Mstmp33 + Mstmp123*Mstmp33 + Mstmp127*Mstmp27 + Mstmp127*M[12] + Mstmp145*Mstmp20 + Mstmp146*Mstmp20 + Mstmp149*Mstmp20 + Mstmp2*Mstmp292 + Mstmp20*M[42] + Mstmp212*z + Mstmp213*z + Mstmp214*z + Mstmp217*y + Mstmp225*M[6] + Mstmp292*M[3] + Mstmp33*M[37] + Mstmp55*Mstmp85 + Mstmp55*Mstmp86 + Mstmp55*Mstmp89 + Mstmp55*M[27] + x*M[63] + y*M[60] + z*M[59] + M[91];
#pragma omp atomic
Ms[92] += Mstmp1*Mstmp294 + Mstmp115*Mstmp49 + Mstmp115*Mstmp50 + Mstmp115*Mstmp51 + Mstmp115*M[18] + Mstmp118*Mstmp43 + Mstmp119*Mstmp43 + Mstmp120*Mstmp43 + Mstmp151*Mstmp20 + Mstmp152*Mstmp20 + Mstmp155*Mstmp20 + Mstmp20*M[43] + Mstmp217*z + Mstmp218*z + Mstmp219*z + Mstmp220*y + Mstmp234*Mstmp7 + Mstmp234*M[5] + Mstmp24*Mstmp78 + Mstmp294*M[2] + Mstmp43*M[36] + Mstmp55*Mstmp91 + Mstmp55*Mstmp92 + Mstmp55*Mstmp95 + Mstmp55*M[28] + Mstmp78*M[11] + x*M[64] + y*M[61] + z*M[60] + M[92];
#pragma omp atomic
Ms[93] += Mstmp115*Mstmp52 + Mstmp115*M[19] + Mstmp121*Mstmp43 + Mstmp157*Mstmp20 + Mstmp20*M[44] + Mstmp220*z + Mstmp234*M[6] + Mstmp236*M[1] + Mstmp253*M[4] + Mstmp293*Mstmp295 + Mstmp294*M[3] + Mstmp43*M[37] + Mstmp54*Mstmp99 + Mstmp55*Mstmp97 + Mstmp55*M[29] + Mstmp78*M[12] + Mstmp99*M[20] + x*M[65] + z*M[61] + M[93];
#pragma omp atomic
Ms[94] += Mstmp100*Mstmp55 + Mstmp124*Mstmp33 + Mstmp127*M[13] + Mstmp161*Mstmp20 + Mstmp164*Mstmp19 + Mstmp164*M[10] + Mstmp20*M[45] + Mstmp221*y + Mstmp224*Mstmp280 + Mstmp225*M[7] + Mstmp240*M[5] + Mstmp291*M[2] + Mstmp296*Mstmp297 + Mstmp33*M[38] + Mstmp55*M[30] + Mstmp58*Mstmp83 + Mstmp83*M[21] + x*M[66] + y*M[62] + M[94];
#pragma omp atomic
Ms[95] += Mstmp10*Mstmp240 + Mstmp102*Mstmp55 + Mstmp103*Mstmp55 + Mstmp104*Mstmp55 + Mstmp127*Mstmp36 + Mstmp127*M[14] + Mstmp129*Mstmp33 + Mstmp131*Mstmp33 + Mstmp133*Mstmp33 + Mstmp15*Mstmp225 + Mstmp166*Mstmp20 + Mstmp167*Mstmp20 + Mstmp170*Mstmp20 + Mstmp2*Mstmp291 + Mstmp20*M[46] + Mstmp221*z + Mstmp222*z + Mstmp223*z + Mstmp225*M[8] + Mstmp228*y + Mstmp240*M[6] + Mstmp291*M[3] + Mstmp33*M[39] + Mstmp55*M[31] + Mstmp61*Mstmp83 + Mstmp62*Mstmp83 + Mstmp63*Mstmp83 + Mstmp83*M[22] + x*M[67] + y*M[63] + z*M[62] + M[95];
#pragma omp atomic
Ms[96] += Mstmp105*Mstmp55 + Mstmp106*Mstmp55 + Mstmp107*Mstmp55 + Mstmp12*Mstmp234 + Mstmp124*Mstmp43 + Mstmp125*Mstmp43 + Mstmp126*Mstmp43 + Mstmp127*Mstmp41 + Mstmp127*M[15] + Mstmp13*Mstmp21*Mstmp298 + Mstmp13*Mstmp224*Mstmp248 + Mstmp135*Mstmp33 + Mstmp136*Mstmp33 + Mstmp137*Mstmp33 + Mstmp17*Mstmp225 + Mstmp172*Mstmp20 + Mstmp173*Mstmp20 + Mstmp176*Mstmp20 + Mstmp177*Mstmp19 + Mstmp177*M[10] + Mstmp20*M[47] + Mstmp225*M[9] + Mstmp228*z + Mstmp229*z + Mstmp230*z + Mstmp231*y + Mstmp234*M[7] + Mstmp30*Mstmp78 + Mstmp33*M[40] + Mstmp43*M[38] + Mstmp55*M[32] + Mstmp78*M[13] + x*M[68] + y*M[64] + z*M[63] + M[96];
#pragma omp atomic
Ms[97] += Mstmp1*Mstmp236 + Mstmp109*Mstmp55 + Mstmp110*Mstmp55 + Mstmp111*Mstmp55 + Mstmp129*Mstmp43 + Mstmp130*Mstmp43 + Mstmp132*Mstmp43 + Mstmp14*Mstmp234 + Mstmp179*Mstmp20 + Mstmp180*Mstmp20 + Mstmp183*Mstmp20 + Mstmp20*M[48] + Mstmp231*z + Mstmp232*z + Mstmp233*z + Mstmp234*M[8] + Mstmp235*y + Mstmp236*M[2] + Mstmp253*Mstmp7 + Mstmp253*M[5] + Mstmp35*Mstmp78 + Mstmp43*M[39] + Mstmp55*M[33] + Mstmp58*Mstmp99 + Mstmp59*Mstmp99 + Mstmp60*Mstmp99 + Mstmp78*M[14] + Mstmp99*M[21] + x*M[69] + y*M[65] + z*M[64] + M[97];
#pragma omp atomic
Ms[98] += Mstmp112*Mstmp55 + Mstmp135*Mstmp43 + Mstmp185*Mstmp20 + Mstmp187*Mstmp19 + Mstmp187*M[10] + Mstmp20*M[49] + Mstmp224*Mstmp284 + Mstmp234*M[9] + Mstmp235*z + Mstmp236*M[3] + Mstmp253*M[6] + Mstmp297*Mstmp299 + Mstmp43*M[40] + Mstmp55*M[34] + Mstmp61*Mstmp99 + Mstmp78*M[15] + Mstmp99*M[22] + x*M[70] + z*M[65] + M[98];
#pragma omp atomic
Ms[99] += Mstmp127*M[16] + Mstmp139*Mstmp33 + Mstmp164*Mstmp23 + Mstmp164*M[11] + Mstmp188*Mstmp20 + Mstmp189*Mstmp301 + Mstmp20*M[50] + Mstmp237*y + Mstmp240*M[7] + Mstmp259*Mstmp3 + Mstmp259*M[4] + Mstmp300*M[2] + Mstmp33*M[41] + Mstmp64*Mstmp83 + Mstmp83*M[23] + x*M[71] + y*M[66] + M[99];
#pragma omp atomic
Ms[100] += Mstmp127*Mstmp47 + Mstmp127*M[17] + Mstmp144*Mstmp33 + Mstmp146*Mstmp33 + Mstmp148*Mstmp33 + Mstmp15*Mstmp240 + Mstmp164*Mstmp26 + Mstmp164*Mstmp27 + Mstmp164*Mstmp28 + Mstmp164*M[12] + Mstmp190*Mstmp20 + Mstmp191*Mstmp20 + Mstmp192*Mstmp20 + Mstmp2*Mstmp300 + Mstmp20*M[51] + Mstmp237*z + Mstmp238*z + Mstmp239*z + Mstmp240*M[8] + Mstmp242*y + Mstmp300*M[3] + Mstmp33*M[42] + Mstmp69*Mstmp83 + Mstmp71*Mstmp83 + Mstmp73*Mstmp83 + Mstmp83*M[24] + x*M[72] + y*M[67] + z*M[66] + M[100];
#pragma omp atomic
Ms[101] += Mstmp127*Mstmp50 + Mstmp127*M[18] + Mstmp139*Mstmp43 + Mstmp140*Mstmp43 + Mstmp141*Mstmp43 + Mstmp150*Mstmp33 + Mstmp152*Mstmp33 + Mstmp154*Mstmp33 + Mstmp17*Mstmp240 + Mstmp177*Mstmp23 + Mstmp177*M[11] + Mstmp193*Mstmp20 + Mstmp194*Mstmp20 + Mstmp195*Mstmp20 + Mstmp20*M[52] + Mstmp240*M[9] + Mstmp242*z + Mstmp243*z + Mstmp244*z + Mstmp245*y + Mstmp266*Mstmp3 + Mstmp266*M[4] + Mstmp298*Mstmp4*Mstmp45 + Mstmp302*M[2] + Mstmp33*M[43] + Mstmp43*M[41] + Mstmp44*Mstmp78 + Mstmp75*Mstmp83 + Mstmp76*Mstmp83 + Mstmp77*Mstmp83 + Mstmp78*M[16] + Mstmp83*M[25] + x*M[73] + y*M[68] + z*M[67] + M[101];
#pragma omp atomic
Ms[102] += Mstmp12*Mstmp253 + Mstmp127*Mstmp52 + Mstmp127*M[19] + Mstmp13*Mstmp4*Mstmp53*Mstmp57 + Mstmp144*Mstmp43 + Mstmp145*Mstmp43 + Mstmp147*Mstmp43 + Mstmp156*Mstmp33 + Mstmp157*Mstmp33 + Mstmp158*Mstmp33 + Mstmp177*Mstmp26 + Mstmp177*M[12] + Mstmp196*Mstmp20 + Mstmp197*Mstmp20 + Mstmp198*Mstmp20 + Mstmp20*M[53] + Mstmp245*z + Mstmp246*z + Mstmp247*z + Mstmp250*y + Mstmp253*M[7] + Mstmp270*Mstmp3 + Mstmp270*M[4] + Mstmp302*M[3] + Mstmp33*M[44] + Mstmp43*M[42] + Mstmp46*Mstmp78 + Mstmp64*Mstmp99 + Mstmp65*Mstmp99 + Mstmp66*Mstmp99 + Mstmp78*M[17] + Mstmp99*M[23] + x*M[74] + y*M[69] + z*M[68] + M[102];
#pragma omp atomic
Ms[103] += Mstmp1*Mstmp303 + Mstmp14*Mstmp253 + Mstmp150*Mstmp43 + Mstmp151*Mstmp43 + Mstmp153*Mstmp43 + Mstmp187*Mstmp23 + Mstmp187*Mstmp24 + Mstmp187*Mstmp25 + Mstmp187*M[11] + Mstmp199*Mstmp20 + Mstmp20*Mstmp200 + Mstmp20*Mstmp201 + Mstmp20*M[54] + Mstmp250*z + Mstmp251*z + Mstmp252*z + Mstmp253*M[8] + Mstmp254*y + Mstmp303*M[2] + Mstmp43*M[43] + Mstmp49*Mstmp78 + Mstmp69*Mstmp99 + Mstmp70*Mstmp99 + Mstmp72*Mstmp99 + Mstmp78*M[18] + Mstmp99*M[24] + x*M[75] + y*M[70] + z*M[69] + M[103];
#pragma omp atomic
Ms[104] += Mstmp156*Mstmp43 + Mstmp187*Mstmp26 + Mstmp187*M[12] + Mstmp20*Mstmp202 + Mstmp20*M[55] + Mstmp203*Mstmp301 + Mstmp253*M[9] + Mstmp254*z + Mstmp275*Mstmp3 + Mstmp275*M[4] + Mstmp303*M[3] + Mstmp43*M[44] + Mstmp75*Mstmp99 + Mstmp78*M[19] + Mstmp99*M[25] + x*M[76] + z*M[70] + M[104];
#pragma omp atomic
Ms[105] += Mstmp0*Mstmp305 + Mstmp160*Mstmp33 + Mstmp164*Mstmp29 + Mstmp164*M[13] + Mstmp255*y + Mstmp259*Mstmp6 + Mstmp259*M[5] + Mstmp277*Mstmp304 + Mstmp33*M[45] + Mstmp79*Mstmp83 + Mstmp83*M[26] + x*M[77] + y*M[71] + M[105];
#pragma omp atomic
Ms[106] += Mstmp10*Mstmp259 + Mstmp11*Mstmp259 + Mstmp164*Mstmp34 + Mstmp164*Mstmp36 + Mstmp164*Mstmp38 + Mstmp164*M[14] + Mstmp165*Mstmp33 + Mstmp167*Mstmp33 + Mstmp169*Mstmp33 + Mstmp255*z + Mstmp256*z + Mstmp257*z + Mstmp259*Mstmp9 + Mstmp259*M[6] + Mstmp260*y + Mstmp33*M[46] + Mstmp83*Mstmp84 + Mstmp83*Mstmp86 + Mstmp83*Mstmp88 + Mstmp83*M[27] + x*M[78] + y*M[72] + z*M[71] + M[106];
#pragma omp atomic
Ms[107] += Mstmp0*Mstmp306 + Mstmp160*Mstmp43 + Mstmp161*Mstmp43 + Mstmp162*Mstmp43 + Mstmp164*Mstmp40 + Mstmp164*Mstmp41 + Mstmp164*Mstmp42 + Mstmp164*M[15] + Mstmp171*Mstmp33 + Mstmp173*Mstmp33 + Mstmp175*Mstmp33 + Mstmp177*Mstmp29 + Mstmp177*M[13] + Mstmp260*z + Mstmp261*z + Mstmp262*z + Mstmp263*y + Mstmp266*Mstmp6 + Mstmp266*M[5] + Mstmp306*M[1] + Mstmp33*M[47] + Mstmp43*M[45] + Mstmp83*Mstmp90 + Mstmp83*Mstmp92 + Mstmp83*Mstmp94 + Mstmp83*M[28] + x*M[79] + y*M[73] + z*M[72] + M[107];
#pragma omp atomic
Ms[108] += Mstmp0*Mstmp307 + Mstmp165*Mstmp43 + Mstmp166*Mstmp43 + Mstmp168*Mstmp43 + Mstmp177*Mstmp34 + Mstmp177*M[14] + Mstmp178*Mstmp33 + Mstmp180*Mstmp33 + Mstmp182*Mstmp33 + Mstmp263*z + Mstmp264*z + Mstmp265*z + Mstmp266*Mstmp9 + Mstmp266*M[6] + Mstmp267*y + Mstmp270*Mstmp6 + Mstmp270*M[5] + Mstmp307*M[1] + Mstmp33*M[48] + Mstmp43*M[46] + Mstmp79*Mstmp99 + Mstmp80*Mstmp99 + Mstmp81*Mstmp99 + Mstmp83*Mstmp96 + Mstmp83*Mstmp97 + Mstmp83*Mstmp98 + Mstmp83*M[29] + Mstmp99*M[26] + x*M[80] + y*M[74] + z*M[73] + M[108];
#pragma omp atomic
Ms[109] += Mstmp0*Mstmp308 + Mstmp171*Mstmp43 + Mstmp172*Mstmp43 + Mstmp174*Mstmp43 + Mstmp177*Mstmp40 + Mstmp177*M[15] + Mstmp184*Mstmp33 + Mstmp185*Mstmp33 + Mstmp186*Mstmp33 + Mstmp187*Mstmp29 + Mstmp187*Mstmp30 + Mstmp187*Mstmp31 + Mstmp187*M[13] + Mstmp267*z + Mstmp268*z + Mstmp269*z + Mstmp270*Mstmp9 + Mstmp270*M[6] + Mstmp271*y + Mstmp308*M[1] + Mstmp33*M[49] + Mstmp43*M[47] + Mstmp84*Mstmp99 + Mstmp85*Mstmp99 + Mstmp87*Mstmp99 + Mstmp99*M[27] + x*M[81] + y*M[75] + z*M[74] + M[109];
#pragma omp atomic
Ms[110] += Mstmp178*Mstmp43 + Mstmp179*Mstmp43 + Mstmp181*Mstmp43 + Mstmp187*Mstmp34 + Mstmp187*Mstmp35 + Mstmp187*Mstmp37 + Mstmp187*M[14] + Mstmp271*z + Mstmp272*z + Mstmp273*z + Mstmp274*y + Mstmp275*Mstmp6 + Mstmp275*Mstmp7 + Mstmp275*Mstmp8 + Mstmp275*M[5] + Mstmp43*M[48] + Mstmp90*Mstmp99 + Mstmp91*Mstmp99 + Mstmp93*Mstmp99 + Mstmp99*M[28] + x*M[82] + y*M[76] + z*M[75] + M[110];
#pragma omp atomic
Ms[111] += Mstmp0*Mstmp309 + Mstmp184*Mstmp43 + Mstmp187*Mstmp40 + Mstmp187*M[15] + Mstmp274*z + Mstmp275*Mstmp9 + Mstmp275*M[6] + Mstmp286*Mstmp304 + Mstmp43*M[49] + Mstmp96*Mstmp99 + Mstmp99*M[29] + x*M[83] + z*M[76] + M[111];
#pragma omp atomic
Ms[112] += Mstmp164*M[16] + Mstmp259*M[7] + Mstmp288*(y*y*y*y*y*y*y) + Mstmp305*M[2] + Mstmp33*M[50] + Mstmp83*M[30] + y*M[77] + M[112];
#pragma omp atomic
Ms[113] += Mstmp103*Mstmp83 + Mstmp15*Mstmp259 + Mstmp164*Mstmp47 + Mstmp164*M[17] + Mstmp191*Mstmp33 + Mstmp2*Mstmp305 + Mstmp259*M[8] + Mstmp276*z + Mstmp305*M[3] + Mstmp33*M[51] + Mstmp83*M[31] + y*M[78] + z*M[77] + M[113];
#pragma omp atomic
Ms[114] += Mstmp106*Mstmp83 + (1.0/240.0)*Mstmp138*Mstmp189 + Mstmp164*Mstmp50 + Mstmp164*M[18] + Mstmp17*Mstmp259 + Mstmp177*M[16] + Mstmp188*Mstmp43 + Mstmp194*Mstmp33 + Mstmp259*M[9] + Mstmp266*M[7] + Mstmp278*z + Mstmp306*M[2] + Mstmp33*M[52] + Mstmp43*M[50] + Mstmp83*M[32] + y*M[79] + z*M[78] + M[114];
#pragma omp atomic
Ms[115] += Mstmp100*Mstmp99 + Mstmp110*Mstmp83 + Mstmp164*Mstmp52 + Mstmp164*M[19] + Mstmp177*M[17] + Mstmp190*Mstmp43 + Mstmp197*Mstmp33 + Mstmp266*M[8] + Mstmp270*M[7] + Mstmp279*z + Mstmp295*Mstmp296 + Mstmp306*M[3] + Mstmp307*M[2] + Mstmp33*M[53] + Mstmp43*M[51] + Mstmp83*M[33] + Mstmp99*M[30] + y*M[80] + z*M[79] + M[115];
#pragma omp atomic
Ms[116] += Mstmp102*Mstmp99 + Mstmp112*Mstmp83 + Mstmp177*M[18] + Mstmp187*Mstmp44 + Mstmp187*M[16] + Mstmp193*Mstmp43 + Mstmp200*Mstmp33 + Mstmp226*Mstmp299 + Mstmp266*M[9] + Mstmp270*M[8] + Mstmp281*z + Mstmp307*M[3] + Mstmp308*M[2] + Mstmp33*M[54] + Mstmp43*M[52] + Mstmp83*M[34] + Mstmp99*M[31] + y*M[81] + z*M[80] + M[116];
#pragma omp atomic
Ms[117] += Mstmp105*Mstmp99 + Mstmp12*Mstmp275 + Mstmp177*M[19] + Mstmp187*Mstmp46 + Mstmp187*M[17] + Mstmp196*Mstmp43 + Mstmp202*Mstmp33 + (1.0/240.0)*Mstmp203*Mstmp67 + Mstmp270*M[9] + Mstmp275*M[7] + Mstmp283*z + Mstmp308*M[3] + Mstmp33*M[55] + Mstmp43*M[53] + Mstmp99*M[32] + y*M[82] + z*M[81] + M[117];
#pragma omp atomic
Ms[118] += Mstmp1*Mstmp309 + Mstmp109*Mstmp99 + Mstmp14*Mstmp275 + Mstmp187*Mstmp49 + Mstmp187*M[18] + Mstmp199*Mstmp43 + Mstmp275*M[8] + Mstmp285*z + Mstmp309*M[2] + Mstmp43*M[54] + Mstmp99*M[33] + y*M[83] + z*M[82] + M[118];
#pragma omp atomic
Ms[119] += Mstmp187*M[19] + Mstmp275*M[9] + Mstmp288*(z*z*z*z*z*z*z) + Mstmp309*M[3] + Mstmp43*M[55] + Mstmp99*M[34] + z*M[83] + M[119];

}

void field_m0_M2L_7(double x, double y, double z, double * M, double * L) {
double R = sqrt(x*x + y*y + z*z);
double D[120];
double Dtmp0 = 1.0*pow(R, -3.0);
double Dtmp1 = -Dtmp0;
double Dtmp2 = (x*x);
double Dtmp3 = pow(R, -5.0);
double Dtmp4 = 3.0*Dtmp3;
double Dtmp5 = x*y;
double Dtmp6 = x*z;
double Dtmp7 = (y*y);
double Dtmp8 = y*z;
double Dtmp9 = 9.0*Dtmp3;
double Dtmp10 = -Dtmp9;
double Dtmp11 = pow(R, -7.0);
double Dtmp12 = 15.0*Dtmp11;
double Dtmp13 = Dtmp12*Dtmp2;
double Dtmp14 = -Dtmp4;
double Dtmp15 = Dtmp13 + Dtmp14;
double Dtmp16 = Dtmp12*Dtmp7;
double Dtmp17 = Dtmp14 + Dtmp16;
double Dtmp18 = 1.0*x;
double Dtmp19 = Dtmp8*x;
double Dtmp20 = (x*x*x*x);
double Dtmp21 = pow(R, -9.0);
double Dtmp22 = 105.0*Dtmp21;
double Dtmp23 = 90.0*Dtmp11;
double Dtmp24 = 45.0*Dtmp11;
double Dtmp25 = -Dtmp24;
double Dtmp26 = Dtmp2*Dtmp22;
double Dtmp27 = x*(Dtmp25 + Dtmp26);
double Dtmp28 = -Dtmp12;
double Dtmp29 = Dtmp22*Dtmp7;
double Dtmp30 = Dtmp25 + Dtmp29;
double Dtmp31 = Dtmp18*y;
double Dtmp32 = Dtmp18*z;
double Dtmp33 = (y*y*y*y);
double Dtmp34 = 225.0*Dtmp11;
double Dtmp35 = pow(R, -11.0);
double Dtmp36 = 945.0*Dtmp35;
double Dtmp37 = Dtmp20*Dtmp36;
double Dtmp38 = Dtmp2*Dtmp21;
double Dtmp39 = 630.0*Dtmp38;
double Dtmp40 = Dtmp24 + Dtmp37 - Dtmp39;
double Dtmp41 = -Dtmp26;
double Dtmp42 = 315.0*Dtmp21;
double Dtmp43 = Dtmp42*Dtmp7;
double Dtmp44 = Dtmp2*Dtmp36;
double Dtmp45 = Dtmp44*Dtmp7;
double Dtmp46 = Dtmp24 + Dtmp45;
double Dtmp47 = -Dtmp42;
double Dtmp48 = -Dtmp29;
double Dtmp49 = Dtmp2*Dtmp42;
double Dtmp50 = Dtmp33*Dtmp36;
double Dtmp51 = Dtmp21*Dtmp7;
double Dtmp52 = 630.0*Dtmp51;
double Dtmp53 = Dtmp24 + Dtmp50 - Dtmp52;
double Dtmp54 = Dtmp36*Dtmp7;
double Dtmp55 = Dtmp18*Dtmp8;
double Dtmp56 = -Dtmp34;
double Dtmp57 = (x*x*x*x*x*x);
double Dtmp58 = pow(R, -13.0);
double Dtmp59 = 10395.0*Dtmp58;
double Dtmp60 = 14175.0*Dtmp35;
double Dtmp61 = 1575.0*Dtmp21;
double Dtmp62 = Dtmp20*Dtmp59;
double Dtmp63 = Dtmp2*Dtmp35;
double Dtmp64 = 9450.0*Dtmp63;
double Dtmp65 = x*(Dtmp61 + Dtmp62 - Dtmp64);
double Dtmp66 = 5670.0*Dtmp63;
double Dtmp67 = Dtmp25 - Dtmp66*Dtmp7;
double Dtmp68 = 945.0*Dtmp21;
double Dtmp69 = 2835.0*Dtmp63;
double Dtmp70 = -Dtmp69;
double Dtmp71 = Dtmp35*Dtmp7;
double Dtmp72 = 2835.0*Dtmp71;
double Dtmp73 = Dtmp2*Dtmp7;
double Dtmp74 = Dtmp59*Dtmp73;
double Dtmp75 = -Dtmp72 + Dtmp74;
double Dtmp76 = Dtmp33*Dtmp59;
double Dtmp77 = 9450.0*Dtmp71;
double Dtmp78 = Dtmp61 + Dtmp76 - Dtmp77;
double Dtmp79 = 5670.0*Dtmp71;
double Dtmp80 = (y*y*y*y*y*y);
double Dtmp81 = -11025.0*Dtmp21;
double Dtmp82 = 135135.0*pow(R, -15.0);
double Dtmp83 = Dtmp57*Dtmp82;
double Dtmp84 = Dtmp20*Dtmp58;
double Dtmp85 = -Dtmp61;
double Dtmp86 = 42525.0*Dtmp63 + Dtmp83 - 155925.0*Dtmp84 + Dtmp85;
double Dtmp87 = Dtmp20*Dtmp82;
double Dtmp88 = Dtmp7*Dtmp87;
double Dtmp89 = -Dtmp62 + Dtmp88;
double Dtmp90 = Dtmp2*Dtmp58;
double Dtmp91 = 103950.0*Dtmp90;
double Dtmp92 = -Dtmp7*Dtmp91 + Dtmp85;
double Dtmp93 = -Dtmp68;
double Dtmp94 = -62370.0*Dtmp7*Dtmp90;
double Dtmp95 = Dtmp72 + Dtmp94;
double Dtmp96 = 31185.0*Dtmp58;
double Dtmp97 = Dtmp33*Dtmp82;
double Dtmp98 = Dtmp2*Dtmp97;
double Dtmp99 = Dtmp69 + Dtmp94 + Dtmp98;
double Dtmp100 = -Dtmp76;
double Dtmp101 = Dtmp80*Dtmp82;
double Dtmp102 = Dtmp33*Dtmp58;
double Dtmp103 = Dtmp101 - 155925.0*Dtmp102 + 42525.0*Dtmp71 + Dtmp85;
D[0] = 1.0/R;
D[1] = -Dtmp0*x;
D[2] = -Dtmp0*y;
D[3] = -Dtmp0*z;
D[4] = Dtmp1 + Dtmp2*Dtmp4;
D[5] = Dtmp4*Dtmp5;
D[6] = Dtmp4*Dtmp6;
D[7] = Dtmp1 + Dtmp4*Dtmp7;
D[8] = Dtmp4*Dtmp8;
D[9] = -D[4] - D[7];
D[10] = -x*(Dtmp10 + Dtmp13);
D[11] = -Dtmp15*y;
D[12] = -Dtmp15*z;
D[13] = -Dtmp17*Dtmp18;
D[14] = -Dtmp12*Dtmp19;
D[15] = -D[10] - D[13];
D[16] = -y*(Dtmp10 + Dtmp16);
D[17] = -Dtmp17*z;
D[18] = -D[11] - D[16];
D[19] = -D[12] - D[17];
D[20] = -Dtmp2*Dtmp23 + Dtmp20*Dtmp22 + Dtmp9;
D[21] = Dtmp27*y;
D[22] = Dtmp27*z;
D[23] = -Dtmp13 - Dtmp16 + Dtmp26*Dtmp7 + Dtmp4;
D[24] = Dtmp8*(Dtmp26 + Dtmp28);
D[25] = -D[20] - D[23];
D[26] = Dtmp30*Dtmp31;
D[27] = Dtmp32*(Dtmp28 + Dtmp29);
D[28] = -D[21] - D[26];
D[29] = -D[22] - D[27];
D[30] = Dtmp22*Dtmp33 - Dtmp23*Dtmp7 + Dtmp9;
D[31] = Dtmp30*Dtmp8;
D[32] = -D[23] - D[30];
D[33] = -D[24] - D[31];
D[34] = -D[25] - D[32];
D[35] = -x*(Dtmp34 + Dtmp37 - 1050.0*Dtmp38);
D[36] = -Dtmp40*y;
D[37] = -Dtmp40*z;
D[38] = -x*(Dtmp41 - Dtmp43 + Dtmp46);
D[39] = -Dtmp19*(Dtmp44 + Dtmp47);
D[40] = -D[35] - D[38];
D[41] = -y*(Dtmp46 + Dtmp48 - Dtmp49);
D[42] = -z*(Dtmp12 + Dtmp41 + Dtmp45 + Dtmp48);
D[43] = -D[36] - D[41];
D[44] = -D[37] - D[42];
D[45] = -Dtmp18*Dtmp53;
D[46] = -Dtmp55*(Dtmp47 + Dtmp54);
D[47] = -D[38] - D[45];
D[48] = -D[39] - D[46];
D[49] = -D[40] - D[47];
D[50] = -y*(Dtmp34 + Dtmp50 - 1050.0*Dtmp51);
D[51] = -Dtmp53*z;
D[52] = -D[41] - D[50];
D[53] = -D[42] - D[51];
D[54] = -D[43] - D[52];
D[55] = -D[44] - D[53];
D[56] = -Dtmp20*Dtmp60 + 4725.0*Dtmp38 + Dtmp56 + Dtmp57*Dtmp59;
D[57] = Dtmp65*y;
D[58] = Dtmp65*z;
D[59] = -Dtmp37 + Dtmp39 + Dtmp43 + Dtmp62*Dtmp7 + Dtmp67;
D[60] = Dtmp8*(Dtmp42 + Dtmp62 - Dtmp66);
D[61] = -D[56] - D[59];
D[62] = Dtmp5*(Dtmp68 + Dtmp70 + Dtmp75);
D[63] = Dtmp6*(Dtmp42 - Dtmp44 + Dtmp75);
D[64] = -D[57] - D[62];
D[65] = -D[58] - D[63];
D[66] = Dtmp2*Dtmp76 + Dtmp49 - Dtmp50 + Dtmp52 + Dtmp67;
D[67] = Dtmp8*(Dtmp42 - Dtmp54 + Dtmp70 + Dtmp74);
D[68] = -D[59] - D[66];
D[69] = -D[60] - D[67];
D[70] = -D[61] - D[68];
D[71] = Dtmp31*Dtmp78;
D[72] = Dtmp32*(Dtmp42 + Dtmp76 - Dtmp79);
D[73] = -D[62] - D[71];
D[74] = -D[63] - D[72];
D[75] = -D[64] - D[73];
D[76] = -D[65] - D[74];
D[77] = -Dtmp33*Dtmp60 + 4725.0*Dtmp51 + Dtmp56 + Dtmp59*Dtmp80;
D[78] = Dtmp78*Dtmp8;
D[79] = -D[66] - D[77];
D[80] = -D[67] - D[78];
D[81] = -D[68] - D[79];
D[82] = -D[69] - D[80];
D[83] = -D[70] - D[81];
D[84] = -x*(99225.0*Dtmp63 + Dtmp81 + Dtmp83 - 218295.0*Dtmp84);
D[85] = -Dtmp86*y;
D[86] = -Dtmp86*z;
D[87] = -x*(Dtmp60*Dtmp7 + Dtmp64 + Dtmp89 + Dtmp92);
D[88] = -Dtmp19*(Dtmp60 + Dtmp87 - Dtmp91);
D[89] = -D[84] - D[87];
D[90] = -y*(17010.0*Dtmp63 - 31185.0*Dtmp84 + Dtmp88 + Dtmp93 + Dtmp95);
D[91] = -z*(Dtmp47 + Dtmp66 + Dtmp89 + Dtmp95);
D[92] = -D[85] - D[90];
D[93] = -D[86] - D[91];
D[94] = -x*(-Dtmp33*Dtmp96 + 17010.0*Dtmp71 + Dtmp93 + Dtmp99);
D[95] = -Dtmp19*(8505.0*Dtmp35 - Dtmp7*Dtmp96 + Dtmp73*Dtmp82 - 31185.0*Dtmp90);
D[96] = -D[87] - D[94];
D[97] = -D[88] - D[95];
D[98] = -D[89] - D[96];
D[99] = -y*(Dtmp100 + Dtmp2*Dtmp60 + Dtmp77 + Dtmp92 + Dtmp98);
D[100] = -z*(Dtmp100 + Dtmp47 + Dtmp79 + Dtmp99);
D[101] = -D[90] - D[99];
D[102] = -D[91] - D[100];
D[103] = -D[92] - D[101];
D[104] = -D[93] - D[102];
D[105] = -Dtmp103*Dtmp18;
D[106] = -Dtmp55*(-103950.0*Dtmp58*Dtmp7 + Dtmp60 + Dtmp97);
D[107] = -D[94] - D[105];
D[108] = -D[95] - D[106];
D[109] = -D[96] - D[107];
D[110] = -D[97] - D[108];
D[111] = -D[98] - D[109];
D[112] = -y*(Dtmp101 - 218295.0*Dtmp102 + 99225.0*Dtmp71 + Dtmp81);
D[113] = -Dtmp103*z;
D[114] = -D[99] - D[112];
D[115] = -D[100] - D[113];
D[116] = -D[101] - D[114];
D[117] = -D[102] - D[115];
D[118] = -D[103] - D[116];
D[119] = -D[104] - D[117];
#pragma omp atomic
L[0] += D[0]*M[0] + D[1]*M[1] + D[2]*M[2] + D[3]*M[3] + D[4]*M[4] + D[5]*M[5] + D[6]*M[6] + D[7]*M[7] + D[8]*M[8] + D[9]*M[9] + D[10]*M[10] + D[11]*M[11] + D[12]*M[12] + D[13]*M[13] + D[14]*M[14] + D[15]*M[15] + D[16]*M[16] + D[17]*M[17] + D[18]*M[18] + D[19]*M[19] + D[20]*M[20] + D[21]*M[21] + D[22]*M[22] + D[23]*M[23] + D[24]*M[24] + D[25]*M[25] + D[26]*M[26] + D[27]*M[27] + D[28]*M[28] + D[29]*M[29] + D[30]*M[30] + D[31]*M[31] + D[32]*M[32] + D[33]*M[33] + D[34]*M[34] + D[35]*M[35] + D[36]*M[36] + D[37]*M[37] + D[38]*M[38] + D[39]*M[39] + D[40]*M[40] + D[41]*M[41] + D[42]*M[42] + D[43]*M[43] + D[44]*M[44] + D[45]*M[45] + D[46]*M[46] + D[47]*M[47] + D[48]*M[48] + D[49]*M[49] + D[50]*M[50] + D[51]*M[51] + D[52]*M[52] + D[53]*M[53] + D[54]*M[54] + D[55]*M[55] + D[56]*M[56] + D[57]*M[57] + D[58]*M[58] + D[59]*M[59] + D[60]*M[60] + D[61]*M[61] + D[62]*M[62] + D[63]*M[63] + D[64]*M[64] + D[65]*M[65] + D[66]*M[66] + D[67]*M[67] + D[68]*M[68] + D[69]*M[69] + D[70]*M[70] + D[71]*M[71] + D[72]*M[72] + D[73]*M[73] + D[74]*M[74] + D[75]*M[75] + D[76]*M[76] + D[77]*M[77] + D[78]*M[78] + D[79]*M[79] + D[80]*M[80] + D[81]*M[81] + D[82]*M[82] + D[83]*M[83] + D[84]*M[84] + D[85]*M[85] + D[86]*M[86] + D[87]*M[87] + D[88]*M[88] + D[89]*M[89] + D[90]*M[90] + D[91]*M[91] + D[92]*M[92] + D[93]*M[93] + D[94]*M[94] + D[95]*M[95] + D[96]*M[96] + D[97]*M[97] + D[98]*M[98] + D[99]*M[99] + D[100]*M[100] + D[101]*M[101] + D[102]*M[102] + D[103]*M[103] + D[104]*M[104] + D[105]*M[105] + D[106]*M[106] + D[107]*M[107] + D[108]*M[108] + D[109]*M[109] + D[110]*M[110] + D[111]*M[111] + D[112]*M[112] + D[113]*M[113] + D[114]*M[114] + D[115]*M[115] + D[116]*M[116] + D[117]*M[117] + D[118]*M[118] + D[119]*M[119];
#pragma omp atomic
L[1] += D[1]*M[0] + D[4]*M[1] + D[5]*M[2] + D[6]*M[3] + D[10]*M[4] + D[11]*M[5] + D[12]*M[6] + D[13]*M[7] + D[14]*M[8] + D[15]*M[9] + D[20]*M[10] + D[21]*M[11] + D[22]*M[12] + D[23]*M[13] + D[24]*M[14] + D[25]*M[15] + D[26]*M[16] + D[27]*M[17] + D[28]*M[18] + D[29]*M[19] + D[35]*M[20] + D[36]*M[21] + D[37]*M[22] + D[38]*M[23] + D[39]*M[24] + D[40]*M[25] + D[41]*M[26] + D[42]*M[27] + D[43]*M[28] + D[44]*M[29] + D[45]*M[30] + D[46]*M[31] + D[47]*M[32] + D[48]*M[33] + D[49]*M[34] + D[56]*M[35] + D[57]*M[36] + D[58]*M[37] + D[59]*M[38] + D[60]*M[39] + D[61]*M[40] + D[62]*M[41] + D[63]*M[42] + D[64]*M[43] + D[65]*M[44] + D[66]*M[45] + D[67]*M[46] + D[68]*M[47] + D[69]*M[48] + D[70]*M[49] + D[71]*M[50] + D[72]*M[51] + D[73]*M[52] + D[74]*M[53] + D[75]*M[54] + D[76]*M[55] + D[84]*M[56] + D[85]*M[57] + D[86]*M[58] + D[87]*M[59] + D[88]*M[60] + D[89]*M[61] + D[90]*M[62] + D[91]*M[63] + D[92]*M[64] + D[93]*M[65] + D[94]*M[66] + D[95]*M[67] + D[96]*M[68] + D[97]*M[69] + D[98]*M[70] + D[99]*M[71] + D[100]*M[72] + D[101]*M[73] + D[102]*M[74] + D[103]*M[75] + D[104]*M[76] + D[105]*M[77] + D[106]*M[78] + D[107]*M[79] + D[108]*M[80] + D[109]*M[81] + D[110]*M[82] + D[111]*M[83];
#pragma omp atomic
L[2] += D[2]*M[0] + D[5]*M[1] + D[7]*M[2] + D[8]*M[3] + D[11]*M[4] + D[13]*M[5] + D[14]*M[6] + D[16]*M[7] + D[17]*M[8] + D[18]*M[9] + D[21]*M[10] + D[23]*M[11] + D[24]*M[12] + D[26]*M[13] + D[27]*M[14] + D[28]*M[15] + D[30]*M[16] + D[31]*M[17] + D[32]*M[18] + D[33]*M[19] + D[36]*M[20] + D[38]*M[21] + D[39]*M[22] + D[41]*M[23] + D[42]*M[24] + D[43]*M[25] + D[45]*M[26] + D[46]*M[27] + D[47]*M[28] + D[48]*M[29] + D[50]*M[30] + D[51]*M[31] + D[52]*M[32] + D[53]*M[33] + D[54]*M[34] + D[57]*M[35] + D[59]*M[36] + D[60]*M[37] + D[62]*M[38] + D[63]*M[39] + D[64]*M[40] + D[66]*M[41] + D[67]*M[42] + D[68]*M[43] + D[69]*M[44] + D[71]*M[45] + D[72]*M[46] + D[73]*M[47] + D[74]*M[48] + D[75]*M[49] + D[77]*M[50] + D[78]*M[51] + D[79]*M[52] + D[80]*M[53] + D[81]*M[54] + D[82]*M[55] + D[85]*M[56] + D[87]*M[57] + D[88]*M[58] + D[90]*M[59] + D[91]*M[60] + D[92]*M[61] + D[94]*M[62] + D[95]*M[63] + D[96]*M[64] + D[97]*M[65] + D[99]*M[66] + D[100]*M[67] + D[101]*M[68] + D[102]*M[69] + D[103]*M[70] + D[105]*M[71] + D[106]*M[72] + D[107]*M[73] + D[108]*M[74] + D[109]*M[75] + D[110]*M[76] + D[112]*M[77] + D[113]*M[78] + D[114]*M[79] + D[115]*M[80] + D[116]*M[81] + D[117]*M[82] + D[118]*M[83];
#pragma omp atomic
L[3] += D[3]*M[0] + D[6]*M[1] + D[8]*M[2] + D[9]*M[3] + D[12]*M[4] + D[14]*M[5] + D[15]*M[6] + D[17]*M[7] + D[18]*M[8] + D[19]*M[9] + D[22]*M[10] + D[24]*M[11] + D[25]*M[12] + D[27]*M[13] + D[28]*M[14] + D[29]*M[15] + D[31]*M[16] + D[32]*M[17] + D[33]*M[18] + D[34]*M[19] + D[37]*M[20] + D[39]*M[21] + D[40]*M[22] + D[42]*M[23] + D[43]*M[24] + D[44]*M[25] + D[46]*M[26] + D[47]*M[27] + D[48]*M[28] + D[49]*M[29] + D[51]*M[30] + D[52]*M[31] + D[53]*M[32] + D[54]*M[33] + D[55]*M[34] + D[58]*M[35] + D[60]*M[36] + D[61]*M[37] + D[63]*M[38] + D[64]*M[39] + D[65]*M[40] + D[67]*M[41] + D[68]*M[42] + D[69]*M[43] + D[70]*M[44] + D[72]*M[45] + D[73]*M[46] + D[74]*M[47] + D[75]*M[48] + D[76]*M[49] + D[78]*M[50] + D[79]*M[51] + D[80]*M[52] + D[81]*M[53] + D[82]*M[54] + D[83]*M[55] + D[86]*M[56] + D[88]*M[57] + D[89]*M[58] + D[91]*M[59] + D[92]*M[60] + D[93]*M[61] + D[95]*M[62] + D[96]*M[63] + D[97]*M[64] + D[98]*M[65] + D[100]*M[66] + D[101]*M[67] + D[102]*M[68] + D[103]*M[69] + D[104]*M[70] + D[106]*M[71] + D[107]*M[72] + D[108]*M[73] + D[109]*M[74] + D[110]*M[75] + D[111]*M[76] + D[113]*M[77] + D[114]*M[78] + D[115]*M[79] + D[116]*M[80] + D[117]*M[81] + D[118]*M[82] + D[119]*M[83];
#pragma omp atomic
L[4] += D[4]*M[0] + D[10]*M[1] + D[11]*M[2] + D[12]*M[3] + D[20]*M[4] + D[21]*M[5] + D[22]*M[6] + D[23]*M[7] + D[24]*M[8] + D[25]*M[9] + D[35]*M[10] + D[36]*M[11] + D[37]*M[12] + D[38]*M[13] + D[39]*M[14] + D[40]*M[15] + D[41]*M[16] + D[42]*M[17] + D[43]*M[18] + D[44]*M[19] + D[56]*M[20] + D[57]*M[21] + D[58]*M[22] + D[59]*M[23] + D[60]*M[24] + D[61]*M[25] + D[62]*M[26] + D[63]*M[27] + D[64]*M[28] + D[65]*M[29] + D[66]*M[30] + D[67]*M[31] + D[68]*M[32] + D[69]*M[33] + D[70]*M[34] + D[84]*M[35] + D[85]*M[36] + D[86]*M[37] + D[87]*M[38] + D[88]*M[39] + D[89]*M[40] + D[90]*M[41] + D[91]*M[42] + D[92]*M[43] + D[93]*M[44] + D[94]*M[45] + D[95]*M[46] + D[96]*M[47] + D[97]*M[48] + D[98]*M[49] + D[99]*M[50] + D[100]*M[51] + D[101]*M[52] + D[102]*M[53] + D[103]*M[54] + D[104]*M[55];
#pragma omp atomic
L[5] += D[5]*M[0] + D[11]*M[1] + D[13]*M[2] + D[14]*M[3] + D[21]*M[4] + D[23]*M[5] + D[24]*M[6] + D[26]*M[7] + D[27]*M[8] + D[28]*M[9] + D[36]*M[10] + D[38]*M[11] + D[39]*M[12] + D[41]*M[13] + D[42]*M[14] + D[43]*M[15] + D[45]*M[16] + D[46]*M[17] + D[47]*M[18] + D[48]*M[19] + D[57]*M[20] + D[59]*M[21] + D[60]*M[22] + D[62]*M[23] + D[63]*M[24] + D[64]*M[25] + D[66]*M[26] + D[67]*M[27] + D[68]*M[28] + D[69]*M[29] + D[71]*M[30] + D[72]*M[31] + D[73]*M[32] + D[74]*M[33] + D[75]*M[34] + D[85]*M[35] + D[87]*M[36] + D[88]*M[37] + D[90]*M[38] + D[91]*M[39] + D[92]*M[40] + D[94]*M[41] + D[95]*M[42] + D[96]*M[43] + D[97]*M[44] + D[99]*M[45] + D[100]*M[46] + D[101]*M[47] + D[102]*M[48] + D[103]*M[49] + D[105]*M[50] + D[106]*M[51] + D[107]*M[52] + D[108]*M[53] + D[109]*M[54] + D[110]*M[55];
#pragma omp atomic
L[6] += D[6]*M[0] + D[12]*M[1] + D[14]*M[2] + D[15]*M[3] + D[22]*M[4] + D[24]*M[5] + D[25]*M[6] + D[27]*M[7] + D[28]*M[8] + D[29]*M[9] + D[37]*M[10] + D[39]*M[11] + D[40]*M[12] + D[42]*M[13] + D[43]*M[14] + D[44]*M[15] + D[46]*M[16] + D[47]*M[17] + D[48]*M[18] + D[49]*M[19] + D[58]*M[20] + D[60]*M[21] + D[61]*M[22] + D[63]*M[23] + D[64]*M[24] + D[65]*M[25] + D[67]*M[26] + D[68]*M[27] + D[69]*M[28] + D[70]*M[29] + D[72]*M[30] + D[73]*M[31] + D[74]*M[32] + D[75]*M[33] + D[76]*M[34] + D[86]*M[35] + D[88]*M[36] + D[89]*M[37] + D[91]*M[38] + D[92]*M[39] + D[93]*M[40] + D[95]*M[41] + D[96]*M[42] + D[97]*M[43] + D[98]*M[44] + D[100]*M[45] + D[101]*M[46] + D[102]*M[47] + D[103]*M[48] + D[104]*M[49] + D[106]*M[50] + D[107]*M[51] + D[108]*M[52] + D[109]*M[53] + D[110]*M[54] + D[111]*M[55];
#pragma omp atomic
L[7] += D[7]*M[0] + D[13]*M[1] + D[16]*M[2] + D[17]*M[3] + D[23]*M[4] + D[26]*M[5] + D[27]*M[6] + D[30]*M[7] + D[31]*M[8] + D[32]*M[9] + D[38]*M[10] + D[41]*M[11] + D[42]*M[12] + D[45]*M[13] + D[46]*M[14] + D[47]*M[15] + D[50]*M[16] + D[51]*M[17] + D[52]*M[18] + D[53]*M[19] + D[59]*M[20] + D[62]*M[21] + D[63]*M[22] + D[66]*M[23] + D[67]*M[24] + D[68]*M[25] + D[71]*M[26] + D[72]*M[27] + D[73]*M[28] + D[74]*M[29] + D[77]*M[30] + D[78]*M[31] + D[79]*M[32] + D[80]*M[33] + D[81]*M[34] + D[87]*M[35] + D[90]*M[36] + D[91]*M[37] + D[94]*M[38] + D[95]*M[39] + D[96]*M[40] + D[99]*M[41] + D[100]*M[42] + D[101]*M[43] + D[102]*M[44] + D[105]*M[45] + D[106]*M[46] + D[107]*M[47] + D[108]*M[48] + D[109]*M[49] + D[112]*M[50] + D[113]*M[51] + D[114]*M[52] + D[115]*M[53] + D[116]*M[54] + D[117]*M[55];
#pragma omp atomic
L[8] += D[8]*M[0] + D[14]*M[1] + D[17]*M[2] + D[18]*M[3] + D[24]*M[4] + D[27]*M[5] + D[28]*M[6] + D[31]*M[7] + D[32]*M[8] + D[33]*M[9] + D[39]*M[10] + D[42]*M[11] + D[43]*M[12] + D[46]*M[13] + D[47]*M[14] + D[48]*M[15] + D[51]*M[16] + D[52]*M[17] + D[53]*M[18] + D[54]*M[19] + D[60]*M[20] + D[63]*M[21] + D[64]*M[22] + D[67]*M[23] + D[68]*M[24] + D[69]*M[25] + D[72]*M[26] + D[73]*M[27] + D[74]*M[28] + D[75]*M[29] + D[78]*M[30] + D[79]*M[31] + D[80]*M[32] + D[81]*M[33] + D[82]*M[34] + D[88]*M[35] + D[91]*M[36] + D[92]*M[37] + D[95]*M[38] + D[96]*M[39] + D[97]*M[40] + D[100]*M[41] + D[101]*M[42] + D[102]*M[43] + D[103]*M[44] + D[106]*M[45] + D[107]*M[46] + D[108]*M[47] + D[109]*M[48] + D[110]*M[49] + D[113]*M[50] + D[114]*M[51] + D[115]*M[52] + D[116]*M[53] + D[117]*M[54] + D[118]*M[55];
#pragma omp atomic
L[9] += D[9]*M[0] + D[15]*M[1] + D[18]*M[2] + D[19]*M[3] + D[25]*M[4] + D[28]*M[5] + D[29]*M[6] + D[32]*M[7] + D[33]*M[8] + D[34]*M[9] + D[40]*M[10] + D[43]*M[11] + D[44]*M[12] + D[47]*M[13] + D[48]*M[14] + D[49]*M[15] + D[52]*M[16] + D[53]*M[17] + D[54]*M[18] + D[55]*M[19] + D[61]*M[20] + D[64]*M[21] + D[65]*M[22] + D[68]*M[23] + D[69]*M[24] + D[70]*M[25] + D[73]*M[26] + D[74]*M[27] + D[75]*M[28] + D[76]*M[29] + D[79]*M[30] + D[80]*M[31] + D[81]*M[32] + D[82]*M[33] + D[83]*M[34] + D[89]*M[35] + D[92]*M[36] + D[93]*M[37] + D[96]*M[38] + D[97]*M[39] + D[98]*M[40] + D[101]*M[41] + D[102]*M[42] + D[103]*M[43] + D[104]*M[44] + D[107]*M[45] + D[108]*M[46] + D[109]*M[47] + D[110]*M[48] + D[111]*M[49] + D[114]*M[50] + D[115]*M[51] + D[116]*M[52] + D[117]*M[53] + D[118]*M[54] + D[119]*M[55];
#pragma omp atomic
L[10] += D[10]*M[0] + D[20]*M[1] + D[21]*M[2] + D[22]*M[3] + D[35]*M[4] + D[36]*M[5] + D[37]*M[6] + D[38]*M[7] + D[39]*M[8] + D[40]*M[9] + D[56]*M[10] + D[57]*M[11] + D[58]*M[12] + D[59]*M[13] + D[60]*M[14] + D[61]*M[15] + D[62]*M[16] + D[63]*M[17] + D[64]*M[18] + D[65]*M[19] + D[84]*M[20] + D[85]*M[21] + D[86]*M[22] + D[87]*M[23] + D[88]*M[24] + D[89]*M[25] + D[90]*M[26] + D[91]*M[27] + D[92]*M[28] + D[93]*M[29] + D[94]*M[30] + D[95]*M[31] + D[96]*M[32] + D[97]*M[33] + D[98]*M[34];
#pragma omp atomic
L[11] += D[11]*M[0] + D[21]*M[1] + D[23]*M[2] + D[24]*M[3] + D[36]*M[4] + D[38]*M[5] + D[39]*M[6] + D[41]*M[7] + D[42]*M[8] + D[43]*M[9] + D[57]*M[10] + D[59]*M[11] + D[60]*M[12] + D[62]*M[13] + D[63]*M[14] + D[64]*M[15] + D[66]*M[16] + D[67]*M[17] + D[68]*M[18] + D[69]*M[19] + D[85]*M[20] + D[87]*M[21] + D[88]*M[22] + D[90]*M[23] + D[91]*M[24] + D[92]*M[25] + D[94]*M[26] + D[95]*M[27] + D[96]*M[28] + D[97]*M[29] + D[99]*M[30] + D[100]*M[31] + D[101]*M[32] + D[102]*M[33] + D[103]*M[34];
#pragma omp atomic
L[12] += D[12]*M[0] + D[22]*M[1] + D[24]*M[2] + D[25]*M[3] + D[37]*M[4] + D[39]*M[5] + D[40]*M[6] + D[42]*M[7] + D[43]*M[8] + D[44]*M[9] + D[58]*M[10] + D[60]*M[11] + D[61]*M[12] + D[63]*M[13] + D[64]*M[14] + D[65]*M[15] + D[67]*M[16] + D[68]*M[17] + D[69]*M[18] + D[70]*M[19] + D[86]*M[20] + D[88]*M[21] + D[89]*M[22] + D[91]*M[23] + D[92]*M[24] + D[93]*M[25] + D[95]*M[26] + D[96]*M[27] + D[97]*M[28] + D[98]*M[29] + D[100]*M[30] + D[101]*M[31] + D[102]*M[32] + D[103]*M[33] + D[104]*M[34];
#pragma omp atomic
L[13] += D[13]*M[0] + D[23]*M[1] + D[26]*M[2] + D[27]*M[3] + D[38]*M[4] + D[41]*M[5] + D[42]*M[6] + D[45]*M[7] + D[46]*M[8] + D[47]*M[9] + D[59]*M[10] + D[62]*M[11] + D[63]*M[12] + D[66]*M[13] + D[67]*M[14] + D[68]*M[15] + D[71]*M[16] + D[72]*M[17] + D[73]*M[18] + D[74]*M[19] + D[87]*M[20] + D[90]*M[21] + D[91]*M[22] + D[94]*M[23] + D[95]*M[24] + D[96]*M[25] + D[99]*M[26] + D[100]*M[27] + D[101]*M[28] + D[102]*M[29] + D[105]*M[30] + D[106]*M[31] + D[107]*M[32] + D[108]*M[33] + D[109]*M[34];
#pragma omp atomic
L[14] += D[14]*M[0] + D[24]*M[1] + D[27]*M[2] + D[28]*M[3] + D[39]*M[4] + D[42]*M[5] + D[43]*M[6] + D[46]*M[7] + D[47]*M[8] + D[48]*M[9] + D[60]*M[10] + D[63]*M[11] + D[64]*M[12] + D[67]*M[13] + D[68]*M[14] + D[69]*M[15] + D[72]*M[16] + D[73]*M[17] + D[74]*M[18] + D[75]*M[19] + D[88]*M[20] + D[91]*M[21] + D[92]*M[22] + D[95]*M[23] + D[96]*M[24] + D[97]*M[25] + D[100]*M[26] + D[101]*M[27] + D[102]*M[28] + D[103]*M[29] + D[106]*M[30] + D[107]*M[31] + D[108]*M[32] + D[109]*M[33] + D[110]*M[34];
#pragma omp atomic
L[15] += D[15]*M[0] + D[25]*M[1] + D[28]*M[2] + D[29]*M[3] + D[40]*M[4] + D[43]*M[5] + D[44]*M[6] + D[47]*M[7] + D[48]*M[8] + D[49]*M[9] + D[61]*M[10] + D[64]*M[11] + D[65]*M[12] + D[68]*M[13] + D[69]*M[14] + D[70]*M[15] + D[73]*M[16] + D[74]*M[17] + D[75]*M[18] + D[76]*M[19] + D[89]*M[20] + D[92]*M[21] + D[93]*M[22] + D[96]*M[23] + D[97]*M[24] + D[98]*M[25] + D[101]*M[26] + D[102]*M[27] + D[103]*M[28] + D[104]*M[29] + D[107]*M[30] + D[108]*M[31] + D[109]*M[32] + D[110]*M[33] + D[111]*M[34];
#pragma omp atomic
L[16] += D[16]*M[0] + D[26]*M[1] + D[30]*M[2] + D[31]*M[3] + D[41]*M[4] + D[45]*M[5] + D[46]*M[6] + D[50]*M[7] + D[51]*M[8] + D[52]*M[9] + D[62]*M[10] + D[66]*M[11] + D[67]*M[12] + D[71]*M[13] + D[72]*M[14] + D[73]*M[15] + D[77]*M[16] + D[78]*M[17] + D[79]*M[18] + D[80]*M[19] + D[90]*M[20] + D[94]*M[21] + D[95]*M[22] + D[99]*M[23] + D[100]*M[24] + D[101]*M[25] + D[105]*M[26] + D[106]*M[27] + D[107]*M[28] + D[108]*M[29] + D[112]*M[30] + D[113]*M[31] + D[114]*M[32] + D[115]*M[33] + D[116]*M[34];
#pragma omp atomic
L[17] += D[17]*M[0] + D[27]*M[1] + D[31]*M[2] + D[32]*M[3] + D[42]*M[4] + D[46]*M[5] + D[47]*M[6] + D[51]*M[7] + D[52]*M[8] + D[53]*M[9] + D[63]*M[10] + D[67]*M[11] + D[68]*M[12] + D[72]*M[13] + D[73]*M[14] + D[74]*M[15] + D[78]*M[16] + D[79]*M[17] + D[80]*M[18] + D[81]*M[19] + D[91]*M[20] + D[95]*M[21] + D[96]*M[22] + D[100]*M[23] + D[101]*M[24] + D[102]*M[25] + D[106]*M[26] + D[107]*M[27] + D[108]*M[28] + D[109]*M[29] + D[113]*M[30] + D[114]*M[31] + D[115]*M[32] + D[116]*M[33] + D[117]*M[34];
#pragma omp atomic
L[18] += D[18]*M[0] + D[28]*M[1] + D[32]*M[2] + D[33]*M[3] + D[43]*M[4] + D[47]*M[5] + D[48]*M[6] + D[52]*M[7] + D[53]*M[8] + D[54]*M[9] + D[64]*M[10] + D[68]*M[11] + D[69]*M[12] + D[73]*M[13] + D[74]*M[14] + D[75]*M[15] + D[79]*M[16] + D[80]*M[17] + D[81]*M[18] + D[82]*M[19] + D[92]*M[20] + D[96]*M[21] + D[97]*M[22] + D[101]*M[23] + D[102]*M[24] + D[103]*M[25] + D[107]*M[26] + D[108]*M[27] + D[109]*M[28] + D[110]*M[29] + D[114]*M[30] + D[115]*M[31] + D[116]*M[32] + D[117]*M[33] + D[118]*M[34];
#pragma omp atomic
L[19] += D[19]*M[0] + D[29]*M[1] + D[33]*M[2] + D[34]*M[3] + D[44]*M[4] + D[48]*M[5] + D[49]*M[6] + D[53]*M[7] + D[54]*M[8] + D[55]*M[9] + D[65]*M[10] + D[69]*M[11] + D[70]*M[12] + D[74]*M[13] + D[75]*M[14] + D[76]*M[15] + D[80]*M[16] + D[81]*M[17] + D[82]*M[18] + D[83]*M[19] + D[93]*M[20] + D[97]*M[21] + D[98]*M[22] + D[102]*M[23] + D[103]*M[24] + D[104]*M[25] + D[108]*M[26] + D[109]*M[27] + D[110]*M[28] + D[111]*M[29] + D[115]*M[30] + D[116]*M[31] + D[117]*M[32] + D[118]*M[33] + D[119]*M[34];
#pragma omp atomic
L[20] += D[20]*M[0] + D[35]*M[1] + D[36]*M[2] + D[37]*M[3] + D[56]*M[4] + D[57]*M[5] + D[58]*M[6] + D[59]*M[7] + D[60]*M[8] + D[61]*M[9] + D[84]*M[10] + D[85]*M[11] + D[86]*M[12] + D[87]*M[13] + D[88]*M[14] + D[89]*M[15] + D[90]*M[16] + D[91]*M[17] + D[92]*M[18] + D[93]*M[19];
#pragma omp atomic
L[21] += D[21]*M[0] + D[36]*M[1] + D[38]*M[2] + D[39]*M[3] + D[57]*M[4] + D[59]*M[5] + D[60]*M[6] + D[62]*M[7] + D[63]*M[8] + D[64]*M[9] + D[85]*M[10] + D[87]*M[11] + D[88]*M[12] + D[90]*M[13] + D[91]*M[14] + D[92]*M[15] + D[94]*M[16] + D[95]*M[17] + D[96]*M[18] + D[97]*M[19];
#pragma omp atomic
L[22] += D[22]*M[0] + D[37]*M[1] + D[39]*M[2] + D[40]*M[3] + D[58]*M[4] + D[60]*M[5] + D[61]*M[6] + D[63]*M[7] + D[64]*M[8] + D[65]*M[9] + D[86]*M[10] + D[88]*M[11] + D[89]*M[12] + D[91]*M[13] + D[92]*M[14] + D[93]*M[15] + D[95]*M[16] + D[96]*M[17] + D[97]*M[18] + D[98]*M[19];
#pragma omp atomic
L[23] += D[23]*M[0] + D[38]*M[1] + D[41]*M[2] + D[42]*M[3] + D[59]*M[4] + D[62]*M[5] + D[63]*M[6] + D[66]*M[7] + D[67]*M[8] + D[68]*M[9] + D[87]*M[10] + D[90]*M[11] + D[91]*M[12] + D[94]*M[13] + D[95]*M[14] + D[96]*M[15] + D[99]*M[16] + D[100]*M[17] + D[101]*M[18] + D[102]*M[19];
#pragma omp atomic
L[24] += D[24]*M[0] + D[39]*M[1] + D[42]*M[2] + D[43]*M[3] + D[60]*M[4] + D[63]*M[5] + D[64]*M[6] + D[67]*M[7] + D[68]*M[8] + D[69]*M[9] + D[88]*M[10] + D[91]*M[11] + D[92]*M[12] + D[95]*M[13] + D[96]*M[14] + D[97]*M[15] + D[100]*M[16] + D[101]*M[17] + D[102]*M[18] + D[103]*M[19];
#pragma omp atomic
L[25] += D[25]*M[0] + D[40]*M[1] + D[43]*M[2] + D[44]*M[3] + D[61]*M[4] + D[64]*M[5] + D[65]*M[6] + D[68]*M[7] + D[69]*M[8] + D[70]*M[9] + D[89]*M[10] + D[92]*M[11] + D[93]*M[12] + D[96]*M[13] + D[97]*M[14] + D[98]*M[15] + D[101]*M[16] + D[102]*M[17] + D[103]*M[18] + D[104]*M[19];
#pragma omp atomic
L[26] += D[26]*M[0] + D[41]*M[1] + D[45]*M[2] + D[46]*M[3] + D[62]*M[4] + D[66]*M[5] + D[67]*M[6] + D[71]*M[7] + D[72]*M[8] + D[73]*M[9] + D[90]*M[10] + D[94]*M[11] + D[95]*M[12] + D[99]*M[13] + D[100]*M[14] + D[101]*M[15] + D[105]*M[16] + D[106]*M[17] + D[107]*M[18] + D[108]*M[19];
#pragma omp atomic
L[27] += D[27]*M[0] + D[42]*M[1] + D[46]*M[2] + D[47]*M[3] + D[63]*M[4] + D[67]*M[5] + D[68]*M[6] + D[72]*M[7] + D[73]*M[8] + D[74]*M[9] + D[91]*M[10] + D[95]*M[11] + D[96]*M[12] + D[100]*M[13] + D[101]*M[14] + D[102]*M[15] + D[106]*M[16] + D[107]*M[17] + D[108]*M[18] + D[109]*M[19];
#pragma omp atomic
L[28] += D[28]*M[0] + D[43]*M[1] + D[47]*M[2] + D[48]*M[3] + D[64]*M[4] + D[68]*M[5] + D[69]*M[6] + D[73]*M[7] + D[74]*M[8] + D[75]*M[9] + D[92]*M[10] + D[96]*M[11] + D[97]*M[12] + D[101]*M[13] + D[102]*M[14] + D[103]*M[15] + D[107]*M[16] + D[108]*M[17] + D[109]*M[18] + D[110]*M[19];
#pragma omp atomic
L[29] += D[29]*M[0] + D[44]*M[1] + D[48]*M[2] + D[49]*M[3] + D[65]*M[4] + D[69]*M[5] + D[70]*M[6] + D[74]*M[7] + D[75]*M[8] + D[76]*M[9] + D[93]*M[10] + D[97]*M[11] + D[98]*M[12] + D[102]*M[13] + D[103]*M[14] + D[104]*M[15] + D[108]*M[16] + D[109]*M[17] + D[110]*M[18] + D[111]*M[19];
#pragma omp atomic
L[30] += D[30]*M[0] + D[45]*M[1] + D[50]*M[2] + D[51]*M[3] + D[66]*M[4] + D[71]*M[5] + D[72]*M[6] + D[77]*M[7] + D[78]*M[8] + D[79]*M[9] + D[94]*M[10] + D[99]*M[11] + D[100]*M[12] + D[105]*M[13] + D[106]*M[14] + D[107]*M[15] + D[112]*M[16] + D[113]*M[17] + D[114]*M[18] + D[115]*M[19];
#pragma omp atomic
L[31] += D[31]*M[0] + D[46]*M[1] + D[51]*M[2] + D[52]*M[3] + D[67]*M[4] + D[72]*M[5] + D[73]*M[6] + D[78]*M[7] + D[79]*M[8] + D[80]*M[9] + D[95]*M[10] + D[100]*M[11] + D[101]*M[12] + D[106]*M[13] + D[107]*M[14] + D[108]*M[15] + D[113]*M[16] + D[114]*M[17] + D[115]*M[18] + D[116]*M[19];
#pragma omp atomic
L[32] += D[32]*M[0] + D[47]*M[1] + D[52]*M[2] + D[53]*M[3] + D[68]*M[4] + D[73]*M[5] + D[74]*M[6] + D[79]*M[7] + D[80]*M[8] + D[81]*M[9] + D[96]*M[10] + D[101]*M[11] + D[102]*M[12] + D[107]*M[13] + D[108]*M[14] + D[109]*M[15] + D[114]*M[16] + D[115]*M[17] + D[116]*M[18] + D[117]*M[19];
#pragma omp atomic
L[33] += D[33]*M[0] + D[48]*M[1] + D[53]*M[2] + D[54]*M[3] + D[69]*M[4] + D[74]*M[5] + D[75]*M[6] + D[80]*M[7] + D[81]*M[8] + D[82]*M[9] + D[97]*M[10] + D[102]*M[11] + D[103]*M[12] + D[108]*M[13] + D[109]*M[14] + D[110]*M[15] + D[115]*M[16] + D[116]*M[17] + D[117]*M[18] + D[118]*M[19];
#pragma omp atomic
L[34] += D[34]*M[0] + D[49]*M[1] + D[54]*M[2] + D[55]*M[3] + D[70]*M[4] + D[75]*M[5] + D[76]*M[6] + D[81]*M[7] + D[82]*M[8] + D[83]*M[9] + D[98]*M[10] + D[103]*M[11] + D[104]*M[12] + D[109]*M[13] + D[110]*M[14] + D[111]*M[15] + D[116]*M[16] + D[117]*M[17] + D[118]*M[18] + D[119]*M[19];
#pragma omp atomic
L[35] += D[35]*M[0] + D[56]*M[1] + D[57]*M[2] + D[58]*M[3] + D[84]*M[4] + D[85]*M[5] + D[86]*M[6] + D[87]*M[7] + D[88]*M[8] + D[89]*M[9];
#pragma omp atomic
L[36] += D[36]*M[0] + D[57]*M[1] + D[59]*M[2] + D[60]*M[3] + D[85]*M[4] + D[87]*M[5] + D[88]*M[6] + D[90]*M[7] + D[91]*M[8] + D[92]*M[9];
#pragma omp atomic
L[37] += D[37]*M[0] + D[58]*M[1] + D[60]*M[2] + D[61]*M[3] + D[86]*M[4] + D[88]*M[5] + D[89]*M[6] + D[91]*M[7] + D[92]*M[8] + D[93]*M[9];
#pragma omp atomic
L[38] += D[38]*M[0] + D[59]*M[1] + D[62]*M[2] + D[63]*M[3] + D[87]*M[4] + D[90]*M[5] + D[91]*M[6] + D[94]*M[7] + D[95]*M[8] + D[96]*M[9];
#pragma omp atomic
L[39] += D[39]*M[0] + D[60]*M[1] + D[63]*M[2] + D[64]*M[3] + D[88]*M[4] + D[91]*M[5] + D[92]*M[6] + D[95]*M[7] + D[96]*M[8] + D[97]*M[9];
#pragma omp atomic
L[40] += D[40]*M[0] + D[61]*M[1] + D[64]*M[2] + D[65]*M[3] + D[89]*M[4] + D[92]*M[5] + D[93]*M[6] + D[96]*M[7] + D[97]*M[8] + D[98]*M[9];
#pragma omp atomic
L[41] += D[41]*M[0] + D[62]*M[1] + D[66]*M[2] + D[67]*M[3] + D[90]*M[4] + D[94]*M[5] + D[95]*M[6] + D[99]*M[7] + D[100]*M[8] + D[101]*M[9];
#pragma omp atomic
L[42] += D[42]*M[0] + D[63]*M[1] + D[67]*M[2] + D[68]*M[3] + D[91]*M[4] + D[95]*M[5] + D[96]*M[6] + D[100]*M[7] + D[101]*M[8] + D[102]*M[9];
#pragma omp atomic
L[43] += D[43]*M[0] + D[64]*M[1] + D[68]*M[2] + D[69]*M[3] + D[92]*M[4] + D[96]*M[5] + D[97]*M[6] + D[101]*M[7] + D[102]*M[8] + D[103]*M[9];
#pragma omp atomic
L[44] += D[44]*M[0] + D[65]*M[1] + D[69]*M[2] + D[70]*M[3] + D[93]*M[4] + D[97]*M[5] + D[98]*M[6] + D[102]*M[7] + D[103]*M[8] + D[104]*M[9];
#pragma omp atomic
L[45] += D[45]*M[0] + D[66]*M[1] + D[71]*M[2] + D[72]*M[3] + D[94]*M[4] + D[99]*M[5] + D[100]*M[6] + D[105]*M[7] + D[106]*M[8] + D[107]*M[9];
#pragma omp atomic
L[46] += D[46]*M[0] + D[67]*M[1] + D[72]*M[2] + D[73]*M[3] + D[95]*M[4] + D[100]*M[5] + D[101]*M[6] + D[106]*M[7] + D[107]*M[8] + D[108]*M[9];
#pragma omp atomic
L[47] += D[47]*M[0] + D[68]*M[1] + D[73]*M[2] + D[74]*M[3] + D[96]*M[4] + D[101]*M[5] + D[102]*M[6] + D[107]*M[7] + D[108]*M[8] + D[109]*M[9];
#pragma omp atomic
L[48] += D[48]*M[0] + D[69]*M[1] + D[74]*M[2] + D[75]*M[3] + D[97]*M[4] + D[102]*M[5] + D[103]*M[6] + D[108]*M[7] + D[109]*M[8] + D[110]*M[9];
#pragma omp atomic
L[49] += D[49]*M[0] + D[70]*M[1] + D[75]*M[2] + D[76]*M[3] + D[98]*M[4] + D[103]*M[5] + D[104]*M[6] + D[109]*M[7] + D[110]*M[8] + D[111]*M[9];
#pragma omp atomic
L[50] += D[50]*M[0] + D[71]*M[1] + D[77]*M[2] + D[78]*M[3] + D[99]*M[4] + D[105]*M[5] + D[106]*M[6] + D[112]*M[7] + D[113]*M[8] + D[114]*M[9];
#pragma omp atomic
L[51] += D[51]*M[0] + D[72]*M[1] + D[78]*M[2] + D[79]*M[3] + D[100]*M[4] + D[106]*M[5] + D[107]*M[6] + D[113]*M[7] + D[114]*M[8] + D[115]*M[9];
#pragma omp atomic
L[52] += D[52]*M[0] + D[73]*M[1] + D[79]*M[2] + D[80]*M[3] + D[101]*M[4] + D[107]*M[5] + D[108]*M[6] + D[114]*M[7] + D[115]*M[8] + D[116]*M[9];
#pragma omp atomic
L[53] += D[53]*M[0] + D[74]*M[1] + D[80]*M[2] + D[81]*M[3] + D[102]*M[4] + D[108]*M[5] + D[109]*M[6] + D[115]*M[7] + D[116]*M[8] + D[117]*M[9];
#pragma omp atomic
L[54] += D[54]*M[0] + D[75]*M[1] + D[81]*M[2] + D[82]*M[3] + D[103]*M[4] + D[109]*M[5] + D[110]*M[6] + D[116]*M[7] + D[117]*M[8] + D[118]*M[9];
#pragma omp atomic
L[55] += D[55]*M[0] + D[76]*M[1] + D[82]*M[2] + D[83]*M[3] + D[104]*M[4] + D[110]*M[5] + D[111]*M[6] + D[117]*M[7] + D[118]*M[8] + D[119]*M[9];
#pragma omp atomic
L[56] += D[56]*M[0] + D[84]*M[1] + D[85]*M[2] + D[86]*M[3];
#pragma omp atomic
L[57] += D[57]*M[0] + D[85]*M[1] + D[87]*M[2] + D[88]*M[3];
#pragma omp atomic
L[58] += D[58]*M[0] + D[86]*M[1] + D[88]*M[2] + D[89]*M[3];
#pragma omp atomic
L[59] += D[59]*M[0] + D[87]*M[1] + D[90]*M[2] + D[91]*M[3];
#pragma omp atomic
L[60] += D[60]*M[0] + D[88]*M[1] + D[91]*M[2] + D[92]*M[3];
#pragma omp atomic
L[61] += D[61]*M[0] + D[89]*M[1] + D[92]*M[2] + D[93]*M[3];
#pragma omp atomic
L[62] += D[62]*M[0] + D[90]*M[1] + D[94]*M[2] + D[95]*M[3];
#pragma omp atomic
L[63] += D[63]*M[0] + D[91]*M[1] + D[95]*M[2] + D[96]*M[3];
#pragma omp atomic
L[64] += D[64]*M[0] + D[92]*M[1] + D[96]*M[2] + D[97]*M[3];
#pragma omp atomic
L[65] += D[65]*M[0] + D[93]*M[1] + D[97]*M[2] + D[98]*M[3];
#pragma omp atomic
L[66] += D[66]*M[0] + D[94]*M[1] + D[99]*M[2] + D[100]*M[3];
#pragma omp atomic
L[67] += D[67]*M[0] + D[95]*M[1] + D[100]*M[2] + D[101]*M[3];
#pragma omp atomic
L[68] += D[68]*M[0] + D[96]*M[1] + D[101]*M[2] + D[102]*M[3];
#pragma omp atomic
L[69] += D[69]*M[0] + D[97]*M[1] + D[102]*M[2] + D[103]*M[3];
#pragma omp atomic
L[70] += D[70]*M[0] + D[98]*M[1] + D[103]*M[2] + D[104]*M[3];
#pragma omp atomic
L[71] += D[71]*M[0] + D[99]*M[1] + D[105]*M[2] + D[106]*M[3];
#pragma omp atomic
L[72] += D[72]*M[0] + D[100]*M[1] + D[106]*M[2] + D[107]*M[3];
#pragma omp atomic
L[73] += D[73]*M[0] + D[101]*M[1] + D[107]*M[2] + D[108]*M[3];
#pragma omp atomic
L[74] += D[74]*M[0] + D[102]*M[1] + D[108]*M[2] + D[109]*M[3];
#pragma omp atomic
L[75] += D[75]*M[0] + D[103]*M[1] + D[109]*M[2] + D[110]*M[3];
#pragma omp atomic
L[76] += D[76]*M[0] + D[104]*M[1] + D[110]*M[2] + D[111]*M[3];
#pragma omp atomic
L[77] += D[77]*M[0] + D[105]*M[1] + D[112]*M[2] + D[113]*M[3];
#pragma omp atomic
L[78] += D[78]*M[0] + D[106]*M[1] + D[113]*M[2] + D[114]*M[3];
#pragma omp atomic
L[79] += D[79]*M[0] + D[107]*M[1] + D[114]*M[2] + D[115]*M[3];
#pragma omp atomic
L[80] += D[80]*M[0] + D[108]*M[1] + D[115]*M[2] + D[116]*M[3];
#pragma omp atomic
L[81] += D[81]*M[0] + D[109]*M[1] + D[116]*M[2] + D[117]*M[3];
#pragma omp atomic
L[82] += D[82]*M[0] + D[110]*M[1] + D[117]*M[2] + D[118]*M[3];
#pragma omp atomic
L[83] += D[83]*M[0] + D[111]*M[1] + D[118]*M[2] + D[119]*M[3];
#pragma omp atomic
L[84] += D[84]*M[0];
#pragma omp atomic
L[85] += D[85]*M[0];
#pragma omp atomic
L[86] += D[86]*M[0];
#pragma omp atomic
L[87] += D[87]*M[0];
#pragma omp atomic
L[88] += D[88]*M[0];
#pragma omp atomic
L[89] += D[89]*M[0];
#pragma omp atomic
L[90] += D[90]*M[0];
#pragma omp atomic
L[91] += D[91]*M[0];
#pragma omp atomic
L[92] += D[92]*M[0];
#pragma omp atomic
L[93] += D[93]*M[0];
#pragma omp atomic
L[94] += D[94]*M[0];
#pragma omp atomic
L[95] += D[95]*M[0];
#pragma omp atomic
L[96] += D[96]*M[0];
#pragma omp atomic
L[97] += D[97]*M[0];
#pragma omp atomic
L[98] += D[98]*M[0];
#pragma omp atomic
L[99] += D[99]*M[0];
#pragma omp atomic
L[100] += D[100]*M[0];
#pragma omp atomic
L[101] += D[101]*M[0];
#pragma omp atomic
L[102] += D[102]*M[0];
#pragma omp atomic
L[103] += D[103]*M[0];
#pragma omp atomic
L[104] += D[104]*M[0];
#pragma omp atomic
L[105] += D[105]*M[0];
#pragma omp atomic
L[106] += D[106]*M[0];
#pragma omp atomic
L[107] += D[107]*M[0];
#pragma omp atomic
L[108] += D[108]*M[0];
#pragma omp atomic
L[109] += D[109]*M[0];
#pragma omp atomic
L[110] += D[110]*M[0];
#pragma omp atomic
L[111] += D[111]*M[0];
#pragma omp atomic
L[112] += D[112]*M[0];
#pragma omp atomic
L[113] += D[113]*M[0];
#pragma omp atomic
L[114] += D[114]*M[0];
#pragma omp atomic
L[115] += D[115]*M[0];
#pragma omp atomic
L[116] += D[116]*M[0];
#pragma omp atomic
L[117] += D[117]*M[0];
#pragma omp atomic
L[118] += D[118]*M[0];
#pragma omp atomic
L[119] += D[119]*M[0];

}

void field_m0_L2L_7(double x, double y, double z, double * L, double * Ls) {
double Lstmp0 = y*L[5];
double Lstmp1 = z*L[6];
double Lstmp2 = z*L[8];
double Lstmp3 = z*L[14];
double Lstmp4 = Lstmp3*y;
double Lstmp5 = (x*x);
double Lstmp6 = (1.0/2.0)*Lstmp5;
double Lstmp7 = (x*x*x);
double Lstmp8 = (1.0/6.0)*Lstmp7;
double Lstmp9 = (x*x*x*x);
double Lstmp10 = (1.0/24.0)*Lstmp9;
double Lstmp11 = (x*x*x*x*x);
double Lstmp12 = (1.0/120.0)*Lstmp11;
double Lstmp13 = (1.0/720.0)*(x*x*x*x*x*x);
double Lstmp14 = (y*y);
double Lstmp15 = (1.0/2.0)*Lstmp14;
double Lstmp16 = (y*y*y);
double Lstmp17 = (1.0/6.0)*Lstmp16;
double Lstmp18 = (y*y*y*y);
double Lstmp19 = (1.0/24.0)*Lstmp18;
double Lstmp20 = (y*y*y*y*y);
double Lstmp21 = (1.0/120.0)*Lstmp20;
double Lstmp22 = (1.0/720.0)*(y*y*y*y*y*y);
double Lstmp23 = (z*z);
double Lstmp24 = (1.0/2.0)*Lstmp23;
double Lstmp25 = (z*z*z);
double Lstmp26 = (1.0/6.0)*Lstmp25;
double Lstmp27 = (z*z*z*z);
double Lstmp28 = (1.0/24.0)*Lstmp27;
double Lstmp29 = (z*z*z*z*z);
double Lstmp30 = (1.0/120.0)*Lstmp29;
double Lstmp31 = (1.0/720.0)*(z*z*z*z*z*z);
double Lstmp32 = x*L[13];
double Lstmp33 = x*L[26];
double Lstmp34 = x*L[45];
double Lstmp35 = x*L[71];
double Lstmp36 = x*L[105];
double Lstmp37 = x*L[15];
double Lstmp38 = x*L[29];
double Lstmp39 = x*L[49];
double Lstmp40 = x*L[76];
double Lstmp41 = x*L[111];
double Lstmp42 = y*L[11];
double Lstmp43 = z*L[12];
double Lstmp44 = y*L[21];
double Lstmp45 = z*L[22];
double Lstmp46 = y*L[36];
double Lstmp47 = z*L[37];
double Lstmp48 = y*L[57];
double Lstmp49 = z*L[58];
double Lstmp50 = y*L[85];
double Lstmp51 = z*L[86];
double Lstmp52 = y*L[18];
double Lstmp53 = y*L[33];
double Lstmp54 = y*L[54];
double Lstmp55 = y*L[82];
double Lstmp56 = y*L[118];
double Lstmp57 = z*L[17];
double Lstmp58 = z*L[31];
double Lstmp59 = z*L[51];
double Lstmp60 = z*L[78];
double Lstmp61 = z*L[113];
double Lstmp62 = y*L[28];
double Lstmp63 = Lstmp62*x;
double Lstmp64 = y*L[48];
double Lstmp65 = Lstmp64*x;
double Lstmp66 = y*L[75];
double Lstmp67 = Lstmp66*x;
double Lstmp68 = y*L[110];
double Lstmp69 = Lstmp68*x;
double Lstmp70 = z*L[27];
double Lstmp71 = Lstmp70*x;
double Lstmp72 = z*L[46];
double Lstmp73 = Lstmp72*x;
double Lstmp74 = z*L[72];
double Lstmp75 = Lstmp74*x;
double Lstmp76 = z*L[106];
double Lstmp77 = Lstmp76*x;
double Lstmp78 = z*L[24];
double Lstmp79 = Lstmp78*y;
double Lstmp80 = z*L[39];
double Lstmp81 = Lstmp80*y;
double Lstmp82 = z*L[60];
double Lstmp83 = Lstmp82*y;
double Lstmp84 = z*L[88];
double Lstmp85 = Lstmp84*y;
double Lstmp86 = (1.0/4.0)*Lstmp5;
double Lstmp87 = Lstmp14*Lstmp86;
double Lstmp88 = (1.0/12.0)*Lstmp5;
double Lstmp89 = Lstmp16*Lstmp88;
double Lstmp90 = (1.0/48.0)*Lstmp5;
double Lstmp91 = Lstmp18*Lstmp90;
double Lstmp92 = (1.0/240.0)*Lstmp5;
double Lstmp93 = Lstmp23*Lstmp86;
double Lstmp94 = Lstmp25*Lstmp88;
double Lstmp95 = Lstmp27*Lstmp90;
double Lstmp96 = (1.0/12.0)*Lstmp7;
double Lstmp97 = Lstmp14*Lstmp96;
double Lstmp98 = (1.0/36.0)*Lstmp7;
double Lstmp99 = Lstmp16*Lstmp98;
double Lstmp100 = (1.0/144.0)*Lstmp7;
double Lstmp101 = Lstmp23*Lstmp96;
double Lstmp102 = Lstmp25*Lstmp98;
double Lstmp103 = (1.0/48.0)*Lstmp9;
double Lstmp104 = Lstmp103*Lstmp14;
double Lstmp105 = (1.0/144.0)*Lstmp9;
double Lstmp106 = Lstmp103*Lstmp23;
double Lstmp107 = (1.0/240.0)*Lstmp11;
double Lstmp108 = Lstmp14*Lstmp23;
double Lstmp109 = (1.0/4.0)*Lstmp108;
double Lstmp110 = Lstmp14*Lstmp25;
double Lstmp111 = (1.0/12.0)*Lstmp110;
double Lstmp112 = (1.0/48.0)*Lstmp14*Lstmp27;
double Lstmp113 = Lstmp16*Lstmp23;
double Lstmp114 = (1.0/12.0)*Lstmp113;
double Lstmp115 = (1.0/36.0)*Lstmp16*Lstmp25;
double Lstmp116 = (1.0/48.0)*Lstmp18*Lstmp23;
double Lstmp117 = x*L[47];
double Lstmp118 = x*L[74];
double Lstmp119 = x*L[109];
double Lstmp120 = x*L[73];
double Lstmp121 = x*L[108];
double Lstmp122 = x*L[107];
double Lstmp123 = y*L[43];
double Lstmp124 = y*L[69];
double Lstmp125 = y*L[103];
double Lstmp126 = z*L[42];
double Lstmp127 = z*L[67];
double Lstmp128 = z*L[100];
double Lstmp129 = y*L[64];
double Lstmp130 = y*L[97];
double Lstmp131 = z*L[63];
double Lstmp132 = z*L[95];
double Lstmp133 = y*L[92];
double Lstmp134 = z*L[91];
double Lstmp135 = (1.0/8.0)*Lstmp108*Lstmp5;
double Lstmp136 = (1.0/24.0)*Lstmp5;
double Lstmp137 = x*L[23];
double Lstmp138 = x*L[41];
double Lstmp139 = x*L[66];
double Lstmp140 = x*L[99];
double Lstmp141 = x*L[25];
double Lstmp142 = x*L[44];
double Lstmp143 = x*L[70];
double Lstmp144 = x*L[104];
double Lstmp145 = Lstmp123*x;
double Lstmp146 = Lstmp124*x;
double Lstmp147 = Lstmp125*x;
double Lstmp148 = Lstmp126*x;
double Lstmp149 = Lstmp127*x;
double Lstmp150 = Lstmp128*x;
double Lstmp151 = x*L[68];
double Lstmp152 = x*L[102];
double Lstmp153 = x*L[101];
double Lstmp154 = y*L[13];
double Lstmp155 = Lstmp70*y;
double Lstmp156 = x*L[28];
double Lstmp157 = x*L[48];
double Lstmp158 = x*L[75];
double Lstmp159 = x*L[110];
double Lstmp160 = y*L[23];
double Lstmp161 = y*L[38];
double Lstmp162 = y*L[59];
double Lstmp163 = y*L[87];
double Lstmp164 = y*L[32];
double Lstmp165 = y*L[53];
double Lstmp166 = y*L[81];
double Lstmp167 = y*L[117];
double Lstmp168 = y*L[47];
double Lstmp169 = Lstmp168*x;
double Lstmp170 = y*L[74];
double Lstmp171 = Lstmp170*x;
double Lstmp172 = y*L[109];
double Lstmp173 = Lstmp172*x;
double Lstmp174 = Lstmp126*y;
double Lstmp175 = Lstmp131*y;
double Lstmp176 = Lstmp134*y;
double Lstmp177 = y*L[68];
double Lstmp178 = y*L[102];
double Lstmp179 = y*L[96];
double Lstmp180 = y*L[14];
double Lstmp181 = z*L[15];
double Lstmp182 = z*L[18];
double Lstmp183 = z*L[28];
double Lstmp184 = Lstmp183*y;
double Lstmp185 = x*L[27];
double Lstmp186 = x*L[46];
double Lstmp187 = x*L[72];
double Lstmp188 = x*L[106];
double Lstmp189 = y*L[24];
double Lstmp190 = z*L[25];
double Lstmp191 = y*L[39];
double Lstmp192 = z*L[40];
double Lstmp193 = y*L[60];
double Lstmp194 = z*L[61];
double Lstmp195 = y*L[88];
double Lstmp196 = z*L[89];
double Lstmp197 = z*L[32];
double Lstmp198 = z*L[52];
double Lstmp199 = z*L[79];
double Lstmp200 = z*L[114];
double Lstmp201 = z*L[47];
double Lstmp202 = Lstmp201*x;
double Lstmp203 = z*L[73];
double Lstmp204 = Lstmp203*x;
double Lstmp205 = z*L[107];
double Lstmp206 = Lstmp205*x;
double Lstmp207 = z*L[43];
double Lstmp208 = Lstmp207*y;
double Lstmp209 = z*L[64];
double Lstmp210 = Lstmp209*y;
double Lstmp211 = z*L[92];
double Lstmp212 = Lstmp211*y;
double Lstmp213 = z*L[68];
double Lstmp214 = z*L[101];
double Lstmp215 = z*L[96];
double Lstmp216 = x*L[38];
double Lstmp217 = x*L[62];
double Lstmp218 = x*L[94];
double Lstmp219 = x*L[40];
double Lstmp220 = x*L[65];
double Lstmp221 = x*L[98];
double Lstmp222 = Lstmp129*x;
double Lstmp223 = Lstmp130*x;
double Lstmp224 = Lstmp131*x;
double Lstmp225 = Lstmp132*x;
double Lstmp226 = x*L[96];
double Lstmp227 = x*L[43];
double Lstmp228 = x*L[69];
double Lstmp229 = x*L[103];
double Lstmp230 = Lstmp177*x;
double Lstmp231 = Lstmp178*x;
double Lstmp232 = x*L[42];
double Lstmp233 = x*L[67];
double Lstmp234 = x*L[100];
double Lstmp235 = Lstmp213*x;
double Lstmp236 = Lstmp214*x;
double Lstmp237 = y*L[26];
double Lstmp238 = Lstmp72*y;
double Lstmp239 = y*L[41];
double Lstmp240 = y*L[62];
double Lstmp241 = y*L[90];
double Lstmp242 = y*L[52];
double Lstmp243 = y*L[80];
double Lstmp244 = y*L[116];
double Lstmp245 = y*L[73];
double Lstmp246 = Lstmp245*x;
double Lstmp247 = y*L[108];
double Lstmp248 = Lstmp247*x;
double Lstmp249 = Lstmp127*y;
double Lstmp250 = Lstmp132*y;
double Lstmp251 = y*L[101];
double Lstmp252 = y*L[27];
double Lstmp253 = Lstmp201*y;
double Lstmp254 = y*L[42];
double Lstmp255 = y*L[63];
double Lstmp256 = y*L[91];
double Lstmp257 = Lstmp213*y;
double Lstmp258 = Lstmp215*y;
double Lstmp259 = z*L[29];
double Lstmp260 = z*L[33];
double Lstmp261 = z*L[48];
double Lstmp262 = Lstmp261*y;
double Lstmp263 = z*L[44];
double Lstmp264 = z*L[65];
double Lstmp265 = z*L[93];
double Lstmp266 = z*L[53];
double Lstmp267 = z*L[80];
double Lstmp268 = z*L[115];
double Lstmp269 = z*L[74];
double Lstmp270 = Lstmp269*x;
double Lstmp271 = z*L[108];
double Lstmp272 = Lstmp271*x;
double Lstmp273 = z*L[69];
double Lstmp274 = Lstmp273*y;
double Lstmp275 = z*L[97];
double Lstmp276 = Lstmp275*y;
double Lstmp277 = z*L[102];
double Lstmp278 = x*L[59];
double Lstmp279 = x*L[90];
double Lstmp280 = x*L[61];
double Lstmp281 = x*L[93];
double Lstmp282 = Lstmp133*x;
double Lstmp283 = Lstmp134*x;
double Lstmp284 = x*L[64];
double Lstmp285 = x*L[97];
double Lstmp286 = Lstmp179*x;
double Lstmp287 = x*L[63];
double Lstmp288 = x*L[95];
double Lstmp289 = Lstmp215*x;
double Lstmp290 = Lstmp251*x;
double Lstmp291 = Lstmp277*x;
double Lstmp292 = y*L[45];
double Lstmp293 = Lstmp74*y;
double Lstmp294 = y*L[66];
double Lstmp295 = y*L[94];
double Lstmp296 = y*L[79];
double Lstmp297 = y*L[115];
double Lstmp298 = y*L[107];
double Lstmp299 = Lstmp298*x;
double Lstmp300 = Lstmp128*y;
double Lstmp301 = y*L[46];
double Lstmp302 = Lstmp203*y;
double Lstmp303 = y*L[67];
double Lstmp304 = y*L[95];
double Lstmp305 = Lstmp214*y;
double Lstmp306 = Lstmp269*y;
double Lstmp307 = Lstmp277*y;
double Lstmp308 = z*L[49];
double Lstmp309 = z*L[54];
double Lstmp310 = z*L[75];
double Lstmp311 = Lstmp310*y;
double Lstmp312 = z*L[70];
double Lstmp313 = z*L[98];
double Lstmp314 = z*L[81];
double Lstmp315 = z*L[116];
double Lstmp316 = z*L[109];
double Lstmp317 = Lstmp316*x;
double Lstmp318 = z*L[103];
double Lstmp319 = Lstmp318*y;
double Lstmp320 = x*L[87];
double Lstmp321 = x*L[89];
double Lstmp322 = x*L[92];
double Lstmp323 = x*L[91];
double Lstmp324 = y*L[71];
double Lstmp325 = Lstmp76*y;
double Lstmp326 = y*L[99];
double Lstmp327 = y*L[114];
double Lstmp328 = y*L[72];
double Lstmp329 = Lstmp205*y;
double Lstmp330 = y*L[100];
double Lstmp331 = Lstmp271*y;
double Lstmp332 = Lstmp316*y;
double Lstmp333 = z*L[76];
double Lstmp334 = z*L[82];
double Lstmp335 = z*L[110];
double Lstmp336 = Lstmp335*y;
double Lstmp337 = z*L[104];
double Lstmp338 = z*L[117];
double Lstmp339 = y*L[105];
double Lstmp340 = y*L[106];
double Lstmp341 = z*L[111];
double Lstmp342 = z*L[118];
#pragma omp atomic
Ls[0] += Lstmp0*x + Lstmp1*x + Lstmp10*Lstmp46 + Lstmp10*Lstmp47 + Lstmp10*Lstmp83 + Lstmp10*L[20] + Lstmp100*Lstmp18*L[94] + Lstmp100*Lstmp27*L[98] + Lstmp101*Lstmp129 + Lstmp101*L[40] + Lstmp102*Lstmp130 + Lstmp102*L[65] + Lstmp104*Lstmp134 + Lstmp104*L[59] + Lstmp105*Lstmp16*L[90] + Lstmp105*Lstmp25*L[93] + Lstmp106*Lstmp133 + Lstmp106*L[61] + Lstmp107*Lstmp14*L[87] + Lstmp107*Lstmp23*L[89] + (1.0/24.0)*Lstmp108*Lstmp7*L[96] + Lstmp109*Lstmp117 + Lstmp109*L[32] + Lstmp110*Lstmp136*L[102] + Lstmp111*Lstmp118 + Lstmp111*L[53] + Lstmp112*Lstmp119 + Lstmp112*L[81] + Lstmp113*Lstmp136*L[101] + Lstmp114*Lstmp120 + Lstmp114*L[52] + Lstmp115*Lstmp121 + Lstmp115*L[80] + Lstmp116*Lstmp122 + Lstmp116*L[79] + Lstmp12*Lstmp48 + Lstmp12*Lstmp49 + Lstmp12*Lstmp85 + Lstmp12*L[35] + Lstmp123*Lstmp93 + Lstmp124*Lstmp94 + Lstmp125*Lstmp95 + Lstmp126*Lstmp87 + Lstmp127*Lstmp89 + Lstmp128*Lstmp91 + Lstmp13*Lstmp50 + Lstmp13*Lstmp51 + Lstmp13*L[56] + Lstmp131*Lstmp97 + Lstmp132*Lstmp99 + Lstmp135*L[68] + (1.0/240.0)*Lstmp14*Lstmp29*L[117] + Lstmp15*Lstmp32 + Lstmp15*Lstmp57 + Lstmp15*Lstmp71 + Lstmp15*L[7] + (1.0/144.0)*Lstmp16*Lstmp27*L[116] + Lstmp17*Lstmp33 + Lstmp17*Lstmp58 + Lstmp17*Lstmp73 + Lstmp17*L[16] + (1.0/144.0)*Lstmp18*Lstmp25*L[115] + Lstmp19*Lstmp34 + Lstmp19*Lstmp59 + Lstmp19*Lstmp75 + Lstmp19*L[30] + Lstmp2*y + (1.0/240.0)*Lstmp20*Lstmp23*L[114] + Lstmp20*Lstmp92*L[99] + Lstmp21*Lstmp35 + Lstmp21*Lstmp60 + Lstmp21*Lstmp77 + Lstmp21*L[50] + Lstmp22*Lstmp36 + Lstmp22*Lstmp61 + Lstmp22*L[77] + Lstmp24*Lstmp37 + Lstmp24*Lstmp52 + Lstmp24*Lstmp63 + Lstmp24*L[9] + Lstmp26*Lstmp38 + Lstmp26*Lstmp53 + Lstmp26*Lstmp65 + Lstmp26*L[19] + Lstmp28*Lstmp39 + Lstmp28*Lstmp54 + Lstmp28*Lstmp67 + Lstmp28*L[34] + Lstmp29*Lstmp92*L[104] + Lstmp30*Lstmp40 + Lstmp30*Lstmp55 + Lstmp30*Lstmp69 + Lstmp30*L[55] + Lstmp31*Lstmp41 + Lstmp31*Lstmp56 + Lstmp31*L[83] + Lstmp4*x + Lstmp42*Lstmp6 + Lstmp43*Lstmp6 + Lstmp44*Lstmp8 + Lstmp45*Lstmp8 + Lstmp6*Lstmp79 + Lstmp6*L[4] + Lstmp8*Lstmp81 + Lstmp8*L[10] + Lstmp87*L[23] + Lstmp89*L[41] + Lstmp91*L[66] + Lstmp93*L[25] + Lstmp94*L[44] + Lstmp95*L[70] + Lstmp97*L[38] + Lstmp99*L[62] + (1.0/5040.0)*(x*x*x*x*x*x*x)*L[84] + x*L[1] + (1.0/5040.0)*(y*y*y*y*y*y*y)*L[112] + y*L[2] + (1.0/5040.0)*(z*z*z*z*z*z*z)*L[119] + z*L[3] + L[0];
#pragma omp atomic
Ls[1] += Lstmp0 + Lstmp1 + Lstmp10*Lstmp48 + Lstmp10*Lstmp49 + Lstmp10*Lstmp85 + Lstmp10*L[35] + Lstmp101*Lstmp133 + Lstmp101*L[61] + Lstmp102*L[93] + Lstmp104*L[87] + Lstmp106*L[89] + Lstmp109*Lstmp151 + Lstmp109*L[47] + Lstmp111*Lstmp152 + Lstmp111*L[74] + Lstmp112*L[109] + Lstmp114*Lstmp153 + Lstmp114*L[73] + Lstmp115*L[108] + Lstmp116*L[107] + Lstmp12*Lstmp50 + Lstmp12*Lstmp51 + Lstmp12*L[56] + Lstmp129*Lstmp93 + Lstmp13*L[84] + Lstmp130*Lstmp94 + Lstmp131*Lstmp87 + Lstmp132*Lstmp89 + Lstmp134*Lstmp97 + Lstmp135*L[96] + Lstmp137*Lstmp15 + Lstmp138*Lstmp17 + Lstmp139*Lstmp19 + Lstmp140*Lstmp21 + Lstmp141*Lstmp24 + Lstmp142*Lstmp26 + Lstmp143*Lstmp28 + Lstmp144*Lstmp30 + Lstmp145*Lstmp24 + Lstmp146*Lstmp26 + Lstmp147*Lstmp28 + Lstmp148*Lstmp15 + Lstmp149*Lstmp17 + Lstmp15*Lstmp70 + Lstmp15*L[13] + Lstmp150*Lstmp19 + Lstmp17*Lstmp72 + Lstmp17*L[26] + Lstmp19*Lstmp74 + Lstmp19*L[45] + Lstmp21*Lstmp76 + Lstmp21*L[71] + Lstmp22*L[105] + Lstmp24*Lstmp62 + Lstmp24*L[15] + Lstmp26*Lstmp64 + Lstmp26*L[29] + Lstmp28*Lstmp66 + Lstmp28*L[49] + Lstmp30*Lstmp68 + Lstmp30*L[76] + Lstmp31*L[111] + Lstmp4 + Lstmp42*x + Lstmp43*x + Lstmp44*Lstmp6 + Lstmp45*Lstmp6 + Lstmp46*Lstmp8 + Lstmp47*Lstmp8 + Lstmp6*Lstmp81 + Lstmp6*L[10] + Lstmp79*x + Lstmp8*Lstmp83 + Lstmp8*L[20] + Lstmp87*L[38] + Lstmp89*L[62] + Lstmp91*L[94] + Lstmp93*L[40] + Lstmp94*L[65] + Lstmp95*L[98] + Lstmp97*L[59] + Lstmp99*L[90] + x*L[4] + L[1];
#pragma omp atomic
Ls[2] += Lstmp10*Lstmp162 + Lstmp10*Lstmp176 + Lstmp10*Lstmp82 + Lstmp10*L[36] + Lstmp101*Lstmp179 + Lstmp101*L[64] + Lstmp102*L[97] + Lstmp104*L[90] + Lstmp106*L[92] + Lstmp109*Lstmp120 + Lstmp109*L[52] + Lstmp111*Lstmp121 + Lstmp111*L[80] + Lstmp112*L[116] + Lstmp114*Lstmp122 + Lstmp114*L[79] + Lstmp115*L[115] + Lstmp116*L[114] + Lstmp12*Lstmp163 + Lstmp12*Lstmp84 + Lstmp12*L[57] + Lstmp127*Lstmp87 + Lstmp128*Lstmp89 + Lstmp13*L[85] + Lstmp132*Lstmp97 + Lstmp135*L[101] + Lstmp15*Lstmp33 + Lstmp15*Lstmp58 + Lstmp15*Lstmp73 + Lstmp15*L[16] + Lstmp154*x + Lstmp155*x + Lstmp156*Lstmp24 + Lstmp157*Lstmp26 + Lstmp158*Lstmp28 + Lstmp159*Lstmp30 + Lstmp160*Lstmp6 + Lstmp161*Lstmp8 + Lstmp164*Lstmp24 + Lstmp165*Lstmp26 + Lstmp166*Lstmp28 + Lstmp167*Lstmp30 + Lstmp169*Lstmp24 + Lstmp17*Lstmp34 + Lstmp17*Lstmp59 + Lstmp17*Lstmp75 + Lstmp17*L[30] + Lstmp171*Lstmp26 + Lstmp173*Lstmp28 + Lstmp174*Lstmp6 + Lstmp175*Lstmp8 + Lstmp177*Lstmp93 + Lstmp178*Lstmp94 + Lstmp19*Lstmp35 + Lstmp19*Lstmp60 + Lstmp19*Lstmp77 + Lstmp19*L[50] + Lstmp2 + Lstmp21*Lstmp36 + Lstmp21*Lstmp61 + Lstmp21*L[77] + Lstmp22*L[112] + Lstmp24*L[18] + Lstmp26*L[33] + Lstmp28*L[54] + Lstmp3*x + Lstmp30*L[82] + Lstmp31*L[118] + Lstmp57*y + Lstmp6*Lstmp78 + Lstmp6*L[11] + Lstmp8*Lstmp80 + Lstmp8*L[21] + Lstmp87*L[41] + Lstmp89*L[66] + Lstmp91*L[99] + Lstmp93*L[43] + Lstmp94*L[69] + Lstmp95*L[103] + Lstmp97*L[62] + Lstmp99*L[94] + x*L[5] + y*L[7] + L[2];
#pragma omp atomic
Ls[3] += Lstmp10*Lstmp193 + Lstmp10*Lstmp194 + Lstmp10*Lstmp212 + Lstmp10*L[37] + Lstmp101*Lstmp130 + Lstmp101*L[65] + Lstmp102*L[98] + Lstmp104*L[91] + Lstmp106*L[93] + Lstmp109*Lstmp118 + Lstmp109*L[53] + Lstmp111*Lstmp119 + Lstmp111*L[81] + Lstmp112*L[117] + Lstmp114*Lstmp121 + Lstmp114*L[80] + Lstmp115*L[116] + Lstmp116*L[115] + Lstmp12*Lstmp195 + Lstmp12*Lstmp196 + Lstmp12*L[58] + Lstmp124*Lstmp93 + Lstmp125*Lstmp94 + Lstmp13*L[86] + Lstmp135*L[102] + Lstmp15*Lstmp185 + Lstmp15*Lstmp197 + Lstmp15*Lstmp202 + Lstmp15*L[17] + Lstmp17*Lstmp186 + Lstmp17*Lstmp198 + Lstmp17*Lstmp204 + Lstmp17*L[31] + Lstmp180*x + Lstmp181*x + Lstmp182*y + Lstmp184*x + Lstmp187*Lstmp19 + Lstmp188*Lstmp21 + Lstmp189*Lstmp6 + Lstmp19*Lstmp199 + Lstmp19*Lstmp206 + Lstmp19*L[51] + Lstmp190*Lstmp6 + Lstmp191*Lstmp8 + Lstmp192*Lstmp8 + Lstmp200*Lstmp21 + Lstmp208*Lstmp6 + Lstmp21*L[78] + Lstmp210*Lstmp8 + Lstmp213*Lstmp87 + Lstmp214*Lstmp89 + Lstmp215*Lstmp97 + Lstmp22*L[113] + Lstmp24*Lstmp38 + Lstmp24*Lstmp53 + Lstmp24*Lstmp65 + Lstmp24*L[19] + Lstmp26*Lstmp39 + Lstmp26*Lstmp54 + Lstmp26*Lstmp67 + Lstmp26*L[34] + Lstmp28*Lstmp40 + Lstmp28*Lstmp55 + Lstmp28*Lstmp69 + Lstmp28*L[55] + Lstmp30*Lstmp41 + Lstmp30*Lstmp56 + Lstmp30*L[83] + Lstmp31*L[119] + Lstmp6*L[12] + Lstmp8*L[22] + Lstmp87*L[42] + Lstmp89*L[67] + Lstmp91*L[100] + Lstmp93*L[44] + Lstmp94*L[70] + Lstmp95*L[104] + Lstmp97*L[63] + Lstmp99*L[95] + x*L[6] + y*L[8] + z*L[9] + L[3];
#pragma omp atomic
Ls[4] += Lstmp10*Lstmp50 + Lstmp10*Lstmp51 + Lstmp10*L[56] + Lstmp101*L[89] + Lstmp109*Lstmp226 + Lstmp109*L[68] + Lstmp111*L[102] + Lstmp114*L[101] + Lstmp12*L[84] + Lstmp123*Lstmp24 + Lstmp124*Lstmp26 + Lstmp125*Lstmp28 + Lstmp126*Lstmp15 + Lstmp127*Lstmp17 + Lstmp128*Lstmp19 + Lstmp133*Lstmp93 + Lstmp134*Lstmp87 + Lstmp15*Lstmp216 + Lstmp15*Lstmp224 + Lstmp15*L[23] + Lstmp17*Lstmp217 + Lstmp17*Lstmp225 + Lstmp17*L[41] + Lstmp19*Lstmp218 + Lstmp19*L[66] + Lstmp21*L[99] + Lstmp219*Lstmp24 + Lstmp220*Lstmp26 + Lstmp221*Lstmp28 + Lstmp222*Lstmp24 + Lstmp223*Lstmp26 + Lstmp24*L[25] + Lstmp26*L[44] + Lstmp28*L[70] + Lstmp30*L[104] + Lstmp42 + Lstmp43 + Lstmp44*x + Lstmp45*x + Lstmp46*Lstmp6 + Lstmp47*Lstmp6 + Lstmp48*Lstmp8 + Lstmp49*Lstmp8 + Lstmp6*Lstmp83 + Lstmp6*L[20] + Lstmp79 + Lstmp8*Lstmp85 + Lstmp8*L[35] + Lstmp81*x + Lstmp87*L[59] + Lstmp89*L[90] + Lstmp93*L[61] + Lstmp94*L[93] + Lstmp97*L[87] + x*L[10] + L[4];
#pragma omp atomic
Ls[5] += Lstmp10*Lstmp163 + Lstmp10*Lstmp84 + Lstmp10*L[57] + Lstmp101*L[92] + Lstmp109*Lstmp153 + Lstmp109*L[73] + Lstmp111*L[108] + Lstmp114*L[107] + Lstmp12*L[85] + Lstmp132*Lstmp87 + Lstmp138*Lstmp15 + Lstmp139*Lstmp17 + Lstmp140*Lstmp19 + Lstmp149*Lstmp15 + Lstmp15*Lstmp72 + Lstmp15*L[26] + Lstmp150*Lstmp17 + Lstmp154 + Lstmp155 + Lstmp160*x + Lstmp161*Lstmp6 + Lstmp162*Lstmp8 + Lstmp168*Lstmp24 + Lstmp17*Lstmp74 + Lstmp17*L[45] + Lstmp170*Lstmp26 + Lstmp172*Lstmp28 + Lstmp174*x + Lstmp175*Lstmp6 + Lstmp176*Lstmp8 + Lstmp179*Lstmp93 + Lstmp19*Lstmp76 + Lstmp19*L[71] + Lstmp21*L[105] + Lstmp227*Lstmp24 + Lstmp228*Lstmp26 + Lstmp229*Lstmp28 + Lstmp230*Lstmp24 + Lstmp231*Lstmp26 + Lstmp24*L[28] + Lstmp26*L[48] + Lstmp28*L[75] + Lstmp3 + Lstmp30*L[110] + Lstmp6*Lstmp80 + Lstmp6*L[21] + Lstmp78*x + Lstmp8*Lstmp82 + Lstmp8*L[36] + Lstmp87*L[62] + Lstmp89*L[94] + Lstmp93*L[64] + Lstmp94*L[97] + Lstmp97*L[90] + x*L[11] + L[5];
#pragma omp atomic
Ls[6] += Lstmp10*Lstmp195 + Lstmp10*Lstmp196 + Lstmp10*L[58] + Lstmp101*L[93] + Lstmp109*Lstmp152 + Lstmp109*L[74] + Lstmp111*L[109] + Lstmp114*L[108] + Lstmp12*L[86] + Lstmp130*Lstmp93 + Lstmp142*Lstmp24 + Lstmp143*Lstmp26 + Lstmp144*Lstmp28 + Lstmp146*Lstmp24 + Lstmp147*Lstmp26 + Lstmp15*Lstmp201 + Lstmp15*Lstmp232 + Lstmp15*Lstmp235 + Lstmp15*L[27] + Lstmp17*Lstmp203 + Lstmp17*Lstmp233 + Lstmp17*Lstmp236 + Lstmp17*L[46] + Lstmp180 + Lstmp181 + Lstmp184 + Lstmp189*x + Lstmp19*Lstmp205 + Lstmp19*Lstmp234 + Lstmp19*L[72] + Lstmp190*x + Lstmp191*Lstmp6 + Lstmp192*Lstmp6 + Lstmp193*Lstmp8 + Lstmp194*Lstmp8 + Lstmp208*x + Lstmp21*L[106] + Lstmp210*Lstmp6 + Lstmp212*Lstmp8 + Lstmp215*Lstmp87 + Lstmp24*Lstmp64 + Lstmp24*L[29] + Lstmp26*Lstmp66 + Lstmp26*L[49] + Lstmp28*Lstmp68 + Lstmp28*L[76] + Lstmp30*L[111] + Lstmp6*L[22] + Lstmp8*L[37] + Lstmp87*L[63] + Lstmp89*L[95] + Lstmp93*L[65] + Lstmp94*L[98] + Lstmp97*L[91] + x*L[12] + L[6];
#pragma omp atomic
Ls[7] += Lstmp10*Lstmp134 + Lstmp10*Lstmp241 + Lstmp10*L[59] + Lstmp101*L[96] + Lstmp109*Lstmp122 + Lstmp109*L[79] + Lstmp111*L[115] + Lstmp114*L[114] + Lstmp117*Lstmp24 + Lstmp118*Lstmp26 + Lstmp119*Lstmp28 + Lstmp12*L[87] + Lstmp126*Lstmp6 + Lstmp128*Lstmp87 + Lstmp131*Lstmp8 + Lstmp15*Lstmp34 + Lstmp15*Lstmp59 + Lstmp15*Lstmp75 + Lstmp15*L[30] + Lstmp17*Lstmp35 + Lstmp17*Lstmp60 + Lstmp17*Lstmp77 + Lstmp17*L[50] + Lstmp19*Lstmp36 + Lstmp19*Lstmp61 + Lstmp19*L[77] + Lstmp21*L[112] + Lstmp237*x + Lstmp238*x + Lstmp239*Lstmp6 + Lstmp24*Lstmp242 + Lstmp24*Lstmp246 + Lstmp24*L[32] + Lstmp240*Lstmp8 + Lstmp243*Lstmp26 + Lstmp244*Lstmp28 + Lstmp248*Lstmp26 + Lstmp249*Lstmp6 + Lstmp250*Lstmp8 + Lstmp251*Lstmp93 + Lstmp26*L[53] + Lstmp28*L[81] + Lstmp30*L[117] + Lstmp32 + Lstmp57 + Lstmp58*y + Lstmp6*L[23] + Lstmp71 + Lstmp8*L[38] + Lstmp87*L[66] + Lstmp89*L[99] + Lstmp93*L[68] + Lstmp94*L[102] + Lstmp97*L[94] + y*L[16] + L[7];
#pragma omp atomic
Ls[8] += Lstmp10*Lstmp211 + Lstmp10*Lstmp256 + Lstmp10*L[60] + Lstmp101*L[97] + Lstmp109*Lstmp121 + Lstmp109*L[80] + Lstmp111*L[116] + Lstmp114*L[115] + Lstmp12*L[88] + Lstmp15*Lstmp186 + Lstmp15*Lstmp198 + Lstmp15*Lstmp204 + Lstmp15*L[31] + Lstmp157*Lstmp24 + Lstmp158*Lstmp26 + Lstmp159*Lstmp28 + Lstmp165*Lstmp24 + Lstmp166*Lstmp26 + Lstmp167*Lstmp28 + Lstmp17*Lstmp187 + Lstmp17*Lstmp199 + Lstmp17*Lstmp206 + Lstmp17*L[51] + Lstmp171*Lstmp24 + Lstmp173*Lstmp26 + Lstmp178*Lstmp93 + Lstmp182 + Lstmp183*x + Lstmp188*Lstmp19 + Lstmp19*Lstmp200 + Lstmp19*L[78] + Lstmp197*y + Lstmp207*Lstmp6 + Lstmp209*Lstmp8 + Lstmp21*L[113] + Lstmp214*Lstmp87 + Lstmp24*L[33] + Lstmp252*x + Lstmp253*x + Lstmp254*Lstmp6 + Lstmp255*Lstmp8 + Lstmp257*Lstmp6 + Lstmp258*Lstmp8 + Lstmp26*L[54] + Lstmp28*L[82] + Lstmp30*L[118] + Lstmp6*L[24] + Lstmp8*L[39] + Lstmp87*L[67] + Lstmp89*L[100] + Lstmp93*L[69] + Lstmp94*L[103] + Lstmp97*L[95] + x*L[14] + y*L[17] + L[8];
#pragma omp atomic
Ls[9] += Lstmp10*Lstmp133 + Lstmp10*Lstmp265 + Lstmp10*L[61] + Lstmp101*L[98] + Lstmp109*Lstmp119 + Lstmp109*L[81] + Lstmp111*L[117] + Lstmp114*L[116] + Lstmp117*Lstmp15 + Lstmp12*L[89] + Lstmp120*Lstmp17 + Lstmp122*Lstmp19 + Lstmp123*Lstmp6 + Lstmp125*Lstmp93 + Lstmp129*Lstmp8 + Lstmp15*Lstmp266 + Lstmp15*Lstmp270 + Lstmp15*L[32] + Lstmp17*Lstmp267 + Lstmp17*Lstmp272 + Lstmp17*L[52] + Lstmp19*Lstmp268 + Lstmp19*L[79] + Lstmp21*L[114] + Lstmp24*Lstmp39 + Lstmp24*Lstmp54 + Lstmp24*Lstmp67 + Lstmp24*L[34] + Lstmp259*x + Lstmp26*Lstmp40 + Lstmp26*Lstmp55 + Lstmp26*Lstmp69 + Lstmp26*L[55] + Lstmp260*y + Lstmp262*x + Lstmp263*Lstmp6 + Lstmp264*Lstmp8 + Lstmp274*Lstmp6 + Lstmp276*Lstmp8 + Lstmp277*Lstmp87 + Lstmp28*Lstmp41 + Lstmp28*Lstmp56 + Lstmp28*L[83] + Lstmp30*L[119] + Lstmp37 + Lstmp52 + Lstmp6*L[25] + Lstmp63 + Lstmp8*L[40] + Lstmp87*L[68] + Lstmp89*L[101] + Lstmp93*L[70] + Lstmp94*L[104] + Lstmp97*L[96] + z*L[19] + L[9];
#pragma omp atomic
Ls[10] += Lstmp10*L[84] + Lstmp109*L[96] + Lstmp129*Lstmp24 + Lstmp130*Lstmp26 + Lstmp131*Lstmp15 + Lstmp132*Lstmp17 + Lstmp15*Lstmp278 + Lstmp15*Lstmp283 + Lstmp15*L[38] + Lstmp17*Lstmp279 + Lstmp17*L[62] + Lstmp19*L[94] + Lstmp24*Lstmp280 + Lstmp24*Lstmp282 + Lstmp24*L[40] + Lstmp26*Lstmp281 + Lstmp26*L[65] + Lstmp28*L[98] + Lstmp44 + Lstmp45 + Lstmp46*x + Lstmp47*x + Lstmp48*Lstmp6 + Lstmp49*Lstmp6 + Lstmp50*Lstmp8 + Lstmp51*Lstmp8 + Lstmp6*Lstmp85 + Lstmp6*L[35] + Lstmp8*L[56] + Lstmp81 + Lstmp83*x + Lstmp87*L[87] + Lstmp93*L[89] + x*L[20] + L[10];
#pragma omp atomic
Ls[11] += Lstmp10*L[85] + Lstmp109*L[101] + Lstmp127*Lstmp15 + Lstmp128*Lstmp17 + Lstmp15*Lstmp217 + Lstmp15*Lstmp225 + Lstmp15*L[41] + Lstmp160 + Lstmp161*x + Lstmp162*Lstmp6 + Lstmp163*Lstmp8 + Lstmp17*Lstmp218 + Lstmp17*L[66] + Lstmp174 + Lstmp175*x + Lstmp176*Lstmp6 + Lstmp177*Lstmp24 + Lstmp178*Lstmp26 + Lstmp19*L[99] + Lstmp24*Lstmp284 + Lstmp24*Lstmp286 + Lstmp24*L[43] + Lstmp26*Lstmp285 + Lstmp26*L[69] + Lstmp28*L[103] + Lstmp6*Lstmp82 + Lstmp6*L[36] + Lstmp78 + Lstmp8*Lstmp84 + Lstmp8*L[57] + Lstmp80*x + Lstmp87*L[90] + Lstmp93*L[92] + x*L[21] + L[11];
#pragma omp atomic
Ls[12] += Lstmp10*L[86] + Lstmp109*L[102] + Lstmp124*Lstmp24 + Lstmp125*Lstmp26 + Lstmp15*Lstmp213 + Lstmp15*Lstmp287 + Lstmp15*Lstmp289 + Lstmp15*L[42] + Lstmp17*Lstmp214 + Lstmp17*Lstmp288 + Lstmp17*L[67] + Lstmp189 + Lstmp19*L[100] + Lstmp190 + Lstmp191*x + Lstmp192*x + Lstmp193*Lstmp6 + Lstmp194*Lstmp6 + Lstmp195*Lstmp8 + Lstmp196*Lstmp8 + Lstmp208 + Lstmp210*x + Lstmp212*Lstmp6 + Lstmp220*Lstmp24 + Lstmp221*Lstmp26 + Lstmp223*Lstmp24 + Lstmp24*L[44] + Lstmp26*L[70] + Lstmp28*L[104] + Lstmp6*L[37] + Lstmp8*L[58] + Lstmp87*L[91] + Lstmp93*L[93] + x*L[22] + L[12];
#pragma omp atomic
Ls[13] += Lstmp10*L[87] + Lstmp109*L[107] + Lstmp131*Lstmp6 + Lstmp134*Lstmp8 + Lstmp137 + Lstmp139*Lstmp15 + Lstmp140*Lstmp17 + Lstmp148 + Lstmp15*Lstmp150 + Lstmp15*Lstmp74 + Lstmp15*L[45] + Lstmp151*Lstmp24 + Lstmp152*Lstmp26 + Lstmp17*Lstmp76 + Lstmp17*L[71] + Lstmp19*L[105] + Lstmp237 + Lstmp238 + Lstmp239*x + Lstmp24*Lstmp245 + Lstmp24*Lstmp290 + Lstmp24*L[47] + Lstmp240*Lstmp6 + Lstmp241*Lstmp8 + Lstmp247*Lstmp26 + Lstmp249*x + Lstmp250*Lstmp6 + Lstmp26*L[74] + Lstmp28*L[109] + Lstmp6*L[38] + Lstmp70 + Lstmp8*L[59] + Lstmp87*L[94] + Lstmp93*L[96] + L[13];
#pragma omp atomic
Ls[14] += Lstmp10*L[88] + Lstmp109*L[108] + Lstmp15*Lstmp203 + Lstmp15*Lstmp233 + Lstmp15*Lstmp236 + Lstmp15*L[46] + Lstmp17*Lstmp205 + Lstmp17*Lstmp234 + Lstmp17*L[72] + Lstmp170*Lstmp24 + Lstmp172*Lstmp26 + Lstmp183 + Lstmp19*L[106] + Lstmp207*x + Lstmp209*Lstmp6 + Lstmp211*Lstmp8 + Lstmp228*Lstmp24 + Lstmp229*Lstmp26 + Lstmp231*Lstmp24 + Lstmp24*L[48] + Lstmp252 + Lstmp253 + Lstmp254*x + Lstmp255*Lstmp6 + Lstmp256*Lstmp8 + Lstmp257*x + Lstmp258*Lstmp6 + Lstmp26*L[75] + Lstmp28*L[110] + Lstmp6*L[39] + Lstmp8*L[60] + Lstmp87*L[95] + Lstmp93*L[97] + x*L[24] + L[14];
#pragma omp atomic
Ls[15] += Lstmp10*L[89] + Lstmp109*L[109] + Lstmp129*Lstmp6 + Lstmp133*Lstmp8 + Lstmp141 + Lstmp143*Lstmp24 + Lstmp144*Lstmp26 + Lstmp145 + Lstmp147*Lstmp24 + Lstmp15*Lstmp151 + Lstmp15*Lstmp269 + Lstmp15*Lstmp291 + Lstmp15*L[47] + Lstmp153*Lstmp17 + Lstmp17*Lstmp271 + Lstmp17*L[73] + Lstmp19*L[107] + Lstmp24*Lstmp66 + Lstmp24*L[49] + Lstmp259 + Lstmp26*Lstmp68 + Lstmp26*L[76] + Lstmp262 + Lstmp263*x + Lstmp264*Lstmp6 + Lstmp265*Lstmp8 + Lstmp274*x + Lstmp276*Lstmp6 + Lstmp28*L[111] + Lstmp6*L[40] + Lstmp62 + Lstmp8*L[61] + Lstmp87*L[96] + Lstmp93*L[98] + L[15];
#pragma omp atomic
Ls[16] += Lstmp10*L[90] + Lstmp109*L[114] + Lstmp120*Lstmp24 + Lstmp121*Lstmp26 + Lstmp127*Lstmp6 + Lstmp132*Lstmp8 + Lstmp15*Lstmp35 + Lstmp15*Lstmp60 + Lstmp15*Lstmp77 + Lstmp15*L[50] + Lstmp17*Lstmp36 + Lstmp17*Lstmp61 + Lstmp17*L[77] + Lstmp19*L[112] + Lstmp24*Lstmp296 + Lstmp24*Lstmp299 + Lstmp24*L[52] + Lstmp26*Lstmp297 + Lstmp26*L[80] + Lstmp28*L[116] + Lstmp292*x + Lstmp293*x + Lstmp294*Lstmp6 + Lstmp295*Lstmp8 + Lstmp300*Lstmp6 + Lstmp33 + Lstmp58 + Lstmp59*y + Lstmp6*L[41] + Lstmp73 + Lstmp8*L[62] + Lstmp87*L[99] + Lstmp93*L[101] + y*L[30] + L[16];
#pragma omp atomic
Ls[17] += Lstmp10*L[91] + Lstmp109*L[115] + Lstmp118*Lstmp24 + Lstmp119*Lstmp26 + Lstmp15*Lstmp187 + Lstmp15*Lstmp199 + Lstmp15*Lstmp206 + Lstmp15*L[51] + Lstmp17*Lstmp188 + Lstmp17*Lstmp200 + Lstmp17*L[78] + Lstmp185 + Lstmp19*L[113] + Lstmp197 + Lstmp198*y + Lstmp202 + Lstmp213*Lstmp6 + Lstmp215*Lstmp8 + Lstmp24*Lstmp243 + Lstmp24*Lstmp248 + Lstmp24*L[53] + Lstmp244*Lstmp26 + Lstmp26*L[81] + Lstmp28*L[117] + Lstmp301*x + Lstmp302*x + Lstmp303*Lstmp6 + Lstmp304*Lstmp8 + Lstmp305*Lstmp6 + Lstmp6*L[42] + Lstmp8*L[63] + Lstmp87*L[100] + Lstmp93*L[102] + y*L[31] + L[17];
#pragma omp atomic
Ls[18] += Lstmp10*L[92] + Lstmp109*L[116] + Lstmp120*Lstmp15 + Lstmp122*Lstmp17 + Lstmp15*Lstmp267 + Lstmp15*Lstmp272 + Lstmp15*L[52] + Lstmp156 + Lstmp158*Lstmp24 + Lstmp159*Lstmp26 + Lstmp164 + Lstmp166*Lstmp24 + Lstmp167*Lstmp26 + Lstmp169 + Lstmp17*Lstmp268 + Lstmp17*L[79] + Lstmp173*Lstmp24 + Lstmp177*Lstmp6 + Lstmp179*Lstmp8 + Lstmp19*L[114] + Lstmp24*L[54] + Lstmp26*L[82] + Lstmp260 + Lstmp261*x + Lstmp266*y + Lstmp273*Lstmp6 + Lstmp275*Lstmp8 + Lstmp28*L[118] + Lstmp306*x + Lstmp307*Lstmp6 + Lstmp6*L[43] + Lstmp8*L[64] + Lstmp87*L[101] + Lstmp93*L[103] + L[18];
#pragma omp atomic
Ls[19] += Lstmp10*L[93] + Lstmp109*L[117] + Lstmp118*Lstmp15 + Lstmp121*Lstmp17 + Lstmp124*Lstmp6 + Lstmp130*Lstmp8 + Lstmp15*Lstmp314 + Lstmp15*Lstmp317 + Lstmp15*L[53] + Lstmp17*Lstmp315 + Lstmp17*L[80] + Lstmp19*L[115] + Lstmp24*Lstmp40 + Lstmp24*Lstmp55 + Lstmp24*Lstmp69 + Lstmp24*L[55] + Lstmp26*Lstmp41 + Lstmp26*Lstmp56 + Lstmp26*L[83] + Lstmp28*L[119] + Lstmp308*x + Lstmp309*y + Lstmp311*x + Lstmp312*Lstmp6 + Lstmp313*Lstmp8 + Lstmp319*Lstmp6 + Lstmp38 + Lstmp53 + Lstmp6*L[44] + Lstmp65 + Lstmp8*L[65] + Lstmp87*L[102] + Lstmp93*L[104] + z*L[34] + L[19];
#pragma omp atomic
Ls[20] += Lstmp133*Lstmp24 + Lstmp134*Lstmp15 + Lstmp15*Lstmp320 + Lstmp15*L[59] + Lstmp17*L[90] + Lstmp24*Lstmp321 + Lstmp24*L[61] + Lstmp26*L[93] + Lstmp46 + Lstmp47 + Lstmp48*x + Lstmp49*x + Lstmp50*Lstmp6 + Lstmp51*Lstmp6 + Lstmp6*L[56] + Lstmp8*L[84] + Lstmp83 + Lstmp85*x + x*L[35] + L[20];
#pragma omp atomic
Ls[21] += Lstmp132*Lstmp15 + Lstmp15*Lstmp279 + Lstmp15*L[62] + Lstmp161 + Lstmp162*x + Lstmp163*Lstmp6 + Lstmp17*L[94] + Lstmp175 + Lstmp176*x + Lstmp179*Lstmp24 + Lstmp24*Lstmp322 + Lstmp24*L[64] + Lstmp26*L[97] + Lstmp6*Lstmp84 + Lstmp6*L[57] + Lstmp8*L[85] + Lstmp80 + Lstmp82*x + x*L[36] + L[21];
#pragma omp atomic
Ls[22] += Lstmp130*Lstmp24 + Lstmp15*Lstmp215 + Lstmp15*Lstmp323 + Lstmp15*L[63] + Lstmp17*L[95] + Lstmp191 + Lstmp192 + Lstmp193*x + Lstmp194*x + Lstmp195*Lstmp6 + Lstmp196*Lstmp6 + Lstmp210 + Lstmp212*x + Lstmp24*Lstmp281 + Lstmp24*L[65] + Lstmp26*L[98] + Lstmp6*L[58] + Lstmp8*L[86] + x*L[37] + L[22];
#pragma omp atomic
Ls[23] += Lstmp126 + Lstmp128*Lstmp15 + Lstmp134*Lstmp6 + Lstmp15*Lstmp218 + Lstmp15*L[66] + Lstmp17*L[99] + Lstmp216 + Lstmp224 + Lstmp226*Lstmp24 + Lstmp239 + Lstmp24*Lstmp251 + Lstmp24*L[68] + Lstmp240*x + Lstmp241*Lstmp6 + Lstmp249 + Lstmp250*x + Lstmp26*L[102] + Lstmp6*L[59] + Lstmp8*L[87] + L[23];
#pragma omp atomic
Ls[24] += Lstmp15*Lstmp214 + Lstmp15*Lstmp288 + Lstmp15*L[67] + Lstmp17*L[100] + Lstmp178*Lstmp24 + Lstmp207 + Lstmp209*x + Lstmp211*Lstmp6 + Lstmp24*Lstmp285 + Lstmp24*L[69] + Lstmp254 + Lstmp255*x + Lstmp256*Lstmp6 + Lstmp257 + Lstmp258*x + Lstmp26*L[103] + Lstmp6*L[60] + Lstmp8*L[88] + x*L[39] + L[24];
#pragma omp atomic
Ls[25] += Lstmp123 + Lstmp125*Lstmp24 + Lstmp133*Lstmp6 + Lstmp15*Lstmp226 + Lstmp15*Lstmp277 + Lstmp15*L[68] + Lstmp17*L[101] + Lstmp219 + Lstmp221*Lstmp24 + Lstmp222 + Lstmp24*L[70] + Lstmp26*L[104] + Lstmp263 + Lstmp264*x + Lstmp265*Lstmp6 + Lstmp274 + Lstmp276*x + Lstmp6*L[61] + Lstmp8*L[89] + L[25];
#pragma omp atomic
Ls[26] += Lstmp132*Lstmp6 + Lstmp138 + Lstmp140*Lstmp15 + Lstmp149 + Lstmp15*Lstmp76 + Lstmp15*L[71] + Lstmp153*Lstmp24 + Lstmp17*L[105] + Lstmp24*Lstmp298 + Lstmp24*L[73] + Lstmp26*L[108] + Lstmp292 + Lstmp293 + Lstmp294*x + Lstmp295*Lstmp6 + Lstmp300*x + Lstmp6*L[62] + Lstmp72 + Lstmp8*L[90] + L[26];
#pragma omp atomic
Ls[27] += Lstmp15*Lstmp205 + Lstmp15*Lstmp234 + Lstmp15*L[72] + Lstmp152*Lstmp24 + Lstmp17*L[106] + Lstmp201 + Lstmp215*Lstmp6 + Lstmp232 + Lstmp235 + Lstmp24*Lstmp247 + Lstmp24*L[74] + Lstmp26*L[109] + Lstmp301 + Lstmp302 + Lstmp303*x + Lstmp304*Lstmp6 + Lstmp305*x + Lstmp6*L[63] + Lstmp8*L[91] + L[27];
#pragma omp atomic
Ls[28] += Lstmp15*Lstmp153 + Lstmp15*Lstmp271 + Lstmp15*L[73] + Lstmp168 + Lstmp17*L[107] + Lstmp172*Lstmp24 + Lstmp179*Lstmp6 + Lstmp227 + Lstmp229*Lstmp24 + Lstmp230 + Lstmp24*L[75] + Lstmp26*L[110] + Lstmp261 + Lstmp273*x + Lstmp275*Lstmp6 + Lstmp306 + Lstmp307*x + Lstmp6*L[64] + Lstmp8*L[92] + L[28];
#pragma omp atomic
Ls[29] += Lstmp130*Lstmp6 + Lstmp142 + Lstmp144*Lstmp24 + Lstmp146 + Lstmp15*Lstmp152 + Lstmp15*Lstmp316 + Lstmp15*L[74] + Lstmp17*L[108] + Lstmp24*Lstmp68 + Lstmp24*L[76] + Lstmp26*L[111] + Lstmp308 + Lstmp311 + Lstmp312*x + Lstmp313*Lstmp6 + Lstmp319*x + Lstmp6*L[65] + Lstmp64 + Lstmp8*L[93] + L[29];
#pragma omp atomic
Ls[30] += Lstmp122*Lstmp24 + Lstmp128*Lstmp6 + Lstmp15*Lstmp36 + Lstmp15*Lstmp61 + Lstmp15*L[77] + Lstmp17*L[112] + Lstmp24*Lstmp327 + Lstmp24*L[79] + Lstmp26*L[115] + Lstmp324*x + Lstmp325*x + Lstmp326*Lstmp6 + Lstmp34 + Lstmp59 + Lstmp6*L[66] + Lstmp60*y + Lstmp75 + Lstmp8*L[94] + y*L[50] + L[30];
#pragma omp atomic
Ls[31] += Lstmp121*Lstmp24 + Lstmp15*Lstmp188 + Lstmp15*Lstmp200 + Lstmp15*L[78] + Lstmp17*L[113] + Lstmp186 + Lstmp198 + Lstmp199*y + Lstmp204 + Lstmp214*Lstmp6 + Lstmp24*Lstmp297 + Lstmp24*L[80] + Lstmp26*L[116] + Lstmp328*x + Lstmp329*x + Lstmp330*Lstmp6 + Lstmp6*L[67] + Lstmp8*L[95] + y*L[51] + L[31];
#pragma omp atomic
Ls[32] += Lstmp117 + Lstmp119*Lstmp24 + Lstmp122*Lstmp15 + Lstmp15*Lstmp268 + Lstmp15*L[79] + Lstmp17*L[114] + Lstmp24*Lstmp244 + Lstmp24*L[81] + Lstmp242 + Lstmp246 + Lstmp251*Lstmp6 + Lstmp26*L[117] + Lstmp266 + Lstmp267*y + Lstmp270 + Lstmp277*Lstmp6 + Lstmp331*x + Lstmp6*L[68] + Lstmp8*L[96] + L[32];
#pragma omp atomic
Ls[33] += Lstmp121*Lstmp15 + Lstmp15*Lstmp315 + Lstmp15*L[80] + Lstmp157 + Lstmp159*Lstmp24 + Lstmp165 + Lstmp167*Lstmp24 + Lstmp17*L[115] + Lstmp171 + Lstmp178*Lstmp6 + Lstmp24*L[82] + Lstmp26*L[118] + Lstmp309 + Lstmp310*x + Lstmp314*y + Lstmp318*Lstmp6 + Lstmp332*x + Lstmp6*L[69] + Lstmp8*L[97] + L[33];
#pragma omp atomic
Ls[34] += Lstmp119*Lstmp15 + Lstmp125*Lstmp6 + Lstmp15*Lstmp338 + Lstmp15*L[81] + Lstmp17*L[116] + Lstmp24*Lstmp41 + Lstmp24*Lstmp56 + Lstmp24*L[83] + Lstmp26*L[119] + Lstmp333*x + Lstmp334*y + Lstmp336*x + Lstmp337*Lstmp6 + Lstmp39 + Lstmp54 + Lstmp6*L[70] + Lstmp67 + Lstmp8*L[98] + z*L[55] + L[34];
#pragma omp atomic
Ls[35] += Lstmp15*L[87] + Lstmp24*L[89] + Lstmp48 + Lstmp49 + Lstmp50*x + Lstmp51*x + Lstmp6*L[84] + Lstmp85 + x*L[56] + L[35];
#pragma omp atomic
Ls[36] += Lstmp15*L[90] + Lstmp162 + Lstmp163*x + Lstmp176 + Lstmp24*L[92] + Lstmp6*L[85] + Lstmp82 + Lstmp84*x + x*L[57] + L[36];
#pragma omp atomic
Ls[37] += Lstmp15*L[91] + Lstmp193 + Lstmp194 + Lstmp195*x + Lstmp196*x + Lstmp212 + Lstmp24*L[93] + Lstmp6*L[86] + x*L[58] + L[37];
#pragma omp atomic
Ls[38] += Lstmp131 + Lstmp15*L[94] + Lstmp24*L[96] + Lstmp240 + Lstmp241*x + Lstmp250 + Lstmp278 + Lstmp283 + Lstmp6*L[87] + L[38];
#pragma omp atomic
Ls[39] += Lstmp15*L[95] + Lstmp209 + Lstmp211*x + Lstmp24*L[97] + Lstmp255 + Lstmp256*x + Lstmp258 + Lstmp6*L[88] + x*L[60] + L[39];
#pragma omp atomic
Ls[40] += Lstmp129 + Lstmp15*L[96] + Lstmp24*L[98] + Lstmp264 + Lstmp265*x + Lstmp276 + Lstmp280 + Lstmp282 + Lstmp6*L[89] + L[40];
#pragma omp atomic
Ls[41] += Lstmp127 + Lstmp15*L[99] + Lstmp217 + Lstmp225 + Lstmp24*L[101] + Lstmp294 + Lstmp295*x + Lstmp300 + Lstmp6*L[90] + L[41];
#pragma omp atomic
Ls[42] += Lstmp15*L[100] + Lstmp213 + Lstmp24*L[102] + Lstmp287 + Lstmp289 + Lstmp303 + Lstmp304*x + Lstmp305 + Lstmp6*L[91] + L[42];
#pragma omp atomic
Ls[43] += Lstmp15*L[101] + Lstmp177 + Lstmp24*L[103] + Lstmp273 + Lstmp275*x + Lstmp284 + Lstmp286 + Lstmp307 + Lstmp6*L[92] + L[43];
#pragma omp atomic
Ls[44] += Lstmp124 + Lstmp15*L[102] + Lstmp220 + Lstmp223 + Lstmp24*L[104] + Lstmp312 + Lstmp313*x + Lstmp319 + Lstmp6*L[93] + L[44];
#pragma omp atomic
Ls[45] += Lstmp139 + Lstmp15*L[105] + Lstmp150 + Lstmp24*L[107] + Lstmp324 + Lstmp325 + Lstmp326*x + Lstmp6*L[94] + Lstmp74 + L[45];
#pragma omp atomic
Ls[46] += Lstmp15*L[106] + Lstmp203 + Lstmp233 + Lstmp236 + Lstmp24*L[108] + Lstmp328 + Lstmp329 + Lstmp330*x + Lstmp6*L[95] + L[46];
#pragma omp atomic
Ls[47] += Lstmp15*L[107] + Lstmp151 + Lstmp24*L[109] + Lstmp245 + Lstmp269 + Lstmp290 + Lstmp291 + Lstmp331 + Lstmp6*L[96] + L[47];
#pragma omp atomic
Ls[48] += Lstmp15*L[108] + Lstmp170 + Lstmp228 + Lstmp231 + Lstmp24*L[110] + Lstmp310 + Lstmp318*x + Lstmp332 + Lstmp6*L[97] + L[48];
#pragma omp atomic
Ls[49] += Lstmp143 + Lstmp147 + Lstmp15*L[109] + Lstmp24*L[111] + Lstmp333 + Lstmp336 + Lstmp337*x + Lstmp6*L[98] + Lstmp66 + L[49];
#pragma omp atomic
Ls[50] += Lstmp15*L[112] + Lstmp24*L[114] + Lstmp339*x + Lstmp35 + Lstmp6*L[99] + Lstmp60 + Lstmp61*y + Lstmp77 + y*L[77] + L[50];
#pragma omp atomic
Ls[51] += Lstmp15*L[113] + Lstmp187 + Lstmp199 + Lstmp200*y + Lstmp206 + Lstmp24*L[115] + Lstmp340*x + Lstmp6*L[100] + y*L[78] + L[51];
#pragma omp atomic
Ls[52] += Lstmp120 + Lstmp15*L[114] + Lstmp24*L[116] + Lstmp267 + Lstmp268*y + Lstmp272 + Lstmp296 + Lstmp299 + Lstmp6*L[101] + L[52];
#pragma omp atomic
Ls[53] += Lstmp118 + Lstmp15*L[115] + Lstmp24*L[117] + Lstmp243 + Lstmp248 + Lstmp314 + Lstmp315*y + Lstmp317 + Lstmp6*L[102] + L[53];
#pragma omp atomic
Ls[54] += Lstmp15*L[116] + Lstmp158 + Lstmp166 + Lstmp173 + Lstmp24*L[118] + Lstmp334 + Lstmp335*x + Lstmp338*y + Lstmp6*L[103] + L[54];
#pragma omp atomic
Ls[55] += Lstmp15*L[117] + Lstmp24*L[119] + Lstmp341*x + Lstmp342*y + Lstmp40 + Lstmp55 + Lstmp6*L[104] + Lstmp69 + z*L[83] + L[55];
#pragma omp atomic
Ls[56] += Lstmp50 + Lstmp51 + x*L[84] + L[56];
#pragma omp atomic
Ls[57] += Lstmp163 + Lstmp84 + x*L[85] + L[57];
#pragma omp atomic
Ls[58] += Lstmp195 + Lstmp196 + x*L[86] + L[58];
#pragma omp atomic
Ls[59] += Lstmp134 + Lstmp241 + Lstmp320 + L[59];
#pragma omp atomic
Ls[60] += Lstmp211 + Lstmp256 + x*L[88] + L[60];
#pragma omp atomic
Ls[61] += Lstmp133 + Lstmp265 + Lstmp321 + L[61];
#pragma omp atomic
Ls[62] += Lstmp132 + Lstmp279 + Lstmp295 + L[62];
#pragma omp atomic
Ls[63] += Lstmp215 + Lstmp304 + Lstmp323 + L[63];
#pragma omp atomic
Ls[64] += Lstmp179 + Lstmp275 + Lstmp322 + L[64];
#pragma omp atomic
Ls[65] += Lstmp130 + Lstmp281 + Lstmp313 + L[65];
#pragma omp atomic
Ls[66] += Lstmp128 + Lstmp218 + Lstmp326 + L[66];
#pragma omp atomic
Ls[67] += Lstmp214 + Lstmp288 + Lstmp330 + L[67];
#pragma omp atomic
Ls[68] += Lstmp226 + Lstmp251 + Lstmp277 + L[68];
#pragma omp atomic
Ls[69] += Lstmp178 + Lstmp285 + Lstmp318 + L[69];
#pragma omp atomic
Ls[70] += Lstmp125 + Lstmp221 + Lstmp337 + L[70];
#pragma omp atomic
Ls[71] += Lstmp140 + Lstmp339 + Lstmp76 + L[71];
#pragma omp atomic
Ls[72] += Lstmp205 + Lstmp234 + Lstmp340 + L[72];
#pragma omp atomic
Ls[73] += Lstmp153 + Lstmp271 + Lstmp298 + L[73];
#pragma omp atomic
Ls[74] += Lstmp152 + Lstmp247 + Lstmp316 + L[74];
#pragma omp atomic
Ls[75] += Lstmp172 + Lstmp229 + Lstmp335 + L[75];
#pragma omp atomic
Ls[76] += Lstmp144 + Lstmp341 + Lstmp68 + L[76];
#pragma omp atomic
Ls[77] += Lstmp36 + Lstmp61 + y*L[112] + L[77];
#pragma omp atomic
Ls[78] += Lstmp188 + Lstmp200 + y*L[113] + L[78];
#pragma omp atomic
Ls[79] += Lstmp122 + Lstmp268 + Lstmp327 + L[79];
#pragma omp atomic
Ls[80] += Lstmp121 + Lstmp297 + Lstmp315 + L[80];
#pragma omp atomic
Ls[81] += Lstmp119 + Lstmp244 + Lstmp338 + L[81];
#pragma omp atomic
Ls[82] += Lstmp159 + Lstmp167 + Lstmp342 + L[82];
#pragma omp atomic
Ls[83] += Lstmp41 + Lstmp56 + z*L[119] + L[83];
#pragma omp atomic
Ls[84] += L[84];
#pragma omp atomic
Ls[85] += L[85];
#pragma omp atomic
Ls[86] += L[86];
#pragma omp atomic
Ls[87] += L[87];
#pragma omp atomic
Ls[88] += L[88];
#pragma omp atomic
Ls[89] += L[89];
#pragma omp atomic
Ls[90] += L[90];
#pragma omp atomic
Ls[91] += L[91];
#pragma omp atomic
Ls[92] += L[92];
#pragma omp atomic
Ls[93] += L[93];
#pragma omp atomic
Ls[94] += L[94];
#pragma omp atomic
Ls[95] += L[95];
#pragma omp atomic
Ls[96] += L[96];
#pragma omp atomic
Ls[97] += L[97];
#pragma omp atomic
Ls[98] += L[98];
#pragma omp atomic
Ls[99] += L[99];
#pragma omp atomic
Ls[100] += L[100];
#pragma omp atomic
Ls[101] += L[101];
#pragma omp atomic
Ls[102] += L[102];
#pragma omp atomic
Ls[103] += L[103];
#pragma omp atomic
Ls[104] += L[104];
#pragma omp atomic
Ls[105] += L[105];
#pragma omp atomic
Ls[106] += L[106];
#pragma omp atomic
Ls[107] += L[107];
#pragma omp atomic
Ls[108] += L[108];
#pragma omp atomic
Ls[109] += L[109];
#pragma omp atomic
Ls[110] += L[110];
#pragma omp atomic
Ls[111] += L[111];
#pragma omp atomic
Ls[112] += L[112];
#pragma omp atomic
Ls[113] += L[113];
#pragma omp atomic
Ls[114] += L[114];
#pragma omp atomic
Ls[115] += L[115];
#pragma omp atomic
Ls[116] += L[116];
#pragma omp atomic
Ls[117] += L[117];
#pragma omp atomic
Ls[118] += L[118];
#pragma omp atomic
Ls[119] += L[119];

}

void field_m0_L2P_7(double x, double y, double z, double * L, double * F) {
double Ftmp0 = x*y;
double Ftmp1 = x*z;
double Ftmp2 = y*z;
double Ftmp3 = Ftmp0*z;
double Ftmp4 = (x*x);
double Ftmp5 = (1.0/2.0)*Ftmp4;
double Ftmp6 = (x*x*x);
double Ftmp7 = (1.0/6.0)*Ftmp6;
double Ftmp8 = (x*x*x*x);
double Ftmp9 = (1.0/24.0)*Ftmp8;
double Ftmp10 = (1.0/120.0)*(x*x*x*x*x);
double Ftmp11 = (1.0/720.0)*(x*x*x*x*x*x);
double Ftmp12 = (y*y);
double Ftmp13 = (1.0/2.0)*Ftmp12;
double Ftmp14 = (y*y*y);
double Ftmp15 = (1.0/6.0)*Ftmp14;
double Ftmp16 = (y*y*y*y);
double Ftmp17 = (1.0/24.0)*Ftmp16;
double Ftmp18 = (1.0/120.0)*(y*y*y*y*y);
double Ftmp19 = (1.0/720.0)*(y*y*y*y*y*y);
double Ftmp20 = (z*z);
double Ftmp21 = (1.0/2.0)*Ftmp20;
double Ftmp22 = (z*z*z);
double Ftmp23 = (1.0/6.0)*Ftmp22;
double Ftmp24 = (z*z*z*z);
double Ftmp25 = (1.0/24.0)*Ftmp24;
double Ftmp26 = (1.0/120.0)*(z*z*z*z*z);
double Ftmp27 = (1.0/720.0)*(z*z*z*z*z*z);
double Ftmp28 = Ftmp13*x;
double Ftmp29 = Ftmp15*x;
double Ftmp30 = Ftmp17*x;
double Ftmp31 = Ftmp18*x;
double Ftmp32 = Ftmp21*x;
double Ftmp33 = Ftmp23*x;
double Ftmp34 = Ftmp25*x;
double Ftmp35 = Ftmp26*x;
double Ftmp36 = Ftmp5*y;
double Ftmp37 = Ftmp5*z;
double Ftmp38 = Ftmp7*y;
double Ftmp39 = Ftmp7*z;
double Ftmp40 = Ftmp9*y;
double Ftmp41 = Ftmp9*z;
double Ftmp42 = Ftmp10*y;
double Ftmp43 = Ftmp10*z;
double Ftmp44 = Ftmp21*y;
double Ftmp45 = Ftmp23*y;
double Ftmp46 = Ftmp25*y;
double Ftmp47 = Ftmp26*y;
double Ftmp48 = Ftmp13*z;
double Ftmp49 = Ftmp15*z;
double Ftmp50 = Ftmp17*z;
double Ftmp51 = Ftmp18*z;
double Ftmp52 = Ftmp0*Ftmp21;
double Ftmp53 = Ftmp0*Ftmp23;
double Ftmp54 = Ftmp0*Ftmp25;
double Ftmp55 = Ftmp1*Ftmp13;
double Ftmp56 = Ftmp1*Ftmp15;
double Ftmp57 = Ftmp1*Ftmp17;
double Ftmp58 = Ftmp2*Ftmp5;
double Ftmp59 = Ftmp2*Ftmp7;
double Ftmp60 = Ftmp2*Ftmp9;
double Ftmp61 = (1.0/4.0)*Ftmp4;
double Ftmp62 = Ftmp12*Ftmp61;
double Ftmp63 = (1.0/12.0)*Ftmp4;
double Ftmp64 = Ftmp14*Ftmp63;
double Ftmp65 = (1.0/48.0)*Ftmp4;
double Ftmp66 = Ftmp16*Ftmp65;
double Ftmp67 = Ftmp20*Ftmp61;
double Ftmp68 = Ftmp22*Ftmp63;
double Ftmp69 = Ftmp24*Ftmp65;
double Ftmp70 = (1.0/12.0)*Ftmp6;
double Ftmp71 = Ftmp12*Ftmp70;
double Ftmp72 = (1.0/36.0)*Ftmp6;
double Ftmp73 = Ftmp14*Ftmp72;
double Ftmp74 = Ftmp20*Ftmp70;
double Ftmp75 = Ftmp22*Ftmp72;
double Ftmp76 = (1.0/48.0)*Ftmp8;
double Ftmp77 = Ftmp12*Ftmp76;
double Ftmp78 = Ftmp20*Ftmp76;
double Ftmp79 = Ftmp12*Ftmp20;
double Ftmp80 = (1.0/4.0)*Ftmp79;
double Ftmp81 = (1.0/12.0)*Ftmp12*Ftmp22;
double Ftmp82 = (1.0/48.0)*Ftmp12*Ftmp24;
double Ftmp83 = (1.0/12.0)*Ftmp14*Ftmp20;
double Ftmp84 = (1.0/36.0)*Ftmp14*Ftmp22;
double Ftmp85 = (1.0/48.0)*Ftmp16*Ftmp20;
double Ftmp86 = Ftmp80*x;
double Ftmp87 = Ftmp81*x;
double Ftmp88 = Ftmp83*x;
double Ftmp89 = Ftmp67*y;
double Ftmp90 = Ftmp68*y;
double Ftmp91 = Ftmp62*z;
double Ftmp92 = Ftmp64*z;
double Ftmp93 = Ftmp74*y;
double Ftmp94 = Ftmp71*z;
double Ftmp95 = (1.0/8.0)*Ftmp4*Ftmp79;
#pragma omp atomic
F[0] += -Ftmp0*L[11] - Ftmp1*L[12] - Ftmp10*L[56] - Ftmp11*L[84] - Ftmp13*L[13] - Ftmp15*L[26] - Ftmp17*L[45] - Ftmp18*L[71] - Ftmp19*L[105] - Ftmp2*L[14] - Ftmp21*L[15] - Ftmp23*L[29] - Ftmp25*L[49] - Ftmp26*L[76] - Ftmp27*L[111] - Ftmp28*L[23] - Ftmp29*L[41] - Ftmp3*L[24] - Ftmp30*L[66] - Ftmp31*L[99] - Ftmp32*L[25] - Ftmp33*L[44] - Ftmp34*L[70] - Ftmp35*L[104] - Ftmp36*L[21] - Ftmp37*L[22] - Ftmp38*L[36] - Ftmp39*L[37] - Ftmp40*L[57] - Ftmp41*L[58] - Ftmp42*L[85] - Ftmp43*L[86] - Ftmp44*L[28] - Ftmp45*L[48] - Ftmp46*L[75] - Ftmp47*L[110] - Ftmp48*L[27] - Ftmp49*L[46] - Ftmp5*L[10] - Ftmp50*L[72] - Ftmp51*L[106] - Ftmp52*L[43] - Ftmp53*L[69] - Ftmp54*L[103] - Ftmp55*L[42] - Ftmp56*L[67] - Ftmp57*L[100] - Ftmp58*L[39] - Ftmp59*L[60] - Ftmp60*L[88] - Ftmp62*L[38] - Ftmp64*L[62] - Ftmp66*L[94] - Ftmp67*L[40] - Ftmp68*L[65] - Ftmp69*L[98] - Ftmp7*L[20] - Ftmp71*L[59] - Ftmp73*L[90] - Ftmp74*L[61] - Ftmp75*L[93] - Ftmp77*L[87] - Ftmp78*L[89] - Ftmp80*L[47] - Ftmp81*L[74] - Ftmp82*L[109] - Ftmp83*L[73] - Ftmp84*L[108] - Ftmp85*L[107] - Ftmp86*L[68] - Ftmp87*L[102] - Ftmp88*L[101] - Ftmp89*L[64] - Ftmp9*L[35] - Ftmp90*L[97] - Ftmp91*L[63] - Ftmp92*L[95] - Ftmp93*L[92] - Ftmp94*L[91] - Ftmp95*L[96] - x*L[4] - y*L[5] - z*L[6] - L[1];
#pragma omp atomic
F[1] += -Ftmp0*L[13] - Ftmp1*L[14] - Ftmp10*L[57] - Ftmp11*L[85] - Ftmp13*L[16] - Ftmp15*L[30] - Ftmp17*L[50] - Ftmp18*L[77] - Ftmp19*L[112] - Ftmp2*L[17] - Ftmp21*L[18] - Ftmp23*L[33] - Ftmp25*L[54] - Ftmp26*L[82] - Ftmp27*L[118] - Ftmp28*L[26] - Ftmp29*L[45] - Ftmp3*L[27] - Ftmp30*L[71] - Ftmp31*L[105] - Ftmp32*L[28] - Ftmp33*L[48] - Ftmp34*L[75] - Ftmp35*L[110] - Ftmp36*L[23] - Ftmp37*L[24] - Ftmp38*L[38] - Ftmp39*L[39] - Ftmp40*L[59] - Ftmp41*L[60] - Ftmp42*L[87] - Ftmp43*L[88] - Ftmp44*L[32] - Ftmp45*L[53] - Ftmp46*L[81] - Ftmp47*L[117] - Ftmp48*L[31] - Ftmp49*L[51] - Ftmp5*L[11] - Ftmp50*L[78] - Ftmp51*L[113] - Ftmp52*L[47] - Ftmp53*L[74] - Ftmp54*L[109] - Ftmp55*L[46] - Ftmp56*L[72] - Ftmp57*L[106] - Ftmp58*L[42] - Ftmp59*L[63] - Ftmp60*L[91] - Ftmp62*L[41] - Ftmp64*L[66] - Ftmp66*L[99] - Ftmp67*L[43] - Ftmp68*L[69] - Ftmp69*L[103] - Ftmp7*L[21] - Ftmp71*L[62] - Ftmp73*L[94] - Ftmp74*L[64] - Ftmp75*L[97] - Ftmp77*L[90] - Ftmp78*L[92] - Ftmp80*L[52] - Ftmp81*L[80] - Ftmp82*L[116] - Ftmp83*L[79] - Ftmp84*L[115] - Ftmp85*L[114] - Ftmp86*L[73] - Ftmp87*L[108] - Ftmp88*L[107] - Ftmp89*L[68] - Ftmp9*L[36] - Ftmp90*L[102] - Ftmp91*L[67] - Ftmp92*L[100] - Ftmp93*L[96] - Ftmp94*L[95] - Ftmp95*L[101] - x*L[5] - y*L[7] - z*L[8] - L[2];
#pragma omp atomic
F[2] += -Ftmp0*L[14] - Ftmp1*L[15] - Ftmp10*L[58] - Ftmp11*L[86] - Ftmp13*L[17] - Ftmp15*L[31] - Ftmp17*L[51] - Ftmp18*L[78] - Ftmp19*L[113] - Ftmp2*L[18] - Ftmp21*L[19] - Ftmp23*L[34] - Ftmp25*L[55] - Ftmp26*L[83] - Ftmp27*L[119] - Ftmp28*L[27] - Ftmp29*L[46] - Ftmp3*L[28] - Ftmp30*L[72] - Ftmp31*L[106] - Ftmp32*L[29] - Ftmp33*L[49] - Ftmp34*L[76] - Ftmp35*L[111] - Ftmp36*L[24] - Ftmp37*L[25] - Ftmp38*L[39] - Ftmp39*L[40] - Ftmp40*L[60] - Ftmp41*L[61] - Ftmp42*L[88] - Ftmp43*L[89] - Ftmp44*L[33] - Ftmp45*L[54] - Ftmp46*L[82] - Ftmp47*L[118] - Ftmp48*L[32] - Ftmp49*L[52] - Ftmp5*L[12] - Ftmp50*L[79] - Ftmp51*L[114] - Ftmp52*L[48] - Ftmp53*L[75] - Ftmp54*L[110] - Ftmp55*L[47] - Ftmp56*L[73] - Ftmp57*L[107] - Ftmp58*L[43] - Ftmp59*L[64] - Ftmp60*L[92] - Ftmp62*L[42] - Ftmp64*L[67] - Ftmp66*L[100] - Ftmp67*L[44] - Ftmp68*L[70] - Ftmp69*L[104] - Ftmp7*L[22] - Ftmp71*L[63] - Ftmp73*L[95] - Ftmp74*L[65] - Ftmp75*L[98] - Ftmp77*L[91] - Ftmp78*L[93] - Ftmp80*L[53] - Ftmp81*L[81] - Ftmp82*L[117] - Ftmp83*L[80] - Ftmp84*L[116] - Ftmp85*L[115] - Ftmp86*L[74] - Ftmp87*L[109] - Ftmp88*L[108] - Ftmp89*L[69] - Ftmp9*L[37] - Ftmp90*L[103] - Ftmp91*L[68] - Ftmp92*L[101] - Ftmp93*L[97] - Ftmp94*L[96] - Ftmp95*L[102] - x*L[6] - y*L[8] - z*L[9] - L[3];

}

void field_m0_M2P_7(double x, double y, double z, double * M, double * F) {
double R = sqrt(x*x + y*y + z*z);
double Ftmp0 = 1.0*pow(R, -3.0);
double Ftmp1 = pow(R, -5.0);
double Ftmp2 = 3.0*Ftmp1;
double Ftmp3 = Ftmp2*M[5];
double Ftmp4 = Ftmp2*z;
double Ftmp5 = Ftmp0*M[0];
double Ftmp6 = y*z;
double Ftmp7 = pow(R, -7.0);
double Ftmp8 = 15.0*Ftmp7;
double Ftmp9 = Ftmp8*M[14];
double Ftmp10 = x*y;
double Ftmp11 = Ftmp2*M[2];
double Ftmp12 = Ftmp4*M[3];
double Ftmp13 = (x*x);
double Ftmp14 = Ftmp2*M[1];
double Ftmp15 = Ftmp6*x;
double Ftmp16 = Ftmp15*Ftmp8;
double Ftmp17 = Ftmp13*Ftmp8;
double Ftmp18 = pow(R, -9.0);
double Ftmp19 = 105.0*Ftmp18;
double Ftmp20 = Ftmp13*Ftmp19;
double Ftmp21 = -9.0*Ftmp1;
double Ftmp22 = Ftmp17 + Ftmp21;
double Ftmp23 = -Ftmp2;
double Ftmp24 = (y*y);
double Ftmp25 = Ftmp24*Ftmp8;
double Ftmp26 = Ftmp23 + Ftmp25;
double Ftmp27 = (z*z);
double Ftmp28 = Ftmp27*Ftmp8;
double Ftmp29 = Ftmp23 + Ftmp28;
double Ftmp30 = Ftmp26*M[7];
double Ftmp31 = Ftmp29*M[9];
double Ftmp32 = 45.0*Ftmp7;
double Ftmp33 = -Ftmp32;
double Ftmp34 = Ftmp20 + Ftmp33;
double Ftmp35 = Ftmp34*M[21];
double Ftmp36 = 1.0*y;
double Ftmp37 = Ftmp19*Ftmp24;
double Ftmp38 = Ftmp33 + Ftmp37;
double Ftmp39 = Ftmp38*M[26];
double Ftmp40 = 3.0*y;
double Ftmp41 = 35.0*Ftmp18;
double Ftmp42 = (Ftmp27*Ftmp41 - 5.0*Ftmp7)*M[28];
double Ftmp43 = Ftmp34*M[22];
double Ftmp44 = 1.0*z;
double Ftmp45 = -Ftmp8;
double Ftmp46 = Ftmp37 + Ftmp45;
double Ftmp47 = Ftmp46*M[27];
double Ftmp48 = Ftmp19*Ftmp27;
double Ftmp49 = Ftmp33 + Ftmp48;
double Ftmp50 = Ftmp44*Ftmp49;
double Ftmp51 = 315.0*Ftmp18;
double Ftmp52 = -Ftmp51;
double Ftmp53 = pow(R, -11.0);
double Ftmp54 = 945.0*Ftmp53;
double Ftmp55 = Ftmp13*Ftmp54;
double Ftmp56 = Ftmp52 + Ftmp55;
double Ftmp57 = Ftmp56*M[39];
double Ftmp58 = Ftmp10*Ftmp34;
double Ftmp59 = Ftmp38*M[16];
double Ftmp60 = Ftmp45 + Ftmp48;
double Ftmp61 = Ftmp60*M[18];
double Ftmp62 = Ftmp36*x;
double Ftmp63 = x*z;
double Ftmp64 = Ftmp34*Ftmp63;
double Ftmp65 = Ftmp46*M[17];
double Ftmp66 = Ftmp49*M[19];
double Ftmp67 = Ftmp36*z;
double Ftmp68 = Ftmp24*Ftmp54;
double Ftmp69 = Ftmp52 + Ftmp68;
double Ftmp70 = Ftmp69*M[46];
double Ftmp71 = Ftmp40*z;
double Ftmp72 = -Ftmp19;
double Ftmp73 = Ftmp27*Ftmp53;
double Ftmp74 = 315.0*Ftmp73;
double Ftmp75 = Ftmp72 + Ftmp74;
double Ftmp76 = Ftmp75*M[48];
double Ftmp77 = Ftmp69*M[31];
double Ftmp78 = -75.0*Ftmp7;
double Ftmp79 = 1.0*Ftmp13;
double Ftmp80 = Ftmp46*M[13];
double Ftmp81 = Ftmp60*M[15];
double Ftmp82 = 525.0*Ftmp18;
double Ftmp83 = -Ftmp82;
double Ftmp84 = Ftmp55 + Ftmp83;
double Ftmp85 = Ftmp13*y;
double Ftmp86 = Ftmp13*z;
double Ftmp87 = Ftmp27*Ftmp54;
double Ftmp88 = Ftmp52 + Ftmp87;
double Ftmp89 = Ftmp36*M[33];
double Ftmp90 = Ftmp13*Ftmp36;
double Ftmp91 = Ftmp69*M[26];
double Ftmp92 = Ftmp13*Ftmp40;
double Ftmp93 = (-Ftmp41 + Ftmp74)*M[28];
double Ftmp94 = Ftmp13*Ftmp44;
double Ftmp95 = (Ftmp68 + Ftmp72)*M[27];
double Ftmp96 = Ftmp88*M[29];
double Ftmp97 = 4725.0*Ftmp53;
double Ftmp98 = -Ftmp97;
double Ftmp99 = pow(R, -13.0);
double Ftmp100 = 10395.0*Ftmp99;
double Ftmp101 = Ftmp100*Ftmp13;
double Ftmp102 = Ftmp101 + Ftmp98;
double Ftmp103 = Ftmp13*Ftmp6;
double Ftmp104 = Ftmp36*Ftmp86;
double Ftmp105 = 2835.0*Ftmp53;
double Ftmp106 = -Ftmp105;
double Ftmp107 = Ftmp100*Ftmp24;
double Ftmp108 = (Ftmp106 + Ftmp107)*M[46];
double Ftmp109 = Ftmp40*Ftmp86;
double Ftmp110 = 3465.0*Ftmp99;
double Ftmp111 = Ftmp110*Ftmp27;
double Ftmp112 = (Ftmp111 - Ftmp54)*M[48];
double Ftmp113 = 225.0*Ftmp7;
double Ftmp114 = (x*x*x*x);
double Ftmp115 = Ftmp114*Ftmp54;
double Ftmp116 = 1050.0*Ftmp18;
double Ftmp117 = Ftmp113 + Ftmp115 - Ftmp116*Ftmp13;
double Ftmp118 = (y*y*y*y);
double Ftmp119 = Ftmp118*Ftmp54;
double Ftmp120 = 630.0*Ftmp18;
double Ftmp121 = Ftmp119 - Ftmp120*Ftmp24 + Ftmp32;
double Ftmp122 = (z*z*z*z);
double Ftmp123 = Ftmp122*Ftmp54;
double Ftmp124 = -Ftmp120*Ftmp27 + Ftmp123 + Ftmp32;
double Ftmp125 = Ftmp121*M[30];
double Ftmp126 = Ftmp124*M[34];
double Ftmp127 = 1575.0*Ftmp18;
double Ftmp128 = Ftmp114*Ftmp99;
double Ftmp129 = 10395.0*Ftmp128;
double Ftmp130 = Ftmp13*Ftmp53;
double Ftmp131 = 9450.0*Ftmp130;
double Ftmp132 = Ftmp127 + Ftmp129 - Ftmp131;
double Ftmp133 = Ftmp132*M[57];
double Ftmp134 = Ftmp100*Ftmp118;
double Ftmp135 = 9450.0*Ftmp53;
double Ftmp136 = Ftmp135*Ftmp24;
double Ftmp137 = Ftmp127 + Ftmp134 - Ftmp136;
double Ftmp138 = Ftmp137*M[71];
double Ftmp139 = (Ftmp110*Ftmp122 + Ftmp19 - 1890.0*Ftmp73)*M[75];
double Ftmp140 = Ftmp132*M[58];
double Ftmp141 = 5670.0*Ftmp53;
double Ftmp142 = Ftmp141*Ftmp24;
double Ftmp143 = Ftmp134 - Ftmp142 + Ftmp51;
double Ftmp144 = Ftmp143*M[72];
double Ftmp145 = Ftmp100*Ftmp122;
double Ftmp146 = Ftmp135*Ftmp27;
double Ftmp147 = Ftmp127 + Ftmp145 - Ftmp146;
double Ftmp148 = Ftmp147*Ftmp44;
double Ftmp149 = 14175.0*Ftmp53;
double Ftmp150 = pow(R, -15.0);
double Ftmp151 = 135135.0*Ftmp150;
double Ftmp152 = Ftmp114*Ftmp151;
double Ftmp153 = Ftmp13*Ftmp99;
double Ftmp154 = 103950.0*Ftmp153;
double Ftmp155 = Ftmp149 + Ftmp152 - Ftmp154;
double Ftmp156 = Ftmp155*M[88];
double Ftmp157 = Ftmp10*Ftmp132;
double Ftmp158 = Ftmp137*M[50];
double Ftmp159 = Ftmp141*Ftmp27;
double Ftmp160 = Ftmp145 - Ftmp159 + Ftmp51;
double Ftmp161 = Ftmp160*M[54];
double Ftmp162 = Ftmp132*Ftmp63;
double Ftmp163 = Ftmp143*M[51];
double Ftmp164 = Ftmp147*M[55];
double Ftmp165 = Ftmp118*Ftmp151;
double Ftmp166 = Ftmp24*Ftmp99;
double Ftmp167 = 103950.0*Ftmp166;
double Ftmp168 = Ftmp149 + Ftmp165 - Ftmp167;
double Ftmp169 = Ftmp168*M[106];
double Ftmp170 = Ftmp122*Ftmp150;
double Ftmp171 = 45045.0*Ftmp170;
double Ftmp172 = Ftmp27*Ftmp99;
double Ftmp173 = Ftmp171 - 34650.0*Ftmp172 + Ftmp97;
double Ftmp174 = Ftmp173*M[110];
double Ftmp175 = Ftmp168*M[78];
double Ftmp176 = 3675.0*Ftmp18;
double Ftmp177 = Ftmp143*M[45];
double Ftmp178 = Ftmp160*M[49];
double Ftmp179 = 33075.0*Ftmp53;
double Ftmp180 = 145530.0*Ftmp153;
double Ftmp181 = Ftmp152 + Ftmp179 - Ftmp180;
double Ftmp182 = Ftmp122*Ftmp151;
double Ftmp183 = 103950.0*Ftmp172;
double Ftmp184 = Ftmp149 + Ftmp182 - Ftmp183;
double Ftmp185 = Ftmp36*M[82];
double Ftmp186 = Ftmp168*M[71];
double Ftmp187 = (Ftmp171 - 20790.0*Ftmp172 + Ftmp54)*M[75];
double Ftmp188 = 62370.0*Ftmp99;
double Ftmp189 = Ftmp188*Ftmp24;
double Ftmp190 = (Ftmp105 + Ftmp165 - Ftmp189)*M[72];
double Ftmp191 = Ftmp184*M[76];
double Ftmp192 = 363825.0*Ftmp99;
double Ftmp193 = pow(R, -17.0);
double Ftmp194 = 2027025.0*Ftmp193;
double Ftmp195 = Ftmp114*Ftmp194;
double Ftmp196 = Ftmp13*Ftmp150;
double Ftmp197 = 1891890.0*Ftmp196;
double Ftmp198 = 155925.0*Ftmp99;
double Ftmp199 = Ftmp118*Ftmp194;
double Ftmp200 = Ftmp150*Ftmp24;
double Ftmp201 = 1351350.0*Ftmp200;
double Ftmp202 = (Ftmp198 + Ftmp199 - Ftmp201)*M[106];
double Ftmp203 = 51975.0*Ftmp99;
double Ftmp204 = 675675.0*Ftmp122*Ftmp193;
double Ftmp205 = Ftmp150*Ftmp27;
double Ftmp206 = (Ftmp203 + Ftmp204 - 450450.0*Ftmp205)*M[110];
double Ftmp207 = -11025.0*Ftmp18;
double Ftmp208 = (x*x*x*x*x*x);
double Ftmp209 = Ftmp151*Ftmp208;
double Ftmp210 = 99225.0*Ftmp53;
double Ftmp211 = -218295.0*Ftmp128 + Ftmp13*Ftmp210 + Ftmp207 + Ftmp209;
double Ftmp212 = -Ftmp127;
double Ftmp213 = (y*y*y*y*y*y);
double Ftmp214 = Ftmp151*Ftmp213;
double Ftmp215 = 42525.0*Ftmp53;
double Ftmp216 = -Ftmp118*Ftmp198 + Ftmp212 + Ftmp214 + Ftmp215*Ftmp24;
double Ftmp217 = (z*z*z*z*z*z);
double Ftmp218 = Ftmp151*Ftmp217;
double Ftmp219 = -Ftmp122*Ftmp198 + Ftmp212 + Ftmp215*Ftmp27 + Ftmp218;
double Ftmp220 = Ftmp216*M[77];
double Ftmp221 = Ftmp219*M[83];
double Ftmp222 = -Ftmp24*Ftmp51;
double Ftmp223 = Ftmp24*Ftmp55;
double Ftmp224 = -Ftmp20;
double Ftmp225 = Ftmp224 + Ftmp32;
double Ftmp226 = Ftmp222 + Ftmp223 + Ftmp225;
double Ftmp227 = -Ftmp27*Ftmp51;
double Ftmp228 = Ftmp27*Ftmp55;
double Ftmp229 = Ftmp225 + Ftmp227 + Ftmp228;
double Ftmp230 = -Ftmp48;
double Ftmp231 = Ftmp230 + Ftmp8;
double Ftmp232 = -Ftmp37;
double Ftmp233 = Ftmp27*Ftmp68;
double Ftmp234 = Ftmp232 + Ftmp233;
double Ftmp235 = Ftmp231 + Ftmp234;
double Ftmp236 = -Ftmp210;
double Ftmp237 = Ftmp194*Ftmp208;
double Ftmp238 = Ftmp114*Ftmp150;
double Ftmp239 = 1091475.0*Ftmp153 + Ftmp236 + Ftmp237 - 2837835.0*Ftmp238;
double Ftmp240 = Ftmp10*Ftmp239;
double Ftmp241 = Ftmp194*Ftmp213;
double Ftmp242 = Ftmp118*Ftmp150;
double Ftmp243 = 1091475.0*Ftmp166 + Ftmp236 + Ftmp241 - 2837835.0*Ftmp242;
double Ftmp244 = Ftmp243*M[112];
double Ftmp245 = -Ftmp149;
double Ftmp246 = Ftmp194*Ftmp217;
double Ftmp247 = 2027025.0*Ftmp150;
double Ftmp248 = 467775.0*Ftmp99;
double Ftmp249 = -Ftmp122*Ftmp247 + Ftmp245 + Ftmp246 + Ftmp248*Ftmp27;
double Ftmp250 = Ftmp249*M[118];
double Ftmp251 = Ftmp239*Ftmp63;
double Ftmp252 = -Ftmp118*Ftmp247 + Ftmp24*Ftmp248 + Ftmp241 + Ftmp245;
double Ftmp253 = Ftmp252*M[113];
double Ftmp254 = -2837835.0*Ftmp170 + 1091475.0*Ftmp172 + Ftmp236 + Ftmp246;
double Ftmp255 = Ftmp254*M[119];
double Ftmp256 = -297675.0*Ftmp53;
double Ftmp257 = Ftmp252*M[105];
double Ftmp258 = Ftmp249*M[111];
double Ftmp259 = Ftmp105*Ftmp24;
double Ftmp260 = -Ftmp259;
double Ftmp261 = Ftmp101*Ftmp24;
double Ftmp262 = Ftmp260 + Ftmp261;
double Ftmp263 = 945.0*Ftmp18;
double Ftmp264 = Ftmp105*Ftmp13;
double Ftmp265 = -Ftmp264;
double Ftmp266 = Ftmp263 + Ftmp265;
double Ftmp267 = Ftmp262 + Ftmp266;
double Ftmp268 = Ftmp267*M[62];
double Ftmp269 = -Ftmp55;
double Ftmp270 = Ftmp269 + Ftmp51;
double Ftmp271 = Ftmp105*Ftmp27;
double Ftmp272 = -Ftmp271;
double Ftmp273 = Ftmp101*Ftmp27;
double Ftmp274 = Ftmp272 + Ftmp273;
double Ftmp275 = Ftmp270 + Ftmp274;
double Ftmp276 = Ftmp275*M[64];
double Ftmp277 = -Ftmp68;
double Ftmp278 = Ftmp107*Ftmp27;
double Ftmp279 = Ftmp278 + Ftmp51;
double Ftmp280 = Ftmp272 + Ftmp277 + Ftmp279;
double Ftmp281 = Ftmp280*M[73];
double Ftmp282 = Ftmp262 + Ftmp270;
double Ftmp283 = Ftmp282*M[63];
double Ftmp284 = Ftmp266 + Ftmp274;
double Ftmp285 = Ftmp284*M[65];
double Ftmp286 = -Ftmp87;
double Ftmp287 = Ftmp260 + Ftmp279 + Ftmp286;
double Ftmp288 = Ftmp287*M[74];
double Ftmp289 = 8505.0*Ftmp53;
double Ftmp290 = 31185.0*Ftmp99;
double Ftmp291 = Ftmp13*Ftmp290;
double Ftmp292 = -Ftmp291;
double Ftmp293 = Ftmp289 + Ftmp292;
double Ftmp294 = Ftmp24*Ftmp290;
double Ftmp295 = -Ftmp294;
double Ftmp296 = Ftmp13*Ftmp151;
double Ftmp297 = Ftmp24*Ftmp296;
double Ftmp298 = Ftmp295 + Ftmp297;
double Ftmp299 = Ftmp293 + Ftmp298;
double Ftmp300 = Ftmp299*M[95];
double Ftmp301 = Ftmp27*Ftmp290;
double Ftmp302 = -Ftmp301;
double Ftmp303 = Ftmp27*Ftmp296;
double Ftmp304 = Ftmp302 + Ftmp303;
double Ftmp305 = Ftmp293 + Ftmp304;
double Ftmp306 = Ftmp305*M[97];
double Ftmp307 = Ftmp10*Ftmp267;
double Ftmp308 = Ftmp10*Ftmp275;
double Ftmp309 = Ftmp282*Ftmp63;
double Ftmp310 = Ftmp284*Ftmp63;
double Ftmp311 = Ftmp24*Ftmp27;
double Ftmp312 = Ftmp151*Ftmp311;
double Ftmp313 = Ftmp302 + Ftmp312;
double Ftmp314 = Ftmp289 + Ftmp295 + Ftmp313;
double Ftmp315 = Ftmp314*M[108];
double Ftmp316 = Ftmp15*Ftmp299;
double Ftmp317 = Ftmp15*Ftmp305;
double Ftmp318 = -Ftmp24*Ftmp97;
double Ftmp319 = Ftmp269 + Ftmp82;
double Ftmp320 = -Ftmp27*Ftmp97;
double Ftmp321 = Ftmp19 + Ftmp286;
double Ftmp322 = Ftmp277 + Ftmp278;
double Ftmp323 = Ftmp203*Ftmp24;
double Ftmp324 = -Ftmp323;
double Ftmp325 = Ftmp297 + Ftmp324;
double Ftmp326 = Ftmp149 + Ftmp292;
double Ftmp327 = -Ftmp101;
double Ftmp328 = Ftmp327 + Ftmp97;
double Ftmp329 = Ftmp203*Ftmp27;
double Ftmp330 = -Ftmp329;
double Ftmp331 = Ftmp303 + Ftmp330;
double Ftmp332 = -Ftmp107;
double Ftmp333 = Ftmp105 + Ftmp332;
double Ftmp334 = Ftmp100*Ftmp27;
double Ftmp335 = -Ftmp334;
double Ftmp336 = Ftmp105 + Ftmp335;
double Ftmp337 = Ftmp295 + Ftmp312;
double Ftmp338 = 675675.0*Ftmp150;
double Ftmp339 = Ftmp24*Ftmp338;
double Ftmp340 = -Ftmp339;
double Ftmp341 = Ftmp13*Ftmp194;
double Ftmp342 = Ftmp24*Ftmp341;
double Ftmp343 = 405405.0*Ftmp196;
double Ftmp344 = -Ftmp343;
double Ftmp345 = Ftmp198 + Ftmp344;
double Ftmp346 = Ftmp27*Ftmp338;
double Ftmp347 = -Ftmp346;
double Ftmp348 = Ftmp27*Ftmp341;
double Ftmp349 = 93555.0*Ftmp99;
double Ftmp350 = -405405.0*Ftmp205;
double Ftmp351 = Ftmp349 + Ftmp350;
double Ftmp352 = 405405.0*Ftmp200;
double Ftmp353 = -Ftmp352;
double Ftmp354 = Ftmp194*Ftmp311;
double Ftmp355 = Ftmp353 + Ftmp354;
double Ftmp356 = -Ftmp263;
double Ftmp357 = Ftmp264 + Ftmp356;
double Ftmp358 = Ftmp13*Ftmp165;
double Ftmp359 = 62370.0*Ftmp153;
double Ftmp360 = -Ftmp24*Ftmp359;
double Ftmp361 = Ftmp358 + Ftmp360;
double Ftmp362 = 17010.0*Ftmp53;
double Ftmp363 = -Ftmp118*Ftmp290 + Ftmp24*Ftmp362;
double Ftmp364 = Ftmp357 + Ftmp361 + Ftmp363;
double Ftmp365 = Ftmp13*Ftmp182;
double Ftmp366 = -Ftmp27*Ftmp359;
double Ftmp367 = Ftmp365 + Ftmp366;
double Ftmp368 = -Ftmp122*Ftmp290 + Ftmp27*Ftmp362;
double Ftmp369 = Ftmp357 + Ftmp367 + Ftmp368;
double Ftmp370 = Ftmp149*Ftmp24;
double Ftmp371 = -Ftmp154*Ftmp24 + Ftmp212;
double Ftmp372 = -Ftmp129;
double Ftmp373 = Ftmp152*Ftmp24;
double Ftmp374 = Ftmp372 + Ftmp373;
double Ftmp375 = Ftmp131 + Ftmp370 + Ftmp371 + Ftmp374;
double Ftmp376 = -Ftmp154*Ftmp27;
double Ftmp377 = Ftmp152*Ftmp27;
double Ftmp378 = Ftmp372 + Ftmp377;
double Ftmp379 = Ftmp149*Ftmp27 + Ftmp212;
double Ftmp380 = Ftmp131 + Ftmp376 + Ftmp378 + Ftmp379;
double Ftmp381 = Ftmp188*Ftmp27;
double Ftmp382 = -Ftmp24*Ftmp381;
double Ftmp383 = Ftmp382 + Ftmp52;
double Ftmp384 = -Ftmp145;
double Ftmp385 = Ftmp182*Ftmp24;
double Ftmp386 = Ftmp384 + Ftmp385;
double Ftmp387 = Ftmp159 + Ftmp259 + Ftmp383 + Ftmp386;
double Ftmp388 = -Ftmp134;
double Ftmp389 = Ftmp165*Ftmp27;
double Ftmp390 = Ftmp388 + Ftmp389;
double Ftmp391 = Ftmp142 + Ftmp271 + Ftmp383 + Ftmp390;
double Ftmp392 = Ftmp13*Ftmp198;
double Ftmp393 = -405405.0*Ftmp242;
double Ftmp394 = 311850.0*Ftmp99;
double Ftmp395 = Ftmp24*Ftmp394;
double Ftmp396 = Ftmp13*Ftmp199;
double Ftmp397 = Ftmp395 + Ftmp396;
double Ftmp398 = -Ftmp215;
double Ftmp399 = 1351350.0*Ftmp196;
double Ftmp400 = -Ftmp24*Ftmp399;
double Ftmp401 = Ftmp398 + Ftmp400;
double Ftmp402 = Ftmp10*(Ftmp392 + Ftmp393 + Ftmp397 + Ftmp401);
double Ftmp403 = Ftmp122*Ftmp194;
double Ftmp404 = Ftmp13*Ftmp403;
double Ftmp405 = 810810.0*Ftmp196;
double Ftmp406 = -Ftmp27*Ftmp405;
double Ftmp407 = Ftmp404 + Ftmp406;
double Ftmp408 = -Ftmp289;
double Ftmp409 = Ftmp291 + Ftmp408;
double Ftmp410 = -405405.0*Ftmp170;
double Ftmp411 = 187110.0*Ftmp172 + Ftmp410;
double Ftmp412 = Ftmp10*(Ftmp407 + Ftmp409 + Ftmp411);
double Ftmp413 = Ftmp195*Ftmp24;
double Ftmp414 = Ftmp198*Ftmp24;
double Ftmp415 = 311850.0*Ftmp153;
double Ftmp416 = -405405.0*Ftmp238;
double Ftmp417 = Ftmp415 + Ftmp416;
double Ftmp418 = Ftmp10*(Ftmp401 + Ftmp413 + Ftmp414 + Ftmp417);
double Ftmp419 = -Ftmp27*Ftmp399;
double Ftmp420 = -Ftmp152;
double Ftmp421 = Ftmp195*Ftmp27;
double Ftmp422 = Ftmp420 + Ftmp421;
double Ftmp423 = Ftmp198*Ftmp27;
double Ftmp424 = Ftmp245 + Ftmp423;
double Ftmp425 = Ftmp10*(Ftmp154 + Ftmp419 + Ftmp422 + Ftmp424);
double Ftmp426 = -810810.0*Ftmp200*Ftmp27;
double Ftmp427 = Ftmp408 + Ftmp426;
double Ftmp428 = Ftmp24*Ftmp403;
double Ftmp429 = Ftmp294 + Ftmp428;
double Ftmp430 = Ftmp411 + Ftmp427 + Ftmp429;
double Ftmp431 = -Ftmp201*Ftmp27;
double Ftmp432 = -Ftmp165;
double Ftmp433 = Ftmp199*Ftmp27;
double Ftmp434 = Ftmp432 + Ftmp433;
double Ftmp435 = Ftmp167 + Ftmp424 + Ftmp431 + Ftmp434;
double Ftmp436 = 187110.0*Ftmp166 + Ftmp393;
double Ftmp437 = -Ftmp24*Ftmp405;
double Ftmp438 = Ftmp396 + Ftmp437;
double Ftmp439 = Ftmp63*(Ftmp409 + Ftmp436 + Ftmp438);
double Ftmp440 = Ftmp398 + Ftmp419;
double Ftmp441 = Ftmp392 + Ftmp404;
double Ftmp442 = Ftmp27*Ftmp394;
double Ftmp443 = Ftmp410 + Ftmp442;
double Ftmp444 = Ftmp63*(Ftmp440 + Ftmp441 + Ftmp443);
double Ftmp445 = Ftmp413 + Ftmp420;
double Ftmp446 = Ftmp245 + Ftmp414;
double Ftmp447 = Ftmp63*(Ftmp154 + Ftmp400 + Ftmp445 + Ftmp446);
double Ftmp448 = Ftmp63*(Ftmp417 + Ftmp421 + Ftmp423 + Ftmp440);
double Ftmp449 = -Ftmp182;
double Ftmp450 = Ftmp428 + Ftmp449;
double Ftmp451 = Ftmp183 + Ftmp431 + Ftmp446 + Ftmp450;
double Ftmp452 = Ftmp301 + Ftmp433;
double Ftmp453 = Ftmp427 + Ftmp436 + Ftmp452;
double Ftmp454 = -Ftmp118*Ftmp338;
double Ftmp455 = Ftmp245 + Ftmp291;
double Ftmp456 = -Ftmp122*Ftmp338 + Ftmp442;
double Ftmp457 = Ftmp192*Ftmp24;
double Ftmp458 = -Ftmp179;
double Ftmp459 = -Ftmp197*Ftmp24 + Ftmp458;
double Ftmp460 = -Ftmp197*Ftmp27;
double Ftmp461 = Ftmp192*Ftmp27 + Ftmp458;
double Ftmp462 = Ftmp106 + Ftmp426;
double Ftmp463 = Ftmp27*Ftmp297;
double Ftmp464 = -Ftmp261 + Ftmp271 + Ftmp463;
double Ftmp465 = Ftmp259 - Ftmp273;
double Ftmp466 = -Ftmp27*Ftmp294 + Ftmp464 + Ftmp465 + Ftmp56;
double Ftmp467 = Ftmp311*Ftmp341;
double Ftmp468 = -Ftmp297 + Ftmp467;
double Ftmp469 = -Ftmp27*Ftmp352 + Ftmp409;
double Ftmp470 = -Ftmp27*Ftmp343 + Ftmp294;
double Ftmp471 = Ftmp10*(Ftmp27*Ftmp349 + Ftmp468 + Ftmp469 + Ftmp470);
double Ftmp472 = -Ftmp303;
double Ftmp473 = -Ftmp24*Ftmp343 + Ftmp301 + Ftmp467;
double Ftmp474 = Ftmp63*(Ftmp24*Ftmp349 + Ftmp469 + Ftmp472 + Ftmp473);
double Ftmp475 = Ftmp329 + Ftmp468;
double Ftmp476 = Ftmp323 + Ftmp472;
double Ftmp477 = x*M[6];
double Ftmp478 = Ftmp17 + Ftmp23;
double Ftmp479 = Ftmp21 + Ftmp25;
double Ftmp480 = Ftmp478*M[4];
double Ftmp481 = 1.0*x;
double Ftmp482 = 3.0*x;
double Ftmp483 = Ftmp20 + Ftmp45;
double Ftmp484 = Ftmp483*M[24];
double Ftmp485 = Ftmp38*M[31];
double Ftmp486 = Ftmp44*x;
double Ftmp487 = 3.0*Ftmp63;
double Ftmp488 = Ftmp483*M[12];
double Ftmp489 = Ftmp56*M[22];
double Ftmp490 = Ftmp483*M[11];
double Ftmp491 = 1.0*Ftmp24;
double Ftmp492 = Ftmp24*x;
double Ftmp493 = Ftmp56*M[21];
double Ftmp494 = Ftmp24*z;
double Ftmp495 = (Ftmp55 + Ftmp72)*M[24];
double Ftmp496 = Ftmp68 + Ftmp83;
double Ftmp497 = Ftmp36*Ftmp63;
double Ftmp498 = Ftmp24*Ftmp481;
double Ftmp499 = Ftmp24*Ftmp482;
double Ftmp500 = Ftmp24*Ftmp44;
double Ftmp501 = Ftmp24*Ftmp63;
double Ftmp502 = (Ftmp101 + Ftmp106)*M[39];
double Ftmp503 = Ftmp107 + Ftmp98;
double Ftmp504 = Ftmp44*Ftmp492;
double Ftmp505 = Ftmp24*Ftmp487;
double Ftmp506 = Ftmp115 - Ftmp120*Ftmp13 + Ftmp32;
double Ftmp507 = Ftmp113 - Ftmp116*Ftmp24 + Ftmp119;
double Ftmp508 = Ftmp506*M[20];
double Ftmp509 = Ftmp13*Ftmp141;
double Ftmp510 = Ftmp129 - Ftmp509 + Ftmp51;
double Ftmp511 = Ftmp510*M[60];
double Ftmp512 = Ftmp137*M[78];
double Ftmp513 = Ftmp510*M[37];
double Ftmp514 = Ftmp155*M[58];
double Ftmp515 = Ftmp510*M[36];
double Ftmp516 = Ftmp155*M[57];
double Ftmp517 = (Ftmp105 + Ftmp152 - Ftmp359)*M[60];
double Ftmp518 = 145530.0*Ftmp166;
double Ftmp519 = Ftmp165 + Ftmp179 - Ftmp518;
double Ftmp520 = (Ftmp195 + Ftmp198 - Ftmp399)*M[88];
double Ftmp521 = 1891890.0*Ftmp200;
double Ftmp522 = -Ftmp114*Ftmp198 + Ftmp13*Ftmp215 + Ftmp209 + Ftmp212;
double Ftmp523 = 218295.0*Ftmp99;
double Ftmp524 = -Ftmp118*Ftmp523 + Ftmp207 + Ftmp210*Ftmp24 + Ftmp214;
double Ftmp525 = Ftmp522*M[56];
double Ftmp526 = Ftmp223 + Ftmp232;
double Ftmp527 = -Ftmp13*Ftmp51 + Ftmp32;
double Ftmp528 = Ftmp526 + Ftmp527;
double Ftmp529 = Ftmp224 + Ftmp228 + Ftmp231;
double Ftmp530 = Ftmp227 + Ftmp234 + Ftmp32;
double Ftmp531 = 467775.0*Ftmp153 + Ftmp237 - 2027025.0*Ftmp238 + Ftmp245;
double Ftmp532 = Ftmp531*M[86];
double Ftmp533 = Ftmp531*M[85];
double Ftmp534 = Ftmp265 + Ftmp51;
double Ftmp535 = Ftmp261 + Ftmp277;
double Ftmp536 = Ftmp534 + Ftmp535;
double Ftmp537 = Ftmp536*M[67];
double Ftmp538 = Ftmp273 + Ftmp286;
double Ftmp539 = Ftmp534 + Ftmp538;
double Ftmp540 = Ftmp539*M[69];
double Ftmp541 = Ftmp260 + Ftmp263 + Ftmp272 + Ftmp278;
double Ftmp542 = Ftmp541*M[80];
double Ftmp543 = Ftmp536*Ftmp6;
double Ftmp544 = Ftmp539*Ftmp6;
double Ftmp545 = Ftmp541*Ftmp6;
double Ftmp546 = -Ftmp13*Ftmp97 + Ftmp82;
double Ftmp547 = Ftmp13*Ftmp203;
double Ftmp548 = -Ftmp547;
double Ftmp549 = Ftmp149 + Ftmp548;
double Ftmp550 = Ftmp105 + Ftmp327;
double Ftmp551 = Ftmp332 + Ftmp97;
double Ftmp552 = Ftmp314*Ftmp497;
double Ftmp553 = Ftmp342 + Ftmp353;
double Ftmp554 = -Ftmp13*Ftmp338 + Ftmp198;
double Ftmp555 = Ftmp13*Ftmp149;
double Ftmp556 = Ftmp136 + Ftmp358 + Ftmp371 + Ftmp388 + Ftmp555;
double Ftmp557 = Ftmp264 + Ftmp52;
double Ftmp558 = Ftmp159 + Ftmp367 + Ftmp384 + Ftmp557;
double Ftmp559 = Ftmp259 + Ftmp356;
double Ftmp560 = -31185.0*Ftmp128 + Ftmp13*Ftmp362;
double Ftmp561 = Ftmp360 + Ftmp373 + Ftmp559 + Ftmp560;
double Ftmp562 = Ftmp509 + Ftmp52;
double Ftmp563 = Ftmp271 + Ftmp366;
double Ftmp564 = Ftmp378 + Ftmp562 + Ftmp563;
double Ftmp565 = Ftmp368 + Ftmp382 + Ftmp385 + Ftmp559;
double Ftmp566 = -Ftmp167*Ftmp27;
double Ftmp567 = Ftmp136 + Ftmp379 + Ftmp390 + Ftmp566;
double Ftmp568 = Ftmp396 + Ftmp432;
double Ftmp569 = Ftmp6*(Ftmp167 + Ftmp245 + Ftmp392 + Ftmp400 + Ftmp568);
double Ftmp570 = Ftmp6*(Ftmp183 + Ftmp245 + Ftmp419 + Ftmp441 + Ftmp449);
double Ftmp571 = Ftmp294 + Ftmp437;
double Ftmp572 = Ftmp413 + Ftmp571;
double Ftmp573 = 187110.0*Ftmp153 + Ftmp408 + Ftmp416;
double Ftmp574 = Ftmp6*(Ftmp572 + Ftmp573);
double Ftmp575 = Ftmp301 + Ftmp406;
double Ftmp576 = Ftmp421 + Ftmp575;
double Ftmp577 = Ftmp6*(Ftmp573 + Ftmp576);
double Ftmp578 = Ftmp398 + Ftmp431;
double Ftmp579 = Ftmp6*(Ftmp414 + Ftmp428 + Ftmp443 + Ftmp578);
double Ftmp580 = Ftmp6*(Ftmp393 + Ftmp395 + Ftmp423 + Ftmp433 + Ftmp578);
double Ftmp581 = Ftmp13*Ftmp192;
double Ftmp582 = Ftmp106 + Ftmp291;
double Ftmp583 = -675675.0*Ftmp238 + Ftmp245 + Ftmp415;
double Ftmp584 = Ftmp106 + Ftmp359;
double Ftmp585 = Ftmp245 + Ftmp426;
double Ftmp586 = -Ftmp27*Ftmp521;
double Ftmp587 = Ftmp264 - Ftmp278;
double Ftmp588 = -Ftmp27*Ftmp291 + Ftmp464 + Ftmp587 + Ftmp69;
double Ftmp589 = -Ftmp312;
double Ftmp590 = Ftmp6*(Ftmp13*Ftmp349 + Ftmp408 + Ftmp470 + Ftmp473 + Ftmp589);
double Ftmp591 = Ftmp547 + Ftmp589;
double Ftmp592 = y*M[8];
double Ftmp593 = Ftmp21 + Ftmp28;
double Ftmp594 = Ftmp481*M[29];
double Ftmp595 = Ftmp40*x;
double Ftmp596 = Ftmp27*x;
double Ftmp597 = Ftmp27*y;
double Ftmp598 = Ftmp40*Ftmp63;
double Ftmp599 = Ftmp27*Ftmp481;
double Ftmp600 = Ftmp27*(Ftmp83 + Ftmp87);
double Ftmp601 = Ftmp10*Ftmp27;
double Ftmp602 = Ftmp36*Ftmp596;
double Ftmp603 = Ftmp40*Ftmp596;
double Ftmp604 = Ftmp113 - Ftmp116*Ftmp27 + Ftmp123;
double Ftmp605 = Ftmp481*M[76];
double Ftmp606 = 145530.0*Ftmp172;
double Ftmp607 = Ftmp27*(Ftmp179 + Ftmp182 - Ftmp606);
double Ftmp608 = -Ftmp122*Ftmp523 + Ftmp207 + Ftmp210*Ftmp27 + Ftmp218;
double Ftmp609 = Ftmp224 + Ftmp526 + Ftmp8;
double Ftmp610 = Ftmp228 + Ftmp230 + Ftmp527;
double Ftmp611 = Ftmp222 + Ftmp230 + Ftmp233 + Ftmp32;
double Ftmp612 = Ftmp335 + Ftmp97;
double Ftmp613 = Ftmp142 + Ftmp361 + Ftmp388 + Ftmp557;
double Ftmp614 = Ftmp146 + Ftmp212;
double Ftmp615 = Ftmp365 + Ftmp376 + Ftmp384 + Ftmp555 + Ftmp614;
double Ftmp616 = Ftmp259 + Ftmp360 + Ftmp374 + Ftmp562;
double Ftmp617 = Ftmp356 + Ftmp377 + Ftmp560 + Ftmp563;
double Ftmp618 = Ftmp370 + Ftmp386 + Ftmp566 + Ftmp614;
double Ftmp619 = Ftmp271 + Ftmp356 + Ftmp363 + Ftmp382 + Ftmp389;
double Ftmp620 = Ftmp458 + Ftmp606;
double Ftmp621 = -Ftmp24*Ftmp291 + Ftmp463 + Ftmp465 + Ftmp587 + Ftmp88;
#pragma omp atomic
F[0] += Ftmp0*M[1] - Ftmp10*Ftmp11 - Ftmp10*Ftmp158 - Ftmp10*Ftmp244 - Ftmp10*Ftmp280*M[52] - Ftmp10*Ftmp430*M[116] - Ftmp10*Ftmp435*M[114] - Ftmp10*Ftmp59 - Ftmp102*Ftmp103*M[39] - Ftmp103*(Ftmp192 + Ftmp195 - Ftmp197)*M[88] - Ftmp103*(Ftmp340 + Ftmp342 + Ftmp345)*M[95] - Ftmp103*(Ftmp345 + Ftmp347 + Ftmp348)*M[97] - Ftmp104*Ftmp108 - Ftmp104*Ftmp202 - Ftmp104*(Ftmp351 + Ftmp355)*M[108] - Ftmp109*Ftmp112 - Ftmp109*Ftmp206 + Ftmp117*x*M[20] + Ftmp117*M[35] - Ftmp12*x + Ftmp121*M[45] + Ftmp124*M[49] + Ftmp125*x + Ftmp126*x - Ftmp13*Ftmp14 - Ftmp13*(Ftmp20 + Ftmp78)*M[10] - Ftmp13*(Ftmp129 - 13230.0*Ftmp130 + Ftmp176)*M[35] - Ftmp13*(Ftmp261 + Ftmp318 + Ftmp319)*M[38] - Ftmp13*(Ftmp273 + Ftmp319 + Ftmp320)*M[40] - Ftmp13*(Ftmp407 + Ftmp455 + Ftmp456)*M[98] - Ftmp13*(Ftmp102 - Ftmp27*Ftmp339 + Ftmp475 + Ftmp476)*M[96] - Ftmp13*(1964655.0*Ftmp153 + Ftmp237 - 3648645.0*Ftmp238 + Ftmp256)*M[84] - Ftmp13*(Ftmp180 + Ftmp422 + Ftmp460 + Ftmp461)*M[89] - Ftmp13*(Ftmp180 + Ftmp445 + Ftmp457 + Ftmp459)*M[87] - Ftmp13*(Ftmp397 + Ftmp437 + Ftmp454 + Ftmp455)*M[94] - Ftmp133*y - Ftmp138*Ftmp36 - Ftmp139*Ftmp40 - Ftmp140*z - Ftmp144*Ftmp44 - Ftmp148*M[76] + Ftmp15*Ftmp155*M[60] + Ftmp15*Ftmp175 + Ftmp15*Ftmp314*M[80] + Ftmp15*Ftmp56*M[24] + Ftmp15*Ftmp77 + Ftmp156*Ftmp6 - Ftmp157*M[36] + Ftmp16*M[8] - Ftmp161*Ftmp62 - Ftmp162*M[37] - Ftmp163*Ftmp63 - Ftmp164*Ftmp63 + Ftmp169*Ftmp67 + Ftmp17*y*M[5] + Ftmp17*z*M[6] + Ftmp174*Ftmp71 - Ftmp177*Ftmp79 - Ftmp178*Ftmp79 + Ftmp181*Ftmp85*M[57] + Ftmp181*Ftmp86*M[58] + Ftmp184*Ftmp185*Ftmp63 + Ftmp186*Ftmp90 + Ftmp187*Ftmp92 + Ftmp190*Ftmp94 + Ftmp191*Ftmp94 - Ftmp20*Ftmp6*M[14] + Ftmp211*x*M[56] + Ftmp211*M[84] + Ftmp216*M[105] + Ftmp219*M[111] + Ftmp22*x*M[4] + Ftmp22*M[10] + Ftmp220*x + Ftmp221*x + Ftmp226*x*M[23] + Ftmp226*M[38] + Ftmp229*x*M[25] + Ftmp229*M[40] + Ftmp235*x*M[32] + Ftmp235*M[47] - Ftmp240*M[85] - Ftmp250*Ftmp62 - Ftmp251*M[86] - Ftmp253*Ftmp63 - Ftmp255*Ftmp63 - Ftmp257*Ftmp79 - Ftmp258*Ftmp79 + Ftmp26*M[13] - Ftmp268*y - Ftmp276*y - Ftmp281*Ftmp36 - Ftmp283*z - Ftmp285*z - Ftmp287*Ftmp63*M[53] - Ftmp288*Ftmp44 + Ftmp29*M[15] - Ftmp3*y + Ftmp30*x + Ftmp300*Ftmp6 + Ftmp306*Ftmp6 - Ftmp307*M[41] - Ftmp308*M[43] - Ftmp309*M[42] + Ftmp31*x - Ftmp310*M[44] + Ftmp315*Ftmp67 + Ftmp316*M[67] + Ftmp317*M[69] - Ftmp35*y - Ftmp36*Ftmp39 + Ftmp364*x*M[66] + Ftmp364*M[94] + Ftmp369*x*M[70] + Ftmp369*M[98] + Ftmp375*x*M[59] + Ftmp375*M[87] + Ftmp380*x*M[61] + Ftmp380*M[89] + Ftmp387*x*M[81] + Ftmp387*M[109] + Ftmp391*x*M[79] + Ftmp391*M[107] - Ftmp4*M[6] - Ftmp40*Ftmp42 - Ftmp402*M[99] - Ftmp412*M[103] - Ftmp418*M[90] - Ftmp425*M[92] - Ftmp43*z - Ftmp439*M[100] - Ftmp44*Ftmp47 - Ftmp444*M[104] - Ftmp447*M[91] - Ftmp448*M[93] - Ftmp451*Ftmp63*M[117] - Ftmp453*Ftmp63*M[115] + Ftmp466*x*M[68] + Ftmp466*M[96] - Ftmp471*M[101] - Ftmp474*M[102] + Ftmp5*x - Ftmp50*M[29] + Ftmp57*Ftmp6 - Ftmp58*M[11] + Ftmp6*Ftmp9 - Ftmp61*Ftmp62 - Ftmp63*Ftmp65 - Ftmp63*Ftmp66 + Ftmp63*Ftmp88*Ftmp89 - Ftmp64*M[12] + Ftmp67*Ftmp70 + Ftmp71*Ftmp76 - Ftmp79*Ftmp80 - Ftmp79*Ftmp81 - Ftmp79*(Ftmp321 + Ftmp322)*M[47] - Ftmp79*(Ftmp189 + Ftmp301 + Ftmp434 + Ftmp462)*M[107] - Ftmp79*(Ftmp294 + Ftmp381 + Ftmp450 + Ftmp462)*M[109] + Ftmp84*Ftmp85*M[21] + Ftmp84*Ftmp86*M[22] + Ftmp85*(Ftmp325 + Ftmp326)*M[62] + Ftmp85*(Ftmp328 + Ftmp331)*M[64] + Ftmp86*(Ftmp325 + Ftmp328)*M[63] + Ftmp86*(Ftmp326 + Ftmp331)*M[65] + Ftmp90*Ftmp91 + Ftmp90*(Ftmp313 + Ftmp333)*M[73] + Ftmp92*Ftmp93 + Ftmp94*Ftmp95 + Ftmp94*Ftmp96 + Ftmp94*(Ftmp336 + Ftmp337)*M[74];
#pragma omp atomic
F[1] += Ftmp0*M[2] - Ftmp10*Ftmp14 - Ftmp11*Ftmp24 - Ftmp112*Ftmp505 - Ftmp12*y + Ftmp124*M[54] + Ftmp126*y - Ftmp133*x - Ftmp137*Ftmp6*M[51] - Ftmp137*Ftmp62*M[45] - Ftmp138*Ftmp481 - Ftmp139*Ftmp482 - Ftmp148*M[82] + Ftmp15*Ftmp489 + Ftmp15*Ftmp514 + Ftmp156*Ftmp63 - Ftmp157*M[35] - Ftmp161*Ftmp491 - Ftmp164*Ftmp6 + Ftmp168*Ftmp497*M[72] + Ftmp169*Ftmp486 + Ftmp174*Ftmp487 - Ftmp178*Ftmp62 + Ftmp184*Ftmp500*M[82] + Ftmp187*Ftmp499 + Ftmp191*Ftmp497 - Ftmp206*Ftmp505 + Ftmp219*M[118] + Ftmp221*y - Ftmp24*Ftmp490 - Ftmp24*Ftmp515 - Ftmp24*Ftmp533 - Ftmp24*(Ftmp37 + Ftmp78)*M[16] - Ftmp24*(Ftmp535 + Ftmp546)*M[41] - Ftmp24*(Ftmp572 + Ftmp583)*M[90] - Ftmp24*(Ftmp134 + Ftmp176 - 13230.0*Ftmp24*Ftmp53)*M[50] - Ftmp24*(Ftmp269 + Ftmp273 + Ftmp321)*M[43] - Ftmp24*(Ftmp320 + Ftmp322 + Ftmp82)*M[52] - Ftmp24*(Ftmp422 + Ftmp575 + Ftmp584)*M[92] - Ftmp24*(Ftmp429 + Ftmp456 + Ftmp585)*M[116] - Ftmp24*(1964655.0*Ftmp166 + Ftmp241 - 3648645.0*Ftmp242 + Ftmp256)*M[112] - Ftmp24*(Ftmp381 + Ftmp407 + Ftmp449 + Ftmp582)*M[103] - Ftmp24*(Ftmp434 + Ftmp461 + Ftmp518 + Ftmp586)*M[114] - Ftmp24*(Ftmp459 + Ftmp518 + Ftmp568 + Ftmp581)*M[99] - Ftmp24*(-Ftmp13*Ftmp346 + Ftmp475 + Ftmp503 + Ftmp591)*M[101] - Ftmp240*M[84] - Ftmp243*Ftmp6*M[113] - Ftmp243*Ftmp62*M[105] + Ftmp25*x*M[5] + Ftmp25*z*M[8] - Ftmp250*Ftmp491 - Ftmp255*Ftmp6 - Ftmp258*Ftmp62 - Ftmp268*x - Ftmp276*x - Ftmp280*Ftmp62*M[47] - Ftmp281*Ftmp481 + Ftmp29*M[18] - Ftmp3*x + Ftmp300*Ftmp63 + Ftmp306*Ftmp63 - Ftmp307*M[38] - Ftmp308*M[40] + Ftmp31*y + Ftmp315*Ftmp486 + Ftmp316*M[63] + Ftmp317*M[65] - Ftmp35*x - Ftmp37*Ftmp63*M[14] - Ftmp38*Ftmp6*M[17] - Ftmp38*Ftmp62*M[13] - Ftmp39*Ftmp481 - Ftmp4*M[8] - Ftmp402*M[94] - Ftmp412*M[98] - Ftmp418*M[87] - Ftmp42*Ftmp482 - Ftmp425*M[89] - Ftmp430*Ftmp62*M[109] - Ftmp435*Ftmp62*M[107] - Ftmp471*M[96] + Ftmp477*Ftmp6*Ftmp8 + Ftmp478*M[11] + Ftmp479*y*M[7] + Ftmp479*M[16] + Ftmp480*y - Ftmp484*z - Ftmp485*z + Ftmp486*Ftmp70 + Ftmp487*Ftmp76 - Ftmp488*Ftmp6 - Ftmp491*Ftmp61 + Ftmp492*Ftmp493 + Ftmp492*Ftmp516 + Ftmp492*(Ftmp298 + Ftmp549)*M[62] + Ftmp492*(Ftmp304 + Ftmp550)*M[64] + Ftmp494*Ftmp495 + Ftmp494*Ftmp496*M[31] + Ftmp494*Ftmp517 + Ftmp494*Ftmp519*M[78] + Ftmp494*(Ftmp149 + Ftmp330 + Ftmp337)*M[80] + Ftmp494*(Ftmp292 + Ftmp303 + Ftmp336)*M[69] + Ftmp494*(Ftmp297 + Ftmp548 + Ftmp551)*M[67] + Ftmp496*Ftmp498*M[26] + Ftmp497*Ftmp69*M[27] + Ftmp497*Ftmp96 + Ftmp498*Ftmp519*M[71] + Ftmp498*(Ftmp312 + Ftmp330 + Ftmp551)*M[73] + Ftmp499*Ftmp93 + Ftmp5*y - Ftmp50*M[33] + Ftmp500*Ftmp88*M[33] - Ftmp501*Ftmp502 - Ftmp501*Ftmp520 - Ftmp501*(Ftmp553 + Ftmp554)*M[95] - Ftmp501*(Ftmp344 + Ftmp348 + Ftmp351)*M[97] - Ftmp503*Ftmp504*M[46] - Ftmp504*(Ftmp192 + Ftmp199 - Ftmp521)*M[106] - Ftmp504*(Ftmp198 + Ftmp347 + Ftmp355)*M[108] + Ftmp506*M[36] + Ftmp507*y*M[30] + Ftmp507*M[50] + Ftmp508*y - Ftmp511*z - Ftmp512*z - Ftmp513*Ftmp6 + Ftmp522*M[85] + Ftmp524*y*M[77] + Ftmp524*M[112] + Ftmp525*y + Ftmp528*y*M[23] + Ftmp528*M[41] + Ftmp529*y*M[25] + Ftmp529*M[43] + Ftmp530*y*M[32] + Ftmp530*M[52] - Ftmp532*Ftmp6 - Ftmp537*z - Ftmp540*z - Ftmp542*z - Ftmp543*M[42] - Ftmp544*M[44] - Ftmp545*M[53] + Ftmp552*M[74] + Ftmp556*y*M[66] + Ftmp556*M[99] + Ftmp558*y*M[70] + Ftmp558*M[103] + Ftmp561*y*M[59] + Ftmp561*M[90] + Ftmp564*y*M[61] + Ftmp564*M[92] + Ftmp565*y*M[81] + Ftmp565*M[116] + Ftmp567*y*M[79] + Ftmp567*M[114] - Ftmp569*M[100] + Ftmp57*Ftmp63 - Ftmp570*M[104] - Ftmp574*M[91] - Ftmp577*M[93] - Ftmp579*M[117] - Ftmp58*M[10] - Ftmp580*M[115] + Ftmp588*y*M[68] + Ftmp588*M[101] - Ftmp590*M[102] - Ftmp6*Ftmp66 - Ftmp62*Ftmp81 + Ftmp63*Ftmp9;
#pragma omp atomic
F[2] += Ftmp0*M[3] + Ftmp10*Ftmp156 + Ftmp10*Ftmp300 + Ftmp10*Ftmp306 - Ftmp10*Ftmp48*M[14] + Ftmp10*Ftmp57 + Ftmp10*Ftmp9 - Ftmp108*Ftmp602 + Ftmp121*M[51] + Ftmp125*z - Ftmp140*x - Ftmp144*Ftmp481 - Ftmp147*Ftmp185 - Ftmp147*Ftmp605 - Ftmp147*Ftmp67*M[54] - Ftmp148*x*M[49] + Ftmp15*Ftmp493 + Ftmp15*Ftmp516 - Ftmp158*Ftmp6 + Ftmp16*M[5] - Ftmp162*M[35] - Ftmp163*Ftmp27 + Ftmp169*Ftmp62 + Ftmp173*Ftmp598*M[75] + Ftmp174*Ftmp595 + Ftmp175*Ftmp597 - Ftmp177*Ftmp486 + Ftmp185*Ftmp607 + Ftmp186*Ftmp497 + Ftmp190*Ftmp599 - Ftmp2*Ftmp27*M[3] - Ftmp2*Ftmp477 - Ftmp2*Ftmp592 - Ftmp202*Ftmp602 + Ftmp216*M[113] + Ftmp220*z - Ftmp244*Ftmp6 - Ftmp251*M[84] - Ftmp253*Ftmp27 - Ftmp254*Ftmp486*M[111] - Ftmp254*Ftmp67*M[118] - Ftmp257*Ftmp486 + Ftmp26*M[17] - Ftmp27*Ftmp488 - Ftmp27*Ftmp513 - Ftmp27*Ftmp532 - Ftmp27*Ftmp65 - Ftmp27*(Ftmp48 + Ftmp78)*M[19] - Ftmp27*(Ftmp538 + Ftmp546)*M[44] - Ftmp27*(Ftmp576 + Ftmp583)*M[93] - Ftmp27*(Ftmp145 + Ftmp176 - 13230.0*Ftmp73)*M[55] - Ftmp27*(Ftmp19 + Ftmp269 + Ftmp535)*M[42] - Ftmp27*(Ftmp445 + Ftmp571 + Ftmp584)*M[91] - Ftmp27*(-3648645.0*Ftmp170 + 1964655.0*Ftmp172 + Ftmp246 + Ftmp256)*M[119] - Ftmp27*(Ftmp189 + Ftmp432 + Ftmp438 + Ftmp582)*M[100] - Ftmp27*(Ftmp278 + Ftmp286 + Ftmp318 + Ftmp82)*M[53] - Ftmp27*(Ftmp395 + Ftmp452 + Ftmp454 + Ftmp585)*M[115] - Ftmp27*(Ftmp450 + Ftmp457 + Ftmp586 + Ftmp620)*M[117] - Ftmp27*(Ftmp404 + Ftmp449 + Ftmp460 + Ftmp581 + Ftmp620)*M[104] - Ftmp27*(-Ftmp13*Ftmp339 + Ftmp334 + Ftmp467 + Ftmp476 + Ftmp591 + Ftmp98)*M[102] + Ftmp28*Ftmp477 + Ftmp28*Ftmp592 - Ftmp283*x - Ftmp285*x - Ftmp287*Ftmp486*M[47] - Ftmp288*Ftmp481 + Ftmp30*z - Ftmp309*M[38] - Ftmp310*M[40] + Ftmp315*Ftmp62 + Ftmp316*M[62] + Ftmp317*M[64] - Ftmp4*x*M[1] - Ftmp4*y*M[2] - Ftmp43*x - Ftmp439*M[94] - Ftmp444*M[98] - Ftmp447*M[87] - Ftmp448*M[89] - Ftmp451*Ftmp486*M[109] - Ftmp453*Ftmp486*M[107] - Ftmp47*Ftmp481 - Ftmp474*M[96] + Ftmp478*M[12] + Ftmp480*z - Ftmp484*y - Ftmp485*y - Ftmp486*Ftmp80 + Ftmp489*Ftmp596 - Ftmp49*Ftmp594 - Ftmp49*Ftmp67*M[18] - Ftmp49*Ftmp89 - Ftmp490*Ftmp6 + Ftmp495*Ftmp597 + Ftmp497*Ftmp91 + Ftmp5*z - Ftmp50*x*M[15] - Ftmp502*Ftmp601 + Ftmp506*M[37] + Ftmp508*z - Ftmp511*y - Ftmp512*y + Ftmp514*Ftmp596 - Ftmp515*Ftmp6 + Ftmp517*Ftmp597 - Ftmp520*Ftmp601 + Ftmp522*M[86] + Ftmp525*z - Ftmp533*Ftmp6 - Ftmp537*y - Ftmp540*y - Ftmp542*y - Ftmp543*M[41] - Ftmp544*M[43] - Ftmp545*M[52] + Ftmp552*M[73] - Ftmp569*M[99] - Ftmp570*M[103] - Ftmp574*M[90] - Ftmp577*M[92] - Ftmp579*M[116] - Ftmp580*M[114] - Ftmp59*Ftmp6 - Ftmp590*M[101] + Ftmp593*z*M[9] + Ftmp593*M[19] + Ftmp594*Ftmp600 + Ftmp595*Ftmp76 + Ftmp596*(Ftmp298 + Ftmp550)*M[63] + Ftmp596*(Ftmp304 + Ftmp549)*M[65] + Ftmp597*Ftmp77 + Ftmp597*(Ftmp149 + Ftmp313 + Ftmp324)*M[80] + Ftmp597*(Ftmp292 + Ftmp297 + Ftmp333)*M[67] + Ftmp597*(Ftmp303 + Ftmp548 + Ftmp612)*M[69] + Ftmp598*Ftmp75*M[28] + Ftmp599*Ftmp95 + Ftmp599*(Ftmp312 + Ftmp324 + Ftmp612)*M[74] + Ftmp600*Ftmp89 - Ftmp601*(Ftmp344 + Ftmp349 + Ftmp553)*M[95] - Ftmp601*(Ftmp348 + Ftmp350 + Ftmp554)*M[97] - Ftmp602*(Ftmp198 + Ftmp340 + Ftmp350 + Ftmp354)*M[108] - Ftmp603*(Ftmp111 - 1575.0*Ftmp53)*M[48] - Ftmp603*(Ftmp204 - 630630.0*Ftmp205 + 121275.0*Ftmp99)*M[110] + Ftmp604*z*M[34] + Ftmp604*M[55] + Ftmp605*Ftmp607 + Ftmp608*z*M[83] + Ftmp608*M[119] + Ftmp609*z*M[23] + Ftmp609*M[42] + Ftmp610*z*M[25] + Ftmp610*M[44] + Ftmp611*z*M[32] + Ftmp611*M[53] + Ftmp613*z*M[66] + Ftmp613*M[100] + Ftmp615*z*M[70] + Ftmp615*M[104] + Ftmp616*z*M[59] + Ftmp616*M[91] + Ftmp617*z*M[61] + Ftmp617*M[93] + Ftmp618*z*M[81] + Ftmp618*M[117] + Ftmp619*z*M[79] + Ftmp619*M[115] + Ftmp62*Ftmp70 + Ftmp621*z*M[68] + Ftmp621*M[102] - Ftmp64*M[10];

}

template<>
void P2M<0, 3>(double x, double y, double z, double q, double * M, int order) {
switch (order) {
  case 1:
    field_m0_P2M_1(x, y, z, q, M);
    break;
  case 2:
    field_m0_P2M_2(x, y, z, q, M);
    break;
  case 3:
    field_m0_P2M_3(x, y, z, q, M);
    break;
  case 4:
    field_m0_P2M_4(x, y, z, q, M);
    break;
  case 5:
    field_m0_P2M_5(x, y, z, q, M);
    break;
  case 6:
    field_m0_P2M_6(x, y, z, q, M);
    break;
  case 7:
    field_m0_P2M_7(x, y, z, q, M);
    break;
  }
}
template<>
void M2M<0, 3>(double x, double y, double z, double * M, double * Ms, int order) {
switch (order) {
  case 1:
    field_m0_M2M_1(x, y, z, M, Ms);
    break;
  case 2:
    field_m0_M2M_2(x, y, z, M, Ms);
    break;
  case 3:
    field_m0_M2M_3(x, y, z, M, Ms);
    break;
  case 4:
    field_m0_M2M_4(x, y, z, M, Ms);
    break;
  case 5:
    field_m0_M2M_5(x, y, z, M, Ms);
    break;
  case 6:
    field_m0_M2M_6(x, y, z, M, Ms);
    break;
  case 7:
    field_m0_M2M_7(x, y, z, M, Ms);
    break;
  }
}
template<>
void M2L<0, 3>(double x, double y, double z, double * M, double * L, int order) {
switch (order) {
  case 1:
    field_m0_M2L_1(x, y, z, M, L);
    break;
  case 2:
    field_m0_M2L_2(x, y, z, M, L);
    break;
  case 3:
    field_m0_M2L_3(x, y, z, M, L);
    break;
  case 4:
    field_m0_M2L_4(x, y, z, M, L);
    break;
  case 5:
    field_m0_M2L_5(x, y, z, M, L);
    break;
  case 6:
    field_m0_M2L_6(x, y, z, M, L);
    break;
  case 7:
    field_m0_M2L_7(x, y, z, M, L);
    break;
  }
}
template<>
void L2L<0, 3>(double x, double y, double z, double * L, double * Ls, int order) {
switch (order) {
  case 1:
    field_m0_L2L_1(x, y, z, L, Ls);
    break;
  case 2:
    field_m0_L2L_2(x, y, z, L, Ls);
    break;
  case 3:
    field_m0_L2L_3(x, y, z, L, Ls);
    break;
  case 4:
    field_m0_L2L_4(x, y, z, L, Ls);
    break;
  case 5:
    field_m0_L2L_5(x, y, z, L, Ls);
    break;
  case 6:
    field_m0_L2L_6(x, y, z, L, Ls);
    break;
  case 7:
    field_m0_L2L_7(x, y, z, L, Ls);
    break;
  }
}
template<>
void L2P<0, 3>(double x, double y, double z, double * L, double * F, int order) {
switch (order) {
  case 1:
    field_m0_L2P_1(x, y, z, L, F);
    break;
  case 2:
    field_m0_L2P_2(x, y, z, L, F);
    break;
  case 3:
    field_m0_L2P_3(x, y, z, L, F);
    break;
  case 4:
    field_m0_L2P_4(x, y, z, L, F);
    break;
  case 5:
    field_m0_L2P_5(x, y, z, L, F);
    break;
  case 6:
    field_m0_L2P_6(x, y, z, L, F);
    break;
  case 7:
    field_m0_L2P_7(x, y, z, L, F);
    break;
  }
}
template<>
void M2P<0, 3>(double x, double y, double z, double * M, double * F, int order) {
switch (order) {
  case 1:
    field_m0_M2P_1(x, y, z, M, F);
    break;
  case 2:
    field_m0_M2P_2(x, y, z, M, F);
    break;
  case 3:
    field_m0_M2P_3(x, y, z, M, F);
    break;
  case 4:
    field_m0_M2P_4(x, y, z, M, F);
    break;
  case 5:
    field_m0_M2P_5(x, y, z, M, F);
    break;
  case 6:
    field_m0_M2P_6(x, y, z, M, F);
    break;
  case 7:
    field_m0_M2P_7(x, y, z, M, F);
    break;
  }
}

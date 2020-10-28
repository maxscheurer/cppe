#include "field_m1.hh"
#include <cmath>
void field_m1_P2M_2(double x, double y, double z, double q, double* M) {
  double Mtmp0 = q * x;
  double Mtmp1 = q * y;
  double Mtmp2 = (1.0 / 2.0) * q;
  M[0] += -Mtmp0;
  M[1] += -Mtmp1;
  M[2] += -q * z;
  M[3] += Mtmp2 * (x * x);
  M[4] += Mtmp0 * y;
  M[5] += Mtmp0 * z;
  M[6] += Mtmp2 * (y * y);
  M[7] += Mtmp1 * z;
  M[8] += Mtmp2 * (z * z);
}
void field_m1_M2M_2(double x, double y, double z, double* M, double* Ms) {
#pragma omp atomic
  Ms[0] += M[0];
#pragma omp atomic
  Ms[1] += M[1];
#pragma omp atomic
  Ms[2] += M[2];
#pragma omp atomic
  Ms[3] += x * M[0] + M[3];
#pragma omp atomic
  Ms[4] += x * M[1] + y * M[0] + M[4];
#pragma omp atomic
  Ms[5] += x * M[2] + z * M[0] + M[5];
#pragma omp atomic
  Ms[6] += y * M[1] + M[6];
#pragma omp atomic
  Ms[7] += y * M[2] + z * M[1] + M[7];
#pragma omp atomic
  Ms[8] += z * M[2] + M[8];
}

void field_m1_M2L_2(double x, double y, double z, double* M, double* L) {
  double R = sqrt(x * x + y * y + z * z);
  double D[9];
  double Dtmp0 = 1.0 * pow(R, -3.0);
  double Dtmp1 = -Dtmp0;
  double Dtmp2 = 3.0 * pow(R, -5.0);
  double Dtmp3 = Dtmp2 * x;
  D[0]         = -Dtmp0 * x;
  D[1]         = -Dtmp0 * y;
  D[2]         = -Dtmp0 * z;
  D[3]         = Dtmp1 + Dtmp2 * (x * x);
  D[4]         = Dtmp3 * y;
  D[5]         = Dtmp3 * z;
  D[6]         = Dtmp1 + Dtmp2 * (y * y);
  D[7]         = Dtmp2 * y * z;
  D[8]         = -D[3] - D[6];
#pragma omp atomic
  L[0] += D[0] * M[0] + D[1] * M[1] + D[2] * M[2] + D[3] * M[3] + D[4] * M[4] +
          D[5] * M[5] + D[6] * M[6] + D[7] * M[7] + D[8] * M[8];
#pragma omp atomic
  L[1] += D[3] * M[0] + D[4] * M[1] + D[5] * M[2];
#pragma omp atomic
  L[2] += D[4] * M[0] + D[6] * M[1] + D[7] * M[2];
#pragma omp atomic
  L[3] += D[5] * M[0] + D[7] * M[1] + D[8] * M[2];
}

void field_m1_L2L_2(double x, double y, double z, double* L, double* Ls) {
#pragma omp atomic
  Ls[0] += x * L[1] + y * L[2] + z * L[3] + L[0];
#pragma omp atomic
  Ls[1] += L[1];
#pragma omp atomic
  Ls[2] += L[2];
#pragma omp atomic
  Ls[3] += L[3];
}

void field_m1_L2P_2(double x, double y, double z, double* L, double* F) {
#pragma omp atomic
  F[0] += -L[1];
#pragma omp atomic
  F[1] += -L[2];
#pragma omp atomic
  F[2] += -L[3];
}

void field_m1_M2P_2(double x, double y, double z, double* M, double* F) {
  double R      = sqrt(x * x + y * y + z * z);
  double Ftmp0  = 1.0 * pow(R, -3.0);
  double Ftmp1  = pow(R, -5.0);
  double Ftmp2  = 3.0 * Ftmp1;
  double Ftmp3  = Ftmp2 * M[4];
  double Ftmp4  = Ftmp2 * z;
  double Ftmp5  = Ftmp2 * M[1];
  double Ftmp6  = x * y;
  double Ftmp7  = Ftmp4 * M[2];
  double Ftmp8  = (x * x);
  double Ftmp9  = Ftmp2 * M[0];
  double Ftmp10 = y * M[7];
  double Ftmp11 = 15.0 * pow(R, -7.0);
  double Ftmp12 = Ftmp11 * z;
  double Ftmp13 = Ftmp12 * x;
  double Ftmp14 = Ftmp11 * Ftmp8;
  double Ftmp15 = y * M[4];
  double Ftmp16 = -9.0 * Ftmp1;
  double Ftmp17 = -Ftmp2;
  double Ftmp18 = (y * y);
  double Ftmp19 = Ftmp11 * Ftmp18;
  double Ftmp20 = (Ftmp17 + Ftmp19) * M[6];
  double Ftmp21 = (z * z);
  double Ftmp22 = Ftmp11 * Ftmp21;
  double Ftmp23 = (Ftmp17 + Ftmp22) * M[8];
  double Ftmp24 = x * M[5];
  double Ftmp25 = (Ftmp14 + Ftmp17) * M[3];
#pragma omp atomic
  F[0] += Ftmp0 * M[0] + Ftmp10 * Ftmp13 + Ftmp14 * Ftmp15 + Ftmp14 * z * M[5] +
          Ftmp20 * x + Ftmp23 * x - Ftmp3 * y - Ftmp4 * M[5] - Ftmp5 * Ftmp6 - Ftmp7 * x -
          Ftmp8 * Ftmp9 + x * (Ftmp14 + Ftmp16) * M[3];
#pragma omp atomic
  F[1] += Ftmp0 * M[1] + Ftmp12 * Ftmp24 * y - Ftmp18 * Ftmp5 + Ftmp19 * x * M[4] +
          Ftmp19 * z * M[7] + Ftmp23 * y + Ftmp25 * y - Ftmp3 * x - Ftmp4 * M[7] -
          Ftmp6 * Ftmp9 - Ftmp7 * y + y * (Ftmp16 + Ftmp19) * M[6];
#pragma omp atomic
  F[2] += Ftmp0 * M[2] - Ftmp10 * Ftmp2 + Ftmp10 * Ftmp22 + Ftmp13 * Ftmp15 -
          Ftmp2 * Ftmp21 * M[2] - Ftmp2 * Ftmp24 + Ftmp20 * z + Ftmp22 * Ftmp24 +
          Ftmp25 * z - Ftmp4 * x * M[0] - Ftmp4 * y * M[1] + z * (Ftmp16 + Ftmp22) * M[8];
}

template <>
void P2P<1, 3>(double x, double y, double z, double* S, double* F) {
  double R     = sqrt(x * x + y * y + z * z);
  double Ftmp0 = 1.0 * pow(R, -3.0);
  double Ftmp1 = 3.0 * pow(R, -5.0);
  double Ftmp2 = Ftmp1 * S[1];
  double Ftmp3 = x * y;
  double Ftmp4 = Ftmp1 * S[2];
  double Ftmp5 = Ftmp4 * z;
  double Ftmp6 = Ftmp1 * S[0];
#pragma omp atomic
  F[0] += Ftmp0 * S[0] - Ftmp2 * Ftmp3 - Ftmp5 * x - Ftmp6 * (x * x);
#pragma omp atomic
  F[1] += Ftmp0 * S[1] - Ftmp2 * (y * y) - Ftmp3 * Ftmp6 - Ftmp5 * y;
#pragma omp atomic
  F[2] += Ftmp0 * S[2] - Ftmp2 * y * z - Ftmp4 * (z * z) - Ftmp6 * x * z;
}

void field_m1_P2M_3(double x, double y, double z, double q, double* M) {
  double Mtmp0  = q * x;
  double Mtmp1  = q * y;
  double Mtmp2  = q * z;
  double Mtmp3  = (x * x);
  double Mtmp4  = (1.0 / 2.0) * q;
  double Mtmp5  = Mtmp0 * y;
  double Mtmp6  = (y * y);
  double Mtmp7  = (z * z);
  double Mtmp8  = (1.0 / 6.0) * q;
  double Mtmp9  = (1.0 / 2.0) * Mtmp3;
  double Mtmp10 = (1.0 / 2.0) * Mtmp0;
  M[0] += -Mtmp0;
  M[1] += -Mtmp1;
  M[2] += -Mtmp2;
  M[3] += Mtmp3 * Mtmp4;
  M[4] += Mtmp5;
  M[5] += Mtmp0 * z;
  M[6] += Mtmp4 * Mtmp6;
  M[7] += Mtmp1 * z;
  M[8] += Mtmp4 * Mtmp7;
  M[9] += -Mtmp8 * (x * x * x);
  M[10] += -Mtmp1 * Mtmp9;
  M[11] += -Mtmp2 * Mtmp9;
  M[12] += -Mtmp10 * Mtmp6;
  M[13] += -Mtmp5 * z;
  M[14] += -Mtmp10 * Mtmp7;
  M[15] += -Mtmp8 * (y * y * y);
  M[16] += -1.0 / 2.0 * Mtmp2 * Mtmp6;
  M[17] += -1.0 / 2.0 * Mtmp1 * Mtmp7;
  M[18] += -Mtmp8 * (z * z * z);
}
void field_m1_M2M_3(double x, double y, double z, double* M, double* Ms) {
  double Mstmp0  = x * M[0];
  double Mstmp1  = x * M[1];
  double Mstmp2  = y * M[0];
  double Mstmp3  = x * M[2];
  double Mstmp4  = y * M[1];
  double Mstmp5  = y * M[2];
  double Mstmp6  = (1.0 / 2.0) * (x * x);
  double Mstmp7  = (y * y);
  double Mstmp8  = (1.0 / 2.0) * M[0];
  double Mstmp9  = (z * z);
  double Mstmp10 = (1.0 / 2.0) * Mstmp7;
  double Mstmp11 = (1.0 / 2.0) * Mstmp9;
#pragma omp atomic
  Ms[0] += M[0];
#pragma omp atomic
  Ms[1] += M[1];
#pragma omp atomic
  Ms[2] += M[2];
#pragma omp atomic
  Ms[3] += Mstmp0 + M[3];
#pragma omp atomic
  Ms[4] += Mstmp1 + Mstmp2 + M[4];
#pragma omp atomic
  Ms[5] += Mstmp3 + z * M[0] + M[5];
#pragma omp atomic
  Ms[6] += Mstmp4 + M[6];
#pragma omp atomic
  Ms[7] += Mstmp5 + z * M[1] + M[7];
#pragma omp atomic
  Ms[8] += z * M[2] + M[8];
#pragma omp atomic
  Ms[9] += Mstmp6 * M[0] + x * M[3] + M[9];
#pragma omp atomic
  Ms[10] += Mstmp0 * y + Mstmp6 * M[1] + x * M[4] + y * M[3] + M[10];
#pragma omp atomic
  Ms[11] += Mstmp0 * z + Mstmp6 * M[2] + x * M[5] + z * M[3] + M[11];
#pragma omp atomic
  Ms[12] += Mstmp1 * y + Mstmp7 * Mstmp8 + x * M[6] + y * M[4] + M[12];
#pragma omp atomic
  Ms[13] += Mstmp1 * z + Mstmp2 * z + Mstmp3 * y + x * M[7] + y * M[5] + z * M[4] + M[13];
#pragma omp atomic
  Ms[14] += Mstmp3 * z + Mstmp8 * Mstmp9 + x * M[8] + z * M[5] + M[14];
#pragma omp atomic
  Ms[15] += Mstmp10 * M[1] + y * M[6] + M[15];
#pragma omp atomic
  Ms[16] += Mstmp10 * M[2] + Mstmp4 * z + y * M[7] + z * M[6] + M[16];
#pragma omp atomic
  Ms[17] += Mstmp11 * M[1] + Mstmp5 * z + y * M[8] + z * M[7] + M[17];
#pragma omp atomic
  Ms[18] += Mstmp11 * M[2] + z * M[8] + M[18];
}

void field_m1_M2L_3(double x, double y, double z, double* M, double* L) {
  double R = sqrt(x * x + y * y + z * z);
  double D[19];
  double Dtmp0  = 1.0 * pow(R, -3.0);
  double Dtmp1  = -Dtmp0;
  double Dtmp2  = (x * x);
  double Dtmp3  = pow(R, -5.0);
  double Dtmp4  = 3.0 * Dtmp3;
  double Dtmp5  = Dtmp4 * x;
  double Dtmp6  = (y * y);
  double Dtmp7  = y * z;
  double Dtmp8  = -9.0 * Dtmp3;
  double Dtmp9  = 15.0 * pow(R, -7.0);
  double Dtmp10 = Dtmp2 * Dtmp9;
  double Dtmp11 = -Dtmp4;
  double Dtmp12 = Dtmp10 + Dtmp11;
  double Dtmp13 = Dtmp6 * Dtmp9;
  double Dtmp14 = Dtmp11 + Dtmp13;
  D[0]          = -Dtmp0 * x;
  D[1]          = -Dtmp0 * y;
  D[2]          = -Dtmp0 * z;
  D[3]          = Dtmp1 + Dtmp2 * Dtmp4;
  D[4]          = Dtmp5 * y;
  D[5]          = Dtmp5 * z;
  D[6]          = Dtmp1 + Dtmp4 * Dtmp6;
  D[7]          = Dtmp4 * Dtmp7;
  D[8]          = -D[3] - D[6];
  D[9]          = -x * (Dtmp10 + Dtmp8);
  D[10]         = -Dtmp12 * y;
  D[11]         = -Dtmp12 * z;
  D[12]         = -1.0 * Dtmp14 * x;
  D[13]         = -Dtmp7 * Dtmp9 * x;
  D[14]         = -D[9] - D[12];
  D[15]         = -y * (Dtmp13 + Dtmp8);
  D[16]         = -Dtmp14 * z;
  D[17]         = -D[10] - D[15];
  D[18]         = -D[11] - D[16];
#pragma omp atomic
  L[0] += D[0] * M[0] + D[1] * M[1] + D[2] * M[2] + D[3] * M[3] + D[4] * M[4] +
          D[5] * M[5] + D[6] * M[6] + D[7] * M[7] + D[8] * M[8] + D[9] * M[9] +
          D[10] * M[10] + D[11] * M[11] + D[12] * M[12] + D[13] * M[13] + D[14] * M[14] +
          D[15] * M[15] + D[16] * M[16] + D[17] * M[17] + D[18] * M[18];
#pragma omp atomic
  L[1] += D[3] * M[0] + D[4] * M[1] + D[5] * M[2] + D[9] * M[3] + D[10] * M[4] +
          D[11] * M[5] + D[12] * M[6] + D[13] * M[7] + D[14] * M[8];
#pragma omp atomic
  L[2] += D[4] * M[0] + D[6] * M[1] + D[7] * M[2] + D[10] * M[3] + D[12] * M[4] +
          D[13] * M[5] + D[15] * M[6] + D[16] * M[7] + D[17] * M[8];
#pragma omp atomic
  L[3] += D[5] * M[0] + D[7] * M[1] + D[8] * M[2] + D[11] * M[3] + D[13] * M[4] +
          D[14] * M[5] + D[16] * M[6] + D[17] * M[7] + D[18] * M[8];
#pragma omp atomic
  L[4] += D[9] * M[0] + D[10] * M[1] + D[11] * M[2];
#pragma omp atomic
  L[5] += D[10] * M[0] + D[12] * M[1] + D[13] * M[2];
#pragma omp atomic
  L[6] += D[11] * M[0] + D[13] * M[1] + D[14] * M[2];
#pragma omp atomic
  L[7] += D[12] * M[0] + D[15] * M[1] + D[16] * M[2];
#pragma omp atomic
  L[8] += D[13] * M[0] + D[16] * M[1] + D[17] * M[2];
#pragma omp atomic
  L[9] += D[14] * M[0] + D[17] * M[1] + D[18] * M[2];
}

void field_m1_L2L_3(double x, double y, double z, double* L, double* Ls) {
  double Lstmp0 = y * L[5];
  double Lstmp1 = z * L[6];
  double Lstmp2 = z * L[8];
#pragma omp atomic
  Ls[0] += Lstmp0 * x + Lstmp1 * x + Lstmp2 * y + (1.0 / 2.0) * (x * x) * L[4] +
           x * L[1] + (1.0 / 2.0) * (y * y) * L[7] + y * L[2] +
           (1.0 / 2.0) * (z * z) * L[9] + z * L[3] + L[0];
#pragma omp atomic
  Ls[1] += Lstmp0 + Lstmp1 + x * L[4] + L[1];
#pragma omp atomic
  Ls[2] += Lstmp2 + x * L[5] + y * L[7] + L[2];
#pragma omp atomic
  Ls[3] += x * L[6] + y * L[8] + z * L[9] + L[3];
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

void field_m1_L2P_3(double x, double y, double z, double* L, double* F) {
#pragma omp atomic
  F[0] += -x * L[4] - y * L[5] - z * L[6] - L[1];
#pragma omp atomic
  F[1] += -x * L[5] - y * L[7] - z * L[8] - L[2];
#pragma omp atomic
  F[2] += -x * L[6] - y * L[8] - z * L[9] - L[3];
}

void field_m1_M2P_3(double x, double y, double z, double* M, double* F) {
  double R      = sqrt(x * x + y * y + z * z);
  double Ftmp0  = 1.0 * pow(R, -3.0);
  double Ftmp1  = pow(R, -5.0);
  double Ftmp2  = 3.0 * Ftmp1;
  double Ftmp3  = Ftmp2 * M[4];
  double Ftmp4  = Ftmp2 * z;
  double Ftmp5  = y * z;
  double Ftmp6  = pow(R, -7.0);
  double Ftmp7  = 15.0 * Ftmp6;
  double Ftmp8  = Ftmp7 * M[13];
  double Ftmp9  = x * y;
  double Ftmp10 = Ftmp2 * M[1];
  double Ftmp11 = Ftmp4 * M[2];
  double Ftmp12 = (x * x);
  double Ftmp13 = Ftmp2 * M[0];
  double Ftmp14 = z * M[7];
  double Ftmp15 = Ftmp7 * Ftmp9;
  double Ftmp16 = Ftmp12 * Ftmp7;
  double Ftmp17 = z * M[5];
  double Ftmp18 = 105.0 * pow(R, -9.0);
  double Ftmp19 = Ftmp12 * Ftmp18;
  double Ftmp20 = -9.0 * Ftmp1;
  double Ftmp21 = Ftmp16 + Ftmp20;
  double Ftmp22 = -Ftmp2;
  double Ftmp23 = (y * y);
  double Ftmp24 = Ftmp23 * Ftmp7;
  double Ftmp25 = Ftmp22 + Ftmp24;
  double Ftmp26 = (z * z);
  double Ftmp27 = Ftmp26 * Ftmp7;
  double Ftmp28 = Ftmp22 + Ftmp27;
  double Ftmp29 = Ftmp25 * M[6];
  double Ftmp30 = Ftmp28 * M[8];
  double Ftmp31 = -45.0 * Ftmp6;
  double Ftmp32 = Ftmp19 + Ftmp31;
  double Ftmp33 = Ftmp32 * Ftmp9;
  double Ftmp34 = Ftmp18 * Ftmp23;
  double Ftmp35 = Ftmp31 + Ftmp34;
  double Ftmp36 = Ftmp35 * M[15];
  double Ftmp37 = -Ftmp7;
  double Ftmp38 = Ftmp18 * Ftmp26;
  double Ftmp39 = Ftmp37 + Ftmp38;
  double Ftmp40 = 1.0 * M[17];
  double Ftmp41 = Ftmp39 * Ftmp40;
  double Ftmp42 = x * z;
  double Ftmp43 = Ftmp32 * Ftmp42;
  double Ftmp44 = Ftmp34 + Ftmp37;
  double Ftmp45 = Ftmp44 * M[16];
  double Ftmp46 = Ftmp31 + Ftmp38;
  double Ftmp47 = Ftmp46 * M[18];
  double Ftmp48 = -75.0 * Ftmp6;
  double Ftmp49 = 1.0 * Ftmp12;
  double Ftmp50 = Ftmp44 * M[12];
  double Ftmp51 = Ftmp39 * M[14];
  double Ftmp52 = Ftmp16 + Ftmp22;
  double Ftmp53 = Ftmp20 + Ftmp24;
  double Ftmp54 = Ftmp52 * M[3];
  double Ftmp55 = 1.0 * Ftmp9;
  double Ftmp56 = Ftmp19 + Ftmp37;
  double Ftmp57 = Ftmp56 * M[11];
  double Ftmp58 = Ftmp56 * M[10];
  double Ftmp59 = x * M[5];
  double Ftmp60 = y * M[7];
  double Ftmp61 = Ftmp20 + Ftmp27;
  double Ftmp62 = 1.0 * Ftmp42;
#pragma omp atomic
  F[0] += Ftmp0 * M[0] - Ftmp10 * Ftmp9 - Ftmp11 * x - Ftmp12 * Ftmp13 -
          Ftmp12 * (Ftmp19 + Ftmp48) * M[9] + Ftmp14 * Ftmp15 + Ftmp16 * Ftmp17 +
          Ftmp16 * y * M[4] - Ftmp19 * Ftmp5 * M[13] + Ftmp21 * x * M[3] + Ftmp21 * M[9] +
          Ftmp25 * M[12] + Ftmp28 * M[14] + Ftmp29 * x - Ftmp3 * y + Ftmp30 * x -
          Ftmp33 * M[10] - Ftmp36 * Ftmp9 - Ftmp4 * M[5] - Ftmp41 * Ftmp9 -
          Ftmp42 * Ftmp45 - Ftmp42 * Ftmp47 - Ftmp43 * M[11] - Ftmp49 * Ftmp50 -
          Ftmp49 * Ftmp51 + Ftmp5 * Ftmp8;
#pragma omp atomic
  F[1] += Ftmp0 * M[1] - Ftmp10 * Ftmp23 - Ftmp11 * y - Ftmp13 * Ftmp9 + Ftmp14 * Ftmp24 +
          Ftmp15 * Ftmp17 - Ftmp23 * Ftmp41 - Ftmp23 * Ftmp58 -
          Ftmp23 * (Ftmp34 + Ftmp48) * M[15] + Ftmp24 * x * M[4] + Ftmp28 * M[17] -
          Ftmp3 * x + Ftmp30 * y - Ftmp33 * M[9] - Ftmp34 * Ftmp42 * M[13] -
          Ftmp35 * Ftmp5 * M[16] - Ftmp35 * Ftmp55 * M[12] - Ftmp4 * M[7] +
          Ftmp42 * Ftmp8 - Ftmp47 * Ftmp5 - Ftmp5 * Ftmp57 - Ftmp51 * Ftmp55 +
          Ftmp52 * M[10] + Ftmp53 * y * M[6] + Ftmp53 * M[15] + Ftmp54 * y;
#pragma omp atomic
  F[2] += Ftmp0 * M[2] + Ftmp15 * z * M[4] - Ftmp2 * Ftmp26 * M[2] - Ftmp2 * Ftmp59 -
          Ftmp2 * Ftmp60 + Ftmp25 * M[16] - Ftmp26 * Ftmp45 - Ftmp26 * Ftmp57 -
          Ftmp26 * (Ftmp38 + Ftmp48) * M[18] + Ftmp27 * Ftmp59 + Ftmp27 * Ftmp60 +
          Ftmp29 * z - Ftmp36 * Ftmp5 - Ftmp38 * Ftmp9 * M[13] - Ftmp4 * x * M[0] -
          Ftmp4 * y * M[1] - Ftmp40 * Ftmp46 * Ftmp5 - Ftmp43 * M[9] -
          Ftmp46 * Ftmp62 * M[14] - Ftmp5 * Ftmp58 - Ftmp50 * Ftmp62 + Ftmp52 * M[11] +
          Ftmp54 * z + Ftmp61 * z * M[8] + Ftmp61 * M[18] + Ftmp8 * Ftmp9;
}

void field_m1_P2M_4(double x, double y, double z, double q, double* M) {
  double Mtmp0  = q * x;
  double Mtmp1  = q * y;
  double Mtmp2  = q * z;
  double Mtmp3  = (x * x);
  double Mtmp4  = (1.0 / 2.0) * q;
  double Mtmp5  = Mtmp0 * y;
  double Mtmp6  = Mtmp0 * z;
  double Mtmp7  = (y * y);
  double Mtmp8  = Mtmp1 * z;
  double Mtmp9  = (z * z);
  double Mtmp10 = (x * x * x);
  double Mtmp11 = (1.0 / 6.0) * q;
  double Mtmp12 = (1.0 / 2.0) * Mtmp3;
  double Mtmp13 = (1.0 / 2.0) * Mtmp0;
  double Mtmp14 = (y * y * y);
  double Mtmp15 = (1.0 / 2.0) * Mtmp7;
  double Mtmp16 = (1.0 / 2.0) * Mtmp9;
  double Mtmp17 = (z * z * z);
  double Mtmp18 = (1.0 / 24.0) * q;
  double Mtmp19 = (1.0 / 6.0) * Mtmp10;
  double Mtmp20 = (1.0 / 4.0) * Mtmp3 * q;
  double Mtmp21 = (1.0 / 6.0) * Mtmp0;
  M[0] += -Mtmp0;
  M[1] += -Mtmp1;
  M[2] += -Mtmp2;
  M[3] += Mtmp3 * Mtmp4;
  M[4] += Mtmp5;
  M[5] += Mtmp6;
  M[6] += Mtmp4 * Mtmp7;
  M[7] += Mtmp8;
  M[8] += Mtmp4 * Mtmp9;
  M[9] += -Mtmp10 * Mtmp11;
  M[10] += -Mtmp1 * Mtmp12;
  M[11] += -Mtmp12 * Mtmp2;
  M[12] += -Mtmp13 * Mtmp7;
  M[13] += -Mtmp5 * z;
  M[14] += -Mtmp13 * Mtmp9;
  M[15] += -Mtmp11 * Mtmp14;
  M[16] += -Mtmp15 * Mtmp2;
  M[17] += -Mtmp1 * Mtmp16;
  M[18] += -Mtmp11 * Mtmp17;
  M[19] += Mtmp18 * (x * x * x * x);
  M[20] += Mtmp1 * Mtmp19;
  M[21] += Mtmp19 * Mtmp2;
  M[22] += Mtmp20 * Mtmp7;
  M[23] += Mtmp12 * Mtmp8;
  M[24] += Mtmp20 * Mtmp9;
  M[25] += Mtmp14 * Mtmp21;
  M[26] += Mtmp15 * Mtmp6;
  M[27] += Mtmp16 * Mtmp5;
  M[28] += Mtmp17 * Mtmp21;
  M[29] += Mtmp18 * (y * y * y * y);
  M[30] += (1.0 / 6.0) * Mtmp14 * Mtmp2;
  M[31] += (1.0 / 4.0) * Mtmp7 * Mtmp9 * q;
  M[32] += (1.0 / 6.0) * Mtmp1 * Mtmp17;
  M[33] += Mtmp18 * (z * z * z * z);
}
void field_m1_M2M_4(double x, double y, double z, double* M, double* Ms) {
  double Mstmp0  = x * M[0];
  double Mstmp1  = x * M[1];
  double Mstmp2  = y * M[0];
  double Mstmp3  = x * M[2];
  double Mstmp4  = z * M[0];
  double Mstmp5  = y * M[1];
  double Mstmp6  = y * M[2];
  double Mstmp7  = z * M[1];
  double Mstmp8  = z * M[2];
  double Mstmp9  = x * M[3];
  double Mstmp10 = (1.0 / 2.0) * (x * x);
  double Mstmp11 = x * M[4];
  double Mstmp12 = y * M[3];
  double Mstmp13 = Mstmp0 * y;
  double Mstmp14 = x * M[5];
  double Mstmp15 = x * M[6];
  double Mstmp16 = y * M[4];
  double Mstmp17 = Mstmp1 * y;
  double Mstmp18 = (y * y);
  double Mstmp19 = (1.0 / 2.0) * M[0];
  double Mstmp20 = x * M[7];
  double Mstmp21 = y * M[5];
  double Mstmp22 = Mstmp3 * y;
  double Mstmp23 = x * M[8];
  double Mstmp24 = (z * z);
  double Mstmp25 = y * M[6];
  double Mstmp26 = (1.0 / 2.0) * Mstmp18;
  double Mstmp27 = y * M[7];
  double Mstmp28 = y * M[8];
  double Mstmp29 = (1.0 / 2.0) * Mstmp24;
  double Mstmp30 = (1.0 / 6.0) * (x * x * x);
  double Mstmp31 = (y * y * y);
  double Mstmp32 = (1.0 / 6.0) * M[0];
  double Mstmp33 = (z * z * z);
  double Mstmp34 = (1.0 / 6.0) * Mstmp31;
  double Mstmp35 = (1.0 / 6.0) * Mstmp33;
#pragma omp atomic
  Ms[0] += M[0];
#pragma omp atomic
  Ms[1] += M[1];
#pragma omp atomic
  Ms[2] += M[2];
#pragma omp atomic
  Ms[3] += Mstmp0 + M[3];
#pragma omp atomic
  Ms[4] += Mstmp1 + Mstmp2 + M[4];
#pragma omp atomic
  Ms[5] += Mstmp3 + Mstmp4 + M[5];
#pragma omp atomic
  Ms[6] += Mstmp5 + M[6];
#pragma omp atomic
  Ms[7] += Mstmp6 + Mstmp7 + M[7];
#pragma omp atomic
  Ms[8] += Mstmp8 + M[8];
#pragma omp atomic
  Ms[9] += Mstmp10 * M[0] + Mstmp9 + M[9];
#pragma omp atomic
  Ms[10] += Mstmp10 * M[1] + Mstmp11 + Mstmp12 + Mstmp13 + M[10];
#pragma omp atomic
  Ms[11] += Mstmp0 * z + Mstmp10 * M[2] + Mstmp14 + z * M[3] + M[11];
#pragma omp atomic
  Ms[12] += Mstmp15 + Mstmp16 + Mstmp17 + Mstmp18 * Mstmp19 + M[12];
#pragma omp atomic
  Ms[13] += Mstmp1 * z + Mstmp2 * z + Mstmp20 + Mstmp21 + Mstmp22 + z * M[4] + M[13];
#pragma omp atomic
  Ms[14] += Mstmp19 * Mstmp24 + Mstmp23 + Mstmp3 * z + z * M[5] + M[14];
#pragma omp atomic
  Ms[15] += Mstmp25 + Mstmp26 * M[1] + M[15];
#pragma omp atomic
  Ms[16] += Mstmp26 * M[2] + Mstmp27 + Mstmp5 * z + z * M[6] + M[16];
#pragma omp atomic
  Ms[17] += Mstmp28 + Mstmp29 * M[1] + Mstmp6 * z + z * M[7] + M[17];
#pragma omp atomic
  Ms[18] += Mstmp29 * M[2] + z * M[8] + M[18];
#pragma omp atomic
  Ms[19] += Mstmp10 * M[3] + Mstmp30 * M[0] + x * M[9] + M[19];
#pragma omp atomic
  Ms[20] += Mstmp10 * Mstmp2 + Mstmp10 * M[4] + Mstmp30 * M[1] + Mstmp9 * y + x * M[10] +
            y * M[9] + M[20];
#pragma omp atomic
  Ms[21] += Mstmp10 * Mstmp4 + Mstmp10 * M[5] + Mstmp30 * M[2] + Mstmp9 * z + x * M[11] +
            z * M[9] + M[21];
#pragma omp atomic
  Ms[22] += Mstmp0 * Mstmp26 + Mstmp10 * Mstmp5 + Mstmp10 * M[6] + Mstmp11 * y +
            Mstmp26 * M[3] + x * M[12] + y * M[10] + M[22];
#pragma omp atomic
  Ms[23] += Mstmp10 * Mstmp6 + Mstmp10 * Mstmp7 + Mstmp10 * M[7] + Mstmp11 * z +
            Mstmp12 * z + Mstmp13 * z + Mstmp14 * y + x * M[13] + y * M[11] + z * M[10] +
            M[23];
#pragma omp atomic
  Ms[24] += Mstmp0 * Mstmp29 + Mstmp10 * Mstmp8 + Mstmp10 * M[8] + Mstmp14 * z +
            Mstmp29 * M[3] + x * M[14] + z * M[11] + M[24];
#pragma omp atomic
  Ms[25] += Mstmp1 * Mstmp26 + Mstmp15 * y + Mstmp26 * M[4] + Mstmp31 * Mstmp32 +
            x * M[15] + y * M[12] + M[25];
#pragma omp atomic
  Ms[26] += Mstmp15 * z + Mstmp16 * z + Mstmp17 * z + Mstmp20 * y + Mstmp26 * Mstmp3 +
            Mstmp26 * Mstmp4 + Mstmp26 * M[5] + x * M[16] + y * M[13] + z * M[12] + M[26];
#pragma omp atomic
  Ms[27] += Mstmp1 * Mstmp29 + Mstmp2 * Mstmp29 + Mstmp20 * z + Mstmp21 * z +
            Mstmp22 * z + Mstmp23 * y + Mstmp29 * M[4] + x * M[17] + y * M[14] +
            z * M[13] + M[27];
#pragma omp atomic
  Ms[28] += Mstmp23 * z + Mstmp29 * Mstmp3 + Mstmp29 * M[5] + Mstmp32 * Mstmp33 +
            x * M[18] + z * M[14] + M[28];
#pragma omp atomic
  Ms[29] += Mstmp26 * M[6] + Mstmp34 * M[1] + y * M[15] + M[29];
#pragma omp atomic
  Ms[30] += Mstmp25 * z + Mstmp26 * Mstmp7 + Mstmp26 * M[7] + Mstmp34 * M[2] + y * M[16] +
            z * M[15] + M[30];
#pragma omp atomic
  Ms[31] += Mstmp26 * Mstmp8 + Mstmp26 * M[8] + Mstmp27 * z + Mstmp29 * Mstmp5 +
            Mstmp29 * M[6] + y * M[17] + z * M[16] + M[31];
#pragma omp atomic
  Ms[32] += Mstmp28 * z + Mstmp29 * Mstmp6 + Mstmp29 * M[7] + Mstmp35 * M[1] + y * M[18] +
            z * M[17] + M[32];
#pragma omp atomic
  Ms[33] += Mstmp29 * M[8] + Mstmp35 * M[2] + z * M[18] + M[33];
}

void field_m1_M2L_4(double x, double y, double z, double* M, double* L) {
  double R = sqrt(x * x + y * y + z * z);
  double D[34];
  double Dtmp0  = 1.0 * pow(R, -3.0);
  double Dtmp1  = -Dtmp0;
  double Dtmp2  = (x * x);
  double Dtmp3  = pow(R, -5.0);
  double Dtmp4  = 3.0 * Dtmp3;
  double Dtmp5  = Dtmp4 * x;
  double Dtmp6  = (y * y);
  double Dtmp7  = y * z;
  double Dtmp8  = 9.0 * Dtmp3;
  double Dtmp9  = -Dtmp8;
  double Dtmp10 = pow(R, -7.0);
  double Dtmp11 = 15.0 * Dtmp10;
  double Dtmp12 = Dtmp11 * Dtmp2;
  double Dtmp13 = -Dtmp4;
  double Dtmp14 = Dtmp12 + Dtmp13;
  double Dtmp15 = Dtmp11 * Dtmp6;
  double Dtmp16 = Dtmp13 + Dtmp15;
  double Dtmp17 = 1.0 * x;
  double Dtmp18 = 105.0 * pow(R, -9.0);
  double Dtmp19 = 90.0 * Dtmp10;
  double Dtmp20 = -45.0 * Dtmp10;
  double Dtmp21 = Dtmp18 * Dtmp2;
  double Dtmp22 = x * (Dtmp20 + Dtmp21);
  double Dtmp23 = -Dtmp11;
  double Dtmp24 = Dtmp18 * Dtmp6;
  double Dtmp25 = Dtmp20 + Dtmp24;
  D[0]          = -Dtmp0 * x;
  D[1]          = -Dtmp0 * y;
  D[2]          = -Dtmp0 * z;
  D[3]          = Dtmp1 + Dtmp2 * Dtmp4;
  D[4]          = Dtmp5 * y;
  D[5]          = Dtmp5 * z;
  D[6]          = Dtmp1 + Dtmp4 * Dtmp6;
  D[7]          = Dtmp4 * Dtmp7;
  D[8]          = -D[3] - D[6];
  D[9]          = -x * (Dtmp12 + Dtmp9);
  D[10]         = -Dtmp14 * y;
  D[11]         = -Dtmp14 * z;
  D[12]         = -Dtmp16 * Dtmp17;
  D[13]         = -Dtmp11 * Dtmp7 * x;
  D[14]         = -D[9] - D[12];
  D[15]         = -y * (Dtmp15 + Dtmp9);
  D[16]         = -Dtmp16 * z;
  D[17]         = -D[10] - D[15];
  D[18]         = -D[11] - D[16];
  D[19]         = Dtmp18 * (x * x * x * x) - Dtmp19 * Dtmp2 + Dtmp8;
  D[20]         = Dtmp22 * y;
  D[21]         = Dtmp22 * z;
  D[22]         = -Dtmp12 - Dtmp15 + Dtmp21 * Dtmp6 + Dtmp4;
  D[23]         = Dtmp7 * (Dtmp21 + Dtmp23);
  D[24]         = -D[19] - D[22];
  D[25]         = Dtmp17 * Dtmp25 * y;
  D[26]         = Dtmp17 * z * (Dtmp23 + Dtmp24);
  D[27]         = -D[20] - D[25];
  D[28]         = -D[21] - D[26];
  D[29]         = Dtmp18 * (y * y * y * y) - Dtmp19 * Dtmp6 + Dtmp8;
  D[30]         = Dtmp25 * Dtmp7;
  D[31]         = -D[22] - D[29];
  D[32]         = -D[23] - D[30];
  D[33]         = -D[24] - D[31];
#pragma omp atomic
  L[0] += D[0] * M[0] + D[1] * M[1] + D[2] * M[2] + D[3] * M[3] + D[4] * M[4] +
          D[5] * M[5] + D[6] * M[6] + D[7] * M[7] + D[8] * M[8] + D[9] * M[9] +
          D[10] * M[10] + D[11] * M[11] + D[12] * M[12] + D[13] * M[13] + D[14] * M[14] +
          D[15] * M[15] + D[16] * M[16] + D[17] * M[17] + D[18] * M[18] + D[19] * M[19] +
          D[20] * M[20] + D[21] * M[21] + D[22] * M[22] + D[23] * M[23] + D[24] * M[24] +
          D[25] * M[25] + D[26] * M[26] + D[27] * M[27] + D[28] * M[28] + D[29] * M[29] +
          D[30] * M[30] + D[31] * M[31] + D[32] * M[32] + D[33] * M[33];
#pragma omp atomic
  L[1] += D[3] * M[0] + D[4] * M[1] + D[5] * M[2] + D[9] * M[3] + D[10] * M[4] +
          D[11] * M[5] + D[12] * M[6] + D[13] * M[7] + D[14] * M[8] + D[19] * M[9] +
          D[20] * M[10] + D[21] * M[11] + D[22] * M[12] + D[23] * M[13] + D[24] * M[14] +
          D[25] * M[15] + D[26] * M[16] + D[27] * M[17] + D[28] * M[18];
#pragma omp atomic
  L[2] += D[4] * M[0] + D[6] * M[1] + D[7] * M[2] + D[10] * M[3] + D[12] * M[4] +
          D[13] * M[5] + D[15] * M[6] + D[16] * M[7] + D[17] * M[8] + D[20] * M[9] +
          D[22] * M[10] + D[23] * M[11] + D[25] * M[12] + D[26] * M[13] + D[27] * M[14] +
          D[29] * M[15] + D[30] * M[16] + D[31] * M[17] + D[32] * M[18];
#pragma omp atomic
  L[3] += D[5] * M[0] + D[7] * M[1] + D[8] * M[2] + D[11] * M[3] + D[13] * M[4] +
          D[14] * M[5] + D[16] * M[6] + D[17] * M[7] + D[18] * M[8] + D[21] * M[9] +
          D[23] * M[10] + D[24] * M[11] + D[26] * M[12] + D[27] * M[13] + D[28] * M[14] +
          D[30] * M[15] + D[31] * M[16] + D[32] * M[17] + D[33] * M[18];
#pragma omp atomic
  L[4] += D[9] * M[0] + D[10] * M[1] + D[11] * M[2] + D[19] * M[3] + D[20] * M[4] +
          D[21] * M[5] + D[22] * M[6] + D[23] * M[7] + D[24] * M[8];
#pragma omp atomic
  L[5] += D[10] * M[0] + D[12] * M[1] + D[13] * M[2] + D[20] * M[3] + D[22] * M[4] +
          D[23] * M[5] + D[25] * M[6] + D[26] * M[7] + D[27] * M[8];
#pragma omp atomic
  L[6] += D[11] * M[0] + D[13] * M[1] + D[14] * M[2] + D[21] * M[3] + D[23] * M[4] +
          D[24] * M[5] + D[26] * M[6] + D[27] * M[7] + D[28] * M[8];
#pragma omp atomic
  L[7] += D[12] * M[0] + D[15] * M[1] + D[16] * M[2] + D[22] * M[3] + D[25] * M[4] +
          D[26] * M[5] + D[29] * M[6] + D[30] * M[7] + D[31] * M[8];
#pragma omp atomic
  L[8] += D[13] * M[0] + D[16] * M[1] + D[17] * M[2] + D[23] * M[3] + D[26] * M[4] +
          D[27] * M[5] + D[30] * M[6] + D[31] * M[7] + D[32] * M[8];
#pragma omp atomic
  L[9] += D[14] * M[0] + D[17] * M[1] + D[18] * M[2] + D[24] * M[3] + D[27] * M[4] +
          D[28] * M[5] + D[31] * M[6] + D[32] * M[7] + D[33] * M[8];
#pragma omp atomic
  L[10] += D[19] * M[0] + D[20] * M[1] + D[21] * M[2];
#pragma omp atomic
  L[11] += D[20] * M[0] + D[22] * M[1] + D[23] * M[2];
#pragma omp atomic
  L[12] += D[21] * M[0] + D[23] * M[1] + D[24] * M[2];
#pragma omp atomic
  L[13] += D[22] * M[0] + D[25] * M[1] + D[26] * M[2];
#pragma omp atomic
  L[14] += D[23] * M[0] + D[26] * M[1] + D[27] * M[2];
#pragma omp atomic
  L[15] += D[24] * M[0] + D[27] * M[1] + D[28] * M[2];
#pragma omp atomic
  L[16] += D[25] * M[0] + D[29] * M[1] + D[30] * M[2];
#pragma omp atomic
  L[17] += D[26] * M[0] + D[30] * M[1] + D[31] * M[2];
#pragma omp atomic
  L[18] += D[27] * M[0] + D[31] * M[1] + D[32] * M[2];
#pragma omp atomic
  L[19] += D[28] * M[0] + D[32] * M[1] + D[33] * M[2];
}

void field_m1_L2L_4(double x, double y, double z, double* L, double* Ls) {
  double Lstmp0  = y * L[5];
  double Lstmp1  = z * L[6];
  double Lstmp2  = z * L[8];
  double Lstmp3  = z * L[14];
  double Lstmp4  = Lstmp3 * y;
  double Lstmp5  = (1.0 / 2.0) * (x * x);
  double Lstmp6  = (1.0 / 2.0) * (y * y);
  double Lstmp7  = (1.0 / 2.0) * (z * z);
  double Lstmp8  = x * L[13];
  double Lstmp9  = x * L[15];
  double Lstmp10 = y * L[11];
  double Lstmp11 = z * L[12];
  double Lstmp12 = y * L[18];
  double Lstmp13 = z * L[17];
  double Lstmp14 = y * L[13];
  double Lstmp15 = y * L[14];
  double Lstmp16 = z * L[15];
  double Lstmp17 = z * L[18];
#pragma omp atomic
  Ls[0] += Lstmp0 * x + Lstmp1 * x + Lstmp10 * Lstmp5 + Lstmp11 * Lstmp5 +
           Lstmp12 * Lstmp7 + Lstmp13 * Lstmp6 + Lstmp2 * y + Lstmp4 * x + Lstmp5 * L[4] +
           Lstmp6 * Lstmp8 + Lstmp6 * L[7] + Lstmp7 * Lstmp9 + Lstmp7 * L[9] +
           (1.0 / 6.0) * (x * x * x) * L[10] + x * L[1] +
           (1.0 / 6.0) * (y * y * y) * L[16] + y * L[2] +
           (1.0 / 6.0) * (z * z * z) * L[19] + z * L[3] + L[0];
#pragma omp atomic
  Ls[1] += Lstmp0 + Lstmp1 + Lstmp10 * x + Lstmp11 * x + Lstmp4 + Lstmp5 * L[10] +
           Lstmp6 * L[13] + Lstmp7 * L[15] + x * L[4] + L[1];
#pragma omp atomic
  Ls[2] += Lstmp13 * y + Lstmp14 * x + Lstmp2 + Lstmp3 * x + Lstmp5 * L[11] +
           Lstmp6 * L[16] + Lstmp7 * L[18] + x * L[5] + y * L[7] + L[2];
#pragma omp atomic
  Ls[3] += Lstmp15 * x + Lstmp16 * x + Lstmp17 * y + Lstmp5 * L[12] + Lstmp6 * L[17] +
           Lstmp7 * L[19] + x * L[6] + y * L[8] + z * L[9] + L[3];
#pragma omp atomic
  Ls[4] += Lstmp10 + Lstmp11 + x * L[10] + L[4];
#pragma omp atomic
  Ls[5] += Lstmp14 + Lstmp3 + x * L[11] + L[5];
#pragma omp atomic
  Ls[6] += Lstmp15 + Lstmp16 + x * L[12] + L[6];
#pragma omp atomic
  Ls[7] += Lstmp13 + Lstmp8 + y * L[16] + L[7];
#pragma omp atomic
  Ls[8] += Lstmp17 + x * L[14] + y * L[17] + L[8];
#pragma omp atomic
  Ls[9] += Lstmp12 + Lstmp9 + z * L[19] + L[9];
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

void field_m1_L2P_4(double x, double y, double z, double* L, double* F) {
  double Ftmp0 = x * y;
  double Ftmp1 = x * z;
  double Ftmp2 = y * z;
  double Ftmp3 = (1.0 / 2.0) * (x * x);
  double Ftmp4 = (1.0 / 2.0) * (y * y);
  double Ftmp5 = (1.0 / 2.0) * (z * z);
#pragma omp atomic
  F[0] += -Ftmp0 * L[11] - Ftmp1 * L[12] - Ftmp2 * L[14] - Ftmp3 * L[10] - Ftmp4 * L[13] -
          Ftmp5 * L[15] - x * L[4] - y * L[5] - z * L[6] - L[1];
#pragma omp atomic
  F[1] += -Ftmp0 * L[13] - Ftmp1 * L[14] - Ftmp2 * L[17] - Ftmp3 * L[11] - Ftmp4 * L[16] -
          Ftmp5 * L[18] - x * L[5] - y * L[7] - z * L[8] - L[2];
#pragma omp atomic
  F[2] += -Ftmp0 * L[14] - Ftmp1 * L[15] - Ftmp2 * L[18] - Ftmp3 * L[12] - Ftmp4 * L[17] -
          Ftmp5 * L[19] - x * L[6] - y * L[8] - z * L[9] - L[3];
}

void field_m1_M2P_4(double x, double y, double z, double* M, double* F) {
  double R       = sqrt(x * x + y * y + z * z);
  double Ftmp0   = 1.0 * pow(R, -3.0);
  double Ftmp1   = pow(R, -5.0);
  double Ftmp2   = 3.0 * Ftmp1;
  double Ftmp3   = Ftmp2 * M[4];
  double Ftmp4   = Ftmp2 * z;
  double Ftmp5   = y * z;
  double Ftmp6   = pow(R, -7.0);
  double Ftmp7   = 15.0 * Ftmp6;
  double Ftmp8   = Ftmp7 * M[13];
  double Ftmp9   = Ftmp2 * M[1];
  double Ftmp10  = x * y;
  double Ftmp11  = Ftmp4 * M[2];
  double Ftmp12  = (x * x);
  double Ftmp13  = Ftmp2 * M[0];
  double Ftmp14  = y * M[7];
  double Ftmp15  = x * z;
  double Ftmp16  = Ftmp15 * Ftmp7;
  double Ftmp17  = Ftmp12 * Ftmp7;
  double Ftmp18  = y * M[4];
  double Ftmp19  = pow(R, -9.0);
  double Ftmp20  = 105.0 * Ftmp19;
  double Ftmp21  = Ftmp12 * Ftmp20;
  double Ftmp22  = -9.0 * Ftmp1;
  double Ftmp23  = Ftmp17 + Ftmp22;
  double Ftmp24  = -Ftmp2;
  double Ftmp25  = (y * y);
  double Ftmp26  = Ftmp25 * Ftmp7;
  double Ftmp27  = Ftmp24 + Ftmp26;
  double Ftmp28  = (z * z);
  double Ftmp29  = Ftmp28 * Ftmp7;
  double Ftmp30  = Ftmp24 + Ftmp29;
  double Ftmp31  = Ftmp27 * M[6];
  double Ftmp32  = Ftmp30 * M[8];
  double Ftmp33  = 45.0 * Ftmp6;
  double Ftmp34  = -Ftmp33;
  double Ftmp35  = Ftmp21 + Ftmp34;
  double Ftmp36  = Ftmp35 * M[20];
  double Ftmp37  = Ftmp20 * Ftmp25;
  double Ftmp38  = Ftmp34 + Ftmp37;
  double Ftmp39  = Ftmp38 * y;
  double Ftmp40  = 1.0 * M[25];
  double Ftmp41  = 35.0 * Ftmp19;
  double Ftmp42  = 3.0 * M[27];
  double Ftmp43  = Ftmp42 * (Ftmp28 * Ftmp41 - 5.0 * Ftmp6);
  double Ftmp44  = Ftmp35 * M[21];
  double Ftmp45  = 1.0 * z;
  double Ftmp46  = -Ftmp7;
  double Ftmp47  = Ftmp37 + Ftmp46;
  double Ftmp48  = Ftmp47 * M[26];
  double Ftmp49  = Ftmp20 * Ftmp28;
  double Ftmp50  = Ftmp34 + Ftmp49;
  double Ftmp51  = Ftmp45 * Ftmp50;
  double Ftmp52  = Ftmp35 * x;
  double Ftmp53  = Ftmp52 * y;
  double Ftmp54  = Ftmp39 * M[15];
  double Ftmp55  = Ftmp46 + Ftmp49;
  double Ftmp56  = 1.0 * Ftmp55 * M[17];
  double Ftmp57  = Ftmp52 * z;
  double Ftmp58  = Ftmp47 * M[16];
  double Ftmp59  = Ftmp50 * M[18];
  double Ftmp60  = 315.0 * Ftmp19;
  double Ftmp61  = -Ftmp60;
  double Ftmp62  = pow(R, -11.0);
  double Ftmp63  = 945.0 * Ftmp62;
  double Ftmp64  = Ftmp12 * Ftmp63;
  double Ftmp65  = Ftmp61 + Ftmp64;
  double Ftmp66  = Ftmp15 * y;
  double Ftmp67  = Ftmp25 * Ftmp63;
  double Ftmp68  = Ftmp61 + Ftmp67;
  double Ftmp69  = Ftmp68 * M[30];
  double Ftmp70  = -75.0 * Ftmp6;
  double Ftmp71  = 1.0 * Ftmp12;
  double Ftmp72  = Ftmp47 * M[12];
  double Ftmp73  = Ftmp55 * M[14];
  double Ftmp74  = -525.0 * Ftmp19;
  double Ftmp75  = Ftmp12 * (Ftmp64 + Ftmp74);
  double Ftmp76  = y * M[20];
  double Ftmp77  = z * M[21];
  double Ftmp78  = y * M[32];
  double Ftmp79  = Ftmp28 * Ftmp63;
  double Ftmp80  = Ftmp61 + Ftmp79;
  double Ftmp81  = Ftmp45 * Ftmp80;
  double Ftmp82  = Ftmp12 * y;
  double Ftmp83  = Ftmp40 * Ftmp68;
  double Ftmp84  = 315.0 * Ftmp28 * Ftmp62;
  double Ftmp85  = Ftmp42 * (-Ftmp41 + Ftmp84);
  double Ftmp86  = Ftmp12 * Ftmp45;
  double Ftmp87  = -Ftmp20;
  double Ftmp88  = (Ftmp67 + Ftmp87) * M[26];
  double Ftmp89  = Ftmp80 * M[28];
  double Ftmp90  = 225.0 * Ftmp6;
  double Ftmp91  = Ftmp63 * (x * x * x * x);
  double Ftmp92  = 1050.0 * Ftmp19;
  double Ftmp93  = Ftmp63 * (y * y * y * y);
  double Ftmp94  = 630.0 * Ftmp19;
  double Ftmp95  = (-Ftmp25 * Ftmp94 + Ftmp33 + Ftmp93) * M[29];
  double Ftmp96  = Ftmp63 * (z * z * z * z);
  double Ftmp97  = (-Ftmp28 * Ftmp94 + Ftmp33 + Ftmp96) * M[33];
  double Ftmp98  = -Ftmp25 * Ftmp60;
  double Ftmp99  = Ftmp25 * Ftmp64;
  double Ftmp100 = -Ftmp21;
  double Ftmp101 = Ftmp100 + Ftmp33;
  double Ftmp102 = -Ftmp28 * Ftmp60;
  double Ftmp103 = Ftmp28 * Ftmp64;
  double Ftmp104 = -Ftmp49;
  double Ftmp105 = Ftmp104 + Ftmp7;
  double Ftmp106 = -Ftmp37;
  double Ftmp107 = Ftmp28 * Ftmp67;
  double Ftmp108 = Ftmp106 + Ftmp107;
  double Ftmp109 = Ftmp17 + Ftmp24;
  double Ftmp110 = Ftmp22 + Ftmp26;
  double Ftmp111 = Ftmp109 * M[3];
  double Ftmp112 = Ftmp40 * x;
  double Ftmp113 = Ftmp21 + Ftmp46;
  double Ftmp114 = Ftmp113 * M[23];
  double Ftmp115 = Ftmp38 * M[30];
  double Ftmp116 = 1.0 * x;
  double Ftmp117 = Ftmp113 * M[11];
  double Ftmp118 = Ftmp65 * x;
  double Ftmp119 = Ftmp113 * M[10];
  double Ftmp120 = Ftmp25 * z;
  double Ftmp121 = (Ftmp64 + Ftmp87) * M[23];
  double Ftmp122 = Ftmp67 + Ftmp74;
  double Ftmp123 = Ftmp10 * Ftmp45;
  double Ftmp124 = (-Ftmp12 * Ftmp94 + Ftmp33 + Ftmp91) * M[19];
  double Ftmp125 = Ftmp106 + Ftmp99;
  double Ftmp126 = -Ftmp12 * Ftmp60 + Ftmp33;
  double Ftmp127 = x * M[5];
  double Ftmp128 = Ftmp22 + Ftmp29;
  double Ftmp129 = Ftmp116 * M[28];
  double Ftmp130 = 1.0 * Ftmp78;
  double Ftmp131 = Ftmp28 * y;
  double Ftmp132 = Ftmp28 * (Ftmp74 + Ftmp79);
#pragma omp atomic
  F[0] += Ftmp0 * M[0] - Ftmp10 * Ftmp56 - Ftmp10 * Ftmp9 - Ftmp11 * x - Ftmp12 * Ftmp13 -
          Ftmp12 * (Ftmp21 + Ftmp70) * M[9] + Ftmp14 * Ftmp16 - Ftmp15 * Ftmp58 -
          Ftmp15 * Ftmp59 + Ftmp17 * Ftmp18 + Ftmp17 * z * M[5] - Ftmp21 * Ftmp5 * M[13] +
          Ftmp23 * x * M[3] + Ftmp23 * M[9] + Ftmp27 * M[12] - Ftmp3 * y +
          Ftmp30 * M[14] + Ftmp31 * x + Ftmp32 * x - Ftmp36 * y - Ftmp39 * Ftmp40 -
          Ftmp4 * M[5] - Ftmp43 * y - Ftmp44 * z - Ftmp45 * Ftmp48 + Ftmp5 * Ftmp8 -
          Ftmp51 * M[28] - Ftmp53 * M[10] - Ftmp54 * x - Ftmp57 * M[11] +
          Ftmp65 * Ftmp66 * M[23] + Ftmp66 * Ftmp69 - Ftmp71 * Ftmp72 - Ftmp71 * Ftmp73 +
          Ftmp75 * Ftmp76 + Ftmp75 * Ftmp77 + Ftmp78 * Ftmp81 * x + Ftmp82 * Ftmp83 +
          Ftmp82 * Ftmp85 + Ftmp86 * Ftmp88 + Ftmp86 * Ftmp89 + Ftmp95 * x + Ftmp97 * x +
          x * (Ftmp105 + Ftmp108) * M[31] + x * (Ftmp101 + Ftmp102 + Ftmp103) * M[24] +
          x * (Ftmp101 + Ftmp98 + Ftmp99) * M[22] +
          x * (-Ftmp12 * Ftmp92 + Ftmp90 + Ftmp91) * M[19];
#pragma omp atomic
  F[1] += Ftmp0 * M[1] - Ftmp10 * Ftmp13 + Ftmp109 * M[10] - Ftmp11 * y +
          Ftmp110 * y * M[6] + Ftmp110 * M[15] + Ftmp111 * y +
          Ftmp112 * Ftmp122 * Ftmp25 - Ftmp112 * Ftmp38 - Ftmp114 * z - Ftmp115 * z -
          Ftmp116 * Ftmp39 * M[12] - Ftmp116 * Ftmp73 * y - Ftmp117 * Ftmp5 +
          Ftmp118 * Ftmp25 * M[20] + Ftmp118 * Ftmp77 * y - Ftmp119 * Ftmp25 +
          Ftmp120 * Ftmp121 + Ftmp120 * Ftmp122 * M[30] + Ftmp123 * Ftmp68 * M[26] +
          Ftmp123 * Ftmp89 + Ftmp124 * y - Ftmp15 * Ftmp37 * M[13] + Ftmp15 * Ftmp8 -
          Ftmp25 * Ftmp56 + Ftmp25 * Ftmp81 * M[32] + Ftmp25 * Ftmp85 * x -
          Ftmp25 * Ftmp9 - Ftmp25 * (Ftmp37 + Ftmp70) * M[15] + Ftmp26 * x * M[4] +
          Ftmp26 * z * M[7] - Ftmp3 * x + Ftmp30 * M[17] + Ftmp32 * y - Ftmp36 * x -
          Ftmp39 * z * M[16] - Ftmp4 * M[7] - Ftmp43 * x - Ftmp5 * Ftmp59 -
          Ftmp51 * M[32] - Ftmp53 * M[9] + Ftmp66 * Ftmp7 * M[5] + Ftmp97 * y +
          y * (Ftmp125 + Ftmp126) * M[22] + y * (Ftmp100 + Ftmp103 + Ftmp105) * M[24] +
          y * (Ftmp102 + Ftmp108 + Ftmp33) * M[31] +
          y * (-Ftmp25 * Ftmp92 + Ftmp90 + Ftmp93) * M[29];
#pragma omp atomic
  F[2] += Ftmp0 * M[2] - Ftmp10 * Ftmp49 * M[13] + Ftmp10 * Ftmp8 + Ftmp109 * M[11] +
          Ftmp111 * z - Ftmp114 * y - Ftmp115 * y + Ftmp116 * Ftmp28 * Ftmp88 -
          Ftmp116 * Ftmp48 - Ftmp117 * Ftmp28 + Ftmp118 * Ftmp28 * M[21] -
          Ftmp119 * Ftmp5 + Ftmp121 * Ftmp131 + Ftmp124 * z - Ftmp127 * Ftmp2 +
          Ftmp127 * Ftmp29 + Ftmp128 * z * M[8] + Ftmp128 * M[18] + Ftmp129 * Ftmp132 -
          Ftmp129 * Ftmp50 + Ftmp130 * Ftmp132 - Ftmp130 * Ftmp50 + Ftmp131 * Ftmp69 -
          Ftmp14 * Ftmp2 + Ftmp14 * Ftmp29 + Ftmp15 * Ftmp65 * Ftmp76 + Ftmp16 * Ftmp18 -
          Ftmp2 * Ftmp28 * M[2] + Ftmp27 * M[16] - Ftmp28 * Ftmp58 -
          Ftmp28 * (Ftmp49 + Ftmp70) * M[18] + Ftmp31 * z - Ftmp4 * x * M[0] -
          Ftmp4 * y * M[1] + Ftmp42 * Ftmp66 * (Ftmp84 + Ftmp87) - Ftmp44 * x -
          Ftmp45 * Ftmp72 * x - Ftmp51 * x * M[14] - Ftmp51 * y * M[17] - Ftmp54 * z -
          Ftmp57 * M[9] + Ftmp66 * Ftmp83 + Ftmp95 * z +
          z * (Ftmp100 + Ftmp125 + Ftmp7) * M[22] +
          z * (Ftmp103 + Ftmp104 + Ftmp126) * M[24] +
          z * (-Ftmp28 * Ftmp92 + Ftmp90 + Ftmp96) * M[33] +
          z * (Ftmp104 + Ftmp107 + Ftmp33 + Ftmp98) * M[31];
}

void field_m1_P2M_5(double x, double y, double z, double q, double* M) {
  double Mtmp0  = q * x;
  double Mtmp1  = q * y;
  double Mtmp2  = q * z;
  double Mtmp3  = (x * x);
  double Mtmp4  = (1.0 / 2.0) * q;
  double Mtmp5  = Mtmp0 * y;
  double Mtmp6  = Mtmp0 * z;
  double Mtmp7  = (y * y);
  double Mtmp8  = Mtmp1 * z;
  double Mtmp9  = (z * z);
  double Mtmp10 = (x * x * x);
  double Mtmp11 = (1.0 / 6.0) * q;
  double Mtmp12 = (1.0 / 2.0) * Mtmp3;
  double Mtmp13 = (1.0 / 2.0) * Mtmp0;
  double Mtmp14 = (y * y * y);
  double Mtmp15 = (1.0 / 2.0) * Mtmp7;
  double Mtmp16 = (1.0 / 2.0) * Mtmp9;
  double Mtmp17 = (z * z * z);
  double Mtmp18 = (x * x * x * x);
  double Mtmp19 = (1.0 / 24.0) * q;
  double Mtmp20 = (1.0 / 6.0) * Mtmp10;
  double Mtmp21 = Mtmp7 * q;
  double Mtmp22 = (1.0 / 4.0) * Mtmp3;
  double Mtmp23 = Mtmp9 * q;
  double Mtmp24 = (1.0 / 6.0) * Mtmp0;
  double Mtmp25 = (y * y * y * y);
  double Mtmp26 = (1.0 / 6.0) * Mtmp14;
  double Mtmp27 = (1.0 / 4.0) * Mtmp9;
  double Mtmp28 = (1.0 / 6.0) * Mtmp17;
  double Mtmp29 = (z * z * z * z);
  double Mtmp30 = (1.0 / 120.0) * q;
  double Mtmp31 = (1.0 / 24.0) * Mtmp18;
  double Mtmp32 = (1.0 / 12.0) * Mtmp10;
  double Mtmp33 = (1.0 / 12.0) * Mtmp14;
  double Mtmp34 = Mtmp3 * q;
  double Mtmp35 = (1.0 / 12.0) * Mtmp17;
  double Mtmp36 = (1.0 / 24.0) * Mtmp0;
  M[0] += -Mtmp0;
  M[1] += -Mtmp1;
  M[2] += -Mtmp2;
  M[3] += Mtmp3 * Mtmp4;
  M[4] += Mtmp5;
  M[5] += Mtmp6;
  M[6] += Mtmp4 * Mtmp7;
  M[7] += Mtmp8;
  M[8] += Mtmp4 * Mtmp9;
  M[9] += -Mtmp10 * Mtmp11;
  M[10] += -Mtmp1 * Mtmp12;
  M[11] += -Mtmp12 * Mtmp2;
  M[12] += -Mtmp13 * Mtmp7;
  M[13] += -Mtmp5 * z;
  M[14] += -Mtmp13 * Mtmp9;
  M[15] += -Mtmp11 * Mtmp14;
  M[16] += -Mtmp15 * Mtmp2;
  M[17] += -Mtmp1 * Mtmp16;
  M[18] += -Mtmp11 * Mtmp17;
  M[19] += Mtmp18 * Mtmp19;
  M[20] += Mtmp1 * Mtmp20;
  M[21] += Mtmp2 * Mtmp20;
  M[22] += Mtmp21 * Mtmp22;
  M[23] += Mtmp12 * Mtmp8;
  M[24] += Mtmp22 * Mtmp23;
  M[25] += Mtmp14 * Mtmp24;
  M[26] += Mtmp15 * Mtmp6;
  M[27] += Mtmp16 * Mtmp5;
  M[28] += Mtmp17 * Mtmp24;
  M[29] += Mtmp19 * Mtmp25;
  M[30] += Mtmp2 * Mtmp26;
  M[31] += Mtmp21 * Mtmp27;
  M[32] += Mtmp1 * Mtmp28;
  M[33] += Mtmp19 * Mtmp29;
  M[34] += -Mtmp30 * (x * x * x * x * x);
  M[35] += -Mtmp1 * Mtmp31;
  M[36] += -Mtmp2 * Mtmp31;
  M[37] += -Mtmp21 * Mtmp32;
  M[38] += -Mtmp20 * Mtmp8;
  M[39] += -Mtmp23 * Mtmp32;
  M[40] += -Mtmp33 * Mtmp34;
  M[41] += -Mtmp2 * Mtmp22 * Mtmp7;
  M[42] += -Mtmp1 * Mtmp22 * Mtmp9;
  M[43] += -Mtmp34 * Mtmp35;
  M[44] += -Mtmp25 * Mtmp36;
  M[45] += -Mtmp26 * Mtmp6;
  M[46] += -Mtmp0 * Mtmp27 * Mtmp7;
  M[47] += -Mtmp28 * Mtmp5;
  M[48] += -Mtmp29 * Mtmp36;
  M[49] += -Mtmp30 * (y * y * y * y * y);
  M[50] += -1.0 / 24.0 * Mtmp2 * Mtmp25;
  M[51] += -Mtmp23 * Mtmp33;
  M[52] += -Mtmp21 * Mtmp35;
  M[53] += -1.0 / 24.0 * Mtmp1 * Mtmp29;
  M[54] += -Mtmp30 * (z * z * z * z * z);
}
void field_m1_M2M_5(double x, double y, double z, double* M, double* Ms) {
  double Mstmp0  = x * M[0];
  double Mstmp1  = x * M[1];
  double Mstmp2  = y * M[0];
  double Mstmp3  = x * M[2];
  double Mstmp4  = z * M[0];
  double Mstmp5  = y * M[1];
  double Mstmp6  = y * M[2];
  double Mstmp7  = z * M[1];
  double Mstmp8  = z * M[2];
  double Mstmp9  = x * M[3];
  double Mstmp10 = (x * x);
  double Mstmp11 = (1.0 / 2.0) * Mstmp10;
  double Mstmp12 = x * M[4];
  double Mstmp13 = y * M[3];
  double Mstmp14 = Mstmp0 * y;
  double Mstmp15 = x * M[5];
  double Mstmp16 = z * M[3];
  double Mstmp17 = Mstmp0 * z;
  double Mstmp18 = x * M[6];
  double Mstmp19 = y * M[4];
  double Mstmp20 = Mstmp1 * y;
  double Mstmp21 = (y * y);
  double Mstmp22 = (1.0 / 2.0) * M[0];
  double Mstmp23 = x * M[7];
  double Mstmp24 = y * M[5];
  double Mstmp25 = z * M[4];
  double Mstmp26 = Mstmp3 * y;
  double Mstmp27 = Mstmp1 * z;
  double Mstmp28 = Mstmp2 * z;
  double Mstmp29 = x * M[8];
  double Mstmp30 = z * M[5];
  double Mstmp31 = Mstmp3 * z;
  double Mstmp32 = (z * z);
  double Mstmp33 = y * M[6];
  double Mstmp34 = (1.0 / 2.0) * Mstmp21;
  double Mstmp35 = y * M[7];
  double Mstmp36 = z * M[6];
  double Mstmp37 = Mstmp5 * z;
  double Mstmp38 = y * M[8];
  double Mstmp39 = z * M[7];
  double Mstmp40 = Mstmp6 * z;
  double Mstmp41 = (1.0 / 2.0) * Mstmp32;
  double Mstmp42 = z * M[8];
  double Mstmp43 = x * M[9];
  double Mstmp44 = (1.0 / 6.0) * (x * x * x);
  double Mstmp45 = x * M[10];
  double Mstmp46 = y * M[9];
  double Mstmp47 = Mstmp9 * y;
  double Mstmp48 = x * M[11];
  double Mstmp49 = x * M[12];
  double Mstmp50 = y * M[10];
  double Mstmp51 = Mstmp12 * y;
  double Mstmp52 = x * M[13];
  double Mstmp53 = y * M[11];
  double Mstmp54 = Mstmp15 * y;
  double Mstmp55 = x * M[14];
  double Mstmp56 = x * M[15];
  double Mstmp57 = y * M[12];
  double Mstmp58 = Mstmp18 * y;
  double Mstmp59 = (y * y * y);
  double Mstmp60 = (1.0 / 6.0) * M[0];
  double Mstmp61 = x * M[16];
  double Mstmp62 = y * M[13];
  double Mstmp63 = Mstmp23 * y;
  double Mstmp64 = x * M[17];
  double Mstmp65 = y * M[14];
  double Mstmp66 = Mstmp29 * y;
  double Mstmp67 = x * M[18];
  double Mstmp68 = (z * z * z);
  double Mstmp69 = y * M[15];
  double Mstmp70 = (1.0 / 6.0) * Mstmp59;
  double Mstmp71 = y * M[16];
  double Mstmp72 = y * M[17];
  double Mstmp73 = y * M[18];
  double Mstmp74 = (1.0 / 6.0) * Mstmp68;
  double Mstmp75 = (1.0 / 24.0) * (x * x * x * x);
  double Mstmp76 = (1.0 / 4.0) * Mstmp10;
  double Mstmp77 = Mstmp76 * M[0];
  double Mstmp78 = Mstmp21 * Mstmp76;
  double Mstmp79 = Mstmp32 * Mstmp76;
  double Mstmp80 = (y * y * y * y);
  double Mstmp81 = (1.0 / 24.0) * M[0];
  double Mstmp82 = (1.0 / 4.0) * Mstmp21 * Mstmp32;
  double Mstmp83 = (z * z * z * z);
  double Mstmp84 = (1.0 / 24.0) * Mstmp80;
  double Mstmp85 = (1.0 / 24.0) * Mstmp83;
#pragma omp atomic
  Ms[0] += M[0];
#pragma omp atomic
  Ms[1] += M[1];
#pragma omp atomic
  Ms[2] += M[2];
#pragma omp atomic
  Ms[3] += Mstmp0 + M[3];
#pragma omp atomic
  Ms[4] += Mstmp1 + Mstmp2 + M[4];
#pragma omp atomic
  Ms[5] += Mstmp3 + Mstmp4 + M[5];
#pragma omp atomic
  Ms[6] += Mstmp5 + M[6];
#pragma omp atomic
  Ms[7] += Mstmp6 + Mstmp7 + M[7];
#pragma omp atomic
  Ms[8] += Mstmp8 + M[8];
#pragma omp atomic
  Ms[9] += Mstmp11 * M[0] + Mstmp9 + M[9];
#pragma omp atomic
  Ms[10] += Mstmp11 * M[1] + Mstmp12 + Mstmp13 + Mstmp14 + M[10];
#pragma omp atomic
  Ms[11] += Mstmp11 * M[2] + Mstmp15 + Mstmp16 + Mstmp17 + M[11];
#pragma omp atomic
  Ms[12] += Mstmp18 + Mstmp19 + Mstmp20 + Mstmp21 * Mstmp22 + M[12];
#pragma omp atomic
  Ms[13] += Mstmp23 + Mstmp24 + Mstmp25 + Mstmp26 + Mstmp27 + Mstmp28 + M[13];
#pragma omp atomic
  Ms[14] += Mstmp22 * Mstmp32 + Mstmp29 + Mstmp30 + Mstmp31 + M[14];
#pragma omp atomic
  Ms[15] += Mstmp33 + Mstmp34 * M[1] + M[15];
#pragma omp atomic
  Ms[16] += Mstmp34 * M[2] + Mstmp35 + Mstmp36 + Mstmp37 + M[16];
#pragma omp atomic
  Ms[17] += Mstmp38 + Mstmp39 + Mstmp40 + Mstmp41 * M[1] + M[17];
#pragma omp atomic
  Ms[18] += Mstmp41 * M[2] + Mstmp42 + M[18];
#pragma omp atomic
  Ms[19] += Mstmp11 * M[3] + Mstmp43 + Mstmp44 * M[0] + M[19];
#pragma omp atomic
  Ms[20] += Mstmp11 * Mstmp2 + Mstmp11 * M[4] + Mstmp44 * M[1] + Mstmp45 + Mstmp46 +
            Mstmp47 + M[20];
#pragma omp atomic
  Ms[21] += Mstmp11 * Mstmp4 + Mstmp11 * M[5] + Mstmp44 * M[2] + Mstmp48 + Mstmp9 * z +
            z * M[9] + M[21];
#pragma omp atomic
  Ms[22] += Mstmp0 * Mstmp34 + Mstmp11 * Mstmp5 + Mstmp11 * M[6] + Mstmp34 * M[3] +
            Mstmp49 + Mstmp50 + Mstmp51 + M[22];
#pragma omp atomic
  Ms[23] += Mstmp11 * Mstmp6 + Mstmp11 * Mstmp7 + Mstmp11 * M[7] + Mstmp12 * z +
            Mstmp13 * z + Mstmp14 * z + Mstmp52 + Mstmp53 + Mstmp54 + z * M[10] + M[23];
#pragma omp atomic
  Ms[24] += Mstmp0 * Mstmp41 + Mstmp11 * Mstmp8 + Mstmp11 * M[8] + Mstmp15 * z +
            Mstmp41 * M[3] + Mstmp55 + z * M[11] + M[24];
#pragma omp atomic
  Ms[25] += Mstmp1 * Mstmp34 + Mstmp34 * M[4] + Mstmp56 + Mstmp57 + Mstmp58 +
            Mstmp59 * Mstmp60 + M[25];
#pragma omp atomic
  Ms[26] += Mstmp18 * z + Mstmp19 * z + Mstmp20 * z + Mstmp3 * Mstmp34 +
            Mstmp34 * Mstmp4 + Mstmp34 * M[5] + Mstmp61 + Mstmp62 + Mstmp63 + z * M[12] +
            M[26];
#pragma omp atomic
  Ms[27] += Mstmp1 * Mstmp41 + Mstmp2 * Mstmp41 + Mstmp23 * z + Mstmp24 * z +
            Mstmp26 * z + Mstmp41 * M[4] + Mstmp64 + Mstmp65 + Mstmp66 + z * M[13] +
            M[27];
#pragma omp atomic
  Ms[28] += Mstmp29 * z + Mstmp3 * Mstmp41 + Mstmp41 * M[5] + Mstmp60 * Mstmp68 +
            Mstmp67 + z * M[14] + M[28];
#pragma omp atomic
  Ms[29] += Mstmp34 * M[6] + Mstmp69 + Mstmp70 * M[1] + M[29];
#pragma omp atomic
  Ms[30] += Mstmp33 * z + Mstmp34 * Mstmp7 + Mstmp34 * M[7] + Mstmp70 * M[2] + Mstmp71 +
            z * M[15] + M[30];
#pragma omp atomic
  Ms[31] += Mstmp34 * Mstmp8 + Mstmp34 * M[8] + Mstmp35 * z + Mstmp41 * Mstmp5 +
            Mstmp41 * M[6] + Mstmp72 + z * M[16] + M[31];
#pragma omp atomic
  Ms[32] += Mstmp38 * z + Mstmp41 * Mstmp6 + Mstmp41 * M[7] + Mstmp73 + Mstmp74 * M[1] +
            z * M[17] + M[32];
#pragma omp atomic
  Ms[33] += Mstmp41 * M[8] + Mstmp74 * M[2] + z * M[18] + M[33];
#pragma omp atomic
  Ms[34] += Mstmp11 * M[9] + Mstmp44 * M[3] + Mstmp75 * M[0] + x * M[19] + M[34];
#pragma omp atomic
  Ms[35] += Mstmp11 * Mstmp13 + Mstmp11 * M[10] + Mstmp2 * Mstmp44 + Mstmp43 * y +
            Mstmp44 * M[4] + Mstmp75 * M[1] + x * M[20] + y * M[19] + M[35];
#pragma omp atomic
  Ms[36] += Mstmp11 * Mstmp16 + Mstmp11 * M[11] + Mstmp4 * Mstmp44 + Mstmp43 * z +
            Mstmp44 * M[5] + Mstmp75 * M[2] + x * M[21] + z * M[19] + M[36];
#pragma omp atomic
  Ms[37] += Mstmp11 * Mstmp19 + Mstmp11 * M[12] + Mstmp21 * Mstmp77 + Mstmp34 * Mstmp9 +
            Mstmp34 * M[9] + Mstmp44 * Mstmp5 + Mstmp44 * M[6] + Mstmp45 * y + x * M[22] +
            y * M[20] + M[37];
#pragma omp atomic
  Ms[38] += Mstmp11 * Mstmp24 + Mstmp11 * Mstmp25 + Mstmp11 * Mstmp28 + Mstmp11 * M[13] +
            Mstmp44 * Mstmp6 + Mstmp44 * Mstmp7 + Mstmp44 * M[7] + Mstmp45 * z +
            Mstmp46 * z + Mstmp47 * z + Mstmp48 * y + x * M[23] + y * M[21] + z * M[20] +
            M[38];
#pragma omp atomic
  Ms[39] += Mstmp11 * Mstmp30 + Mstmp11 * M[14] + Mstmp32 * Mstmp77 + Mstmp41 * Mstmp9 +
            Mstmp41 * M[9] + Mstmp44 * Mstmp8 + Mstmp44 * M[8] + Mstmp48 * z + x * M[24] +
            z * M[21] + M[39];
#pragma omp atomic
  Ms[40] += Mstmp0 * Mstmp70 + Mstmp11 * Mstmp33 + Mstmp11 * M[15] + Mstmp12 * Mstmp34 +
            Mstmp34 * M[10] + Mstmp49 * y + Mstmp70 * M[3] + Mstmp78 * M[1] + x * M[25] +
            y * M[22] + M[40];
#pragma omp atomic
  Ms[41] += Mstmp11 * Mstmp35 + Mstmp11 * Mstmp36 + Mstmp11 * Mstmp37 + Mstmp11 * M[16] +
            Mstmp15 * Mstmp34 + Mstmp16 * Mstmp34 + Mstmp17 * Mstmp34 + Mstmp34 * M[11] +
            Mstmp49 * z + Mstmp50 * z + Mstmp51 * z + Mstmp52 * y + Mstmp78 * M[2] +
            x * M[26] + y * M[23] + z * M[22] + M[41];
#pragma omp atomic
  Ms[42] += Mstmp11 * Mstmp38 + Mstmp11 * Mstmp39 + Mstmp11 * Mstmp40 + Mstmp11 * M[17] +
            Mstmp12 * Mstmp41 + Mstmp13 * Mstmp41 + Mstmp14 * Mstmp41 + Mstmp41 * M[10] +
            Mstmp52 * z + Mstmp53 * z + Mstmp54 * z + Mstmp55 * y + Mstmp79 * M[1] +
            x * M[27] + y * M[24] + z * M[23] + M[42];
#pragma omp atomic
  Ms[43] += Mstmp0 * Mstmp74 + Mstmp11 * Mstmp42 + Mstmp11 * M[18] + Mstmp15 * Mstmp41 +
            Mstmp41 * M[11] + Mstmp55 * z + Mstmp74 * M[3] + Mstmp79 * M[2] + x * M[28] +
            z * M[24] + M[43];
#pragma omp atomic
  Ms[44] += Mstmp1 * Mstmp70 + Mstmp18 * Mstmp34 + Mstmp34 * M[12] + Mstmp56 * y +
            Mstmp70 * M[4] + Mstmp80 * Mstmp81 + x * M[29] + y * M[25] + M[44];
#pragma omp atomic
  Ms[45] += Mstmp23 * Mstmp34 + Mstmp25 * Mstmp34 + Mstmp27 * Mstmp34 + Mstmp3 * Mstmp70 +
            Mstmp34 * M[13] + Mstmp4 * Mstmp70 + Mstmp56 * z + Mstmp57 * z + Mstmp58 * z +
            Mstmp61 * y + Mstmp70 * M[5] + x * M[30] + y * M[26] + z * M[25] + M[45];
#pragma omp atomic
  Ms[46] += Mstmp18 * Mstmp41 + Mstmp19 * Mstmp41 + Mstmp20 * Mstmp41 +
            Mstmp29 * Mstmp34 + Mstmp30 * Mstmp34 + Mstmp31 * Mstmp34 + Mstmp34 * M[14] +
            Mstmp41 * M[12] + Mstmp61 * z + Mstmp62 * z + Mstmp63 * z + Mstmp64 * y +
            Mstmp82 * M[0] + x * M[31] + y * M[27] + z * M[26] + M[46];
#pragma omp atomic
  Ms[47] += Mstmp1 * Mstmp74 + Mstmp2 * Mstmp74 + Mstmp23 * Mstmp41 + Mstmp24 * Mstmp41 +
            Mstmp26 * Mstmp41 + Mstmp41 * M[13] + Mstmp64 * z + Mstmp65 * z +
            Mstmp66 * z + Mstmp67 * y + Mstmp74 * M[4] + x * M[32] + y * M[28] +
            z * M[27] + M[47];
#pragma omp atomic
  Ms[48] += Mstmp29 * Mstmp41 + Mstmp3 * Mstmp74 + Mstmp41 * M[14] + Mstmp67 * z +
            Mstmp74 * M[5] + Mstmp81 * Mstmp83 + x * M[33] + z * M[28] + M[48];
#pragma omp atomic
  Ms[49] += Mstmp34 * M[15] + Mstmp70 * M[6] + Mstmp84 * M[1] + y * M[29] + M[49];
#pragma omp atomic
  Ms[50] += Mstmp34 * Mstmp36 + Mstmp34 * M[16] + Mstmp69 * z + Mstmp7 * Mstmp70 +
            Mstmp70 * M[7] + Mstmp84 * M[2] + y * M[30] + z * M[29] + M[50];
#pragma omp atomic
  Ms[51] += Mstmp33 * Mstmp41 + Mstmp34 * Mstmp39 + Mstmp34 * M[17] + Mstmp41 * M[15] +
            Mstmp70 * Mstmp8 + Mstmp70 * M[8] + Mstmp71 * z + Mstmp82 * M[1] + y * M[31] +
            z * M[30] + M[51];
#pragma omp atomic
  Ms[52] += Mstmp34 * Mstmp42 + Mstmp34 * M[18] + Mstmp35 * Mstmp41 + Mstmp41 * M[16] +
            Mstmp5 * Mstmp74 + Mstmp72 * z + Mstmp74 * M[6] + Mstmp82 * M[2] + y * M[32] +
            z * M[31] + M[52];
#pragma omp atomic
  Ms[53] += Mstmp38 * Mstmp41 + Mstmp41 * M[17] + Mstmp6 * Mstmp74 + Mstmp73 * z +
            Mstmp74 * M[7] + Mstmp85 * M[1] + y * M[33] + z * M[32] + M[53];
#pragma omp atomic
  Ms[54] += Mstmp41 * M[18] + Mstmp74 * M[8] + Mstmp85 * M[2] + z * M[33] + M[54];
}

void field_m1_M2L_5(double x, double y, double z, double* M, double* L) {
  double R = sqrt(x * x + y * y + z * z);
  double D[55];
  double Dtmp0  = 1.0 * pow(R, -3.0);
  double Dtmp1  = -Dtmp0;
  double Dtmp2  = (x * x);
  double Dtmp3  = pow(R, -5.0);
  double Dtmp4  = 3.0 * Dtmp3;
  double Dtmp5  = Dtmp4 * x;
  double Dtmp6  = (y * y);
  double Dtmp7  = y * z;
  double Dtmp8  = 9.0 * Dtmp3;
  double Dtmp9  = -Dtmp8;
  double Dtmp10 = pow(R, -7.0);
  double Dtmp11 = 15.0 * Dtmp10;
  double Dtmp12 = Dtmp11 * Dtmp2;
  double Dtmp13 = -Dtmp4;
  double Dtmp14 = Dtmp12 + Dtmp13;
  double Dtmp15 = Dtmp11 * Dtmp6;
  double Dtmp16 = Dtmp13 + Dtmp15;
  double Dtmp17 = 1.0 * x;
  double Dtmp18 = Dtmp7 * x;
  double Dtmp19 = (x * x * x * x);
  double Dtmp20 = pow(R, -9.0);
  double Dtmp21 = 105.0 * Dtmp20;
  double Dtmp22 = 90.0 * Dtmp10;
  double Dtmp23 = 45.0 * Dtmp10;
  double Dtmp24 = -Dtmp23;
  double Dtmp25 = Dtmp2 * Dtmp21;
  double Dtmp26 = x * (Dtmp24 + Dtmp25);
  double Dtmp27 = -Dtmp11;
  double Dtmp28 = Dtmp21 * Dtmp6;
  double Dtmp29 = Dtmp24 + Dtmp28;
  double Dtmp30 = (y * y * y * y);
  double Dtmp31 = 225.0 * Dtmp10;
  double Dtmp32 = 945.0 * pow(R, -11.0);
  double Dtmp33 = Dtmp19 * Dtmp32;
  double Dtmp34 = Dtmp2 * Dtmp20;
  double Dtmp35 = Dtmp23 + Dtmp33 - 630.0 * Dtmp34;
  double Dtmp36 = -Dtmp25;
  double Dtmp37 = 315.0 * Dtmp20;
  double Dtmp38 = Dtmp2 * Dtmp32;
  double Dtmp39 = Dtmp38 * Dtmp6;
  double Dtmp40 = Dtmp23 + Dtmp39;
  double Dtmp41 = -Dtmp37;
  double Dtmp42 = -Dtmp28;
  double Dtmp43 = Dtmp30 * Dtmp32;
  double Dtmp44 = Dtmp20 * Dtmp6;
  double Dtmp45 = Dtmp23 + Dtmp43 - 630.0 * Dtmp44;
  D[0]          = -Dtmp0 * x;
  D[1]          = -Dtmp0 * y;
  D[2]          = -Dtmp0 * z;
  D[3]          = Dtmp1 + Dtmp2 * Dtmp4;
  D[4]          = Dtmp5 * y;
  D[5]          = Dtmp5 * z;
  D[6]          = Dtmp1 + Dtmp4 * Dtmp6;
  D[7]          = Dtmp4 * Dtmp7;
  D[8]          = -D[3] - D[6];
  D[9]          = -x * (Dtmp12 + Dtmp9);
  D[10]         = -Dtmp14 * y;
  D[11]         = -Dtmp14 * z;
  D[12]         = -Dtmp16 * Dtmp17;
  D[13]         = -Dtmp11 * Dtmp18;
  D[14]         = -D[9] - D[12];
  D[15]         = -y * (Dtmp15 + Dtmp9);
  D[16]         = -Dtmp16 * z;
  D[17]         = -D[10] - D[15];
  D[18]         = -D[11] - D[16];
  D[19]         = Dtmp19 * Dtmp21 - Dtmp2 * Dtmp22 + Dtmp8;
  D[20]         = Dtmp26 * y;
  D[21]         = Dtmp26 * z;
  D[22]         = -Dtmp12 - Dtmp15 + Dtmp25 * Dtmp6 + Dtmp4;
  D[23]         = Dtmp7 * (Dtmp25 + Dtmp27);
  D[24]         = -D[19] - D[22];
  D[25]         = Dtmp17 * Dtmp29 * y;
  D[26]         = Dtmp17 * z * (Dtmp27 + Dtmp28);
  D[27]         = -D[20] - D[25];
  D[28]         = -D[21] - D[26];
  D[29]         = Dtmp21 * Dtmp30 - Dtmp22 * Dtmp6 + Dtmp8;
  D[30]         = Dtmp29 * Dtmp7;
  D[31]         = -D[22] - D[29];
  D[32]         = -D[23] - D[30];
  D[33]         = -D[24] - D[31];
  D[34]         = -x * (Dtmp31 + Dtmp33 - 1050.0 * Dtmp34);
  D[35]         = -Dtmp35 * y;
  D[36]         = -Dtmp35 * z;
  D[37]         = -x * (Dtmp36 - Dtmp37 * Dtmp6 + Dtmp40);
  D[38]         = -Dtmp18 * (Dtmp38 + Dtmp41);
  D[39]         = -D[34] - D[37];
  D[40]         = -y * (-Dtmp2 * Dtmp37 + Dtmp40 + Dtmp42);
  D[41]         = -z * (Dtmp11 + Dtmp36 + Dtmp39 + Dtmp42);
  D[42]         = -D[35] - D[40];
  D[43]         = -D[36] - D[41];
  D[44]         = -Dtmp17 * Dtmp45;
  D[45]         = -Dtmp17 * Dtmp7 * (Dtmp32 * Dtmp6 + Dtmp41);
  D[46]         = -D[37] - D[44];
  D[47]         = -D[38] - D[45];
  D[48]         = -D[39] - D[46];
  D[49]         = -y * (Dtmp31 + Dtmp43 - 1050.0 * Dtmp44);
  D[50]         = -Dtmp45 * z;
  D[51]         = -D[40] - D[49];
  D[52]         = -D[41] - D[50];
  D[53]         = -D[42] - D[51];
  D[54]         = -D[43] - D[52];
#pragma omp atomic
  L[0] += D[0] * M[0] + D[1] * M[1] + D[2] * M[2] + D[3] * M[3] + D[4] * M[4] +
          D[5] * M[5] + D[6] * M[6] + D[7] * M[7] + D[8] * M[8] + D[9] * M[9] +
          D[10] * M[10] + D[11] * M[11] + D[12] * M[12] + D[13] * M[13] + D[14] * M[14] +
          D[15] * M[15] + D[16] * M[16] + D[17] * M[17] + D[18] * M[18] + D[19] * M[19] +
          D[20] * M[20] + D[21] * M[21] + D[22] * M[22] + D[23] * M[23] + D[24] * M[24] +
          D[25] * M[25] + D[26] * M[26] + D[27] * M[27] + D[28] * M[28] + D[29] * M[29] +
          D[30] * M[30] + D[31] * M[31] + D[32] * M[32] + D[33] * M[33] + D[34] * M[34] +
          D[35] * M[35] + D[36] * M[36] + D[37] * M[37] + D[38] * M[38] + D[39] * M[39] +
          D[40] * M[40] + D[41] * M[41] + D[42] * M[42] + D[43] * M[43] + D[44] * M[44] +
          D[45] * M[45] + D[46] * M[46] + D[47] * M[47] + D[48] * M[48] + D[49] * M[49] +
          D[50] * M[50] + D[51] * M[51] + D[52] * M[52] + D[53] * M[53] + D[54] * M[54];
#pragma omp atomic
  L[1] += D[3] * M[0] + D[4] * M[1] + D[5] * M[2] + D[9] * M[3] + D[10] * M[4] +
          D[11] * M[5] + D[12] * M[6] + D[13] * M[7] + D[14] * M[8] + D[19] * M[9] +
          D[20] * M[10] + D[21] * M[11] + D[22] * M[12] + D[23] * M[13] + D[24] * M[14] +
          D[25] * M[15] + D[26] * M[16] + D[27] * M[17] + D[28] * M[18] + D[34] * M[19] +
          D[35] * M[20] + D[36] * M[21] + D[37] * M[22] + D[38] * M[23] + D[39] * M[24] +
          D[40] * M[25] + D[41] * M[26] + D[42] * M[27] + D[43] * M[28] + D[44] * M[29] +
          D[45] * M[30] + D[46] * M[31] + D[47] * M[32] + D[48] * M[33];
#pragma omp atomic
  L[2] += D[4] * M[0] + D[6] * M[1] + D[7] * M[2] + D[10] * M[3] + D[12] * M[4] +
          D[13] * M[5] + D[15] * M[6] + D[16] * M[7] + D[17] * M[8] + D[20] * M[9] +
          D[22] * M[10] + D[23] * M[11] + D[25] * M[12] + D[26] * M[13] + D[27] * M[14] +
          D[29] * M[15] + D[30] * M[16] + D[31] * M[17] + D[32] * M[18] + D[35] * M[19] +
          D[37] * M[20] + D[38] * M[21] + D[40] * M[22] + D[41] * M[23] + D[42] * M[24] +
          D[44] * M[25] + D[45] * M[26] + D[46] * M[27] + D[47] * M[28] + D[49] * M[29] +
          D[50] * M[30] + D[51] * M[31] + D[52] * M[32] + D[53] * M[33];
#pragma omp atomic
  L[3] += D[5] * M[0] + D[7] * M[1] + D[8] * M[2] + D[11] * M[3] + D[13] * M[4] +
          D[14] * M[5] + D[16] * M[6] + D[17] * M[7] + D[18] * M[8] + D[21] * M[9] +
          D[23] * M[10] + D[24] * M[11] + D[26] * M[12] + D[27] * M[13] + D[28] * M[14] +
          D[30] * M[15] + D[31] * M[16] + D[32] * M[17] + D[33] * M[18] + D[36] * M[19] +
          D[38] * M[20] + D[39] * M[21] + D[41] * M[22] + D[42] * M[23] + D[43] * M[24] +
          D[45] * M[25] + D[46] * M[26] + D[47] * M[27] + D[48] * M[28] + D[50] * M[29] +
          D[51] * M[30] + D[52] * M[31] + D[53] * M[32] + D[54] * M[33];
#pragma omp atomic
  L[4] += D[9] * M[0] + D[10] * M[1] + D[11] * M[2] + D[19] * M[3] + D[20] * M[4] +
          D[21] * M[5] + D[22] * M[6] + D[23] * M[7] + D[24] * M[8] + D[34] * M[9] +
          D[35] * M[10] + D[36] * M[11] + D[37] * M[12] + D[38] * M[13] + D[39] * M[14] +
          D[40] * M[15] + D[41] * M[16] + D[42] * M[17] + D[43] * M[18];
#pragma omp atomic
  L[5] += D[10] * M[0] + D[12] * M[1] + D[13] * M[2] + D[20] * M[3] + D[22] * M[4] +
          D[23] * M[5] + D[25] * M[6] + D[26] * M[7] + D[27] * M[8] + D[35] * M[9] +
          D[37] * M[10] + D[38] * M[11] + D[40] * M[12] + D[41] * M[13] + D[42] * M[14] +
          D[44] * M[15] + D[45] * M[16] + D[46] * M[17] + D[47] * M[18];
#pragma omp atomic
  L[6] += D[11] * M[0] + D[13] * M[1] + D[14] * M[2] + D[21] * M[3] + D[23] * M[4] +
          D[24] * M[5] + D[26] * M[6] + D[27] * M[7] + D[28] * M[8] + D[36] * M[9] +
          D[38] * M[10] + D[39] * M[11] + D[41] * M[12] + D[42] * M[13] + D[43] * M[14] +
          D[45] * M[15] + D[46] * M[16] + D[47] * M[17] + D[48] * M[18];
#pragma omp atomic
  L[7] += D[12] * M[0] + D[15] * M[1] + D[16] * M[2] + D[22] * M[3] + D[25] * M[4] +
          D[26] * M[5] + D[29] * M[6] + D[30] * M[7] + D[31] * M[8] + D[37] * M[9] +
          D[40] * M[10] + D[41] * M[11] + D[44] * M[12] + D[45] * M[13] + D[46] * M[14] +
          D[49] * M[15] + D[50] * M[16] + D[51] * M[17] + D[52] * M[18];
#pragma omp atomic
  L[8] += D[13] * M[0] + D[16] * M[1] + D[17] * M[2] + D[23] * M[3] + D[26] * M[4] +
          D[27] * M[5] + D[30] * M[6] + D[31] * M[7] + D[32] * M[8] + D[38] * M[9] +
          D[41] * M[10] + D[42] * M[11] + D[45] * M[12] + D[46] * M[13] + D[47] * M[14] +
          D[50] * M[15] + D[51] * M[16] + D[52] * M[17] + D[53] * M[18];
#pragma omp atomic
  L[9] += D[14] * M[0] + D[17] * M[1] + D[18] * M[2] + D[24] * M[3] + D[27] * M[4] +
          D[28] * M[5] + D[31] * M[6] + D[32] * M[7] + D[33] * M[8] + D[39] * M[9] +
          D[42] * M[10] + D[43] * M[11] + D[46] * M[12] + D[47] * M[13] + D[48] * M[14] +
          D[51] * M[15] + D[52] * M[16] + D[53] * M[17] + D[54] * M[18];
#pragma omp atomic
  L[10] += D[19] * M[0] + D[20] * M[1] + D[21] * M[2] + D[34] * M[3] + D[35] * M[4] +
           D[36] * M[5] + D[37] * M[6] + D[38] * M[7] + D[39] * M[8];
#pragma omp atomic
  L[11] += D[20] * M[0] + D[22] * M[1] + D[23] * M[2] + D[35] * M[3] + D[37] * M[4] +
           D[38] * M[5] + D[40] * M[6] + D[41] * M[7] + D[42] * M[8];
#pragma omp atomic
  L[12] += D[21] * M[0] + D[23] * M[1] + D[24] * M[2] + D[36] * M[3] + D[38] * M[4] +
           D[39] * M[5] + D[41] * M[6] + D[42] * M[7] + D[43] * M[8];
#pragma omp atomic
  L[13] += D[22] * M[0] + D[25] * M[1] + D[26] * M[2] + D[37] * M[3] + D[40] * M[4] +
           D[41] * M[5] + D[44] * M[6] + D[45] * M[7] + D[46] * M[8];
#pragma omp atomic
  L[14] += D[23] * M[0] + D[26] * M[1] + D[27] * M[2] + D[38] * M[3] + D[41] * M[4] +
           D[42] * M[5] + D[45] * M[6] + D[46] * M[7] + D[47] * M[8];
#pragma omp atomic
  L[15] += D[24] * M[0] + D[27] * M[1] + D[28] * M[2] + D[39] * M[3] + D[42] * M[4] +
           D[43] * M[5] + D[46] * M[6] + D[47] * M[7] + D[48] * M[8];
#pragma omp atomic
  L[16] += D[25] * M[0] + D[29] * M[1] + D[30] * M[2] + D[40] * M[3] + D[44] * M[4] +
           D[45] * M[5] + D[49] * M[6] + D[50] * M[7] + D[51] * M[8];
#pragma omp atomic
  L[17] += D[26] * M[0] + D[30] * M[1] + D[31] * M[2] + D[41] * M[3] + D[45] * M[4] +
           D[46] * M[5] + D[50] * M[6] + D[51] * M[7] + D[52] * M[8];
#pragma omp atomic
  L[18] += D[27] * M[0] + D[31] * M[1] + D[32] * M[2] + D[42] * M[3] + D[46] * M[4] +
           D[47] * M[5] + D[51] * M[6] + D[52] * M[7] + D[53] * M[8];
#pragma omp atomic
  L[19] += D[28] * M[0] + D[32] * M[1] + D[33] * M[2] + D[43] * M[3] + D[47] * M[4] +
           D[48] * M[5] + D[52] * M[6] + D[53] * M[7] + D[54] * M[8];
#pragma omp atomic
  L[20] += D[34] * M[0] + D[35] * M[1] + D[36] * M[2];
#pragma omp atomic
  L[21] += D[35] * M[0] + D[37] * M[1] + D[38] * M[2];
#pragma omp atomic
  L[22] += D[36] * M[0] + D[38] * M[1] + D[39] * M[2];
#pragma omp atomic
  L[23] += D[37] * M[0] + D[40] * M[1] + D[41] * M[2];
#pragma omp atomic
  L[24] += D[38] * M[0] + D[41] * M[1] + D[42] * M[2];
#pragma omp atomic
  L[25] += D[39] * M[0] + D[42] * M[1] + D[43] * M[2];
#pragma omp atomic
  L[26] += D[40] * M[0] + D[44] * M[1] + D[45] * M[2];
#pragma omp atomic
  L[27] += D[41] * M[0] + D[45] * M[1] + D[46] * M[2];
#pragma omp atomic
  L[28] += D[42] * M[0] + D[46] * M[1] + D[47] * M[2];
#pragma omp atomic
  L[29] += D[43] * M[0] + D[47] * M[1] + D[48] * M[2];
#pragma omp atomic
  L[30] += D[44] * M[0] + D[49] * M[1] + D[50] * M[2];
#pragma omp atomic
  L[31] += D[45] * M[0] + D[50] * M[1] + D[51] * M[2];
#pragma omp atomic
  L[32] += D[46] * M[0] + D[51] * M[1] + D[52] * M[2];
#pragma omp atomic
  L[33] += D[47] * M[0] + D[52] * M[1] + D[53] * M[2];
#pragma omp atomic
  L[34] += D[48] * M[0] + D[53] * M[1] + D[54] * M[2];
}

void field_m1_L2L_5(double x, double y, double z, double* L, double* Ls) {
  double Lstmp0  = y * L[5];
  double Lstmp1  = z * L[6];
  double Lstmp2  = z * L[8];
  double Lstmp3  = z * L[14];
  double Lstmp4  = Lstmp3 * y;
  double Lstmp5  = (x * x);
  double Lstmp6  = (1.0 / 2.0) * Lstmp5;
  double Lstmp7  = (1.0 / 6.0) * (x * x * x);
  double Lstmp8  = (y * y);
  double Lstmp9  = (1.0 / 2.0) * Lstmp8;
  double Lstmp10 = (1.0 / 6.0) * (y * y * y);
  double Lstmp11 = (z * z);
  double Lstmp12 = (1.0 / 2.0) * Lstmp11;
  double Lstmp13 = (1.0 / 6.0) * (z * z * z);
  double Lstmp14 = x * L[13];
  double Lstmp15 = x * L[26];
  double Lstmp16 = x * L[15];
  double Lstmp17 = x * L[29];
  double Lstmp18 = y * L[11];
  double Lstmp19 = z * L[12];
  double Lstmp20 = y * L[21];
  double Lstmp21 = z * L[22];
  double Lstmp22 = y * L[18];
  double Lstmp23 = y * L[33];
  double Lstmp24 = z * L[17];
  double Lstmp25 = z * L[31];
  double Lstmp26 = y * L[28];
  double Lstmp27 = Lstmp26 * x;
  double Lstmp28 = z * L[27];
  double Lstmp29 = Lstmp28 * x;
  double Lstmp30 = z * L[24];
  double Lstmp31 = Lstmp30 * y;
  double Lstmp32 = (1.0 / 4.0) * Lstmp5;
  double Lstmp33 = x * L[23];
  double Lstmp34 = x * L[25];
  double Lstmp35 = y * L[13];
  double Lstmp36 = Lstmp28 * y;
  double Lstmp37 = x * L[28];
  double Lstmp38 = y * L[23];
  double Lstmp39 = y * L[32];
  double Lstmp40 = y * L[14];
  double Lstmp41 = z * L[15];
  double Lstmp42 = z * L[18];
  double Lstmp43 = z * L[28];
  double Lstmp44 = Lstmp43 * y;
  double Lstmp45 = x * L[27];
  double Lstmp46 = y * L[24];
  double Lstmp47 = z * L[25];
  double Lstmp48 = z * L[32];
  double Lstmp49 = y * L[26];
  double Lstmp50 = y * L[27];
  double Lstmp51 = z * L[29];
  double Lstmp52 = z * L[33];
#pragma omp atomic
  Ls[0] += Lstmp0 * x + Lstmp1 * x + Lstmp10 * Lstmp15 + Lstmp10 * Lstmp25 +
           Lstmp10 * L[16] + Lstmp11 * Lstmp32 * L[25] +
           (1.0 / 4.0) * Lstmp11 * Lstmp8 * L[32] + Lstmp12 * Lstmp16 +
           Lstmp12 * Lstmp22 + Lstmp12 * Lstmp27 + Lstmp12 * L[9] + Lstmp13 * Lstmp17 +
           Lstmp13 * Lstmp23 + Lstmp13 * L[19] + Lstmp14 * Lstmp9 + Lstmp18 * Lstmp6 +
           Lstmp19 * Lstmp6 + Lstmp2 * y + Lstmp20 * Lstmp7 + Lstmp21 * Lstmp7 +
           Lstmp24 * Lstmp9 + Lstmp29 * Lstmp9 + Lstmp31 * Lstmp6 +
           Lstmp32 * Lstmp8 * L[23] + Lstmp4 * x + Lstmp6 * L[4] + Lstmp7 * L[10] +
           Lstmp9 * L[7] + (1.0 / 24.0) * (x * x * x * x) * L[20] + x * L[1] +
           (1.0 / 24.0) * (y * y * y * y) * L[30] + y * L[2] +
           (1.0 / 24.0) * (z * z * z * z) * L[34] + z * L[3] + L[0];
#pragma omp atomic
  Ls[1] += Lstmp0 + Lstmp1 + Lstmp10 * L[26] + Lstmp12 * Lstmp26 + Lstmp12 * Lstmp34 +
           Lstmp12 * L[15] + Lstmp13 * L[29] + Lstmp18 * x + Lstmp19 * x +
           Lstmp20 * Lstmp6 + Lstmp21 * Lstmp6 + Lstmp28 * Lstmp9 + Lstmp31 * x +
           Lstmp33 * Lstmp9 + Lstmp4 + Lstmp6 * L[10] + Lstmp7 * L[20] + Lstmp9 * L[13] +
           x * L[4] + L[1];
#pragma omp atomic
  Ls[2] += Lstmp10 * L[30] + Lstmp12 * Lstmp37 + Lstmp12 * Lstmp39 + Lstmp12 * L[18] +
           Lstmp13 * L[33] + Lstmp15 * Lstmp9 + Lstmp2 + Lstmp24 * y + Lstmp25 * Lstmp9 +
           Lstmp3 * x + Lstmp30 * Lstmp6 + Lstmp35 * x + Lstmp36 * x + Lstmp38 * Lstmp6 +
           Lstmp6 * L[11] + Lstmp7 * L[21] + Lstmp9 * L[16] + x * L[5] + y * L[7] + L[2];
#pragma omp atomic
  Ls[3] += Lstmp10 * L[31] + Lstmp12 * Lstmp17 + Lstmp12 * Lstmp23 + Lstmp12 * L[19] +
           Lstmp13 * L[34] + Lstmp40 * x + Lstmp41 * x + Lstmp42 * y + Lstmp44 * x +
           Lstmp45 * Lstmp9 + Lstmp46 * Lstmp6 + Lstmp47 * Lstmp6 + Lstmp48 * Lstmp9 +
           Lstmp6 * L[12] + Lstmp7 * L[22] + Lstmp9 * L[17] + x * L[6] + y * L[8] +
           z * L[9] + L[3];
#pragma omp atomic
  Ls[4] += Lstmp12 * L[25] + Lstmp18 + Lstmp19 + Lstmp20 * x + Lstmp21 * x + Lstmp31 +
           Lstmp6 * L[20] + Lstmp9 * L[23] + x * L[10] + L[4];
#pragma omp atomic
  Ls[5] += Lstmp12 * L[28] + Lstmp3 + Lstmp30 * x + Lstmp35 + Lstmp36 + Lstmp38 * x +
           Lstmp6 * L[21] + Lstmp9 * L[26] + x * L[11] + L[5];
#pragma omp atomic
  Ls[6] += Lstmp12 * L[29] + Lstmp40 + Lstmp41 + Lstmp44 + Lstmp46 * x + Lstmp47 * x +
           Lstmp6 * L[22] + Lstmp9 * L[27] + x * L[12] + L[6];
#pragma omp atomic
  Ls[7] += Lstmp12 * L[32] + Lstmp14 + Lstmp24 + Lstmp25 * y + Lstmp29 + Lstmp49 * x +
           Lstmp6 * L[23] + Lstmp9 * L[30] + y * L[16] + L[7];
#pragma omp atomic
  Ls[8] += Lstmp12 * L[33] + Lstmp42 + Lstmp43 * x + Lstmp48 * y + Lstmp50 * x +
           Lstmp6 * L[24] + Lstmp9 * L[31] + x * L[14] + y * L[17] + L[8];
#pragma omp atomic
  Ls[9] += Lstmp12 * L[34] + Lstmp16 + Lstmp22 + Lstmp27 + Lstmp51 * x + Lstmp52 * y +
           Lstmp6 * L[25] + Lstmp9 * L[32] + z * L[19] + L[9];
#pragma omp atomic
  Ls[10] += Lstmp20 + Lstmp21 + x * L[20] + L[10];
#pragma omp atomic
  Ls[11] += Lstmp30 + Lstmp38 + x * L[21] + L[11];
#pragma omp atomic
  Ls[12] += Lstmp46 + Lstmp47 + x * L[22] + L[12];
#pragma omp atomic
  Ls[13] += Lstmp28 + Lstmp33 + Lstmp49 + L[13];
#pragma omp atomic
  Ls[14] += Lstmp43 + Lstmp50 + x * L[24] + L[14];
#pragma omp atomic
  Ls[15] += Lstmp26 + Lstmp34 + Lstmp51 + L[15];
#pragma omp atomic
  Ls[16] += Lstmp15 + Lstmp25 + y * L[30] + L[16];
#pragma omp atomic
  Ls[17] += Lstmp45 + Lstmp48 + y * L[31] + L[17];
#pragma omp atomic
  Ls[18] += Lstmp37 + Lstmp39 + Lstmp52 + L[18];
#pragma omp atomic
  Ls[19] += Lstmp17 + Lstmp23 + z * L[34] + L[19];
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

void field_m1_L2P_5(double x, double y, double z, double* L, double* F) {
  double Ftmp0  = x * y;
  double Ftmp1  = x * z;
  double Ftmp2  = y * z;
  double Ftmp3  = Ftmp0 * z;
  double Ftmp4  = (1.0 / 2.0) * (x * x);
  double Ftmp5  = (1.0 / 6.0) * (x * x * x);
  double Ftmp6  = (1.0 / 2.0) * (y * y);
  double Ftmp7  = (1.0 / 6.0) * (y * y * y);
  double Ftmp8  = (1.0 / 2.0) * (z * z);
  double Ftmp9  = (1.0 / 6.0) * (z * z * z);
  double Ftmp10 = Ftmp6 * x;
  double Ftmp11 = Ftmp8 * x;
  double Ftmp12 = Ftmp4 * y;
  double Ftmp13 = Ftmp4 * z;
  double Ftmp14 = Ftmp8 * y;
  double Ftmp15 = Ftmp6 * z;
#pragma omp atomic
  F[0] += -Ftmp0 * L[11] - Ftmp1 * L[12] - Ftmp10 * L[23] - Ftmp11 * L[25] -
          Ftmp12 * L[21] - Ftmp13 * L[22] - Ftmp14 * L[28] - Ftmp15 * L[27] -
          Ftmp2 * L[14] - Ftmp3 * L[24] - Ftmp4 * L[10] - Ftmp5 * L[20] - Ftmp6 * L[13] -
          Ftmp7 * L[26] - Ftmp8 * L[15] - Ftmp9 * L[29] - x * L[4] - y * L[5] - z * L[6] -
          L[1];
#pragma omp atomic
  F[1] += -Ftmp0 * L[13] - Ftmp1 * L[14] - Ftmp10 * L[26] - Ftmp11 * L[28] -
          Ftmp12 * L[23] - Ftmp13 * L[24] - Ftmp14 * L[32] - Ftmp15 * L[31] -
          Ftmp2 * L[17] - Ftmp3 * L[27] - Ftmp4 * L[11] - Ftmp5 * L[21] - Ftmp6 * L[16] -
          Ftmp7 * L[30] - Ftmp8 * L[18] - Ftmp9 * L[33] - x * L[5] - y * L[7] - z * L[8] -
          L[2];
#pragma omp atomic
  F[2] += -Ftmp0 * L[14] - Ftmp1 * L[15] - Ftmp10 * L[27] - Ftmp11 * L[29] -
          Ftmp12 * L[24] - Ftmp13 * L[25] - Ftmp14 * L[33] - Ftmp15 * L[32] -
          Ftmp2 * L[18] - Ftmp3 * L[28] - Ftmp4 * L[12] - Ftmp5 * L[22] - Ftmp6 * L[17] -
          Ftmp7 * L[31] - Ftmp8 * L[19] - Ftmp9 * L[34] - x * L[6] - y * L[8] - z * L[9] -
          L[3];
}

void field_m1_M2P_5(double x, double y, double z, double* M, double* F) {
  double R       = sqrt(x * x + y * y + z * z);
  double Ftmp0   = 1.0 * pow(R, -3.0);
  double Ftmp1   = pow(R, -5.0);
  double Ftmp2   = 3.0 * Ftmp1;
  double Ftmp3   = Ftmp2 * M[4];
  double Ftmp4   = Ftmp2 * z;
  double Ftmp5   = y * z;
  double Ftmp6   = pow(R, -7.0);
  double Ftmp7   = 15.0 * Ftmp6;
  double Ftmp8   = Ftmp7 * M[13];
  double Ftmp9   = x * y;
  double Ftmp10  = Ftmp2 * M[1];
  double Ftmp11  = Ftmp4 * M[2];
  double Ftmp12  = (x * x);
  double Ftmp13  = Ftmp2 * M[0];
  double Ftmp14  = Ftmp5 * x;
  double Ftmp15  = Ftmp14 * Ftmp7;
  double Ftmp16  = Ftmp12 * Ftmp7;
  double Ftmp17  = pow(R, -9.0);
  double Ftmp18  = 105.0 * Ftmp17;
  double Ftmp19  = Ftmp12 * Ftmp18;
  double Ftmp20  = -9.0 * Ftmp1;
  double Ftmp21  = Ftmp16 + Ftmp20;
  double Ftmp22  = -Ftmp2;
  double Ftmp23  = (y * y);
  double Ftmp24  = Ftmp23 * Ftmp7;
  double Ftmp25  = Ftmp22 + Ftmp24;
  double Ftmp26  = (z * z);
  double Ftmp27  = Ftmp26 * Ftmp7;
  double Ftmp28  = Ftmp22 + Ftmp27;
  double Ftmp29  = Ftmp25 * M[6];
  double Ftmp30  = Ftmp28 * M[8];
  double Ftmp31  = 45.0 * Ftmp6;
  double Ftmp32  = -Ftmp31;
  double Ftmp33  = Ftmp19 + Ftmp32;
  double Ftmp34  = Ftmp33 * M[20];
  double Ftmp35  = Ftmp18 * Ftmp23;
  double Ftmp36  = Ftmp32 + Ftmp35;
  double Ftmp37  = Ftmp36 * y;
  double Ftmp38  = 1.0 * M[25];
  double Ftmp39  = 3.0 * y;
  double Ftmp40  = 35.0 * Ftmp17;
  double Ftmp41  = (Ftmp26 * Ftmp40 - 5.0 * Ftmp6) * M[27];
  double Ftmp42  = Ftmp33 * M[21];
  double Ftmp43  = 1.0 * z;
  double Ftmp44  = -Ftmp7;
  double Ftmp45  = Ftmp35 + Ftmp44;
  double Ftmp46  = Ftmp45 * M[26];
  double Ftmp47  = Ftmp18 * Ftmp26;
  double Ftmp48  = Ftmp32 + Ftmp47;
  double Ftmp49  = Ftmp43 * Ftmp48;
  double Ftmp50  = 315.0 * Ftmp17;
  double Ftmp51  = -Ftmp50;
  double Ftmp52  = pow(R, -11.0);
  double Ftmp53  = 945.0 * Ftmp52;
  double Ftmp54  = Ftmp12 * Ftmp53;
  double Ftmp55  = Ftmp51 + Ftmp54;
  double Ftmp56  = Ftmp55 * M[38];
  double Ftmp57  = Ftmp33 * Ftmp9;
  double Ftmp58  = Ftmp37 * M[15];
  double Ftmp59  = Ftmp44 + Ftmp47;
  double Ftmp60  = Ftmp59 * M[17];
  double Ftmp61  = 1.0 * Ftmp9;
  double Ftmp62  = x * z;
  double Ftmp63  = Ftmp33 * Ftmp62;
  double Ftmp64  = Ftmp45 * M[16];
  double Ftmp65  = Ftmp48 * M[18];
  double Ftmp66  = Ftmp23 * Ftmp53;
  double Ftmp67  = Ftmp51 + Ftmp66;
  double Ftmp68  = Ftmp67 * y;
  double Ftmp69  = Ftmp43 * M[45];
  double Ftmp70  = -Ftmp18;
  double Ftmp71  = Ftmp26 * Ftmp52;
  double Ftmp72  = 315.0 * Ftmp71;
  double Ftmp73  = Ftmp70 + Ftmp72;
  double Ftmp74  = Ftmp73 * M[47];
  double Ftmp75  = Ftmp39 * Ftmp74;
  double Ftmp76  = Ftmp55 * x;
  double Ftmp77  = Ftmp67 * M[30];
  double Ftmp78  = -75.0 * Ftmp6;
  double Ftmp79  = 1.0 * Ftmp12;
  double Ftmp80  = Ftmp45 * M[12];
  double Ftmp81  = Ftmp59 * M[14];
  double Ftmp82  = 525.0 * Ftmp17;
  double Ftmp83  = -Ftmp82;
  double Ftmp84  = Ftmp12 * (Ftmp54 + Ftmp83);
  double Ftmp85  = Ftmp26 * Ftmp53;
  double Ftmp86  = Ftmp51 + Ftmp85;
  double Ftmp87  = Ftmp43 * Ftmp86 * M[32];
  double Ftmp88  = (-Ftmp40 + Ftmp72) * M[27];
  double Ftmp89  = Ftmp12 * Ftmp39;
  double Ftmp90  = Ftmp12 * Ftmp43;
  double Ftmp91  = (Ftmp66 + Ftmp70) * M[26];
  double Ftmp92  = Ftmp86 * M[28];
  double Ftmp93  = 4725.0 * Ftmp52;
  double Ftmp94  = -Ftmp93;
  double Ftmp95  = pow(R, -13.0);
  double Ftmp96  = 10395.0 * Ftmp95;
  double Ftmp97  = Ftmp12 * Ftmp96;
  double Ftmp98  = 2835.0 * Ftmp52;
  double Ftmp99  = -Ftmp98;
  double Ftmp100 = Ftmp23 * Ftmp96;
  double Ftmp101 = Ftmp100 + Ftmp99;
  double Ftmp102 = 3465.0 * Ftmp26 * Ftmp95;
  double Ftmp103 = (Ftmp102 - Ftmp53) * M[47];
  double Ftmp104 = 225.0 * Ftmp6;
  double Ftmp105 = (x * x * x * x);
  double Ftmp106 = Ftmp105 * Ftmp53;
  double Ftmp107 = 1050.0 * Ftmp17;
  double Ftmp108 = Ftmp104 + Ftmp106 - Ftmp107 * Ftmp12;
  double Ftmp109 = (y * y * y * y);
  double Ftmp110 = Ftmp109 * Ftmp53;
  double Ftmp111 = 630.0 * Ftmp17;
  double Ftmp112 = Ftmp110 - Ftmp111 * Ftmp23 + Ftmp31;
  double Ftmp113 = (z * z * z * z);
  double Ftmp114 = Ftmp113 * Ftmp53;
  double Ftmp115 = -Ftmp111 * Ftmp26 + Ftmp114 + Ftmp31;
  double Ftmp116 = Ftmp112 * M[29];
  double Ftmp117 = Ftmp115 * M[33];
  double Ftmp118 = 1575.0 * Ftmp17;
  double Ftmp119 = Ftmp105 * Ftmp96;
  double Ftmp120 = Ftmp12 * Ftmp52;
  double Ftmp121 = Ftmp118 + Ftmp119 - 9450.0 * Ftmp120;
  double Ftmp122 = Ftmp121 * Ftmp9;
  double Ftmp123 = Ftmp109 * Ftmp96;
  double Ftmp124 = 9450.0 * Ftmp52;
  double Ftmp125 = Ftmp118 + Ftmp123 - Ftmp124 * Ftmp23;
  double Ftmp126 = Ftmp125 * M[49];
  double Ftmp127 = Ftmp113 * Ftmp96;
  double Ftmp128 = 5670.0 * Ftmp52;
  double Ftmp129 = Ftmp127 - Ftmp128 * Ftmp26 + Ftmp50;
  double Ftmp130 = Ftmp129 * M[53];
  double Ftmp131 = Ftmp121 * Ftmp62;
  double Ftmp132 = Ftmp123 - Ftmp128 * Ftmp23 + Ftmp50;
  double Ftmp133 = Ftmp132 * M[50];
  double Ftmp134 = Ftmp118 - Ftmp124 * Ftmp26 + Ftmp127;
  double Ftmp135 = Ftmp134 * M[54];
  double Ftmp136 = 3675.0 * Ftmp17;
  double Ftmp137 = Ftmp132 * M[44];
  double Ftmp138 = Ftmp129 * M[48];
  double Ftmp139 = -Ftmp23 * Ftmp50;
  double Ftmp140 = Ftmp23 * Ftmp54;
  double Ftmp141 = -Ftmp19;
  double Ftmp142 = Ftmp141 + Ftmp31;
  double Ftmp143 = Ftmp139 + Ftmp140 + Ftmp142;
  double Ftmp144 = -Ftmp26 * Ftmp50;
  double Ftmp145 = Ftmp26 * Ftmp54;
  double Ftmp146 = Ftmp142 + Ftmp144 + Ftmp145;
  double Ftmp147 = -Ftmp47;
  double Ftmp148 = Ftmp147 + Ftmp7;
  double Ftmp149 = -Ftmp35;
  double Ftmp150 = Ftmp26 * Ftmp66;
  double Ftmp151 = Ftmp149 + Ftmp150;
  double Ftmp152 = Ftmp148 + Ftmp151;
  double Ftmp153 = Ftmp23 * Ftmp97;
  double Ftmp154 = -Ftmp23 * Ftmp98;
  double Ftmp155 = Ftmp153 + Ftmp154;
  double Ftmp156 = 945.0 * Ftmp17;
  double Ftmp157 = -Ftmp12 * Ftmp98;
  double Ftmp158 = Ftmp156 + Ftmp157;
  double Ftmp159 = Ftmp9 * (Ftmp155 + Ftmp158);
  double Ftmp160 = -Ftmp26 * Ftmp98;
  double Ftmp161 = Ftmp160 + Ftmp50;
  double Ftmp162 = -Ftmp54;
  double Ftmp163 = Ftmp26 * Ftmp97;
  double Ftmp164 = Ftmp162 + Ftmp163;
  double Ftmp165 = Ftmp9 * (Ftmp161 + Ftmp164);
  double Ftmp166 = -Ftmp66;
  double Ftmp167 = Ftmp100 * Ftmp26;
  double Ftmp168 = Ftmp166 + Ftmp167;
  double Ftmp169 = Ftmp9 * (Ftmp161 + Ftmp168);
  double Ftmp170 = Ftmp62 * (Ftmp155 + Ftmp162 + Ftmp50);
  double Ftmp171 = Ftmp62 * (Ftmp158 + Ftmp160 + Ftmp163);
  double Ftmp172 = -Ftmp85;
  double Ftmp173 = Ftmp172 + Ftmp50;
  double Ftmp174 = Ftmp154 + Ftmp167;
  double Ftmp175 = Ftmp173 + Ftmp174;
  double Ftmp176 = -Ftmp23 * Ftmp93;
  double Ftmp177 = Ftmp162 + Ftmp82;
  double Ftmp178 = -Ftmp26 * Ftmp93;
  double Ftmp179 = Ftmp172 + Ftmp18;
  double Ftmp180 = x * M[5];
  double Ftmp181 = Ftmp16 + Ftmp22;
  double Ftmp182 = Ftmp20 + Ftmp24;
  double Ftmp183 = Ftmp181 * M[3];
  double Ftmp184 = Ftmp38 * x;
  double Ftmp185 = 3.0 * x;
  double Ftmp186 = Ftmp19 + Ftmp44;
  double Ftmp187 = Ftmp186 * M[23];
  double Ftmp188 = Ftmp36 * M[30];
  double Ftmp189 = 1.0 * x;
  double Ftmp190 = Ftmp69 * x;
  double Ftmp191 = 3.0 * Ftmp62;
  double Ftmp192 = Ftmp186 * M[11];
  double Ftmp193 = Ftmp76 * M[21];
  double Ftmp194 = Ftmp186 * M[10];
  double Ftmp195 = 1.0 * Ftmp23;
  double Ftmp196 = Ftmp76 * M[20];
  double Ftmp197 = Ftmp23 * z;
  double Ftmp198 = (Ftmp54 + Ftmp70) * M[23];
  double Ftmp199 = Ftmp66 + Ftmp83;
  double Ftmp200 = Ftmp43 * Ftmp9;
  double Ftmp201 = (Ftmp97 + Ftmp99) * M[38];
  double Ftmp202 = Ftmp106 - Ftmp111 * Ftmp12 + Ftmp31;
  double Ftmp203 = Ftmp104 - Ftmp107 * Ftmp23 + Ftmp110;
  double Ftmp204 = Ftmp202 * M[19];
  double Ftmp205 = Ftmp119 - 5670.0 * Ftmp120 + Ftmp50;
  double Ftmp206 = Ftmp205 * M[36];
  double Ftmp207 = Ftmp205 * M[35];
  double Ftmp208 = Ftmp140 + Ftmp149;
  double Ftmp209 = -Ftmp12 * Ftmp50 + Ftmp31;
  double Ftmp210 = Ftmp208 + Ftmp209;
  double Ftmp211 = Ftmp141 + Ftmp145 + Ftmp148;
  double Ftmp212 = Ftmp144 + Ftmp151 + Ftmp31;
  double Ftmp213 = Ftmp153 + Ftmp166;
  double Ftmp214 = Ftmp5 * (Ftmp157 + Ftmp213 + Ftmp50);
  double Ftmp215 = Ftmp5 * (Ftmp157 + Ftmp163 + Ftmp173);
  double Ftmp216 = Ftmp5 * (Ftmp156 + Ftmp160 + Ftmp174);
  double Ftmp217 = -Ftmp12 * Ftmp93 + Ftmp82;
  double Ftmp218 = y * M[7];
  double Ftmp219 = Ftmp20 + Ftmp27;
  double Ftmp220 = Ftmp189 * M[28];
  double Ftmp221 = 1.0 * y * M[32];
  double Ftmp222 = Ftmp61 * M[45];
  double Ftmp223 = Ftmp43 * x;
  double Ftmp224 = Ftmp26 * y;
  double Ftmp225 = Ftmp26 * (Ftmp83 + Ftmp85);
  double Ftmp226 = Ftmp104 - Ftmp107 * Ftmp26 + Ftmp114;
  double Ftmp227 = Ftmp141 + Ftmp208 + Ftmp7;
  double Ftmp228 = Ftmp145 + Ftmp147 + Ftmp209;
  double Ftmp229 = Ftmp139 + Ftmp147 + Ftmp150 + Ftmp31;
#pragma omp atomic
  F[0] += Ftmp0 * M[0] - Ftmp10 * Ftmp9 - Ftmp101 * Ftmp12 * Ftmp69 * y -
          Ftmp103 * Ftmp89 * z + Ftmp108 * x * M[19] + Ftmp108 * M[34] - Ftmp11 * x +
          Ftmp112 * M[44] + Ftmp115 * M[48] + Ftmp116 * x + Ftmp117 * x -
          Ftmp12 * Ftmp13 + Ftmp12 * Ftmp38 * Ftmp68 -
          Ftmp12 * Ftmp5 * (Ftmp94 + Ftmp97) * M[38] - Ftmp12 * (Ftmp19 + Ftmp78) * M[9] -
          Ftmp12 * (Ftmp119 - 13230.0 * Ftmp120 + Ftmp136) * M[34] -
          Ftmp12 * (Ftmp153 + Ftmp176 + Ftmp177) * M[37] -
          Ftmp12 * (Ftmp163 + Ftmp177 + Ftmp178) * M[39] - Ftmp122 * M[35] -
          Ftmp126 * Ftmp9 - Ftmp130 * Ftmp61 - Ftmp131 * M[36] - Ftmp133 * Ftmp62 -
          Ftmp135 * Ftmp62 - Ftmp137 * Ftmp79 - Ftmp138 * Ftmp79 + Ftmp14 * Ftmp77 +
          Ftmp143 * x * M[22] + Ftmp143 * M[37] + Ftmp146 * x * M[24] + Ftmp146 * M[39] +
          Ftmp15 * M[7] + Ftmp152 * x * M[31] + Ftmp152 * M[46] - Ftmp159 * M[40] +
          Ftmp16 * y * M[4] + Ftmp16 * z * M[5] - Ftmp165 * M[42] - Ftmp169 * M[51] -
          Ftmp170 * M[41] - Ftmp171 * M[43] - Ftmp175 * Ftmp62 * M[52] -
          Ftmp19 * Ftmp5 * M[13] + Ftmp21 * x * M[3] + Ftmp21 * M[9] + Ftmp25 * M[12] +
          Ftmp28 * M[14] + Ftmp29 * x - Ftmp3 * y + Ftmp30 * x - Ftmp34 * y -
          Ftmp37 * Ftmp38 - Ftmp39 * Ftmp41 - Ftmp4 * M[5] - Ftmp42 * z -
          Ftmp43 * Ftmp46 - Ftmp49 * M[28] + Ftmp5 * Ftmp56 + Ftmp5 * Ftmp76 * M[23] +
          Ftmp5 * Ftmp8 - Ftmp57 * M[10] - Ftmp58 * x - Ftmp60 * Ftmp61 -
          Ftmp62 * Ftmp64 - Ftmp62 * Ftmp65 - Ftmp63 * M[11] + Ftmp68 * Ftmp69 +
          Ftmp75 * z - Ftmp79 * Ftmp80 - Ftmp79 * Ftmp81 -
          Ftmp79 * (Ftmp168 + Ftmp179) * M[46] + Ftmp84 * y * M[20] + Ftmp84 * z * M[21] +
          Ftmp87 * Ftmp9 + Ftmp88 * Ftmp89 + Ftmp90 * Ftmp91 + Ftmp90 * Ftmp92;
#pragma omp atomic
  F[1] += Ftmp0 * M[1] - Ftmp10 * Ftmp23 - Ftmp103 * Ftmp191 * Ftmp23 - Ftmp11 * y +
          Ftmp115 * M[53] + Ftmp117 * y - Ftmp122 * M[34] - Ftmp125 * Ftmp5 * M[50] -
          Ftmp125 * Ftmp61 * M[44] - Ftmp13 * Ftmp9 - Ftmp130 * Ftmp195 -
          Ftmp135 * Ftmp5 - Ftmp138 * Ftmp61 - Ftmp159 * M[37] - Ftmp165 * M[39] -
          1.0 * Ftmp169 * M[46] + Ftmp180 * Ftmp5 * Ftmp7 + Ftmp181 * M[10] +
          Ftmp182 * y * M[6] + Ftmp182 * M[15] + Ftmp183 * y +
          Ftmp184 * Ftmp199 * Ftmp23 - Ftmp184 * Ftmp36 + Ftmp185 * Ftmp23 * Ftmp88 -
          Ftmp185 * Ftmp41 - Ftmp187 * z - Ftmp188 * z - Ftmp189 * Ftmp37 * M[12] -
          Ftmp190 * Ftmp23 * (Ftmp100 + Ftmp94) + Ftmp190 * Ftmp67 + Ftmp191 * Ftmp74 -
          Ftmp192 * Ftmp5 + Ftmp193 * Ftmp5 - Ftmp194 * Ftmp23 - Ftmp195 * Ftmp60 +
          Ftmp196 * Ftmp23 + Ftmp197 * Ftmp198 + Ftmp197 * Ftmp199 * M[30] +
          Ftmp200 * Ftmp67 * M[26] + Ftmp200 * Ftmp92 - Ftmp201 * Ftmp23 * Ftmp62 +
          Ftmp202 * M[35] + Ftmp203 * y * M[29] + Ftmp203 * M[49] + Ftmp204 * y -
          Ftmp206 * Ftmp5 - Ftmp207 * Ftmp23 + Ftmp210 * y * M[22] + Ftmp210 * M[40] +
          Ftmp211 * y * M[24] + Ftmp211 * M[42] + Ftmp212 * y * M[31] + Ftmp212 * M[51] -
          Ftmp214 * M[41] - Ftmp215 * M[43] - Ftmp216 * M[52] + Ftmp23 * Ftmp87 -
          Ftmp23 * (Ftmp164 + Ftmp179) * M[42] - Ftmp23 * (Ftmp213 + Ftmp217) * M[40] -
          Ftmp23 * (Ftmp35 + Ftmp78) * M[15] -
          Ftmp23 * (Ftmp123 + Ftmp136 - 13230.0 * Ftmp23 * Ftmp52) * M[49] -
          Ftmp23 * (Ftmp168 + Ftmp178 + Ftmp82) * M[51] + Ftmp24 * x * M[4] +
          Ftmp24 * z * M[7] + Ftmp28 * M[17] - Ftmp3 * x + Ftmp30 * y - Ftmp34 * x -
          Ftmp35 * Ftmp62 * M[13] - Ftmp37 * z * M[16] - Ftmp4 * M[7] - Ftmp49 * M[32] -
          Ftmp5 * Ftmp65 + Ftmp56 * Ftmp62 - Ftmp57 * M[9] - Ftmp61 * Ftmp81 +
          Ftmp62 * Ftmp8;
#pragma omp atomic
  F[2] += Ftmp0 * M[2] - Ftmp101 * Ftmp222 * Ftmp26 + Ftmp112 * M[50] + Ftmp116 * z -
          Ftmp126 * Ftmp5 - Ftmp131 * M[34] - Ftmp133 * Ftmp26 -
          Ftmp134 * Ftmp223 * M[48] - Ftmp134 * Ftmp43 * y * M[53] - Ftmp137 * Ftmp223 +
          Ftmp15 * M[4] - Ftmp170 * M[37] - Ftmp171 * M[39] - Ftmp175 * Ftmp223 * M[46] -
          Ftmp180 * Ftmp2 + Ftmp180 * Ftmp27 + Ftmp181 * M[11] + Ftmp183 * z +
          Ftmp184 * Ftmp5 * Ftmp67 - Ftmp187 * y - Ftmp188 * y +
          Ftmp189 * Ftmp26 * Ftmp91 - Ftmp189 * Ftmp46 - Ftmp192 * Ftmp26 +
          Ftmp193 * Ftmp26 - Ftmp194 * Ftmp5 + Ftmp196 * Ftmp5 + Ftmp198 * Ftmp224 -
          Ftmp2 * Ftmp218 - Ftmp2 * Ftmp26 * M[2] - Ftmp201 * Ftmp26 * Ftmp9 +
          Ftmp202 * M[36] + Ftmp204 * z - Ftmp206 * Ftmp26 - Ftmp207 * Ftmp5 -
          Ftmp214 * M[40] - Ftmp215 * M[42] - Ftmp216 * M[51] + Ftmp218 * Ftmp27 +
          Ftmp219 * z * M[8] + Ftmp219 * M[18] + Ftmp220 * Ftmp225 - Ftmp220 * Ftmp48 +
          Ftmp221 * Ftmp225 - Ftmp221 * Ftmp48 + Ftmp222 * Ftmp67 - Ftmp223 * Ftmp80 +
          Ftmp224 * Ftmp77 + Ftmp226 * z * M[33] + Ftmp226 * M[54] + Ftmp227 * z * M[22] +
          Ftmp227 * M[41] + Ftmp228 * z * M[24] + Ftmp228 * M[43] + Ftmp229 * z * M[31] +
          Ftmp229 * M[52] + Ftmp25 * M[16] -
          Ftmp26 * Ftmp39 * x * (Ftmp102 - 1575.0 * Ftmp52) * M[47] - Ftmp26 * Ftmp64 -
          Ftmp26 * (Ftmp47 + Ftmp78) * M[18] -
          Ftmp26 * (Ftmp127 + Ftmp136 - 13230.0 * Ftmp71) * M[54] -
          Ftmp26 * (Ftmp162 + Ftmp18 + Ftmp213) * M[41] -
          Ftmp26 * (Ftmp163 + Ftmp172 + Ftmp217) * M[43] -
          Ftmp26 * (Ftmp167 + Ftmp172 + Ftmp176 + Ftmp82) * M[52] + Ftmp29 * z +
          Ftmp39 * Ftmp62 * Ftmp73 * M[27] - Ftmp4 * x * M[0] - Ftmp4 * y * M[1] -
          Ftmp42 * x - Ftmp47 * Ftmp9 * M[13] - Ftmp49 * x * M[14] - Ftmp49 * y * M[17] +
          Ftmp56 * Ftmp9 - Ftmp58 * z - Ftmp63 * M[9] + Ftmp75 * x + Ftmp8 * Ftmp9;
}

void field_m1_P2M_6(double x, double y, double z, double q, double* M) {
  double Mtmp0  = q * x;
  double Mtmp1  = q * y;
  double Mtmp2  = q * z;
  double Mtmp3  = (x * x);
  double Mtmp4  = (1.0 / 2.0) * q;
  double Mtmp5  = Mtmp0 * y;
  double Mtmp6  = Mtmp0 * z;
  double Mtmp7  = (y * y);
  double Mtmp8  = Mtmp1 * z;
  double Mtmp9  = (z * z);
  double Mtmp10 = (x * x * x);
  double Mtmp11 = (1.0 / 6.0) * q;
  double Mtmp12 = (1.0 / 2.0) * Mtmp3;
  double Mtmp13 = (1.0 / 2.0) * Mtmp0;
  double Mtmp14 = (y * y * y);
  double Mtmp15 = (1.0 / 2.0) * Mtmp7;
  double Mtmp16 = (1.0 / 2.0) * Mtmp9;
  double Mtmp17 = (z * z * z);
  double Mtmp18 = (x * x * x * x);
  double Mtmp19 = (1.0 / 24.0) * q;
  double Mtmp20 = (1.0 / 6.0) * Mtmp10;
  double Mtmp21 = Mtmp7 * q;
  double Mtmp22 = (1.0 / 4.0) * Mtmp3;
  double Mtmp23 = Mtmp9 * q;
  double Mtmp24 = (1.0 / 6.0) * Mtmp0;
  double Mtmp25 = (y * y * y * y);
  double Mtmp26 = (1.0 / 6.0) * Mtmp14;
  double Mtmp27 = (1.0 / 4.0) * Mtmp9;
  double Mtmp28 = (1.0 / 6.0) * Mtmp17;
  double Mtmp29 = (z * z * z * z);
  double Mtmp30 = (x * x * x * x * x);
  double Mtmp31 = (1.0 / 120.0) * q;
  double Mtmp32 = (1.0 / 24.0) * Mtmp18;
  double Mtmp33 = (1.0 / 12.0) * Mtmp10;
  double Mtmp34 = (1.0 / 12.0) * Mtmp14;
  double Mtmp35 = Mtmp3 * q;
  double Mtmp36 = Mtmp2 * Mtmp7;
  double Mtmp37 = Mtmp1 * Mtmp9;
  double Mtmp38 = (1.0 / 12.0) * Mtmp17;
  double Mtmp39 = (1.0 / 24.0) * Mtmp0;
  double Mtmp40 = Mtmp0 * Mtmp7;
  double Mtmp41 = (y * y * y * y * y);
  double Mtmp42 = (1.0 / 24.0) * Mtmp25;
  double Mtmp43 = (1.0 / 24.0) * Mtmp29;
  double Mtmp44 = (z * z * z * z * z);
  double Mtmp45 = (1.0 / 720.0) * q;
  double Mtmp46 = (1.0 / 120.0) * Mtmp30;
  double Mtmp47 = (1.0 / 48.0) * Mtmp18;
  double Mtmp48 = (1.0 / 36.0) * Mtmp10 * q;
  double Mtmp49 = (1.0 / 48.0) * Mtmp35;
  double Mtmp50 = (1.0 / 120.0) * Mtmp0;
  M[0] += -Mtmp0;
  M[1] += -Mtmp1;
  M[2] += -Mtmp2;
  M[3] += Mtmp3 * Mtmp4;
  M[4] += Mtmp5;
  M[5] += Mtmp6;
  M[6] += Mtmp4 * Mtmp7;
  M[7] += Mtmp8;
  M[8] += Mtmp4 * Mtmp9;
  M[9] += -Mtmp10 * Mtmp11;
  M[10] += -Mtmp1 * Mtmp12;
  M[11] += -Mtmp12 * Mtmp2;
  M[12] += -Mtmp13 * Mtmp7;
  M[13] += -Mtmp5 * z;
  M[14] += -Mtmp13 * Mtmp9;
  M[15] += -Mtmp11 * Mtmp14;
  M[16] += -Mtmp15 * Mtmp2;
  M[17] += -Mtmp1 * Mtmp16;
  M[18] += -Mtmp11 * Mtmp17;
  M[19] += Mtmp18 * Mtmp19;
  M[20] += Mtmp1 * Mtmp20;
  M[21] += Mtmp2 * Mtmp20;
  M[22] += Mtmp21 * Mtmp22;
  M[23] += Mtmp12 * Mtmp8;
  M[24] += Mtmp22 * Mtmp23;
  M[25] += Mtmp14 * Mtmp24;
  M[26] += Mtmp15 * Mtmp6;
  M[27] += Mtmp16 * Mtmp5;
  M[28] += Mtmp17 * Mtmp24;
  M[29] += Mtmp19 * Mtmp25;
  M[30] += Mtmp2 * Mtmp26;
  M[31] += Mtmp21 * Mtmp27;
  M[32] += Mtmp1 * Mtmp28;
  M[33] += Mtmp19 * Mtmp29;
  M[34] += -Mtmp30 * Mtmp31;
  M[35] += -Mtmp1 * Mtmp32;
  M[36] += -Mtmp2 * Mtmp32;
  M[37] += -Mtmp21 * Mtmp33;
  M[38] += -Mtmp20 * Mtmp8;
  M[39] += -Mtmp23 * Mtmp33;
  M[40] += -Mtmp34 * Mtmp35;
  M[41] += -Mtmp22 * Mtmp36;
  M[42] += -Mtmp22 * Mtmp37;
  M[43] += -Mtmp35 * Mtmp38;
  M[44] += -Mtmp25 * Mtmp39;
  M[45] += -Mtmp26 * Mtmp6;
  M[46] += -Mtmp27 * Mtmp40;
  M[47] += -Mtmp28 * Mtmp5;
  M[48] += -Mtmp29 * Mtmp39;
  M[49] += -Mtmp31 * Mtmp41;
  M[50] += -Mtmp2 * Mtmp42;
  M[51] += -Mtmp23 * Mtmp34;
  M[52] += -Mtmp21 * Mtmp38;
  M[53] += -Mtmp1 * Mtmp43;
  M[54] += -Mtmp31 * Mtmp44;
  M[55] += Mtmp45 * (x * x * x * x * x * x);
  M[56] += Mtmp1 * Mtmp46;
  M[57] += Mtmp2 * Mtmp46;
  M[58] += Mtmp21 * Mtmp47;
  M[59] += Mtmp32 * Mtmp8;
  M[60] += Mtmp23 * Mtmp47;
  M[61] += Mtmp14 * Mtmp48;
  M[62] += Mtmp33 * Mtmp36;
  M[63] += Mtmp33 * Mtmp37;
  M[64] += Mtmp17 * Mtmp48;
  M[65] += Mtmp25 * Mtmp49;
  M[66] += Mtmp2 * Mtmp3 * Mtmp34;
  M[67] += (1.0 / 8.0) * Mtmp21 * Mtmp3 * Mtmp9;
  M[68] += Mtmp1 * Mtmp3 * Mtmp38;
  M[69] += Mtmp29 * Mtmp49;
  M[70] += Mtmp41 * Mtmp50;
  M[71] += Mtmp42 * Mtmp6;
  M[72] += Mtmp0 * Mtmp34 * Mtmp9;
  M[73] += Mtmp38 * Mtmp40;
  M[74] += Mtmp43 * Mtmp5;
  M[75] += Mtmp44 * Mtmp50;
  M[76] += Mtmp45 * (y * y * y * y * y * y);
  M[77] += (1.0 / 120.0) * Mtmp2 * Mtmp41;
  M[78] += (1.0 / 48.0) * Mtmp23 * Mtmp25;
  M[79] += (1.0 / 36.0) * Mtmp14 * Mtmp17 * q;
  M[80] += (1.0 / 48.0) * Mtmp21 * Mtmp29;
  M[81] += (1.0 / 120.0) * Mtmp1 * Mtmp44;
  M[82] += Mtmp45 * (z * z * z * z * z * z);
}
void field_m1_M2M_6(double x, double y, double z, double* M, double* Ms) {
  double Mstmp0   = x * M[0];
  double Mstmp1   = x * M[1];
  double Mstmp2   = y * M[0];
  double Mstmp3   = x * M[2];
  double Mstmp4   = z * M[0];
  double Mstmp5   = y * M[1];
  double Mstmp6   = y * M[2];
  double Mstmp7   = z * M[1];
  double Mstmp8   = z * M[2];
  double Mstmp9   = x * M[3];
  double Mstmp10  = (x * x);
  double Mstmp11  = (1.0 / 2.0) * Mstmp10;
  double Mstmp12  = x * M[4];
  double Mstmp13  = y * M[3];
  double Mstmp14  = Mstmp0 * y;
  double Mstmp15  = x * M[5];
  double Mstmp16  = z * M[3];
  double Mstmp17  = Mstmp0 * z;
  double Mstmp18  = x * M[6];
  double Mstmp19  = y * M[4];
  double Mstmp20  = Mstmp1 * y;
  double Mstmp21  = (y * y);
  double Mstmp22  = (1.0 / 2.0) * M[0];
  double Mstmp23  = x * M[7];
  double Mstmp24  = y * M[5];
  double Mstmp25  = z * M[4];
  double Mstmp26  = Mstmp3 * y;
  double Mstmp27  = Mstmp1 * z;
  double Mstmp28  = Mstmp2 * z;
  double Mstmp29  = x * M[8];
  double Mstmp30  = z * M[5];
  double Mstmp31  = Mstmp3 * z;
  double Mstmp32  = (z * z);
  double Mstmp33  = y * M[6];
  double Mstmp34  = (1.0 / 2.0) * Mstmp21;
  double Mstmp35  = y * M[7];
  double Mstmp36  = z * M[6];
  double Mstmp37  = Mstmp5 * z;
  double Mstmp38  = y * M[8];
  double Mstmp39  = z * M[7];
  double Mstmp40  = Mstmp6 * z;
  double Mstmp41  = (1.0 / 2.0) * Mstmp32;
  double Mstmp42  = z * M[8];
  double Mstmp43  = x * M[9];
  double Mstmp44  = (x * x * x);
  double Mstmp45  = (1.0 / 6.0) * Mstmp44;
  double Mstmp46  = x * M[10];
  double Mstmp47  = y * M[9];
  double Mstmp48  = Mstmp9 * y;
  double Mstmp49  = x * M[11];
  double Mstmp50  = z * M[9];
  double Mstmp51  = Mstmp9 * z;
  double Mstmp52  = x * M[12];
  double Mstmp53  = y * M[10];
  double Mstmp54  = Mstmp12 * y;
  double Mstmp55  = x * M[13];
  double Mstmp56  = y * M[11];
  double Mstmp57  = z * M[10];
  double Mstmp58  = Mstmp15 * y;
  double Mstmp59  = Mstmp12 * z;
  double Mstmp60  = Mstmp13 * z;
  double Mstmp61  = x * M[14];
  double Mstmp62  = z * M[11];
  double Mstmp63  = Mstmp15 * z;
  double Mstmp64  = x * M[15];
  double Mstmp65  = y * M[12];
  double Mstmp66  = Mstmp18 * y;
  double Mstmp67  = (y * y * y);
  double Mstmp68  = (1.0 / 6.0) * M[0];
  double Mstmp69  = x * M[16];
  double Mstmp70  = y * M[13];
  double Mstmp71  = z * M[12];
  double Mstmp72  = Mstmp23 * y;
  double Mstmp73  = Mstmp18 * z;
  double Mstmp74  = Mstmp19 * z;
  double Mstmp75  = x * M[17];
  double Mstmp76  = y * M[14];
  double Mstmp77  = z * M[13];
  double Mstmp78  = Mstmp29 * y;
  double Mstmp79  = Mstmp23 * z;
  double Mstmp80  = Mstmp24 * z;
  double Mstmp81  = x * M[18];
  double Mstmp82  = z * M[14];
  double Mstmp83  = Mstmp29 * z;
  double Mstmp84  = (z * z * z);
  double Mstmp85  = y * M[15];
  double Mstmp86  = (1.0 / 6.0) * Mstmp67;
  double Mstmp87  = y * M[16];
  double Mstmp88  = z * M[15];
  double Mstmp89  = Mstmp33 * z;
  double Mstmp90  = y * M[17];
  double Mstmp91  = z * M[16];
  double Mstmp92  = Mstmp35 * z;
  double Mstmp93  = y * M[18];
  double Mstmp94  = z * M[17];
  double Mstmp95  = Mstmp38 * z;
  double Mstmp96  = (1.0 / 6.0) * Mstmp84;
  double Mstmp97  = z * M[18];
  double Mstmp98  = x * M[19];
  double Mstmp99  = (1.0 / 24.0) * (x * x * x * x);
  double Mstmp100 = x * M[20];
  double Mstmp101 = y * M[19];
  double Mstmp102 = Mstmp43 * y;
  double Mstmp103 = x * M[21];
  double Mstmp104 = x * M[22];
  double Mstmp105 = y * M[20];
  double Mstmp106 = Mstmp46 * y;
  double Mstmp107 = (1.0 / 4.0) * Mstmp10;
  double Mstmp108 = Mstmp21 * M[0];
  double Mstmp109 = x * M[23];
  double Mstmp110 = y * M[21];
  double Mstmp111 = Mstmp49 * y;
  double Mstmp112 = x * M[24];
  double Mstmp113 = Mstmp107 * Mstmp32;
  double Mstmp114 = x * M[25];
  double Mstmp115 = y * M[22];
  double Mstmp116 = Mstmp52 * y;
  double Mstmp117 = Mstmp107 * Mstmp21;
  double Mstmp118 = x * M[26];
  double Mstmp119 = y * M[23];
  double Mstmp120 = Mstmp55 * y;
  double Mstmp121 = x * M[27];
  double Mstmp122 = y * M[24];
  double Mstmp123 = Mstmp61 * y;
  double Mstmp124 = x * M[28];
  double Mstmp125 = x * M[29];
  double Mstmp126 = y * M[25];
  double Mstmp127 = Mstmp64 * y;
  double Mstmp128 = (y * y * y * y);
  double Mstmp129 = (1.0 / 24.0) * M[0];
  double Mstmp130 = x * M[30];
  double Mstmp131 = y * M[26];
  double Mstmp132 = Mstmp69 * y;
  double Mstmp133 = x * M[31];
  double Mstmp134 = y * M[27];
  double Mstmp135 = Mstmp75 * y;
  double Mstmp136 = (1.0 / 4.0) * Mstmp32;
  double Mstmp137 = x * M[32];
  double Mstmp138 = y * M[28];
  double Mstmp139 = Mstmp81 * y;
  double Mstmp140 = x * M[33];
  double Mstmp141 = (z * z * z * z);
  double Mstmp142 = y * M[29];
  double Mstmp143 = (1.0 / 24.0) * Mstmp128;
  double Mstmp144 = y * M[30];
  double Mstmp145 = y * M[31];
  double Mstmp146 = Mstmp136 * Mstmp21;
  double Mstmp147 = y * M[32];
  double Mstmp148 = y * M[33];
  double Mstmp149 = (1.0 / 24.0) * Mstmp141;
  double Mstmp150 = (1.0 / 120.0) * (x * x * x * x * x);
  double Mstmp151 = (1.0 / 12.0) * Mstmp44;
  double Mstmp152 = Mstmp151 * Mstmp32;
  double Mstmp153 = (1.0 / 12.0) * Mstmp10;
  double Mstmp154 = Mstmp153 * M[0];
  double Mstmp155 = Mstmp151 * Mstmp21;
  double Mstmp156 = Mstmp153 * Mstmp67;
  double Mstmp157 = Mstmp153 * Mstmp84;
  double Mstmp158 = (y * y * y * y * y);
  double Mstmp159 = (1.0 / 120.0) * M[0];
  double Mstmp160 = (1.0 / 12.0) * Mstmp32 * Mstmp67;
  double Mstmp161 = (1.0 / 12.0) * Mstmp84;
  double Mstmp162 = (z * z * z * z * z);
  double Mstmp163 = (1.0 / 120.0) * Mstmp158;
  double Mstmp164 = Mstmp161 * Mstmp21;
  double Mstmp165 = (1.0 / 120.0) * Mstmp162;
#pragma omp atomic
  Ms[0] += M[0];
#pragma omp atomic
  Ms[1] += M[1];
#pragma omp atomic
  Ms[2] += M[2];
#pragma omp atomic
  Ms[3] += Mstmp0 + M[3];
#pragma omp atomic
  Ms[4] += Mstmp1 + Mstmp2 + M[4];
#pragma omp atomic
  Ms[5] += Mstmp3 + Mstmp4 + M[5];
#pragma omp atomic
  Ms[6] += Mstmp5 + M[6];
#pragma omp atomic
  Ms[7] += Mstmp6 + Mstmp7 + M[7];
#pragma omp atomic
  Ms[8] += Mstmp8 + M[8];
#pragma omp atomic
  Ms[9] += Mstmp11 * M[0] + Mstmp9 + M[9];
#pragma omp atomic
  Ms[10] += Mstmp11 * M[1] + Mstmp12 + Mstmp13 + Mstmp14 + M[10];
#pragma omp atomic
  Ms[11] += Mstmp11 * M[2] + Mstmp15 + Mstmp16 + Mstmp17 + M[11];
#pragma omp atomic
  Ms[12] += Mstmp18 + Mstmp19 + Mstmp20 + Mstmp21 * Mstmp22 + M[12];
#pragma omp atomic
  Ms[13] += Mstmp23 + Mstmp24 + Mstmp25 + Mstmp26 + Mstmp27 + Mstmp28 + M[13];
#pragma omp atomic
  Ms[14] += Mstmp22 * Mstmp32 + Mstmp29 + Mstmp30 + Mstmp31 + M[14];
#pragma omp atomic
  Ms[15] += Mstmp33 + Mstmp34 * M[1] + M[15];
#pragma omp atomic
  Ms[16] += Mstmp34 * M[2] + Mstmp35 + Mstmp36 + Mstmp37 + M[16];
#pragma omp atomic
  Ms[17] += Mstmp38 + Mstmp39 + Mstmp40 + Mstmp41 * M[1] + M[17];
#pragma omp atomic
  Ms[18] += Mstmp41 * M[2] + Mstmp42 + M[18];
#pragma omp atomic
  Ms[19] += Mstmp11 * M[3] + Mstmp43 + Mstmp45 * M[0] + M[19];
#pragma omp atomic
  Ms[20] += Mstmp11 * Mstmp2 + Mstmp11 * M[4] + Mstmp45 * M[1] + Mstmp46 + Mstmp47 +
            Mstmp48 + M[20];
#pragma omp atomic
  Ms[21] += Mstmp11 * Mstmp4 + Mstmp11 * M[5] + Mstmp45 * M[2] + Mstmp49 + Mstmp50 +
            Mstmp51 + M[21];
#pragma omp atomic
  Ms[22] += Mstmp0 * Mstmp34 + Mstmp11 * Mstmp5 + Mstmp11 * M[6] + Mstmp34 * M[3] +
            Mstmp52 + Mstmp53 + Mstmp54 + M[22];
#pragma omp atomic
  Ms[23] += Mstmp11 * Mstmp6 + Mstmp11 * Mstmp7 + Mstmp11 * M[7] + Mstmp14 * z + Mstmp55 +
            Mstmp56 + Mstmp57 + Mstmp58 + Mstmp59 + Mstmp60 + M[23];
#pragma omp atomic
  Ms[24] += Mstmp0 * Mstmp41 + Mstmp11 * Mstmp8 + Mstmp11 * M[8] + Mstmp41 * M[3] +
            Mstmp61 + Mstmp62 + Mstmp63 + M[24];
#pragma omp atomic
  Ms[25] += Mstmp1 * Mstmp34 + Mstmp34 * M[4] + Mstmp64 + Mstmp65 + Mstmp66 +
            Mstmp67 * Mstmp68 + M[25];
#pragma omp atomic
  Ms[26] += Mstmp20 * z + Mstmp3 * Mstmp34 + Mstmp34 * Mstmp4 + Mstmp34 * M[5] + Mstmp69 +
            Mstmp70 + Mstmp71 + Mstmp72 + Mstmp73 + Mstmp74 + M[26];
#pragma omp atomic
  Ms[27] += Mstmp1 * Mstmp41 + Mstmp2 * Mstmp41 + Mstmp26 * z + Mstmp41 * M[4] + Mstmp75 +
            Mstmp76 + Mstmp77 + Mstmp78 + Mstmp79 + Mstmp80 + M[27];
#pragma omp atomic
  Ms[28] += Mstmp3 * Mstmp41 + Mstmp41 * M[5] + Mstmp68 * Mstmp84 + Mstmp81 + Mstmp82 +
            Mstmp83 + M[28];
#pragma omp atomic
  Ms[29] += Mstmp34 * M[6] + Mstmp85 + Mstmp86 * M[1] + M[29];
#pragma omp atomic
  Ms[30] += Mstmp34 * Mstmp7 + Mstmp34 * M[7] + Mstmp86 * M[2] + Mstmp87 + Mstmp88 +
            Mstmp89 + M[30];
#pragma omp atomic
  Ms[31] += Mstmp34 * Mstmp8 + Mstmp34 * M[8] + Mstmp41 * Mstmp5 + Mstmp41 * M[6] +
            Mstmp90 + Mstmp91 + Mstmp92 + M[31];
#pragma omp atomic
  Ms[32] += Mstmp41 * Mstmp6 + Mstmp41 * M[7] + Mstmp93 + Mstmp94 + Mstmp95 +
            Mstmp96 * M[1] + M[32];
#pragma omp atomic
  Ms[33] += Mstmp41 * M[8] + Mstmp96 * M[2] + Mstmp97 + M[33];
#pragma omp atomic
  Ms[34] += Mstmp11 * M[9] + Mstmp45 * M[3] + Mstmp98 + Mstmp99 * M[0] + M[34];
#pragma omp atomic
  Ms[35] += Mstmp100 + Mstmp101 + Mstmp102 + Mstmp11 * Mstmp13 + Mstmp11 * M[10] +
            Mstmp2 * Mstmp45 + Mstmp45 * M[4] + Mstmp99 * M[1] + M[35];
#pragma omp atomic
  Ms[36] += Mstmp103 + Mstmp11 * Mstmp16 + Mstmp11 * M[11] + Mstmp4 * Mstmp45 +
            Mstmp43 * z + Mstmp45 * M[5] + Mstmp99 * M[2] + z * M[19] + M[36];
#pragma omp atomic
  Ms[37] += Mstmp104 + Mstmp105 + Mstmp106 + Mstmp107 * Mstmp108 + Mstmp11 * Mstmp19 +
            Mstmp11 * M[12] + Mstmp34 * Mstmp9 + Mstmp34 * M[9] + Mstmp45 * Mstmp5 +
            Mstmp45 * M[6] + M[37];
#pragma omp atomic
  Ms[38] += Mstmp109 + Mstmp11 * Mstmp24 + Mstmp11 * Mstmp25 + Mstmp11 * Mstmp28 +
            Mstmp11 * M[13] + Mstmp110 + Mstmp111 + Mstmp45 * Mstmp6 + Mstmp45 * Mstmp7 +
            Mstmp45 * M[7] + Mstmp46 * z + Mstmp47 * z + Mstmp48 * z + z * M[20] + M[38];
#pragma omp atomic
  Ms[39] += Mstmp11 * Mstmp30 + Mstmp11 * M[14] + Mstmp112 + Mstmp113 * M[0] +
            Mstmp41 * Mstmp9 + Mstmp41 * M[9] + Mstmp45 * Mstmp8 + Mstmp45 * M[8] +
            Mstmp49 * z + z * M[21] + M[39];
#pragma omp atomic
  Ms[40] += Mstmp0 * Mstmp86 + Mstmp11 * Mstmp33 + Mstmp11 * M[15] + Mstmp114 + Mstmp115 +
            Mstmp116 + Mstmp117 * M[1] + Mstmp12 * Mstmp34 + Mstmp34 * M[10] +
            Mstmp86 * M[3] + M[40];
#pragma omp atomic
  Ms[41] += Mstmp11 * Mstmp35 + Mstmp11 * Mstmp36 + Mstmp11 * Mstmp37 + Mstmp11 * M[16] +
            Mstmp117 * M[2] + Mstmp118 + Mstmp119 + Mstmp120 + Mstmp15 * Mstmp34 +
            Mstmp16 * Mstmp34 + Mstmp17 * Mstmp34 + Mstmp34 * M[11] + Mstmp52 * z +
            Mstmp53 * z + Mstmp54 * z + z * M[22] + M[41];
#pragma omp atomic
  Ms[42] += Mstmp11 * Mstmp38 + Mstmp11 * Mstmp39 + Mstmp11 * Mstmp40 + Mstmp11 * M[17] +
            Mstmp113 * M[1] + Mstmp12 * Mstmp41 + Mstmp121 + Mstmp122 + Mstmp123 +
            Mstmp13 * Mstmp41 + Mstmp14 * Mstmp41 + Mstmp41 * M[10] + Mstmp55 * z +
            Mstmp56 * z + Mstmp58 * z + z * M[23] + M[42];
#pragma omp atomic
  Ms[43] += Mstmp0 * Mstmp96 + Mstmp11 * Mstmp42 + Mstmp11 * M[18] + Mstmp113 * M[2] +
            Mstmp124 + Mstmp15 * Mstmp41 + Mstmp41 * M[11] + Mstmp61 * z +
            Mstmp96 * M[3] + z * M[24] + M[43];
#pragma omp atomic
  Ms[44] += Mstmp1 * Mstmp86 + Mstmp125 + Mstmp126 + Mstmp127 + Mstmp128 * Mstmp129 +
            Mstmp18 * Mstmp34 + Mstmp34 * M[12] + Mstmp86 * M[4] + M[44];
#pragma omp atomic
  Ms[45] += Mstmp130 + Mstmp131 + Mstmp132 + Mstmp23 * Mstmp34 + Mstmp25 * Mstmp34 +
            Mstmp27 * Mstmp34 + Mstmp3 * Mstmp86 + Mstmp34 * M[13] + Mstmp4 * Mstmp86 +
            Mstmp64 * z + Mstmp65 * z + Mstmp66 * z + Mstmp86 * M[5] + z * M[25] + M[45];
#pragma omp atomic
  Ms[46] += Mstmp108 * Mstmp136 + Mstmp133 + Mstmp134 + Mstmp135 + Mstmp18 * Mstmp41 +
            Mstmp19 * Mstmp41 + Mstmp20 * Mstmp41 + Mstmp29 * Mstmp34 +
            Mstmp30 * Mstmp34 + Mstmp31 * Mstmp34 + Mstmp34 * M[14] + Mstmp41 * M[12] +
            Mstmp69 * z + Mstmp70 * z + Mstmp72 * z + z * M[26] + M[46];
#pragma omp atomic
  Ms[47] += Mstmp1 * Mstmp96 + Mstmp137 + Mstmp138 + Mstmp139 + Mstmp2 * Mstmp96 +
            Mstmp23 * Mstmp41 + Mstmp24 * Mstmp41 + Mstmp26 * Mstmp41 + Mstmp41 * M[13] +
            Mstmp75 * z + Mstmp76 * z + Mstmp78 * z + Mstmp96 * M[4] + z * M[27] + M[47];
#pragma omp atomic
  Ms[48] += Mstmp129 * Mstmp141 + Mstmp140 + Mstmp29 * Mstmp41 + Mstmp3 * Mstmp96 +
            Mstmp41 * M[14] + Mstmp81 * z + Mstmp96 * M[5] + z * M[28] + M[48];
#pragma omp atomic
  Ms[49] += Mstmp142 + Mstmp143 * M[1] + Mstmp34 * M[15] + Mstmp86 * M[6] + M[49];
#pragma omp atomic
  Ms[50] += Mstmp143 * M[2] + Mstmp144 + Mstmp34 * Mstmp36 + Mstmp34 * M[16] +
            Mstmp7 * Mstmp86 + Mstmp85 * z + Mstmp86 * M[7] + z * M[29] + M[50];
#pragma omp atomic
  Ms[51] += Mstmp145 + Mstmp146 * M[1] + Mstmp33 * Mstmp41 + Mstmp34 * Mstmp39 +
            Mstmp34 * M[17] + Mstmp41 * M[15] + Mstmp8 * Mstmp86 + Mstmp86 * M[8] +
            Mstmp87 * z + z * M[30] + M[51];
#pragma omp atomic
  Ms[52] += Mstmp146 * M[2] + Mstmp147 + Mstmp34 * Mstmp42 + Mstmp34 * M[18] +
            Mstmp35 * Mstmp41 + Mstmp41 * M[16] + Mstmp5 * Mstmp96 + Mstmp90 * z +
            Mstmp96 * M[6] + z * M[31] + M[52];
#pragma omp atomic
  Ms[53] += Mstmp148 + Mstmp149 * M[1] + Mstmp38 * Mstmp41 + Mstmp41 * M[17] +
            Mstmp6 * Mstmp96 + Mstmp93 * z + Mstmp96 * M[7] + z * M[32] + M[53];
#pragma omp atomic
  Ms[54] += Mstmp149 * M[2] + Mstmp41 * M[18] + Mstmp96 * M[8] + z * M[33] + M[54];
#pragma omp atomic
  Ms[55] += Mstmp11 * M[19] + Mstmp150 * M[0] + Mstmp45 * M[9] + Mstmp99 * M[3] +
            x * M[34] + M[55];
#pragma omp atomic
  Ms[56] += Mstmp11 * Mstmp47 + Mstmp11 * M[20] + Mstmp13 * Mstmp45 + Mstmp150 * M[1] +
            Mstmp2 * Mstmp99 + Mstmp45 * M[10] + Mstmp98 * y + Mstmp99 * M[4] +
            x * M[35] + y * M[34] + M[56];
#pragma omp atomic
  Ms[57] += Mstmp11 * Mstmp50 + Mstmp11 * M[21] + Mstmp150 * M[2] + Mstmp16 * Mstmp45 +
            Mstmp4 * Mstmp99 + Mstmp45 * M[11] + Mstmp98 * z + Mstmp99 * M[5] +
            x * M[36] + z * M[34] + M[57];
#pragma omp atomic
  Ms[58] += Mstmp100 * y + Mstmp108 * Mstmp151 + Mstmp11 * Mstmp53 + Mstmp11 * M[22] +
            Mstmp117 * M[3] + Mstmp19 * Mstmp45 + Mstmp34 * Mstmp43 + Mstmp34 * M[19] +
            Mstmp45 * M[12] + Mstmp5 * Mstmp99 + Mstmp99 * M[6] + x * M[37] + y * M[35] +
            M[58];
#pragma omp atomic
  Ms[59] += Mstmp100 * z + Mstmp101 * z + Mstmp102 * z + Mstmp103 * y +
            Mstmp11 * Mstmp56 + Mstmp11 * Mstmp57 + Mstmp11 * Mstmp60 + Mstmp11 * M[23] +
            Mstmp24 * Mstmp45 + Mstmp25 * Mstmp45 + Mstmp28 * Mstmp45 + Mstmp45 * M[13] +
            Mstmp6 * Mstmp99 + Mstmp7 * Mstmp99 + Mstmp99 * M[7] + x * M[38] + y * M[36] +
            z * M[35] + M[59];
#pragma omp atomic
  Ms[60] += Mstmp103 * z + Mstmp11 * Mstmp62 + Mstmp11 * M[24] + Mstmp113 * M[3] +
            Mstmp152 * M[0] + Mstmp30 * Mstmp45 + Mstmp41 * Mstmp43 + Mstmp41 * M[19] +
            Mstmp45 * M[14] + Mstmp8 * Mstmp99 + Mstmp99 * M[8] + x * M[39] + z * M[36] +
            M[60];
#pragma omp atomic
  Ms[61] += Mstmp104 * y + Mstmp11 * Mstmp65 + Mstmp11 * M[25] + Mstmp117 * M[4] +
            Mstmp154 * Mstmp67 + Mstmp155 * M[1] + Mstmp33 * Mstmp45 + Mstmp34 * Mstmp46 +
            Mstmp34 * M[20] + Mstmp45 * M[15] + Mstmp86 * Mstmp9 + Mstmp86 * M[9] +
            x * M[40] + y * M[37] + M[61];
#pragma omp atomic
  Ms[62] += Mstmp104 * z + Mstmp105 * z + Mstmp106 * z + Mstmp109 * y +
            Mstmp11 * Mstmp70 + Mstmp11 * Mstmp71 + Mstmp11 * Mstmp74 + Mstmp11 * M[26] +
            Mstmp117 * Mstmp4 + Mstmp117 * M[5] + Mstmp155 * M[2] + Mstmp34 * Mstmp49 +
            Mstmp34 * Mstmp50 + Mstmp34 * Mstmp51 + Mstmp34 * M[21] + Mstmp35 * Mstmp45 +
            Mstmp36 * Mstmp45 + Mstmp37 * Mstmp45 + Mstmp45 * M[16] + x * M[41] +
            y * M[38] + z * M[37] + M[62];
#pragma omp atomic
  Ms[63] += Mstmp109 * z + Mstmp11 * Mstmp76 + Mstmp11 * Mstmp77 + Mstmp11 * Mstmp80 +
            Mstmp11 * M[27] + Mstmp110 * z + Mstmp111 * z + Mstmp112 * y +
            Mstmp113 * Mstmp2 + Mstmp113 * M[4] + Mstmp152 * M[1] + Mstmp38 * Mstmp45 +
            Mstmp39 * Mstmp45 + Mstmp40 * Mstmp45 + Mstmp41 * Mstmp46 +
            Mstmp41 * Mstmp47 + Mstmp41 * Mstmp48 + Mstmp41 * M[20] + Mstmp45 * M[17] +
            x * M[42] + y * M[39] + z * M[38] + M[63];
#pragma omp atomic
  Ms[64] += Mstmp11 * Mstmp82 + Mstmp11 * M[28] + Mstmp112 * z + Mstmp113 * M[5] +
            Mstmp152 * M[2] + Mstmp154 * Mstmp84 + Mstmp41 * Mstmp49 + Mstmp41 * M[21] +
            Mstmp42 * Mstmp45 + Mstmp45 * M[18] + Mstmp9 * Mstmp96 + Mstmp96 * M[9] +
            x * M[43] + z * M[39] + M[64];
#pragma omp atomic
  Ms[65] += Mstmp0 * Mstmp143 + Mstmp11 * Mstmp85 + Mstmp11 * M[29] + Mstmp114 * y +
            Mstmp117 * M[6] + Mstmp12 * Mstmp86 + Mstmp143 * M[3] + Mstmp156 * M[1] +
            Mstmp34 * Mstmp52 + Mstmp34 * M[22] + Mstmp86 * M[10] + x * M[44] +
            y * M[40] + M[65];
#pragma omp atomic
  Ms[66] += Mstmp11 * Mstmp87 + Mstmp11 * Mstmp88 + Mstmp11 * Mstmp89 + Mstmp11 * M[30] +
            Mstmp114 * z + Mstmp115 * z + Mstmp116 * z + Mstmp117 * Mstmp7 +
            Mstmp117 * M[7] + Mstmp118 * y + Mstmp15 * Mstmp86 + Mstmp156 * M[2] +
            Mstmp16 * Mstmp86 + Mstmp17 * Mstmp86 + Mstmp34 * Mstmp55 +
            Mstmp34 * Mstmp57 + Mstmp34 * Mstmp59 + Mstmp34 * M[23] + Mstmp86 * M[11] +
            x * M[45] + y * M[41] + z * M[40] + M[66];
#pragma omp atomic
  Ms[67] += Mstmp0 * Mstmp146 + Mstmp11 * Mstmp90 + Mstmp11 * Mstmp91 +
            Mstmp11 * Mstmp92 + Mstmp11 * M[31] + Mstmp113 * Mstmp5 + Mstmp113 * M[6] +
            Mstmp117 * Mstmp8 + Mstmp117 * M[8] + Mstmp118 * z + Mstmp119 * z +
            Mstmp120 * z + Mstmp121 * y + Mstmp146 * M[3] + Mstmp34 * Mstmp61 +
            Mstmp34 * Mstmp62 + Mstmp34 * Mstmp63 + Mstmp34 * M[24] + Mstmp41 * Mstmp52 +
            Mstmp41 * Mstmp53 + Mstmp41 * Mstmp54 + Mstmp41 * M[22] + x * M[46] +
            y * M[42] + z * M[41] + M[67];
#pragma omp atomic
  Ms[68] += Mstmp11 * Mstmp93 + Mstmp11 * Mstmp94 + Mstmp11 * Mstmp95 + Mstmp11 * M[32] +
            Mstmp113 * Mstmp6 + Mstmp113 * M[7] + Mstmp12 * Mstmp96 + Mstmp121 * z +
            Mstmp122 * z + Mstmp123 * z + Mstmp124 * y + Mstmp13 * Mstmp96 +
            Mstmp14 * Mstmp96 + Mstmp157 * M[1] + Mstmp41 * Mstmp55 + Mstmp41 * Mstmp56 +
            Mstmp41 * Mstmp58 + Mstmp41 * M[23] + Mstmp96 * M[10] + x * M[47] +
            y * M[43] + z * M[42] + M[68];
#pragma omp atomic
  Ms[69] += Mstmp0 * Mstmp149 + Mstmp11 * Mstmp97 + Mstmp11 * M[33] + Mstmp113 * M[8] +
            Mstmp124 * z + Mstmp149 * M[3] + Mstmp15 * Mstmp96 + Mstmp157 * M[2] +
            Mstmp41 * Mstmp61 + Mstmp41 * M[24] + Mstmp96 * M[11] + x * M[48] +
            z * M[43] + M[69];
#pragma omp atomic
  Ms[70] += Mstmp1 * Mstmp143 + Mstmp125 * y + Mstmp143 * M[4] + Mstmp158 * Mstmp159 +
            Mstmp18 * Mstmp86 + Mstmp34 * Mstmp64 + Mstmp34 * M[25] + Mstmp86 * M[12] +
            x * M[49] + y * M[44] + M[70];
#pragma omp atomic
  Ms[71] += Mstmp125 * z + Mstmp126 * z + Mstmp127 * z + Mstmp130 * y +
            Mstmp143 * Mstmp3 + Mstmp143 * Mstmp4 + Mstmp143 * M[5] + Mstmp23 * Mstmp86 +
            Mstmp25 * Mstmp86 + Mstmp27 * Mstmp86 + Mstmp34 * Mstmp69 +
            Mstmp34 * Mstmp71 + Mstmp34 * Mstmp73 + Mstmp34 * M[26] + Mstmp86 * M[13] +
            x * M[50] + y * M[45] + z * M[44] + M[71];
#pragma omp atomic
  Ms[72] += Mstmp1 * Mstmp146 + Mstmp130 * z + Mstmp131 * z + Mstmp132 * z +
            Mstmp133 * y + Mstmp146 * M[4] + Mstmp160 * M[0] + Mstmp29 * Mstmp86 +
            Mstmp30 * Mstmp86 + Mstmp31 * Mstmp86 + Mstmp34 * Mstmp75 +
            Mstmp34 * Mstmp77 + Mstmp34 * Mstmp79 + Mstmp34 * M[27] + Mstmp41 * Mstmp64 +
            Mstmp41 * Mstmp65 + Mstmp41 * Mstmp66 + Mstmp41 * M[25] + Mstmp86 * M[14] +
            x * M[51] + y * M[46] + z * M[45] + M[72];
#pragma omp atomic
  Ms[73] += Mstmp108 * Mstmp161 + Mstmp133 * z + Mstmp134 * z + Mstmp135 * z +
            Mstmp137 * y + Mstmp146 * Mstmp3 + Mstmp146 * M[5] + Mstmp18 * Mstmp96 +
            Mstmp19 * Mstmp96 + Mstmp20 * Mstmp96 + Mstmp34 * Mstmp81 +
            Mstmp34 * Mstmp82 + Mstmp34 * Mstmp83 + Mstmp34 * M[28] + Mstmp41 * Mstmp69 +
            Mstmp41 * Mstmp70 + Mstmp41 * Mstmp72 + Mstmp41 * M[26] + Mstmp96 * M[12] +
            x * M[52] + y * M[47] + z * M[46] + M[73];
#pragma omp atomic
  Ms[74] += Mstmp1 * Mstmp149 + Mstmp137 * z + Mstmp138 * z + Mstmp139 * z +
            Mstmp140 * y + Mstmp149 * Mstmp2 + Mstmp149 * M[4] + Mstmp23 * Mstmp96 +
            Mstmp24 * Mstmp96 + Mstmp26 * Mstmp96 + Mstmp41 * Mstmp75 +
            Mstmp41 * Mstmp76 + Mstmp41 * Mstmp78 + Mstmp41 * M[27] + Mstmp96 * M[13] +
            x * M[53] + y * M[48] + z * M[47] + M[74];
#pragma omp atomic
  Ms[75] += Mstmp140 * z + Mstmp149 * Mstmp3 + Mstmp149 * M[5] + Mstmp159 * Mstmp162 +
            Mstmp29 * Mstmp96 + Mstmp41 * Mstmp81 + Mstmp41 * M[28] + Mstmp96 * M[14] +
            x * M[54] + z * M[48] + M[75];
#pragma omp atomic
  Ms[76] += Mstmp143 * M[6] + Mstmp163 * M[1] + Mstmp34 * M[29] + Mstmp86 * M[15] +
            y * M[49] + M[76];
#pragma omp atomic
  Ms[77] += Mstmp142 * z + Mstmp143 * Mstmp7 + Mstmp143 * M[7] + Mstmp163 * M[2] +
            Mstmp34 * Mstmp88 + Mstmp34 * M[30] + Mstmp36 * Mstmp86 + Mstmp86 * M[16] +
            y * M[50] + z * M[49] + M[77];
#pragma omp atomic
  Ms[78] += Mstmp143 * Mstmp8 + Mstmp143 * M[8] + Mstmp144 * z + Mstmp146 * M[6] +
            Mstmp160 * M[1] + Mstmp34 * Mstmp91 + Mstmp34 * M[31] + Mstmp39 * Mstmp86 +
            Mstmp41 * Mstmp85 + Mstmp41 * M[29] + Mstmp86 * M[17] + y * M[51] +
            z * M[50] + M[78];
#pragma omp atomic
  Ms[79] += Mstmp145 * z + Mstmp146 * M[7] + Mstmp160 * M[2] + Mstmp164 * M[1] +
            Mstmp33 * Mstmp96 + Mstmp34 * Mstmp94 + Mstmp34 * M[32] + Mstmp41 * Mstmp87 +
            Mstmp41 * M[30] + Mstmp42 * Mstmp86 + Mstmp86 * M[18] + Mstmp96 * M[15] +
            y * M[52] + z * M[51] + M[79];
#pragma omp atomic
  Ms[80] += Mstmp146 * M[8] + Mstmp147 * z + Mstmp149 * Mstmp5 + Mstmp149 * M[6] +
            Mstmp164 * M[2] + Mstmp34 * Mstmp97 + Mstmp34 * M[33] + Mstmp35 * Mstmp96 +
            Mstmp41 * Mstmp90 + Mstmp41 * M[31] + Mstmp96 * M[16] + y * M[53] +
            z * M[52] + M[80];
#pragma omp atomic
  Ms[81] += Mstmp148 * z + Mstmp149 * Mstmp6 + Mstmp149 * M[7] + Mstmp165 * M[1] +
            Mstmp38 * Mstmp96 + Mstmp41 * Mstmp93 + Mstmp41 * M[32] + Mstmp96 * M[17] +
            y * M[54] + z * M[53] + M[81];
#pragma omp atomic
  Ms[82] += Mstmp149 * M[8] + Mstmp165 * M[2] + Mstmp41 * M[33] + Mstmp96 * M[18] +
            z * M[54] + M[82];
}

void field_m1_M2L_6(double x, double y, double z, double* M, double* L) {
  double R = sqrt(x * x + y * y + z * z);
  double D[83];
  double Dtmp0  = 1.0 * pow(R, -3.0);
  double Dtmp1  = -Dtmp0;
  double Dtmp2  = (x * x);
  double Dtmp3  = pow(R, -5.0);
  double Dtmp4  = 3.0 * Dtmp3;
  double Dtmp5  = x * y;
  double Dtmp6  = x * z;
  double Dtmp7  = (y * y);
  double Dtmp8  = y * z;
  double Dtmp9  = 9.0 * Dtmp3;
  double Dtmp10 = -Dtmp9;
  double Dtmp11 = pow(R, -7.0);
  double Dtmp12 = 15.0 * Dtmp11;
  double Dtmp13 = Dtmp12 * Dtmp2;
  double Dtmp14 = -Dtmp4;
  double Dtmp15 = Dtmp13 + Dtmp14;
  double Dtmp16 = Dtmp12 * Dtmp7;
  double Dtmp17 = Dtmp14 + Dtmp16;
  double Dtmp18 = 1.0 * x;
  double Dtmp19 = Dtmp8 * x;
  double Dtmp20 = (x * x * x * x);
  double Dtmp21 = pow(R, -9.0);
  double Dtmp22 = 105.0 * Dtmp21;
  double Dtmp23 = 90.0 * Dtmp11;
  double Dtmp24 = 45.0 * Dtmp11;
  double Dtmp25 = -Dtmp24;
  double Dtmp26 = Dtmp2 * Dtmp22;
  double Dtmp27 = x * (Dtmp25 + Dtmp26);
  double Dtmp28 = -Dtmp12;
  double Dtmp29 = Dtmp22 * Dtmp7;
  double Dtmp30 = Dtmp25 + Dtmp29;
  double Dtmp31 = Dtmp18 * y;
  double Dtmp32 = Dtmp18 * z;
  double Dtmp33 = (y * y * y * y);
  double Dtmp34 = 225.0 * Dtmp11;
  double Dtmp35 = pow(R, -11.0);
  double Dtmp36 = 945.0 * Dtmp35;
  double Dtmp37 = Dtmp20 * Dtmp36;
  double Dtmp38 = Dtmp2 * Dtmp21;
  double Dtmp39 = 630.0 * Dtmp38;
  double Dtmp40 = Dtmp24 + Dtmp37 - Dtmp39;
  double Dtmp41 = -Dtmp26;
  double Dtmp42 = 315.0 * Dtmp21;
  double Dtmp43 = Dtmp42 * Dtmp7;
  double Dtmp44 = Dtmp2 * Dtmp36;
  double Dtmp45 = Dtmp44 * Dtmp7;
  double Dtmp46 = Dtmp24 + Dtmp45;
  double Dtmp47 = -Dtmp42;
  double Dtmp48 = -Dtmp29;
  double Dtmp49 = Dtmp2 * Dtmp42;
  double Dtmp50 = Dtmp33 * Dtmp36;
  double Dtmp51 = Dtmp21 * Dtmp7;
  double Dtmp52 = 630.0 * Dtmp51;
  double Dtmp53 = Dtmp24 + Dtmp50 - Dtmp52;
  double Dtmp54 = Dtmp36 * Dtmp7;
  double Dtmp55 = -Dtmp34;
  double Dtmp56 = 10395.0 * pow(R, -13.0);
  double Dtmp57 = 14175.0 * Dtmp35;
  double Dtmp58 = 1575.0 * Dtmp21;
  double Dtmp59 = Dtmp20 * Dtmp56;
  double Dtmp60 = Dtmp2 * Dtmp35;
  double Dtmp61 = x * (Dtmp58 + Dtmp59 - 9450.0 * Dtmp60);
  double Dtmp62 = 5670.0 * Dtmp60;
  double Dtmp63 = Dtmp25 - Dtmp62 * Dtmp7;
  double Dtmp64 = -2835.0 * Dtmp60;
  double Dtmp65 = Dtmp35 * Dtmp7;
  double Dtmp66 = Dtmp2 * Dtmp56 * Dtmp7;
  double Dtmp67 = -2835.0 * Dtmp65 + Dtmp66;
  double Dtmp68 = Dtmp33 * Dtmp56;
  double Dtmp69 = Dtmp58 - 9450.0 * Dtmp65 + Dtmp68;
  D[0]          = -Dtmp0 * x;
  D[1]          = -Dtmp0 * y;
  D[2]          = -Dtmp0 * z;
  D[3]          = Dtmp1 + Dtmp2 * Dtmp4;
  D[4]          = Dtmp4 * Dtmp5;
  D[5]          = Dtmp4 * Dtmp6;
  D[6]          = Dtmp1 + Dtmp4 * Dtmp7;
  D[7]          = Dtmp4 * Dtmp8;
  D[8]          = -D[3] - D[6];
  D[9]          = -x * (Dtmp10 + Dtmp13);
  D[10]         = -Dtmp15 * y;
  D[11]         = -Dtmp15 * z;
  D[12]         = -Dtmp17 * Dtmp18;
  D[13]         = -Dtmp12 * Dtmp19;
  D[14]         = -D[9] - D[12];
  D[15]         = -y * (Dtmp10 + Dtmp16);
  D[16]         = -Dtmp17 * z;
  D[17]         = -D[10] - D[15];
  D[18]         = -D[11] - D[16];
  D[19]         = -Dtmp2 * Dtmp23 + Dtmp20 * Dtmp22 + Dtmp9;
  D[20]         = Dtmp27 * y;
  D[21]         = Dtmp27 * z;
  D[22]         = -Dtmp13 - Dtmp16 + Dtmp26 * Dtmp7 + Dtmp4;
  D[23]         = Dtmp8 * (Dtmp26 + Dtmp28);
  D[24]         = -D[19] - D[22];
  D[25]         = Dtmp30 * Dtmp31;
  D[26]         = Dtmp32 * (Dtmp28 + Dtmp29);
  D[27]         = -D[20] - D[25];
  D[28]         = -D[21] - D[26];
  D[29]         = Dtmp22 * Dtmp33 - Dtmp23 * Dtmp7 + Dtmp9;
  D[30]         = Dtmp30 * Dtmp8;
  D[31]         = -D[22] - D[29];
  D[32]         = -D[23] - D[30];
  D[33]         = -D[24] - D[31];
  D[34]         = -x * (Dtmp34 + Dtmp37 - 1050.0 * Dtmp38);
  D[35]         = -Dtmp40 * y;
  D[36]         = -Dtmp40 * z;
  D[37]         = -x * (Dtmp41 - Dtmp43 + Dtmp46);
  D[38]         = -Dtmp19 * (Dtmp44 + Dtmp47);
  D[39]         = -D[34] - D[37];
  D[40]         = -y * (Dtmp46 + Dtmp48 - Dtmp49);
  D[41]         = -z * (Dtmp12 + Dtmp41 + Dtmp45 + Dtmp48);
  D[42]         = -D[35] - D[40];
  D[43]         = -D[36] - D[41];
  D[44]         = -Dtmp18 * Dtmp53;
  D[45]         = -Dtmp18 * Dtmp8 * (Dtmp47 + Dtmp54);
  D[46]         = -D[37] - D[44];
  D[47]         = -D[38] - D[45];
  D[48]         = -D[39] - D[46];
  D[49]         = -y * (Dtmp34 + Dtmp50 - 1050.0 * Dtmp51);
  D[50]         = -Dtmp53 * z;
  D[51]         = -D[40] - D[49];
  D[52]         = -D[41] - D[50];
  D[53]         = -D[42] - D[51];
  D[54]         = -D[43] - D[52];
  D[55] = -Dtmp20 * Dtmp57 + 4725.0 * Dtmp38 + Dtmp55 + Dtmp56 * (x * x * x * x * x * x);
  D[56] = Dtmp61 * y;
  D[57] = Dtmp61 * z;
  D[58] = -Dtmp37 + Dtmp39 + Dtmp43 + Dtmp59 * Dtmp7 + Dtmp63;
  D[59] = Dtmp8 * (Dtmp42 + Dtmp59 - Dtmp62);
  D[60] = -D[55] - D[58];
  D[61] = Dtmp5 * (945.0 * Dtmp21 + Dtmp64 + Dtmp67);
  D[62] = Dtmp6 * (Dtmp42 - Dtmp44 + Dtmp67);
  D[63] = -D[56] - D[61];
  D[64] = -D[57] - D[62];
  D[65] = Dtmp2 * Dtmp68 + Dtmp49 - Dtmp50 + Dtmp52 + Dtmp63;
  D[66] = Dtmp8 * (Dtmp42 - Dtmp54 + Dtmp64 + Dtmp66);
  D[67] = -D[58] - D[65];
  D[68] = -D[59] - D[66];
  D[69] = -D[60] - D[67];
  D[70] = Dtmp31 * Dtmp69;
  D[71] = Dtmp32 * (Dtmp42 - 5670.0 * Dtmp65 + Dtmp68);
  D[72] = -D[61] - D[70];
  D[73] = -D[62] - D[71];
  D[74] = -D[63] - D[72];
  D[75] = -D[64] - D[73];
  D[76] = -Dtmp33 * Dtmp57 + 4725.0 * Dtmp51 + Dtmp55 + Dtmp56 * (y * y * y * y * y * y);
  D[77] = Dtmp69 * Dtmp8;
  D[78] = -D[65] - D[76];
  D[79] = -D[66] - D[77];
  D[80] = -D[67] - D[78];
  D[81] = -D[68] - D[79];
  D[82] = -D[69] - D[80];
#pragma omp atomic
  L[0] += D[0] * M[0] + D[1] * M[1] + D[2] * M[2] + D[3] * M[3] + D[4] * M[4] +
          D[5] * M[5] + D[6] * M[6] + D[7] * M[7] + D[8] * M[8] + D[9] * M[9] +
          D[10] * M[10] + D[11] * M[11] + D[12] * M[12] + D[13] * M[13] + D[14] * M[14] +
          D[15] * M[15] + D[16] * M[16] + D[17] * M[17] + D[18] * M[18] + D[19] * M[19] +
          D[20] * M[20] + D[21] * M[21] + D[22] * M[22] + D[23] * M[23] + D[24] * M[24] +
          D[25] * M[25] + D[26] * M[26] + D[27] * M[27] + D[28] * M[28] + D[29] * M[29] +
          D[30] * M[30] + D[31] * M[31] + D[32] * M[32] + D[33] * M[33] + D[34] * M[34] +
          D[35] * M[35] + D[36] * M[36] + D[37] * M[37] + D[38] * M[38] + D[39] * M[39] +
          D[40] * M[40] + D[41] * M[41] + D[42] * M[42] + D[43] * M[43] + D[44] * M[44] +
          D[45] * M[45] + D[46] * M[46] + D[47] * M[47] + D[48] * M[48] + D[49] * M[49] +
          D[50] * M[50] + D[51] * M[51] + D[52] * M[52] + D[53] * M[53] + D[54] * M[54] +
          D[55] * M[55] + D[56] * M[56] + D[57] * M[57] + D[58] * M[58] + D[59] * M[59] +
          D[60] * M[60] + D[61] * M[61] + D[62] * M[62] + D[63] * M[63] + D[64] * M[64] +
          D[65] * M[65] + D[66] * M[66] + D[67] * M[67] + D[68] * M[68] + D[69] * M[69] +
          D[70] * M[70] + D[71] * M[71] + D[72] * M[72] + D[73] * M[73] + D[74] * M[74] +
          D[75] * M[75] + D[76] * M[76] + D[77] * M[77] + D[78] * M[78] + D[79] * M[79] +
          D[80] * M[80] + D[81] * M[81] + D[82] * M[82];
#pragma omp atomic
  L[1] += D[3] * M[0] + D[4] * M[1] + D[5] * M[2] + D[9] * M[3] + D[10] * M[4] +
          D[11] * M[5] + D[12] * M[6] + D[13] * M[7] + D[14] * M[8] + D[19] * M[9] +
          D[20] * M[10] + D[21] * M[11] + D[22] * M[12] + D[23] * M[13] + D[24] * M[14] +
          D[25] * M[15] + D[26] * M[16] + D[27] * M[17] + D[28] * M[18] + D[34] * M[19] +
          D[35] * M[20] + D[36] * M[21] + D[37] * M[22] + D[38] * M[23] + D[39] * M[24] +
          D[40] * M[25] + D[41] * M[26] + D[42] * M[27] + D[43] * M[28] + D[44] * M[29] +
          D[45] * M[30] + D[46] * M[31] + D[47] * M[32] + D[48] * M[33] + D[55] * M[34] +
          D[56] * M[35] + D[57] * M[36] + D[58] * M[37] + D[59] * M[38] + D[60] * M[39] +
          D[61] * M[40] + D[62] * M[41] + D[63] * M[42] + D[64] * M[43] + D[65] * M[44] +
          D[66] * M[45] + D[67] * M[46] + D[68] * M[47] + D[69] * M[48] + D[70] * M[49] +
          D[71] * M[50] + D[72] * M[51] + D[73] * M[52] + D[74] * M[53] + D[75] * M[54];
#pragma omp atomic
  L[2] += D[4] * M[0] + D[6] * M[1] + D[7] * M[2] + D[10] * M[3] + D[12] * M[4] +
          D[13] * M[5] + D[15] * M[6] + D[16] * M[7] + D[17] * M[8] + D[20] * M[9] +
          D[22] * M[10] + D[23] * M[11] + D[25] * M[12] + D[26] * M[13] + D[27] * M[14] +
          D[29] * M[15] + D[30] * M[16] + D[31] * M[17] + D[32] * M[18] + D[35] * M[19] +
          D[37] * M[20] + D[38] * M[21] + D[40] * M[22] + D[41] * M[23] + D[42] * M[24] +
          D[44] * M[25] + D[45] * M[26] + D[46] * M[27] + D[47] * M[28] + D[49] * M[29] +
          D[50] * M[30] + D[51] * M[31] + D[52] * M[32] + D[53] * M[33] + D[56] * M[34] +
          D[58] * M[35] + D[59] * M[36] + D[61] * M[37] + D[62] * M[38] + D[63] * M[39] +
          D[65] * M[40] + D[66] * M[41] + D[67] * M[42] + D[68] * M[43] + D[70] * M[44] +
          D[71] * M[45] + D[72] * M[46] + D[73] * M[47] + D[74] * M[48] + D[76] * M[49] +
          D[77] * M[50] + D[78] * M[51] + D[79] * M[52] + D[80] * M[53] + D[81] * M[54];
#pragma omp atomic
  L[3] += D[5] * M[0] + D[7] * M[1] + D[8] * M[2] + D[11] * M[3] + D[13] * M[4] +
          D[14] * M[5] + D[16] * M[6] + D[17] * M[7] + D[18] * M[8] + D[21] * M[9] +
          D[23] * M[10] + D[24] * M[11] + D[26] * M[12] + D[27] * M[13] + D[28] * M[14] +
          D[30] * M[15] + D[31] * M[16] + D[32] * M[17] + D[33] * M[18] + D[36] * M[19] +
          D[38] * M[20] + D[39] * M[21] + D[41] * M[22] + D[42] * M[23] + D[43] * M[24] +
          D[45] * M[25] + D[46] * M[26] + D[47] * M[27] + D[48] * M[28] + D[50] * M[29] +
          D[51] * M[30] + D[52] * M[31] + D[53] * M[32] + D[54] * M[33] + D[57] * M[34] +
          D[59] * M[35] + D[60] * M[36] + D[62] * M[37] + D[63] * M[38] + D[64] * M[39] +
          D[66] * M[40] + D[67] * M[41] + D[68] * M[42] + D[69] * M[43] + D[71] * M[44] +
          D[72] * M[45] + D[73] * M[46] + D[74] * M[47] + D[75] * M[48] + D[77] * M[49] +
          D[78] * M[50] + D[79] * M[51] + D[80] * M[52] + D[81] * M[53] + D[82] * M[54];
#pragma omp atomic
  L[4] += D[9] * M[0] + D[10] * M[1] + D[11] * M[2] + D[19] * M[3] + D[20] * M[4] +
          D[21] * M[5] + D[22] * M[6] + D[23] * M[7] + D[24] * M[8] + D[34] * M[9] +
          D[35] * M[10] + D[36] * M[11] + D[37] * M[12] + D[38] * M[13] + D[39] * M[14] +
          D[40] * M[15] + D[41] * M[16] + D[42] * M[17] + D[43] * M[18] + D[55] * M[19] +
          D[56] * M[20] + D[57] * M[21] + D[58] * M[22] + D[59] * M[23] + D[60] * M[24] +
          D[61] * M[25] + D[62] * M[26] + D[63] * M[27] + D[64] * M[28] + D[65] * M[29] +
          D[66] * M[30] + D[67] * M[31] + D[68] * M[32] + D[69] * M[33];
#pragma omp atomic
  L[5] += D[10] * M[0] + D[12] * M[1] + D[13] * M[2] + D[20] * M[3] + D[22] * M[4] +
          D[23] * M[5] + D[25] * M[6] + D[26] * M[7] + D[27] * M[8] + D[35] * M[9] +
          D[37] * M[10] + D[38] * M[11] + D[40] * M[12] + D[41] * M[13] + D[42] * M[14] +
          D[44] * M[15] + D[45] * M[16] + D[46] * M[17] + D[47] * M[18] + D[56] * M[19] +
          D[58] * M[20] + D[59] * M[21] + D[61] * M[22] + D[62] * M[23] + D[63] * M[24] +
          D[65] * M[25] + D[66] * M[26] + D[67] * M[27] + D[68] * M[28] + D[70] * M[29] +
          D[71] * M[30] + D[72] * M[31] + D[73] * M[32] + D[74] * M[33];
#pragma omp atomic
  L[6] += D[11] * M[0] + D[13] * M[1] + D[14] * M[2] + D[21] * M[3] + D[23] * M[4] +
          D[24] * M[5] + D[26] * M[6] + D[27] * M[7] + D[28] * M[8] + D[36] * M[9] +
          D[38] * M[10] + D[39] * M[11] + D[41] * M[12] + D[42] * M[13] + D[43] * M[14] +
          D[45] * M[15] + D[46] * M[16] + D[47] * M[17] + D[48] * M[18] + D[57] * M[19] +
          D[59] * M[20] + D[60] * M[21] + D[62] * M[22] + D[63] * M[23] + D[64] * M[24] +
          D[66] * M[25] + D[67] * M[26] + D[68] * M[27] + D[69] * M[28] + D[71] * M[29] +
          D[72] * M[30] + D[73] * M[31] + D[74] * M[32] + D[75] * M[33];
#pragma omp atomic
  L[7] += D[12] * M[0] + D[15] * M[1] + D[16] * M[2] + D[22] * M[3] + D[25] * M[4] +
          D[26] * M[5] + D[29] * M[6] + D[30] * M[7] + D[31] * M[8] + D[37] * M[9] +
          D[40] * M[10] + D[41] * M[11] + D[44] * M[12] + D[45] * M[13] + D[46] * M[14] +
          D[49] * M[15] + D[50] * M[16] + D[51] * M[17] + D[52] * M[18] + D[58] * M[19] +
          D[61] * M[20] + D[62] * M[21] + D[65] * M[22] + D[66] * M[23] + D[67] * M[24] +
          D[70] * M[25] + D[71] * M[26] + D[72] * M[27] + D[73] * M[28] + D[76] * M[29] +
          D[77] * M[30] + D[78] * M[31] + D[79] * M[32] + D[80] * M[33];
#pragma omp atomic
  L[8] += D[13] * M[0] + D[16] * M[1] + D[17] * M[2] + D[23] * M[3] + D[26] * M[4] +
          D[27] * M[5] + D[30] * M[6] + D[31] * M[7] + D[32] * M[8] + D[38] * M[9] +
          D[41] * M[10] + D[42] * M[11] + D[45] * M[12] + D[46] * M[13] + D[47] * M[14] +
          D[50] * M[15] + D[51] * M[16] + D[52] * M[17] + D[53] * M[18] + D[59] * M[19] +
          D[62] * M[20] + D[63] * M[21] + D[66] * M[22] + D[67] * M[23] + D[68] * M[24] +
          D[71] * M[25] + D[72] * M[26] + D[73] * M[27] + D[74] * M[28] + D[77] * M[29] +
          D[78] * M[30] + D[79] * M[31] + D[80] * M[32] + D[81] * M[33];
#pragma omp atomic
  L[9] += D[14] * M[0] + D[17] * M[1] + D[18] * M[2] + D[24] * M[3] + D[27] * M[4] +
          D[28] * M[5] + D[31] * M[6] + D[32] * M[7] + D[33] * M[8] + D[39] * M[9] +
          D[42] * M[10] + D[43] * M[11] + D[46] * M[12] + D[47] * M[13] + D[48] * M[14] +
          D[51] * M[15] + D[52] * M[16] + D[53] * M[17] + D[54] * M[18] + D[60] * M[19] +
          D[63] * M[20] + D[64] * M[21] + D[67] * M[22] + D[68] * M[23] + D[69] * M[24] +
          D[72] * M[25] + D[73] * M[26] + D[74] * M[27] + D[75] * M[28] + D[78] * M[29] +
          D[79] * M[30] + D[80] * M[31] + D[81] * M[32] + D[82] * M[33];
#pragma omp atomic
  L[10] += D[19] * M[0] + D[20] * M[1] + D[21] * M[2] + D[34] * M[3] + D[35] * M[4] +
           D[36] * M[5] + D[37] * M[6] + D[38] * M[7] + D[39] * M[8] + D[55] * M[9] +
           D[56] * M[10] + D[57] * M[11] + D[58] * M[12] + D[59] * M[13] + D[60] * M[14] +
           D[61] * M[15] + D[62] * M[16] + D[63] * M[17] + D[64] * M[18];
#pragma omp atomic
  L[11] += D[20] * M[0] + D[22] * M[1] + D[23] * M[2] + D[35] * M[3] + D[37] * M[4] +
           D[38] * M[5] + D[40] * M[6] + D[41] * M[7] + D[42] * M[8] + D[56] * M[9] +
           D[58] * M[10] + D[59] * M[11] + D[61] * M[12] + D[62] * M[13] + D[63] * M[14] +
           D[65] * M[15] + D[66] * M[16] + D[67] * M[17] + D[68] * M[18];
#pragma omp atomic
  L[12] += D[21] * M[0] + D[23] * M[1] + D[24] * M[2] + D[36] * M[3] + D[38] * M[4] +
           D[39] * M[5] + D[41] * M[6] + D[42] * M[7] + D[43] * M[8] + D[57] * M[9] +
           D[59] * M[10] + D[60] * M[11] + D[62] * M[12] + D[63] * M[13] + D[64] * M[14] +
           D[66] * M[15] + D[67] * M[16] + D[68] * M[17] + D[69] * M[18];
#pragma omp atomic
  L[13] += D[22] * M[0] + D[25] * M[1] + D[26] * M[2] + D[37] * M[3] + D[40] * M[4] +
           D[41] * M[5] + D[44] * M[6] + D[45] * M[7] + D[46] * M[8] + D[58] * M[9] +
           D[61] * M[10] + D[62] * M[11] + D[65] * M[12] + D[66] * M[13] + D[67] * M[14] +
           D[70] * M[15] + D[71] * M[16] + D[72] * M[17] + D[73] * M[18];
#pragma omp atomic
  L[14] += D[23] * M[0] + D[26] * M[1] + D[27] * M[2] + D[38] * M[3] + D[41] * M[4] +
           D[42] * M[5] + D[45] * M[6] + D[46] * M[7] + D[47] * M[8] + D[59] * M[9] +
           D[62] * M[10] + D[63] * M[11] + D[66] * M[12] + D[67] * M[13] + D[68] * M[14] +
           D[71] * M[15] + D[72] * M[16] + D[73] * M[17] + D[74] * M[18];
#pragma omp atomic
  L[15] += D[24] * M[0] + D[27] * M[1] + D[28] * M[2] + D[39] * M[3] + D[42] * M[4] +
           D[43] * M[5] + D[46] * M[6] + D[47] * M[7] + D[48] * M[8] + D[60] * M[9] +
           D[63] * M[10] + D[64] * M[11] + D[67] * M[12] + D[68] * M[13] + D[69] * M[14] +
           D[72] * M[15] + D[73] * M[16] + D[74] * M[17] + D[75] * M[18];
#pragma omp atomic
  L[16] += D[25] * M[0] + D[29] * M[1] + D[30] * M[2] + D[40] * M[3] + D[44] * M[4] +
           D[45] * M[5] + D[49] * M[6] + D[50] * M[7] + D[51] * M[8] + D[61] * M[9] +
           D[65] * M[10] + D[66] * M[11] + D[70] * M[12] + D[71] * M[13] + D[72] * M[14] +
           D[76] * M[15] + D[77] * M[16] + D[78] * M[17] + D[79] * M[18];
#pragma omp atomic
  L[17] += D[26] * M[0] + D[30] * M[1] + D[31] * M[2] + D[41] * M[3] + D[45] * M[4] +
           D[46] * M[5] + D[50] * M[6] + D[51] * M[7] + D[52] * M[8] + D[62] * M[9] +
           D[66] * M[10] + D[67] * M[11] + D[71] * M[12] + D[72] * M[13] + D[73] * M[14] +
           D[77] * M[15] + D[78] * M[16] + D[79] * M[17] + D[80] * M[18];
#pragma omp atomic
  L[18] += D[27] * M[0] + D[31] * M[1] + D[32] * M[2] + D[42] * M[3] + D[46] * M[4] +
           D[47] * M[5] + D[51] * M[6] + D[52] * M[7] + D[53] * M[8] + D[63] * M[9] +
           D[67] * M[10] + D[68] * M[11] + D[72] * M[12] + D[73] * M[13] + D[74] * M[14] +
           D[78] * M[15] + D[79] * M[16] + D[80] * M[17] + D[81] * M[18];
#pragma omp atomic
  L[19] += D[28] * M[0] + D[32] * M[1] + D[33] * M[2] + D[43] * M[3] + D[47] * M[4] +
           D[48] * M[5] + D[52] * M[6] + D[53] * M[7] + D[54] * M[8] + D[64] * M[9] +
           D[68] * M[10] + D[69] * M[11] + D[73] * M[12] + D[74] * M[13] + D[75] * M[14] +
           D[79] * M[15] + D[80] * M[16] + D[81] * M[17] + D[82] * M[18];
#pragma omp atomic
  L[20] += D[34] * M[0] + D[35] * M[1] + D[36] * M[2] + D[55] * M[3] + D[56] * M[4] +
           D[57] * M[5] + D[58] * M[6] + D[59] * M[7] + D[60] * M[8];
#pragma omp atomic
  L[21] += D[35] * M[0] + D[37] * M[1] + D[38] * M[2] + D[56] * M[3] + D[58] * M[4] +
           D[59] * M[5] + D[61] * M[6] + D[62] * M[7] + D[63] * M[8];
#pragma omp atomic
  L[22] += D[36] * M[0] + D[38] * M[1] + D[39] * M[2] + D[57] * M[3] + D[59] * M[4] +
           D[60] * M[5] + D[62] * M[6] + D[63] * M[7] + D[64] * M[8];
#pragma omp atomic
  L[23] += D[37] * M[0] + D[40] * M[1] + D[41] * M[2] + D[58] * M[3] + D[61] * M[4] +
           D[62] * M[5] + D[65] * M[6] + D[66] * M[7] + D[67] * M[8];
#pragma omp atomic
  L[24] += D[38] * M[0] + D[41] * M[1] + D[42] * M[2] + D[59] * M[3] + D[62] * M[4] +
           D[63] * M[5] + D[66] * M[6] + D[67] * M[7] + D[68] * M[8];
#pragma omp atomic
  L[25] += D[39] * M[0] + D[42] * M[1] + D[43] * M[2] + D[60] * M[3] + D[63] * M[4] +
           D[64] * M[5] + D[67] * M[6] + D[68] * M[7] + D[69] * M[8];
#pragma omp atomic
  L[26] += D[40] * M[0] + D[44] * M[1] + D[45] * M[2] + D[61] * M[3] + D[65] * M[4] +
           D[66] * M[5] + D[70] * M[6] + D[71] * M[7] + D[72] * M[8];
#pragma omp atomic
  L[27] += D[41] * M[0] + D[45] * M[1] + D[46] * M[2] + D[62] * M[3] + D[66] * M[4] +
           D[67] * M[5] + D[71] * M[6] + D[72] * M[7] + D[73] * M[8];
#pragma omp atomic
  L[28] += D[42] * M[0] + D[46] * M[1] + D[47] * M[2] + D[63] * M[3] + D[67] * M[4] +
           D[68] * M[5] + D[72] * M[6] + D[73] * M[7] + D[74] * M[8];
#pragma omp atomic
  L[29] += D[43] * M[0] + D[47] * M[1] + D[48] * M[2] + D[64] * M[3] + D[68] * M[4] +
           D[69] * M[5] + D[73] * M[6] + D[74] * M[7] + D[75] * M[8];
#pragma omp atomic
  L[30] += D[44] * M[0] + D[49] * M[1] + D[50] * M[2] + D[65] * M[3] + D[70] * M[4] +
           D[71] * M[5] + D[76] * M[6] + D[77] * M[7] + D[78] * M[8];
#pragma omp atomic
  L[31] += D[45] * M[0] + D[50] * M[1] + D[51] * M[2] + D[66] * M[3] + D[71] * M[4] +
           D[72] * M[5] + D[77] * M[6] + D[78] * M[7] + D[79] * M[8];
#pragma omp atomic
  L[32] += D[46] * M[0] + D[51] * M[1] + D[52] * M[2] + D[67] * M[3] + D[72] * M[4] +
           D[73] * M[5] + D[78] * M[6] + D[79] * M[7] + D[80] * M[8];
#pragma omp atomic
  L[33] += D[47] * M[0] + D[52] * M[1] + D[53] * M[2] + D[68] * M[3] + D[73] * M[4] +
           D[74] * M[5] + D[79] * M[6] + D[80] * M[7] + D[81] * M[8];
#pragma omp atomic
  L[34] += D[48] * M[0] + D[53] * M[1] + D[54] * M[2] + D[69] * M[3] + D[74] * M[4] +
           D[75] * M[5] + D[80] * M[6] + D[81] * M[7] + D[82] * M[8];
#pragma omp atomic
  L[35] += D[55] * M[0] + D[56] * M[1] + D[57] * M[2];
#pragma omp atomic
  L[36] += D[56] * M[0] + D[58] * M[1] + D[59] * M[2];
#pragma omp atomic
  L[37] += D[57] * M[0] + D[59] * M[1] + D[60] * M[2];
#pragma omp atomic
  L[38] += D[58] * M[0] + D[61] * M[1] + D[62] * M[2];
#pragma omp atomic
  L[39] += D[59] * M[0] + D[62] * M[1] + D[63] * M[2];
#pragma omp atomic
  L[40] += D[60] * M[0] + D[63] * M[1] + D[64] * M[2];
#pragma omp atomic
  L[41] += D[61] * M[0] + D[65] * M[1] + D[66] * M[2];
#pragma omp atomic
  L[42] += D[62] * M[0] + D[66] * M[1] + D[67] * M[2];
#pragma omp atomic
  L[43] += D[63] * M[0] + D[67] * M[1] + D[68] * M[2];
#pragma omp atomic
  L[44] += D[64] * M[0] + D[68] * M[1] + D[69] * M[2];
#pragma omp atomic
  L[45] += D[65] * M[0] + D[70] * M[1] + D[71] * M[2];
#pragma omp atomic
  L[46] += D[66] * M[0] + D[71] * M[1] + D[72] * M[2];
#pragma omp atomic
  L[47] += D[67] * M[0] + D[72] * M[1] + D[73] * M[2];
#pragma omp atomic
  L[48] += D[68] * M[0] + D[73] * M[1] + D[74] * M[2];
#pragma omp atomic
  L[49] += D[69] * M[0] + D[74] * M[1] + D[75] * M[2];
#pragma omp atomic
  L[50] += D[70] * M[0] + D[76] * M[1] + D[77] * M[2];
#pragma omp atomic
  L[51] += D[71] * M[0] + D[77] * M[1] + D[78] * M[2];
#pragma omp atomic
  L[52] += D[72] * M[0] + D[78] * M[1] + D[79] * M[2];
#pragma omp atomic
  L[53] += D[73] * M[0] + D[79] * M[1] + D[80] * M[2];
#pragma omp atomic
  L[54] += D[74] * M[0] + D[80] * M[1] + D[81] * M[2];
#pragma omp atomic
  L[55] += D[75] * M[0] + D[81] * M[1] + D[82] * M[2];
}

void field_m1_L2L_6(double x, double y, double z, double* L, double* Ls) {
  double Lstmp0   = y * L[5];
  double Lstmp1   = z * L[6];
  double Lstmp2   = z * L[8];
  double Lstmp3   = z * L[14];
  double Lstmp4   = Lstmp3 * y;
  double Lstmp5   = (x * x);
  double Lstmp6   = (1.0 / 2.0) * Lstmp5;
  double Lstmp7   = (x * x * x);
  double Lstmp8   = (1.0 / 6.0) * Lstmp7;
  double Lstmp9   = (1.0 / 24.0) * (x * x * x * x);
  double Lstmp10  = (y * y);
  double Lstmp11  = (1.0 / 2.0) * Lstmp10;
  double Lstmp12  = (y * y * y);
  double Lstmp13  = (1.0 / 6.0) * Lstmp12;
  double Lstmp14  = (1.0 / 24.0) * (y * y * y * y);
  double Lstmp15  = (z * z);
  double Lstmp16  = (1.0 / 2.0) * Lstmp15;
  double Lstmp17  = (z * z * z);
  double Lstmp18  = (1.0 / 6.0) * Lstmp17;
  double Lstmp19  = (1.0 / 24.0) * (z * z * z * z);
  double Lstmp20  = x * L[13];
  double Lstmp21  = x * L[26];
  double Lstmp22  = x * L[45];
  double Lstmp23  = x * L[15];
  double Lstmp24  = x * L[29];
  double Lstmp25  = x * L[49];
  double Lstmp26  = y * L[11];
  double Lstmp27  = z * L[12];
  double Lstmp28  = y * L[21];
  double Lstmp29  = z * L[22];
  double Lstmp30  = y * L[36];
  double Lstmp31  = z * L[37];
  double Lstmp32  = y * L[18];
  double Lstmp33  = y * L[33];
  double Lstmp34  = y * L[54];
  double Lstmp35  = z * L[17];
  double Lstmp36  = z * L[31];
  double Lstmp37  = z * L[51];
  double Lstmp38  = y * L[28];
  double Lstmp39  = Lstmp38 * x;
  double Lstmp40  = y * L[48];
  double Lstmp41  = Lstmp40 * x;
  double Lstmp42  = z * L[27];
  double Lstmp43  = Lstmp42 * x;
  double Lstmp44  = z * L[46];
  double Lstmp45  = Lstmp44 * x;
  double Lstmp46  = z * L[24];
  double Lstmp47  = Lstmp46 * y;
  double Lstmp48  = z * L[39];
  double Lstmp49  = Lstmp48 * y;
  double Lstmp50  = (1.0 / 4.0) * Lstmp5;
  double Lstmp51  = Lstmp10 * Lstmp50;
  double Lstmp52  = (1.0 / 12.0) * Lstmp5;
  double Lstmp53  = Lstmp15 * Lstmp50;
  double Lstmp54  = (1.0 / 12.0) * Lstmp7;
  double Lstmp55  = (1.0 / 4.0) * Lstmp10 * Lstmp15;
  double Lstmp56  = x * L[47];
  double Lstmp57  = y * L[43];
  double Lstmp58  = z * L[42];
  double Lstmp59  = x * L[23];
  double Lstmp60  = x * L[41];
  double Lstmp61  = x * L[25];
  double Lstmp62  = x * L[44];
  double Lstmp63  = Lstmp57 * x;
  double Lstmp64  = Lstmp58 * x;
  double Lstmp65  = y * L[13];
  double Lstmp66  = Lstmp42 * y;
  double Lstmp67  = x * L[28];
  double Lstmp68  = x * L[48];
  double Lstmp69  = y * L[23];
  double Lstmp70  = y * L[38];
  double Lstmp71  = y * L[32];
  double Lstmp72  = y * L[53];
  double Lstmp73  = y * L[47];
  double Lstmp74  = Lstmp73 * x;
  double Lstmp75  = Lstmp58 * y;
  double Lstmp76  = y * L[14];
  double Lstmp77  = z * L[15];
  double Lstmp78  = z * L[18];
  double Lstmp79  = z * L[28];
  double Lstmp80  = Lstmp79 * y;
  double Lstmp81  = x * L[27];
  double Lstmp82  = x * L[46];
  double Lstmp83  = y * L[24];
  double Lstmp84  = z * L[25];
  double Lstmp85  = y * L[39];
  double Lstmp86  = z * L[40];
  double Lstmp87  = z * L[32];
  double Lstmp88  = z * L[52];
  double Lstmp89  = z * L[47];
  double Lstmp90  = Lstmp89 * x;
  double Lstmp91  = z * L[43];
  double Lstmp92  = Lstmp91 * y;
  double Lstmp93  = x * L[38];
  double Lstmp94  = x * L[40];
  double Lstmp95  = x * L[43];
  double Lstmp96  = x * L[42];
  double Lstmp97  = y * L[26];
  double Lstmp98  = Lstmp44 * y;
  double Lstmp99  = y * L[41];
  double Lstmp100 = y * L[52];
  double Lstmp101 = y * L[27];
  double Lstmp102 = Lstmp89 * y;
  double Lstmp103 = y * L[42];
  double Lstmp104 = z * L[29];
  double Lstmp105 = z * L[33];
  double Lstmp106 = z * L[48];
  double Lstmp107 = Lstmp106 * y;
  double Lstmp108 = z * L[44];
  double Lstmp109 = z * L[53];
  double Lstmp110 = y * L[45];
  double Lstmp111 = y * L[46];
  double Lstmp112 = z * L[49];
  double Lstmp113 = z * L[54];
#pragma omp atomic
  Ls[0] += Lstmp0 * x + Lstmp1 * x + (1.0 / 12.0) * Lstmp10 * Lstmp17 * L[53] +
           Lstmp10 * Lstmp54 * L[38] + Lstmp11 * Lstmp20 + Lstmp11 * Lstmp35 +
           Lstmp11 * Lstmp43 + Lstmp11 * L[7] + (1.0 / 12.0) * Lstmp12 * Lstmp15 * L[52] +
           Lstmp12 * Lstmp52 * L[41] + Lstmp13 * Lstmp21 + Lstmp13 * Lstmp36 +
           Lstmp13 * Lstmp45 + Lstmp13 * L[16] + Lstmp14 * Lstmp22 + Lstmp14 * Lstmp37 +
           Lstmp14 * L[30] + Lstmp15 * Lstmp54 * L[40] + Lstmp16 * Lstmp23 +
           Lstmp16 * Lstmp32 + Lstmp16 * Lstmp39 + Lstmp16 * L[9] +
           Lstmp17 * Lstmp52 * L[44] + Lstmp18 * Lstmp24 + Lstmp18 * Lstmp33 +
           Lstmp18 * Lstmp41 + Lstmp18 * L[19] + Lstmp19 * Lstmp25 + Lstmp19 * Lstmp34 +
           Lstmp19 * L[34] + Lstmp2 * y + Lstmp26 * Lstmp6 + Lstmp27 * Lstmp6 +
           Lstmp28 * Lstmp8 + Lstmp29 * Lstmp8 + Lstmp30 * Lstmp9 + Lstmp31 * Lstmp9 +
           Lstmp4 * x + Lstmp47 * Lstmp6 + Lstmp49 * Lstmp8 + Lstmp51 * Lstmp58 +
           Lstmp51 * L[23] + Lstmp53 * Lstmp57 + Lstmp53 * L[25] + Lstmp55 * Lstmp56 +
           Lstmp55 * L[32] + Lstmp6 * L[4] + Lstmp8 * L[10] + Lstmp9 * L[20] +
           (1.0 / 120.0) * (x * x * x * x * x) * L[35] + x * L[1] +
           (1.0 / 120.0) * (y * y * y * y * y) * L[50] + y * L[2] +
           (1.0 / 120.0) * (z * z * z * z * z) * L[55] + z * L[3] + L[0];
#pragma omp atomic
  Ls[1] += Lstmp0 + Lstmp1 + Lstmp11 * Lstmp42 + Lstmp11 * Lstmp59 + Lstmp11 * Lstmp64 +
           Lstmp11 * L[13] + Lstmp13 * Lstmp44 + Lstmp13 * Lstmp60 + Lstmp13 * L[26] +
           Lstmp14 * L[45] + Lstmp16 * Lstmp38 + Lstmp16 * Lstmp61 + Lstmp16 * Lstmp63 +
           Lstmp16 * L[15] + Lstmp18 * Lstmp40 + Lstmp18 * Lstmp62 + Lstmp18 * L[29] +
           Lstmp19 * L[49] + Lstmp26 * x + Lstmp27 * x + Lstmp28 * Lstmp6 +
           Lstmp29 * Lstmp6 + Lstmp30 * Lstmp8 + Lstmp31 * Lstmp8 + Lstmp4 + Lstmp47 * x +
           Lstmp49 * Lstmp6 + Lstmp51 * L[38] + Lstmp53 * L[40] + Lstmp55 * L[47] +
           Lstmp6 * L[10] + Lstmp8 * L[20] + Lstmp9 * L[35] + x * L[4] + L[1];
#pragma omp atomic
  Ls[2] += Lstmp11 * Lstmp21 + Lstmp11 * Lstmp36 + Lstmp11 * Lstmp45 + Lstmp11 * L[16] +
           Lstmp13 * Lstmp22 + Lstmp13 * Lstmp37 + Lstmp13 * L[30] + Lstmp14 * L[50] +
           Lstmp16 * Lstmp67 + Lstmp16 * Lstmp71 + Lstmp16 * Lstmp74 + Lstmp16 * L[18] +
           Lstmp18 * Lstmp68 + Lstmp18 * Lstmp72 + Lstmp18 * L[33] + Lstmp19 * L[54] +
           Lstmp2 + Lstmp3 * x + Lstmp35 * y + Lstmp46 * Lstmp6 + Lstmp48 * Lstmp8 +
           Lstmp51 * L[41] + Lstmp53 * L[43] + Lstmp55 * L[52] + Lstmp6 * Lstmp69 +
           Lstmp6 * Lstmp75 + Lstmp6 * L[11] + Lstmp65 * x + Lstmp66 * x +
           Lstmp70 * Lstmp8 + Lstmp8 * L[21] + Lstmp9 * L[36] + x * L[5] + y * L[7] +
           L[2];
#pragma omp atomic
  Ls[3] += Lstmp11 * Lstmp81 + Lstmp11 * Lstmp87 + Lstmp11 * Lstmp90 + Lstmp11 * L[17] +
           Lstmp13 * Lstmp82 + Lstmp13 * Lstmp88 + Lstmp13 * L[31] + Lstmp14 * L[51] +
           Lstmp16 * Lstmp24 + Lstmp16 * Lstmp33 + Lstmp16 * Lstmp41 + Lstmp16 * L[19] +
           Lstmp18 * Lstmp25 + Lstmp18 * Lstmp34 + Lstmp18 * L[34] + Lstmp19 * L[55] +
           Lstmp51 * L[42] + Lstmp53 * L[44] + Lstmp55 * L[53] + Lstmp6 * Lstmp83 +
           Lstmp6 * Lstmp84 + Lstmp6 * Lstmp92 + Lstmp6 * L[12] + Lstmp76 * x +
           Lstmp77 * x + Lstmp78 * y + Lstmp8 * Lstmp85 + Lstmp8 * Lstmp86 +
           Lstmp8 * L[22] + Lstmp80 * x + Lstmp9 * L[37] + x * L[6] + y * L[8] +
           z * L[9] + L[3];
#pragma omp atomic
  Ls[4] += Lstmp11 * Lstmp58 + Lstmp11 * Lstmp93 + Lstmp11 * L[23] + Lstmp13 * L[41] +
           Lstmp16 * Lstmp57 + Lstmp16 * Lstmp94 + Lstmp16 * L[25] + Lstmp18 * L[44] +
           Lstmp26 + Lstmp27 + Lstmp28 * x + Lstmp29 * x + Lstmp30 * Lstmp6 +
           Lstmp31 * Lstmp6 + Lstmp47 + Lstmp49 * x + Lstmp6 * L[20] + Lstmp8 * L[35] +
           x * L[10] + L[4];
#pragma omp atomic
  Ls[5] += Lstmp11 * Lstmp44 + Lstmp11 * Lstmp60 + Lstmp11 * L[26] + Lstmp13 * L[45] +
           Lstmp16 * Lstmp73 + Lstmp16 * Lstmp95 + Lstmp16 * L[28] + Lstmp18 * L[48] +
           Lstmp3 + Lstmp46 * x + Lstmp48 * Lstmp6 + Lstmp6 * Lstmp70 + Lstmp6 * L[21] +
           Lstmp65 + Lstmp66 + Lstmp69 * x + Lstmp75 * x + Lstmp8 * L[36] + x * L[11] +
           L[5];
#pragma omp atomic
  Ls[6] += Lstmp11 * Lstmp89 + Lstmp11 * Lstmp96 + Lstmp11 * L[27] + Lstmp13 * L[46] +
           Lstmp16 * Lstmp40 + Lstmp16 * Lstmp62 + Lstmp16 * L[29] + Lstmp18 * L[49] +
           Lstmp6 * Lstmp85 + Lstmp6 * Lstmp86 + Lstmp6 * L[22] + Lstmp76 + Lstmp77 +
           Lstmp8 * L[37] + Lstmp80 + Lstmp83 * x + Lstmp84 * x + Lstmp92 * x +
           x * L[12] + L[6];
#pragma omp atomic
  Ls[7] += Lstmp100 * Lstmp16 + Lstmp11 * Lstmp22 + Lstmp11 * Lstmp37 + Lstmp11 * L[30] +
           Lstmp13 * L[50] + Lstmp16 * Lstmp56 + Lstmp16 * L[32] + Lstmp18 * L[53] +
           Lstmp20 + Lstmp35 + Lstmp36 * y + Lstmp43 + Lstmp58 * Lstmp6 +
           Lstmp6 * Lstmp99 + Lstmp6 * L[23] + Lstmp8 * L[38] + Lstmp97 * x +
           Lstmp98 * x + y * L[16] + L[7];
#pragma omp atomic
  Ls[8] += Lstmp101 * x + Lstmp102 * x + Lstmp103 * Lstmp6 + Lstmp11 * Lstmp82 +
           Lstmp11 * Lstmp88 + Lstmp11 * L[31] + Lstmp13 * L[51] + Lstmp16 * Lstmp68 +
           Lstmp16 * Lstmp72 + Lstmp16 * L[33] + Lstmp18 * L[54] + Lstmp6 * Lstmp91 +
           Lstmp6 * L[24] + Lstmp78 + Lstmp79 * x + Lstmp8 * L[39] + Lstmp87 * y +
           x * L[14] + y * L[17] + L[8];
#pragma omp atomic
  Ls[9] += Lstmp104 * x + Lstmp105 * y + Lstmp107 * x + Lstmp108 * Lstmp6 +
           Lstmp109 * Lstmp11 + Lstmp11 * Lstmp56 + Lstmp11 * L[32] + Lstmp13 * L[52] +
           Lstmp16 * Lstmp25 + Lstmp16 * Lstmp34 + Lstmp16 * L[34] + Lstmp18 * L[55] +
           Lstmp23 + Lstmp32 + Lstmp39 + Lstmp57 * Lstmp6 + Lstmp6 * L[25] +
           Lstmp8 * L[40] + z * L[19] + L[9];
#pragma omp atomic
  Ls[10] += Lstmp11 * L[38] + Lstmp16 * L[40] + Lstmp28 + Lstmp29 + Lstmp30 * x +
            Lstmp31 * x + Lstmp49 + Lstmp6 * L[35] + x * L[20] + L[10];
#pragma omp atomic
  Ls[11] += Lstmp11 * L[41] + Lstmp16 * L[43] + Lstmp46 + Lstmp48 * x + Lstmp6 * L[36] +
            Lstmp69 + Lstmp70 * x + Lstmp75 + x * L[21] + L[11];
#pragma omp atomic
  Ls[12] += Lstmp11 * L[42] + Lstmp16 * L[44] + Lstmp6 * L[37] + Lstmp83 + Lstmp84 +
            Lstmp85 * x + Lstmp86 * x + Lstmp92 + x * L[22] + L[12];
#pragma omp atomic
  Ls[13] += Lstmp11 * L[45] + Lstmp16 * L[47] + Lstmp42 + Lstmp59 + Lstmp6 * L[38] +
            Lstmp64 + Lstmp97 + Lstmp98 + Lstmp99 * x + L[13];
#pragma omp atomic
  Ls[14] += Lstmp101 + Lstmp102 + Lstmp103 * x + Lstmp11 * L[46] + Lstmp16 * L[48] +
            Lstmp6 * L[39] + Lstmp79 + Lstmp91 * x + x * L[24] + L[14];
#pragma omp atomic
  Ls[15] += Lstmp104 + Lstmp107 + Lstmp108 * x + Lstmp11 * L[47] + Lstmp16 * L[49] +
            Lstmp38 + Lstmp6 * L[40] + Lstmp61 + Lstmp63 + L[15];
#pragma omp atomic
  Ls[16] += Lstmp11 * L[50] + Lstmp110 * x + Lstmp16 * L[52] + Lstmp21 + Lstmp36 +
            Lstmp37 * y + Lstmp45 + Lstmp6 * L[41] + y * L[30] + L[16];
#pragma omp atomic
  Ls[17] += Lstmp11 * L[51] + Lstmp111 * x + Lstmp16 * L[53] + Lstmp6 * L[42] + Lstmp81 +
            Lstmp87 + Lstmp88 * y + Lstmp90 + y * L[31] + L[17];
#pragma omp atomic
  Ls[18] += Lstmp105 + Lstmp106 * x + Lstmp109 * y + Lstmp11 * L[52] + Lstmp16 * L[54] +
            Lstmp6 * L[43] + Lstmp67 + Lstmp71 + Lstmp74 + L[18];
#pragma omp atomic
  Ls[19] += Lstmp11 * L[53] + Lstmp112 * x + Lstmp113 * y + Lstmp16 * L[55] + Lstmp24 +
            Lstmp33 + Lstmp41 + Lstmp6 * L[44] + z * L[34] + L[19];
#pragma omp atomic
  Ls[20] += Lstmp30 + Lstmp31 + x * L[35] + L[20];
#pragma omp atomic
  Ls[21] += Lstmp48 + Lstmp70 + x * L[36] + L[21];
#pragma omp atomic
  Ls[22] += Lstmp85 + Lstmp86 + x * L[37] + L[22];
#pragma omp atomic
  Ls[23] += Lstmp58 + Lstmp93 + Lstmp99 + L[23];
#pragma omp atomic
  Ls[24] += Lstmp103 + Lstmp91 + x * L[39] + L[24];
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
  Ls[30] += Lstmp22 + Lstmp37 + y * L[50] + L[30];
#pragma omp atomic
  Ls[31] += Lstmp82 + Lstmp88 + y * L[51] + L[31];
#pragma omp atomic
  Ls[32] += Lstmp100 + Lstmp109 + Lstmp56 + L[32];
#pragma omp atomic
  Ls[33] += Lstmp113 + Lstmp68 + Lstmp72 + L[33];
#pragma omp atomic
  Ls[34] += Lstmp25 + Lstmp34 + z * L[55] + L[34];
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

void field_m1_L2P_6(double x, double y, double z, double* L, double* F) {
  double Ftmp0  = x * y;
  double Ftmp1  = x * z;
  double Ftmp2  = y * z;
  double Ftmp3  = Ftmp0 * z;
  double Ftmp4  = (x * x);
  double Ftmp5  = (1.0 / 2.0) * Ftmp4;
  double Ftmp6  = (1.0 / 6.0) * (x * x * x);
  double Ftmp7  = (1.0 / 24.0) * (x * x * x * x);
  double Ftmp8  = (y * y);
  double Ftmp9  = (1.0 / 2.0) * Ftmp8;
  double Ftmp10 = (1.0 / 6.0) * (y * y * y);
  double Ftmp11 = (1.0 / 24.0) * (y * y * y * y);
  double Ftmp12 = (z * z);
  double Ftmp13 = (1.0 / 2.0) * Ftmp12;
  double Ftmp14 = (1.0 / 6.0) * (z * z * z);
  double Ftmp15 = (1.0 / 24.0) * (z * z * z * z);
  double Ftmp16 = Ftmp9 * x;
  double Ftmp17 = Ftmp10 * x;
  double Ftmp18 = Ftmp13 * x;
  double Ftmp19 = Ftmp14 * x;
  double Ftmp20 = Ftmp5 * y;
  double Ftmp21 = Ftmp5 * z;
  double Ftmp22 = Ftmp6 * y;
  double Ftmp23 = Ftmp6 * z;
  double Ftmp24 = Ftmp13 * y;
  double Ftmp25 = Ftmp14 * y;
  double Ftmp26 = Ftmp9 * z;
  double Ftmp27 = Ftmp10 * z;
  double Ftmp28 = Ftmp0 * Ftmp13;
  double Ftmp29 = Ftmp1 * Ftmp9;
  double Ftmp30 = Ftmp2 * Ftmp5;
  double Ftmp31 = (1.0 / 4.0) * Ftmp4;
  double Ftmp32 = Ftmp31 * Ftmp8;
  double Ftmp33 = Ftmp12 * Ftmp31;
  double Ftmp34 = (1.0 / 4.0) * Ftmp12 * Ftmp8;
#pragma omp atomic
  F[0] += -Ftmp0 * L[11] - Ftmp1 * L[12] - Ftmp10 * L[26] - Ftmp11 * L[45] -
          Ftmp13 * L[15] - Ftmp14 * L[29] - Ftmp15 * L[49] - Ftmp16 * L[23] -
          Ftmp17 * L[41] - Ftmp18 * L[25] - Ftmp19 * L[44] - Ftmp2 * L[14] -
          Ftmp20 * L[21] - Ftmp21 * L[22] - Ftmp22 * L[36] - Ftmp23 * L[37] -
          Ftmp24 * L[28] - Ftmp25 * L[48] - Ftmp26 * L[27] - Ftmp27 * L[46] -
          Ftmp28 * L[43] - Ftmp29 * L[42] - Ftmp3 * L[24] - Ftmp30 * L[39] -
          Ftmp32 * L[38] - Ftmp33 * L[40] - Ftmp34 * L[47] - Ftmp5 * L[10] -
          Ftmp6 * L[20] - Ftmp7 * L[35] - Ftmp9 * L[13] - x * L[4] - y * L[5] - z * L[6] -
          L[1];
#pragma omp atomic
  F[1] += -Ftmp0 * L[13] - Ftmp1 * L[14] - Ftmp10 * L[30] - Ftmp11 * L[50] -
          Ftmp13 * L[18] - Ftmp14 * L[33] - Ftmp15 * L[54] - Ftmp16 * L[26] -
          Ftmp17 * L[45] - Ftmp18 * L[28] - Ftmp19 * L[48] - Ftmp2 * L[17] -
          Ftmp20 * L[23] - Ftmp21 * L[24] - Ftmp22 * L[38] - Ftmp23 * L[39] -
          Ftmp24 * L[32] - Ftmp25 * L[53] - Ftmp26 * L[31] - Ftmp27 * L[51] -
          Ftmp28 * L[47] - Ftmp29 * L[46] - Ftmp3 * L[27] - Ftmp30 * L[42] -
          Ftmp32 * L[41] - Ftmp33 * L[43] - Ftmp34 * L[52] - Ftmp5 * L[11] -
          Ftmp6 * L[21] - Ftmp7 * L[36] - Ftmp9 * L[16] - x * L[5] - y * L[7] - z * L[8] -
          L[2];
#pragma omp atomic
  F[2] += -Ftmp0 * L[14] - Ftmp1 * L[15] - Ftmp10 * L[31] - Ftmp11 * L[51] -
          Ftmp13 * L[19] - Ftmp14 * L[34] - Ftmp15 * L[55] - Ftmp16 * L[27] -
          Ftmp17 * L[46] - Ftmp18 * L[29] - Ftmp19 * L[49] - Ftmp2 * L[18] -
          Ftmp20 * L[24] - Ftmp21 * L[25] - Ftmp22 * L[39] - Ftmp23 * L[40] -
          Ftmp24 * L[33] - Ftmp25 * L[54] - Ftmp26 * L[32] - Ftmp27 * L[52] -
          Ftmp28 * L[48] - Ftmp29 * L[47] - Ftmp3 * L[28] - Ftmp30 * L[43] -
          Ftmp32 * L[42] - Ftmp33 * L[44] - Ftmp34 * L[53] - Ftmp5 * L[12] -
          Ftmp6 * L[22] - Ftmp7 * L[37] - Ftmp9 * L[17] - x * L[6] - y * L[8] - z * L[9] -
          L[3];
}

void field_m1_M2P_6(double x, double y, double z, double* M, double* F) {
  double R       = sqrt(x * x + y * y + z * z);
  double Ftmp0   = 1.0 * pow(R, -3.0);
  double Ftmp1   = pow(R, -5.0);
  double Ftmp2   = 3.0 * Ftmp1;
  double Ftmp3   = Ftmp2 * M[4];
  double Ftmp4   = Ftmp2 * z;
  double Ftmp5   = y * z;
  double Ftmp6   = pow(R, -7.0);
  double Ftmp7   = 15.0 * Ftmp6;
  double Ftmp8   = Ftmp7 * M[13];
  double Ftmp9   = x * y;
  double Ftmp10  = Ftmp2 * M[1];
  double Ftmp11  = Ftmp4 * M[2];
  double Ftmp12  = (x * x);
  double Ftmp13  = Ftmp2 * M[0];
  double Ftmp14  = Ftmp5 * x;
  double Ftmp15  = Ftmp14 * Ftmp7;
  double Ftmp16  = Ftmp12 * Ftmp7;
  double Ftmp17  = pow(R, -9.0);
  double Ftmp18  = 105.0 * Ftmp17;
  double Ftmp19  = Ftmp12 * Ftmp18;
  double Ftmp20  = -9.0 * Ftmp1;
  double Ftmp21  = Ftmp16 + Ftmp20;
  double Ftmp22  = -Ftmp2;
  double Ftmp23  = (y * y);
  double Ftmp24  = Ftmp23 * Ftmp7;
  double Ftmp25  = Ftmp22 + Ftmp24;
  double Ftmp26  = (z * z);
  double Ftmp27  = Ftmp26 * Ftmp7;
  double Ftmp28  = Ftmp22 + Ftmp27;
  double Ftmp29  = Ftmp25 * M[6];
  double Ftmp30  = Ftmp28 * M[8];
  double Ftmp31  = 45.0 * Ftmp6;
  double Ftmp32  = -Ftmp31;
  double Ftmp33  = Ftmp19 + Ftmp32;
  double Ftmp34  = Ftmp33 * M[20];
  double Ftmp35  = 1.0 * y;
  double Ftmp36  = Ftmp18 * Ftmp23;
  double Ftmp37  = Ftmp32 + Ftmp36;
  double Ftmp38  = Ftmp37 * M[25];
  double Ftmp39  = 3.0 * y;
  double Ftmp40  = 35.0 * Ftmp17;
  double Ftmp41  = (Ftmp26 * Ftmp40 - 5.0 * Ftmp6) * M[27];
  double Ftmp42  = Ftmp33 * M[21];
  double Ftmp43  = 1.0 * z;
  double Ftmp44  = -Ftmp7;
  double Ftmp45  = Ftmp36 + Ftmp44;
  double Ftmp46  = Ftmp45 * M[26];
  double Ftmp47  = Ftmp18 * Ftmp26;
  double Ftmp48  = Ftmp32 + Ftmp47;
  double Ftmp49  = Ftmp43 * Ftmp48;
  double Ftmp50  = 315.0 * Ftmp17;
  double Ftmp51  = -Ftmp50;
  double Ftmp52  = pow(R, -11.0);
  double Ftmp53  = 945.0 * Ftmp52;
  double Ftmp54  = Ftmp12 * Ftmp53;
  double Ftmp55  = Ftmp51 + Ftmp54;
  double Ftmp56  = Ftmp55 * M[38];
  double Ftmp57  = Ftmp33 * Ftmp9;
  double Ftmp58  = Ftmp37 * M[15];
  double Ftmp59  = Ftmp44 + Ftmp47;
  double Ftmp60  = Ftmp59 * M[17];
  double Ftmp61  = Ftmp35 * x;
  double Ftmp62  = x * z;
  double Ftmp63  = Ftmp33 * Ftmp62;
  double Ftmp64  = Ftmp45 * M[16];
  double Ftmp65  = Ftmp48 * M[18];
  double Ftmp66  = Ftmp23 * Ftmp53;
  double Ftmp67  = Ftmp51 + Ftmp66;
  double Ftmp68  = Ftmp35 * Ftmp67;
  double Ftmp69  = Ftmp68 * M[45];
  double Ftmp70  = -Ftmp18;
  double Ftmp71  = Ftmp26 * Ftmp52;
  double Ftmp72  = 315.0 * Ftmp71;
  double Ftmp73  = Ftmp70 + Ftmp72;
  double Ftmp74  = Ftmp73 * M[47];
  double Ftmp75  = Ftmp39 * Ftmp74;
  double Ftmp76  = Ftmp67 * M[30];
  double Ftmp77  = -75.0 * Ftmp6;
  double Ftmp78  = 1.0 * Ftmp12;
  double Ftmp79  = Ftmp45 * M[12];
  double Ftmp80  = Ftmp59 * M[14];
  double Ftmp81  = 525.0 * Ftmp17;
  double Ftmp82  = -Ftmp81;
  double Ftmp83  = Ftmp54 + Ftmp82;
  double Ftmp84  = Ftmp12 * y;
  double Ftmp85  = Ftmp12 * z;
  double Ftmp86  = Ftmp26 * Ftmp53;
  double Ftmp87  = Ftmp51 + Ftmp86;
  double Ftmp88  = Ftmp35 * M[32];
  double Ftmp89  = Ftmp68 * M[25];
  double Ftmp90  = Ftmp12 * Ftmp39;
  double Ftmp91  = (-Ftmp40 + Ftmp72) * M[27];
  double Ftmp92  = Ftmp12 * Ftmp43;
  double Ftmp93  = (Ftmp66 + Ftmp70) * M[26];
  double Ftmp94  = Ftmp87 * M[28];
  double Ftmp95  = 4725.0 * Ftmp52;
  double Ftmp96  = -Ftmp95;
  double Ftmp97  = pow(R, -13.0);
  double Ftmp98  = 10395.0 * Ftmp97;
  double Ftmp99  = Ftmp12 * Ftmp98;
  double Ftmp100 = 2835.0 * Ftmp52;
  double Ftmp101 = -Ftmp100;
  double Ftmp102 = Ftmp23 * Ftmp98;
  double Ftmp103 = Ftmp35 * (Ftmp101 + Ftmp102) * M[45];
  double Ftmp104 = 3465.0 * Ftmp97;
  double Ftmp105 = Ftmp104 * Ftmp26;
  double Ftmp106 = (Ftmp105 - Ftmp53) * M[47];
  double Ftmp107 = 225.0 * Ftmp6;
  double Ftmp108 = (x * x * x * x);
  double Ftmp109 = Ftmp108 * Ftmp53;
  double Ftmp110 = 1050.0 * Ftmp17;
  double Ftmp111 = Ftmp107 + Ftmp109 - Ftmp110 * Ftmp12;
  double Ftmp112 = (y * y * y * y);
  double Ftmp113 = Ftmp112 * Ftmp53;
  double Ftmp114 = 630.0 * Ftmp17;
  double Ftmp115 = Ftmp113 - Ftmp114 * Ftmp23 + Ftmp31;
  double Ftmp116 = (z * z * z * z);
  double Ftmp117 = Ftmp116 * Ftmp53;
  double Ftmp118 = -Ftmp114 * Ftmp26 + Ftmp117 + Ftmp31;
  double Ftmp119 = Ftmp115 * M[29];
  double Ftmp120 = Ftmp118 * M[33];
  double Ftmp121 = 1575.0 * Ftmp17;
  double Ftmp122 = Ftmp108 * Ftmp97;
  double Ftmp123 = 10395.0 * Ftmp122;
  double Ftmp124 = Ftmp12 * Ftmp52;
  double Ftmp125 = 9450.0 * Ftmp124;
  double Ftmp126 = Ftmp121 + Ftmp123 - Ftmp125;
  double Ftmp127 = Ftmp126 * M[56];
  double Ftmp128 = Ftmp112 * Ftmp98;
  double Ftmp129 = 9450.0 * Ftmp52;
  double Ftmp130 = Ftmp129 * Ftmp23;
  double Ftmp131 = Ftmp121 + Ftmp128 - Ftmp130;
  double Ftmp132 = Ftmp131 * M[70];
  double Ftmp133 = (Ftmp104 * Ftmp116 + Ftmp18 - 1890.0 * Ftmp71) * M[74];
  double Ftmp134 = Ftmp126 * M[57];
  double Ftmp135 = 5670.0 * Ftmp52;
  double Ftmp136 = Ftmp135 * Ftmp23;
  double Ftmp137 = Ftmp128 - Ftmp136 + Ftmp50;
  double Ftmp138 = Ftmp137 * M[71];
  double Ftmp139 = Ftmp116 * Ftmp98;
  double Ftmp140 = Ftmp129 * Ftmp26;
  double Ftmp141 = Ftmp121 + Ftmp139 - Ftmp140;
  double Ftmp142 = Ftmp141 * Ftmp43;
  double Ftmp143 = Ftmp126 * Ftmp9;
  double Ftmp144 = Ftmp131 * M[49];
  double Ftmp145 = Ftmp135 * Ftmp26;
  double Ftmp146 = Ftmp139 - Ftmp145 + Ftmp50;
  double Ftmp147 = Ftmp146 * M[53];
  double Ftmp148 = Ftmp126 * Ftmp62;
  double Ftmp149 = Ftmp137 * M[50];
  double Ftmp150 = Ftmp141 * M[54];
  double Ftmp151 = 14175.0 * Ftmp52;
  double Ftmp152 = pow(R, -15.0);
  double Ftmp153 = 135135.0 * Ftmp152;
  double Ftmp154 = Ftmp108 * Ftmp153;
  double Ftmp155 = Ftmp12 * Ftmp97;
  double Ftmp156 = 103950.0 * Ftmp155;
  double Ftmp157 = Ftmp151 + Ftmp154 - Ftmp156;
  double Ftmp158 = Ftmp112 * Ftmp153;
  double Ftmp159 = Ftmp23 * Ftmp97;
  double Ftmp160 = 103950.0 * Ftmp159;
  double Ftmp161 = Ftmp151 + Ftmp158 - Ftmp160;
  double Ftmp162 = Ftmp161 * M[77];
  double Ftmp163 = 3675.0 * Ftmp17;
  double Ftmp164 = Ftmp137 * M[44];
  double Ftmp165 = Ftmp146 * M[48];
  double Ftmp166 = 33075.0 * Ftmp52;
  double Ftmp167 = Ftmp154 - 145530.0 * Ftmp155 + Ftmp166;
  double Ftmp168 = Ftmp116 * Ftmp153;
  double Ftmp169 = Ftmp26 * Ftmp97;
  double Ftmp170 = Ftmp151 + Ftmp168 - 103950.0 * Ftmp169;
  double Ftmp171 = Ftmp35 * M[81];
  double Ftmp172 = Ftmp12 * Ftmp35;
  double Ftmp173 = Ftmp161 * M[70];
  double Ftmp174 = 45045.0 * Ftmp116 * Ftmp152;
  double Ftmp175 = (-20790.0 * Ftmp169 + Ftmp174 + Ftmp53) * M[74];
  double Ftmp176 = 62370.0 * Ftmp159;
  double Ftmp177 = (Ftmp100 + Ftmp158 - Ftmp176) * M[71];
  double Ftmp178 = Ftmp170 * M[75];
  double Ftmp179 = -11025.0 * Ftmp17;
  double Ftmp180 = Ftmp153 * (x * x * x * x * x * x);
  double Ftmp181 = -Ftmp121;
  double Ftmp182 = Ftmp153 * (y * y * y * y * y * y);
  double Ftmp183 = 155925.0 * Ftmp97;
  double Ftmp184 = 42525.0 * Ftmp52;
  double Ftmp185 = (-Ftmp112 * Ftmp183 + Ftmp181 + Ftmp182 + Ftmp184 * Ftmp23) * M[76];
  double Ftmp186 = Ftmp153 * (z * z * z * z * z * z);
  double Ftmp187 = (-Ftmp116 * Ftmp183 + Ftmp181 + Ftmp184 * Ftmp26 + Ftmp186) * M[82];
  double Ftmp188 = -Ftmp23 * Ftmp50;
  double Ftmp189 = Ftmp23 * Ftmp54;
  double Ftmp190 = -Ftmp19;
  double Ftmp191 = Ftmp190 + Ftmp31;
  double Ftmp192 = Ftmp188 + Ftmp189 + Ftmp191;
  double Ftmp193 = -Ftmp26 * Ftmp50;
  double Ftmp194 = Ftmp26 * Ftmp54;
  double Ftmp195 = Ftmp191 + Ftmp193 + Ftmp194;
  double Ftmp196 = -Ftmp47;
  double Ftmp197 = Ftmp196 + Ftmp7;
  double Ftmp198 = -Ftmp36;
  double Ftmp199 = Ftmp26 * Ftmp66;
  double Ftmp200 = Ftmp198 + Ftmp199;
  double Ftmp201 = Ftmp197 + Ftmp200;
  double Ftmp202 = Ftmp100 * Ftmp23;
  double Ftmp203 = -Ftmp202;
  double Ftmp204 = Ftmp23 * Ftmp99;
  double Ftmp205 = Ftmp203 + Ftmp204;
  double Ftmp206 = 945.0 * Ftmp17;
  double Ftmp207 = Ftmp100 * Ftmp12;
  double Ftmp208 = -Ftmp207;
  double Ftmp209 = Ftmp206 + Ftmp208;
  double Ftmp210 = Ftmp205 + Ftmp209;
  double Ftmp211 = Ftmp210 * M[61];
  double Ftmp212 = -Ftmp54;
  double Ftmp213 = Ftmp212 + Ftmp50;
  double Ftmp214 = Ftmp100 * Ftmp26;
  double Ftmp215 = -Ftmp214;
  double Ftmp216 = Ftmp26 * Ftmp99;
  double Ftmp217 = Ftmp215 + Ftmp216;
  double Ftmp218 = Ftmp213 + Ftmp217;
  double Ftmp219 = Ftmp218 * M[63];
  double Ftmp220 = -Ftmp66;
  double Ftmp221 = Ftmp102 * Ftmp26;
  double Ftmp222 = Ftmp221 + Ftmp50;
  double Ftmp223 = Ftmp215 + Ftmp220 + Ftmp222;
  double Ftmp224 = Ftmp223 * M[72];
  double Ftmp225 = Ftmp205 + Ftmp213;
  double Ftmp226 = Ftmp225 * M[62];
  double Ftmp227 = Ftmp209 + Ftmp217;
  double Ftmp228 = Ftmp227 * M[64];
  double Ftmp229 = -Ftmp86;
  double Ftmp230 = Ftmp203 + Ftmp222 + Ftmp229;
  double Ftmp231 = Ftmp230 * M[73];
  double Ftmp232 = Ftmp210 * Ftmp9;
  double Ftmp233 = Ftmp218 * Ftmp9;
  double Ftmp234 = Ftmp225 * Ftmp62;
  double Ftmp235 = Ftmp227 * Ftmp62;
  double Ftmp236 = 31185.0 * Ftmp97;
  double Ftmp237 = Ftmp12 * Ftmp236;
  double Ftmp238 = -Ftmp237;
  double Ftmp239 = 8505.0 * Ftmp52;
  double Ftmp240 = Ftmp238 + Ftmp239;
  double Ftmp241 = Ftmp12 * Ftmp153;
  double Ftmp242 = Ftmp23 * Ftmp241;
  double Ftmp243 = -Ftmp23 * Ftmp236;
  double Ftmp244 = Ftmp242 + Ftmp243;
  double Ftmp245 = Ftmp14 * (Ftmp240 + Ftmp244);
  double Ftmp246 = Ftmp241 * Ftmp26;
  double Ftmp247 = Ftmp236 * Ftmp26;
  double Ftmp248 = -Ftmp247;
  double Ftmp249 = Ftmp246 + Ftmp248;
  double Ftmp250 = Ftmp14 * (Ftmp240 + Ftmp249);
  double Ftmp251 = Ftmp153 * Ftmp23 * Ftmp26;
  double Ftmp252 = Ftmp243 + Ftmp251;
  double Ftmp253 = Ftmp239 + Ftmp248 + Ftmp252;
  double Ftmp254 = -Ftmp23 * Ftmp95;
  double Ftmp255 = Ftmp212 + Ftmp81;
  double Ftmp256 = -Ftmp26 * Ftmp95;
  double Ftmp257 = Ftmp18 + Ftmp229;
  double Ftmp258 = Ftmp220 + Ftmp221;
  double Ftmp259 = 51975.0 * Ftmp97;
  double Ftmp260 = -Ftmp23 * Ftmp259;
  double Ftmp261 = Ftmp242 + Ftmp260;
  double Ftmp262 = Ftmp151 + Ftmp238;
  double Ftmp263 = -Ftmp99;
  double Ftmp264 = Ftmp263 + Ftmp95;
  double Ftmp265 = -Ftmp259 * Ftmp26;
  double Ftmp266 = Ftmp246 + Ftmp265;
  double Ftmp267 = -Ftmp102;
  double Ftmp268 = Ftmp100 + Ftmp251;
  double Ftmp269 = -Ftmp26 * Ftmp98;
  double Ftmp270 = -Ftmp206;
  double Ftmp271 = Ftmp207 + Ftmp270;
  double Ftmp272 = Ftmp12 * Ftmp158;
  double Ftmp273 = 62370.0 * Ftmp155;
  double Ftmp274 = -Ftmp23 * Ftmp273;
  double Ftmp275 = Ftmp272 + Ftmp274;
  double Ftmp276 = 17010.0 * Ftmp52;
  double Ftmp277 = -Ftmp112 * Ftmp236 + Ftmp23 * Ftmp276;
  double Ftmp278 = Ftmp12 * Ftmp168;
  double Ftmp279 = -Ftmp26 * Ftmp273;
  double Ftmp280 = Ftmp278 + Ftmp279;
  double Ftmp281 = -Ftmp116 * Ftmp236 + Ftmp26 * Ftmp276;
  double Ftmp282 = Ftmp151 * Ftmp23;
  double Ftmp283 = -Ftmp156 * Ftmp23 + Ftmp181;
  double Ftmp284 = -Ftmp123;
  double Ftmp285 = Ftmp154 * Ftmp23;
  double Ftmp286 = Ftmp284 + Ftmp285;
  double Ftmp287 = -Ftmp156 * Ftmp26;
  double Ftmp288 = Ftmp154 * Ftmp26;
  double Ftmp289 = Ftmp284 + Ftmp288;
  double Ftmp290 = Ftmp151 * Ftmp26 + Ftmp181;
  double Ftmp291 = -Ftmp176 * Ftmp26;
  double Ftmp292 = Ftmp291 + Ftmp51;
  double Ftmp293 = -Ftmp139;
  double Ftmp294 = Ftmp168 * Ftmp23;
  double Ftmp295 = Ftmp293 + Ftmp294;
  double Ftmp296 = -Ftmp128;
  double Ftmp297 = Ftmp158 * Ftmp26;
  double Ftmp298 = Ftmp296 + Ftmp297;
  double Ftmp299 = Ftmp23 * Ftmp246;
  double Ftmp300 = -Ftmp204 + Ftmp214 + Ftmp299;
  double Ftmp301 = Ftmp202 - Ftmp216;
  double Ftmp302 = x * M[5];
  double Ftmp303 = Ftmp16 + Ftmp22;
  double Ftmp304 = Ftmp20 + Ftmp24;
  double Ftmp305 = Ftmp303 * M[3];
  double Ftmp306 = 1.0 * x;
  double Ftmp307 = 3.0 * x;
  double Ftmp308 = Ftmp19 + Ftmp44;
  double Ftmp309 = Ftmp308 * M[23];
  double Ftmp310 = Ftmp37 * M[30];
  double Ftmp311 = Ftmp43 * x;
  double Ftmp312 = 3.0 * Ftmp62;
  double Ftmp313 = Ftmp308 * M[11];
  double Ftmp314 = Ftmp55 * M[21];
  double Ftmp315 = Ftmp308 * M[10];
  double Ftmp316 = 1.0 * Ftmp23;
  double Ftmp317 = Ftmp23 * x;
  double Ftmp318 = Ftmp55 * M[20];
  double Ftmp319 = Ftmp23 * z;
  double Ftmp320 = (Ftmp54 + Ftmp70) * M[23];
  double Ftmp321 = Ftmp66 + Ftmp82;
  double Ftmp322 = Ftmp35 * Ftmp62;
  double Ftmp323 = Ftmp23 * Ftmp306;
  double Ftmp324 = Ftmp23 * Ftmp307;
  double Ftmp325 = Ftmp23 * Ftmp43;
  double Ftmp326 = (Ftmp101 + Ftmp99) * M[38];
  double Ftmp327 = Ftmp109 - Ftmp114 * Ftmp12 + Ftmp31;
  double Ftmp328 = Ftmp107 - Ftmp110 * Ftmp23 + Ftmp113;
  double Ftmp329 = Ftmp327 * M[19];
  double Ftmp330 = 5670.0 * Ftmp124;
  double Ftmp331 = Ftmp123 - Ftmp330 + Ftmp50;
  double Ftmp332 = Ftmp331 * M[59];
  double Ftmp333 = Ftmp131 * M[77];
  double Ftmp334 = Ftmp331 * M[36];
  double Ftmp335 = Ftmp157 * M[57];
  double Ftmp336 = Ftmp331 * M[35];
  double Ftmp337 = Ftmp23 * Ftmp52;
  double Ftmp338 = Ftmp157 * M[56];
  double Ftmp339 = (Ftmp100 + Ftmp154 - Ftmp273) * M[59];
  double Ftmp340 = Ftmp158 - 145530.0 * Ftmp159 + Ftmp166;
  double Ftmp341 = (-Ftmp108 * Ftmp183 + Ftmp12 * Ftmp184 + Ftmp180 + Ftmp181) * M[55];
  double Ftmp342 = 218295.0 * Ftmp97;
  double Ftmp343 = Ftmp189 + Ftmp198;
  double Ftmp344 = -Ftmp12 * Ftmp50 + Ftmp31;
  double Ftmp345 = Ftmp343 + Ftmp344;
  double Ftmp346 = Ftmp190 + Ftmp194 + Ftmp197;
  double Ftmp347 = Ftmp193 + Ftmp200 + Ftmp31;
  double Ftmp348 = Ftmp208 + Ftmp50;
  double Ftmp349 = Ftmp204 + Ftmp220;
  double Ftmp350 = Ftmp348 + Ftmp349;
  double Ftmp351 = Ftmp350 * M[66];
  double Ftmp352 = Ftmp216 + Ftmp229;
  double Ftmp353 = Ftmp348 + Ftmp352;
  double Ftmp354 = Ftmp353 * M[68];
  double Ftmp355 = Ftmp203 + Ftmp206 + Ftmp215 + Ftmp221;
  double Ftmp356 = Ftmp355 * M[79];
  double Ftmp357 = Ftmp350 * Ftmp5;
  double Ftmp358 = Ftmp353 * Ftmp5;
  double Ftmp359 = Ftmp355 * Ftmp5;
  double Ftmp360 = -Ftmp12 * Ftmp95 + Ftmp81;
  double Ftmp361 = -51975.0 * Ftmp155;
  double Ftmp362 = Ftmp151 + Ftmp361;
  double Ftmp363 = Ftmp100 + Ftmp263;
  double Ftmp364 = Ftmp267 + Ftmp95;
  double Ftmp365 = Ftmp100 + Ftmp238;
  double Ftmp366 = Ftmp246 + Ftmp269;
  double Ftmp367 = Ftmp253 * Ftmp322;
  double Ftmp368 = Ftmp12 * Ftmp151;
  double Ftmp369 = Ftmp207 + Ftmp51;
  double Ftmp370 = Ftmp202 + Ftmp270;
  double Ftmp371 = -31185.0 * Ftmp122 + 17010.0 * Ftmp124;
  double Ftmp372 = Ftmp330 + Ftmp51;
  double Ftmp373 = Ftmp214 + Ftmp279;
  double Ftmp374 = -Ftmp160 * Ftmp26;
  double Ftmp375 = Ftmp207 - Ftmp221;
  double Ftmp376 = y * M[7];
  double Ftmp377 = Ftmp20 + Ftmp27;
  double Ftmp378 = Ftmp306 * M[28];
  double Ftmp379 = Ftmp35 * z;
  double Ftmp380 = Ftmp26 * x;
  double Ftmp381 = Ftmp26 * y;
  double Ftmp382 = Ftmp39 * Ftmp62;
  double Ftmp383 = Ftmp26 * Ftmp306;
  double Ftmp384 = Ftmp26 * (Ftmp82 + Ftmp86);
  double Ftmp385 = Ftmp107 - Ftmp110 * Ftmp26 + Ftmp117;
  double Ftmp386 = Ftmp306 * M[75];
  double Ftmp387 = Ftmp26 * (Ftmp166 + Ftmp168 - 145530.0 * Ftmp169);
  double Ftmp388 = Ftmp190 + Ftmp343 + Ftmp7;
  double Ftmp389 = Ftmp194 + Ftmp196 + Ftmp344;
  double Ftmp390 = Ftmp188 + Ftmp196 + Ftmp199 + Ftmp31;
  double Ftmp391 = Ftmp251 + Ftmp260;
  double Ftmp392 = Ftmp140 + Ftmp181;
#pragma omp atomic
  F[0] += Ftmp0 * M[0] - Ftmp10 * Ftmp9 - Ftmp103 * Ftmp85 - Ftmp106 * Ftmp39 * Ftmp85 -
          Ftmp11 * x + Ftmp111 * x * M[19] + Ftmp111 * M[34] + Ftmp115 * M[44] +
          Ftmp118 * M[48] + Ftmp119 * x - Ftmp12 * Ftmp13 -
          Ftmp12 * Ftmp5 * (Ftmp96 + Ftmp99) * M[38] + Ftmp12 * Ftmp89 -
          Ftmp12 * (Ftmp19 + Ftmp77) * M[9] -
          Ftmp12 * (Ftmp123 - 13230.0 * Ftmp124 + Ftmp163) * M[34] -
          Ftmp12 * (Ftmp204 + Ftmp254 + Ftmp255) * M[37] -
          Ftmp12 * (Ftmp216 + Ftmp255 + Ftmp256) * M[39] + Ftmp120 * x - Ftmp127 * y -
          Ftmp132 * Ftmp35 - Ftmp133 * Ftmp39 - Ftmp134 * z - Ftmp138 * Ftmp43 +
          Ftmp14 * Ftmp157 * M[59] + Ftmp14 * Ftmp162 + Ftmp14 * Ftmp253 * M[79] +
          Ftmp14 * Ftmp55 * M[23] + Ftmp14 * Ftmp76 - Ftmp142 * M[75] - Ftmp143 * M[35] -
          Ftmp144 * Ftmp9 - Ftmp147 * Ftmp61 - Ftmp148 * M[36] - Ftmp149 * Ftmp62 +
          Ftmp15 * M[7] - Ftmp150 * Ftmp62 + Ftmp16 * y * M[4] + Ftmp16 * z * M[5] -
          Ftmp164 * Ftmp78 - Ftmp165 * Ftmp78 + Ftmp167 * Ftmp84 * M[56] +
          Ftmp167 * Ftmp85 * M[57] + Ftmp170 * Ftmp171 * Ftmp62 + Ftmp172 * Ftmp173 +
          Ftmp172 * (Ftmp248 + Ftmp267 + Ftmp268) * M[72] + Ftmp175 * Ftmp90 +
          Ftmp177 * Ftmp92 + Ftmp178 * Ftmp92 + Ftmp185 * x + Ftmp187 * x -
          Ftmp19 * Ftmp5 * M[13] + Ftmp192 * x * M[22] + Ftmp192 * M[37] +
          Ftmp195 * x * M[24] + Ftmp195 * M[39] + Ftmp201 * x * M[31] + Ftmp201 * M[46] +
          Ftmp21 * x * M[3] + Ftmp21 * M[9] - Ftmp211 * y - Ftmp219 * y -
          Ftmp223 * Ftmp9 * M[51] - Ftmp224 * Ftmp35 - Ftmp226 * z - Ftmp228 * z -
          Ftmp230 * Ftmp62 * M[52] - Ftmp231 * Ftmp43 - Ftmp232 * M[40] -
          Ftmp233 * M[42] - Ftmp234 * M[41] - Ftmp235 * M[43] + Ftmp245 * M[66] +
          Ftmp25 * M[12] + Ftmp250 * M[68] + Ftmp28 * M[14] + Ftmp29 * x - Ftmp3 * y +
          Ftmp30 * x - Ftmp34 * y - Ftmp35 * Ftmp38 - Ftmp39 * Ftmp41 - Ftmp4 * M[5] -
          Ftmp42 * z - Ftmp43 * Ftmp46 - Ftmp49 * M[28] + Ftmp5 * Ftmp56 + Ftmp5 * Ftmp8 -
          Ftmp57 * M[10] - Ftmp58 * Ftmp9 - Ftmp60 * Ftmp61 - Ftmp62 * Ftmp64 -
          Ftmp62 * Ftmp65 + Ftmp62 * Ftmp87 * Ftmp88 - Ftmp63 * M[11] + Ftmp69 * z +
          Ftmp75 * z - Ftmp78 * Ftmp79 - Ftmp78 * Ftmp80 -
          Ftmp78 * (Ftmp257 + Ftmp258) * M[46] + Ftmp83 * Ftmp84 * M[20] +
          Ftmp83 * Ftmp85 * M[21] + Ftmp84 * (Ftmp261 + Ftmp262) * M[61] +
          Ftmp84 * (Ftmp264 + Ftmp266) * M[63] + Ftmp85 * (Ftmp261 + Ftmp264) * M[62] +
          Ftmp85 * (Ftmp262 + Ftmp266) * M[64] + Ftmp90 * Ftmp91 + Ftmp92 * Ftmp93 +
          Ftmp92 * Ftmp94 + Ftmp92 * (Ftmp243 + Ftmp268 + Ftmp269) * M[73] +
          x * (Ftmp271 + Ftmp275 + Ftmp277) * M[65] +
          x * (Ftmp271 + Ftmp280 + Ftmp281) * M[69] +
          x * (-218295.0 * Ftmp122 + 99225.0 * Ftmp124 + Ftmp179 + Ftmp180) * M[55] +
          x * (Ftmp125 + Ftmp282 + Ftmp283 + Ftmp286) * M[58] +
          x * (Ftmp125 + Ftmp287 + Ftmp289 + Ftmp290) * M[60] +
          x * (Ftmp136 + Ftmp214 + Ftmp292 + Ftmp298) * M[78] +
          x * (Ftmp145 + Ftmp202 + Ftmp292 + Ftmp295) * M[80] +
          x * (-Ftmp23 * Ftmp247 + Ftmp300 + Ftmp301 + Ftmp55) * M[67];
#pragma omp atomic
  F[1] += Ftmp0 * M[1] - Ftmp10 * Ftmp23 - Ftmp106 * Ftmp23 * Ftmp312 - Ftmp11 * y +
          Ftmp118 * M[53] + Ftmp120 * y - Ftmp127 * x - Ftmp13 * Ftmp9 -
          Ftmp131 * Ftmp5 * M[50] - Ftmp131 * Ftmp61 * M[44] - Ftmp132 * Ftmp306 -
          Ftmp133 * Ftmp307 + Ftmp14 * Ftmp314 + Ftmp14 * Ftmp335 - Ftmp142 * M[81] -
          Ftmp143 * M[34] - Ftmp147 * Ftmp316 - Ftmp150 * Ftmp5 +
          Ftmp161 * Ftmp322 * M[71] - Ftmp165 * Ftmp61 + Ftmp170 * Ftmp325 * M[81] +
          Ftmp175 * Ftmp324 + Ftmp178 * Ftmp322 + Ftmp187 * y - Ftmp211 * x -
          Ftmp219 * x - Ftmp223 * Ftmp61 * M[46] - Ftmp224 * Ftmp306 - Ftmp23 * Ftmp315 -
          Ftmp23 * Ftmp326 * Ftmp62 - Ftmp23 * Ftmp336 -
          Ftmp23 * (Ftmp349 + Ftmp360) * M[40] - Ftmp23 * (Ftmp36 + Ftmp77) * M[15] -
          Ftmp23 * (Ftmp128 + Ftmp163 - 13230.0 * Ftmp337) * M[49] -
          Ftmp23 * (Ftmp212 + Ftmp216 + Ftmp257) * M[42] -
          Ftmp23 * (Ftmp256 + Ftmp258 + Ftmp81) * M[51] - Ftmp232 * M[37] -
          Ftmp233 * M[39] + Ftmp24 * x * M[4] + Ftmp24 * z * M[7] + Ftmp245 * M[62] +
          Ftmp250 * M[64] + Ftmp28 * M[17] - Ftmp3 * x + Ftmp30 * y +
          Ftmp302 * Ftmp5 * Ftmp7 + Ftmp303 * M[10] + Ftmp304 * y * M[6] +
          Ftmp304 * M[15] + Ftmp305 * y - Ftmp306 * Ftmp38 - Ftmp307 * Ftmp41 -
          Ftmp309 * z - Ftmp310 * z + Ftmp311 * Ftmp67 * M[45] + Ftmp312 * Ftmp74 -
          Ftmp313 * Ftmp5 - Ftmp316 * Ftmp60 + Ftmp317 * Ftmp318 + Ftmp317 * Ftmp338 -
          Ftmp317 * Ftmp43 * (Ftmp102 + Ftmp96) * M[45] +
          Ftmp317 * (Ftmp244 + Ftmp362) * M[61] + Ftmp317 * (Ftmp249 + Ftmp363) * M[63] +
          Ftmp319 * Ftmp320 + Ftmp319 * Ftmp321 * M[30] + Ftmp319 * Ftmp339 +
          Ftmp319 * Ftmp340 * M[77] + Ftmp319 * (Ftmp365 + Ftmp366) * M[68] +
          Ftmp319 * (Ftmp151 + Ftmp252 + Ftmp265) * M[79] +
          Ftmp319 * (Ftmp242 + Ftmp361 + Ftmp364) * M[66] + Ftmp321 * Ftmp323 * M[25] +
          Ftmp322 * Ftmp94 + Ftmp323 * Ftmp340 * M[70] +
          Ftmp323 * (Ftmp251 + Ftmp265 + Ftmp364) * M[72] + Ftmp324 * Ftmp91 +
          Ftmp325 * Ftmp87 * M[32] + Ftmp327 * M[35] + Ftmp328 * y * M[29] +
          Ftmp328 * M[49] + Ftmp329 * y - Ftmp332 * z - Ftmp333 * z - Ftmp334 * Ftmp5 -
          Ftmp34 * x + Ftmp341 * y + Ftmp345 * y * M[22] + Ftmp345 * M[40] +
          Ftmp346 * y * M[24] + Ftmp346 * M[42] + Ftmp347 * y * M[31] + Ftmp347 * M[51] -
          Ftmp351 * z - Ftmp354 * z - Ftmp356 * z - Ftmp357 * M[41] - Ftmp358 * M[43] -
          Ftmp359 * M[52] - Ftmp36 * Ftmp62 * M[13] + Ftmp367 * M[73] -
          Ftmp37 * Ftmp5 * M[16] - Ftmp37 * Ftmp61 * M[12] - Ftmp4 * M[7] -
          Ftmp49 * M[32] - Ftmp5 * Ftmp65 + Ftmp56 * Ftmp62 - Ftmp57 * M[9] -
          Ftmp61 * Ftmp80 + Ftmp62 * Ftmp68 * M[26] + Ftmp62 * Ftmp8 +
          y * (Ftmp289 + Ftmp372 + Ftmp373) * M[60] +
          y * (Ftmp130 + Ftmp290 + Ftmp298 + Ftmp374) * M[78] +
          y * (Ftmp145 + Ftmp280 + Ftmp293 + Ftmp369) * M[69] +
          y * (Ftmp274 + Ftmp285 + Ftmp370 + Ftmp371) * M[58] +
          y * (Ftmp281 + Ftmp291 + Ftmp294 + Ftmp370) * M[80] +
          y * (-Ftmp112 * Ftmp342 + Ftmp179 + Ftmp182 + 99225.0 * Ftmp337) * M[76] +
          y * (-Ftmp237 * Ftmp26 + Ftmp300 + Ftmp375 + Ftmp67) * M[67] +
          y * (Ftmp130 + Ftmp272 + Ftmp283 + Ftmp296 + Ftmp368) * M[65];
#pragma omp atomic
  F[2] += Ftmp0 * M[2] - Ftmp103 * Ftmp380 + Ftmp115 * M[50] + Ftmp119 * z - Ftmp134 * x -
          Ftmp138 * Ftmp306 + Ftmp14 * Ftmp318 + Ftmp14 * Ftmp338 - Ftmp141 * Ftmp171 -
          Ftmp141 * Ftmp379 * M[53] - Ftmp141 * Ftmp386 - Ftmp142 * x * M[48] -
          Ftmp144 * Ftmp5 - Ftmp148 * M[34] - Ftmp149 * Ftmp26 + Ftmp15 * M[4] +
          Ftmp162 * Ftmp381 - Ftmp164 * Ftmp311 + Ftmp171 * Ftmp387 + Ftmp173 * Ftmp322 +
          Ftmp177 * Ftmp383 + Ftmp185 * z - Ftmp2 * Ftmp26 * M[2] - Ftmp2 * Ftmp302 -
          Ftmp2 * Ftmp376 - Ftmp226 * x - Ftmp228 * x - Ftmp230 * Ftmp311 * M[46] -
          Ftmp231 * Ftmp306 - Ftmp234 * M[37] - Ftmp235 * M[39] + Ftmp245 * M[61] +
          Ftmp25 * M[16] + Ftmp250 * M[63] - Ftmp26 * Ftmp313 - Ftmp26 * Ftmp326 * Ftmp9 -
          Ftmp26 * Ftmp334 - Ftmp26 * Ftmp64 - Ftmp26 * (Ftmp352 + Ftmp360) * M[43] -
          Ftmp26 * (Ftmp47 + Ftmp77) * M[18] -
          Ftmp26 * (Ftmp139 + Ftmp163 - 13230.0 * Ftmp71) * M[54] -
          Ftmp26 * (Ftmp18 + Ftmp212 + Ftmp349) * M[41] -
          Ftmp26 * (Ftmp221 + Ftmp229 + Ftmp254 + Ftmp81) * M[52] + Ftmp27 * Ftmp302 +
          Ftmp27 * Ftmp376 + Ftmp29 * z + Ftmp303 * M[11] + Ftmp305 * z -
          Ftmp306 * Ftmp46 - Ftmp309 * y - Ftmp310 * y - Ftmp311 * Ftmp79 +
          Ftmp314 * Ftmp380 - Ftmp315 * Ftmp5 + Ftmp320 * Ftmp381 + Ftmp327 * M[36] +
          Ftmp329 * z - Ftmp332 * y - Ftmp333 * y + Ftmp335 * Ftmp380 - Ftmp336 * Ftmp5 +
          Ftmp339 * Ftmp381 + Ftmp341 * z - Ftmp351 * y - Ftmp354 * y - Ftmp356 * y -
          Ftmp357 * M[40] - Ftmp358 * M[42] - Ftmp359 * M[51] + Ftmp367 * M[72] +
          Ftmp377 * z * M[8] + Ftmp377 * M[18] + Ftmp378 * Ftmp384 - Ftmp378 * Ftmp48 -
          Ftmp379 * Ftmp48 * M[17] -
          Ftmp380 * Ftmp39 * (Ftmp105 - 1575.0 * Ftmp52) * M[47] +
          Ftmp380 * (Ftmp244 + Ftmp363) * M[62] + Ftmp380 * (Ftmp249 + Ftmp362) * M[64] +
          Ftmp381 * Ftmp76 + Ftmp381 * (Ftmp151 + Ftmp248 + Ftmp391) * M[79] +
          Ftmp381 * (Ftmp242 + Ftmp267 + Ftmp365) * M[66] +
          Ftmp381 * (Ftmp361 + Ftmp366 + Ftmp95) * M[68] + Ftmp382 * Ftmp73 * M[27] +
          Ftmp382 * (-34650.0 * Ftmp169 + Ftmp174 + Ftmp95) * M[74] + Ftmp383 * Ftmp93 +
          Ftmp383 * (Ftmp269 + Ftmp391 + Ftmp95) * M[73] + Ftmp384 * Ftmp88 +
          Ftmp385 * z * M[33] + Ftmp385 * M[54] + Ftmp386 * Ftmp387 +
          Ftmp388 * z * M[22] + Ftmp388 * M[41] + Ftmp389 * z * M[24] + Ftmp389 * M[43] +
          Ftmp390 * z * M[31] + Ftmp390 * M[52] - Ftmp4 * x * M[0] - Ftmp4 * y * M[1] -
          Ftmp42 * x - Ftmp47 * Ftmp9 * M[13] - Ftmp48 * Ftmp88 - Ftmp49 * x * M[14] -
          Ftmp5 * Ftmp58 + Ftmp56 * Ftmp9 + Ftmp62 * Ftmp89 - Ftmp63 * M[9] + Ftmp69 * x +
          Ftmp75 * x + Ftmp8 * Ftmp9 +
          z * (Ftmp136 + Ftmp275 + Ftmp296 + Ftmp369) * M[65] +
          z * (Ftmp202 + Ftmp274 + Ftmp286 + Ftmp372) * M[58] +
          z * (Ftmp270 + Ftmp288 + Ftmp371 + Ftmp373) * M[60] +
          z * (Ftmp282 + Ftmp295 + Ftmp374 + Ftmp392) * M[80] +
          z * (-Ftmp116 * Ftmp342 + Ftmp179 + Ftmp186 + 99225.0 * Ftmp71) * M[82] +
          z * (Ftmp214 + Ftmp270 + Ftmp277 + Ftmp291 + Ftmp297) * M[78] +
          z * (Ftmp278 + Ftmp287 + Ftmp293 + Ftmp368 + Ftmp392) * M[69] +
          z * (-Ftmp23 * Ftmp237 + Ftmp299 + Ftmp301 + Ftmp375 + Ftmp87) * M[67];
}

void field_m1_P2M_7(double x, double y, double z, double q, double* M) {
  double Mtmp0  = q * x;
  double Mtmp1  = q * y;
  double Mtmp2  = q * z;
  double Mtmp3  = (x * x);
  double Mtmp4  = (1.0 / 2.0) * q;
  double Mtmp5  = Mtmp0 * y;
  double Mtmp6  = Mtmp0 * z;
  double Mtmp7  = (y * y);
  double Mtmp8  = Mtmp1 * z;
  double Mtmp9  = (z * z);
  double Mtmp10 = (x * x * x);
  double Mtmp11 = (1.0 / 6.0) * q;
  double Mtmp12 = (1.0 / 2.0) * Mtmp3;
  double Mtmp13 = (1.0 / 2.0) * Mtmp0;
  double Mtmp14 = (y * y * y);
  double Mtmp15 = (1.0 / 2.0) * Mtmp7;
  double Mtmp16 = (1.0 / 2.0) * Mtmp9;
  double Mtmp17 = (z * z * z);
  double Mtmp18 = (x * x * x * x);
  double Mtmp19 = (1.0 / 24.0) * q;
  double Mtmp20 = (1.0 / 6.0) * Mtmp10;
  double Mtmp21 = Mtmp7 * q;
  double Mtmp22 = (1.0 / 4.0) * Mtmp3;
  double Mtmp23 = Mtmp9 * q;
  double Mtmp24 = (1.0 / 6.0) * Mtmp0;
  double Mtmp25 = (y * y * y * y);
  double Mtmp26 = (1.0 / 6.0) * Mtmp14;
  double Mtmp27 = (1.0 / 4.0) * Mtmp9;
  double Mtmp28 = (1.0 / 6.0) * Mtmp17;
  double Mtmp29 = (z * z * z * z);
  double Mtmp30 = (x * x * x * x * x);
  double Mtmp31 = (1.0 / 120.0) * q;
  double Mtmp32 = (1.0 / 24.0) * Mtmp18;
  double Mtmp33 = (1.0 / 12.0) * Mtmp10;
  double Mtmp34 = (1.0 / 12.0) * Mtmp14;
  double Mtmp35 = Mtmp3 * q;
  double Mtmp36 = Mtmp2 * Mtmp7;
  double Mtmp37 = Mtmp1 * Mtmp9;
  double Mtmp38 = (1.0 / 12.0) * Mtmp17;
  double Mtmp39 = (1.0 / 24.0) * Mtmp0;
  double Mtmp40 = Mtmp0 * Mtmp7;
  double Mtmp41 = (y * y * y * y * y);
  double Mtmp42 = (1.0 / 24.0) * Mtmp25;
  double Mtmp43 = (1.0 / 24.0) * Mtmp29;
  double Mtmp44 = (z * z * z * z * z);
  double Mtmp45 = (x * x * x * x * x * x);
  double Mtmp46 = (1.0 / 720.0) * q;
  double Mtmp47 = (1.0 / 120.0) * Mtmp30;
  double Mtmp48 = (1.0 / 48.0) * Mtmp18;
  double Mtmp49 = Mtmp14 * q;
  double Mtmp50 = (1.0 / 36.0) * Mtmp10;
  double Mtmp51 = Mtmp17 * q;
  double Mtmp52 = (1.0 / 48.0) * Mtmp35;
  double Mtmp53 = Mtmp2 * Mtmp3;
  double Mtmp54 = Mtmp3 * Mtmp9;
  double Mtmp55 = Mtmp1 * Mtmp3;
  double Mtmp56 = (1.0 / 120.0) * Mtmp0;
  double Mtmp57 = Mtmp0 * Mtmp9;
  double Mtmp58 = (y * y * y * y * y * y);
  double Mtmp59 = (1.0 / 120.0) * Mtmp41;
  double Mtmp60 = (1.0 / 48.0) * Mtmp25;
  double Mtmp61 = (1.0 / 36.0) * Mtmp17;
  double Mtmp62 = (1.0 / 48.0) * Mtmp29;
  double Mtmp63 = (1.0 / 120.0) * Mtmp44;
  double Mtmp64 = (z * z * z * z * z * z);
  double Mtmp65 = (1.0 / 5040.0) * q;
  double Mtmp66 = (1.0 / 720.0) * Mtmp45;
  double Mtmp67 = (1.0 / 240.0) * Mtmp30;
  double Mtmp68 = (1.0 / 144.0) * Mtmp18;
  double Mtmp69 = (1.0 / 144.0) * Mtmp25;
  double Mtmp70 = Mtmp10 * q;
  double Mtmp71 = Mtmp19 * Mtmp7;
  double Mtmp72 = (1.0 / 144.0) * Mtmp29;
  double Mtmp73 = (1.0 / 240.0) * Mtmp35;
  double Mtmp74 = (1.0 / 720.0) * Mtmp0;
  M[0] += -Mtmp0;
  M[1] += -Mtmp1;
  M[2] += -Mtmp2;
  M[3] += Mtmp3 * Mtmp4;
  M[4] += Mtmp5;
  M[5] += Mtmp6;
  M[6] += Mtmp4 * Mtmp7;
  M[7] += Mtmp8;
  M[8] += Mtmp4 * Mtmp9;
  M[9] += -Mtmp10 * Mtmp11;
  M[10] += -Mtmp1 * Mtmp12;
  M[11] += -Mtmp12 * Mtmp2;
  M[12] += -Mtmp13 * Mtmp7;
  M[13] += -Mtmp5 * z;
  M[14] += -Mtmp13 * Mtmp9;
  M[15] += -Mtmp11 * Mtmp14;
  M[16] += -Mtmp15 * Mtmp2;
  M[17] += -Mtmp1 * Mtmp16;
  M[18] += -Mtmp11 * Mtmp17;
  M[19] += Mtmp18 * Mtmp19;
  M[20] += Mtmp1 * Mtmp20;
  M[21] += Mtmp2 * Mtmp20;
  M[22] += Mtmp21 * Mtmp22;
  M[23] += Mtmp12 * Mtmp8;
  M[24] += Mtmp22 * Mtmp23;
  M[25] += Mtmp14 * Mtmp24;
  M[26] += Mtmp15 * Mtmp6;
  M[27] += Mtmp16 * Mtmp5;
  M[28] += Mtmp17 * Mtmp24;
  M[29] += Mtmp19 * Mtmp25;
  M[30] += Mtmp2 * Mtmp26;
  M[31] += Mtmp21 * Mtmp27;
  M[32] += Mtmp1 * Mtmp28;
  M[33] += Mtmp19 * Mtmp29;
  M[34] += -Mtmp30 * Mtmp31;
  M[35] += -Mtmp1 * Mtmp32;
  M[36] += -Mtmp2 * Mtmp32;
  M[37] += -Mtmp21 * Mtmp33;
  M[38] += -Mtmp20 * Mtmp8;
  M[39] += -Mtmp23 * Mtmp33;
  M[40] += -Mtmp34 * Mtmp35;
  M[41] += -Mtmp22 * Mtmp36;
  M[42] += -Mtmp22 * Mtmp37;
  M[43] += -Mtmp35 * Mtmp38;
  M[44] += -Mtmp25 * Mtmp39;
  M[45] += -Mtmp26 * Mtmp6;
  M[46] += -Mtmp27 * Mtmp40;
  M[47] += -Mtmp28 * Mtmp5;
  M[48] += -Mtmp29 * Mtmp39;
  M[49] += -Mtmp31 * Mtmp41;
  M[50] += -Mtmp2 * Mtmp42;
  M[51] += -Mtmp23 * Mtmp34;
  M[52] += -Mtmp21 * Mtmp38;
  M[53] += -Mtmp1 * Mtmp43;
  M[54] += -Mtmp31 * Mtmp44;
  M[55] += Mtmp45 * Mtmp46;
  M[56] += Mtmp1 * Mtmp47;
  M[57] += Mtmp2 * Mtmp47;
  M[58] += Mtmp21 * Mtmp48;
  M[59] += Mtmp32 * Mtmp8;
  M[60] += Mtmp23 * Mtmp48;
  M[61] += Mtmp49 * Mtmp50;
  M[62] += Mtmp33 * Mtmp36;
  M[63] += Mtmp33 * Mtmp37;
  M[64] += Mtmp50 * Mtmp51;
  M[65] += Mtmp25 * Mtmp52;
  M[66] += Mtmp34 * Mtmp53;
  M[67] += (1.0 / 8.0) * Mtmp21 * Mtmp54;
  M[68] += Mtmp38 * Mtmp55;
  M[69] += Mtmp29 * Mtmp52;
  M[70] += Mtmp41 * Mtmp56;
  M[71] += Mtmp42 * Mtmp6;
  M[72] += Mtmp34 * Mtmp57;
  M[73] += Mtmp38 * Mtmp40;
  M[74] += Mtmp43 * Mtmp5;
  M[75] += Mtmp44 * Mtmp56;
  M[76] += Mtmp46 * Mtmp58;
  M[77] += Mtmp2 * Mtmp59;
  M[78] += Mtmp23 * Mtmp60;
  M[79] += Mtmp49 * Mtmp61;
  M[80] += Mtmp21 * Mtmp62;
  M[81] += Mtmp1 * Mtmp63;
  M[82] += Mtmp46 * Mtmp64;
  M[83] += -Mtmp65 * (x * x * x * x * x * x * x);
  M[84] += -Mtmp1 * Mtmp66;
  M[85] += -Mtmp2 * Mtmp66;
  M[86] += -Mtmp21 * Mtmp67;
  M[87] += -Mtmp47 * Mtmp8;
  M[88] += -Mtmp23 * Mtmp67;
  M[89] += -Mtmp49 * Mtmp68;
  M[90] += -Mtmp36 * Mtmp48;
  M[91] += -Mtmp37 * Mtmp48;
  M[92] += -Mtmp51 * Mtmp68;
  M[93] += -Mtmp69 * Mtmp70;
  M[94] += -Mtmp14 * Mtmp2 * Mtmp50;
  M[95] += -Mtmp10 * Mtmp71 * Mtmp9;
  M[96] += -Mtmp1 * Mtmp17 * Mtmp50;
  M[97] += -Mtmp70 * Mtmp72;
  M[98] += -Mtmp41 * Mtmp73;
  M[99] += -Mtmp53 * Mtmp60;
  M[100] += -Mtmp14 * Mtmp19 * Mtmp54;
  M[101] += -Mtmp17 * Mtmp3 * Mtmp71;
  M[102] += -Mtmp55 * Mtmp62;
  M[103] += -Mtmp44 * Mtmp73;
  M[104] += -Mtmp58 * Mtmp74;
  M[105] += -Mtmp59 * Mtmp6;
  M[106] += -Mtmp57 * Mtmp60;
  M[107] += -Mtmp0 * Mtmp14 * Mtmp61;
  M[108] += -Mtmp40 * Mtmp62;
  M[109] += -Mtmp5 * Mtmp63;
  M[110] += -Mtmp64 * Mtmp74;
  M[111] += -Mtmp65 * (y * y * y * y * y * y * y);
  M[112] += -1.0 / 720.0 * Mtmp2 * Mtmp58;
  M[113] += -1.0 / 240.0 * Mtmp23 * Mtmp41;
  M[114] += -Mtmp51 * Mtmp69;
  M[115] += -Mtmp49 * Mtmp72;
  M[116] += -1.0 / 240.0 * Mtmp21 * Mtmp44;
  M[117] += -1.0 / 720.0 * Mtmp1 * Mtmp64;
  M[118] += -Mtmp65 * (z * z * z * z * z * z * z);
}
void field_m1_M2M_7(double x, double y, double z, double* M, double* Ms) {
  double Mstmp0   = x * M[0];
  double Mstmp1   = x * M[1];
  double Mstmp2   = y * M[0];
  double Mstmp3   = x * M[2];
  double Mstmp4   = z * M[0];
  double Mstmp5   = y * M[1];
  double Mstmp6   = y * M[2];
  double Mstmp7   = z * M[1];
  double Mstmp8   = z * M[2];
  double Mstmp9   = x * M[3];
  double Mstmp10  = (x * x);
  double Mstmp11  = (1.0 / 2.0) * Mstmp10;
  double Mstmp12  = x * M[4];
  double Mstmp13  = y * M[3];
  double Mstmp14  = Mstmp0 * y;
  double Mstmp15  = x * M[5];
  double Mstmp16  = z * M[3];
  double Mstmp17  = Mstmp0 * z;
  double Mstmp18  = x * M[6];
  double Mstmp19  = y * M[4];
  double Mstmp20  = Mstmp1 * y;
  double Mstmp21  = (y * y);
  double Mstmp22  = (1.0 / 2.0) * M[0];
  double Mstmp23  = x * M[7];
  double Mstmp24  = y * M[5];
  double Mstmp25  = z * M[4];
  double Mstmp26  = Mstmp3 * y;
  double Mstmp27  = Mstmp1 * z;
  double Mstmp28  = Mstmp2 * z;
  double Mstmp29  = x * M[8];
  double Mstmp30  = z * M[5];
  double Mstmp31  = Mstmp3 * z;
  double Mstmp32  = (z * z);
  double Mstmp33  = y * M[6];
  double Mstmp34  = (1.0 / 2.0) * Mstmp21;
  double Mstmp35  = y * M[7];
  double Mstmp36  = z * M[6];
  double Mstmp37  = Mstmp5 * z;
  double Mstmp38  = y * M[8];
  double Mstmp39  = z * M[7];
  double Mstmp40  = Mstmp6 * z;
  double Mstmp41  = (1.0 / 2.0) * Mstmp32;
  double Mstmp42  = z * M[8];
  double Mstmp43  = x * M[9];
  double Mstmp44  = (x * x * x);
  double Mstmp45  = (1.0 / 6.0) * Mstmp44;
  double Mstmp46  = x * M[10];
  double Mstmp47  = y * M[9];
  double Mstmp48  = Mstmp9 * y;
  double Mstmp49  = x * M[11];
  double Mstmp50  = z * M[9];
  double Mstmp51  = Mstmp9 * z;
  double Mstmp52  = x * M[12];
  double Mstmp53  = y * M[10];
  double Mstmp54  = Mstmp12 * y;
  double Mstmp55  = x * M[13];
  double Mstmp56  = y * M[11];
  double Mstmp57  = z * M[10];
  double Mstmp58  = Mstmp15 * y;
  double Mstmp59  = Mstmp12 * z;
  double Mstmp60  = Mstmp13 * z;
  double Mstmp61  = x * M[14];
  double Mstmp62  = z * M[11];
  double Mstmp63  = Mstmp15 * z;
  double Mstmp64  = x * M[15];
  double Mstmp65  = y * M[12];
  double Mstmp66  = Mstmp18 * y;
  double Mstmp67  = (y * y * y);
  double Mstmp68  = (1.0 / 6.0) * M[0];
  double Mstmp69  = x * M[16];
  double Mstmp70  = y * M[13];
  double Mstmp71  = z * M[12];
  double Mstmp72  = Mstmp23 * y;
  double Mstmp73  = Mstmp18 * z;
  double Mstmp74  = Mstmp19 * z;
  double Mstmp75  = x * M[17];
  double Mstmp76  = y * M[14];
  double Mstmp77  = z * M[13];
  double Mstmp78  = Mstmp29 * y;
  double Mstmp79  = Mstmp23 * z;
  double Mstmp80  = Mstmp24 * z;
  double Mstmp81  = x * M[18];
  double Mstmp82  = z * M[14];
  double Mstmp83  = Mstmp29 * z;
  double Mstmp84  = (z * z * z);
  double Mstmp85  = y * M[15];
  double Mstmp86  = (1.0 / 6.0) * Mstmp67;
  double Mstmp87  = y * M[16];
  double Mstmp88  = z * M[15];
  double Mstmp89  = Mstmp33 * z;
  double Mstmp90  = y * M[17];
  double Mstmp91  = z * M[16];
  double Mstmp92  = Mstmp35 * z;
  double Mstmp93  = y * M[18];
  double Mstmp94  = z * M[17];
  double Mstmp95  = Mstmp38 * z;
  double Mstmp96  = (1.0 / 6.0) * Mstmp84;
  double Mstmp97  = z * M[18];
  double Mstmp98  = x * M[19];
  double Mstmp99  = (x * x * x * x);
  double Mstmp100 = (1.0 / 24.0) * Mstmp99;
  double Mstmp101 = x * M[20];
  double Mstmp102 = y * M[19];
  double Mstmp103 = Mstmp43 * y;
  double Mstmp104 = x * M[21];
  double Mstmp105 = z * M[19];
  double Mstmp106 = Mstmp43 * z;
  double Mstmp107 = x * M[22];
  double Mstmp108 = y * M[20];
  double Mstmp109 = Mstmp46 * y;
  double Mstmp110 = (1.0 / 4.0) * Mstmp10;
  double Mstmp111 = Mstmp21 * M[0];
  double Mstmp112 = x * M[23];
  double Mstmp113 = y * M[21];
  double Mstmp114 = z * M[20];
  double Mstmp115 = Mstmp49 * y;
  double Mstmp116 = Mstmp46 * z;
  double Mstmp117 = Mstmp47 * z;
  double Mstmp118 = x * M[24];
  double Mstmp119 = z * M[21];
  double Mstmp120 = Mstmp49 * z;
  double Mstmp121 = Mstmp110 * Mstmp32;
  double Mstmp122 = x * M[25];
  double Mstmp123 = y * M[22];
  double Mstmp124 = Mstmp52 * y;
  double Mstmp125 = Mstmp110 * Mstmp21;
  double Mstmp126 = x * M[26];
  double Mstmp127 = y * M[23];
  double Mstmp128 = z * M[22];
  double Mstmp129 = Mstmp55 * y;
  double Mstmp130 = Mstmp52 * z;
  double Mstmp131 = Mstmp53 * z;
  double Mstmp132 = x * M[27];
  double Mstmp133 = y * M[24];
  double Mstmp134 = z * M[23];
  double Mstmp135 = Mstmp61 * y;
  double Mstmp136 = Mstmp55 * z;
  double Mstmp137 = Mstmp56 * z;
  double Mstmp138 = x * M[28];
  double Mstmp139 = z * M[24];
  double Mstmp140 = Mstmp61 * z;
  double Mstmp141 = x * M[29];
  double Mstmp142 = y * M[25];
  double Mstmp143 = Mstmp64 * y;
  double Mstmp144 = (y * y * y * y);
  double Mstmp145 = (1.0 / 24.0) * M[0];
  double Mstmp146 = x * M[30];
  double Mstmp147 = y * M[26];
  double Mstmp148 = z * M[25];
  double Mstmp149 = Mstmp69 * y;
  double Mstmp150 = Mstmp64 * z;
  double Mstmp151 = Mstmp65 * z;
  double Mstmp152 = x * M[31];
  double Mstmp153 = y * M[27];
  double Mstmp154 = z * M[26];
  double Mstmp155 = Mstmp75 * y;
  double Mstmp156 = Mstmp69 * z;
  double Mstmp157 = Mstmp70 * z;
  double Mstmp158 = (1.0 / 4.0) * Mstmp32;
  double Mstmp159 = x * M[32];
  double Mstmp160 = y * M[28];
  double Mstmp161 = z * M[27];
  double Mstmp162 = Mstmp81 * y;
  double Mstmp163 = Mstmp75 * z;
  double Mstmp164 = Mstmp76 * z;
  double Mstmp165 = x * M[33];
  double Mstmp166 = z * M[28];
  double Mstmp167 = Mstmp81 * z;
  double Mstmp168 = (z * z * z * z);
  double Mstmp169 = y * M[29];
  double Mstmp170 = (1.0 / 24.0) * Mstmp144;
  double Mstmp171 = y * M[30];
  double Mstmp172 = z * M[29];
  double Mstmp173 = Mstmp85 * z;
  double Mstmp174 = y * M[31];
  double Mstmp175 = z * M[30];
  double Mstmp176 = Mstmp87 * z;
  double Mstmp177 = Mstmp158 * Mstmp21;
  double Mstmp178 = y * M[32];
  double Mstmp179 = z * M[31];
  double Mstmp180 = Mstmp90 * z;
  double Mstmp181 = y * M[33];
  double Mstmp182 = z * M[32];
  double Mstmp183 = Mstmp93 * z;
  double Mstmp184 = (1.0 / 24.0) * Mstmp168;
  double Mstmp185 = z * M[33];
  double Mstmp186 = x * M[34];
  double Mstmp187 = (1.0 / 120.0) * (x * x * x * x * x);
  double Mstmp188 = x * M[35];
  double Mstmp189 = y * M[34];
  double Mstmp190 = Mstmp98 * y;
  double Mstmp191 = x * M[36];
  double Mstmp192 = x * M[37];
  double Mstmp193 = y * M[35];
  double Mstmp194 = Mstmp101 * y;
  double Mstmp195 = (1.0 / 12.0) * Mstmp44;
  double Mstmp196 = x * M[38];
  double Mstmp197 = y * M[36];
  double Mstmp198 = Mstmp104 * y;
  double Mstmp199 = x * M[39];
  double Mstmp200 = Mstmp195 * Mstmp32;
  double Mstmp201 = x * M[40];
  double Mstmp202 = y * M[37];
  double Mstmp203 = Mstmp107 * y;
  double Mstmp204 = (1.0 / 12.0) * Mstmp10;
  double Mstmp205 = Mstmp67 * M[0];
  double Mstmp206 = Mstmp195 * Mstmp21;
  double Mstmp207 = x * M[41];
  double Mstmp208 = y * M[38];
  double Mstmp209 = Mstmp112 * y;
  double Mstmp210 = x * M[42];
  double Mstmp211 = y * M[39];
  double Mstmp212 = Mstmp118 * y;
  double Mstmp213 = x * M[43];
  double Mstmp214 = Mstmp204 * Mstmp84;
  double Mstmp215 = x * M[44];
  double Mstmp216 = y * M[40];
  double Mstmp217 = Mstmp122 * y;
  double Mstmp218 = Mstmp204 * Mstmp67;
  double Mstmp219 = x * M[45];
  double Mstmp220 = y * M[41];
  double Mstmp221 = Mstmp126 * y;
  double Mstmp222 = x * M[46];
  double Mstmp223 = y * M[42];
  double Mstmp224 = Mstmp132 * y;
  double Mstmp225 = x * M[47];
  double Mstmp226 = y * M[43];
  double Mstmp227 = Mstmp138 * y;
  double Mstmp228 = x * M[48];
  double Mstmp229 = x * M[49];
  double Mstmp230 = y * M[44];
  double Mstmp231 = Mstmp141 * y;
  double Mstmp232 = (y * y * y * y * y);
  double Mstmp233 = (1.0 / 120.0) * M[0];
  double Mstmp234 = x * M[50];
  double Mstmp235 = y * M[45];
  double Mstmp236 = Mstmp146 * y;
  double Mstmp237 = x * M[51];
  double Mstmp238 = y * M[46];
  double Mstmp239 = Mstmp152 * y;
  double Mstmp240 = (1.0 / 12.0) * Mstmp32;
  double Mstmp241 = x * M[52];
  double Mstmp242 = y * M[47];
  double Mstmp243 = Mstmp159 * y;
  double Mstmp244 = (1.0 / 12.0) * Mstmp84;
  double Mstmp245 = x * M[53];
  double Mstmp246 = y * M[48];
  double Mstmp247 = Mstmp165 * y;
  double Mstmp248 = x * M[54];
  double Mstmp249 = (z * z * z * z * z);
  double Mstmp250 = y * M[49];
  double Mstmp251 = (1.0 / 120.0) * Mstmp232;
  double Mstmp252 = y * M[50];
  double Mstmp253 = y * M[51];
  double Mstmp254 = Mstmp240 * Mstmp67;
  double Mstmp255 = y * M[52];
  double Mstmp256 = Mstmp21 * Mstmp244;
  double Mstmp257 = y * M[53];
  double Mstmp258 = y * M[54];
  double Mstmp259 = (1.0 / 120.0) * Mstmp249;
  double Mstmp260 = (1.0 / 720.0) * (x * x * x * x * x * x);
  double Mstmp261 = (1.0 / 48.0) * Mstmp99;
  double Mstmp262 = Mstmp261 * Mstmp32;
  double Mstmp263 = (1.0 / 36.0) * Mstmp44;
  double Mstmp264 = Mstmp21 * Mstmp261;
  double Mstmp265 = Mstmp263 * Mstmp84;
  double Mstmp266 = (1.0 / 48.0) * Mstmp10;
  double Mstmp267 = Mstmp266 * M[0];
  double Mstmp268 = Mstmp263 * Mstmp67;
  double Mstmp269 = (1.0 / 8.0) * Mstmp10 * Mstmp32;
  double Mstmp270 = Mstmp144 * Mstmp266;
  double Mstmp271 = Mstmp21 * Mstmp269;
  double Mstmp272 = Mstmp168 * Mstmp266;
  double Mstmp273 = (y * y * y * y * y * y);
  double Mstmp274 = (1.0 / 720.0) * M[0];
  double Mstmp275 = (1.0 / 48.0) * Mstmp144 * Mstmp32;
  double Mstmp276 = (1.0 / 36.0) * Mstmp84;
  double Mstmp277 = (1.0 / 48.0) * Mstmp168;
  double Mstmp278 = (z * z * z * z * z * z);
  double Mstmp279 = (1.0 / 720.0) * Mstmp273;
  double Mstmp280 = Mstmp276 * Mstmp67;
  double Mstmp281 = Mstmp21 * Mstmp277;
  double Mstmp282 = (1.0 / 720.0) * Mstmp278;
#pragma omp atomic
  Ms[0] += M[0];
#pragma omp atomic
  Ms[1] += M[1];
#pragma omp atomic
  Ms[2] += M[2];
#pragma omp atomic
  Ms[3] += Mstmp0 + M[3];
#pragma omp atomic
  Ms[4] += Mstmp1 + Mstmp2 + M[4];
#pragma omp atomic
  Ms[5] += Mstmp3 + Mstmp4 + M[5];
#pragma omp atomic
  Ms[6] += Mstmp5 + M[6];
#pragma omp atomic
  Ms[7] += Mstmp6 + Mstmp7 + M[7];
#pragma omp atomic
  Ms[8] += Mstmp8 + M[8];
#pragma omp atomic
  Ms[9] += Mstmp11 * M[0] + Mstmp9 + M[9];
#pragma omp atomic
  Ms[10] += Mstmp11 * M[1] + Mstmp12 + Mstmp13 + Mstmp14 + M[10];
#pragma omp atomic
  Ms[11] += Mstmp11 * M[2] + Mstmp15 + Mstmp16 + Mstmp17 + M[11];
#pragma omp atomic
  Ms[12] += Mstmp18 + Mstmp19 + Mstmp20 + Mstmp21 * Mstmp22 + M[12];
#pragma omp atomic
  Ms[13] += Mstmp23 + Mstmp24 + Mstmp25 + Mstmp26 + Mstmp27 + Mstmp28 + M[13];
#pragma omp atomic
  Ms[14] += Mstmp22 * Mstmp32 + Mstmp29 + Mstmp30 + Mstmp31 + M[14];
#pragma omp atomic
  Ms[15] += Mstmp33 + Mstmp34 * M[1] + M[15];
#pragma omp atomic
  Ms[16] += Mstmp34 * M[2] + Mstmp35 + Mstmp36 + Mstmp37 + M[16];
#pragma omp atomic
  Ms[17] += Mstmp38 + Mstmp39 + Mstmp40 + Mstmp41 * M[1] + M[17];
#pragma omp atomic
  Ms[18] += Mstmp41 * M[2] + Mstmp42 + M[18];
#pragma omp atomic
  Ms[19] += Mstmp11 * M[3] + Mstmp43 + Mstmp45 * M[0] + M[19];
#pragma omp atomic
  Ms[20] += Mstmp11 * Mstmp2 + Mstmp11 * M[4] + Mstmp45 * M[1] + Mstmp46 + Mstmp47 +
            Mstmp48 + M[20];
#pragma omp atomic
  Ms[21] += Mstmp11 * Mstmp4 + Mstmp11 * M[5] + Mstmp45 * M[2] + Mstmp49 + Mstmp50 +
            Mstmp51 + M[21];
#pragma omp atomic
  Ms[22] += Mstmp0 * Mstmp34 + Mstmp11 * Mstmp5 + Mstmp11 * M[6] + Mstmp34 * M[3] +
            Mstmp52 + Mstmp53 + Mstmp54 + M[22];
#pragma omp atomic
  Ms[23] += Mstmp11 * Mstmp6 + Mstmp11 * Mstmp7 + Mstmp11 * M[7] + Mstmp14 * z + Mstmp55 +
            Mstmp56 + Mstmp57 + Mstmp58 + Mstmp59 + Mstmp60 + M[23];
#pragma omp atomic
  Ms[24] += Mstmp0 * Mstmp41 + Mstmp11 * Mstmp8 + Mstmp11 * M[8] + Mstmp41 * M[3] +
            Mstmp61 + Mstmp62 + Mstmp63 + M[24];
#pragma omp atomic
  Ms[25] += Mstmp1 * Mstmp34 + Mstmp34 * M[4] + Mstmp64 + Mstmp65 + Mstmp66 +
            Mstmp67 * Mstmp68 + M[25];
#pragma omp atomic
  Ms[26] += Mstmp20 * z + Mstmp3 * Mstmp34 + Mstmp34 * Mstmp4 + Mstmp34 * M[5] + Mstmp69 +
            Mstmp70 + Mstmp71 + Mstmp72 + Mstmp73 + Mstmp74 + M[26];
#pragma omp atomic
  Ms[27] += Mstmp1 * Mstmp41 + Mstmp2 * Mstmp41 + Mstmp26 * z + Mstmp41 * M[4] + Mstmp75 +
            Mstmp76 + Mstmp77 + Mstmp78 + Mstmp79 + Mstmp80 + M[27];
#pragma omp atomic
  Ms[28] += Mstmp3 * Mstmp41 + Mstmp41 * M[5] + Mstmp68 * Mstmp84 + Mstmp81 + Mstmp82 +
            Mstmp83 + M[28];
#pragma omp atomic
  Ms[29] += Mstmp34 * M[6] + Mstmp85 + Mstmp86 * M[1] + M[29];
#pragma omp atomic
  Ms[30] += Mstmp34 * Mstmp7 + Mstmp34 * M[7] + Mstmp86 * M[2] + Mstmp87 + Mstmp88 +
            Mstmp89 + M[30];
#pragma omp atomic
  Ms[31] += Mstmp34 * Mstmp8 + Mstmp34 * M[8] + Mstmp41 * Mstmp5 + Mstmp41 * M[6] +
            Mstmp90 + Mstmp91 + Mstmp92 + M[31];
#pragma omp atomic
  Ms[32] += Mstmp41 * Mstmp6 + Mstmp41 * M[7] + Mstmp93 + Mstmp94 + Mstmp95 +
            Mstmp96 * M[1] + M[32];
#pragma omp atomic
  Ms[33] += Mstmp41 * M[8] + Mstmp96 * M[2] + Mstmp97 + M[33];
#pragma omp atomic
  Ms[34] += Mstmp100 * M[0] + Mstmp11 * M[9] + Mstmp45 * M[3] + Mstmp98 + M[34];
#pragma omp atomic
  Ms[35] += Mstmp100 * M[1] + Mstmp101 + Mstmp102 + Mstmp103 + Mstmp11 * Mstmp13 +
            Mstmp11 * M[10] + Mstmp2 * Mstmp45 + Mstmp45 * M[4] + M[35];
#pragma omp atomic
  Ms[36] += Mstmp100 * M[2] + Mstmp104 + Mstmp105 + Mstmp106 + Mstmp11 * Mstmp16 +
            Mstmp11 * M[11] + Mstmp4 * Mstmp45 + Mstmp45 * M[5] + M[36];
#pragma omp atomic
  Ms[37] += Mstmp107 + Mstmp108 + Mstmp109 + Mstmp11 * Mstmp19 + Mstmp11 * M[12] +
            Mstmp110 * Mstmp111 + Mstmp34 * Mstmp9 + Mstmp34 * M[9] + Mstmp45 * Mstmp5 +
            Mstmp45 * M[6] + M[37];
#pragma omp atomic
  Ms[38] += Mstmp11 * Mstmp24 + Mstmp11 * Mstmp25 + Mstmp11 * Mstmp28 + Mstmp11 * M[13] +
            Mstmp112 + Mstmp113 + Mstmp114 + Mstmp115 + Mstmp116 + Mstmp117 +
            Mstmp45 * Mstmp6 + Mstmp45 * Mstmp7 + Mstmp45 * M[7] + Mstmp48 * z + M[38];
#pragma omp atomic
  Ms[39] += Mstmp11 * Mstmp30 + Mstmp11 * M[14] + Mstmp118 + Mstmp119 + Mstmp120 +
            Mstmp121 * M[0] + Mstmp41 * Mstmp9 + Mstmp41 * M[9] + Mstmp45 * Mstmp8 +
            Mstmp45 * M[8] + M[39];
#pragma omp atomic
  Ms[40] += Mstmp0 * Mstmp86 + Mstmp11 * Mstmp33 + Mstmp11 * M[15] + Mstmp12 * Mstmp34 +
            Mstmp122 + Mstmp123 + Mstmp124 + Mstmp125 * M[1] + Mstmp34 * M[10] +
            Mstmp86 * M[3] + M[40];
#pragma omp atomic
  Ms[41] += Mstmp11 * Mstmp35 + Mstmp11 * Mstmp36 + Mstmp11 * Mstmp37 + Mstmp11 * M[16] +
            Mstmp125 * M[2] + Mstmp126 + Mstmp127 + Mstmp128 + Mstmp129 + Mstmp130 +
            Mstmp131 + Mstmp15 * Mstmp34 + Mstmp16 * Mstmp34 + Mstmp17 * Mstmp34 +
            Mstmp34 * M[11] + Mstmp54 * z + M[41];
#pragma omp atomic
  Ms[42] += Mstmp11 * Mstmp38 + Mstmp11 * Mstmp39 + Mstmp11 * Mstmp40 + Mstmp11 * M[17] +
            Mstmp12 * Mstmp41 + Mstmp121 * M[1] + Mstmp13 * Mstmp41 + Mstmp132 +
            Mstmp133 + Mstmp134 + Mstmp135 + Mstmp136 + Mstmp137 + Mstmp14 * Mstmp41 +
            Mstmp41 * M[10] + Mstmp58 * z + M[42];
#pragma omp atomic
  Ms[43] += Mstmp0 * Mstmp96 + Mstmp11 * Mstmp42 + Mstmp11 * M[18] + Mstmp121 * M[2] +
            Mstmp138 + Mstmp139 + Mstmp140 + Mstmp15 * Mstmp41 + Mstmp41 * M[11] +
            Mstmp96 * M[3] + M[43];
#pragma omp atomic
  Ms[44] += Mstmp1 * Mstmp86 + Mstmp141 + Mstmp142 + Mstmp143 + Mstmp144 * Mstmp145 +
            Mstmp18 * Mstmp34 + Mstmp34 * M[12] + Mstmp86 * M[4] + M[44];
#pragma omp atomic
  Ms[45] += Mstmp146 + Mstmp147 + Mstmp148 + Mstmp149 + Mstmp150 + Mstmp151 +
            Mstmp23 * Mstmp34 + Mstmp25 * Mstmp34 + Mstmp27 * Mstmp34 + Mstmp3 * Mstmp86 +
            Mstmp34 * M[13] + Mstmp4 * Mstmp86 + Mstmp66 * z + Mstmp86 * M[5] + M[45];
#pragma omp atomic
  Ms[46] += Mstmp111 * Mstmp158 + Mstmp152 + Mstmp153 + Mstmp154 + Mstmp155 + Mstmp156 +
            Mstmp157 + Mstmp18 * Mstmp41 + Mstmp19 * Mstmp41 + Mstmp20 * Mstmp41 +
            Mstmp29 * Mstmp34 + Mstmp30 * Mstmp34 + Mstmp31 * Mstmp34 + Mstmp34 * M[14] +
            Mstmp41 * M[12] + Mstmp72 * z + M[46];
#pragma omp atomic
  Ms[47] += Mstmp1 * Mstmp96 + Mstmp159 + Mstmp160 + Mstmp161 + Mstmp162 + Mstmp163 +
            Mstmp164 + Mstmp2 * Mstmp96 + Mstmp23 * Mstmp41 + Mstmp24 * Mstmp41 +
            Mstmp26 * Mstmp41 + Mstmp41 * M[13] + Mstmp78 * z + Mstmp96 * M[4] + M[47];
#pragma omp atomic
  Ms[48] += Mstmp145 * Mstmp168 + Mstmp165 + Mstmp166 + Mstmp167 + Mstmp29 * Mstmp41 +
            Mstmp3 * Mstmp96 + Mstmp41 * M[14] + Mstmp96 * M[5] + M[48];
#pragma omp atomic
  Ms[49] += Mstmp169 + Mstmp170 * M[1] + Mstmp34 * M[15] + Mstmp86 * M[6] + M[49];
#pragma omp atomic
  Ms[50] += Mstmp170 * M[2] + Mstmp171 + Mstmp172 + Mstmp173 + Mstmp34 * Mstmp36 +
            Mstmp34 * M[16] + Mstmp7 * Mstmp86 + Mstmp86 * M[7] + M[50];
#pragma omp atomic
  Ms[51] += Mstmp174 + Mstmp175 + Mstmp176 + Mstmp177 * M[1] + Mstmp33 * Mstmp41 +
            Mstmp34 * Mstmp39 + Mstmp34 * M[17] + Mstmp41 * M[15] + Mstmp8 * Mstmp86 +
            Mstmp86 * M[8] + M[51];
#pragma omp atomic
  Ms[52] += Mstmp177 * M[2] + Mstmp178 + Mstmp179 + Mstmp180 + Mstmp34 * Mstmp42 +
            Mstmp34 * M[18] + Mstmp35 * Mstmp41 + Mstmp41 * M[16] + Mstmp5 * Mstmp96 +
            Mstmp96 * M[6] + M[52];
#pragma omp atomic
  Ms[53] += Mstmp181 + Mstmp182 + Mstmp183 + Mstmp184 * M[1] + Mstmp38 * Mstmp41 +
            Mstmp41 * M[17] + Mstmp6 * Mstmp96 + Mstmp96 * M[7] + M[53];
#pragma omp atomic
  Ms[54] += Mstmp184 * M[2] + Mstmp185 + Mstmp41 * M[18] + Mstmp96 * M[8] + M[54];
#pragma omp atomic
  Ms[55] += Mstmp100 * M[3] + Mstmp11 * M[19] + Mstmp186 + Mstmp187 * M[0] +
            Mstmp45 * M[9] + M[55];
#pragma omp atomic
  Ms[56] += Mstmp100 * Mstmp2 + Mstmp100 * M[4] + Mstmp11 * Mstmp47 + Mstmp11 * M[20] +
            Mstmp13 * Mstmp45 + Mstmp187 * M[1] + Mstmp188 + Mstmp189 + Mstmp190 +
            Mstmp45 * M[10] + M[56];
#pragma omp atomic
  Ms[57] += Mstmp100 * Mstmp4 + Mstmp100 * M[5] + Mstmp11 * Mstmp50 + Mstmp11 * M[21] +
            Mstmp16 * Mstmp45 + Mstmp187 * M[2] + Mstmp191 + Mstmp45 * M[11] +
            Mstmp98 * z + z * M[34] + M[57];
#pragma omp atomic
  Ms[58] += Mstmp100 * Mstmp5 + Mstmp100 * M[6] + Mstmp11 * Mstmp53 + Mstmp11 * M[22] +
            Mstmp111 * Mstmp195 + Mstmp125 * M[3] + Mstmp19 * Mstmp45 + Mstmp192 +
            Mstmp193 + Mstmp194 + Mstmp34 * Mstmp43 + Mstmp34 * M[19] + Mstmp45 * M[12] +
            M[58];
#pragma omp atomic
  Ms[59] += Mstmp100 * Mstmp6 + Mstmp100 * Mstmp7 + Mstmp100 * M[7] + Mstmp101 * z +
            Mstmp102 * z + Mstmp103 * z + Mstmp11 * Mstmp56 + Mstmp11 * Mstmp57 +
            Mstmp11 * Mstmp60 + Mstmp11 * M[23] + Mstmp196 + Mstmp197 + Mstmp198 +
            Mstmp24 * Mstmp45 + Mstmp25 * Mstmp45 + Mstmp28 * Mstmp45 + Mstmp45 * M[13] +
            z * M[35] + M[59];
#pragma omp atomic
  Ms[60] += Mstmp100 * Mstmp8 + Mstmp100 * M[8] + Mstmp104 * z + Mstmp11 * Mstmp62 +
            Mstmp11 * M[24] + Mstmp121 * M[3] + Mstmp199 + Mstmp200 * M[0] +
            Mstmp30 * Mstmp45 + Mstmp41 * Mstmp43 + Mstmp41 * M[19] + Mstmp45 * M[14] +
            z * M[36] + M[60];
#pragma omp atomic
  Ms[61] += Mstmp11 * Mstmp65 + Mstmp11 * M[25] + Mstmp125 * M[4] + Mstmp201 + Mstmp202 +
            Mstmp203 + Mstmp204 * Mstmp205 + Mstmp206 * M[1] + Mstmp33 * Mstmp45 +
            Mstmp34 * Mstmp46 + Mstmp34 * M[20] + Mstmp45 * M[15] + Mstmp86 * Mstmp9 +
            Mstmp86 * M[9] + M[61];
#pragma omp atomic
  Ms[62] += Mstmp107 * z + Mstmp108 * z + Mstmp109 * z + Mstmp11 * Mstmp70 +
            Mstmp11 * Mstmp71 + Mstmp11 * Mstmp74 + Mstmp11 * M[26] + Mstmp125 * Mstmp4 +
            Mstmp125 * M[5] + Mstmp206 * M[2] + Mstmp207 + Mstmp208 + Mstmp209 +
            Mstmp34 * Mstmp49 + Mstmp34 * Mstmp50 + Mstmp34 * Mstmp51 + Mstmp34 * M[21] +
            Mstmp35 * Mstmp45 + Mstmp36 * Mstmp45 + Mstmp37 * Mstmp45 + Mstmp45 * M[16] +
            z * M[37] + M[62];
#pragma omp atomic
  Ms[63] += Mstmp11 * Mstmp76 + Mstmp11 * Mstmp77 + Mstmp11 * Mstmp80 + Mstmp11 * M[27] +
            Mstmp112 * z + Mstmp113 * z + Mstmp115 * z + Mstmp121 * Mstmp2 +
            Mstmp121 * M[4] + Mstmp200 * M[1] + Mstmp210 + Mstmp211 + Mstmp212 +
            Mstmp38 * Mstmp45 + Mstmp39 * Mstmp45 + Mstmp40 * Mstmp45 +
            Mstmp41 * Mstmp46 + Mstmp41 * Mstmp47 + Mstmp41 * Mstmp48 + Mstmp41 * M[20] +
            Mstmp45 * M[17] + z * M[38] + M[63];
#pragma omp atomic
  Ms[64] += Mstmp11 * Mstmp82 + Mstmp11 * M[28] + Mstmp118 * z + Mstmp121 * M[5] +
            Mstmp200 * M[2] + Mstmp213 + Mstmp214 * M[0] + Mstmp41 * Mstmp49 +
            Mstmp41 * M[21] + Mstmp42 * Mstmp45 + Mstmp45 * M[18] + Mstmp9 * Mstmp96 +
            Mstmp96 * M[9] + z * M[39] + M[64];
#pragma omp atomic
  Ms[65] += Mstmp0 * Mstmp170 + Mstmp11 * Mstmp85 + Mstmp11 * M[29] + Mstmp12 * Mstmp86 +
            Mstmp125 * M[6] + Mstmp170 * M[3] + Mstmp215 + Mstmp216 + Mstmp217 +
            Mstmp218 * M[1] + Mstmp34 * Mstmp52 + Mstmp34 * M[22] + Mstmp86 * M[10] +
            M[65];
#pragma omp atomic
  Ms[66] += Mstmp11 * Mstmp87 + Mstmp11 * Mstmp88 + Mstmp11 * Mstmp89 + Mstmp11 * M[30] +
            Mstmp122 * z + Mstmp123 * z + Mstmp124 * z + Mstmp125 * Mstmp7 +
            Mstmp125 * M[7] + Mstmp15 * Mstmp86 + Mstmp16 * Mstmp86 + Mstmp17 * Mstmp86 +
            Mstmp218 * M[2] + Mstmp219 + Mstmp220 + Mstmp221 + Mstmp34 * Mstmp55 +
            Mstmp34 * Mstmp57 + Mstmp34 * Mstmp59 + Mstmp34 * M[23] + Mstmp86 * M[11] +
            z * M[40] + M[66];
#pragma omp atomic
  Ms[67] += Mstmp0 * Mstmp177 + Mstmp11 * Mstmp90 + Mstmp11 * Mstmp91 +
            Mstmp11 * Mstmp92 + Mstmp11 * M[31] + Mstmp121 * Mstmp5 + Mstmp121 * M[6] +
            Mstmp125 * Mstmp8 + Mstmp125 * M[8] + Mstmp126 * z + Mstmp127 * z +
            Mstmp129 * z + Mstmp177 * M[3] + Mstmp222 + Mstmp223 + Mstmp224 +
            Mstmp34 * Mstmp61 + Mstmp34 * Mstmp62 + Mstmp34 * Mstmp63 + Mstmp34 * M[24] +
            Mstmp41 * Mstmp52 + Mstmp41 * Mstmp53 + Mstmp41 * Mstmp54 + Mstmp41 * M[22] +
            z * M[41] + M[67];
#pragma omp atomic
  Ms[68] += Mstmp11 * Mstmp93 + Mstmp11 * Mstmp94 + Mstmp11 * Mstmp95 + Mstmp11 * M[32] +
            Mstmp12 * Mstmp96 + Mstmp121 * Mstmp6 + Mstmp121 * M[7] + Mstmp13 * Mstmp96 +
            Mstmp132 * z + Mstmp133 * z + Mstmp135 * z + Mstmp14 * Mstmp96 +
            Mstmp214 * M[1] + Mstmp225 + Mstmp226 + Mstmp227 + Mstmp41 * Mstmp55 +
            Mstmp41 * Mstmp56 + Mstmp41 * Mstmp58 + Mstmp41 * M[23] + Mstmp96 * M[10] +
            z * M[42] + M[68];
#pragma omp atomic
  Ms[69] += Mstmp0 * Mstmp184 + Mstmp11 * Mstmp97 + Mstmp11 * M[33] + Mstmp121 * M[8] +
            Mstmp138 * z + Mstmp15 * Mstmp96 + Mstmp184 * M[3] + Mstmp214 * M[2] +
            Mstmp228 + Mstmp41 * Mstmp61 + Mstmp41 * M[24] + Mstmp96 * M[11] + z * M[43] +
            M[69];
#pragma omp atomic
  Ms[70] += Mstmp1 * Mstmp170 + Mstmp170 * M[4] + Mstmp18 * Mstmp86 + Mstmp229 +
            Mstmp230 + Mstmp231 + Mstmp232 * Mstmp233 + Mstmp34 * Mstmp64 +
            Mstmp34 * M[25] + Mstmp86 * M[12] + M[70];
#pragma omp atomic
  Ms[71] += Mstmp141 * z + Mstmp142 * z + Mstmp143 * z + Mstmp170 * Mstmp3 +
            Mstmp170 * Mstmp4 + Mstmp170 * M[5] + Mstmp23 * Mstmp86 + Mstmp234 +
            Mstmp235 + Mstmp236 + Mstmp25 * Mstmp86 + Mstmp27 * Mstmp86 +
            Mstmp34 * Mstmp69 + Mstmp34 * Mstmp71 + Mstmp34 * Mstmp73 + Mstmp34 * M[26] +
            Mstmp86 * M[13] + z * M[44] + M[71];
#pragma omp atomic
  Ms[72] += Mstmp1 * Mstmp177 + Mstmp146 * z + Mstmp147 * z + Mstmp149 * z +
            Mstmp177 * M[4] + Mstmp205 * Mstmp240 + Mstmp237 + Mstmp238 + Mstmp239 +
            Mstmp29 * Mstmp86 + Mstmp30 * Mstmp86 + Mstmp31 * Mstmp86 +
            Mstmp34 * Mstmp75 + Mstmp34 * Mstmp77 + Mstmp34 * Mstmp79 + Mstmp34 * M[27] +
            Mstmp41 * Mstmp64 + Mstmp41 * Mstmp65 + Mstmp41 * Mstmp66 + Mstmp41 * M[25] +
            Mstmp86 * M[14] + z * M[45] + M[72];
#pragma omp atomic
  Ms[73] += Mstmp111 * Mstmp244 + Mstmp152 * z + Mstmp153 * z + Mstmp155 * z +
            Mstmp177 * Mstmp3 + Mstmp177 * M[5] + Mstmp18 * Mstmp96 + Mstmp19 * Mstmp96 +
            Mstmp20 * Mstmp96 + Mstmp241 + Mstmp242 + Mstmp243 + Mstmp34 * Mstmp81 +
            Mstmp34 * Mstmp82 + Mstmp34 * Mstmp83 + Mstmp34 * M[28] + Mstmp41 * Mstmp69 +
            Mstmp41 * Mstmp70 + Mstmp41 * Mstmp72 + Mstmp41 * M[26] + Mstmp96 * M[12] +
            z * M[46] + M[73];
#pragma omp atomic
  Ms[74] += Mstmp1 * Mstmp184 + Mstmp159 * z + Mstmp160 * z + Mstmp162 * z +
            Mstmp184 * Mstmp2 + Mstmp184 * M[4] + Mstmp23 * Mstmp96 + Mstmp24 * Mstmp96 +
            Mstmp245 + Mstmp246 + Mstmp247 + Mstmp26 * Mstmp96 + Mstmp41 * Mstmp75 +
            Mstmp41 * Mstmp76 + Mstmp41 * Mstmp78 + Mstmp41 * M[27] + Mstmp96 * M[13] +
            z * M[47] + M[74];
#pragma omp atomic
  Ms[75] += Mstmp165 * z + Mstmp184 * Mstmp3 + Mstmp184 * M[5] + Mstmp233 * Mstmp249 +
            Mstmp248 + Mstmp29 * Mstmp96 + Mstmp41 * Mstmp81 + Mstmp41 * M[28] +
            Mstmp96 * M[14] + z * M[48] + M[75];
#pragma omp atomic
  Ms[76] += Mstmp170 * M[6] + Mstmp250 + Mstmp251 * M[1] + Mstmp34 * M[29] +
            Mstmp86 * M[15] + M[76];
#pragma omp atomic
  Ms[77] += Mstmp169 * z + Mstmp170 * Mstmp7 + Mstmp170 * M[7] + Mstmp251 * M[2] +
            Mstmp252 + Mstmp34 * Mstmp88 + Mstmp34 * M[30] + Mstmp36 * Mstmp86 +
            Mstmp86 * M[16] + z * M[49] + M[77];
#pragma omp atomic
  Ms[78] += Mstmp170 * Mstmp8 + Mstmp170 * M[8] + Mstmp171 * z + Mstmp177 * M[6] +
            Mstmp253 + Mstmp254 * M[1] + Mstmp34 * Mstmp91 + Mstmp34 * M[31] +
            Mstmp39 * Mstmp86 + Mstmp41 * Mstmp85 + Mstmp41 * M[29] + Mstmp86 * M[17] +
            z * M[50] + M[78];
#pragma omp atomic
  Ms[79] += Mstmp174 * z + Mstmp177 * M[7] + Mstmp254 * M[2] + Mstmp255 +
            Mstmp256 * M[1] + Mstmp33 * Mstmp96 + Mstmp34 * Mstmp94 + Mstmp34 * M[32] +
            Mstmp41 * Mstmp87 + Mstmp41 * M[30] + Mstmp42 * Mstmp86 + Mstmp86 * M[18] +
            Mstmp96 * M[15] + z * M[51] + M[79];
#pragma omp atomic
  Ms[80] += Mstmp177 * M[8] + Mstmp178 * z + Mstmp184 * Mstmp5 + Mstmp184 * M[6] +
            Mstmp256 * M[2] + Mstmp257 + Mstmp34 * Mstmp97 + Mstmp34 * M[33] +
            Mstmp35 * Mstmp96 + Mstmp41 * Mstmp90 + Mstmp41 * M[31] + Mstmp96 * M[16] +
            z * M[52] + M[80];
#pragma omp atomic
  Ms[81] += Mstmp181 * z + Mstmp184 * Mstmp6 + Mstmp184 * M[7] + Mstmp258 +
            Mstmp259 * M[1] + Mstmp38 * Mstmp96 + Mstmp41 * Mstmp93 + Mstmp41 * M[32] +
            Mstmp96 * M[17] + z * M[53] + M[81];
#pragma omp atomic
  Ms[82] += Mstmp184 * M[8] + Mstmp259 * M[2] + Mstmp41 * M[33] + Mstmp96 * M[18] +
            z * M[54] + M[82];
#pragma omp atomic
  Ms[83] += Mstmp100 * M[9] + Mstmp11 * M[34] + Mstmp187 * M[3] + Mstmp260 * M[0] +
            Mstmp45 * M[19] + x * M[55] + M[83];
#pragma omp atomic
  Ms[84] += Mstmp100 * Mstmp13 + Mstmp100 * M[10] + Mstmp102 * Mstmp11 + Mstmp11 * M[35] +
            Mstmp186 * y + Mstmp187 * Mstmp2 + Mstmp187 * M[4] + Mstmp260 * M[1] +
            Mstmp45 * Mstmp47 + Mstmp45 * M[20] + x * M[56] + y * M[55] + M[84];
#pragma omp atomic
  Ms[85] += Mstmp100 * Mstmp16 + Mstmp100 * M[11] + Mstmp105 * Mstmp11 + Mstmp11 * M[36] +
            Mstmp186 * z + Mstmp187 * Mstmp4 + Mstmp187 * M[5] + Mstmp260 * M[2] +
            Mstmp45 * Mstmp50 + Mstmp45 * M[21] + x * M[57] + z * M[55] + M[85];
#pragma omp atomic
  Ms[86] += Mstmp100 * Mstmp19 + Mstmp100 * M[12] + Mstmp108 * Mstmp11 + Mstmp11 * M[37] +
            Mstmp111 * Mstmp261 + Mstmp125 * M[9] + Mstmp187 * Mstmp5 + Mstmp187 * M[6] +
            Mstmp188 * y + Mstmp206 * M[3] + Mstmp34 * Mstmp98 + Mstmp34 * M[34] +
            Mstmp45 * Mstmp53 + Mstmp45 * M[22] + x * M[58] + y * M[56] + M[86];
#pragma omp atomic
  Ms[87] += Mstmp100 * Mstmp24 + Mstmp100 * Mstmp25 + Mstmp100 * Mstmp28 +
            Mstmp100 * M[13] + Mstmp11 * Mstmp113 + Mstmp11 * Mstmp114 +
            Mstmp11 * Mstmp117 + Mstmp11 * M[38] + Mstmp187 * Mstmp6 + Mstmp187 * Mstmp7 +
            Mstmp187 * M[7] + Mstmp188 * z + Mstmp189 * z + Mstmp190 * z + Mstmp191 * y +
            Mstmp45 * Mstmp56 + Mstmp45 * Mstmp57 + Mstmp45 * Mstmp60 + Mstmp45 * M[23] +
            x * M[59] + y * M[57] + z * M[56] + M[87];
#pragma omp atomic
  Ms[88] += Mstmp100 * Mstmp30 + Mstmp100 * M[14] + Mstmp11 * Mstmp119 + Mstmp11 * M[39] +
            Mstmp121 * M[9] + Mstmp187 * Mstmp8 + Mstmp187 * M[8] + Mstmp191 * z +
            Mstmp200 * M[3] + Mstmp262 * M[0] + Mstmp41 * Mstmp98 + Mstmp41 * M[34] +
            Mstmp45 * Mstmp62 + Mstmp45 * M[24] + x * M[60] + z * M[57] + M[88];
#pragma omp atomic
  Ms[89] += Mstmp100 * Mstmp33 + Mstmp100 * M[15] + Mstmp101 * Mstmp34 +
            Mstmp11 * Mstmp123 + Mstmp11 * M[40] + Mstmp125 * M[10] + Mstmp192 * y +
            Mstmp205 * Mstmp263 + Mstmp206 * M[4] + Mstmp218 * M[3] + Mstmp264 * M[1] +
            Mstmp34 * M[35] + Mstmp43 * Mstmp86 + Mstmp45 * Mstmp65 + Mstmp45 * M[25] +
            Mstmp86 * M[19] + x * M[61] + y * M[58] + M[89];
#pragma omp atomic
  Ms[90] += Mstmp100 * Mstmp35 + Mstmp100 * Mstmp36 + Mstmp100 * Mstmp37 +
            Mstmp100 * M[16] + Mstmp104 * Mstmp34 + Mstmp105 * Mstmp34 +
            Mstmp106 * Mstmp34 + Mstmp11 * Mstmp127 + Mstmp11 * Mstmp128 +
            Mstmp11 * Mstmp131 + Mstmp11 * M[41] + Mstmp125 * Mstmp16 + Mstmp125 * M[11] +
            Mstmp192 * z + Mstmp193 * z + Mstmp194 * z + Mstmp196 * y +
            Mstmp206 * Mstmp4 + Mstmp206 * M[5] + Mstmp264 * M[2] + Mstmp34 * M[36] +
            Mstmp45 * Mstmp70 + Mstmp45 * Mstmp71 + Mstmp45 * Mstmp74 + Mstmp45 * M[26] +
            x * M[62] + y * M[59] + z * M[58] + M[90];
#pragma omp atomic
  Ms[91] += Mstmp100 * Mstmp38 + Mstmp100 * Mstmp39 + Mstmp100 * Mstmp40 +
            Mstmp100 * M[17] + Mstmp101 * Mstmp41 + Mstmp102 * Mstmp41 +
            Mstmp103 * Mstmp41 + Mstmp11 * Mstmp133 + Mstmp11 * Mstmp134 +
            Mstmp11 * Mstmp137 + Mstmp11 * M[42] + Mstmp121 * Mstmp13 + Mstmp121 * M[10] +
            Mstmp196 * z + Mstmp197 * z + Mstmp198 * z + Mstmp199 * y +
            Mstmp2 * Mstmp200 + Mstmp200 * M[4] + Mstmp262 * M[1] + Mstmp41 * M[35] +
            Mstmp45 * Mstmp76 + Mstmp45 * Mstmp77 + Mstmp45 * Mstmp80 + Mstmp45 * M[27] +
            x * M[63] + y * M[60] + z * M[59] + M[91];
#pragma omp atomic
  Ms[92] += Mstmp100 * Mstmp42 + Mstmp100 * M[18] + Mstmp104 * Mstmp41 +
            Mstmp11 * Mstmp139 + Mstmp11 * M[43] + Mstmp121 * M[11] + Mstmp199 * z +
            Mstmp200 * M[5] + Mstmp214 * M[3] + Mstmp262 * M[2] + Mstmp265 * M[0] +
            Mstmp41 * M[36] + Mstmp43 * Mstmp96 + Mstmp45 * Mstmp82 + Mstmp45 * M[28] +
            Mstmp96 * M[19] + x * M[64] + z * M[60] + M[92];
#pragma omp atomic
  Ms[93] += Mstmp107 * Mstmp34 + Mstmp11 * Mstmp142 + Mstmp11 * M[44] + Mstmp125 * M[12] +
            Mstmp144 * Mstmp267 + Mstmp170 * Mstmp9 + Mstmp170 * M[9] + Mstmp201 * y +
            Mstmp206 * M[6] + Mstmp218 * M[4] + Mstmp268 * M[1] + Mstmp34 * M[37] +
            Mstmp45 * Mstmp85 + Mstmp45 * M[29] + Mstmp46 * Mstmp86 + Mstmp86 * M[20] +
            x * M[65] + y * M[61] + M[93];
#pragma omp atomic
  Ms[94] += Mstmp11 * Mstmp147 + Mstmp11 * Mstmp148 + Mstmp11 * Mstmp151 +
            Mstmp11 * M[45] + Mstmp112 * Mstmp34 + Mstmp114 * Mstmp34 +
            Mstmp116 * Mstmp34 + Mstmp125 * Mstmp25 + Mstmp125 * M[13] + Mstmp201 * z +
            Mstmp202 * z + Mstmp203 * z + Mstmp206 * Mstmp7 + Mstmp206 * M[7] +
            Mstmp207 * y + Mstmp218 * Mstmp4 + Mstmp218 * M[5] + Mstmp268 * M[2] +
            Mstmp34 * M[38] + Mstmp45 * Mstmp87 + Mstmp45 * Mstmp88 + Mstmp45 * Mstmp89 +
            Mstmp45 * M[30] + Mstmp49 * Mstmp86 + Mstmp50 * Mstmp86 + Mstmp51 * Mstmp86 +
            Mstmp86 * M[21] + x * M[66] + y * M[62] + z * M[61] + M[94];
#pragma omp atomic
  Ms[95] += Mstmp107 * Mstmp41 + Mstmp108 * Mstmp41 + Mstmp109 * Mstmp41 +
            Mstmp11 * Mstmp153 + Mstmp11 * Mstmp154 + Mstmp11 * Mstmp157 +
            Mstmp11 * M[46] + Mstmp111 * Mstmp269 + Mstmp118 * Mstmp34 +
            Mstmp119 * Mstmp34 + Mstmp120 * Mstmp34 + Mstmp121 * Mstmp19 +
            Mstmp121 * M[12] + Mstmp125 * Mstmp30 + Mstmp125 * M[14] + Mstmp177 * Mstmp9 +
            Mstmp177 * M[9] + Mstmp200 * Mstmp5 + Mstmp200 * M[6] + Mstmp206 * Mstmp8 +
            Mstmp206 * M[8] + Mstmp207 * z + Mstmp208 * z + Mstmp209 * z + Mstmp210 * y +
            Mstmp34 * M[39] + Mstmp41 * M[37] + Mstmp45 * Mstmp90 + Mstmp45 * Mstmp91 +
            Mstmp45 * Mstmp92 + Mstmp45 * M[31] + x * M[67] + y * M[63] + z * M[62] +
            M[95];
#pragma omp atomic
  Ms[96] += Mstmp11 * Mstmp160 + Mstmp11 * Mstmp161 + Mstmp11 * Mstmp164 +
            Mstmp11 * M[47] + Mstmp112 * Mstmp41 + Mstmp113 * Mstmp41 +
            Mstmp115 * Mstmp41 + Mstmp121 * Mstmp24 + Mstmp121 * M[13] +
            Mstmp2 * Mstmp214 + Mstmp200 * Mstmp6 + Mstmp200 * M[7] + Mstmp210 * z +
            Mstmp211 * z + Mstmp212 * z + Mstmp213 * y + Mstmp214 * M[4] +
            Mstmp265 * M[1] + Mstmp41 * M[38] + Mstmp45 * Mstmp93 + Mstmp45 * Mstmp94 +
            Mstmp45 * Mstmp95 + Mstmp45 * M[32] + Mstmp46 * Mstmp96 + Mstmp47 * Mstmp96 +
            Mstmp48 * Mstmp96 + Mstmp96 * M[20] + x * M[68] + y * M[64] + z * M[63] +
            M[96];
#pragma omp atomic
  Ms[97] += Mstmp11 * Mstmp166 + Mstmp11 * M[48] + Mstmp118 * Mstmp41 + Mstmp121 * M[14] +
            Mstmp168 * Mstmp267 + Mstmp184 * Mstmp9 + Mstmp184 * M[9] + Mstmp200 * M[8] +
            Mstmp213 * z + Mstmp214 * M[5] + Mstmp265 * M[2] + Mstmp41 * M[39] +
            Mstmp45 * Mstmp97 + Mstmp45 * M[33] + Mstmp49 * Mstmp96 + Mstmp96 * M[21] +
            x * M[69] + z * M[64] + M[97];
#pragma omp atomic
  Ms[98] += Mstmp0 * Mstmp251 + Mstmp11 * Mstmp169 + Mstmp11 * M[49] +
            Mstmp12 * Mstmp170 + Mstmp122 * Mstmp34 + Mstmp125 * M[15] +
            Mstmp170 * M[10] + Mstmp215 * y + Mstmp218 * M[6] + Mstmp251 * M[3] +
            Mstmp270 * M[1] + Mstmp34 * M[40] + Mstmp52 * Mstmp86 + Mstmp86 * M[22] +
            x * M[70] + y * M[65] + M[98];
#pragma omp atomic
  Ms[99] += Mstmp11 * Mstmp171 + Mstmp11 * Mstmp172 + Mstmp11 * Mstmp173 +
            Mstmp11 * M[50] + Mstmp125 * Mstmp36 + Mstmp125 * M[16] + Mstmp126 * Mstmp34 +
            Mstmp128 * Mstmp34 + Mstmp130 * Mstmp34 + Mstmp15 * Mstmp170 +
            Mstmp16 * Mstmp170 + Mstmp17 * Mstmp170 + Mstmp170 * M[11] + Mstmp215 * z +
            Mstmp216 * z + Mstmp217 * z + Mstmp218 * Mstmp7 + Mstmp218 * M[7] +
            Mstmp219 * y + Mstmp270 * M[2] + Mstmp34 * M[41] + Mstmp55 * Mstmp86 +
            Mstmp57 * Mstmp86 + Mstmp59 * Mstmp86 + Mstmp86 * M[23] + x * M[71] +
            y * M[66] + z * M[65] + M[99];
#pragma omp atomic
  Ms[100] += Mstmp0 * Mstmp254 + Mstmp11 * Mstmp174 + Mstmp11 * Mstmp175 +
             Mstmp11 * Mstmp176 + Mstmp11 * M[51] + Mstmp12 * Mstmp177 +
             Mstmp121 * Mstmp33 + Mstmp121 * M[15] + Mstmp122 * Mstmp41 +
             Mstmp123 * Mstmp41 + Mstmp124 * Mstmp41 + Mstmp125 * Mstmp39 +
             Mstmp125 * M[17] + Mstmp132 * Mstmp34 + Mstmp134 * Mstmp34 +
             Mstmp136 * Mstmp34 + Mstmp177 * M[10] + Mstmp218 * Mstmp8 + Mstmp218 * M[8] +
             Mstmp219 * z + Mstmp220 * z + Mstmp221 * z + Mstmp222 * y + Mstmp254 * M[3] +
             Mstmp271 * M[1] + Mstmp34 * M[42] + Mstmp41 * M[40] + Mstmp61 * Mstmp86 +
             Mstmp62 * Mstmp86 + Mstmp63 * Mstmp86 + Mstmp86 * M[24] + x * M[72] +
             y * M[67] + z * M[66] + M[100];
#pragma omp atomic
  Ms[101] += Mstmp0 * Mstmp256 + Mstmp11 * Mstmp178 + Mstmp11 * Mstmp179 +
             Mstmp11 * Mstmp180 + Mstmp11 * M[52] + Mstmp121 * Mstmp35 +
             Mstmp121 * M[16] + Mstmp125 * Mstmp42 + Mstmp125 * M[18] +
             Mstmp126 * Mstmp41 + Mstmp127 * Mstmp41 + Mstmp129 * Mstmp41 +
             Mstmp138 * Mstmp34 + Mstmp139 * Mstmp34 + Mstmp140 * Mstmp34 +
             Mstmp15 * Mstmp177 + Mstmp177 * M[11] + Mstmp214 * Mstmp5 + Mstmp214 * M[6] +
             Mstmp222 * z + Mstmp223 * z + Mstmp224 * z + Mstmp225 * y + Mstmp256 * M[3] +
             Mstmp271 * M[2] + Mstmp34 * M[43] + Mstmp41 * M[41] + Mstmp52 * Mstmp96 +
             Mstmp53 * Mstmp96 + Mstmp54 * Mstmp96 + Mstmp96 * M[22] + x * M[73] +
             y * M[68] + z * M[67] + M[101];
#pragma omp atomic
  Ms[102] += Mstmp11 * Mstmp181 + Mstmp11 * Mstmp182 + Mstmp11 * Mstmp183 +
             Mstmp11 * M[53] + Mstmp12 * Mstmp184 + Mstmp121 * Mstmp38 +
             Mstmp121 * M[17] + Mstmp13 * Mstmp184 + Mstmp132 * Mstmp41 +
             Mstmp133 * Mstmp41 + Mstmp135 * Mstmp41 + Mstmp14 * Mstmp184 +
             Mstmp184 * M[10] + Mstmp214 * Mstmp6 + Mstmp214 * M[7] + Mstmp225 * z +
             Mstmp226 * z + Mstmp227 * z + Mstmp228 * y + Mstmp272 * M[1] +
             Mstmp41 * M[42] + Mstmp55 * Mstmp96 + Mstmp56 * Mstmp96 + Mstmp58 * Mstmp96 +
             Mstmp96 * M[23] + x * M[74] + y * M[69] + z * M[68] + M[102];
#pragma omp atomic
  Ms[103] += Mstmp0 * Mstmp259 + Mstmp11 * Mstmp185 + Mstmp11 * M[54] + Mstmp121 * M[18] +
             Mstmp138 * Mstmp41 + Mstmp15 * Mstmp184 + Mstmp184 * M[11] +
             Mstmp214 * M[8] + Mstmp228 * z + Mstmp259 * M[3] + Mstmp272 * M[2] +
             Mstmp41 * M[43] + Mstmp61 * Mstmp96 + Mstmp96 * M[24] + x * M[75] +
             z * M[69] + M[103];
#pragma omp atomic
  Ms[104] += Mstmp1 * Mstmp251 + Mstmp141 * Mstmp34 + Mstmp170 * Mstmp18 +
             Mstmp170 * M[12] + Mstmp229 * y + Mstmp251 * M[4] + Mstmp273 * Mstmp274 +
             Mstmp34 * M[44] + Mstmp64 * Mstmp86 + Mstmp86 * M[25] + x * M[76] +
             y * M[70] + M[104];
#pragma omp atomic
  Ms[105] += Mstmp146 * Mstmp34 + Mstmp148 * Mstmp34 + Mstmp150 * Mstmp34 +
             Mstmp170 * Mstmp23 + Mstmp170 * Mstmp25 + Mstmp170 * Mstmp27 +
             Mstmp170 * M[13] + Mstmp229 * z + Mstmp230 * z + Mstmp231 * z +
             Mstmp234 * y + Mstmp251 * Mstmp3 + Mstmp251 * Mstmp4 + Mstmp251 * M[5] +
             Mstmp34 * M[45] + Mstmp69 * Mstmp86 + Mstmp71 * Mstmp86 + Mstmp73 * Mstmp86 +
             Mstmp86 * M[26] + x * M[77] + y * M[71] + z * M[70] + M[105];
#pragma omp atomic
  Ms[106] += Mstmp1 * Mstmp254 + Mstmp141 * Mstmp41 + Mstmp142 * Mstmp41 +
             Mstmp143 * Mstmp41 + Mstmp152 * Mstmp34 + Mstmp154 * Mstmp34 +
             Mstmp156 * Mstmp34 + Mstmp170 * Mstmp29 + Mstmp170 * Mstmp30 +
             Mstmp170 * Mstmp31 + Mstmp170 * M[14] + Mstmp177 * Mstmp18 +
             Mstmp177 * M[12] + Mstmp234 * z + Mstmp235 * z + Mstmp236 * z +
             Mstmp237 * y + Mstmp254 * M[4] + Mstmp275 * M[0] + Mstmp34 * M[46] +
             Mstmp41 * M[44] + Mstmp75 * Mstmp86 + Mstmp77 * Mstmp86 + Mstmp79 * Mstmp86 +
             Mstmp86 * M[27] + x * M[78] + y * M[72] + z * M[71] + M[106];
#pragma omp atomic
  Ms[107] += Mstmp1 * Mstmp256 + Mstmp146 * Mstmp41 + Mstmp147 * Mstmp41 +
             Mstmp149 * Mstmp41 + Mstmp159 * Mstmp34 + Mstmp161 * Mstmp34 +
             Mstmp163 * Mstmp34 + Mstmp177 * Mstmp23 + Mstmp177 * M[13] +
             Mstmp205 * Mstmp276 + Mstmp237 * z + Mstmp238 * z + Mstmp239 * z +
             Mstmp241 * y + Mstmp254 * Mstmp3 + Mstmp254 * M[5] + Mstmp256 * M[4] +
             Mstmp34 * M[47] + Mstmp41 * M[45] + Mstmp64 * Mstmp96 + Mstmp65 * Mstmp96 +
             Mstmp66 * Mstmp96 + Mstmp81 * Mstmp86 + Mstmp82 * Mstmp86 +
             Mstmp83 * Mstmp86 + Mstmp86 * M[28] + Mstmp96 * M[25] + x * M[79] +
             y * M[73] + z * M[72] + M[107];
#pragma omp atomic
  Ms[108] += Mstmp111 * Mstmp277 + Mstmp152 * Mstmp41 + Mstmp153 * Mstmp41 +
             Mstmp155 * Mstmp41 + Mstmp165 * Mstmp34 + Mstmp166 * Mstmp34 +
             Mstmp167 * Mstmp34 + Mstmp177 * Mstmp29 + Mstmp177 * M[14] +
             Mstmp18 * Mstmp184 + Mstmp184 * Mstmp19 + Mstmp184 * Mstmp20 +
             Mstmp184 * M[12] + Mstmp241 * z + Mstmp242 * z + Mstmp243 * z +
             Mstmp245 * y + Mstmp256 * Mstmp3 + Mstmp256 * M[5] + Mstmp34 * M[48] +
             Mstmp41 * M[46] + Mstmp69 * Mstmp96 + Mstmp70 * Mstmp96 + Mstmp72 * Mstmp96 +
             Mstmp96 * M[26] + x * M[80] + y * M[74] + z * M[73] + M[108];
#pragma omp atomic
  Ms[109] += Mstmp1 * Mstmp259 + Mstmp159 * Mstmp41 + Mstmp160 * Mstmp41 +
             Mstmp162 * Mstmp41 + Mstmp184 * Mstmp23 + Mstmp184 * Mstmp24 +
             Mstmp184 * Mstmp26 + Mstmp184 * M[13] + Mstmp2 * Mstmp259 + Mstmp245 * z +
             Mstmp246 * z + Mstmp247 * z + Mstmp248 * y + Mstmp259 * M[4] +
             Mstmp41 * M[47] + Mstmp75 * Mstmp96 + Mstmp76 * Mstmp96 + Mstmp78 * Mstmp96 +
             Mstmp96 * M[27] + x * M[81] + y * M[75] + z * M[74] + M[109];
#pragma omp atomic
  Ms[110] += Mstmp165 * Mstmp41 + Mstmp184 * Mstmp29 + Mstmp184 * M[14] + Mstmp248 * z +
             Mstmp259 * Mstmp3 + Mstmp259 * M[5] + Mstmp274 * Mstmp278 + Mstmp41 * M[48] +
             Mstmp81 * Mstmp96 + Mstmp96 * M[28] + x * M[82] + z * M[75] + M[110];
#pragma omp atomic
  Ms[111] += Mstmp170 * M[15] + Mstmp251 * M[6] + Mstmp279 * M[1] + Mstmp34 * M[49] +
             Mstmp86 * M[29] + y * M[76] + M[111];
#pragma omp atomic
  Ms[112] += Mstmp170 * Mstmp36 + Mstmp170 * M[16] + Mstmp172 * Mstmp34 + Mstmp250 * z +
             Mstmp251 * Mstmp7 + Mstmp251 * M[7] + Mstmp279 * M[2] + Mstmp34 * M[50] +
             Mstmp86 * Mstmp88 + Mstmp86 * M[30] + y * M[77] + z * M[76] + M[112];
#pragma omp atomic
  Ms[113] += Mstmp169 * Mstmp41 + Mstmp170 * Mstmp39 + Mstmp170 * M[17] +
             Mstmp175 * Mstmp34 + Mstmp177 * M[15] + Mstmp251 * Mstmp8 + Mstmp251 * M[8] +
             Mstmp252 * z + Mstmp254 * M[6] + Mstmp275 * M[1] + Mstmp34 * M[51] +
             Mstmp41 * M[49] + Mstmp86 * Mstmp91 + Mstmp86 * M[31] + y * M[78] +
             z * M[77] + M[113];
#pragma omp atomic
  Ms[114] += Mstmp170 * Mstmp42 + Mstmp170 * M[18] + Mstmp171 * Mstmp41 +
             Mstmp177 * M[16] + Mstmp179 * Mstmp34 + Mstmp253 * z + Mstmp254 * M[7] +
             Mstmp256 * M[6] + Mstmp275 * M[2] + Mstmp280 * M[1] + Mstmp34 * M[52] +
             Mstmp41 * M[50] + Mstmp85 * Mstmp96 + Mstmp86 * Mstmp94 + Mstmp86 * M[32] +
             Mstmp96 * M[29] + y * M[79] + z * M[78] + M[114];
#pragma omp atomic
  Ms[115] += Mstmp174 * Mstmp41 + Mstmp177 * M[17] + Mstmp182 * Mstmp34 +
             Mstmp184 * Mstmp33 + Mstmp184 * M[15] + Mstmp254 * M[8] + Mstmp255 * z +
             Mstmp256 * M[7] + Mstmp280 * M[2] + Mstmp281 * M[1] + Mstmp34 * M[53] +
             Mstmp41 * M[51] + Mstmp86 * Mstmp97 + Mstmp86 * M[33] + Mstmp87 * Mstmp96 +
             Mstmp96 * M[30] + y * M[80] + z * M[79] + M[115];
#pragma omp atomic
  Ms[116] += Mstmp177 * M[18] + Mstmp178 * Mstmp41 + Mstmp184 * Mstmp35 +
             Mstmp184 * M[16] + Mstmp185 * Mstmp34 + Mstmp256 * M[8] + Mstmp257 * z +
             Mstmp259 * Mstmp5 + Mstmp259 * M[6] + Mstmp281 * M[2] + Mstmp34 * M[54] +
             Mstmp41 * M[52] + Mstmp90 * Mstmp96 + Mstmp96 * M[31] + y * M[81] +
             z * M[80] + M[116];
#pragma omp atomic
  Ms[117] += Mstmp181 * Mstmp41 + Mstmp184 * Mstmp38 + Mstmp184 * M[17] + Mstmp258 * z +
             Mstmp259 * Mstmp6 + Mstmp259 * M[7] + Mstmp282 * M[1] + Mstmp41 * M[53] +
             Mstmp93 * Mstmp96 + Mstmp96 * M[32] + y * M[82] + z * M[81] + M[117];
#pragma omp atomic
  Ms[118] += Mstmp184 * M[18] + Mstmp259 * M[8] + Mstmp282 * M[2] + Mstmp41 * M[54] +
             Mstmp96 * M[33] + z * M[82] + M[118];
}

void field_m1_M2L_7(double x, double y, double z, double* M, double* L) {
  double R = sqrt(x * x + y * y + z * z);
  double D[119];
  double Dtmp0   = 1.0 * pow(R, -3.0);
  double Dtmp1   = -Dtmp0;
  double Dtmp2   = (x * x);
  double Dtmp3   = pow(R, -5.0);
  double Dtmp4   = 3.0 * Dtmp3;
  double Dtmp5   = x * y;
  double Dtmp6   = x * z;
  double Dtmp7   = (y * y);
  double Dtmp8   = y * z;
  double Dtmp9   = 9.0 * Dtmp3;
  double Dtmp10  = -Dtmp9;
  double Dtmp11  = pow(R, -7.0);
  double Dtmp12  = 15.0 * Dtmp11;
  double Dtmp13  = Dtmp12 * Dtmp2;
  double Dtmp14  = -Dtmp4;
  double Dtmp15  = Dtmp13 + Dtmp14;
  double Dtmp16  = Dtmp12 * Dtmp7;
  double Dtmp17  = Dtmp14 + Dtmp16;
  double Dtmp18  = 1.0 * x;
  double Dtmp19  = Dtmp8 * x;
  double Dtmp20  = (x * x * x * x);
  double Dtmp21  = pow(R, -9.0);
  double Dtmp22  = 105.0 * Dtmp21;
  double Dtmp23  = 90.0 * Dtmp11;
  double Dtmp24  = 45.0 * Dtmp11;
  double Dtmp25  = -Dtmp24;
  double Dtmp26  = Dtmp2 * Dtmp22;
  double Dtmp27  = x * (Dtmp25 + Dtmp26);
  double Dtmp28  = -Dtmp12;
  double Dtmp29  = Dtmp22 * Dtmp7;
  double Dtmp30  = Dtmp25 + Dtmp29;
  double Dtmp31  = Dtmp18 * y;
  double Dtmp32  = Dtmp18 * z;
  double Dtmp33  = (y * y * y * y);
  double Dtmp34  = 225.0 * Dtmp11;
  double Dtmp35  = pow(R, -11.0);
  double Dtmp36  = 945.0 * Dtmp35;
  double Dtmp37  = Dtmp20 * Dtmp36;
  double Dtmp38  = Dtmp2 * Dtmp21;
  double Dtmp39  = 630.0 * Dtmp38;
  double Dtmp40  = Dtmp24 + Dtmp37 - Dtmp39;
  double Dtmp41  = -Dtmp26;
  double Dtmp42  = 315.0 * Dtmp21;
  double Dtmp43  = Dtmp42 * Dtmp7;
  double Dtmp44  = Dtmp2 * Dtmp36;
  double Dtmp45  = Dtmp44 * Dtmp7;
  double Dtmp46  = Dtmp24 + Dtmp45;
  double Dtmp47  = -Dtmp42;
  double Dtmp48  = -Dtmp29;
  double Dtmp49  = Dtmp2 * Dtmp42;
  double Dtmp50  = Dtmp33 * Dtmp36;
  double Dtmp51  = Dtmp21 * Dtmp7;
  double Dtmp52  = 630.0 * Dtmp51;
  double Dtmp53  = Dtmp24 + Dtmp50 - Dtmp52;
  double Dtmp54  = Dtmp36 * Dtmp7;
  double Dtmp55  = Dtmp18 * Dtmp8;
  double Dtmp56  = -Dtmp34;
  double Dtmp57  = (x * x * x * x * x * x);
  double Dtmp58  = pow(R, -13.0);
  double Dtmp59  = 10395.0 * Dtmp58;
  double Dtmp60  = 14175.0 * Dtmp35;
  double Dtmp61  = 1575.0 * Dtmp21;
  double Dtmp62  = Dtmp20 * Dtmp59;
  double Dtmp63  = Dtmp2 * Dtmp35;
  double Dtmp64  = 9450.0 * Dtmp63;
  double Dtmp65  = x * (Dtmp61 + Dtmp62 - Dtmp64);
  double Dtmp66  = 5670.0 * Dtmp63;
  double Dtmp67  = Dtmp25 - Dtmp66 * Dtmp7;
  double Dtmp68  = 945.0 * Dtmp21;
  double Dtmp69  = 2835.0 * Dtmp63;
  double Dtmp70  = -Dtmp69;
  double Dtmp71  = Dtmp35 * Dtmp7;
  double Dtmp72  = 2835.0 * Dtmp71;
  double Dtmp73  = Dtmp2 * Dtmp7;
  double Dtmp74  = Dtmp59 * Dtmp73;
  double Dtmp75  = -Dtmp72 + Dtmp74;
  double Dtmp76  = Dtmp33 * Dtmp59;
  double Dtmp77  = 9450.0 * Dtmp71;
  double Dtmp78  = Dtmp61 + Dtmp76 - Dtmp77;
  double Dtmp79  = 5670.0 * Dtmp71;
  double Dtmp80  = (y * y * y * y * y * y);
  double Dtmp81  = -11025.0 * Dtmp21;
  double Dtmp82  = 135135.0 * pow(R, -15.0);
  double Dtmp83  = Dtmp57 * Dtmp82;
  double Dtmp84  = Dtmp20 * Dtmp58;
  double Dtmp85  = -Dtmp61;
  double Dtmp86  = 42525.0 * Dtmp63 + Dtmp83 - 155925.0 * Dtmp84 + Dtmp85;
  double Dtmp87  = Dtmp20 * Dtmp82;
  double Dtmp88  = Dtmp7 * Dtmp87;
  double Dtmp89  = -Dtmp62 + Dtmp88;
  double Dtmp90  = Dtmp2 * Dtmp58;
  double Dtmp91  = 103950.0 * Dtmp90;
  double Dtmp92  = -Dtmp7 * Dtmp91 + Dtmp85;
  double Dtmp93  = -Dtmp68;
  double Dtmp94  = -62370.0 * Dtmp7 * Dtmp90;
  double Dtmp95  = Dtmp72 + Dtmp94;
  double Dtmp96  = 31185.0 * Dtmp58;
  double Dtmp97  = Dtmp33 * Dtmp82;
  double Dtmp98  = Dtmp2 * Dtmp97;
  double Dtmp99  = Dtmp69 + Dtmp94 + Dtmp98;
  double Dtmp100 = -Dtmp76;
  double Dtmp101 = Dtmp80 * Dtmp82;
  double Dtmp102 = Dtmp33 * Dtmp58;
  double Dtmp103 = Dtmp101 - 155925.0 * Dtmp102 + 42525.0 * Dtmp71 + Dtmp85;
  D[0]           = -Dtmp0 * x;
  D[1]           = -Dtmp0 * y;
  D[2]           = -Dtmp0 * z;
  D[3]           = Dtmp1 + Dtmp2 * Dtmp4;
  D[4]           = Dtmp4 * Dtmp5;
  D[5]           = Dtmp4 * Dtmp6;
  D[6]           = Dtmp1 + Dtmp4 * Dtmp7;
  D[7]           = Dtmp4 * Dtmp8;
  D[8]           = -D[3] - D[6];
  D[9]           = -x * (Dtmp10 + Dtmp13);
  D[10]          = -Dtmp15 * y;
  D[11]          = -Dtmp15 * z;
  D[12]          = -Dtmp17 * Dtmp18;
  D[13]          = -Dtmp12 * Dtmp19;
  D[14]          = -D[9] - D[12];
  D[15]          = -y * (Dtmp10 + Dtmp16);
  D[16]          = -Dtmp17 * z;
  D[17]          = -D[10] - D[15];
  D[18]          = -D[11] - D[16];
  D[19]          = -Dtmp2 * Dtmp23 + Dtmp20 * Dtmp22 + Dtmp9;
  D[20]          = Dtmp27 * y;
  D[21]          = Dtmp27 * z;
  D[22]          = -Dtmp13 - Dtmp16 + Dtmp26 * Dtmp7 + Dtmp4;
  D[23]          = Dtmp8 * (Dtmp26 + Dtmp28);
  D[24]          = -D[19] - D[22];
  D[25]          = Dtmp30 * Dtmp31;
  D[26]          = Dtmp32 * (Dtmp28 + Dtmp29);
  D[27]          = -D[20] - D[25];
  D[28]          = -D[21] - D[26];
  D[29]          = Dtmp22 * Dtmp33 - Dtmp23 * Dtmp7 + Dtmp9;
  D[30]          = Dtmp30 * Dtmp8;
  D[31]          = -D[22] - D[29];
  D[32]          = -D[23] - D[30];
  D[33]          = -D[24] - D[31];
  D[34]          = -x * (Dtmp34 + Dtmp37 - 1050.0 * Dtmp38);
  D[35]          = -Dtmp40 * y;
  D[36]          = -Dtmp40 * z;
  D[37]          = -x * (Dtmp41 - Dtmp43 + Dtmp46);
  D[38]          = -Dtmp19 * (Dtmp44 + Dtmp47);
  D[39]          = -D[34] - D[37];
  D[40]          = -y * (Dtmp46 + Dtmp48 - Dtmp49);
  D[41]          = -z * (Dtmp12 + Dtmp41 + Dtmp45 + Dtmp48);
  D[42]          = -D[35] - D[40];
  D[43]          = -D[36] - D[41];
  D[44]          = -Dtmp18 * Dtmp53;
  D[45]          = -Dtmp55 * (Dtmp47 + Dtmp54);
  D[46]          = -D[37] - D[44];
  D[47]          = -D[38] - D[45];
  D[48]          = -D[39] - D[46];
  D[49]          = -y * (Dtmp34 + Dtmp50 - 1050.0 * Dtmp51);
  D[50]          = -Dtmp53 * z;
  D[51]          = -D[40] - D[49];
  D[52]          = -D[41] - D[50];
  D[53]          = -D[42] - D[51];
  D[54]          = -D[43] - D[52];
  D[55]          = -Dtmp20 * Dtmp60 + 4725.0 * Dtmp38 + Dtmp56 + Dtmp57 * Dtmp59;
  D[56]          = Dtmp65 * y;
  D[57]          = Dtmp65 * z;
  D[58]          = -Dtmp37 + Dtmp39 + Dtmp43 + Dtmp62 * Dtmp7 + Dtmp67;
  D[59]          = Dtmp8 * (Dtmp42 + Dtmp62 - Dtmp66);
  D[60]          = -D[55] - D[58];
  D[61]          = Dtmp5 * (Dtmp68 + Dtmp70 + Dtmp75);
  D[62]          = Dtmp6 * (Dtmp42 - Dtmp44 + Dtmp75);
  D[63]          = -D[56] - D[61];
  D[64]          = -D[57] - D[62];
  D[65]          = Dtmp2 * Dtmp76 + Dtmp49 - Dtmp50 + Dtmp52 + Dtmp67;
  D[66]          = Dtmp8 * (Dtmp42 - Dtmp54 + Dtmp70 + Dtmp74);
  D[67]          = -D[58] - D[65];
  D[68]          = -D[59] - D[66];
  D[69]          = -D[60] - D[67];
  D[70]          = Dtmp31 * Dtmp78;
  D[71]          = Dtmp32 * (Dtmp42 + Dtmp76 - Dtmp79);
  D[72]          = -D[61] - D[70];
  D[73]          = -D[62] - D[71];
  D[74]          = -D[63] - D[72];
  D[75]          = -D[64] - D[73];
  D[76]          = -Dtmp33 * Dtmp60 + 4725.0 * Dtmp51 + Dtmp56 + Dtmp59 * Dtmp80;
  D[77]          = Dtmp78 * Dtmp8;
  D[78]          = -D[65] - D[76];
  D[79]          = -D[66] - D[77];
  D[80]          = -D[67] - D[78];
  D[81]          = -D[68] - D[79];
  D[82]          = -D[69] - D[80];
  D[83]          = -x * (99225.0 * Dtmp63 + Dtmp81 + Dtmp83 - 218295.0 * Dtmp84);
  D[84]          = -Dtmp86 * y;
  D[85]          = -Dtmp86 * z;
  D[86]          = -x * (Dtmp60 * Dtmp7 + Dtmp64 + Dtmp89 + Dtmp92);
  D[87]          = -Dtmp19 * (Dtmp60 + Dtmp87 - Dtmp91);
  D[88]          = -D[83] - D[86];
  D[89]          = -y * (17010.0 * Dtmp63 - 31185.0 * Dtmp84 + Dtmp88 + Dtmp93 + Dtmp95);
  D[90]          = -z * (Dtmp47 + Dtmp66 + Dtmp89 + Dtmp95);
  D[91]          = -D[84] - D[89];
  D[92]          = -D[85] - D[90];
  D[93]          = -x * (-Dtmp33 * Dtmp96 + 17010.0 * Dtmp71 + Dtmp93 + Dtmp99);
  D[94] =
        -Dtmp19 * (8505.0 * Dtmp35 - Dtmp7 * Dtmp96 + Dtmp73 * Dtmp82 - 31185.0 * Dtmp90);
  D[95]  = -D[86] - D[93];
  D[96]  = -D[87] - D[94];
  D[97]  = -D[88] - D[95];
  D[98]  = -y * (Dtmp100 + Dtmp2 * Dtmp60 + Dtmp77 + Dtmp92 + Dtmp98);
  D[99]  = -z * (Dtmp100 + Dtmp47 + Dtmp79 + Dtmp99);
  D[100] = -D[89] - D[98];
  D[101] = -D[90] - D[99];
  D[102] = -D[91] - D[100];
  D[103] = -D[92] - D[101];
  D[104] = -Dtmp103 * Dtmp18;
  D[105] = -Dtmp55 * (-103950.0 * Dtmp58 * Dtmp7 + Dtmp60 + Dtmp97);
  D[106] = -D[93] - D[104];
  D[107] = -D[94] - D[105];
  D[108] = -D[95] - D[106];
  D[109] = -D[96] - D[107];
  D[110] = -D[97] - D[108];
  D[111] = -y * (Dtmp101 - 218295.0 * Dtmp102 + 99225.0 * Dtmp71 + Dtmp81);
  D[112] = -Dtmp103 * z;
  D[113] = -D[98] - D[111];
  D[114] = -D[99] - D[112];
  D[115] = -D[100] - D[113];
  D[116] = -D[101] - D[114];
  D[117] = -D[102] - D[115];
  D[118] = -D[103] - D[116];
#pragma omp atomic
  L[0] += D[0] * M[0] + D[1] * M[1] + D[2] * M[2] + D[3] * M[3] + D[4] * M[4] +
          D[5] * M[5] + D[6] * M[6] + D[7] * M[7] + D[8] * M[8] + D[9] * M[9] +
          D[10] * M[10] + D[11] * M[11] + D[12] * M[12] + D[13] * M[13] + D[14] * M[14] +
          D[15] * M[15] + D[16] * M[16] + D[17] * M[17] + D[18] * M[18] + D[19] * M[19] +
          D[20] * M[20] + D[21] * M[21] + D[22] * M[22] + D[23] * M[23] + D[24] * M[24] +
          D[25] * M[25] + D[26] * M[26] + D[27] * M[27] + D[28] * M[28] + D[29] * M[29] +
          D[30] * M[30] + D[31] * M[31] + D[32] * M[32] + D[33] * M[33] + D[34] * M[34] +
          D[35] * M[35] + D[36] * M[36] + D[37] * M[37] + D[38] * M[38] + D[39] * M[39] +
          D[40] * M[40] + D[41] * M[41] + D[42] * M[42] + D[43] * M[43] + D[44] * M[44] +
          D[45] * M[45] + D[46] * M[46] + D[47] * M[47] + D[48] * M[48] + D[49] * M[49] +
          D[50] * M[50] + D[51] * M[51] + D[52] * M[52] + D[53] * M[53] + D[54] * M[54] +
          D[55] * M[55] + D[56] * M[56] + D[57] * M[57] + D[58] * M[58] + D[59] * M[59] +
          D[60] * M[60] + D[61] * M[61] + D[62] * M[62] + D[63] * M[63] + D[64] * M[64] +
          D[65] * M[65] + D[66] * M[66] + D[67] * M[67] + D[68] * M[68] + D[69] * M[69] +
          D[70] * M[70] + D[71] * M[71] + D[72] * M[72] + D[73] * M[73] + D[74] * M[74] +
          D[75] * M[75] + D[76] * M[76] + D[77] * M[77] + D[78] * M[78] + D[79] * M[79] +
          D[80] * M[80] + D[81] * M[81] + D[82] * M[82] + D[83] * M[83] + D[84] * M[84] +
          D[85] * M[85] + D[86] * M[86] + D[87] * M[87] + D[88] * M[88] + D[89] * M[89] +
          D[90] * M[90] + D[91] * M[91] + D[92] * M[92] + D[93] * M[93] + D[94] * M[94] +
          D[95] * M[95] + D[96] * M[96] + D[97] * M[97] + D[98] * M[98] + D[99] * M[99] +
          D[100] * M[100] + D[101] * M[101] + D[102] * M[102] + D[103] * M[103] +
          D[104] * M[104] + D[105] * M[105] + D[106] * M[106] + D[107] * M[107] +
          D[108] * M[108] + D[109] * M[109] + D[110] * M[110] + D[111] * M[111] +
          D[112] * M[112] + D[113] * M[113] + D[114] * M[114] + D[115] * M[115] +
          D[116] * M[116] + D[117] * M[117] + D[118] * M[118];
#pragma omp atomic
  L[1] += D[3] * M[0] + D[4] * M[1] + D[5] * M[2] + D[9] * M[3] + D[10] * M[4] +
          D[11] * M[5] + D[12] * M[6] + D[13] * M[7] + D[14] * M[8] + D[19] * M[9] +
          D[20] * M[10] + D[21] * M[11] + D[22] * M[12] + D[23] * M[13] + D[24] * M[14] +
          D[25] * M[15] + D[26] * M[16] + D[27] * M[17] + D[28] * M[18] + D[34] * M[19] +
          D[35] * M[20] + D[36] * M[21] + D[37] * M[22] + D[38] * M[23] + D[39] * M[24] +
          D[40] * M[25] + D[41] * M[26] + D[42] * M[27] + D[43] * M[28] + D[44] * M[29] +
          D[45] * M[30] + D[46] * M[31] + D[47] * M[32] + D[48] * M[33] + D[55] * M[34] +
          D[56] * M[35] + D[57] * M[36] + D[58] * M[37] + D[59] * M[38] + D[60] * M[39] +
          D[61] * M[40] + D[62] * M[41] + D[63] * M[42] + D[64] * M[43] + D[65] * M[44] +
          D[66] * M[45] + D[67] * M[46] + D[68] * M[47] + D[69] * M[48] + D[70] * M[49] +
          D[71] * M[50] + D[72] * M[51] + D[73] * M[52] + D[74] * M[53] + D[75] * M[54] +
          D[83] * M[55] + D[84] * M[56] + D[85] * M[57] + D[86] * M[58] + D[87] * M[59] +
          D[88] * M[60] + D[89] * M[61] + D[90] * M[62] + D[91] * M[63] + D[92] * M[64] +
          D[93] * M[65] + D[94] * M[66] + D[95] * M[67] + D[96] * M[68] + D[97] * M[69] +
          D[98] * M[70] + D[99] * M[71] + D[100] * M[72] + D[101] * M[73] +
          D[102] * M[74] + D[103] * M[75] + D[104] * M[76] + D[105] * M[77] +
          D[106] * M[78] + D[107] * M[79] + D[108] * M[80] + D[109] * M[81] +
          D[110] * M[82];
#pragma omp atomic
  L[2] += D[4] * M[0] + D[6] * M[1] + D[7] * M[2] + D[10] * M[3] + D[12] * M[4] +
          D[13] * M[5] + D[15] * M[6] + D[16] * M[7] + D[17] * M[8] + D[20] * M[9] +
          D[22] * M[10] + D[23] * M[11] + D[25] * M[12] + D[26] * M[13] + D[27] * M[14] +
          D[29] * M[15] + D[30] * M[16] + D[31] * M[17] + D[32] * M[18] + D[35] * M[19] +
          D[37] * M[20] + D[38] * M[21] + D[40] * M[22] + D[41] * M[23] + D[42] * M[24] +
          D[44] * M[25] + D[45] * M[26] + D[46] * M[27] + D[47] * M[28] + D[49] * M[29] +
          D[50] * M[30] + D[51] * M[31] + D[52] * M[32] + D[53] * M[33] + D[56] * M[34] +
          D[58] * M[35] + D[59] * M[36] + D[61] * M[37] + D[62] * M[38] + D[63] * M[39] +
          D[65] * M[40] + D[66] * M[41] + D[67] * M[42] + D[68] * M[43] + D[70] * M[44] +
          D[71] * M[45] + D[72] * M[46] + D[73] * M[47] + D[74] * M[48] + D[76] * M[49] +
          D[77] * M[50] + D[78] * M[51] + D[79] * M[52] + D[80] * M[53] + D[81] * M[54] +
          D[84] * M[55] + D[86] * M[56] + D[87] * M[57] + D[89] * M[58] + D[90] * M[59] +
          D[91] * M[60] + D[93] * M[61] + D[94] * M[62] + D[95] * M[63] + D[96] * M[64] +
          D[98] * M[65] + D[99] * M[66] + D[100] * M[67] + D[101] * M[68] +
          D[102] * M[69] + D[104] * M[70] + D[105] * M[71] + D[106] * M[72] +
          D[107] * M[73] + D[108] * M[74] + D[109] * M[75] + D[111] * M[76] +
          D[112] * M[77] + D[113] * M[78] + D[114] * M[79] + D[115] * M[80] +
          D[116] * M[81] + D[117] * M[82];
#pragma omp atomic
  L[3] += D[5] * M[0] + D[7] * M[1] + D[8] * M[2] + D[11] * M[3] + D[13] * M[4] +
          D[14] * M[5] + D[16] * M[6] + D[17] * M[7] + D[18] * M[8] + D[21] * M[9] +
          D[23] * M[10] + D[24] * M[11] + D[26] * M[12] + D[27] * M[13] + D[28] * M[14] +
          D[30] * M[15] + D[31] * M[16] + D[32] * M[17] + D[33] * M[18] + D[36] * M[19] +
          D[38] * M[20] + D[39] * M[21] + D[41] * M[22] + D[42] * M[23] + D[43] * M[24] +
          D[45] * M[25] + D[46] * M[26] + D[47] * M[27] + D[48] * M[28] + D[50] * M[29] +
          D[51] * M[30] + D[52] * M[31] + D[53] * M[32] + D[54] * M[33] + D[57] * M[34] +
          D[59] * M[35] + D[60] * M[36] + D[62] * M[37] + D[63] * M[38] + D[64] * M[39] +
          D[66] * M[40] + D[67] * M[41] + D[68] * M[42] + D[69] * M[43] + D[71] * M[44] +
          D[72] * M[45] + D[73] * M[46] + D[74] * M[47] + D[75] * M[48] + D[77] * M[49] +
          D[78] * M[50] + D[79] * M[51] + D[80] * M[52] + D[81] * M[53] + D[82] * M[54] +
          D[85] * M[55] + D[87] * M[56] + D[88] * M[57] + D[90] * M[58] + D[91] * M[59] +
          D[92] * M[60] + D[94] * M[61] + D[95] * M[62] + D[96] * M[63] + D[97] * M[64] +
          D[99] * M[65] + D[100] * M[66] + D[101] * M[67] + D[102] * M[68] +
          D[103] * M[69] + D[105] * M[70] + D[106] * M[71] + D[107] * M[72] +
          D[108] * M[73] + D[109] * M[74] + D[110] * M[75] + D[112] * M[76] +
          D[113] * M[77] + D[114] * M[78] + D[115] * M[79] + D[116] * M[80] +
          D[117] * M[81] + D[118] * M[82];
#pragma omp atomic
  L[4] += D[9] * M[0] + D[10] * M[1] + D[11] * M[2] + D[19] * M[3] + D[20] * M[4] +
          D[21] * M[5] + D[22] * M[6] + D[23] * M[7] + D[24] * M[8] + D[34] * M[9] +
          D[35] * M[10] + D[36] * M[11] + D[37] * M[12] + D[38] * M[13] + D[39] * M[14] +
          D[40] * M[15] + D[41] * M[16] + D[42] * M[17] + D[43] * M[18] + D[55] * M[19] +
          D[56] * M[20] + D[57] * M[21] + D[58] * M[22] + D[59] * M[23] + D[60] * M[24] +
          D[61] * M[25] + D[62] * M[26] + D[63] * M[27] + D[64] * M[28] + D[65] * M[29] +
          D[66] * M[30] + D[67] * M[31] + D[68] * M[32] + D[69] * M[33] + D[83] * M[34] +
          D[84] * M[35] + D[85] * M[36] + D[86] * M[37] + D[87] * M[38] + D[88] * M[39] +
          D[89] * M[40] + D[90] * M[41] + D[91] * M[42] + D[92] * M[43] + D[93] * M[44] +
          D[94] * M[45] + D[95] * M[46] + D[96] * M[47] + D[97] * M[48] + D[98] * M[49] +
          D[99] * M[50] + D[100] * M[51] + D[101] * M[52] + D[102] * M[53] +
          D[103] * M[54];
#pragma omp atomic
  L[5] += D[10] * M[0] + D[12] * M[1] + D[13] * M[2] + D[20] * M[3] + D[22] * M[4] +
          D[23] * M[5] + D[25] * M[6] + D[26] * M[7] + D[27] * M[8] + D[35] * M[9] +
          D[37] * M[10] + D[38] * M[11] + D[40] * M[12] + D[41] * M[13] + D[42] * M[14] +
          D[44] * M[15] + D[45] * M[16] + D[46] * M[17] + D[47] * M[18] + D[56] * M[19] +
          D[58] * M[20] + D[59] * M[21] + D[61] * M[22] + D[62] * M[23] + D[63] * M[24] +
          D[65] * M[25] + D[66] * M[26] + D[67] * M[27] + D[68] * M[28] + D[70] * M[29] +
          D[71] * M[30] + D[72] * M[31] + D[73] * M[32] + D[74] * M[33] + D[84] * M[34] +
          D[86] * M[35] + D[87] * M[36] + D[89] * M[37] + D[90] * M[38] + D[91] * M[39] +
          D[93] * M[40] + D[94] * M[41] + D[95] * M[42] + D[96] * M[43] + D[98] * M[44] +
          D[99] * M[45] + D[100] * M[46] + D[101] * M[47] + D[102] * M[48] +
          D[104] * M[49] + D[105] * M[50] + D[106] * M[51] + D[107] * M[52] +
          D[108] * M[53] + D[109] * M[54];
#pragma omp atomic
  L[6] += D[11] * M[0] + D[13] * M[1] + D[14] * M[2] + D[21] * M[3] + D[23] * M[4] +
          D[24] * M[5] + D[26] * M[6] + D[27] * M[7] + D[28] * M[8] + D[36] * M[9] +
          D[38] * M[10] + D[39] * M[11] + D[41] * M[12] + D[42] * M[13] + D[43] * M[14] +
          D[45] * M[15] + D[46] * M[16] + D[47] * M[17] + D[48] * M[18] + D[57] * M[19] +
          D[59] * M[20] + D[60] * M[21] + D[62] * M[22] + D[63] * M[23] + D[64] * M[24] +
          D[66] * M[25] + D[67] * M[26] + D[68] * M[27] + D[69] * M[28] + D[71] * M[29] +
          D[72] * M[30] + D[73] * M[31] + D[74] * M[32] + D[75] * M[33] + D[85] * M[34] +
          D[87] * M[35] + D[88] * M[36] + D[90] * M[37] + D[91] * M[38] + D[92] * M[39] +
          D[94] * M[40] + D[95] * M[41] + D[96] * M[42] + D[97] * M[43] + D[99] * M[44] +
          D[100] * M[45] + D[101] * M[46] + D[102] * M[47] + D[103] * M[48] +
          D[105] * M[49] + D[106] * M[50] + D[107] * M[51] + D[108] * M[52] +
          D[109] * M[53] + D[110] * M[54];
#pragma omp atomic
  L[7] += D[12] * M[0] + D[15] * M[1] + D[16] * M[2] + D[22] * M[3] + D[25] * M[4] +
          D[26] * M[5] + D[29] * M[6] + D[30] * M[7] + D[31] * M[8] + D[37] * M[9] +
          D[40] * M[10] + D[41] * M[11] + D[44] * M[12] + D[45] * M[13] + D[46] * M[14] +
          D[49] * M[15] + D[50] * M[16] + D[51] * M[17] + D[52] * M[18] + D[58] * M[19] +
          D[61] * M[20] + D[62] * M[21] + D[65] * M[22] + D[66] * M[23] + D[67] * M[24] +
          D[70] * M[25] + D[71] * M[26] + D[72] * M[27] + D[73] * M[28] + D[76] * M[29] +
          D[77] * M[30] + D[78] * M[31] + D[79] * M[32] + D[80] * M[33] + D[86] * M[34] +
          D[89] * M[35] + D[90] * M[36] + D[93] * M[37] + D[94] * M[38] + D[95] * M[39] +
          D[98] * M[40] + D[99] * M[41] + D[100] * M[42] + D[101] * M[43] +
          D[104] * M[44] + D[105] * M[45] + D[106] * M[46] + D[107] * M[47] +
          D[108] * M[48] + D[111] * M[49] + D[112] * M[50] + D[113] * M[51] +
          D[114] * M[52] + D[115] * M[53] + D[116] * M[54];
#pragma omp atomic
  L[8] += D[13] * M[0] + D[16] * M[1] + D[17] * M[2] + D[23] * M[3] + D[26] * M[4] +
          D[27] * M[5] + D[30] * M[6] + D[31] * M[7] + D[32] * M[8] + D[38] * M[9] +
          D[41] * M[10] + D[42] * M[11] + D[45] * M[12] + D[46] * M[13] + D[47] * M[14] +
          D[50] * M[15] + D[51] * M[16] + D[52] * M[17] + D[53] * M[18] + D[59] * M[19] +
          D[62] * M[20] + D[63] * M[21] + D[66] * M[22] + D[67] * M[23] + D[68] * M[24] +
          D[71] * M[25] + D[72] * M[26] + D[73] * M[27] + D[74] * M[28] + D[77] * M[29] +
          D[78] * M[30] + D[79] * M[31] + D[80] * M[32] + D[81] * M[33] + D[87] * M[34] +
          D[90] * M[35] + D[91] * M[36] + D[94] * M[37] + D[95] * M[38] + D[96] * M[39] +
          D[99] * M[40] + D[100] * M[41] + D[101] * M[42] + D[102] * M[43] +
          D[105] * M[44] + D[106] * M[45] + D[107] * M[46] + D[108] * M[47] +
          D[109] * M[48] + D[112] * M[49] + D[113] * M[50] + D[114] * M[51] +
          D[115] * M[52] + D[116] * M[53] + D[117] * M[54];
#pragma omp atomic
  L[9] += D[14] * M[0] + D[17] * M[1] + D[18] * M[2] + D[24] * M[3] + D[27] * M[4] +
          D[28] * M[5] + D[31] * M[6] + D[32] * M[7] + D[33] * M[8] + D[39] * M[9] +
          D[42] * M[10] + D[43] * M[11] + D[46] * M[12] + D[47] * M[13] + D[48] * M[14] +
          D[51] * M[15] + D[52] * M[16] + D[53] * M[17] + D[54] * M[18] + D[60] * M[19] +
          D[63] * M[20] + D[64] * M[21] + D[67] * M[22] + D[68] * M[23] + D[69] * M[24] +
          D[72] * M[25] + D[73] * M[26] + D[74] * M[27] + D[75] * M[28] + D[78] * M[29] +
          D[79] * M[30] + D[80] * M[31] + D[81] * M[32] + D[82] * M[33] + D[88] * M[34] +
          D[91] * M[35] + D[92] * M[36] + D[95] * M[37] + D[96] * M[38] + D[97] * M[39] +
          D[100] * M[40] + D[101] * M[41] + D[102] * M[42] + D[103] * M[43] +
          D[106] * M[44] + D[107] * M[45] + D[108] * M[46] + D[109] * M[47] +
          D[110] * M[48] + D[113] * M[49] + D[114] * M[50] + D[115] * M[51] +
          D[116] * M[52] + D[117] * M[53] + D[118] * M[54];
#pragma omp atomic
  L[10] += D[19] * M[0] + D[20] * M[1] + D[21] * M[2] + D[34] * M[3] + D[35] * M[4] +
           D[36] * M[5] + D[37] * M[6] + D[38] * M[7] + D[39] * M[8] + D[55] * M[9] +
           D[56] * M[10] + D[57] * M[11] + D[58] * M[12] + D[59] * M[13] + D[60] * M[14] +
           D[61] * M[15] + D[62] * M[16] + D[63] * M[17] + D[64] * M[18] + D[83] * M[19] +
           D[84] * M[20] + D[85] * M[21] + D[86] * M[22] + D[87] * M[23] + D[88] * M[24] +
           D[89] * M[25] + D[90] * M[26] + D[91] * M[27] + D[92] * M[28] + D[93] * M[29] +
           D[94] * M[30] + D[95] * M[31] + D[96] * M[32] + D[97] * M[33];
#pragma omp atomic
  L[11] += D[20] * M[0] + D[22] * M[1] + D[23] * M[2] + D[35] * M[3] + D[37] * M[4] +
           D[38] * M[5] + D[40] * M[6] + D[41] * M[7] + D[42] * M[8] + D[56] * M[9] +
           D[58] * M[10] + D[59] * M[11] + D[61] * M[12] + D[62] * M[13] + D[63] * M[14] +
           D[65] * M[15] + D[66] * M[16] + D[67] * M[17] + D[68] * M[18] + D[84] * M[19] +
           D[86] * M[20] + D[87] * M[21] + D[89] * M[22] + D[90] * M[23] + D[91] * M[24] +
           D[93] * M[25] + D[94] * M[26] + D[95] * M[27] + D[96] * M[28] + D[98] * M[29] +
           D[99] * M[30] + D[100] * M[31] + D[101] * M[32] + D[102] * M[33];
#pragma omp atomic
  L[12] += D[21] * M[0] + D[23] * M[1] + D[24] * M[2] + D[36] * M[3] + D[38] * M[4] +
           D[39] * M[5] + D[41] * M[6] + D[42] * M[7] + D[43] * M[8] + D[57] * M[9] +
           D[59] * M[10] + D[60] * M[11] + D[62] * M[12] + D[63] * M[13] + D[64] * M[14] +
           D[66] * M[15] + D[67] * M[16] + D[68] * M[17] + D[69] * M[18] + D[85] * M[19] +
           D[87] * M[20] + D[88] * M[21] + D[90] * M[22] + D[91] * M[23] + D[92] * M[24] +
           D[94] * M[25] + D[95] * M[26] + D[96] * M[27] + D[97] * M[28] + D[99] * M[29] +
           D[100] * M[30] + D[101] * M[31] + D[102] * M[32] + D[103] * M[33];
#pragma omp atomic
  L[13] += D[22] * M[0] + D[25] * M[1] + D[26] * M[2] + D[37] * M[3] + D[40] * M[4] +
           D[41] * M[5] + D[44] * M[6] + D[45] * M[7] + D[46] * M[8] + D[58] * M[9] +
           D[61] * M[10] + D[62] * M[11] + D[65] * M[12] + D[66] * M[13] + D[67] * M[14] +
           D[70] * M[15] + D[71] * M[16] + D[72] * M[17] + D[73] * M[18] + D[86] * M[19] +
           D[89] * M[20] + D[90] * M[21] + D[93] * M[22] + D[94] * M[23] + D[95] * M[24] +
           D[98] * M[25] + D[99] * M[26] + D[100] * M[27] + D[101] * M[28] +
           D[104] * M[29] + D[105] * M[30] + D[106] * M[31] + D[107] * M[32] +
           D[108] * M[33];
#pragma omp atomic
  L[14] += D[23] * M[0] + D[26] * M[1] + D[27] * M[2] + D[38] * M[3] + D[41] * M[4] +
           D[42] * M[5] + D[45] * M[6] + D[46] * M[7] + D[47] * M[8] + D[59] * M[9] +
           D[62] * M[10] + D[63] * M[11] + D[66] * M[12] + D[67] * M[13] + D[68] * M[14] +
           D[71] * M[15] + D[72] * M[16] + D[73] * M[17] + D[74] * M[18] + D[87] * M[19] +
           D[90] * M[20] + D[91] * M[21] + D[94] * M[22] + D[95] * M[23] + D[96] * M[24] +
           D[99] * M[25] + D[100] * M[26] + D[101] * M[27] + D[102] * M[28] +
           D[105] * M[29] + D[106] * M[30] + D[107] * M[31] + D[108] * M[32] +
           D[109] * M[33];
#pragma omp atomic
  L[15] += D[24] * M[0] + D[27] * M[1] + D[28] * M[2] + D[39] * M[3] + D[42] * M[4] +
           D[43] * M[5] + D[46] * M[6] + D[47] * M[7] + D[48] * M[8] + D[60] * M[9] +
           D[63] * M[10] + D[64] * M[11] + D[67] * M[12] + D[68] * M[13] + D[69] * M[14] +
           D[72] * M[15] + D[73] * M[16] + D[74] * M[17] + D[75] * M[18] + D[88] * M[19] +
           D[91] * M[20] + D[92] * M[21] + D[95] * M[22] + D[96] * M[23] + D[97] * M[24] +
           D[100] * M[25] + D[101] * M[26] + D[102] * M[27] + D[103] * M[28] +
           D[106] * M[29] + D[107] * M[30] + D[108] * M[31] + D[109] * M[32] +
           D[110] * M[33];
#pragma omp atomic
  L[16] += D[25] * M[0] + D[29] * M[1] + D[30] * M[2] + D[40] * M[3] + D[44] * M[4] +
           D[45] * M[5] + D[49] * M[6] + D[50] * M[7] + D[51] * M[8] + D[61] * M[9] +
           D[65] * M[10] + D[66] * M[11] + D[70] * M[12] + D[71] * M[13] + D[72] * M[14] +
           D[76] * M[15] + D[77] * M[16] + D[78] * M[17] + D[79] * M[18] + D[89] * M[19] +
           D[93] * M[20] + D[94] * M[21] + D[98] * M[22] + D[99] * M[23] +
           D[100] * M[24] + D[104] * M[25] + D[105] * M[26] + D[106] * M[27] +
           D[107] * M[28] + D[111] * M[29] + D[112] * M[30] + D[113] * M[31] +
           D[114] * M[32] + D[115] * M[33];
#pragma omp atomic
  L[17] += D[26] * M[0] + D[30] * M[1] + D[31] * M[2] + D[41] * M[3] + D[45] * M[4] +
           D[46] * M[5] + D[50] * M[6] + D[51] * M[7] + D[52] * M[8] + D[62] * M[9] +
           D[66] * M[10] + D[67] * M[11] + D[71] * M[12] + D[72] * M[13] + D[73] * M[14] +
           D[77] * M[15] + D[78] * M[16] + D[79] * M[17] + D[80] * M[18] + D[90] * M[19] +
           D[94] * M[20] + D[95] * M[21] + D[99] * M[22] + D[100] * M[23] +
           D[101] * M[24] + D[105] * M[25] + D[106] * M[26] + D[107] * M[27] +
           D[108] * M[28] + D[112] * M[29] + D[113] * M[30] + D[114] * M[31] +
           D[115] * M[32] + D[116] * M[33];
#pragma omp atomic
  L[18] += D[27] * M[0] + D[31] * M[1] + D[32] * M[2] + D[42] * M[3] + D[46] * M[4] +
           D[47] * M[5] + D[51] * M[6] + D[52] * M[7] + D[53] * M[8] + D[63] * M[9] +
           D[67] * M[10] + D[68] * M[11] + D[72] * M[12] + D[73] * M[13] + D[74] * M[14] +
           D[78] * M[15] + D[79] * M[16] + D[80] * M[17] + D[81] * M[18] + D[91] * M[19] +
           D[95] * M[20] + D[96] * M[21] + D[100] * M[22] + D[101] * M[23] +
           D[102] * M[24] + D[106] * M[25] + D[107] * M[26] + D[108] * M[27] +
           D[109] * M[28] + D[113] * M[29] + D[114] * M[30] + D[115] * M[31] +
           D[116] * M[32] + D[117] * M[33];
#pragma omp atomic
  L[19] += D[28] * M[0] + D[32] * M[1] + D[33] * M[2] + D[43] * M[3] + D[47] * M[4] +
           D[48] * M[5] + D[52] * M[6] + D[53] * M[7] + D[54] * M[8] + D[64] * M[9] +
           D[68] * M[10] + D[69] * M[11] + D[73] * M[12] + D[74] * M[13] + D[75] * M[14] +
           D[79] * M[15] + D[80] * M[16] + D[81] * M[17] + D[82] * M[18] + D[92] * M[19] +
           D[96] * M[20] + D[97] * M[21] + D[101] * M[22] + D[102] * M[23] +
           D[103] * M[24] + D[107] * M[25] + D[108] * M[26] + D[109] * M[27] +
           D[110] * M[28] + D[114] * M[29] + D[115] * M[30] + D[116] * M[31] +
           D[117] * M[32] + D[118] * M[33];
#pragma omp atomic
  L[20] += D[34] * M[0] + D[35] * M[1] + D[36] * M[2] + D[55] * M[3] + D[56] * M[4] +
           D[57] * M[5] + D[58] * M[6] + D[59] * M[7] + D[60] * M[8] + D[83] * M[9] +
           D[84] * M[10] + D[85] * M[11] + D[86] * M[12] + D[87] * M[13] + D[88] * M[14] +
           D[89] * M[15] + D[90] * M[16] + D[91] * M[17] + D[92] * M[18];
#pragma omp atomic
  L[21] += D[35] * M[0] + D[37] * M[1] + D[38] * M[2] + D[56] * M[3] + D[58] * M[4] +
           D[59] * M[5] + D[61] * M[6] + D[62] * M[7] + D[63] * M[8] + D[84] * M[9] +
           D[86] * M[10] + D[87] * M[11] + D[89] * M[12] + D[90] * M[13] + D[91] * M[14] +
           D[93] * M[15] + D[94] * M[16] + D[95] * M[17] + D[96] * M[18];
#pragma omp atomic
  L[22] += D[36] * M[0] + D[38] * M[1] + D[39] * M[2] + D[57] * M[3] + D[59] * M[4] +
           D[60] * M[5] + D[62] * M[6] + D[63] * M[7] + D[64] * M[8] + D[85] * M[9] +
           D[87] * M[10] + D[88] * M[11] + D[90] * M[12] + D[91] * M[13] + D[92] * M[14] +
           D[94] * M[15] + D[95] * M[16] + D[96] * M[17] + D[97] * M[18];
#pragma omp atomic
  L[23] += D[37] * M[0] + D[40] * M[1] + D[41] * M[2] + D[58] * M[3] + D[61] * M[4] +
           D[62] * M[5] + D[65] * M[6] + D[66] * M[7] + D[67] * M[8] + D[86] * M[9] +
           D[89] * M[10] + D[90] * M[11] + D[93] * M[12] + D[94] * M[13] + D[95] * M[14] +
           D[98] * M[15] + D[99] * M[16] + D[100] * M[17] + D[101] * M[18];
#pragma omp atomic
  L[24] += D[38] * M[0] + D[41] * M[1] + D[42] * M[2] + D[59] * M[3] + D[62] * M[4] +
           D[63] * M[5] + D[66] * M[6] + D[67] * M[7] + D[68] * M[8] + D[87] * M[9] +
           D[90] * M[10] + D[91] * M[11] + D[94] * M[12] + D[95] * M[13] + D[96] * M[14] +
           D[99] * M[15] + D[100] * M[16] + D[101] * M[17] + D[102] * M[18];
#pragma omp atomic
  L[25] += D[39] * M[0] + D[42] * M[1] + D[43] * M[2] + D[60] * M[3] + D[63] * M[4] +
           D[64] * M[5] + D[67] * M[6] + D[68] * M[7] + D[69] * M[8] + D[88] * M[9] +
           D[91] * M[10] + D[92] * M[11] + D[95] * M[12] + D[96] * M[13] + D[97] * M[14] +
           D[100] * M[15] + D[101] * M[16] + D[102] * M[17] + D[103] * M[18];
#pragma omp atomic
  L[26] += D[40] * M[0] + D[44] * M[1] + D[45] * M[2] + D[61] * M[3] + D[65] * M[4] +
           D[66] * M[5] + D[70] * M[6] + D[71] * M[7] + D[72] * M[8] + D[89] * M[9] +
           D[93] * M[10] + D[94] * M[11] + D[98] * M[12] + D[99] * M[13] +
           D[100] * M[14] + D[104] * M[15] + D[105] * M[16] + D[106] * M[17] +
           D[107] * M[18];
#pragma omp atomic
  L[27] += D[41] * M[0] + D[45] * M[1] + D[46] * M[2] + D[62] * M[3] + D[66] * M[4] +
           D[67] * M[5] + D[71] * M[6] + D[72] * M[7] + D[73] * M[8] + D[90] * M[9] +
           D[94] * M[10] + D[95] * M[11] + D[99] * M[12] + D[100] * M[13] +
           D[101] * M[14] + D[105] * M[15] + D[106] * M[16] + D[107] * M[17] +
           D[108] * M[18];
#pragma omp atomic
  L[28] += D[42] * M[0] + D[46] * M[1] + D[47] * M[2] + D[63] * M[3] + D[67] * M[4] +
           D[68] * M[5] + D[72] * M[6] + D[73] * M[7] + D[74] * M[8] + D[91] * M[9] +
           D[95] * M[10] + D[96] * M[11] + D[100] * M[12] + D[101] * M[13] +
           D[102] * M[14] + D[106] * M[15] + D[107] * M[16] + D[108] * M[17] +
           D[109] * M[18];
#pragma omp atomic
  L[29] += D[43] * M[0] + D[47] * M[1] + D[48] * M[2] + D[64] * M[3] + D[68] * M[4] +
           D[69] * M[5] + D[73] * M[6] + D[74] * M[7] + D[75] * M[8] + D[92] * M[9] +
           D[96] * M[10] + D[97] * M[11] + D[101] * M[12] + D[102] * M[13] +
           D[103] * M[14] + D[107] * M[15] + D[108] * M[16] + D[109] * M[17] +
           D[110] * M[18];
#pragma omp atomic
  L[30] += D[44] * M[0] + D[49] * M[1] + D[50] * M[2] + D[65] * M[3] + D[70] * M[4] +
           D[71] * M[5] + D[76] * M[6] + D[77] * M[7] + D[78] * M[8] + D[93] * M[9] +
           D[98] * M[10] + D[99] * M[11] + D[104] * M[12] + D[105] * M[13] +
           D[106] * M[14] + D[111] * M[15] + D[112] * M[16] + D[113] * M[17] +
           D[114] * M[18];
#pragma omp atomic
  L[31] += D[45] * M[0] + D[50] * M[1] + D[51] * M[2] + D[66] * M[3] + D[71] * M[4] +
           D[72] * M[5] + D[77] * M[6] + D[78] * M[7] + D[79] * M[8] + D[94] * M[9] +
           D[99] * M[10] + D[100] * M[11] + D[105] * M[12] + D[106] * M[13] +
           D[107] * M[14] + D[112] * M[15] + D[113] * M[16] + D[114] * M[17] +
           D[115] * M[18];
#pragma omp atomic
  L[32] += D[46] * M[0] + D[51] * M[1] + D[52] * M[2] + D[67] * M[3] + D[72] * M[4] +
           D[73] * M[5] + D[78] * M[6] + D[79] * M[7] + D[80] * M[8] + D[95] * M[9] +
           D[100] * M[10] + D[101] * M[11] + D[106] * M[12] + D[107] * M[13] +
           D[108] * M[14] + D[113] * M[15] + D[114] * M[16] + D[115] * M[17] +
           D[116] * M[18];
#pragma omp atomic
  L[33] += D[47] * M[0] + D[52] * M[1] + D[53] * M[2] + D[68] * M[3] + D[73] * M[4] +
           D[74] * M[5] + D[79] * M[6] + D[80] * M[7] + D[81] * M[8] + D[96] * M[9] +
           D[101] * M[10] + D[102] * M[11] + D[107] * M[12] + D[108] * M[13] +
           D[109] * M[14] + D[114] * M[15] + D[115] * M[16] + D[116] * M[17] +
           D[117] * M[18];
#pragma omp atomic
  L[34] += D[48] * M[0] + D[53] * M[1] + D[54] * M[2] + D[69] * M[3] + D[74] * M[4] +
           D[75] * M[5] + D[80] * M[6] + D[81] * M[7] + D[82] * M[8] + D[97] * M[9] +
           D[102] * M[10] + D[103] * M[11] + D[108] * M[12] + D[109] * M[13] +
           D[110] * M[14] + D[115] * M[15] + D[116] * M[16] + D[117] * M[17] +
           D[118] * M[18];
#pragma omp atomic
  L[35] += D[55] * M[0] + D[56] * M[1] + D[57] * M[2] + D[83] * M[3] + D[84] * M[4] +
           D[85] * M[5] + D[86] * M[6] + D[87] * M[7] + D[88] * M[8];
#pragma omp atomic
  L[36] += D[56] * M[0] + D[58] * M[1] + D[59] * M[2] + D[84] * M[3] + D[86] * M[4] +
           D[87] * M[5] + D[89] * M[6] + D[90] * M[7] + D[91] * M[8];
#pragma omp atomic
  L[37] += D[57] * M[0] + D[59] * M[1] + D[60] * M[2] + D[85] * M[3] + D[87] * M[4] +
           D[88] * M[5] + D[90] * M[6] + D[91] * M[7] + D[92] * M[8];
#pragma omp atomic
  L[38] += D[58] * M[0] + D[61] * M[1] + D[62] * M[2] + D[86] * M[3] + D[89] * M[4] +
           D[90] * M[5] + D[93] * M[6] + D[94] * M[7] + D[95] * M[8];
#pragma omp atomic
  L[39] += D[59] * M[0] + D[62] * M[1] + D[63] * M[2] + D[87] * M[3] + D[90] * M[4] +
           D[91] * M[5] + D[94] * M[6] + D[95] * M[7] + D[96] * M[8];
#pragma omp atomic
  L[40] += D[60] * M[0] + D[63] * M[1] + D[64] * M[2] + D[88] * M[3] + D[91] * M[4] +
           D[92] * M[5] + D[95] * M[6] + D[96] * M[7] + D[97] * M[8];
#pragma omp atomic
  L[41] += D[61] * M[0] + D[65] * M[1] + D[66] * M[2] + D[89] * M[3] + D[93] * M[4] +
           D[94] * M[5] + D[98] * M[6] + D[99] * M[7] + D[100] * M[8];
#pragma omp atomic
  L[42] += D[62] * M[0] + D[66] * M[1] + D[67] * M[2] + D[90] * M[3] + D[94] * M[4] +
           D[95] * M[5] + D[99] * M[6] + D[100] * M[7] + D[101] * M[8];
#pragma omp atomic
  L[43] += D[63] * M[0] + D[67] * M[1] + D[68] * M[2] + D[91] * M[3] + D[95] * M[4] +
           D[96] * M[5] + D[100] * M[6] + D[101] * M[7] + D[102] * M[8];
#pragma omp atomic
  L[44] += D[64] * M[0] + D[68] * M[1] + D[69] * M[2] + D[92] * M[3] + D[96] * M[4] +
           D[97] * M[5] + D[101] * M[6] + D[102] * M[7] + D[103] * M[8];
#pragma omp atomic
  L[45] += D[65] * M[0] + D[70] * M[1] + D[71] * M[2] + D[93] * M[3] + D[98] * M[4] +
           D[99] * M[5] + D[104] * M[6] + D[105] * M[7] + D[106] * M[8];
#pragma omp atomic
  L[46] += D[66] * M[0] + D[71] * M[1] + D[72] * M[2] + D[94] * M[3] + D[99] * M[4] +
           D[100] * M[5] + D[105] * M[6] + D[106] * M[7] + D[107] * M[8];
#pragma omp atomic
  L[47] += D[67] * M[0] + D[72] * M[1] + D[73] * M[2] + D[95] * M[3] + D[100] * M[4] +
           D[101] * M[5] + D[106] * M[6] + D[107] * M[7] + D[108] * M[8];
#pragma omp atomic
  L[48] += D[68] * M[0] + D[73] * M[1] + D[74] * M[2] + D[96] * M[3] + D[101] * M[4] +
           D[102] * M[5] + D[107] * M[6] + D[108] * M[7] + D[109] * M[8];
#pragma omp atomic
  L[49] += D[69] * M[0] + D[74] * M[1] + D[75] * M[2] + D[97] * M[3] + D[102] * M[4] +
           D[103] * M[5] + D[108] * M[6] + D[109] * M[7] + D[110] * M[8];
#pragma omp atomic
  L[50] += D[70] * M[0] + D[76] * M[1] + D[77] * M[2] + D[98] * M[3] + D[104] * M[4] +
           D[105] * M[5] + D[111] * M[6] + D[112] * M[7] + D[113] * M[8];
#pragma omp atomic
  L[51] += D[71] * M[0] + D[77] * M[1] + D[78] * M[2] + D[99] * M[3] + D[105] * M[4] +
           D[106] * M[5] + D[112] * M[6] + D[113] * M[7] + D[114] * M[8];
#pragma omp atomic
  L[52] += D[72] * M[0] + D[78] * M[1] + D[79] * M[2] + D[100] * M[3] + D[106] * M[4] +
           D[107] * M[5] + D[113] * M[6] + D[114] * M[7] + D[115] * M[8];
#pragma omp atomic
  L[53] += D[73] * M[0] + D[79] * M[1] + D[80] * M[2] + D[101] * M[3] + D[107] * M[4] +
           D[108] * M[5] + D[114] * M[6] + D[115] * M[7] + D[116] * M[8];
#pragma omp atomic
  L[54] += D[74] * M[0] + D[80] * M[1] + D[81] * M[2] + D[102] * M[3] + D[108] * M[4] +
           D[109] * M[5] + D[115] * M[6] + D[116] * M[7] + D[117] * M[8];
#pragma omp atomic
  L[55] += D[75] * M[0] + D[81] * M[1] + D[82] * M[2] + D[103] * M[3] + D[109] * M[4] +
           D[110] * M[5] + D[116] * M[6] + D[117] * M[7] + D[118] * M[8];
#pragma omp atomic
  L[56] += D[83] * M[0] + D[84] * M[1] + D[85] * M[2];
#pragma omp atomic
  L[57] += D[84] * M[0] + D[86] * M[1] + D[87] * M[2];
#pragma omp atomic
  L[58] += D[85] * M[0] + D[87] * M[1] + D[88] * M[2];
#pragma omp atomic
  L[59] += D[86] * M[0] + D[89] * M[1] + D[90] * M[2];
#pragma omp atomic
  L[60] += D[87] * M[0] + D[90] * M[1] + D[91] * M[2];
#pragma omp atomic
  L[61] += D[88] * M[0] + D[91] * M[1] + D[92] * M[2];
#pragma omp atomic
  L[62] += D[89] * M[0] + D[93] * M[1] + D[94] * M[2];
#pragma omp atomic
  L[63] += D[90] * M[0] + D[94] * M[1] + D[95] * M[2];
#pragma omp atomic
  L[64] += D[91] * M[0] + D[95] * M[1] + D[96] * M[2];
#pragma omp atomic
  L[65] += D[92] * M[0] + D[96] * M[1] + D[97] * M[2];
#pragma omp atomic
  L[66] += D[93] * M[0] + D[98] * M[1] + D[99] * M[2];
#pragma omp atomic
  L[67] += D[94] * M[0] + D[99] * M[1] + D[100] * M[2];
#pragma omp atomic
  L[68] += D[95] * M[0] + D[100] * M[1] + D[101] * M[2];
#pragma omp atomic
  L[69] += D[96] * M[0] + D[101] * M[1] + D[102] * M[2];
#pragma omp atomic
  L[70] += D[97] * M[0] + D[102] * M[1] + D[103] * M[2];
#pragma omp atomic
  L[71] += D[98] * M[0] + D[104] * M[1] + D[105] * M[2];
#pragma omp atomic
  L[72] += D[99] * M[0] + D[105] * M[1] + D[106] * M[2];
#pragma omp atomic
  L[73] += D[100] * M[0] + D[106] * M[1] + D[107] * M[2];
#pragma omp atomic
  L[74] += D[101] * M[0] + D[107] * M[1] + D[108] * M[2];
#pragma omp atomic
  L[75] += D[102] * M[0] + D[108] * M[1] + D[109] * M[2];
#pragma omp atomic
  L[76] += D[103] * M[0] + D[109] * M[1] + D[110] * M[2];
#pragma omp atomic
  L[77] += D[104] * M[0] + D[111] * M[1] + D[112] * M[2];
#pragma omp atomic
  L[78] += D[105] * M[0] + D[112] * M[1] + D[113] * M[2];
#pragma omp atomic
  L[79] += D[106] * M[0] + D[113] * M[1] + D[114] * M[2];
#pragma omp atomic
  L[80] += D[107] * M[0] + D[114] * M[1] + D[115] * M[2];
#pragma omp atomic
  L[81] += D[108] * M[0] + D[115] * M[1] + D[116] * M[2];
#pragma omp atomic
  L[82] += D[109] * M[0] + D[116] * M[1] + D[117] * M[2];
#pragma omp atomic
  L[83] += D[110] * M[0] + D[117] * M[1] + D[118] * M[2];
}

void field_m1_L2L_7(double x, double y, double z, double* L, double* Ls) {
  double Lstmp0   = y * L[5];
  double Lstmp1   = z * L[6];
  double Lstmp2   = z * L[8];
  double Lstmp3   = z * L[14];
  double Lstmp4   = Lstmp3 * y;
  double Lstmp5   = (x * x);
  double Lstmp6   = (1.0 / 2.0) * Lstmp5;
  double Lstmp7   = (x * x * x);
  double Lstmp8   = (1.0 / 6.0) * Lstmp7;
  double Lstmp9   = (x * x * x * x);
  double Lstmp10  = (1.0 / 24.0) * Lstmp9;
  double Lstmp11  = (1.0 / 120.0) * (x * x * x * x * x);
  double Lstmp12  = (y * y);
  double Lstmp13  = (1.0 / 2.0) * Lstmp12;
  double Lstmp14  = (y * y * y);
  double Lstmp15  = (1.0 / 6.0) * Lstmp14;
  double Lstmp16  = (y * y * y * y);
  double Lstmp17  = (1.0 / 24.0) * Lstmp16;
  double Lstmp18  = (1.0 / 120.0) * (y * y * y * y * y);
  double Lstmp19  = (z * z);
  double Lstmp20  = (1.0 / 2.0) * Lstmp19;
  double Lstmp21  = (z * z * z);
  double Lstmp22  = (1.0 / 6.0) * Lstmp21;
  double Lstmp23  = (z * z * z * z);
  double Lstmp24  = (1.0 / 24.0) * Lstmp23;
  double Lstmp25  = (1.0 / 120.0) * (z * z * z * z * z);
  double Lstmp26  = x * L[13];
  double Lstmp27  = x * L[26];
  double Lstmp28  = x * L[45];
  double Lstmp29  = x * L[71];
  double Lstmp30  = x * L[15];
  double Lstmp31  = x * L[29];
  double Lstmp32  = x * L[49];
  double Lstmp33  = x * L[76];
  double Lstmp34  = y * L[11];
  double Lstmp35  = z * L[12];
  double Lstmp36  = y * L[21];
  double Lstmp37  = z * L[22];
  double Lstmp38  = y * L[36];
  double Lstmp39  = z * L[37];
  double Lstmp40  = y * L[57];
  double Lstmp41  = z * L[58];
  double Lstmp42  = y * L[18];
  double Lstmp43  = y * L[33];
  double Lstmp44  = y * L[54];
  double Lstmp45  = y * L[82];
  double Lstmp46  = z * L[17];
  double Lstmp47  = z * L[31];
  double Lstmp48  = z * L[51];
  double Lstmp49  = z * L[78];
  double Lstmp50  = y * L[28];
  double Lstmp51  = Lstmp50 * x;
  double Lstmp52  = y * L[48];
  double Lstmp53  = Lstmp52 * x;
  double Lstmp54  = y * L[75];
  double Lstmp55  = Lstmp54 * x;
  double Lstmp56  = z * L[27];
  double Lstmp57  = Lstmp56 * x;
  double Lstmp58  = z * L[46];
  double Lstmp59  = Lstmp58 * x;
  double Lstmp60  = z * L[72];
  double Lstmp61  = Lstmp60 * x;
  double Lstmp62  = z * L[24];
  double Lstmp63  = Lstmp62 * y;
  double Lstmp64  = z * L[39];
  double Lstmp65  = Lstmp64 * y;
  double Lstmp66  = z * L[60];
  double Lstmp67  = Lstmp66 * y;
  double Lstmp68  = (1.0 / 4.0) * Lstmp5;
  double Lstmp69  = Lstmp12 * Lstmp68;
  double Lstmp70  = (1.0 / 12.0) * Lstmp5;
  double Lstmp71  = Lstmp14 * Lstmp70;
  double Lstmp72  = (1.0 / 48.0) * Lstmp5;
  double Lstmp73  = Lstmp19 * Lstmp68;
  double Lstmp74  = Lstmp21 * Lstmp70;
  double Lstmp75  = (1.0 / 12.0) * Lstmp7;
  double Lstmp76  = Lstmp12 * Lstmp75;
  double Lstmp77  = (1.0 / 36.0) * Lstmp7;
  double Lstmp78  = Lstmp19 * Lstmp75;
  double Lstmp79  = (1.0 / 48.0) * Lstmp9;
  double Lstmp80  = Lstmp12 * Lstmp19;
  double Lstmp81  = (1.0 / 4.0) * Lstmp80;
  double Lstmp82  = (1.0 / 12.0) * Lstmp12 * Lstmp21;
  double Lstmp83  = (1.0 / 12.0) * Lstmp14 * Lstmp19;
  double Lstmp84  = x * L[47];
  double Lstmp85  = x * L[74];
  double Lstmp86  = x * L[73];
  double Lstmp87  = y * L[43];
  double Lstmp88  = y * L[69];
  double Lstmp89  = z * L[42];
  double Lstmp90  = z * L[67];
  double Lstmp91  = y * L[64];
  double Lstmp92  = z * L[63];
  double Lstmp93  = x * L[23];
  double Lstmp94  = x * L[41];
  double Lstmp95  = x * L[66];
  double Lstmp96  = x * L[25];
  double Lstmp97  = x * L[44];
  double Lstmp98  = x * L[70];
  double Lstmp99  = Lstmp87 * x;
  double Lstmp100 = Lstmp88 * x;
  double Lstmp101 = Lstmp89 * x;
  double Lstmp102 = Lstmp90 * x;
  double Lstmp103 = x * L[68];
  double Lstmp104 = y * L[13];
  double Lstmp105 = Lstmp56 * y;
  double Lstmp106 = x * L[28];
  double Lstmp107 = x * L[48];
  double Lstmp108 = x * L[75];
  double Lstmp109 = y * L[23];
  double Lstmp110 = y * L[38];
  double Lstmp111 = y * L[59];
  double Lstmp112 = y * L[32];
  double Lstmp113 = y * L[53];
  double Lstmp114 = y * L[81];
  double Lstmp115 = y * L[47];
  double Lstmp116 = Lstmp115 * x;
  double Lstmp117 = y * L[74];
  double Lstmp118 = Lstmp117 * x;
  double Lstmp119 = Lstmp89 * y;
  double Lstmp120 = Lstmp92 * y;
  double Lstmp121 = y * L[68];
  double Lstmp122 = y * L[14];
  double Lstmp123 = z * L[15];
  double Lstmp124 = z * L[18];
  double Lstmp125 = z * L[28];
  double Lstmp126 = Lstmp125 * y;
  double Lstmp127 = x * L[27];
  double Lstmp128 = x * L[46];
  double Lstmp129 = x * L[72];
  double Lstmp130 = y * L[24];
  double Lstmp131 = z * L[25];
  double Lstmp132 = y * L[39];
  double Lstmp133 = z * L[40];
  double Lstmp134 = y * L[60];
  double Lstmp135 = z * L[61];
  double Lstmp136 = z * L[32];
  double Lstmp137 = z * L[52];
  double Lstmp138 = z * L[79];
  double Lstmp139 = z * L[47];
  double Lstmp140 = Lstmp139 * x;
  double Lstmp141 = z * L[73];
  double Lstmp142 = Lstmp141 * x;
  double Lstmp143 = z * L[43];
  double Lstmp144 = Lstmp143 * y;
  double Lstmp145 = z * L[64];
  double Lstmp146 = Lstmp145 * y;
  double Lstmp147 = z * L[68];
  double Lstmp148 = x * L[38];
  double Lstmp149 = x * L[62];
  double Lstmp150 = x * L[40];
  double Lstmp151 = x * L[65];
  double Lstmp152 = Lstmp91 * x;
  double Lstmp153 = Lstmp92 * x;
  double Lstmp154 = x * L[43];
  double Lstmp155 = x * L[69];
  double Lstmp156 = Lstmp121 * x;
  double Lstmp157 = x * L[42];
  double Lstmp158 = x * L[67];
  double Lstmp159 = Lstmp147 * x;
  double Lstmp160 = y * L[26];
  double Lstmp161 = Lstmp58 * y;
  double Lstmp162 = y * L[41];
  double Lstmp163 = y * L[62];
  double Lstmp164 = y * L[52];
  double Lstmp165 = y * L[80];
  double Lstmp166 = y * L[73];
  double Lstmp167 = Lstmp166 * x;
  double Lstmp168 = Lstmp90 * y;
  double Lstmp169 = y * L[27];
  double Lstmp170 = Lstmp139 * y;
  double Lstmp171 = y * L[42];
  double Lstmp172 = y * L[63];
  double Lstmp173 = Lstmp147 * y;
  double Lstmp174 = z * L[29];
  double Lstmp175 = z * L[33];
  double Lstmp176 = z * L[48];
  double Lstmp177 = Lstmp176 * y;
  double Lstmp178 = z * L[44];
  double Lstmp179 = z * L[65];
  double Lstmp180 = z * L[53];
  double Lstmp181 = z * L[80];
  double Lstmp182 = z * L[74];
  double Lstmp183 = Lstmp182 * x;
  double Lstmp184 = z * L[69];
  double Lstmp185 = Lstmp184 * y;
  double Lstmp186 = x * L[59];
  double Lstmp187 = x * L[61];
  double Lstmp188 = x * L[64];
  double Lstmp189 = x * L[63];
  double Lstmp190 = y * L[45];
  double Lstmp191 = Lstmp60 * y;
  double Lstmp192 = y * L[66];
  double Lstmp193 = y * L[79];
  double Lstmp194 = y * L[46];
  double Lstmp195 = Lstmp141 * y;
  double Lstmp196 = y * L[67];
  double Lstmp197 = Lstmp182 * y;
  double Lstmp198 = z * L[49];
  double Lstmp199 = z * L[54];
  double Lstmp200 = z * L[75];
  double Lstmp201 = Lstmp200 * y;
  double Lstmp202 = z * L[70];
  double Lstmp203 = z * L[81];
  double Lstmp204 = y * L[71];
  double Lstmp205 = y * L[72];
  double Lstmp206 = z * L[76];
  double Lstmp207 = z * L[82];
#pragma omp atomic
  Ls[0] += Lstmp0 * x + Lstmp1 * x + Lstmp10 * Lstmp38 + Lstmp10 * Lstmp39 +
           Lstmp10 * Lstmp67 + Lstmp10 * L[20] + Lstmp11 * Lstmp40 + Lstmp11 * Lstmp41 +
           Lstmp11 * L[35] + (1.0 / 48.0) * Lstmp12 * Lstmp23 * L[81] +
           Lstmp12 * Lstmp79 * L[59] + Lstmp13 * Lstmp26 + Lstmp13 * Lstmp46 +
           Lstmp13 * Lstmp57 + Lstmp13 * L[7] + (1.0 / 36.0) * Lstmp14 * Lstmp21 * L[80] +
           Lstmp14 * Lstmp77 * L[62] + Lstmp15 * Lstmp27 + Lstmp15 * Lstmp47 +
           Lstmp15 * Lstmp59 + Lstmp15 * L[16] +
           (1.0 / 48.0) * Lstmp16 * Lstmp19 * L[79] + Lstmp16 * Lstmp72 * L[66] +
           Lstmp17 * Lstmp28 + Lstmp17 * Lstmp48 + Lstmp17 * Lstmp61 + Lstmp17 * L[30] +
           Lstmp18 * Lstmp29 + Lstmp18 * Lstmp49 + Lstmp18 * L[50] +
           Lstmp19 * Lstmp79 * L[61] + Lstmp2 * y + Lstmp20 * Lstmp30 +
           Lstmp20 * Lstmp42 + Lstmp20 * Lstmp51 + Lstmp20 * L[9] +
           Lstmp21 * Lstmp77 * L[65] + Lstmp22 * Lstmp31 + Lstmp22 * Lstmp43 +
           Lstmp22 * Lstmp53 + Lstmp22 * L[19] + Lstmp23 * Lstmp72 * L[70] +
           Lstmp24 * Lstmp32 + Lstmp24 * Lstmp44 + Lstmp24 * Lstmp55 + Lstmp24 * L[34] +
           Lstmp25 * Lstmp33 + Lstmp25 * Lstmp45 + Lstmp25 * L[55] + Lstmp34 * Lstmp6 +
           Lstmp35 * Lstmp6 + Lstmp36 * Lstmp8 + Lstmp37 * Lstmp8 + Lstmp4 * x +
           (1.0 / 8.0) * Lstmp5 * Lstmp80 * L[68] + Lstmp6 * Lstmp63 + Lstmp6 * L[4] +
           Lstmp65 * Lstmp8 + Lstmp69 * Lstmp89 + Lstmp69 * L[23] + Lstmp71 * Lstmp90 +
           Lstmp71 * L[41] + Lstmp73 * Lstmp87 + Lstmp73 * L[25] + Lstmp74 * Lstmp88 +
           Lstmp74 * L[44] + Lstmp76 * Lstmp92 + Lstmp76 * L[38] + Lstmp78 * Lstmp91 +
           Lstmp78 * L[40] + Lstmp8 * L[10] + Lstmp81 * Lstmp84 + Lstmp81 * L[32] +
           Lstmp82 * Lstmp85 + Lstmp82 * L[53] + Lstmp83 * Lstmp86 + Lstmp83 * L[52] +
           (1.0 / 720.0) * (x * x * x * x * x * x) * L[56] + x * L[1] +
           (1.0 / 720.0) * (y * y * y * y * y * y) * L[77] + y * L[2] +
           (1.0 / 720.0) * (z * z * z * z * z * z) * L[83] + z * L[3] + L[0];
#pragma omp atomic
  Ls[1] += Lstmp0 + Lstmp1 + Lstmp10 * Lstmp40 + Lstmp10 * Lstmp41 + Lstmp10 * L[35] +
           Lstmp100 * Lstmp22 + Lstmp101 * Lstmp13 + Lstmp102 * Lstmp15 +
           Lstmp103 * Lstmp81 + Lstmp11 * L[56] + Lstmp13 * Lstmp56 + Lstmp13 * Lstmp93 +
           Lstmp13 * L[13] + Lstmp15 * Lstmp58 + Lstmp15 * Lstmp94 + Lstmp15 * L[26] +
           Lstmp17 * Lstmp60 + Lstmp17 * Lstmp95 + Lstmp17 * L[45] + Lstmp18 * L[71] +
           Lstmp20 * Lstmp50 + Lstmp20 * Lstmp96 + Lstmp20 * Lstmp99 + Lstmp20 * L[15] +
           Lstmp22 * Lstmp52 + Lstmp22 * Lstmp97 + Lstmp22 * L[29] + Lstmp24 * Lstmp54 +
           Lstmp24 * Lstmp98 + Lstmp24 * L[49] + Lstmp25 * L[76] + Lstmp34 * x +
           Lstmp35 * x + Lstmp36 * Lstmp6 + Lstmp37 * Lstmp6 + Lstmp38 * Lstmp8 +
           Lstmp39 * Lstmp8 + Lstmp4 + Lstmp6 * Lstmp65 + Lstmp6 * L[10] + Lstmp63 * x +
           Lstmp67 * Lstmp8 + Lstmp69 * Lstmp92 + Lstmp69 * L[38] + Lstmp71 * L[62] +
           Lstmp73 * Lstmp91 + Lstmp73 * L[40] + Lstmp74 * L[65] + Lstmp76 * L[59] +
           Lstmp78 * L[61] + Lstmp8 * L[20] + Lstmp81 * L[47] + Lstmp82 * L[74] +
           Lstmp83 * L[73] + x * L[4] + L[1];
#pragma omp atomic
  Ls[2] += Lstmp10 * Lstmp111 + Lstmp10 * Lstmp66 + Lstmp10 * L[36] + Lstmp104 * x +
           Lstmp105 * x + Lstmp106 * Lstmp20 + Lstmp107 * Lstmp22 + Lstmp108 * Lstmp24 +
           Lstmp109 * Lstmp6 + Lstmp11 * L[57] + Lstmp110 * Lstmp8 + Lstmp112 * Lstmp20 +
           Lstmp113 * Lstmp22 + Lstmp114 * Lstmp24 + Lstmp116 * Lstmp20 +
           Lstmp118 * Lstmp22 + Lstmp119 * Lstmp6 + Lstmp120 * Lstmp8 +
           Lstmp121 * Lstmp73 + Lstmp13 * Lstmp27 + Lstmp13 * Lstmp47 +
           Lstmp13 * Lstmp59 + Lstmp13 * L[16] + Lstmp15 * Lstmp28 + Lstmp15 * Lstmp48 +
           Lstmp15 * Lstmp61 + Lstmp15 * L[30] + Lstmp17 * Lstmp29 + Lstmp17 * Lstmp49 +
           Lstmp17 * L[50] + Lstmp18 * L[77] + Lstmp2 + Lstmp20 * L[18] +
           Lstmp22 * L[33] + Lstmp24 * L[54] + Lstmp25 * L[82] + Lstmp3 * x +
           Lstmp46 * y + Lstmp6 * Lstmp62 + Lstmp6 * L[11] + Lstmp64 * Lstmp8 +
           Lstmp69 * Lstmp90 + Lstmp69 * L[41] + Lstmp71 * L[66] + Lstmp73 * L[43] +
           Lstmp74 * L[69] + Lstmp76 * L[62] + Lstmp78 * L[64] + Lstmp8 * L[21] +
           Lstmp81 * Lstmp86 + Lstmp81 * L[52] + Lstmp82 * L[80] + Lstmp83 * L[79] +
           x * L[5] + y * L[7] + L[2];
#pragma omp atomic
  Ls[3] += Lstmp10 * Lstmp134 + Lstmp10 * Lstmp135 + Lstmp10 * L[37] + Lstmp11 * L[58] +
           Lstmp122 * x + Lstmp123 * x + Lstmp124 * y + Lstmp126 * x +
           Lstmp127 * Lstmp13 + Lstmp128 * Lstmp15 + Lstmp129 * Lstmp17 +
           Lstmp13 * Lstmp136 + Lstmp13 * Lstmp140 + Lstmp13 * L[17] + Lstmp130 * Lstmp6 +
           Lstmp131 * Lstmp6 + Lstmp132 * Lstmp8 + Lstmp133 * Lstmp8 +
           Lstmp137 * Lstmp15 + Lstmp138 * Lstmp17 + Lstmp142 * Lstmp15 +
           Lstmp144 * Lstmp6 + Lstmp146 * Lstmp8 + Lstmp147 * Lstmp69 + Lstmp15 * L[31] +
           Lstmp17 * L[51] + Lstmp18 * L[78] + Lstmp20 * Lstmp31 + Lstmp20 * Lstmp43 +
           Lstmp20 * Lstmp53 + Lstmp20 * L[19] + Lstmp22 * Lstmp32 + Lstmp22 * Lstmp44 +
           Lstmp22 * Lstmp55 + Lstmp22 * L[34] + Lstmp24 * Lstmp33 + Lstmp24 * Lstmp45 +
           Lstmp24 * L[55] + Lstmp25 * L[83] + Lstmp6 * L[12] + Lstmp69 * L[42] +
           Lstmp71 * L[67] + Lstmp73 * Lstmp88 + Lstmp73 * L[44] + Lstmp74 * L[70] +
           Lstmp76 * L[63] + Lstmp78 * L[65] + Lstmp8 * L[22] + Lstmp81 * Lstmp85 +
           Lstmp81 * L[53] + Lstmp82 * L[81] + Lstmp83 * L[80] + x * L[6] + y * L[8] +
           z * L[9] + L[3];
#pragma omp atomic
  Ls[4] += Lstmp10 * L[56] + Lstmp13 * Lstmp148 + Lstmp13 * Lstmp153 + Lstmp13 * Lstmp89 +
           Lstmp13 * L[23] + Lstmp149 * Lstmp15 + Lstmp15 * Lstmp90 + Lstmp15 * L[41] +
           Lstmp150 * Lstmp20 + Lstmp151 * Lstmp22 + Lstmp152 * Lstmp20 +
           Lstmp17 * L[66] + Lstmp20 * Lstmp87 + Lstmp20 * L[25] + Lstmp22 * Lstmp88 +
           Lstmp22 * L[44] + Lstmp24 * L[70] + Lstmp34 + Lstmp35 + Lstmp36 * x +
           Lstmp37 * x + Lstmp38 * Lstmp6 + Lstmp39 * Lstmp6 + Lstmp40 * Lstmp8 +
           Lstmp41 * Lstmp8 + Lstmp6 * Lstmp67 + Lstmp6 * L[20] + Lstmp63 + Lstmp65 * x +
           Lstmp69 * L[59] + Lstmp73 * L[61] + Lstmp8 * L[35] + Lstmp81 * L[68] +
           x * L[10] + L[4];
#pragma omp atomic
  Ls[5] += Lstmp10 * L[57] + Lstmp102 * Lstmp13 + Lstmp104 + Lstmp105 + Lstmp109 * x +
           Lstmp110 * Lstmp6 + Lstmp111 * Lstmp8 + Lstmp115 * Lstmp20 +
           Lstmp117 * Lstmp22 + Lstmp119 * x + Lstmp120 * Lstmp6 + Lstmp13 * Lstmp58 +
           Lstmp13 * Lstmp94 + Lstmp13 * L[26] + Lstmp15 * Lstmp60 + Lstmp15 * Lstmp95 +
           Lstmp15 * L[45] + Lstmp154 * Lstmp20 + Lstmp155 * Lstmp22 +
           Lstmp156 * Lstmp20 + Lstmp17 * L[71] + Lstmp20 * L[28] + Lstmp22 * L[48] +
           Lstmp24 * L[75] + Lstmp3 + Lstmp6 * Lstmp64 + Lstmp6 * L[21] + Lstmp62 * x +
           Lstmp66 * Lstmp8 + Lstmp69 * L[62] + Lstmp73 * L[64] + Lstmp8 * L[36] +
           Lstmp81 * L[73] + x * L[11] + L[5];
#pragma omp atomic
  Ls[6] += Lstmp10 * L[58] + Lstmp100 * Lstmp20 + Lstmp122 + Lstmp123 + Lstmp126 +
           Lstmp13 * Lstmp139 + Lstmp13 * Lstmp157 + Lstmp13 * Lstmp159 +
           Lstmp13 * L[27] + Lstmp130 * x + Lstmp131 * x + Lstmp132 * Lstmp6 +
           Lstmp133 * Lstmp6 + Lstmp134 * Lstmp8 + Lstmp135 * Lstmp8 +
           Lstmp141 * Lstmp15 + Lstmp144 * x + Lstmp146 * Lstmp6 + Lstmp15 * Lstmp158 +
           Lstmp15 * L[46] + Lstmp17 * L[72] + Lstmp20 * Lstmp52 + Lstmp20 * Lstmp97 +
           Lstmp20 * L[29] + Lstmp22 * Lstmp54 + Lstmp22 * Lstmp98 + Lstmp22 * L[49] +
           Lstmp24 * L[76] + Lstmp6 * L[22] + Lstmp69 * L[63] + Lstmp73 * L[65] +
           Lstmp8 * L[37] + Lstmp81 * L[74] + x * L[12] + L[6];
#pragma omp atomic
  Ls[7] += Lstmp10 * L[59] + Lstmp13 * Lstmp28 + Lstmp13 * Lstmp48 + Lstmp13 * Lstmp61 +
           Lstmp13 * L[30] + Lstmp15 * Lstmp29 + Lstmp15 * Lstmp49 + Lstmp15 * L[50] +
           Lstmp160 * x + Lstmp161 * x + Lstmp162 * Lstmp6 + Lstmp163 * Lstmp8 +
           Lstmp164 * Lstmp20 + Lstmp165 * Lstmp22 + Lstmp167 * Lstmp20 +
           Lstmp168 * Lstmp6 + Lstmp17 * L[77] + Lstmp20 * Lstmp84 + Lstmp20 * L[32] +
           Lstmp22 * Lstmp85 + Lstmp22 * L[53] + Lstmp24 * L[81] + Lstmp26 + Lstmp46 +
           Lstmp47 * y + Lstmp57 + Lstmp6 * Lstmp89 + Lstmp6 * L[23] + Lstmp69 * L[66] +
           Lstmp73 * L[68] + Lstmp8 * Lstmp92 + Lstmp8 * L[38] + Lstmp81 * L[79] +
           y * L[16] + L[7];
#pragma omp atomic
  Ls[8] += Lstmp10 * L[60] + Lstmp107 * Lstmp20 + Lstmp108 * Lstmp22 +
           Lstmp113 * Lstmp20 + Lstmp114 * Lstmp22 + Lstmp118 * Lstmp20 + Lstmp124 +
           Lstmp125 * x + Lstmp128 * Lstmp13 + Lstmp129 * Lstmp15 + Lstmp13 * Lstmp137 +
           Lstmp13 * Lstmp142 + Lstmp13 * L[31] + Lstmp136 * y + Lstmp138 * Lstmp15 +
           Lstmp143 * Lstmp6 + Lstmp145 * Lstmp8 + Lstmp15 * L[51] + Lstmp169 * x +
           Lstmp17 * L[78] + Lstmp170 * x + Lstmp171 * Lstmp6 + Lstmp172 * Lstmp8 +
           Lstmp173 * Lstmp6 + Lstmp20 * L[33] + Lstmp22 * L[54] + Lstmp24 * L[82] +
           Lstmp6 * L[24] + Lstmp69 * L[67] + Lstmp73 * L[69] + Lstmp8 * L[39] +
           Lstmp81 * L[80] + x * L[14] + y * L[17] + L[8];
#pragma omp atomic
  Ls[9] += Lstmp10 * L[61] + Lstmp13 * Lstmp180 + Lstmp13 * Lstmp183 + Lstmp13 * Lstmp84 +
           Lstmp13 * L[32] + Lstmp15 * Lstmp181 + Lstmp15 * Lstmp86 + Lstmp15 * L[52] +
           Lstmp17 * L[79] + Lstmp174 * x + Lstmp175 * y + Lstmp177 * x +
           Lstmp178 * Lstmp6 + Lstmp179 * Lstmp8 + Lstmp185 * Lstmp6 + Lstmp20 * Lstmp32 +
           Lstmp20 * Lstmp44 + Lstmp20 * Lstmp55 + Lstmp20 * L[34] + Lstmp22 * Lstmp33 +
           Lstmp22 * Lstmp45 + Lstmp22 * L[55] + Lstmp24 * L[83] + Lstmp30 + Lstmp42 +
           Lstmp51 + Lstmp6 * Lstmp87 + Lstmp6 * L[25] + Lstmp69 * L[68] +
           Lstmp73 * L[70] + Lstmp8 * Lstmp91 + Lstmp8 * L[40] + Lstmp81 * L[81] +
           z * L[19] + L[9];
#pragma omp atomic
  Ls[10] += Lstmp13 * Lstmp186 + Lstmp13 * Lstmp92 + Lstmp13 * L[38] + Lstmp15 * L[62] +
            Lstmp187 * Lstmp20 + Lstmp20 * Lstmp91 + Lstmp20 * L[40] + Lstmp22 * L[65] +
            Lstmp36 + Lstmp37 + Lstmp38 * x + Lstmp39 * x + Lstmp40 * Lstmp6 +
            Lstmp41 * Lstmp6 + Lstmp6 * L[35] + Lstmp65 + Lstmp67 * x + Lstmp8 * L[56] +
            x * L[20] + L[10];
#pragma omp atomic
  Ls[11] += Lstmp109 + Lstmp110 * x + Lstmp111 * Lstmp6 + Lstmp119 + Lstmp120 * x +
            Lstmp121 * Lstmp20 + Lstmp13 * Lstmp149 + Lstmp13 * Lstmp90 +
            Lstmp13 * L[41] + Lstmp15 * L[66] + Lstmp188 * Lstmp20 + Lstmp20 * L[43] +
            Lstmp22 * L[69] + Lstmp6 * Lstmp66 + Lstmp6 * L[36] + Lstmp62 + Lstmp64 * x +
            Lstmp8 * L[57] + x * L[21] + L[11];
#pragma omp atomic
  Ls[12] += Lstmp13 * Lstmp147 + Lstmp13 * Lstmp189 + Lstmp13 * L[42] + Lstmp130 +
            Lstmp131 + Lstmp132 * x + Lstmp133 * x + Lstmp134 * Lstmp6 +
            Lstmp135 * Lstmp6 + Lstmp144 + Lstmp146 * x + Lstmp15 * L[67] +
            Lstmp151 * Lstmp20 + Lstmp20 * Lstmp88 + Lstmp20 * L[44] + Lstmp22 * L[70] +
            Lstmp6 * L[37] + Lstmp8 * L[58] + x * L[22] + L[12];
#pragma omp atomic
  Ls[13] += Lstmp101 + Lstmp103 * Lstmp20 + Lstmp13 * Lstmp60 + Lstmp13 * Lstmp95 +
            Lstmp13 * L[45] + Lstmp15 * L[71] + Lstmp160 + Lstmp161 + Lstmp162 * x +
            Lstmp163 * Lstmp6 + Lstmp166 * Lstmp20 + Lstmp168 * x + Lstmp20 * L[47] +
            Lstmp22 * L[74] + Lstmp56 + Lstmp6 * Lstmp92 + Lstmp6 * L[38] +
            Lstmp8 * L[59] + Lstmp93 + L[13];
#pragma omp atomic
  Ls[14] += Lstmp117 * Lstmp20 + Lstmp125 + Lstmp13 * Lstmp141 + Lstmp13 * Lstmp158 +
            Lstmp13 * L[46] + Lstmp143 * x + Lstmp145 * Lstmp6 + Lstmp15 * L[72] +
            Lstmp155 * Lstmp20 + Lstmp169 + Lstmp170 + Lstmp171 * x + Lstmp172 * Lstmp6 +
            Lstmp173 * x + Lstmp20 * L[48] + Lstmp22 * L[75] + Lstmp6 * L[39] +
            Lstmp8 * L[60] + x * L[24] + L[14];
#pragma omp atomic
  Ls[15] += Lstmp103 * Lstmp13 + Lstmp13 * Lstmp182 + Lstmp13 * L[47] + Lstmp15 * L[73] +
            Lstmp174 + Lstmp177 + Lstmp178 * x + Lstmp179 * Lstmp6 + Lstmp185 * x +
            Lstmp20 * Lstmp54 + Lstmp20 * Lstmp98 + Lstmp20 * L[49] + Lstmp22 * L[76] +
            Lstmp50 + Lstmp6 * Lstmp91 + Lstmp6 * L[40] + Lstmp8 * L[61] + Lstmp96 +
            Lstmp99 + L[15];
#pragma omp atomic
  Ls[16] += Lstmp13 * Lstmp29 + Lstmp13 * Lstmp49 + Lstmp13 * L[50] + Lstmp15 * L[77] +
            Lstmp190 * x + Lstmp191 * x + Lstmp192 * Lstmp6 + Lstmp193 * Lstmp20 +
            Lstmp20 * Lstmp86 + Lstmp20 * L[52] + Lstmp22 * L[80] + Lstmp27 + Lstmp47 +
            Lstmp48 * y + Lstmp59 + Lstmp6 * Lstmp90 + Lstmp6 * L[41] + Lstmp8 * L[62] +
            y * L[30] + L[16];
#pragma omp atomic
  Ls[17] += Lstmp127 + Lstmp129 * Lstmp13 + Lstmp13 * Lstmp138 + Lstmp13 * L[51] +
            Lstmp136 + Lstmp137 * y + Lstmp140 + Lstmp147 * Lstmp6 + Lstmp15 * L[78] +
            Lstmp165 * Lstmp20 + Lstmp194 * x + Lstmp195 * x + Lstmp196 * Lstmp6 +
            Lstmp20 * Lstmp85 + Lstmp20 * L[53] + Lstmp22 * L[81] + Lstmp6 * L[42] +
            Lstmp8 * L[63] + y * L[31] + L[17];
#pragma omp atomic
  Ls[18] += Lstmp106 + Lstmp108 * Lstmp20 + Lstmp112 + Lstmp114 * Lstmp20 + Lstmp116 +
            Lstmp121 * Lstmp6 + Lstmp13 * Lstmp181 + Lstmp13 * Lstmp86 + Lstmp13 * L[52] +
            Lstmp15 * L[79] + Lstmp175 + Lstmp176 * x + Lstmp180 * y + Lstmp184 * Lstmp6 +
            Lstmp197 * x + Lstmp20 * L[54] + Lstmp22 * L[82] + Lstmp6 * L[43] +
            Lstmp8 * L[64] + L[18];
#pragma omp atomic
  Ls[19] += Lstmp13 * Lstmp203 + Lstmp13 * Lstmp85 + Lstmp13 * L[53] + Lstmp15 * L[80] +
            Lstmp198 * x + Lstmp199 * y + Lstmp20 * Lstmp33 + Lstmp20 * Lstmp45 +
            Lstmp20 * L[55] + Lstmp201 * x + Lstmp202 * Lstmp6 + Lstmp22 * L[83] +
            Lstmp31 + Lstmp43 + Lstmp53 + Lstmp6 * Lstmp88 + Lstmp6 * L[44] +
            Lstmp8 * L[65] + z * L[34] + L[19];
#pragma omp atomic
  Ls[20] += Lstmp13 * L[59] + Lstmp20 * L[61] + Lstmp38 + Lstmp39 + Lstmp40 * x +
            Lstmp41 * x + Lstmp6 * L[56] + Lstmp67 + x * L[35] + L[20];
#pragma omp atomic
  Ls[21] += Lstmp110 + Lstmp111 * x + Lstmp120 + Lstmp13 * L[62] + Lstmp20 * L[64] +
            Lstmp6 * L[57] + Lstmp64 + Lstmp66 * x + x * L[36] + L[21];
#pragma omp atomic
  Ls[22] += Lstmp13 * L[63] + Lstmp132 + Lstmp133 + Lstmp134 * x + Lstmp135 * x +
            Lstmp146 + Lstmp20 * L[65] + Lstmp6 * L[58] + x * L[37] + L[22];
#pragma omp atomic
  Ls[23] += Lstmp13 * L[66] + Lstmp148 + Lstmp153 + Lstmp162 + Lstmp163 * x + Lstmp168 +
            Lstmp20 * L[68] + Lstmp6 * L[59] + Lstmp89 + L[23];
#pragma omp atomic
  Ls[24] += Lstmp13 * L[67] + Lstmp143 + Lstmp145 * x + Lstmp171 + Lstmp172 * x +
            Lstmp173 + Lstmp20 * L[69] + Lstmp6 * L[60] + x * L[39] + L[24];
#pragma omp atomic
  Ls[25] += Lstmp13 * L[68] + Lstmp150 + Lstmp152 + Lstmp178 + Lstmp179 * x + Lstmp185 +
            Lstmp20 * L[70] + Lstmp6 * L[61] + Lstmp87 + L[25];
#pragma omp atomic
  Ls[26] += Lstmp102 + Lstmp13 * L[71] + Lstmp190 + Lstmp191 + Lstmp192 * x +
            Lstmp20 * L[73] + Lstmp58 + Lstmp6 * L[62] + Lstmp94 + L[26];
#pragma omp atomic
  Ls[27] += Lstmp13 * L[72] + Lstmp139 + Lstmp157 + Lstmp159 + Lstmp194 + Lstmp195 +
            Lstmp196 * x + Lstmp20 * L[74] + Lstmp6 * L[63] + L[27];
#pragma omp atomic
  Ls[28] += Lstmp115 + Lstmp13 * L[73] + Lstmp154 + Lstmp156 + Lstmp176 + Lstmp184 * x +
            Lstmp197 + Lstmp20 * L[75] + Lstmp6 * L[64] + L[28];
#pragma omp atomic
  Ls[29] += Lstmp100 + Lstmp13 * L[74] + Lstmp198 + Lstmp20 * L[76] + Lstmp201 +
            Lstmp202 * x + Lstmp52 + Lstmp6 * L[65] + Lstmp97 + L[29];
#pragma omp atomic
  Ls[30] += Lstmp13 * L[77] + Lstmp20 * L[79] + Lstmp204 * x + Lstmp28 + Lstmp48 +
            Lstmp49 * y + Lstmp6 * L[66] + Lstmp61 + y * L[50] + L[30];
#pragma omp atomic
  Ls[31] += Lstmp128 + Lstmp13 * L[78] + Lstmp137 + Lstmp138 * y + Lstmp142 +
            Lstmp20 * L[80] + Lstmp205 * x + Lstmp6 * L[67] + y * L[51] + L[31];
#pragma omp atomic
  Ls[32] += Lstmp13 * L[79] + Lstmp164 + Lstmp167 + Lstmp180 + Lstmp181 * y + Lstmp183 +
            Lstmp20 * L[81] + Lstmp6 * L[68] + Lstmp84 + L[32];
#pragma omp atomic
  Ls[33] += Lstmp107 + Lstmp113 + Lstmp118 + Lstmp13 * L[80] + Lstmp199 +
            Lstmp20 * L[82] + Lstmp200 * x + Lstmp203 * y + Lstmp6 * L[69] + L[33];
#pragma omp atomic
  Ls[34] += Lstmp13 * L[81] + Lstmp20 * L[83] + Lstmp206 * x + Lstmp207 * y + Lstmp32 +
            Lstmp44 + Lstmp55 + Lstmp6 * L[70] + z * L[55] + L[34];
#pragma omp atomic
  Ls[35] += Lstmp40 + Lstmp41 + x * L[56] + L[35];
#pragma omp atomic
  Ls[36] += Lstmp111 + Lstmp66 + x * L[57] + L[36];
#pragma omp atomic
  Ls[37] += Lstmp134 + Lstmp135 + x * L[58] + L[37];
#pragma omp atomic
  Ls[38] += Lstmp163 + Lstmp186 + Lstmp92 + L[38];
#pragma omp atomic
  Ls[39] += Lstmp145 + Lstmp172 + x * L[60] + L[39];
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
  Ls[50] += Lstmp29 + Lstmp49 + y * L[77] + L[50];
#pragma omp atomic
  Ls[51] += Lstmp129 + Lstmp138 + y * L[78] + L[51];
#pragma omp atomic
  Ls[52] += Lstmp181 + Lstmp193 + Lstmp86 + L[52];
#pragma omp atomic
  Ls[53] += Lstmp165 + Lstmp203 + Lstmp85 + L[53];
#pragma omp atomic
  Ls[54] += Lstmp108 + Lstmp114 + Lstmp207 + L[54];
#pragma omp atomic
  Ls[55] += Lstmp33 + Lstmp45 + z * L[83] + L[55];
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

void field_m1_L2P_7(double x, double y, double z, double* L, double* F) {
  double Ftmp0  = x * y;
  double Ftmp1  = x * z;
  double Ftmp2  = y * z;
  double Ftmp3  = Ftmp0 * z;
  double Ftmp4  = (x * x);
  double Ftmp5  = (1.0 / 2.0) * Ftmp4;
  double Ftmp6  = (x * x * x);
  double Ftmp7  = (1.0 / 6.0) * Ftmp6;
  double Ftmp8  = (1.0 / 24.0) * (x * x * x * x);
  double Ftmp9  = (1.0 / 120.0) * (x * x * x * x * x);
  double Ftmp10 = (y * y);
  double Ftmp11 = (1.0 / 2.0) * Ftmp10;
  double Ftmp12 = (y * y * y);
  double Ftmp13 = (1.0 / 6.0) * Ftmp12;
  double Ftmp14 = (1.0 / 24.0) * (y * y * y * y);
  double Ftmp15 = (1.0 / 120.0) * (y * y * y * y * y);
  double Ftmp16 = (z * z);
  double Ftmp17 = (1.0 / 2.0) * Ftmp16;
  double Ftmp18 = (z * z * z);
  double Ftmp19 = (1.0 / 6.0) * Ftmp18;
  double Ftmp20 = (1.0 / 24.0) * (z * z * z * z);
  double Ftmp21 = (1.0 / 120.0) * (z * z * z * z * z);
  double Ftmp22 = Ftmp11 * x;
  double Ftmp23 = Ftmp13 * x;
  double Ftmp24 = Ftmp14 * x;
  double Ftmp25 = Ftmp17 * x;
  double Ftmp26 = Ftmp19 * x;
  double Ftmp27 = Ftmp20 * x;
  double Ftmp28 = Ftmp5 * y;
  double Ftmp29 = Ftmp5 * z;
  double Ftmp30 = Ftmp7 * y;
  double Ftmp31 = Ftmp7 * z;
  double Ftmp32 = Ftmp8 * y;
  double Ftmp33 = Ftmp8 * z;
  double Ftmp34 = Ftmp17 * y;
  double Ftmp35 = Ftmp19 * y;
  double Ftmp36 = Ftmp20 * y;
  double Ftmp37 = Ftmp11 * z;
  double Ftmp38 = Ftmp13 * z;
  double Ftmp39 = Ftmp14 * z;
  double Ftmp40 = Ftmp0 * Ftmp17;
  double Ftmp41 = Ftmp0 * Ftmp19;
  double Ftmp42 = Ftmp1 * Ftmp11;
  double Ftmp43 = Ftmp1 * Ftmp13;
  double Ftmp44 = Ftmp2 * Ftmp5;
  double Ftmp45 = Ftmp2 * Ftmp7;
  double Ftmp46 = (1.0 / 4.0) * Ftmp4;
  double Ftmp47 = Ftmp10 * Ftmp46;
  double Ftmp48 = (1.0 / 12.0) * Ftmp4;
  double Ftmp49 = Ftmp12 * Ftmp48;
  double Ftmp50 = Ftmp16 * Ftmp46;
  double Ftmp51 = Ftmp18 * Ftmp48;
  double Ftmp52 = (1.0 / 12.0) * Ftmp6;
  double Ftmp53 = Ftmp10 * Ftmp52;
  double Ftmp54 = Ftmp16 * Ftmp52;
  double Ftmp55 = (1.0 / 4.0) * Ftmp10 * Ftmp16;
  double Ftmp56 = (1.0 / 12.0) * Ftmp10 * Ftmp18;
  double Ftmp57 = (1.0 / 12.0) * Ftmp12 * Ftmp16;
  double Ftmp58 = Ftmp55 * x;
  double Ftmp59 = Ftmp50 * y;
  double Ftmp60 = Ftmp47 * z;
#pragma omp atomic
  F[0] += -Ftmp0 * L[11] - Ftmp1 * L[12] - Ftmp11 * L[13] - Ftmp13 * L[26] -
          Ftmp14 * L[45] - Ftmp15 * L[71] - Ftmp17 * L[15] - Ftmp19 * L[29] -
          Ftmp2 * L[14] - Ftmp20 * L[49] - Ftmp21 * L[76] - Ftmp22 * L[23] -
          Ftmp23 * L[41] - Ftmp24 * L[66] - Ftmp25 * L[25] - Ftmp26 * L[44] -
          Ftmp27 * L[70] - Ftmp28 * L[21] - Ftmp29 * L[22] - Ftmp3 * L[24] -
          Ftmp30 * L[36] - Ftmp31 * L[37] - Ftmp32 * L[57] - Ftmp33 * L[58] -
          Ftmp34 * L[28] - Ftmp35 * L[48] - Ftmp36 * L[75] - Ftmp37 * L[27] -
          Ftmp38 * L[46] - Ftmp39 * L[72] - Ftmp40 * L[43] - Ftmp41 * L[69] -
          Ftmp42 * L[42] - Ftmp43 * L[67] - Ftmp44 * L[39] - Ftmp45 * L[60] -
          Ftmp47 * L[38] - Ftmp49 * L[62] - Ftmp5 * L[10] - Ftmp50 * L[40] -
          Ftmp51 * L[65] - Ftmp53 * L[59] - Ftmp54 * L[61] - Ftmp55 * L[47] -
          Ftmp56 * L[74] - Ftmp57 * L[73] - Ftmp58 * L[68] - Ftmp59 * L[64] -
          Ftmp60 * L[63] - Ftmp7 * L[20] - Ftmp8 * L[35] - Ftmp9 * L[56] - x * L[4] -
          y * L[5] - z * L[6] - L[1];
#pragma omp atomic
  F[1] += -Ftmp0 * L[13] - Ftmp1 * L[14] - Ftmp11 * L[16] - Ftmp13 * L[30] -
          Ftmp14 * L[50] - Ftmp15 * L[77] - Ftmp17 * L[18] - Ftmp19 * L[33] -
          Ftmp2 * L[17] - Ftmp20 * L[54] - Ftmp21 * L[82] - Ftmp22 * L[26] -
          Ftmp23 * L[45] - Ftmp24 * L[71] - Ftmp25 * L[28] - Ftmp26 * L[48] -
          Ftmp27 * L[75] - Ftmp28 * L[23] - Ftmp29 * L[24] - Ftmp3 * L[27] -
          Ftmp30 * L[38] - Ftmp31 * L[39] - Ftmp32 * L[59] - Ftmp33 * L[60] -
          Ftmp34 * L[32] - Ftmp35 * L[53] - Ftmp36 * L[81] - Ftmp37 * L[31] -
          Ftmp38 * L[51] - Ftmp39 * L[78] - Ftmp40 * L[47] - Ftmp41 * L[74] -
          Ftmp42 * L[46] - Ftmp43 * L[72] - Ftmp44 * L[42] - Ftmp45 * L[63] -
          Ftmp47 * L[41] - Ftmp49 * L[66] - Ftmp5 * L[11] - Ftmp50 * L[43] -
          Ftmp51 * L[69] - Ftmp53 * L[62] - Ftmp54 * L[64] - Ftmp55 * L[52] -
          Ftmp56 * L[80] - Ftmp57 * L[79] - Ftmp58 * L[73] - Ftmp59 * L[68] -
          Ftmp60 * L[67] - Ftmp7 * L[21] - Ftmp8 * L[36] - Ftmp9 * L[57] - x * L[5] -
          y * L[7] - z * L[8] - L[2];
#pragma omp atomic
  F[2] += -Ftmp0 * L[14] - Ftmp1 * L[15] - Ftmp11 * L[17] - Ftmp13 * L[31] -
          Ftmp14 * L[51] - Ftmp15 * L[78] - Ftmp17 * L[19] - Ftmp19 * L[34] -
          Ftmp2 * L[18] - Ftmp20 * L[55] - Ftmp21 * L[83] - Ftmp22 * L[27] -
          Ftmp23 * L[46] - Ftmp24 * L[72] - Ftmp25 * L[29] - Ftmp26 * L[49] -
          Ftmp27 * L[76] - Ftmp28 * L[24] - Ftmp29 * L[25] - Ftmp3 * L[28] -
          Ftmp30 * L[39] - Ftmp31 * L[40] - Ftmp32 * L[60] - Ftmp33 * L[61] -
          Ftmp34 * L[33] - Ftmp35 * L[54] - Ftmp36 * L[82] - Ftmp37 * L[32] -
          Ftmp38 * L[52] - Ftmp39 * L[79] - Ftmp40 * L[48] - Ftmp41 * L[75] -
          Ftmp42 * L[47] - Ftmp43 * L[73] - Ftmp44 * L[43] - Ftmp45 * L[64] -
          Ftmp47 * L[42] - Ftmp49 * L[67] - Ftmp5 * L[12] - Ftmp50 * L[44] -
          Ftmp51 * L[70] - Ftmp53 * L[63] - Ftmp54 * L[65] - Ftmp55 * L[53] -
          Ftmp56 * L[81] - Ftmp57 * L[80] - Ftmp58 * L[74] - Ftmp59 * L[69] -
          Ftmp60 * L[68] - Ftmp7 * L[22] - Ftmp8 * L[37] - Ftmp9 * L[58] - x * L[6] -
          y * L[8] - z * L[9] - L[3];
}

void field_m1_M2P_7(double x, double y, double z, double* M, double* F) {
  double R       = sqrt(x * x + y * y + z * z);
  double Ftmp0   = 1.0 * pow(R, -3.0);
  double Ftmp1   = pow(R, -5.0);
  double Ftmp2   = 3.0 * Ftmp1;
  double Ftmp3   = Ftmp2 * M[4];
  double Ftmp4   = Ftmp2 * z;
  double Ftmp5   = y * z;
  double Ftmp6   = pow(R, -7.0);
  double Ftmp7   = 15.0 * Ftmp6;
  double Ftmp8   = Ftmp7 * M[13];
  double Ftmp9   = x * y;
  double Ftmp10  = Ftmp2 * M[1];
  double Ftmp11  = Ftmp4 * M[2];
  double Ftmp12  = (x * x);
  double Ftmp13  = Ftmp2 * M[0];
  double Ftmp14  = Ftmp5 * x;
  double Ftmp15  = Ftmp14 * Ftmp7;
  double Ftmp16  = Ftmp12 * Ftmp7;
  double Ftmp17  = pow(R, -9.0);
  double Ftmp18  = 105.0 * Ftmp17;
  double Ftmp19  = Ftmp12 * Ftmp18;
  double Ftmp20  = -9.0 * Ftmp1;
  double Ftmp21  = Ftmp16 + Ftmp20;
  double Ftmp22  = -Ftmp2;
  double Ftmp23  = (y * y);
  double Ftmp24  = Ftmp23 * Ftmp7;
  double Ftmp25  = Ftmp22 + Ftmp24;
  double Ftmp26  = (z * z);
  double Ftmp27  = Ftmp26 * Ftmp7;
  double Ftmp28  = Ftmp22 + Ftmp27;
  double Ftmp29  = Ftmp25 * M[6];
  double Ftmp30  = Ftmp28 * M[8];
  double Ftmp31  = 45.0 * Ftmp6;
  double Ftmp32  = -Ftmp31;
  double Ftmp33  = Ftmp19 + Ftmp32;
  double Ftmp34  = Ftmp33 * M[20];
  double Ftmp35  = 1.0 * y;
  double Ftmp36  = Ftmp18 * Ftmp23;
  double Ftmp37  = Ftmp32 + Ftmp36;
  double Ftmp38  = Ftmp37 * M[25];
  double Ftmp39  = 3.0 * y;
  double Ftmp40  = 35.0 * Ftmp17;
  double Ftmp41  = (Ftmp26 * Ftmp40 - 5.0 * Ftmp6) * M[27];
  double Ftmp42  = Ftmp33 * M[21];
  double Ftmp43  = 1.0 * z;
  double Ftmp44  = -Ftmp7;
  double Ftmp45  = Ftmp36 + Ftmp44;
  double Ftmp46  = Ftmp45 * M[26];
  double Ftmp47  = Ftmp18 * Ftmp26;
  double Ftmp48  = Ftmp32 + Ftmp47;
  double Ftmp49  = Ftmp43 * Ftmp48;
  double Ftmp50  = 315.0 * Ftmp17;
  double Ftmp51  = -Ftmp50;
  double Ftmp52  = pow(R, -11.0);
  double Ftmp53  = 945.0 * Ftmp52;
  double Ftmp54  = Ftmp12 * Ftmp53;
  double Ftmp55  = Ftmp51 + Ftmp54;
  double Ftmp56  = Ftmp55 * M[38];
  double Ftmp57  = Ftmp33 * Ftmp9;
  double Ftmp58  = Ftmp37 * M[15];
  double Ftmp59  = Ftmp44 + Ftmp47;
  double Ftmp60  = Ftmp59 * M[17];
  double Ftmp61  = Ftmp35 * x;
  double Ftmp62  = x * z;
  double Ftmp63  = Ftmp33 * Ftmp62;
  double Ftmp64  = Ftmp45 * M[16];
  double Ftmp65  = Ftmp48 * M[18];
  double Ftmp66  = Ftmp35 * z;
  double Ftmp67  = Ftmp23 * Ftmp53;
  double Ftmp68  = Ftmp51 + Ftmp67;
  double Ftmp69  = Ftmp68 * M[45];
  double Ftmp70  = Ftmp39 * z;
  double Ftmp71  = -Ftmp18;
  double Ftmp72  = Ftmp26 * Ftmp52;
  double Ftmp73  = 315.0 * Ftmp72;
  double Ftmp74  = Ftmp71 + Ftmp73;
  double Ftmp75  = Ftmp74 * M[47];
  double Ftmp76  = Ftmp68 * M[30];
  double Ftmp77  = -75.0 * Ftmp6;
  double Ftmp78  = 1.0 * Ftmp12;
  double Ftmp79  = Ftmp45 * M[12];
  double Ftmp80  = Ftmp59 * M[14];
  double Ftmp81  = 525.0 * Ftmp17;
  double Ftmp82  = -Ftmp81;
  double Ftmp83  = Ftmp54 + Ftmp82;
  double Ftmp84  = Ftmp12 * y;
  double Ftmp85  = Ftmp12 * z;
  double Ftmp86  = Ftmp26 * Ftmp53;
  double Ftmp87  = Ftmp51 + Ftmp86;
  double Ftmp88  = Ftmp35 * M[32];
  double Ftmp89  = Ftmp12 * Ftmp35;
  double Ftmp90  = Ftmp68 * M[25];
  double Ftmp91  = Ftmp12 * Ftmp39;
  double Ftmp92  = (-Ftmp40 + Ftmp73) * M[27];
  double Ftmp93  = Ftmp12 * Ftmp43;
  double Ftmp94  = (Ftmp67 + Ftmp71) * M[26];
  double Ftmp95  = Ftmp87 * M[28];
  double Ftmp96  = 4725.0 * Ftmp52;
  double Ftmp97  = -Ftmp96;
  double Ftmp98  = pow(R, -13.0);
  double Ftmp99  = 10395.0 * Ftmp98;
  double Ftmp100 = Ftmp12 * Ftmp99;
  double Ftmp101 = Ftmp100 + Ftmp97;
  double Ftmp102 = Ftmp12 * Ftmp5;
  double Ftmp103 = Ftmp35 * Ftmp85;
  double Ftmp104 = 2835.0 * Ftmp52;
  double Ftmp105 = -Ftmp104;
  double Ftmp106 = Ftmp23 * Ftmp99;
  double Ftmp107 = (Ftmp105 + Ftmp106) * M[45];
  double Ftmp108 = Ftmp39 * Ftmp85;
  double Ftmp109 = 3465.0 * Ftmp98;
  double Ftmp110 = Ftmp109 * Ftmp26;
  double Ftmp111 = (Ftmp110 - Ftmp53) * M[47];
  double Ftmp112 = 225.0 * Ftmp6;
  double Ftmp113 = (x * x * x * x);
  double Ftmp114 = Ftmp113 * Ftmp53;
  double Ftmp115 = 1050.0 * Ftmp17;
  double Ftmp116 = Ftmp112 + Ftmp114 - Ftmp115 * Ftmp12;
  double Ftmp117 = (y * y * y * y);
  double Ftmp118 = Ftmp117 * Ftmp53;
  double Ftmp119 = 630.0 * Ftmp17;
  double Ftmp120 = Ftmp118 - Ftmp119 * Ftmp23 + Ftmp31;
  double Ftmp121 = (z * z * z * z);
  double Ftmp122 = Ftmp121 * Ftmp53;
  double Ftmp123 = -Ftmp119 * Ftmp26 + Ftmp122 + Ftmp31;
  double Ftmp124 = Ftmp120 * M[29];
  double Ftmp125 = Ftmp123 * M[33];
  double Ftmp126 = 1575.0 * Ftmp17;
  double Ftmp127 = Ftmp113 * Ftmp98;
  double Ftmp128 = 10395.0 * Ftmp127;
  double Ftmp129 = Ftmp12 * Ftmp52;
  double Ftmp130 = 9450.0 * Ftmp129;
  double Ftmp131 = Ftmp126 + Ftmp128 - Ftmp130;
  double Ftmp132 = Ftmp131 * M[56];
  double Ftmp133 = Ftmp117 * Ftmp99;
  double Ftmp134 = 9450.0 * Ftmp52;
  double Ftmp135 = Ftmp134 * Ftmp23;
  double Ftmp136 = Ftmp126 + Ftmp133 - Ftmp135;
  double Ftmp137 = Ftmp136 * M[70];
  double Ftmp138 = (Ftmp109 * Ftmp121 + Ftmp18 - 1890.0 * Ftmp72) * M[74];
  double Ftmp139 = Ftmp131 * M[57];
  double Ftmp140 = 5670.0 * Ftmp52;
  double Ftmp141 = Ftmp140 * Ftmp23;
  double Ftmp142 = Ftmp133 - Ftmp141 + Ftmp50;
  double Ftmp143 = Ftmp142 * M[71];
  double Ftmp144 = Ftmp121 * Ftmp99;
  double Ftmp145 = Ftmp134 * Ftmp26;
  double Ftmp146 = Ftmp126 + Ftmp144 - Ftmp145;
  double Ftmp147 = Ftmp146 * Ftmp43;
  double Ftmp148 = 14175.0 * Ftmp52;
  double Ftmp149 = pow(R, -15.0);
  double Ftmp150 = 135135.0 * Ftmp149;
  double Ftmp151 = Ftmp113 * Ftmp150;
  double Ftmp152 = Ftmp12 * Ftmp98;
  double Ftmp153 = 103950.0 * Ftmp152;
  double Ftmp154 = Ftmp148 + Ftmp151 - Ftmp153;
  double Ftmp155 = Ftmp154 * M[87];
  double Ftmp156 = Ftmp131 * Ftmp9;
  double Ftmp157 = Ftmp136 * M[49];
  double Ftmp158 = Ftmp140 * Ftmp26;
  double Ftmp159 = Ftmp144 - Ftmp158 + Ftmp50;
  double Ftmp160 = Ftmp159 * M[53];
  double Ftmp161 = Ftmp131 * Ftmp62;
  double Ftmp162 = Ftmp142 * M[50];
  double Ftmp163 = Ftmp146 * M[54];
  double Ftmp164 = Ftmp117 * Ftmp150;
  double Ftmp165 = Ftmp23 * Ftmp98;
  double Ftmp166 = 103950.0 * Ftmp165;
  double Ftmp167 = Ftmp148 + Ftmp164 - Ftmp166;
  double Ftmp168 = Ftmp167 * M[105];
  double Ftmp169 = Ftmp121 * Ftmp149;
  double Ftmp170 = 45045.0 * Ftmp169;
  double Ftmp171 = Ftmp26 * Ftmp98;
  double Ftmp172 = Ftmp170 - 34650.0 * Ftmp171 + Ftmp96;
  double Ftmp173 = Ftmp172 * M[109];
  double Ftmp174 = Ftmp167 * M[77];
  double Ftmp175 = 3675.0 * Ftmp17;
  double Ftmp176 = Ftmp142 * M[44];
  double Ftmp177 = Ftmp159 * M[48];
  double Ftmp178 = 33075.0 * Ftmp52;
  double Ftmp179 = 145530.0 * Ftmp152;
  double Ftmp180 = Ftmp151 + Ftmp178 - Ftmp179;
  double Ftmp181 = Ftmp121 * Ftmp150;
  double Ftmp182 = 103950.0 * Ftmp171;
  double Ftmp183 = Ftmp148 + Ftmp181 - Ftmp182;
  double Ftmp184 = Ftmp35 * M[81];
  double Ftmp185 = Ftmp167 * M[70];
  double Ftmp186 = (Ftmp170 - 20790.0 * Ftmp171 + Ftmp53) * M[74];
  double Ftmp187 = 62370.0 * Ftmp98;
  double Ftmp188 = Ftmp187 * Ftmp23;
  double Ftmp189 = (Ftmp104 + Ftmp164 - Ftmp188) * M[71];
  double Ftmp190 = Ftmp183 * M[75];
  double Ftmp191 = 363825.0 * Ftmp98;
  double Ftmp192 = pow(R, -17.0);
  double Ftmp193 = 2027025.0 * Ftmp192;
  double Ftmp194 = Ftmp113 * Ftmp193;
  double Ftmp195 = Ftmp12 * Ftmp149;
  double Ftmp196 = 1891890.0 * Ftmp195;
  double Ftmp197 = 155925.0 * Ftmp98;
  double Ftmp198 = Ftmp117 * Ftmp193;
  double Ftmp199 = Ftmp149 * Ftmp23;
  double Ftmp200 = 1351350.0 * Ftmp199;
  double Ftmp201 = (Ftmp197 + Ftmp198 - Ftmp200) * M[105];
  double Ftmp202 = 51975.0 * Ftmp98;
  double Ftmp203 = 675675.0 * Ftmp121 * Ftmp192;
  double Ftmp204 = Ftmp149 * Ftmp26;
  double Ftmp205 = (Ftmp202 + Ftmp203 - 450450.0 * Ftmp204) * M[109];
  double Ftmp206 = -11025.0 * Ftmp17;
  double Ftmp207 = (x * x * x * x * x * x);
  double Ftmp208 = Ftmp150 * Ftmp207;
  double Ftmp209 = 99225.0 * Ftmp52;
  double Ftmp210 = Ftmp12 * Ftmp209 - 218295.0 * Ftmp127 + Ftmp206 + Ftmp208;
  double Ftmp211 = -Ftmp126;
  double Ftmp212 = (y * y * y * y * y * y);
  double Ftmp213 = Ftmp150 * Ftmp212;
  double Ftmp214 = 42525.0 * Ftmp52;
  double Ftmp215 = -Ftmp117 * Ftmp197 + Ftmp211 + Ftmp213 + Ftmp214 * Ftmp23;
  double Ftmp216 = (z * z * z * z * z * z);
  double Ftmp217 = Ftmp150 * Ftmp216;
  double Ftmp218 = -Ftmp121 * Ftmp197 + Ftmp211 + Ftmp214 * Ftmp26 + Ftmp217;
  double Ftmp219 = Ftmp215 * M[76];
  double Ftmp220 = Ftmp218 * M[82];
  double Ftmp221 = -Ftmp23 * Ftmp50;
  double Ftmp222 = Ftmp23 * Ftmp54;
  double Ftmp223 = -Ftmp19;
  double Ftmp224 = Ftmp223 + Ftmp31;
  double Ftmp225 = Ftmp221 + Ftmp222 + Ftmp224;
  double Ftmp226 = -Ftmp26 * Ftmp50;
  double Ftmp227 = Ftmp26 * Ftmp54;
  double Ftmp228 = Ftmp224 + Ftmp226 + Ftmp227;
  double Ftmp229 = -Ftmp47;
  double Ftmp230 = Ftmp229 + Ftmp7;
  double Ftmp231 = -Ftmp36;
  double Ftmp232 = Ftmp26 * Ftmp67;
  double Ftmp233 = Ftmp231 + Ftmp232;
  double Ftmp234 = Ftmp230 + Ftmp233;
  double Ftmp235 = -Ftmp209;
  double Ftmp236 = Ftmp193 * Ftmp207;
  double Ftmp237 = Ftmp113 * Ftmp149;
  double Ftmp238 = 1091475.0 * Ftmp152 + Ftmp235 + Ftmp236 - 2837835.0 * Ftmp237;
  double Ftmp239 = Ftmp238 * Ftmp9;
  double Ftmp240 = Ftmp193 * Ftmp212;
  double Ftmp241 = Ftmp117 * Ftmp149;
  double Ftmp242 = 1091475.0 * Ftmp165 + Ftmp235 + Ftmp240 - 2837835.0 * Ftmp241;
  double Ftmp243 = Ftmp242 * M[111];
  double Ftmp244 = -Ftmp148;
  double Ftmp245 = Ftmp193 * Ftmp216;
  double Ftmp246 = 2027025.0 * Ftmp149;
  double Ftmp247 = 467775.0 * Ftmp98;
  double Ftmp248 = -Ftmp121 * Ftmp246 + Ftmp244 + Ftmp245 + Ftmp247 * Ftmp26;
  double Ftmp249 = Ftmp248 * M[117];
  double Ftmp250 = Ftmp238 * Ftmp62;
  double Ftmp251 = -Ftmp117 * Ftmp246 + Ftmp23 * Ftmp247 + Ftmp240 + Ftmp244;
  double Ftmp252 = Ftmp251 * M[112];
  double Ftmp253 = -2837835.0 * Ftmp169 + 1091475.0 * Ftmp171 + Ftmp235 + Ftmp245;
  double Ftmp254 = Ftmp253 * M[118];
  double Ftmp255 = -297675.0 * Ftmp52;
  double Ftmp256 = Ftmp251 * M[104];
  double Ftmp257 = Ftmp248 * M[110];
  double Ftmp258 = Ftmp104 * Ftmp23;
  double Ftmp259 = -Ftmp258;
  double Ftmp260 = Ftmp100 * Ftmp23;
  double Ftmp261 = Ftmp259 + Ftmp260;
  double Ftmp262 = 945.0 * Ftmp17;
  double Ftmp263 = Ftmp104 * Ftmp12;
  double Ftmp264 = -Ftmp263;
  double Ftmp265 = Ftmp262 + Ftmp264;
  double Ftmp266 = Ftmp261 + Ftmp265;
  double Ftmp267 = Ftmp266 * M[61];
  double Ftmp268 = -Ftmp54;
  double Ftmp269 = Ftmp268 + Ftmp50;
  double Ftmp270 = Ftmp104 * Ftmp26;
  double Ftmp271 = -Ftmp270;
  double Ftmp272 = Ftmp100 * Ftmp26;
  double Ftmp273 = Ftmp271 + Ftmp272;
  double Ftmp274 = Ftmp269 + Ftmp273;
  double Ftmp275 = Ftmp274 * M[63];
  double Ftmp276 = -Ftmp67;
  double Ftmp277 = Ftmp106 * Ftmp26;
  double Ftmp278 = Ftmp277 + Ftmp50;
  double Ftmp279 = Ftmp271 + Ftmp276 + Ftmp278;
  double Ftmp280 = Ftmp279 * M[72];
  double Ftmp281 = Ftmp261 + Ftmp269;
  double Ftmp282 = Ftmp281 * M[62];
  double Ftmp283 = Ftmp265 + Ftmp273;
  double Ftmp284 = Ftmp283 * M[64];
  double Ftmp285 = -Ftmp86;
  double Ftmp286 = Ftmp259 + Ftmp278 + Ftmp285;
  double Ftmp287 = Ftmp286 * M[73];
  double Ftmp288 = 8505.0 * Ftmp52;
  double Ftmp289 = 31185.0 * Ftmp98;
  double Ftmp290 = Ftmp12 * Ftmp289;
  double Ftmp291 = -Ftmp290;
  double Ftmp292 = Ftmp288 + Ftmp291;
  double Ftmp293 = Ftmp23 * Ftmp289;
  double Ftmp294 = -Ftmp293;
  double Ftmp295 = Ftmp12 * Ftmp150;
  double Ftmp296 = Ftmp23 * Ftmp295;
  double Ftmp297 = Ftmp294 + Ftmp296;
  double Ftmp298 = Ftmp292 + Ftmp297;
  double Ftmp299 = Ftmp298 * M[94];
  double Ftmp300 = Ftmp26 * Ftmp289;
  double Ftmp301 = -Ftmp300;
  double Ftmp302 = Ftmp26 * Ftmp295;
  double Ftmp303 = Ftmp301 + Ftmp302;
  double Ftmp304 = Ftmp292 + Ftmp303;
  double Ftmp305 = Ftmp304 * M[96];
  double Ftmp306 = Ftmp266 * Ftmp9;
  double Ftmp307 = Ftmp274 * Ftmp9;
  double Ftmp308 = Ftmp281 * Ftmp62;
  double Ftmp309 = Ftmp283 * Ftmp62;
  double Ftmp310 = Ftmp23 * Ftmp26;
  double Ftmp311 = Ftmp150 * Ftmp310;
  double Ftmp312 = Ftmp301 + Ftmp311;
  double Ftmp313 = Ftmp288 + Ftmp294 + Ftmp312;
  double Ftmp314 = Ftmp313 * M[107];
  double Ftmp315 = Ftmp14 * Ftmp298;
  double Ftmp316 = Ftmp14 * Ftmp304;
  double Ftmp317 = -Ftmp23 * Ftmp96;
  double Ftmp318 = Ftmp268 + Ftmp81;
  double Ftmp319 = -Ftmp26 * Ftmp96;
  double Ftmp320 = Ftmp18 + Ftmp285;
  double Ftmp321 = Ftmp276 + Ftmp277;
  double Ftmp322 = Ftmp202 * Ftmp23;
  double Ftmp323 = -Ftmp322;
  double Ftmp324 = Ftmp296 + Ftmp323;
  double Ftmp325 = Ftmp148 + Ftmp291;
  double Ftmp326 = -Ftmp100;
  double Ftmp327 = Ftmp326 + Ftmp96;
  double Ftmp328 = Ftmp202 * Ftmp26;
  double Ftmp329 = -Ftmp328;
  double Ftmp330 = Ftmp302 + Ftmp329;
  double Ftmp331 = -Ftmp106;
  double Ftmp332 = Ftmp104 + Ftmp331;
  double Ftmp333 = Ftmp26 * Ftmp99;
  double Ftmp334 = -Ftmp333;
  double Ftmp335 = Ftmp104 + Ftmp334;
  double Ftmp336 = Ftmp294 + Ftmp311;
  double Ftmp337 = 675675.0 * Ftmp149;
  double Ftmp338 = Ftmp23 * Ftmp337;
  double Ftmp339 = -Ftmp338;
  double Ftmp340 = Ftmp12 * Ftmp193;
  double Ftmp341 = Ftmp23 * Ftmp340;
  double Ftmp342 = 405405.0 * Ftmp195;
  double Ftmp343 = -Ftmp342;
  double Ftmp344 = Ftmp197 + Ftmp343;
  double Ftmp345 = Ftmp26 * Ftmp337;
  double Ftmp346 = -Ftmp345;
  double Ftmp347 = Ftmp26 * Ftmp340;
  double Ftmp348 = 93555.0 * Ftmp98;
  double Ftmp349 = -405405.0 * Ftmp204;
  double Ftmp350 = Ftmp348 + Ftmp349;
  double Ftmp351 = 405405.0 * Ftmp199;
  double Ftmp352 = -Ftmp351;
  double Ftmp353 = Ftmp193 * Ftmp310;
  double Ftmp354 = Ftmp352 + Ftmp353;
  double Ftmp355 = -Ftmp262;
  double Ftmp356 = Ftmp263 + Ftmp355;
  double Ftmp357 = Ftmp12 * Ftmp164;
  double Ftmp358 = 62370.0 * Ftmp152;
  double Ftmp359 = -Ftmp23 * Ftmp358;
  double Ftmp360 = Ftmp357 + Ftmp359;
  double Ftmp361 = 17010.0 * Ftmp52;
  double Ftmp362 = -Ftmp117 * Ftmp289 + Ftmp23 * Ftmp361;
  double Ftmp363 = Ftmp356 + Ftmp360 + Ftmp362;
  double Ftmp364 = Ftmp12 * Ftmp181;
  double Ftmp365 = -Ftmp26 * Ftmp358;
  double Ftmp366 = Ftmp364 + Ftmp365;
  double Ftmp367 = -Ftmp121 * Ftmp289 + Ftmp26 * Ftmp361;
  double Ftmp368 = Ftmp356 + Ftmp366 + Ftmp367;
  double Ftmp369 = Ftmp148 * Ftmp23;
  double Ftmp370 = -Ftmp153 * Ftmp23 + Ftmp211;
  double Ftmp371 = -Ftmp128;
  double Ftmp372 = Ftmp151 * Ftmp23;
  double Ftmp373 = Ftmp371 + Ftmp372;
  double Ftmp374 = Ftmp130 + Ftmp369 + Ftmp370 + Ftmp373;
  double Ftmp375 = -Ftmp153 * Ftmp26;
  double Ftmp376 = Ftmp151 * Ftmp26;
  double Ftmp377 = Ftmp371 + Ftmp376;
  double Ftmp378 = Ftmp148 * Ftmp26 + Ftmp211;
  double Ftmp379 = Ftmp130 + Ftmp375 + Ftmp377 + Ftmp378;
  double Ftmp380 = Ftmp187 * Ftmp26;
  double Ftmp381 = -Ftmp23 * Ftmp380;
  double Ftmp382 = Ftmp381 + Ftmp51;
  double Ftmp383 = -Ftmp144;
  double Ftmp384 = Ftmp181 * Ftmp23;
  double Ftmp385 = Ftmp383 + Ftmp384;
  double Ftmp386 = Ftmp158 + Ftmp258 + Ftmp382 + Ftmp385;
  double Ftmp387 = -Ftmp133;
  double Ftmp388 = Ftmp164 * Ftmp26;
  double Ftmp389 = Ftmp387 + Ftmp388;
  double Ftmp390 = Ftmp141 + Ftmp270 + Ftmp382 + Ftmp389;
  double Ftmp391 = Ftmp12 * Ftmp197;
  double Ftmp392 = -405405.0 * Ftmp241;
  double Ftmp393 = 311850.0 * Ftmp98;
  double Ftmp394 = Ftmp23 * Ftmp393;
  double Ftmp395 = Ftmp12 * Ftmp198;
  double Ftmp396 = Ftmp394 + Ftmp395;
  double Ftmp397 = -Ftmp214;
  double Ftmp398 = 1351350.0 * Ftmp195;
  double Ftmp399 = -Ftmp23 * Ftmp398;
  double Ftmp400 = Ftmp397 + Ftmp399;
  double Ftmp401 = Ftmp9 * (Ftmp391 + Ftmp392 + Ftmp396 + Ftmp400);
  double Ftmp402 = Ftmp121 * Ftmp193;
  double Ftmp403 = Ftmp12 * Ftmp402;
  double Ftmp404 = 810810.0 * Ftmp195;
  double Ftmp405 = -Ftmp26 * Ftmp404;
  double Ftmp406 = Ftmp403 + Ftmp405;
  double Ftmp407 = -Ftmp288;
  double Ftmp408 = Ftmp290 + Ftmp407;
  double Ftmp409 = -405405.0 * Ftmp169;
  double Ftmp410 = 187110.0 * Ftmp171 + Ftmp409;
  double Ftmp411 = Ftmp9 * (Ftmp406 + Ftmp408 + Ftmp410);
  double Ftmp412 = Ftmp194 * Ftmp23;
  double Ftmp413 = Ftmp197 * Ftmp23;
  double Ftmp414 = 311850.0 * Ftmp152;
  double Ftmp415 = -405405.0 * Ftmp237;
  double Ftmp416 = Ftmp414 + Ftmp415;
  double Ftmp417 = Ftmp9 * (Ftmp400 + Ftmp412 + Ftmp413 + Ftmp416);
  double Ftmp418 = -Ftmp26 * Ftmp398;
  double Ftmp419 = -Ftmp151;
  double Ftmp420 = Ftmp194 * Ftmp26;
  double Ftmp421 = Ftmp419 + Ftmp420;
  double Ftmp422 = Ftmp197 * Ftmp26;
  double Ftmp423 = Ftmp244 + Ftmp422;
  double Ftmp424 = Ftmp9 * (Ftmp153 + Ftmp418 + Ftmp421 + Ftmp423);
  double Ftmp425 = -810810.0 * Ftmp199 * Ftmp26;
  double Ftmp426 = Ftmp407 + Ftmp425;
  double Ftmp427 = Ftmp23 * Ftmp402;
  double Ftmp428 = Ftmp293 + Ftmp427;
  double Ftmp429 = Ftmp410 + Ftmp426 + Ftmp428;
  double Ftmp430 = -Ftmp200 * Ftmp26;
  double Ftmp431 = -Ftmp164;
  double Ftmp432 = Ftmp198 * Ftmp26;
  double Ftmp433 = Ftmp431 + Ftmp432;
  double Ftmp434 = Ftmp166 + Ftmp423 + Ftmp430 + Ftmp433;
  double Ftmp435 = 187110.0 * Ftmp165 + Ftmp392;
  double Ftmp436 = -Ftmp23 * Ftmp404;
  double Ftmp437 = Ftmp395 + Ftmp436;
  double Ftmp438 = Ftmp62 * (Ftmp408 + Ftmp435 + Ftmp437);
  double Ftmp439 = Ftmp397 + Ftmp418;
  double Ftmp440 = Ftmp391 + Ftmp403;
  double Ftmp441 = Ftmp26 * Ftmp393;
  double Ftmp442 = Ftmp409 + Ftmp441;
  double Ftmp443 = Ftmp62 * (Ftmp439 + Ftmp440 + Ftmp442);
  double Ftmp444 = Ftmp412 + Ftmp419;
  double Ftmp445 = Ftmp244 + Ftmp413;
  double Ftmp446 = Ftmp62 * (Ftmp153 + Ftmp399 + Ftmp444 + Ftmp445);
  double Ftmp447 = Ftmp62 * (Ftmp416 + Ftmp420 + Ftmp422 + Ftmp439);
  double Ftmp448 = -Ftmp181;
  double Ftmp449 = Ftmp427 + Ftmp448;
  double Ftmp450 = Ftmp182 + Ftmp430 + Ftmp445 + Ftmp449;
  double Ftmp451 = Ftmp300 + Ftmp432;
  double Ftmp452 = Ftmp426 + Ftmp435 + Ftmp451;
  double Ftmp453 = -Ftmp117 * Ftmp337;
  double Ftmp454 = Ftmp244 + Ftmp290;
  double Ftmp455 = -Ftmp121 * Ftmp337 + Ftmp441;
  double Ftmp456 = Ftmp191 * Ftmp23;
  double Ftmp457 = -Ftmp178;
  double Ftmp458 = -Ftmp196 * Ftmp23 + Ftmp457;
  double Ftmp459 = -Ftmp196 * Ftmp26;
  double Ftmp460 = Ftmp191 * Ftmp26 + Ftmp457;
  double Ftmp461 = Ftmp105 + Ftmp425;
  double Ftmp462 = Ftmp26 * Ftmp296;
  double Ftmp463 = -Ftmp260 + Ftmp270 + Ftmp462;
  double Ftmp464 = Ftmp258 - Ftmp272;
  double Ftmp465 = -Ftmp26 * Ftmp293 + Ftmp463 + Ftmp464 + Ftmp55;
  double Ftmp466 = Ftmp310 * Ftmp340;
  double Ftmp467 = -Ftmp296 + Ftmp466;
  double Ftmp468 = -Ftmp26 * Ftmp351 + Ftmp408;
  double Ftmp469 = -Ftmp26 * Ftmp342 + Ftmp293;
  double Ftmp470 = Ftmp9 * (Ftmp26 * Ftmp348 + Ftmp467 + Ftmp468 + Ftmp469);
  double Ftmp471 = -Ftmp302;
  double Ftmp472 = -Ftmp23 * Ftmp342 + Ftmp300 + Ftmp466;
  double Ftmp473 = Ftmp62 * (Ftmp23 * Ftmp348 + Ftmp468 + Ftmp471 + Ftmp472);
  double Ftmp474 = Ftmp328 + Ftmp467;
  double Ftmp475 = Ftmp322 + Ftmp471;
  double Ftmp476 = x * M[5];
  double Ftmp477 = Ftmp16 + Ftmp22;
  double Ftmp478 = Ftmp20 + Ftmp24;
  double Ftmp479 = Ftmp477 * M[3];
  double Ftmp480 = 1.0 * x;
  double Ftmp481 = 3.0 * x;
  double Ftmp482 = Ftmp19 + Ftmp44;
  double Ftmp483 = Ftmp482 * M[23];
  double Ftmp484 = Ftmp37 * M[30];
  double Ftmp485 = Ftmp43 * x;
  double Ftmp486 = 3.0 * Ftmp62;
  double Ftmp487 = Ftmp482 * M[11];
  double Ftmp488 = Ftmp55 * M[21];
  double Ftmp489 = Ftmp482 * M[10];
  double Ftmp490 = 1.0 * Ftmp23;
  double Ftmp491 = Ftmp23 * x;
  double Ftmp492 = Ftmp55 * M[20];
  double Ftmp493 = Ftmp23 * z;
  double Ftmp494 = (Ftmp54 + Ftmp71) * M[23];
  double Ftmp495 = Ftmp67 + Ftmp82;
  double Ftmp496 = Ftmp35 * Ftmp62;
  double Ftmp497 = Ftmp23 * Ftmp480;
  double Ftmp498 = Ftmp23 * Ftmp481;
  double Ftmp499 = Ftmp23 * Ftmp43;
  double Ftmp500 = Ftmp23 * Ftmp62;
  double Ftmp501 = (Ftmp100 + Ftmp105) * M[38];
  double Ftmp502 = Ftmp106 + Ftmp97;
  double Ftmp503 = Ftmp43 * Ftmp491;
  double Ftmp504 = Ftmp23 * Ftmp486;
  double Ftmp505 = Ftmp114 - Ftmp119 * Ftmp12 + Ftmp31;
  double Ftmp506 = Ftmp112 - Ftmp115 * Ftmp23 + Ftmp118;
  double Ftmp507 = Ftmp505 * M[19];
  double Ftmp508 = Ftmp12 * Ftmp140;
  double Ftmp509 = Ftmp128 + Ftmp50 - Ftmp508;
  double Ftmp510 = Ftmp509 * M[59];
  double Ftmp511 = Ftmp136 * M[77];
  double Ftmp512 = Ftmp509 * M[36];
  double Ftmp513 = Ftmp154 * M[57];
  double Ftmp514 = Ftmp509 * M[35];
  double Ftmp515 = Ftmp154 * M[56];
  double Ftmp516 = (Ftmp104 + Ftmp151 - Ftmp358) * M[59];
  double Ftmp517 = 145530.0 * Ftmp165;
  double Ftmp518 = Ftmp164 + Ftmp178 - Ftmp517;
  double Ftmp519 = (Ftmp194 + Ftmp197 - Ftmp398) * M[87];
  double Ftmp520 = 1891890.0 * Ftmp199;
  double Ftmp521 = -Ftmp113 * Ftmp197 + Ftmp12 * Ftmp214 + Ftmp208 + Ftmp211;
  double Ftmp522 = 218295.0 * Ftmp98;
  double Ftmp523 = -Ftmp117 * Ftmp522 + Ftmp206 + Ftmp209 * Ftmp23 + Ftmp213;
  double Ftmp524 = Ftmp521 * M[55];
  double Ftmp525 = Ftmp222 + Ftmp231;
  double Ftmp526 = -Ftmp12 * Ftmp50 + Ftmp31;
  double Ftmp527 = Ftmp525 + Ftmp526;
  double Ftmp528 = Ftmp223 + Ftmp227 + Ftmp230;
  double Ftmp529 = Ftmp226 + Ftmp233 + Ftmp31;
  double Ftmp530 = 467775.0 * Ftmp152 + Ftmp236 - 2027025.0 * Ftmp237 + Ftmp244;
  double Ftmp531 = Ftmp530 * M[85];
  double Ftmp532 = Ftmp530 * M[84];
  double Ftmp533 = Ftmp264 + Ftmp50;
  double Ftmp534 = Ftmp260 + Ftmp276;
  double Ftmp535 = Ftmp533 + Ftmp534;
  double Ftmp536 = Ftmp535 * M[66];
  double Ftmp537 = Ftmp272 + Ftmp285;
  double Ftmp538 = Ftmp533 + Ftmp537;
  double Ftmp539 = Ftmp538 * M[68];
  double Ftmp540 = Ftmp259 + Ftmp262 + Ftmp271 + Ftmp277;
  double Ftmp541 = Ftmp540 * M[79];
  double Ftmp542 = Ftmp5 * Ftmp535;
  double Ftmp543 = Ftmp5 * Ftmp538;
  double Ftmp544 = Ftmp5 * Ftmp540;
  double Ftmp545 = -Ftmp12 * Ftmp96 + Ftmp81;
  double Ftmp546 = Ftmp12 * Ftmp202;
  double Ftmp547 = -Ftmp546;
  double Ftmp548 = Ftmp148 + Ftmp547;
  double Ftmp549 = Ftmp104 + Ftmp326;
  double Ftmp550 = Ftmp331 + Ftmp96;
  double Ftmp551 = Ftmp313 * Ftmp496;
  double Ftmp552 = Ftmp341 + Ftmp352;
  double Ftmp553 = -Ftmp12 * Ftmp337 + Ftmp197;
  double Ftmp554 = Ftmp12 * Ftmp148;
  double Ftmp555 = Ftmp135 + Ftmp357 + Ftmp370 + Ftmp387 + Ftmp554;
  double Ftmp556 = Ftmp263 + Ftmp51;
  double Ftmp557 = Ftmp158 + Ftmp366 + Ftmp383 + Ftmp556;
  double Ftmp558 = Ftmp258 + Ftmp355;
  double Ftmp559 = Ftmp12 * Ftmp361 - 31185.0 * Ftmp127;
  double Ftmp560 = Ftmp359 + Ftmp372 + Ftmp558 + Ftmp559;
  double Ftmp561 = Ftmp508 + Ftmp51;
  double Ftmp562 = Ftmp270 + Ftmp365;
  double Ftmp563 = Ftmp377 + Ftmp561 + Ftmp562;
  double Ftmp564 = Ftmp367 + Ftmp381 + Ftmp384 + Ftmp558;
  double Ftmp565 = -Ftmp166 * Ftmp26;
  double Ftmp566 = Ftmp135 + Ftmp378 + Ftmp389 + Ftmp565;
  double Ftmp567 = Ftmp395 + Ftmp431;
  double Ftmp568 = Ftmp5 * (Ftmp166 + Ftmp244 + Ftmp391 + Ftmp399 + Ftmp567);
  double Ftmp569 = Ftmp5 * (Ftmp182 + Ftmp244 + Ftmp418 + Ftmp440 + Ftmp448);
  double Ftmp570 = Ftmp293 + Ftmp436;
  double Ftmp571 = Ftmp412 + Ftmp570;
  double Ftmp572 = 187110.0 * Ftmp152 + Ftmp407 + Ftmp415;
  double Ftmp573 = Ftmp5 * (Ftmp571 + Ftmp572);
  double Ftmp574 = Ftmp300 + Ftmp405;
  double Ftmp575 = Ftmp420 + Ftmp574;
  double Ftmp576 = Ftmp5 * (Ftmp572 + Ftmp575);
  double Ftmp577 = Ftmp397 + Ftmp430;
  double Ftmp578 = Ftmp5 * (Ftmp413 + Ftmp427 + Ftmp442 + Ftmp577);
  double Ftmp579 = Ftmp5 * (Ftmp392 + Ftmp394 + Ftmp422 + Ftmp432 + Ftmp577);
  double Ftmp580 = Ftmp12 * Ftmp191;
  double Ftmp581 = Ftmp105 + Ftmp290;
  double Ftmp582 = -675675.0 * Ftmp237 + Ftmp244 + Ftmp414;
  double Ftmp583 = Ftmp105 + Ftmp358;
  double Ftmp584 = Ftmp244 + Ftmp425;
  double Ftmp585 = -Ftmp26 * Ftmp520;
  double Ftmp586 = Ftmp263 - Ftmp277;
  double Ftmp587 = -Ftmp26 * Ftmp290 + Ftmp463 + Ftmp586 + Ftmp68;
  double Ftmp588 = -Ftmp311;
  double Ftmp589 = Ftmp5 * (Ftmp12 * Ftmp348 + Ftmp407 + Ftmp469 + Ftmp472 + Ftmp588);
  double Ftmp590 = Ftmp546 + Ftmp588;
  double Ftmp591 = y * M[7];
  double Ftmp592 = Ftmp20 + Ftmp27;
  double Ftmp593 = Ftmp480 * M[28];
  double Ftmp594 = Ftmp39 * x;
  double Ftmp595 = Ftmp26 * x;
  double Ftmp596 = Ftmp26 * y;
  double Ftmp597 = Ftmp39 * Ftmp62;
  double Ftmp598 = Ftmp26 * Ftmp480;
  double Ftmp599 = Ftmp26 * (Ftmp82 + Ftmp86);
  double Ftmp600 = Ftmp26 * Ftmp9;
  double Ftmp601 = Ftmp35 * Ftmp595;
  double Ftmp602 = Ftmp39 * Ftmp595;
  double Ftmp603 = Ftmp112 - Ftmp115 * Ftmp26 + Ftmp122;
  double Ftmp604 = Ftmp480 * M[75];
  double Ftmp605 = 145530.0 * Ftmp171;
  double Ftmp606 = Ftmp26 * (Ftmp178 + Ftmp181 - Ftmp605);
  double Ftmp607 = -Ftmp121 * Ftmp522 + Ftmp206 + Ftmp209 * Ftmp26 + Ftmp217;
  double Ftmp608 = Ftmp223 + Ftmp525 + Ftmp7;
  double Ftmp609 = Ftmp227 + Ftmp229 + Ftmp526;
  double Ftmp610 = Ftmp221 + Ftmp229 + Ftmp232 + Ftmp31;
  double Ftmp611 = Ftmp334 + Ftmp96;
  double Ftmp612 = Ftmp141 + Ftmp360 + Ftmp387 + Ftmp556;
  double Ftmp613 = Ftmp145 + Ftmp211;
  double Ftmp614 = Ftmp364 + Ftmp375 + Ftmp383 + Ftmp554 + Ftmp613;
  double Ftmp615 = Ftmp258 + Ftmp359 + Ftmp373 + Ftmp561;
  double Ftmp616 = Ftmp355 + Ftmp376 + Ftmp559 + Ftmp562;
  double Ftmp617 = Ftmp369 + Ftmp385 + Ftmp565 + Ftmp613;
  double Ftmp618 = Ftmp270 + Ftmp355 + Ftmp362 + Ftmp381 + Ftmp388;
  double Ftmp619 = Ftmp457 + Ftmp605;
  double Ftmp620 = -Ftmp23 * Ftmp290 + Ftmp462 + Ftmp464 + Ftmp586 + Ftmp87;
#pragma omp atomic
  F[0] +=
        Ftmp0 * M[0] - Ftmp10 * Ftmp9 - Ftmp101 * Ftmp102 * M[38] -
        Ftmp102 * (Ftmp191 + Ftmp194 - Ftmp196) * M[87] -
        Ftmp102 * (Ftmp339 + Ftmp341 + Ftmp344) * M[94] -
        Ftmp102 * (Ftmp344 + Ftmp346 + Ftmp347) * M[96] - Ftmp103 * Ftmp107 -
        Ftmp103 * Ftmp201 - Ftmp103 * (Ftmp350 + Ftmp354) * M[107] - Ftmp108 * Ftmp111 -
        Ftmp108 * Ftmp205 - Ftmp11 * x + Ftmp116 * x * M[19] + Ftmp116 * M[34] -
        Ftmp12 * Ftmp13 - Ftmp12 * (Ftmp19 + Ftmp77) * M[9] -
        Ftmp12 * (Ftmp128 - 13230.0 * Ftmp129 + Ftmp175) * M[34] -
        Ftmp12 * (Ftmp260 + Ftmp317 + Ftmp318) * M[37] -
        Ftmp12 * (Ftmp272 + Ftmp318 + Ftmp319) * M[39] -
        Ftmp12 * (Ftmp406 + Ftmp454 + Ftmp455) * M[97] -
        Ftmp12 * (Ftmp101 - Ftmp26 * Ftmp338 + Ftmp474 + Ftmp475) * M[95] -
        Ftmp12 * (1964655.0 * Ftmp152 + Ftmp236 - 3648645.0 * Ftmp237 + Ftmp255) * M[83] -
        Ftmp12 * (Ftmp179 + Ftmp421 + Ftmp459 + Ftmp460) * M[88] -
        Ftmp12 * (Ftmp179 + Ftmp444 + Ftmp456 + Ftmp458) * M[86] -
        Ftmp12 * (Ftmp396 + Ftmp436 + Ftmp453 + Ftmp454) * M[93] + Ftmp120 * M[44] +
        Ftmp123 * M[48] + Ftmp124 * x + Ftmp125 * x - Ftmp132 * y - Ftmp137 * Ftmp35 -
        Ftmp138 * Ftmp39 - Ftmp139 * z + Ftmp14 * Ftmp154 * M[59] + Ftmp14 * Ftmp174 +
        Ftmp14 * Ftmp313 * M[79] + Ftmp14 * Ftmp55 * M[23] + Ftmp14 * Ftmp76 -
        Ftmp143 * Ftmp43 - Ftmp147 * M[75] + Ftmp15 * M[7] + Ftmp155 * Ftmp5 -
        Ftmp156 * M[35] - Ftmp157 * Ftmp9 + Ftmp16 * y * M[4] + Ftmp16 * z * M[5] -
        Ftmp160 * Ftmp61 - Ftmp161 * M[36] - Ftmp162 * Ftmp62 - Ftmp163 * Ftmp62 +
        Ftmp168 * Ftmp66 + Ftmp173 * Ftmp70 - Ftmp176 * Ftmp78 - Ftmp177 * Ftmp78 +
        Ftmp180 * Ftmp84 * M[56] + Ftmp180 * Ftmp85 * M[57] + Ftmp183 * Ftmp184 * Ftmp62 +
        Ftmp185 * Ftmp89 + Ftmp186 * Ftmp91 + Ftmp189 * Ftmp93 - Ftmp19 * Ftmp5 * M[13] +
        Ftmp190 * Ftmp93 + Ftmp21 * x * M[3] + Ftmp21 * M[9] + Ftmp210 * x * M[55] +
        Ftmp210 * M[83] + Ftmp215 * M[104] + Ftmp218 * M[110] + Ftmp219 * x +
        Ftmp220 * x + Ftmp225 * x * M[22] + Ftmp225 * M[37] + Ftmp228 * x * M[24] +
        Ftmp228 * M[39] + Ftmp234 * x * M[31] + Ftmp234 * M[46] - Ftmp239 * M[84] -
        Ftmp243 * Ftmp9 - Ftmp249 * Ftmp61 + Ftmp25 * M[12] - Ftmp250 * M[85] -
        Ftmp252 * Ftmp62 - Ftmp254 * Ftmp62 - Ftmp256 * Ftmp78 - Ftmp257 * Ftmp78 -
        Ftmp267 * y - Ftmp275 * y - Ftmp279 * Ftmp9 * M[51] + Ftmp28 * M[14] -
        Ftmp280 * Ftmp35 - Ftmp282 * z - Ftmp284 * z - Ftmp286 * Ftmp62 * M[52] -
        Ftmp287 * Ftmp43 + Ftmp29 * x + Ftmp299 * Ftmp5 - Ftmp3 * y + Ftmp30 * x +
        Ftmp305 * Ftmp5 - Ftmp306 * M[40] - Ftmp307 * M[42] - Ftmp308 * M[41] -
        Ftmp309 * M[43] + Ftmp314 * Ftmp66 + Ftmp315 * M[66] + Ftmp316 * M[68] -
        Ftmp34 * y - Ftmp35 * Ftmp38 + Ftmp363 * x * M[65] + Ftmp363 * M[93] +
        Ftmp368 * x * M[69] + Ftmp368 * M[97] + Ftmp374 * x * M[58] + Ftmp374 * M[86] +
        Ftmp379 * x * M[60] + Ftmp379 * M[88] + Ftmp386 * x * M[80] + Ftmp386 * M[108] -
        Ftmp39 * Ftmp41 + Ftmp390 * x * M[78] + Ftmp390 * M[106] - Ftmp4 * M[5] -
        Ftmp401 * M[98] - Ftmp411 * M[102] - Ftmp417 * M[89] - Ftmp42 * z -
        Ftmp424 * M[91] - Ftmp429 * Ftmp9 * M[115] - Ftmp43 * Ftmp46 -
        Ftmp434 * Ftmp9 * M[113] - Ftmp438 * M[99] - Ftmp443 * M[103] - Ftmp446 * M[90] -
        Ftmp447 * M[92] - Ftmp450 * Ftmp62 * M[116] - Ftmp452 * Ftmp62 * M[114] +
        Ftmp465 * x * M[67] + Ftmp465 * M[95] - Ftmp470 * M[100] - Ftmp473 * M[101] -
        Ftmp49 * M[28] + Ftmp5 * Ftmp56 + Ftmp5 * Ftmp8 - Ftmp57 * M[10] -
        Ftmp58 * Ftmp9 - Ftmp60 * Ftmp61 - Ftmp62 * Ftmp64 - Ftmp62 * Ftmp65 +
        Ftmp62 * Ftmp87 * Ftmp88 - Ftmp63 * M[11] + Ftmp66 * Ftmp69 + Ftmp70 * Ftmp75 -
        Ftmp78 * Ftmp79 - Ftmp78 * Ftmp80 - Ftmp78 * (Ftmp320 + Ftmp321) * M[46] -
        Ftmp78 * (Ftmp188 + Ftmp300 + Ftmp433 + Ftmp461) * M[106] -
        Ftmp78 * (Ftmp293 + Ftmp380 + Ftmp449 + Ftmp461) * M[108] +
        Ftmp83 * Ftmp84 * M[20] + Ftmp83 * Ftmp85 * M[21] +
        Ftmp84 * (Ftmp324 + Ftmp325) * M[61] + Ftmp84 * (Ftmp327 + Ftmp330) * M[63] +
        Ftmp85 * (Ftmp324 + Ftmp327) * M[62] + Ftmp85 * (Ftmp325 + Ftmp330) * M[64] +
        Ftmp89 * Ftmp90 + Ftmp89 * (Ftmp312 + Ftmp332) * M[72] + Ftmp91 * Ftmp92 +
        Ftmp93 * Ftmp94 + Ftmp93 * Ftmp95 + Ftmp93 * (Ftmp335 + Ftmp336) * M[73];
#pragma omp atomic
  F[1] += Ftmp0 * M[1] - Ftmp10 * Ftmp23 - Ftmp11 * y - Ftmp111 * Ftmp504 +
          Ftmp123 * M[53] + Ftmp125 * y - Ftmp13 * Ftmp9 - Ftmp132 * x -
          Ftmp136 * Ftmp5 * M[50] - Ftmp136 * Ftmp61 * M[44] - Ftmp137 * Ftmp480 -
          Ftmp138 * Ftmp481 + Ftmp14 * Ftmp488 + Ftmp14 * Ftmp513 - Ftmp147 * M[81] +
          Ftmp155 * Ftmp62 - Ftmp156 * M[34] - Ftmp160 * Ftmp490 - Ftmp163 * Ftmp5 +
          Ftmp167 * Ftmp496 * M[71] + Ftmp168 * Ftmp485 + Ftmp173 * Ftmp486 -
          Ftmp177 * Ftmp61 + Ftmp183 * Ftmp499 * M[81] + Ftmp186 * Ftmp498 +
          Ftmp190 * Ftmp496 - Ftmp205 * Ftmp504 + Ftmp218 * M[117] + Ftmp220 * y -
          Ftmp23 * Ftmp489 - Ftmp23 * Ftmp514 - Ftmp23 * Ftmp532 -
          Ftmp23 * (Ftmp36 + Ftmp77) * M[15] - Ftmp23 * (Ftmp534 + Ftmp545) * M[40] -
          Ftmp23 * (Ftmp571 + Ftmp582) * M[89] -
          Ftmp23 * (Ftmp133 + Ftmp175 - 13230.0 * Ftmp23 * Ftmp52) * M[49] -
          Ftmp23 * (Ftmp268 + Ftmp272 + Ftmp320) * M[42] -
          Ftmp23 * (Ftmp319 + Ftmp321 + Ftmp81) * M[51] -
          Ftmp23 * (Ftmp421 + Ftmp574 + Ftmp583) * M[91] -
          Ftmp23 * (Ftmp428 + Ftmp455 + Ftmp584) * M[115] -
          Ftmp23 * (1964655.0 * Ftmp165 + Ftmp240 - 3648645.0 * Ftmp241 + Ftmp255) *
                M[111] -
          Ftmp23 * (Ftmp380 + Ftmp406 + Ftmp448 + Ftmp581) * M[102] -
          Ftmp23 * (Ftmp433 + Ftmp460 + Ftmp517 + Ftmp585) * M[113] -
          Ftmp23 * (Ftmp458 + Ftmp517 + Ftmp567 + Ftmp580) * M[98] -
          Ftmp23 * (-Ftmp12 * Ftmp345 + Ftmp474 + Ftmp502 + Ftmp590) * M[100] -
          Ftmp239 * M[83] + Ftmp24 * x * M[4] + Ftmp24 * z * M[7] -
          Ftmp242 * Ftmp5 * M[112] - Ftmp242 * Ftmp61 * M[104] - Ftmp249 * Ftmp490 -
          Ftmp254 * Ftmp5 - Ftmp257 * Ftmp61 - Ftmp267 * x - Ftmp275 * x -
          Ftmp279 * Ftmp61 * M[46] + Ftmp28 * M[17] - Ftmp280 * Ftmp480 +
          Ftmp299 * Ftmp62 - Ftmp3 * x + Ftmp30 * y + Ftmp305 * Ftmp62 - Ftmp306 * M[37] -
          Ftmp307 * M[39] + Ftmp314 * Ftmp485 + Ftmp315 * M[62] + Ftmp316 * M[64] -
          Ftmp34 * x - Ftmp36 * Ftmp62 * M[13] - Ftmp37 * Ftmp5 * M[16] -
          Ftmp37 * Ftmp61 * M[12] - Ftmp38 * Ftmp480 - Ftmp4 * M[7] - Ftmp401 * M[93] -
          Ftmp41 * Ftmp481 - Ftmp411 * M[97] - Ftmp417 * M[86] - Ftmp424 * M[88] -
          Ftmp429 * Ftmp61 * M[108] - Ftmp434 * Ftmp61 * M[106] - Ftmp470 * M[95] +
          Ftmp476 * Ftmp5 * Ftmp7 + Ftmp477 * M[10] + Ftmp478 * y * M[6] +
          Ftmp478 * M[15] + Ftmp479 * y - Ftmp483 * z - Ftmp484 * z + Ftmp485 * Ftmp69 +
          Ftmp486 * Ftmp75 - Ftmp487 * Ftmp5 - Ftmp49 * M[32] - Ftmp490 * Ftmp60 +
          Ftmp491 * Ftmp492 + Ftmp491 * Ftmp515 + Ftmp491 * (Ftmp297 + Ftmp548) * M[61] +
          Ftmp491 * (Ftmp303 + Ftmp549) * M[63] + Ftmp493 * Ftmp494 +
          Ftmp493 * Ftmp495 * M[30] + Ftmp493 * Ftmp516 + Ftmp493 * Ftmp518 * M[77] +
          Ftmp493 * (Ftmp148 + Ftmp329 + Ftmp336) * M[79] +
          Ftmp493 * (Ftmp291 + Ftmp302 + Ftmp335) * M[68] +
          Ftmp493 * (Ftmp296 + Ftmp547 + Ftmp550) * M[66] + Ftmp495 * Ftmp497 * M[25] +
          Ftmp496 * Ftmp68 * M[26] + Ftmp496 * Ftmp95 + Ftmp497 * Ftmp518 * M[70] +
          Ftmp497 * (Ftmp311 + Ftmp329 + Ftmp550) * M[72] + Ftmp498 * Ftmp92 +
          Ftmp499 * Ftmp87 * M[32] - Ftmp5 * Ftmp512 - Ftmp5 * Ftmp531 - Ftmp5 * Ftmp65 -
          Ftmp500 * Ftmp501 - Ftmp500 * Ftmp519 - Ftmp500 * (Ftmp552 + Ftmp553) * M[94] -
          Ftmp500 * (Ftmp343 + Ftmp347 + Ftmp350) * M[96] - Ftmp502 * Ftmp503 * M[45] -
          Ftmp503 * (Ftmp191 + Ftmp198 - Ftmp520) * M[105] -
          Ftmp503 * (Ftmp197 + Ftmp346 + Ftmp354) * M[107] + Ftmp505 * M[35] +
          Ftmp506 * y * M[29] + Ftmp506 * M[49] + Ftmp507 * y - Ftmp510 * z -
          Ftmp511 * z + Ftmp521 * M[84] + Ftmp523 * y * M[76] + Ftmp523 * M[111] +
          Ftmp524 * y + Ftmp527 * y * M[22] + Ftmp527 * M[40] + Ftmp528 * y * M[24] +
          Ftmp528 * M[42] + Ftmp529 * y * M[31] + Ftmp529 * M[51] - Ftmp536 * z -
          Ftmp539 * z - Ftmp541 * z - Ftmp542 * M[41] - Ftmp543 * M[43] -
          Ftmp544 * M[52] + Ftmp551 * M[73] + Ftmp555 * y * M[65] + Ftmp555 * M[98] +
          Ftmp557 * y * M[69] + Ftmp557 * M[102] + Ftmp56 * Ftmp62 + Ftmp560 * y * M[58] +
          Ftmp560 * M[89] + Ftmp563 * y * M[60] + Ftmp563 * M[91] + Ftmp564 * y * M[80] +
          Ftmp564 * M[115] + Ftmp566 * y * M[78] + Ftmp566 * M[113] - Ftmp568 * M[99] -
          Ftmp569 * M[103] - Ftmp57 * M[9] - Ftmp573 * M[90] - Ftmp576 * M[92] -
          Ftmp578 * M[116] - Ftmp579 * M[114] + Ftmp587 * y * M[67] + Ftmp587 * M[100] -
          Ftmp589 * M[101] - Ftmp61 * Ftmp80 + Ftmp62 * Ftmp8;
#pragma omp atomic
  F[2] += Ftmp0 * M[2] - Ftmp107 * Ftmp601 + Ftmp120 * M[50] + Ftmp124 * z - Ftmp139 * x +
          Ftmp14 * Ftmp492 + Ftmp14 * Ftmp515 - Ftmp143 * Ftmp480 - Ftmp146 * Ftmp184 -
          Ftmp146 * Ftmp604 - Ftmp146 * Ftmp66 * M[53] - Ftmp147 * x * M[48] +
          Ftmp15 * M[4] + Ftmp155 * Ftmp9 - Ftmp157 * Ftmp5 - Ftmp161 * M[34] -
          Ftmp162 * Ftmp26 + Ftmp168 * Ftmp61 + Ftmp172 * Ftmp597 * M[74] +
          Ftmp173 * Ftmp594 + Ftmp174 * Ftmp596 - Ftmp176 * Ftmp485 + Ftmp184 * Ftmp606 +
          Ftmp185 * Ftmp496 + Ftmp189 * Ftmp598 - Ftmp2 * Ftmp26 * M[2] -
          Ftmp2 * Ftmp476 - Ftmp2 * Ftmp591 - Ftmp201 * Ftmp601 + Ftmp215 * M[112] +
          Ftmp219 * z - Ftmp243 * Ftmp5 + Ftmp25 * M[16] - Ftmp250 * M[83] -
          Ftmp252 * Ftmp26 - Ftmp253 * Ftmp485 * M[110] - Ftmp253 * Ftmp66 * M[117] -
          Ftmp256 * Ftmp485 - Ftmp26 * Ftmp487 - Ftmp26 * Ftmp512 - Ftmp26 * Ftmp531 -
          Ftmp26 * Ftmp64 - Ftmp26 * (Ftmp47 + Ftmp77) * M[18] -
          Ftmp26 * (Ftmp537 + Ftmp545) * M[43] - Ftmp26 * (Ftmp575 + Ftmp582) * M[92] -
          Ftmp26 * (Ftmp144 + Ftmp175 - 13230.0 * Ftmp72) * M[54] -
          Ftmp26 * (Ftmp18 + Ftmp268 + Ftmp534) * M[41] -
          Ftmp26 * (Ftmp444 + Ftmp570 + Ftmp583) * M[90] -
          Ftmp26 * (-3648645.0 * Ftmp169 + 1964655.0 * Ftmp171 + Ftmp245 + Ftmp255) *
                M[118] -
          Ftmp26 * (Ftmp188 + Ftmp431 + Ftmp437 + Ftmp581) * M[99] -
          Ftmp26 * (Ftmp277 + Ftmp285 + Ftmp317 + Ftmp81) * M[52] -
          Ftmp26 * (Ftmp394 + Ftmp451 + Ftmp453 + Ftmp584) * M[114] -
          Ftmp26 * (Ftmp449 + Ftmp456 + Ftmp585 + Ftmp619) * M[116] -
          Ftmp26 * (Ftmp403 + Ftmp448 + Ftmp459 + Ftmp580 + Ftmp619) * M[103] -
          Ftmp26 * (-Ftmp12 * Ftmp338 + Ftmp333 + Ftmp466 + Ftmp475 + Ftmp590 + Ftmp97) *
                M[101] +
          Ftmp27 * Ftmp476 + Ftmp27 * Ftmp591 - Ftmp282 * x - Ftmp284 * x -
          Ftmp286 * Ftmp485 * M[46] - Ftmp287 * Ftmp480 + Ftmp29 * z + Ftmp299 * Ftmp9 +
          Ftmp305 * Ftmp9 - Ftmp308 * M[37] - Ftmp309 * M[39] + Ftmp314 * Ftmp61 +
          Ftmp315 * M[61] + Ftmp316 * M[63] - Ftmp4 * x * M[0] - Ftmp4 * y * M[1] -
          Ftmp42 * x - Ftmp438 * M[93] - Ftmp443 * M[97] - Ftmp446 * M[86] -
          Ftmp447 * M[88] - Ftmp450 * Ftmp485 * M[108] - Ftmp452 * Ftmp485 * M[106] -
          Ftmp46 * Ftmp480 - Ftmp47 * Ftmp9 * M[13] - Ftmp473 * M[95] + Ftmp477 * M[11] +
          Ftmp479 * z - Ftmp48 * Ftmp593 - Ftmp48 * Ftmp66 * M[17] - Ftmp48 * Ftmp88 -
          Ftmp483 * y - Ftmp484 * y - Ftmp485 * Ftmp79 + Ftmp488 * Ftmp595 -
          Ftmp489 * Ftmp5 - Ftmp49 * x * M[14] + Ftmp494 * Ftmp596 + Ftmp496 * Ftmp90 -
          Ftmp5 * Ftmp514 - Ftmp5 * Ftmp532 - Ftmp5 * Ftmp58 - Ftmp501 * Ftmp600 +
          Ftmp505 * M[36] + Ftmp507 * z - Ftmp510 * y - Ftmp511 * y + Ftmp513 * Ftmp595 +
          Ftmp516 * Ftmp596 - Ftmp519 * Ftmp600 + Ftmp521 * M[85] + Ftmp524 * z -
          Ftmp536 * y - Ftmp539 * y - Ftmp541 * y - Ftmp542 * M[40] - Ftmp543 * M[42] -
          Ftmp544 * M[51] + Ftmp551 * M[72] + Ftmp56 * Ftmp9 - Ftmp568 * M[98] -
          Ftmp569 * M[102] - Ftmp573 * M[89] - Ftmp576 * M[91] - Ftmp578 * M[115] -
          Ftmp579 * M[113] - Ftmp589 * M[100] + Ftmp592 * z * M[8] + Ftmp592 * M[18] +
          Ftmp593 * Ftmp599 + Ftmp594 * Ftmp75 + Ftmp595 * (Ftmp297 + Ftmp549) * M[62] +
          Ftmp595 * (Ftmp303 + Ftmp548) * M[64] + Ftmp596 * Ftmp76 +
          Ftmp596 * (Ftmp148 + Ftmp312 + Ftmp323) * M[79] +
          Ftmp596 * (Ftmp291 + Ftmp296 + Ftmp332) * M[66] +
          Ftmp596 * (Ftmp302 + Ftmp547 + Ftmp611) * M[68] + Ftmp597 * Ftmp74 * M[27] +
          Ftmp598 * Ftmp94 + Ftmp598 * (Ftmp311 + Ftmp323 + Ftmp611) * M[73] +
          Ftmp599 * Ftmp88 - Ftmp600 * (Ftmp343 + Ftmp348 + Ftmp552) * M[94] -
          Ftmp600 * (Ftmp347 + Ftmp349 + Ftmp553) * M[96] -
          Ftmp601 * (Ftmp197 + Ftmp339 + Ftmp349 + Ftmp353) * M[107] -
          Ftmp602 * (Ftmp110 - 1575.0 * Ftmp52) * M[47] -
          Ftmp602 * (Ftmp203 - 630630.0 * Ftmp204 + 121275.0 * Ftmp98) * M[109] +
          Ftmp603 * z * M[33] + Ftmp603 * M[54] + Ftmp604 * Ftmp606 +
          Ftmp607 * z * M[82] + Ftmp607 * M[118] + Ftmp608 * z * M[22] + Ftmp608 * M[41] +
          Ftmp609 * z * M[24] + Ftmp609 * M[43] + Ftmp61 * Ftmp69 + Ftmp610 * z * M[31] +
          Ftmp610 * M[52] + Ftmp612 * z * M[65] + Ftmp612 * M[99] + Ftmp614 * z * M[69] +
          Ftmp614 * M[103] + Ftmp615 * z * M[58] + Ftmp615 * M[90] + Ftmp616 * z * M[60] +
          Ftmp616 * M[92] + Ftmp617 * z * M[80] + Ftmp617 * M[116] + Ftmp618 * z * M[78] +
          Ftmp618 * M[114] + Ftmp620 * z * M[67] + Ftmp620 * M[101] - Ftmp63 * M[9] +
          Ftmp8 * Ftmp9;
}

template <>
void P2M<1, 3>(double x, double y, double z, double q, double* M, int order) {
  switch (order) {
    case 2:
      field_m1_P2M_2(x, y, z, q, M);
      break;
    case 3:
      field_m1_P2M_3(x, y, z, q, M);
      break;
    case 4:
      field_m1_P2M_4(x, y, z, q, M);
      break;
    case 5:
      field_m1_P2M_5(x, y, z, q, M);
      break;
    case 6:
      field_m1_P2M_6(x, y, z, q, M);
      break;
    case 7:
      field_m1_P2M_7(x, y, z, q, M);
      break;
  }
}
template <>
void M2M<1, 3>(double x, double y, double z, double* M, double* Ms, int order) {
  switch (order) {
    case 2:
      field_m1_M2M_2(x, y, z, M, Ms);
      break;
    case 3:
      field_m1_M2M_3(x, y, z, M, Ms);
      break;
    case 4:
      field_m1_M2M_4(x, y, z, M, Ms);
      break;
    case 5:
      field_m1_M2M_5(x, y, z, M, Ms);
      break;
    case 6:
      field_m1_M2M_6(x, y, z, M, Ms);
      break;
    case 7:
      field_m1_M2M_7(x, y, z, M, Ms);
      break;
  }
}
template <>
void M2L<1, 3>(double x, double y, double z, double* M, double* L, int order) {
  switch (order) {
    case 2:
      field_m1_M2L_2(x, y, z, M, L);
      break;
    case 3:
      field_m1_M2L_3(x, y, z, M, L);
      break;
    case 4:
      field_m1_M2L_4(x, y, z, M, L);
      break;
    case 5:
      field_m1_M2L_5(x, y, z, M, L);
      break;
    case 6:
      field_m1_M2L_6(x, y, z, M, L);
      break;
    case 7:
      field_m1_M2L_7(x, y, z, M, L);
      break;
  }
}
template <>
void L2L<1, 3>(double x, double y, double z, double* L, double* Ls, int order) {
  switch (order) {
    case 2:
      field_m1_L2L_2(x, y, z, L, Ls);
      break;
    case 3:
      field_m1_L2L_3(x, y, z, L, Ls);
      break;
    case 4:
      field_m1_L2L_4(x, y, z, L, Ls);
      break;
    case 5:
      field_m1_L2L_5(x, y, z, L, Ls);
      break;
    case 6:
      field_m1_L2L_6(x, y, z, L, Ls);
      break;
    case 7:
      field_m1_L2L_7(x, y, z, L, Ls);
      break;
  }
}
template <>
void L2P<1, 3>(double x, double y, double z, double* L, double* F, int order) {
  switch (order) {
    case 2:
      field_m1_L2P_2(x, y, z, L, F);
      break;
    case 3:
      field_m1_L2P_3(x, y, z, L, F);
      break;
    case 4:
      field_m1_L2P_4(x, y, z, L, F);
      break;
    case 5:
      field_m1_L2P_5(x, y, z, L, F);
      break;
    case 6:
      field_m1_L2P_6(x, y, z, L, F);
      break;
    case 7:
      field_m1_L2P_7(x, y, z, L, F);
      break;
  }
}
template <>
void M2P<1, 3>(double x, double y, double z, double* M, double* F, int order) {
  switch (order) {
    case 2:
      field_m1_M2P_2(x, y, z, M, F);
      break;
    case 3:
      field_m1_M2P_3(x, y, z, M, F);
      break;
    case 4:
      field_m1_M2P_4(x, y, z, M, F);
      break;
    case 5:
      field_m1_M2P_5(x, y, z, M, F);
      break;
    case 6:
      field_m1_M2P_6(x, y, z, M, F);
      break;
    case 7:
      field_m1_M2P_7(x, y, z, M, F);
      break;
  }
}

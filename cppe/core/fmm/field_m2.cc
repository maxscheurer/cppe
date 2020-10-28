#include "field_m2.hh"
#include <cmath>
void field_m2_P2M_3(double x, double y, double z, double q, double* M) {
  double Mtmp0 = (1.0 / 2.0) * q;
  double Mtmp1 = Mtmp0 * (x * x);
  double Mtmp2 = q * x;
  double Mtmp3 = Mtmp2 * y;
  double Mtmp4 = Mtmp0 * (y * y);
  double Mtmp5 = Mtmp0 * (z * z);
  double Mtmp6 = (1.0 / 6.0) * q;
  M[0] += Mtmp1;
  M[1] += Mtmp3;
  M[2] += Mtmp2 * z;
  M[3] += Mtmp4;
  M[4] += q * y * z;
  M[5] += Mtmp5;
  M[6] += -Mtmp6 * (x * x * x);
  M[7] += -Mtmp1 * y;
  M[8] += -Mtmp1 * z;
  M[9] += -Mtmp4 * x;
  M[10] += -Mtmp3 * z;
  M[11] += -Mtmp5 * x;
  M[12] += -Mtmp6 * (y * y * y);
  M[13] += -Mtmp4 * z;
  M[14] += -Mtmp5 * y;
  M[15] += -Mtmp6 * (z * z * z);
}
void field_m2_M2M_3(double x, double y, double z, double* M, double* Ms) {
#pragma omp atomic
  Ms[0] += M[0];
#pragma omp atomic
  Ms[1] += M[1];
#pragma omp atomic
  Ms[2] += M[2];
#pragma omp atomic
  Ms[3] += M[3];
#pragma omp atomic
  Ms[4] += M[4];
#pragma omp atomic
  Ms[5] += M[5];
#pragma omp atomic
  Ms[6] += x * M[0] + M[6];
#pragma omp atomic
  Ms[7] += x * M[1] + y * M[0] + M[7];
#pragma omp atomic
  Ms[8] += x * M[2] + z * M[0] + M[8];
#pragma omp atomic
  Ms[9] += x * M[3] + y * M[1] + M[9];
#pragma omp atomic
  Ms[10] += x * M[4] + y * M[2] + z * M[1] + M[10];
#pragma omp atomic
  Ms[11] += x * M[5] + z * M[2] + M[11];
#pragma omp atomic
  Ms[12] += y * M[3] + M[12];
#pragma omp atomic
  Ms[13] += y * M[4] + z * M[3] + M[13];
#pragma omp atomic
  Ms[14] += y * M[5] + z * M[4] + M[14];
#pragma omp atomic
  Ms[15] += z * M[5] + M[15];
}

void field_m2_M2L_3(double x, double y, double z, double* M, double* L) {
  double R = sqrt(x * x + y * y + z * z);
  double D[16];
  double Dtmp0  = -1.0 * pow(R, -3.0);
  double Dtmp1  = (x * x);
  double Dtmp2  = pow(R, -5.0);
  double Dtmp3  = 3.0 * Dtmp2;
  double Dtmp4  = Dtmp3 * x;
  double Dtmp5  = (y * y);
  double Dtmp6  = y * z;
  double Dtmp7  = -9.0 * Dtmp2;
  double Dtmp8  = 15.0 * pow(R, -7.0);
  double Dtmp9  = Dtmp1 * Dtmp8;
  double Dtmp10 = -Dtmp3;
  double Dtmp11 = Dtmp10 + Dtmp9;
  double Dtmp12 = Dtmp5 * Dtmp8;
  double Dtmp13 = Dtmp10 + Dtmp12;
  D[0]          = Dtmp0 + Dtmp1 * Dtmp3;
  D[1]          = Dtmp4 * y;
  D[2]          = Dtmp4 * z;
  D[3]          = Dtmp0 + Dtmp3 * Dtmp5;
  D[4]          = Dtmp3 * Dtmp6;
  D[5]          = -D[0] - D[3];
  D[6]          = -x * (Dtmp7 + Dtmp9);
  D[7]          = -Dtmp11 * y;
  D[8]          = -Dtmp11 * z;
  D[9]          = -1.0 * Dtmp13 * x;
  D[10]         = -Dtmp6 * Dtmp8 * x;
  D[11]         = -D[6] - D[9];
  D[12]         = -y * (Dtmp12 + Dtmp7);
  D[13]         = -Dtmp13 * z;
  D[14]         = -D[7] - D[12];
  D[15]         = -D[8] - D[13];
#pragma omp atomic
  L[0] += D[0] * M[0] + D[1] * M[1] + D[2] * M[2] + D[3] * M[3] + D[4] * M[4] +
          D[5] * M[5] + D[6] * M[6] + D[7] * M[7] + D[8] * M[8] + D[9] * M[9] +
          D[10] * M[10] + D[11] * M[11] + D[12] * M[12] + D[13] * M[13] + D[14] * M[14] +
          D[15] * M[15];
#pragma omp atomic
  L[1] += D[6] * M[0] + D[7] * M[1] + D[8] * M[2] + D[9] * M[3] + D[10] * M[4] +
          D[11] * M[5];
#pragma omp atomic
  L[2] += D[7] * M[0] + D[9] * M[1] + D[10] * M[2] + D[12] * M[3] + D[13] * M[4] +
          D[14] * M[5];
#pragma omp atomic
  L[3] += D[8] * M[0] + D[10] * M[1] + D[11] * M[2] + D[13] * M[3] + D[14] * M[4] +
          D[15] * M[5];
}

void field_m2_L2L_3(double x, double y, double z, double* L, double* Ls) {
#pragma omp atomic
  Ls[0] += x * L[1] + y * L[2] + z * L[3] + L[0];
#pragma omp atomic
  Ls[1] += L[1];
#pragma omp atomic
  Ls[2] += L[2];
#pragma omp atomic
  Ls[3] += L[3];
}

void field_m2_L2P_3(double x, double y, double z, double* L, double* F) {
#pragma omp atomic
  F[0] += -L[1];
#pragma omp atomic
  F[1] += -L[2];
#pragma omp atomic
  F[2] += -L[3];
}

void field_m2_M2P_3(double x, double y, double z, double* M, double* F) {
  double R      = sqrt(x * x + y * y + z * z);
  double Ftmp0  = pow(R, -5.0);
  double Ftmp1  = 3.0 * Ftmp0;
  double Ftmp2  = Ftmp1 * M[1];
  double Ftmp3  = Ftmp1 * z;
  double Ftmp4  = y * z;
  double Ftmp5  = pow(R, -7.0);
  double Ftmp6  = 15.0 * Ftmp5;
  double Ftmp7  = Ftmp6 * M[10];
  double Ftmp8  = z * M[4];
  double Ftmp9  = x * y;
  double Ftmp10 = Ftmp6 * Ftmp9;
  double Ftmp11 = (x * x);
  double Ftmp12 = Ftmp11 * Ftmp6;
  double Ftmp13 = z * M[2];
  double Ftmp14 = 105.0 * pow(R, -9.0);
  double Ftmp15 = Ftmp11 * Ftmp14;
  double Ftmp16 = -9.0 * Ftmp0;
  double Ftmp17 = Ftmp12 + Ftmp16;
  double Ftmp18 = -Ftmp1;
  double Ftmp19 = (y * y);
  double Ftmp20 = Ftmp19 * Ftmp6;
  double Ftmp21 = Ftmp18 + Ftmp20;
  double Ftmp22 = (z * z);
  double Ftmp23 = Ftmp22 * Ftmp6;
  double Ftmp24 = Ftmp18 + Ftmp23;
  double Ftmp25 = Ftmp21 * M[3];
  double Ftmp26 = Ftmp24 * M[5];
  double Ftmp27 = -45.0 * Ftmp5;
  double Ftmp28 = Ftmp15 + Ftmp27;
  double Ftmp29 = Ftmp28 * Ftmp9;
  double Ftmp30 = Ftmp14 * Ftmp19;
  double Ftmp31 = Ftmp27 + Ftmp30;
  double Ftmp32 = Ftmp31 * M[12];
  double Ftmp33 = -Ftmp6;
  double Ftmp34 = Ftmp14 * Ftmp22;
  double Ftmp35 = Ftmp33 + Ftmp34;
  double Ftmp36 = 1.0 * M[14];
  double Ftmp37 = Ftmp35 * Ftmp36;
  double Ftmp38 = x * z;
  double Ftmp39 = Ftmp28 * Ftmp38;
  double Ftmp40 = Ftmp30 + Ftmp33;
  double Ftmp41 = Ftmp40 * M[13];
  double Ftmp42 = Ftmp27 + Ftmp34;
  double Ftmp43 = Ftmp42 * M[15];
  double Ftmp44 = -75.0 * Ftmp5;
  double Ftmp45 = 1.0 * Ftmp11;
  double Ftmp46 = Ftmp40 * M[9];
  double Ftmp47 = Ftmp35 * M[11];
  double Ftmp48 = Ftmp12 + Ftmp18;
  double Ftmp49 = Ftmp16 + Ftmp20;
  double Ftmp50 = Ftmp48 * M[0];
  double Ftmp51 = 1.0 * Ftmp9;
  double Ftmp52 = Ftmp15 + Ftmp33;
  double Ftmp53 = Ftmp52 * M[8];
  double Ftmp54 = Ftmp52 * M[7];
  double Ftmp55 = x * M[2];
  double Ftmp56 = y * M[4];
  double Ftmp57 = Ftmp16 + Ftmp23;
  double Ftmp58 = 1.0 * Ftmp38;
#pragma omp atomic
  F[0] += Ftmp10 * Ftmp8 - Ftmp11 * (Ftmp15 + Ftmp44) * M[6] + Ftmp12 * Ftmp13 +
          Ftmp12 * y * M[1] - Ftmp15 * Ftmp4 * M[10] + Ftmp17 * x * M[0] + Ftmp17 * M[6] -
          Ftmp2 * y + Ftmp21 * M[9] + Ftmp24 * M[11] + Ftmp25 * x + Ftmp26 * x -
          Ftmp29 * M[7] - Ftmp3 * M[2] - Ftmp32 * Ftmp9 - Ftmp37 * Ftmp9 -
          Ftmp38 * Ftmp41 - Ftmp38 * Ftmp43 - Ftmp39 * M[8] + Ftmp4 * Ftmp7 -
          Ftmp45 * Ftmp46 - Ftmp45 * Ftmp47;
#pragma omp atomic
  F[1] += Ftmp10 * Ftmp13 - Ftmp19 * Ftmp37 - Ftmp19 * Ftmp54 -
          Ftmp19 * (Ftmp30 + Ftmp44) * M[12] - Ftmp2 * x + Ftmp20 * Ftmp8 +
          Ftmp20 * x * M[1] + Ftmp24 * M[14] + Ftmp26 * y - Ftmp29 * M[6] - Ftmp3 * M[4] -
          Ftmp30 * Ftmp38 * M[10] - Ftmp31 * Ftmp4 * M[13] - Ftmp31 * Ftmp51 * M[9] +
          Ftmp38 * Ftmp7 - Ftmp4 * Ftmp43 - Ftmp4 * Ftmp53 - Ftmp47 * Ftmp51 +
          Ftmp48 * M[7] + Ftmp49 * y * M[3] + Ftmp49 * M[12] + Ftmp50 * y;
#pragma omp atomic
  F[2] += -Ftmp1 * Ftmp55 - Ftmp1 * Ftmp56 + Ftmp10 * z * M[1] + Ftmp21 * M[13] -
          Ftmp22 * Ftmp41 - Ftmp22 * Ftmp53 - Ftmp22 * (Ftmp34 + Ftmp44) * M[15] +
          Ftmp23 * Ftmp55 + Ftmp23 * Ftmp56 + Ftmp25 * z - Ftmp32 * Ftmp4 -
          Ftmp34 * Ftmp9 * M[10] - Ftmp36 * Ftmp4 * Ftmp42 - Ftmp39 * M[6] -
          Ftmp4 * Ftmp54 - Ftmp42 * Ftmp58 * M[11] - Ftmp46 * Ftmp58 + Ftmp48 * M[8] +
          Ftmp50 * z + Ftmp57 * z * M[5] + Ftmp57 * M[15] + Ftmp7 * Ftmp9;
}

template <>
void P2P<2, 3>(double x, double y, double z, double* S, double* F) {
  double R      = sqrt(x * x + y * y + z * z);
  double Ftmp0  = pow(R, -5.0);
  double Ftmp1  = 3.0 * Ftmp0;
  double Ftmp2  = Ftmp1 * S[1];
  double Ftmp3  = Ftmp1 * z;
  double Ftmp4  = y * S[4];
  double Ftmp5  = 15.0 * pow(R, -7.0);
  double Ftmp6  = Ftmp5 * z;
  double Ftmp7  = Ftmp6 * x;
  double Ftmp8  = Ftmp5 * (x * x);
  double Ftmp9  = y * S[1];
  double Ftmp10 = -9.0 * Ftmp0;
  double Ftmp11 = -Ftmp1;
  double Ftmp12 = Ftmp5 * (y * y);
  double Ftmp13 = (Ftmp11 + Ftmp12) * S[3];
  double Ftmp14 = Ftmp5 * (z * z);
  double Ftmp15 = (Ftmp11 + Ftmp14) * S[5];
  double Ftmp16 = x * S[2];
  double Ftmp17 = (Ftmp11 + Ftmp8) * S[0];
#pragma omp atomic
  F[0] += Ftmp13 * x + Ftmp15 * x - Ftmp2 * y - Ftmp3 * S[2] + Ftmp4 * Ftmp7 +
          Ftmp8 * Ftmp9 + Ftmp8 * z * S[2] + x * (Ftmp10 + Ftmp8) * S[0];
#pragma omp atomic
  F[1] += Ftmp12 * x * S[1] + Ftmp12 * z * S[4] + Ftmp15 * y + Ftmp16 * Ftmp6 * y +
          Ftmp17 * y - Ftmp2 * x - Ftmp3 * S[4] + y * (Ftmp10 + Ftmp12) * S[3];
#pragma omp atomic
  F[2] += -Ftmp1 * Ftmp16 - Ftmp1 * Ftmp4 + Ftmp13 * z + Ftmp14 * Ftmp16 +
          Ftmp14 * Ftmp4 + Ftmp17 * z + Ftmp7 * Ftmp9 + z * (Ftmp10 + Ftmp14) * S[5];
}

void field_m2_P2M_4(double x, double y, double z, double q, double* M) {
  double Mtmp0  = (x * x);
  double Mtmp1  = (1.0 / 2.0) * q;
  double Mtmp2  = Mtmp0 * Mtmp1;
  double Mtmp3  = q * x;
  double Mtmp4  = Mtmp3 * y;
  double Mtmp5  = (y * y);
  double Mtmp6  = Mtmp1 * Mtmp5;
  double Mtmp7  = q * y;
  double Mtmp8  = (z * z);
  double Mtmp9  = Mtmp1 * Mtmp8;
  double Mtmp10 = (x * x * x);
  double Mtmp11 = (1.0 / 6.0) * q;
  double Mtmp12 = Mtmp10 * Mtmp11;
  double Mtmp13 = Mtmp2 * y;
  double Mtmp14 = Mtmp6 * x;
  double Mtmp15 = Mtmp9 * x;
  double Mtmp16 = (y * y * y);
  double Mtmp17 = Mtmp11 * Mtmp16;
  double Mtmp18 = (z * z * z);
  double Mtmp19 = (1.0 / 24.0) * q;
  double Mtmp20 = (1.0 / 6.0) * Mtmp7;
  double Mtmp21 = (1.0 / 4.0) * Mtmp0 * q;
  double Mtmp22 = (1.0 / 6.0) * Mtmp3;
  M[0] += Mtmp2;
  M[1] += Mtmp4;
  M[2] += Mtmp3 * z;
  M[3] += Mtmp6;
  M[4] += Mtmp7 * z;
  M[5] += Mtmp9;
  M[6] += -Mtmp12;
  M[7] += -Mtmp13;
  M[8] += -Mtmp2 * z;
  M[9] += -Mtmp14;
  M[10] += -Mtmp4 * z;
  M[11] += -Mtmp15;
  M[12] += -Mtmp17;
  M[13] += -Mtmp6 * z;
  M[14] += -Mtmp9 * y;
  M[15] += -Mtmp11 * Mtmp18;
  M[16] += Mtmp19 * (x * x * x * x);
  M[17] += Mtmp10 * Mtmp20;
  M[18] += Mtmp12 * z;
  M[19] += Mtmp21 * Mtmp5;
  M[20] += Mtmp13 * z;
  M[21] += Mtmp21 * Mtmp8;
  M[22] += Mtmp16 * Mtmp22;
  M[23] += Mtmp14 * z;
  M[24] += Mtmp15 * y;
  M[25] += Mtmp18 * Mtmp22;
  M[26] += Mtmp19 * (y * y * y * y);
  M[27] += Mtmp17 * z;
  M[28] += (1.0 / 4.0) * Mtmp5 * Mtmp8 * q;
  M[29] += Mtmp18 * Mtmp20;
  M[30] += Mtmp19 * (z * z * z * z);
}
void field_m2_M2M_4(double x, double y, double z, double* M, double* Ms) {
  double Mstmp0  = x * M[0];
  double Mstmp1  = x * M[1];
  double Mstmp2  = y * M[0];
  double Mstmp3  = x * M[2];
  double Mstmp4  = x * M[3];
  double Mstmp5  = y * M[1];
  double Mstmp6  = x * M[4];
  double Mstmp7  = y * M[2];
  double Mstmp8  = x * M[5];
  double Mstmp9  = y * M[3];
  double Mstmp10 = y * M[4];
  double Mstmp11 = y * M[5];
  double Mstmp12 = (1.0 / 2.0) * (x * x);
  double Mstmp13 = (y * y);
  double Mstmp14 = (1.0 / 2.0) * M[0];
  double Mstmp15 = (z * z);
  double Mstmp16 = (1.0 / 2.0) * Mstmp13;
  double Mstmp17 = (1.0 / 2.0) * Mstmp15;
#pragma omp atomic
  Ms[0] += M[0];
#pragma omp atomic
  Ms[1] += M[1];
#pragma omp atomic
  Ms[2] += M[2];
#pragma omp atomic
  Ms[3] += M[3];
#pragma omp atomic
  Ms[4] += M[4];
#pragma omp atomic
  Ms[5] += M[5];
#pragma omp atomic
  Ms[6] += Mstmp0 + M[6];
#pragma omp atomic
  Ms[7] += Mstmp1 + Mstmp2 + M[7];
#pragma omp atomic
  Ms[8] += Mstmp3 + z * M[0] + M[8];
#pragma omp atomic
  Ms[9] += Mstmp4 + Mstmp5 + M[9];
#pragma omp atomic
  Ms[10] += Mstmp6 + Mstmp7 + z * M[1] + M[10];
#pragma omp atomic
  Ms[11] += Mstmp8 + z * M[2] + M[11];
#pragma omp atomic
  Ms[12] += Mstmp9 + M[12];
#pragma omp atomic
  Ms[13] += Mstmp10 + z * M[3] + M[13];
#pragma omp atomic
  Ms[14] += Mstmp11 + z * M[4] + M[14];
#pragma omp atomic
  Ms[15] += z * M[5] + M[15];
#pragma omp atomic
  Ms[16] += Mstmp12 * M[0] + x * M[6] + M[16];
#pragma omp atomic
  Ms[17] += Mstmp0 * y + Mstmp12 * M[1] + x * M[7] + y * M[6] + M[17];
#pragma omp atomic
  Ms[18] += Mstmp0 * z + Mstmp12 * M[2] + x * M[8] + z * M[6] + M[18];
#pragma omp atomic
  Ms[19] += Mstmp1 * y + Mstmp12 * M[3] + Mstmp13 * Mstmp14 + x * M[9] + y * M[7] + M[19];
#pragma omp atomic
  Ms[20] += Mstmp1 * z + Mstmp12 * M[4] + Mstmp2 * z + Mstmp3 * y + x * M[10] + y * M[8] +
            z * M[7] + M[20];
#pragma omp atomic
  Ms[21] +=
        Mstmp12 * M[5] + Mstmp14 * Mstmp15 + Mstmp3 * z + x * M[11] + z * M[8] + M[21];
#pragma omp atomic
  Ms[22] += Mstmp16 * M[1] + Mstmp4 * y + x * M[12] + y * M[9] + M[22];
#pragma omp atomic
  Ms[23] += Mstmp16 * M[2] + Mstmp4 * z + Mstmp5 * z + Mstmp6 * y + x * M[13] +
            y * M[10] + z * M[9] + M[23];
#pragma omp atomic
  Ms[24] += Mstmp17 * M[1] + Mstmp6 * z + Mstmp7 * z + Mstmp8 * y + x * M[14] +
            y * M[11] + z * M[10] + M[24];
#pragma omp atomic
  Ms[25] += Mstmp17 * M[2] + Mstmp8 * z + x * M[15] + z * M[11] + M[25];
#pragma omp atomic
  Ms[26] += Mstmp16 * M[3] + y * M[12] + M[26];
#pragma omp atomic
  Ms[27] += Mstmp16 * M[4] + Mstmp9 * z + y * M[13] + z * M[12] + M[27];
#pragma omp atomic
  Ms[28] += Mstmp10 * z + Mstmp16 * M[5] + Mstmp17 * M[3] + y * M[14] + z * M[13] + M[28];
#pragma omp atomic
  Ms[29] += Mstmp11 * z + Mstmp17 * M[4] + y * M[15] + z * M[14] + M[29];
#pragma omp atomic
  Ms[30] += Mstmp17 * M[5] + z * M[15] + M[30];
}

void field_m2_M2L_4(double x, double y, double z, double* M, double* L) {
  double R = sqrt(x * x + y * y + z * z);
  double D[31];
  double Dtmp0  = -1.0 * pow(R, -3.0);
  double Dtmp1  = (x * x);
  double Dtmp2  = pow(R, -5.0);
  double Dtmp3  = 3.0 * Dtmp2;
  double Dtmp4  = Dtmp3 * x;
  double Dtmp5  = (y * y);
  double Dtmp6  = y * z;
  double Dtmp7  = 9.0 * Dtmp2;
  double Dtmp8  = -Dtmp7;
  double Dtmp9  = pow(R, -7.0);
  double Dtmp10 = 15.0 * Dtmp9;
  double Dtmp11 = Dtmp1 * Dtmp10;
  double Dtmp12 = -Dtmp3;
  double Dtmp13 = Dtmp11 + Dtmp12;
  double Dtmp14 = Dtmp10 * Dtmp5;
  double Dtmp15 = Dtmp12 + Dtmp14;
  double Dtmp16 = 1.0 * x;
  double Dtmp17 = 105.0 * pow(R, -9.0);
  double Dtmp18 = 90.0 * Dtmp9;
  double Dtmp19 = -45.0 * Dtmp9;
  double Dtmp20 = Dtmp1 * Dtmp17;
  double Dtmp21 = x * (Dtmp19 + Dtmp20);
  double Dtmp22 = -Dtmp10;
  double Dtmp23 = Dtmp17 * Dtmp5;
  double Dtmp24 = Dtmp19 + Dtmp23;
  D[0]          = Dtmp0 + Dtmp1 * Dtmp3;
  D[1]          = Dtmp4 * y;
  D[2]          = Dtmp4 * z;
  D[3]          = Dtmp0 + Dtmp3 * Dtmp5;
  D[4]          = Dtmp3 * Dtmp6;
  D[5]          = -D[0] - D[3];
  D[6]          = -x * (Dtmp11 + Dtmp8);
  D[7]          = -Dtmp13 * y;
  D[8]          = -Dtmp13 * z;
  D[9]          = -Dtmp15 * Dtmp16;
  D[10]         = -Dtmp10 * Dtmp6 * x;
  D[11]         = -D[6] - D[9];
  D[12]         = -y * (Dtmp14 + Dtmp8);
  D[13]         = -Dtmp15 * z;
  D[14]         = -D[7] - D[12];
  D[15]         = -D[8] - D[13];
  D[16]         = -Dtmp1 * Dtmp18 + Dtmp17 * (x * x * x * x) + Dtmp7;
  D[17]         = Dtmp21 * y;
  D[18]         = Dtmp21 * z;
  D[19]         = -Dtmp11 - Dtmp14 + Dtmp20 * Dtmp5 + Dtmp3;
  D[20]         = Dtmp6 * (Dtmp20 + Dtmp22);
  D[21]         = -D[16] - D[19];
  D[22]         = Dtmp16 * Dtmp24 * y;
  D[23]         = Dtmp16 * z * (Dtmp22 + Dtmp23);
  D[24]         = -D[17] - D[22];
  D[25]         = -D[18] - D[23];
  D[26]         = Dtmp17 * (y * y * y * y) - Dtmp18 * Dtmp5 + Dtmp7;
  D[27]         = Dtmp24 * Dtmp6;
  D[28]         = -D[19] - D[26];
  D[29]         = -D[20] - D[27];
  D[30]         = -D[21] - D[28];
#pragma omp atomic
  L[0] += D[0] * M[0] + D[1] * M[1] + D[2] * M[2] + D[3] * M[3] + D[4] * M[4] +
          D[5] * M[5] + D[6] * M[6] + D[7] * M[7] + D[8] * M[8] + D[9] * M[9] +
          D[10] * M[10] + D[11] * M[11] + D[12] * M[12] + D[13] * M[13] + D[14] * M[14] +
          D[15] * M[15] + D[16] * M[16] + D[17] * M[17] + D[18] * M[18] + D[19] * M[19] +
          D[20] * M[20] + D[21] * M[21] + D[22] * M[22] + D[23] * M[23] + D[24] * M[24] +
          D[25] * M[25] + D[26] * M[26] + D[27] * M[27] + D[28] * M[28] + D[29] * M[29] +
          D[30] * M[30];
#pragma omp atomic
  L[1] += D[6] * M[0] + D[7] * M[1] + D[8] * M[2] + D[9] * M[3] + D[10] * M[4] +
          D[11] * M[5] + D[16] * M[6] + D[17] * M[7] + D[18] * M[8] + D[19] * M[9] +
          D[20] * M[10] + D[21] * M[11] + D[22] * M[12] + D[23] * M[13] + D[24] * M[14] +
          D[25] * M[15];
#pragma omp atomic
  L[2] += D[7] * M[0] + D[9] * M[1] + D[10] * M[2] + D[12] * M[3] + D[13] * M[4] +
          D[14] * M[5] + D[17] * M[6] + D[19] * M[7] + D[20] * M[8] + D[22] * M[9] +
          D[23] * M[10] + D[24] * M[11] + D[26] * M[12] + D[27] * M[13] + D[28] * M[14] +
          D[29] * M[15];
#pragma omp atomic
  L[3] += D[8] * M[0] + D[10] * M[1] + D[11] * M[2] + D[13] * M[3] + D[14] * M[4] +
          D[15] * M[5] + D[18] * M[6] + D[20] * M[7] + D[21] * M[8] + D[23] * M[9] +
          D[24] * M[10] + D[25] * M[11] + D[27] * M[12] + D[28] * M[13] + D[29] * M[14] +
          D[30] * M[15];
#pragma omp atomic
  L[4] += D[16] * M[0] + D[17] * M[1] + D[18] * M[2] + D[19] * M[3] + D[20] * M[4] +
          D[21] * M[5];
#pragma omp atomic
  L[5] += D[17] * M[0] + D[19] * M[1] + D[20] * M[2] + D[22] * M[3] + D[23] * M[4] +
          D[24] * M[5];
#pragma omp atomic
  L[6] += D[18] * M[0] + D[20] * M[1] + D[21] * M[2] + D[23] * M[3] + D[24] * M[4] +
          D[25] * M[5];
#pragma omp atomic
  L[7] += D[19] * M[0] + D[22] * M[1] + D[23] * M[2] + D[26] * M[3] + D[27] * M[4] +
          D[28] * M[5];
#pragma omp atomic
  L[8] += D[20] * M[0] + D[23] * M[1] + D[24] * M[2] + D[27] * M[3] + D[28] * M[4] +
          D[29] * M[5];
#pragma omp atomic
  L[9] += D[21] * M[0] + D[24] * M[1] + D[25] * M[2] + D[28] * M[3] + D[29] * M[4] +
          D[30] * M[5];
}

void field_m2_L2L_4(double x, double y, double z, double* L, double* Ls) {
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

void field_m2_L2P_4(double x, double y, double z, double* L, double* F) {
#pragma omp atomic
  F[0] += -x * L[4] - y * L[5] - z * L[6] - L[1];
#pragma omp atomic
  F[1] += -x * L[5] - y * L[7] - z * L[8] - L[2];
#pragma omp atomic
  F[2] += -x * L[6] - y * L[8] - z * L[9] - L[3];
}

void field_m2_M2P_4(double x, double y, double z, double* M, double* F) {
  double R       = sqrt(x * x + y * y + z * z);
  double Ftmp0   = pow(R, -5.0);
  double Ftmp1   = 3.0 * Ftmp0;
  double Ftmp2   = Ftmp1 * M[1];
  double Ftmp3   = Ftmp1 * z;
  double Ftmp4   = y * z;
  double Ftmp5   = pow(R, -7.0);
  double Ftmp6   = 15.0 * Ftmp5;
  double Ftmp7   = Ftmp6 * M[10];
  double Ftmp8   = y * M[4];
  double Ftmp9   = x * z;
  double Ftmp10  = Ftmp6 * Ftmp9;
  double Ftmp11  = (x * x);
  double Ftmp12  = Ftmp11 * Ftmp6;
  double Ftmp13  = y * M[1];
  double Ftmp14  = pow(R, -9.0);
  double Ftmp15  = 105.0 * Ftmp14;
  double Ftmp16  = Ftmp11 * Ftmp15;
  double Ftmp17  = -9.0 * Ftmp0;
  double Ftmp18  = Ftmp12 + Ftmp17;
  double Ftmp19  = -Ftmp1;
  double Ftmp20  = (y * y);
  double Ftmp21  = Ftmp20 * Ftmp6;
  double Ftmp22  = Ftmp19 + Ftmp21;
  double Ftmp23  = (z * z);
  double Ftmp24  = Ftmp23 * Ftmp6;
  double Ftmp25  = Ftmp19 + Ftmp24;
  double Ftmp26  = Ftmp22 * M[3];
  double Ftmp27  = Ftmp25 * M[5];
  double Ftmp28  = 45.0 * Ftmp5;
  double Ftmp29  = -Ftmp28;
  double Ftmp30  = Ftmp16 + Ftmp29;
  double Ftmp31  = Ftmp30 * M[17];
  double Ftmp32  = Ftmp15 * Ftmp20;
  double Ftmp33  = Ftmp29 + Ftmp32;
  double Ftmp34  = Ftmp33 * y;
  double Ftmp35  = 1.0 * M[22];
  double Ftmp36  = 35.0 * Ftmp14;
  double Ftmp37  = 3.0 * M[24];
  double Ftmp38  = Ftmp37 * (Ftmp23 * Ftmp36 - 5.0 * Ftmp5);
  double Ftmp39  = Ftmp30 * M[18];
  double Ftmp40  = 1.0 * z;
  double Ftmp41  = -Ftmp6;
  double Ftmp42  = Ftmp32 + Ftmp41;
  double Ftmp43  = Ftmp42 * M[23];
  double Ftmp44  = Ftmp15 * Ftmp23;
  double Ftmp45  = Ftmp29 + Ftmp44;
  double Ftmp46  = Ftmp40 * Ftmp45;
  double Ftmp47  = Ftmp30 * x;
  double Ftmp48  = Ftmp47 * y;
  double Ftmp49  = Ftmp34 * M[12];
  double Ftmp50  = Ftmp41 + Ftmp44;
  double Ftmp51  = 1.0 * Ftmp50 * M[14];
  double Ftmp52  = x * y;
  double Ftmp53  = Ftmp47 * z;
  double Ftmp54  = Ftmp42 * M[13];
  double Ftmp55  = Ftmp45 * M[15];
  double Ftmp56  = 315.0 * Ftmp14;
  double Ftmp57  = -Ftmp56;
  double Ftmp58  = pow(R, -11.0);
  double Ftmp59  = 945.0 * Ftmp58;
  double Ftmp60  = Ftmp11 * Ftmp59;
  double Ftmp61  = Ftmp57 + Ftmp60;
  double Ftmp62  = Ftmp9 * y;
  double Ftmp63  = Ftmp20 * Ftmp59;
  double Ftmp64  = Ftmp57 + Ftmp63;
  double Ftmp65  = Ftmp64 * M[27];
  double Ftmp66  = -75.0 * Ftmp5;
  double Ftmp67  = 1.0 * Ftmp11;
  double Ftmp68  = Ftmp42 * M[9];
  double Ftmp69  = Ftmp50 * M[11];
  double Ftmp70  = -525.0 * Ftmp14;
  double Ftmp71  = Ftmp11 * (Ftmp60 + Ftmp70);
  double Ftmp72  = y * M[17];
  double Ftmp73  = z * M[18];
  double Ftmp74  = y * M[29];
  double Ftmp75  = Ftmp23 * Ftmp59;
  double Ftmp76  = Ftmp57 + Ftmp75;
  double Ftmp77  = Ftmp40 * Ftmp76;
  double Ftmp78  = Ftmp11 * y;
  double Ftmp79  = Ftmp35 * Ftmp64;
  double Ftmp80  = 315.0 * Ftmp23 * Ftmp58;
  double Ftmp81  = Ftmp37 * (-Ftmp36 + Ftmp80);
  double Ftmp82  = Ftmp11 * Ftmp40;
  double Ftmp83  = -Ftmp15;
  double Ftmp84  = (Ftmp63 + Ftmp83) * M[23];
  double Ftmp85  = Ftmp76 * M[25];
  double Ftmp86  = 225.0 * Ftmp5;
  double Ftmp87  = Ftmp59 * (x * x * x * x);
  double Ftmp88  = 1050.0 * Ftmp14;
  double Ftmp89  = Ftmp59 * (y * y * y * y);
  double Ftmp90  = 630.0 * Ftmp14;
  double Ftmp91  = (-Ftmp20 * Ftmp90 + Ftmp28 + Ftmp89) * M[26];
  double Ftmp92  = Ftmp59 * (z * z * z * z);
  double Ftmp93  = (-Ftmp23 * Ftmp90 + Ftmp28 + Ftmp92) * M[30];
  double Ftmp94  = -Ftmp20 * Ftmp56;
  double Ftmp95  = Ftmp20 * Ftmp60;
  double Ftmp96  = -Ftmp16;
  double Ftmp97  = Ftmp28 + Ftmp96;
  double Ftmp98  = -Ftmp23 * Ftmp56;
  double Ftmp99  = Ftmp23 * Ftmp60;
  double Ftmp100 = -Ftmp44;
  double Ftmp101 = Ftmp100 + Ftmp6;
  double Ftmp102 = -Ftmp32;
  double Ftmp103 = Ftmp23 * Ftmp63;
  double Ftmp104 = Ftmp102 + Ftmp103;
  double Ftmp105 = Ftmp12 + Ftmp19;
  double Ftmp106 = Ftmp17 + Ftmp21;
  double Ftmp107 = Ftmp105 * M[0];
  double Ftmp108 = Ftmp35 * x;
  double Ftmp109 = Ftmp16 + Ftmp41;
  double Ftmp110 = Ftmp109 * M[20];
  double Ftmp111 = Ftmp33 * M[27];
  double Ftmp112 = 1.0 * x;
  double Ftmp113 = Ftmp109 * M[8];
  double Ftmp114 = Ftmp61 * x;
  double Ftmp115 = Ftmp109 * M[7];
  double Ftmp116 = Ftmp20 * z;
  double Ftmp117 = (Ftmp60 + Ftmp83) * M[20];
  double Ftmp118 = Ftmp63 + Ftmp70;
  double Ftmp119 = Ftmp40 * Ftmp52;
  double Ftmp120 = (-Ftmp11 * Ftmp90 + Ftmp28 + Ftmp87) * M[16];
  double Ftmp121 = Ftmp102 + Ftmp95;
  double Ftmp122 = -Ftmp11 * Ftmp56 + Ftmp28;
  double Ftmp123 = x * M[2];
  double Ftmp124 = Ftmp17 + Ftmp24;
  double Ftmp125 = Ftmp112 * M[25];
  double Ftmp126 = 1.0 * Ftmp74;
  double Ftmp127 = Ftmp23 * y;
  double Ftmp128 = Ftmp23 * (Ftmp70 + Ftmp75);
#pragma omp atomic
  F[0] += Ftmp10 * Ftmp8 - Ftmp11 * (Ftmp16 + Ftmp66) * M[6] + Ftmp12 * Ftmp13 +
          Ftmp12 * z * M[2] - Ftmp16 * Ftmp4 * M[10] + Ftmp18 * x * M[0] + Ftmp18 * M[6] -
          Ftmp2 * y + Ftmp22 * M[9] + Ftmp25 * M[11] + Ftmp26 * x + Ftmp27 * x -
          Ftmp3 * M[2] - Ftmp31 * y - Ftmp34 * Ftmp35 - Ftmp38 * y - Ftmp39 * z +
          Ftmp4 * Ftmp7 - Ftmp40 * Ftmp43 - Ftmp46 * M[25] - Ftmp48 * M[7] - Ftmp49 * x -
          Ftmp51 * Ftmp52 - Ftmp53 * M[8] - Ftmp54 * Ftmp9 - Ftmp55 * Ftmp9 +
          Ftmp61 * Ftmp62 * M[20] + Ftmp62 * Ftmp65 - Ftmp67 * Ftmp68 - Ftmp67 * Ftmp69 +
          Ftmp71 * Ftmp72 + Ftmp71 * Ftmp73 + Ftmp74 * Ftmp77 * x + Ftmp78 * Ftmp79 +
          Ftmp78 * Ftmp81 + Ftmp82 * Ftmp84 + Ftmp82 * Ftmp85 + Ftmp91 * x + Ftmp93 * x +
          x * (Ftmp101 + Ftmp104) * M[28] + x * (Ftmp94 + Ftmp95 + Ftmp97) * M[19] +
          x * (Ftmp97 + Ftmp98 + Ftmp99) * M[21] +
          x * (-Ftmp11 * Ftmp88 + Ftmp86 + Ftmp87) * M[16];
#pragma omp atomic
  F[1] += Ftmp105 * M[7] + Ftmp106 * y * M[3] + Ftmp106 * M[12] + Ftmp107 * y +
          Ftmp108 * Ftmp118 * Ftmp20 - Ftmp108 * Ftmp33 - Ftmp110 * z - Ftmp111 * z -
          Ftmp112 * Ftmp34 * M[9] - Ftmp112 * Ftmp69 * y - Ftmp113 * Ftmp4 +
          Ftmp114 * Ftmp20 * M[17] + Ftmp114 * Ftmp73 * y - Ftmp115 * Ftmp20 +
          Ftmp116 * Ftmp117 + Ftmp116 * Ftmp118 * M[27] + Ftmp119 * Ftmp64 * M[23] +
          Ftmp119 * Ftmp85 + Ftmp120 * y - Ftmp2 * x - Ftmp20 * Ftmp51 +
          Ftmp20 * Ftmp77 * M[29] + Ftmp20 * Ftmp81 * x -
          Ftmp20 * (Ftmp32 + Ftmp66) * M[12] + Ftmp21 * x * M[1] + Ftmp21 * z * M[4] +
          Ftmp25 * M[14] + Ftmp27 * y - Ftmp3 * M[4] - Ftmp31 * x -
          Ftmp32 * Ftmp9 * M[10] - Ftmp34 * z * M[13] - Ftmp38 * x - Ftmp4 * Ftmp55 -
          Ftmp46 * M[29] - Ftmp48 * M[6] + Ftmp6 * Ftmp62 * M[2] + Ftmp7 * Ftmp9 +
          Ftmp93 * y + y * (Ftmp121 + Ftmp122) * M[19] +
          y * (Ftmp101 + Ftmp96 + Ftmp99) * M[21] +
          y * (Ftmp104 + Ftmp28 + Ftmp98) * M[28] +
          y * (-Ftmp20 * Ftmp88 + Ftmp86 + Ftmp89) * M[26];
#pragma omp atomic
  F[2] += -Ftmp1 * Ftmp123 - Ftmp1 * Ftmp8 + Ftmp10 * Ftmp13 + Ftmp105 * M[8] +
          Ftmp107 * z - Ftmp110 * y - Ftmp111 * y + Ftmp112 * Ftmp23 * Ftmp84 -
          Ftmp112 * Ftmp43 - Ftmp113 * Ftmp23 + Ftmp114 * Ftmp23 * M[18] -
          Ftmp115 * Ftmp4 + Ftmp117 * Ftmp127 + Ftmp120 * z + Ftmp123 * Ftmp24 +
          Ftmp124 * z * M[5] + Ftmp124 * M[15] + Ftmp125 * Ftmp128 - Ftmp125 * Ftmp45 +
          Ftmp126 * Ftmp128 - Ftmp126 * Ftmp45 + Ftmp127 * Ftmp65 + Ftmp22 * M[13] -
          Ftmp23 * Ftmp54 - Ftmp23 * (Ftmp44 + Ftmp66) * M[15] + Ftmp24 * Ftmp8 +
          Ftmp26 * z + Ftmp37 * Ftmp62 * (Ftmp80 + Ftmp83) - Ftmp39 * x -
          Ftmp40 * Ftmp68 * x - Ftmp44 * Ftmp52 * M[10] - Ftmp46 * x * M[11] -
          Ftmp46 * y * M[14] - Ftmp49 * z + Ftmp52 * Ftmp7 - Ftmp53 * M[6] +
          Ftmp61 * Ftmp72 * Ftmp9 + Ftmp62 * Ftmp79 + Ftmp91 * z +
          z * (Ftmp100 + Ftmp122 + Ftmp99) * M[21] +
          z * (Ftmp121 + Ftmp6 + Ftmp96) * M[19] +
          z * (-Ftmp23 * Ftmp88 + Ftmp86 + Ftmp92) * M[30] +
          z * (Ftmp100 + Ftmp103 + Ftmp28 + Ftmp94) * M[28];
}

void field_m2_P2M_5(double x, double y, double z, double q, double* M) {
  double Mtmp0  = (x * x);
  double Mtmp1  = (1.0 / 2.0) * q;
  double Mtmp2  = Mtmp0 * Mtmp1;
  double Mtmp3  = q * x;
  double Mtmp4  = Mtmp3 * y;
  double Mtmp5  = Mtmp3 * z;
  double Mtmp6  = (y * y);
  double Mtmp7  = Mtmp1 * Mtmp6;
  double Mtmp8  = q * y;
  double Mtmp9  = Mtmp8 * z;
  double Mtmp10 = (z * z);
  double Mtmp11 = Mtmp1 * Mtmp10;
  double Mtmp12 = (x * x * x);
  double Mtmp13 = (1.0 / 6.0) * q;
  double Mtmp14 = Mtmp12 * Mtmp13;
  double Mtmp15 = Mtmp2 * y;
  double Mtmp16 = Mtmp7 * x;
  double Mtmp17 = Mtmp11 * x;
  double Mtmp18 = (y * y * y);
  double Mtmp19 = Mtmp13 * Mtmp18;
  double Mtmp20 = (z * z * z);
  double Mtmp21 = (x * x * x * x);
  double Mtmp22 = (1.0 / 24.0) * q;
  double Mtmp23 = Mtmp21 * Mtmp22;
  double Mtmp24 = (1.0 / 6.0) * Mtmp8;
  double Mtmp25 = Mtmp6 * q;
  double Mtmp26 = (1.0 / 4.0) * Mtmp0;
  double Mtmp27 = Mtmp25 * Mtmp26;
  double Mtmp28 = Mtmp10 * q;
  double Mtmp29 = (1.0 / 6.0) * Mtmp3;
  double Mtmp30 = (y * y * y * y);
  double Mtmp31 = Mtmp22 * Mtmp30;
  double Mtmp32 = (1.0 / 4.0) * Mtmp10;
  double Mtmp33 = (z * z * z * z);
  double Mtmp34 = (1.0 / 120.0) * q;
  double Mtmp35 = (1.0 / 24.0) * Mtmp8;
  double Mtmp36 = (1.0 / 12.0) * Mtmp12;
  double Mtmp37 = (1.0 / 12.0) * Mtmp18;
  double Mtmp38 = Mtmp0 * q;
  double Mtmp39 = (1.0 / 12.0) * Mtmp20;
  double Mtmp40 = (1.0 / 24.0) * Mtmp3;
  M[0] += Mtmp2;
  M[1] += Mtmp4;
  M[2] += Mtmp5;
  M[3] += Mtmp7;
  M[4] += Mtmp9;
  M[5] += Mtmp11;
  M[6] += -Mtmp14;
  M[7] += -Mtmp15;
  M[8] += -Mtmp2 * z;
  M[9] += -Mtmp16;
  M[10] += -Mtmp4 * z;
  M[11] += -Mtmp17;
  M[12] += -Mtmp19;
  M[13] += -Mtmp7 * z;
  M[14] += -Mtmp11 * y;
  M[15] += -Mtmp13 * Mtmp20;
  M[16] += Mtmp23;
  M[17] += Mtmp12 * Mtmp24;
  M[18] += Mtmp14 * z;
  M[19] += Mtmp27;
  M[20] += Mtmp15 * z;
  M[21] += Mtmp26 * Mtmp28;
  M[22] += Mtmp18 * Mtmp29;
  M[23] += Mtmp16 * z;
  M[24] += Mtmp17 * y;
  M[25] += Mtmp20 * Mtmp29;
  M[26] += Mtmp31;
  M[27] += Mtmp19 * z;
  M[28] += Mtmp25 * Mtmp32;
  M[29] += Mtmp20 * Mtmp24;
  M[30] += Mtmp22 * Mtmp33;
  M[31] += -Mtmp34 * (x * x * x * x * x);
  M[32] += -Mtmp21 * Mtmp35;
  M[33] += -Mtmp23 * z;
  M[34] += -Mtmp25 * Mtmp36;
  M[35] += -1.0 / 6.0 * Mtmp12 * Mtmp9;
  M[36] += -Mtmp28 * Mtmp36;
  M[37] += -Mtmp37 * Mtmp38;
  M[38] += -Mtmp27 * z;
  M[39] += -Mtmp10 * Mtmp26 * Mtmp8;
  M[40] += -Mtmp38 * Mtmp39;
  M[41] += -Mtmp30 * Mtmp40;
  M[42] += -1.0 / 6.0 * Mtmp18 * Mtmp5;
  M[43] += -Mtmp3 * Mtmp32 * Mtmp6;
  M[44] += -1.0 / 6.0 * Mtmp20 * Mtmp4;
  M[45] += -Mtmp33 * Mtmp40;
  M[46] += -Mtmp34 * (y * y * y * y * y);
  M[47] += -Mtmp31 * z;
  M[48] += -Mtmp28 * Mtmp37;
  M[49] += -Mtmp25 * Mtmp39;
  M[50] += -Mtmp33 * Mtmp35;
  M[51] += -Mtmp34 * (z * z * z * z * z);
}
void field_m2_M2M_5(double x, double y, double z, double* M, double* Ms) {
  double Mstmp0  = x * M[0];
  double Mstmp1  = x * M[1];
  double Mstmp2  = y * M[0];
  double Mstmp3  = x * M[2];
  double Mstmp4  = z * M[0];
  double Mstmp5  = x * M[3];
  double Mstmp6  = y * M[1];
  double Mstmp7  = x * M[4];
  double Mstmp8  = y * M[2];
  double Mstmp9  = z * M[1];
  double Mstmp10 = x * M[5];
  double Mstmp11 = z * M[2];
  double Mstmp12 = y * M[3];
  double Mstmp13 = y * M[4];
  double Mstmp14 = z * M[3];
  double Mstmp15 = y * M[5];
  double Mstmp16 = z * M[4];
  double Mstmp17 = z * M[5];
  double Mstmp18 = x * M[6];
  double Mstmp19 = (1.0 / 2.0) * (x * x);
  double Mstmp20 = x * M[7];
  double Mstmp21 = y * M[6];
  double Mstmp22 = Mstmp0 * y;
  double Mstmp23 = x * M[8];
  double Mstmp24 = x * M[9];
  double Mstmp25 = y * M[7];
  double Mstmp26 = Mstmp1 * y;
  double Mstmp27 = (y * y);
  double Mstmp28 = (1.0 / 2.0) * M[0];
  double Mstmp29 = x * M[10];
  double Mstmp30 = y * M[8];
  double Mstmp31 = Mstmp3 * y;
  double Mstmp32 = x * M[11];
  double Mstmp33 = (z * z);
  double Mstmp34 = x * M[12];
  double Mstmp35 = y * M[9];
  double Mstmp36 = Mstmp5 * y;
  double Mstmp37 = (1.0 / 2.0) * Mstmp27;
  double Mstmp38 = x * M[13];
  double Mstmp39 = y * M[10];
  double Mstmp40 = Mstmp7 * y;
  double Mstmp41 = x * M[14];
  double Mstmp42 = y * M[11];
  double Mstmp43 = Mstmp10 * y;
  double Mstmp44 = (1.0 / 2.0) * Mstmp33;
  double Mstmp45 = x * M[15];
  double Mstmp46 = y * M[12];
  double Mstmp47 = y * M[13];
  double Mstmp48 = y * M[14];
  double Mstmp49 = y * M[15];
  double Mstmp50 = (1.0 / 6.0) * (x * x * x);
  double Mstmp51 = (y * y * y);
  double Mstmp52 = (1.0 / 6.0) * M[0];
  double Mstmp53 = (z * z * z);
  double Mstmp54 = (1.0 / 6.0) * Mstmp51;
  double Mstmp55 = (1.0 / 6.0) * Mstmp53;
#pragma omp atomic
  Ms[0] += M[0];
#pragma omp atomic
  Ms[1] += M[1];
#pragma omp atomic
  Ms[2] += M[2];
#pragma omp atomic
  Ms[3] += M[3];
#pragma omp atomic
  Ms[4] += M[4];
#pragma omp atomic
  Ms[5] += M[5];
#pragma omp atomic
  Ms[6] += Mstmp0 + M[6];
#pragma omp atomic
  Ms[7] += Mstmp1 + Mstmp2 + M[7];
#pragma omp atomic
  Ms[8] += Mstmp3 + Mstmp4 + M[8];
#pragma omp atomic
  Ms[9] += Mstmp5 + Mstmp6 + M[9];
#pragma omp atomic
  Ms[10] += Mstmp7 + Mstmp8 + Mstmp9 + M[10];
#pragma omp atomic
  Ms[11] += Mstmp10 + Mstmp11 + M[11];
#pragma omp atomic
  Ms[12] += Mstmp12 + M[12];
#pragma omp atomic
  Ms[13] += Mstmp13 + Mstmp14 + M[13];
#pragma omp atomic
  Ms[14] += Mstmp15 + Mstmp16 + M[14];
#pragma omp atomic
  Ms[15] += Mstmp17 + M[15];
#pragma omp atomic
  Ms[16] += Mstmp18 + Mstmp19 * M[0] + M[16];
#pragma omp atomic
  Ms[17] += Mstmp19 * M[1] + Mstmp20 + Mstmp21 + Mstmp22 + M[17];
#pragma omp atomic
  Ms[18] += Mstmp0 * z + Mstmp19 * M[2] + Mstmp23 + z * M[6] + M[18];
#pragma omp atomic
  Ms[19] += Mstmp19 * M[3] + Mstmp24 + Mstmp25 + Mstmp26 + Mstmp27 * Mstmp28 + M[19];
#pragma omp atomic
  Ms[20] += Mstmp1 * z + Mstmp19 * M[4] + Mstmp2 * z + Mstmp29 + Mstmp30 + Mstmp31 +
            z * M[7] + M[20];
#pragma omp atomic
  Ms[21] += Mstmp19 * M[5] + Mstmp28 * Mstmp33 + Mstmp3 * z + Mstmp32 + z * M[8] + M[21];
#pragma omp atomic
  Ms[22] += Mstmp34 + Mstmp35 + Mstmp36 + Mstmp37 * M[1] + M[22];
#pragma omp atomic
  Ms[23] += Mstmp37 * M[2] + Mstmp38 + Mstmp39 + Mstmp40 + Mstmp5 * z + Mstmp6 * z +
            z * M[9] + M[23];
#pragma omp atomic
  Ms[24] += Mstmp41 + Mstmp42 + Mstmp43 + Mstmp44 * M[1] + Mstmp7 * z + Mstmp8 * z +
            z * M[10] + M[24];
#pragma omp atomic
  Ms[25] += Mstmp10 * z + Mstmp44 * M[2] + Mstmp45 + z * M[11] + M[25];
#pragma omp atomic
  Ms[26] += Mstmp37 * M[3] + Mstmp46 + M[26];
#pragma omp atomic
  Ms[27] += Mstmp12 * z + Mstmp37 * M[4] + Mstmp47 + z * M[12] + M[27];
#pragma omp atomic
  Ms[28] += Mstmp13 * z + Mstmp37 * M[5] + Mstmp44 * M[3] + Mstmp48 + z * M[13] + M[28];
#pragma omp atomic
  Ms[29] += Mstmp15 * z + Mstmp44 * M[4] + Mstmp49 + z * M[14] + M[29];
#pragma omp atomic
  Ms[30] += Mstmp44 * M[5] + z * M[15] + M[30];
#pragma omp atomic
  Ms[31] += Mstmp19 * M[6] + Mstmp50 * M[0] + x * M[16] + M[31];
#pragma omp atomic
  Ms[32] += Mstmp18 * y + Mstmp19 * Mstmp2 + Mstmp19 * M[7] + Mstmp50 * M[1] + x * M[17] +
            y * M[16] + M[32];
#pragma omp atomic
  Ms[33] += Mstmp18 * z + Mstmp19 * Mstmp4 + Mstmp19 * M[8] + Mstmp50 * M[2] + x * M[18] +
            z * M[16] + M[33];
#pragma omp atomic
  Ms[34] += Mstmp0 * Mstmp37 + Mstmp19 * Mstmp6 + Mstmp19 * M[9] + Mstmp20 * y +
            Mstmp37 * M[6] + Mstmp50 * M[3] + x * M[19] + y * M[17] + M[34];
#pragma omp atomic
  Ms[35] += Mstmp19 * Mstmp8 + Mstmp19 * Mstmp9 + Mstmp19 * M[10] + Mstmp20 * z +
            Mstmp21 * z + Mstmp22 * z + Mstmp23 * y + Mstmp50 * M[4] + x * M[20] +
            y * M[18] + z * M[17] + M[35];
#pragma omp atomic
  Ms[36] += Mstmp0 * Mstmp44 + Mstmp11 * Mstmp19 + Mstmp19 * M[11] + Mstmp23 * z +
            Mstmp44 * M[6] + Mstmp50 * M[5] + x * M[21] + z * M[18] + M[36];
#pragma omp atomic
  Ms[37] += Mstmp1 * Mstmp37 + Mstmp12 * Mstmp19 + Mstmp19 * M[12] + Mstmp24 * y +
            Mstmp37 * M[7] + Mstmp51 * Mstmp52 + x * M[22] + y * M[19] + M[37];
#pragma omp atomic
  Ms[38] += Mstmp13 * Mstmp19 + Mstmp14 * Mstmp19 + Mstmp19 * M[13] + Mstmp24 * z +
            Mstmp25 * z + Mstmp26 * z + Mstmp29 * y + Mstmp3 * Mstmp37 +
            Mstmp37 * Mstmp4 + Mstmp37 * M[8] + x * M[23] + y * M[20] + z * M[19] + M[38];
#pragma omp atomic
  Ms[39] += Mstmp1 * Mstmp44 + Mstmp15 * Mstmp19 + Mstmp16 * Mstmp19 + Mstmp19 * M[14] +
            Mstmp2 * Mstmp44 + Mstmp29 * z + Mstmp30 * z + Mstmp31 * z + Mstmp32 * y +
            Mstmp44 * M[7] + x * M[24] + y * M[21] + z * M[20] + M[39];
#pragma omp atomic
  Ms[40] += Mstmp17 * Mstmp19 + Mstmp19 * M[15] + Mstmp3 * Mstmp44 + Mstmp32 * z +
            Mstmp44 * M[8] + Mstmp52 * Mstmp53 + x * M[25] + z * M[21] + M[40];
#pragma omp atomic
  Ms[41] += Mstmp34 * y + Mstmp37 * Mstmp5 + Mstmp37 * M[9] + Mstmp54 * M[1] + x * M[26] +
            y * M[22] + M[41];
#pragma omp atomic
  Ms[42] += Mstmp34 * z + Mstmp35 * z + Mstmp36 * z + Mstmp37 * Mstmp7 +
            Mstmp37 * Mstmp9 + Mstmp37 * M[10] + Mstmp38 * y + Mstmp54 * M[2] +
            x * M[27] + y * M[23] + z * M[22] + M[42];
#pragma omp atomic
  Ms[43] += Mstmp10 * Mstmp37 + Mstmp11 * Mstmp37 + Mstmp37 * M[11] + Mstmp38 * z +
            Mstmp39 * z + Mstmp40 * z + Mstmp41 * y + Mstmp44 * Mstmp5 +
            Mstmp44 * Mstmp6 + Mstmp44 * M[9] + x * M[28] + y * M[24] + z * M[23] + M[43];
#pragma omp atomic
  Ms[44] += Mstmp41 * z + Mstmp42 * z + Mstmp43 * z + Mstmp44 * Mstmp7 +
            Mstmp44 * Mstmp8 + Mstmp44 * M[10] + Mstmp45 * y + Mstmp55 * M[1] +
            x * M[29] + y * M[25] + z * M[24] + M[44];
#pragma omp atomic
  Ms[45] += Mstmp10 * Mstmp44 + Mstmp44 * M[11] + Mstmp45 * z + Mstmp55 * M[2] +
            x * M[30] + z * M[25] + M[45];
#pragma omp atomic
  Ms[46] += Mstmp37 * M[12] + Mstmp54 * M[3] + y * M[26] + M[46];
#pragma omp atomic
  Ms[47] += Mstmp14 * Mstmp37 + Mstmp37 * M[13] + Mstmp46 * z + Mstmp54 * M[4] +
            y * M[27] + z * M[26] + M[47];
#pragma omp atomic
  Ms[48] += Mstmp12 * Mstmp44 + Mstmp16 * Mstmp37 + Mstmp37 * M[14] + Mstmp44 * M[12] +
            Mstmp47 * z + Mstmp54 * M[5] + y * M[28] + z * M[27] + M[48];
#pragma omp atomic
  Ms[49] += Mstmp13 * Mstmp44 + Mstmp17 * Mstmp37 + Mstmp37 * M[15] + Mstmp44 * M[13] +
            Mstmp48 * z + Mstmp55 * M[3] + y * M[29] + z * M[28] + M[49];
#pragma omp atomic
  Ms[50] += Mstmp15 * Mstmp44 + Mstmp44 * M[14] + Mstmp49 * z + Mstmp55 * M[4] +
            y * M[30] + z * M[29] + M[50];
#pragma omp atomic
  Ms[51] += Mstmp44 * M[15] + Mstmp55 * M[5] + z * M[30] + M[51];
}

void field_m2_M2L_5(double x, double y, double z, double* M, double* L) {
  double R = sqrt(x * x + y * y + z * z);
  double D[52];
  double Dtmp0  = -1.0 * pow(R, -3.0);
  double Dtmp1  = (x * x);
  double Dtmp2  = pow(R, -5.0);
  double Dtmp3  = 3.0 * Dtmp2;
  double Dtmp4  = Dtmp3 * x;
  double Dtmp5  = (y * y);
  double Dtmp6  = y * z;
  double Dtmp7  = 9.0 * Dtmp2;
  double Dtmp8  = -Dtmp7;
  double Dtmp9  = pow(R, -7.0);
  double Dtmp10 = 15.0 * Dtmp9;
  double Dtmp11 = Dtmp1 * Dtmp10;
  double Dtmp12 = -Dtmp3;
  double Dtmp13 = Dtmp11 + Dtmp12;
  double Dtmp14 = Dtmp10 * Dtmp5;
  double Dtmp15 = Dtmp12 + Dtmp14;
  double Dtmp16 = 1.0 * x;
  double Dtmp17 = Dtmp6 * x;
  double Dtmp18 = (x * x * x * x);
  double Dtmp19 = pow(R, -9.0);
  double Dtmp20 = 105.0 * Dtmp19;
  double Dtmp21 = 90.0 * Dtmp9;
  double Dtmp22 = 45.0 * Dtmp9;
  double Dtmp23 = -Dtmp22;
  double Dtmp24 = Dtmp1 * Dtmp20;
  double Dtmp25 = x * (Dtmp23 + Dtmp24);
  double Dtmp26 = -Dtmp10;
  double Dtmp27 = Dtmp20 * Dtmp5;
  double Dtmp28 = Dtmp23 + Dtmp27;
  double Dtmp29 = (y * y * y * y);
  double Dtmp30 = 225.0 * Dtmp9;
  double Dtmp31 = 945.0 * pow(R, -11.0);
  double Dtmp32 = Dtmp18 * Dtmp31;
  double Dtmp33 = Dtmp1 * Dtmp19;
  double Dtmp34 = Dtmp22 + Dtmp32 - 630.0 * Dtmp33;
  double Dtmp35 = -Dtmp24;
  double Dtmp36 = 315.0 * Dtmp19;
  double Dtmp37 = Dtmp1 * Dtmp31;
  double Dtmp38 = Dtmp37 * Dtmp5;
  double Dtmp39 = Dtmp22 + Dtmp38;
  double Dtmp40 = -Dtmp36;
  double Dtmp41 = -Dtmp27;
  double Dtmp42 = Dtmp29 * Dtmp31;
  double Dtmp43 = Dtmp19 * Dtmp5;
  double Dtmp44 = Dtmp22 + Dtmp42 - 630.0 * Dtmp43;
  D[0]          = Dtmp0 + Dtmp1 * Dtmp3;
  D[1]          = Dtmp4 * y;
  D[2]          = Dtmp4 * z;
  D[3]          = Dtmp0 + Dtmp3 * Dtmp5;
  D[4]          = Dtmp3 * Dtmp6;
  D[5]          = -D[0] - D[3];
  D[6]          = -x * (Dtmp11 + Dtmp8);
  D[7]          = -Dtmp13 * y;
  D[8]          = -Dtmp13 * z;
  D[9]          = -Dtmp15 * Dtmp16;
  D[10]         = -Dtmp10 * Dtmp17;
  D[11]         = -D[6] - D[9];
  D[12]         = -y * (Dtmp14 + Dtmp8);
  D[13]         = -Dtmp15 * z;
  D[14]         = -D[7] - D[12];
  D[15]         = -D[8] - D[13];
  D[16]         = -Dtmp1 * Dtmp21 + Dtmp18 * Dtmp20 + Dtmp7;
  D[17]         = Dtmp25 * y;
  D[18]         = Dtmp25 * z;
  D[19]         = -Dtmp11 - Dtmp14 + Dtmp24 * Dtmp5 + Dtmp3;
  D[20]         = Dtmp6 * (Dtmp24 + Dtmp26);
  D[21]         = -D[16] - D[19];
  D[22]         = Dtmp16 * Dtmp28 * y;
  D[23]         = Dtmp16 * z * (Dtmp26 + Dtmp27);
  D[24]         = -D[17] - D[22];
  D[25]         = -D[18] - D[23];
  D[26]         = Dtmp20 * Dtmp29 - Dtmp21 * Dtmp5 + Dtmp7;
  D[27]         = Dtmp28 * Dtmp6;
  D[28]         = -D[19] - D[26];
  D[29]         = -D[20] - D[27];
  D[30]         = -D[21] - D[28];
  D[31]         = -x * (Dtmp30 + Dtmp32 - 1050.0 * Dtmp33);
  D[32]         = -Dtmp34 * y;
  D[33]         = -Dtmp34 * z;
  D[34]         = -x * (Dtmp35 - Dtmp36 * Dtmp5 + Dtmp39);
  D[35]         = -Dtmp17 * (Dtmp37 + Dtmp40);
  D[36]         = -D[31] - D[34];
  D[37]         = -y * (-Dtmp1 * Dtmp36 + Dtmp39 + Dtmp41);
  D[38]         = -z * (Dtmp10 + Dtmp35 + Dtmp38 + Dtmp41);
  D[39]         = -D[32] - D[37];
  D[40]         = -D[33] - D[38];
  D[41]         = -Dtmp16 * Dtmp44;
  D[42]         = -Dtmp16 * Dtmp6 * (Dtmp31 * Dtmp5 + Dtmp40);
  D[43]         = -D[34] - D[41];
  D[44]         = -D[35] - D[42];
  D[45]         = -D[36] - D[43];
  D[46]         = -y * (Dtmp30 + Dtmp42 - 1050.0 * Dtmp43);
  D[47]         = -Dtmp44 * z;
  D[48]         = -D[37] - D[46];
  D[49]         = -D[38] - D[47];
  D[50]         = -D[39] - D[48];
  D[51]         = -D[40] - D[49];
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
          D[50] * M[50] + D[51] * M[51];
#pragma omp atomic
  L[1] += D[6] * M[0] + D[7] * M[1] + D[8] * M[2] + D[9] * M[3] + D[10] * M[4] +
          D[11] * M[5] + D[16] * M[6] + D[17] * M[7] + D[18] * M[8] + D[19] * M[9] +
          D[20] * M[10] + D[21] * M[11] + D[22] * M[12] + D[23] * M[13] + D[24] * M[14] +
          D[25] * M[15] + D[31] * M[16] + D[32] * M[17] + D[33] * M[18] + D[34] * M[19] +
          D[35] * M[20] + D[36] * M[21] + D[37] * M[22] + D[38] * M[23] + D[39] * M[24] +
          D[40] * M[25] + D[41] * M[26] + D[42] * M[27] + D[43] * M[28] + D[44] * M[29] +
          D[45] * M[30];
#pragma omp atomic
  L[2] += D[7] * M[0] + D[9] * M[1] + D[10] * M[2] + D[12] * M[3] + D[13] * M[4] +
          D[14] * M[5] + D[17] * M[6] + D[19] * M[7] + D[20] * M[8] + D[22] * M[9] +
          D[23] * M[10] + D[24] * M[11] + D[26] * M[12] + D[27] * M[13] + D[28] * M[14] +
          D[29] * M[15] + D[32] * M[16] + D[34] * M[17] + D[35] * M[18] + D[37] * M[19] +
          D[38] * M[20] + D[39] * M[21] + D[41] * M[22] + D[42] * M[23] + D[43] * M[24] +
          D[44] * M[25] + D[46] * M[26] + D[47] * M[27] + D[48] * M[28] + D[49] * M[29] +
          D[50] * M[30];
#pragma omp atomic
  L[3] += D[8] * M[0] + D[10] * M[1] + D[11] * M[2] + D[13] * M[3] + D[14] * M[4] +
          D[15] * M[5] + D[18] * M[6] + D[20] * M[7] + D[21] * M[8] + D[23] * M[9] +
          D[24] * M[10] + D[25] * M[11] + D[27] * M[12] + D[28] * M[13] + D[29] * M[14] +
          D[30] * M[15] + D[33] * M[16] + D[35] * M[17] + D[36] * M[18] + D[38] * M[19] +
          D[39] * M[20] + D[40] * M[21] + D[42] * M[22] + D[43] * M[23] + D[44] * M[24] +
          D[45] * M[25] + D[47] * M[26] + D[48] * M[27] + D[49] * M[28] + D[50] * M[29] +
          D[51] * M[30];
#pragma omp atomic
  L[4] += D[16] * M[0] + D[17] * M[1] + D[18] * M[2] + D[19] * M[3] + D[20] * M[4] +
          D[21] * M[5] + D[31] * M[6] + D[32] * M[7] + D[33] * M[8] + D[34] * M[9] +
          D[35] * M[10] + D[36] * M[11] + D[37] * M[12] + D[38] * M[13] + D[39] * M[14] +
          D[40] * M[15];
#pragma omp atomic
  L[5] += D[17] * M[0] + D[19] * M[1] + D[20] * M[2] + D[22] * M[3] + D[23] * M[4] +
          D[24] * M[5] + D[32] * M[6] + D[34] * M[7] + D[35] * M[8] + D[37] * M[9] +
          D[38] * M[10] + D[39] * M[11] + D[41] * M[12] + D[42] * M[13] + D[43] * M[14] +
          D[44] * M[15];
#pragma omp atomic
  L[6] += D[18] * M[0] + D[20] * M[1] + D[21] * M[2] + D[23] * M[3] + D[24] * M[4] +
          D[25] * M[5] + D[33] * M[6] + D[35] * M[7] + D[36] * M[8] + D[38] * M[9] +
          D[39] * M[10] + D[40] * M[11] + D[42] * M[12] + D[43] * M[13] + D[44] * M[14] +
          D[45] * M[15];
#pragma omp atomic
  L[7] += D[19] * M[0] + D[22] * M[1] + D[23] * M[2] + D[26] * M[3] + D[27] * M[4] +
          D[28] * M[5] + D[34] * M[6] + D[37] * M[7] + D[38] * M[8] + D[41] * M[9] +
          D[42] * M[10] + D[43] * M[11] + D[46] * M[12] + D[47] * M[13] + D[48] * M[14] +
          D[49] * M[15];
#pragma omp atomic
  L[8] += D[20] * M[0] + D[23] * M[1] + D[24] * M[2] + D[27] * M[3] + D[28] * M[4] +
          D[29] * M[5] + D[35] * M[6] + D[38] * M[7] + D[39] * M[8] + D[42] * M[9] +
          D[43] * M[10] + D[44] * M[11] + D[47] * M[12] + D[48] * M[13] + D[49] * M[14] +
          D[50] * M[15];
#pragma omp atomic
  L[9] += D[21] * M[0] + D[24] * M[1] + D[25] * M[2] + D[28] * M[3] + D[29] * M[4] +
          D[30] * M[5] + D[36] * M[6] + D[39] * M[7] + D[40] * M[8] + D[43] * M[9] +
          D[44] * M[10] + D[45] * M[11] + D[48] * M[12] + D[49] * M[13] + D[50] * M[14] +
          D[51] * M[15];
#pragma omp atomic
  L[10] += D[31] * M[0] + D[32] * M[1] + D[33] * M[2] + D[34] * M[3] + D[35] * M[4] +
           D[36] * M[5];
#pragma omp atomic
  L[11] += D[32] * M[0] + D[34] * M[1] + D[35] * M[2] + D[37] * M[3] + D[38] * M[4] +
           D[39] * M[5];
#pragma omp atomic
  L[12] += D[33] * M[0] + D[35] * M[1] + D[36] * M[2] + D[38] * M[3] + D[39] * M[4] +
           D[40] * M[5];
#pragma omp atomic
  L[13] += D[34] * M[0] + D[37] * M[1] + D[38] * M[2] + D[41] * M[3] + D[42] * M[4] +
           D[43] * M[5];
#pragma omp atomic
  L[14] += D[35] * M[0] + D[38] * M[1] + D[39] * M[2] + D[42] * M[3] + D[43] * M[4] +
           D[44] * M[5];
#pragma omp atomic
  L[15] += D[36] * M[0] + D[39] * M[1] + D[40] * M[2] + D[43] * M[3] + D[44] * M[4] +
           D[45] * M[5];
#pragma omp atomic
  L[16] += D[37] * M[0] + D[41] * M[1] + D[42] * M[2] + D[46] * M[3] + D[47] * M[4] +
           D[48] * M[5];
#pragma omp atomic
  L[17] += D[38] * M[0] + D[42] * M[1] + D[43] * M[2] + D[47] * M[3] + D[48] * M[4] +
           D[49] * M[5];
#pragma omp atomic
  L[18] += D[39] * M[0] + D[43] * M[1] + D[44] * M[2] + D[48] * M[3] + D[49] * M[4] +
           D[50] * M[5];
#pragma omp atomic
  L[19] += D[40] * M[0] + D[44] * M[1] + D[45] * M[2] + D[49] * M[3] + D[50] * M[4] +
           D[51] * M[5];
}

void field_m2_L2L_5(double x, double y, double z, double* L, double* Ls) {
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

void field_m2_L2P_5(double x, double y, double z, double* L, double* F) {
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

void field_m2_M2P_5(double x, double y, double z, double* M, double* F) {
  double R       = sqrt(x * x + y * y + z * z);
  double Ftmp0   = pow(R, -5.0);
  double Ftmp1   = 3.0 * Ftmp0;
  double Ftmp2   = Ftmp1 * M[1];
  double Ftmp3   = Ftmp1 * z;
  double Ftmp4   = y * z;
  double Ftmp5   = pow(R, -7.0);
  double Ftmp6   = 15.0 * Ftmp5;
  double Ftmp7   = Ftmp6 * M[10];
  double Ftmp8   = Ftmp4 * x;
  double Ftmp9   = Ftmp6 * Ftmp8;
  double Ftmp10  = (x * x);
  double Ftmp11  = Ftmp10 * Ftmp6;
  double Ftmp12  = pow(R, -9.0);
  double Ftmp13  = 105.0 * Ftmp12;
  double Ftmp14  = Ftmp10 * Ftmp13;
  double Ftmp15  = -9.0 * Ftmp0;
  double Ftmp16  = Ftmp11 + Ftmp15;
  double Ftmp17  = -Ftmp1;
  double Ftmp18  = (y * y);
  double Ftmp19  = Ftmp18 * Ftmp6;
  double Ftmp20  = Ftmp17 + Ftmp19;
  double Ftmp21  = (z * z);
  double Ftmp22  = Ftmp21 * Ftmp6;
  double Ftmp23  = Ftmp17 + Ftmp22;
  double Ftmp24  = Ftmp20 * M[3];
  double Ftmp25  = Ftmp23 * M[5];
  double Ftmp26  = 45.0 * Ftmp5;
  double Ftmp27  = -Ftmp26;
  double Ftmp28  = Ftmp14 + Ftmp27;
  double Ftmp29  = Ftmp28 * M[17];
  double Ftmp30  = Ftmp13 * Ftmp18;
  double Ftmp31  = Ftmp27 + Ftmp30;
  double Ftmp32  = Ftmp31 * y;
  double Ftmp33  = 1.0 * M[22];
  double Ftmp34  = 3.0 * y;
  double Ftmp35  = 35.0 * Ftmp12;
  double Ftmp36  = (Ftmp21 * Ftmp35 - 5.0 * Ftmp5) * M[24];
  double Ftmp37  = Ftmp28 * M[18];
  double Ftmp38  = 1.0 * z;
  double Ftmp39  = -Ftmp6;
  double Ftmp40  = Ftmp30 + Ftmp39;
  double Ftmp41  = Ftmp40 * M[23];
  double Ftmp42  = Ftmp13 * Ftmp21;
  double Ftmp43  = Ftmp27 + Ftmp42;
  double Ftmp44  = Ftmp38 * Ftmp43;
  double Ftmp45  = 315.0 * Ftmp12;
  double Ftmp46  = -Ftmp45;
  double Ftmp47  = pow(R, -11.0);
  double Ftmp48  = 945.0 * Ftmp47;
  double Ftmp49  = Ftmp10 * Ftmp48;
  double Ftmp50  = Ftmp46 + Ftmp49;
  double Ftmp51  = Ftmp50 * M[35];
  double Ftmp52  = x * y;
  double Ftmp53  = Ftmp28 * Ftmp52;
  double Ftmp54  = Ftmp32 * M[12];
  double Ftmp55  = Ftmp39 + Ftmp42;
  double Ftmp56  = Ftmp55 * M[14];
  double Ftmp57  = 1.0 * Ftmp52;
  double Ftmp58  = x * z;
  double Ftmp59  = Ftmp28 * Ftmp58;
  double Ftmp60  = Ftmp40 * M[13];
  double Ftmp61  = Ftmp43 * M[15];
  double Ftmp62  = Ftmp18 * Ftmp48;
  double Ftmp63  = Ftmp46 + Ftmp62;
  double Ftmp64  = Ftmp63 * y;
  double Ftmp65  = Ftmp38 * M[42];
  double Ftmp66  = -Ftmp13;
  double Ftmp67  = Ftmp21 * Ftmp47;
  double Ftmp68  = 315.0 * Ftmp67;
  double Ftmp69  = Ftmp66 + Ftmp68;
  double Ftmp70  = Ftmp69 * M[44];
  double Ftmp71  = Ftmp34 * Ftmp70;
  double Ftmp72  = Ftmp50 * x;
  double Ftmp73  = Ftmp63 * M[27];
  double Ftmp74  = -75.0 * Ftmp5;
  double Ftmp75  = 1.0 * Ftmp10;
  double Ftmp76  = Ftmp40 * M[9];
  double Ftmp77  = Ftmp55 * M[11];
  double Ftmp78  = 525.0 * Ftmp12;
  double Ftmp79  = -Ftmp78;
  double Ftmp80  = Ftmp10 * (Ftmp49 + Ftmp79);
  double Ftmp81  = Ftmp21 * Ftmp48;
  double Ftmp82  = Ftmp46 + Ftmp81;
  double Ftmp83  = Ftmp38 * Ftmp82 * M[29];
  double Ftmp84  = (-Ftmp35 + Ftmp68) * M[24];
  double Ftmp85  = Ftmp10 * Ftmp34;
  double Ftmp86  = Ftmp10 * Ftmp38;
  double Ftmp87  = (Ftmp62 + Ftmp66) * M[23];
  double Ftmp88  = Ftmp82 * M[25];
  double Ftmp89  = 4725.0 * Ftmp47;
  double Ftmp90  = -Ftmp89;
  double Ftmp91  = pow(R, -13.0);
  double Ftmp92  = 10395.0 * Ftmp91;
  double Ftmp93  = Ftmp10 * Ftmp92;
  double Ftmp94  = 2835.0 * Ftmp47;
  double Ftmp95  = -Ftmp94;
  double Ftmp96  = Ftmp18 * Ftmp92;
  double Ftmp97  = Ftmp95 + Ftmp96;
  double Ftmp98  = 3465.0 * Ftmp21 * Ftmp91;
  double Ftmp99  = (-Ftmp48 + Ftmp98) * M[44];
  double Ftmp100 = 225.0 * Ftmp5;
  double Ftmp101 = (x * x * x * x);
  double Ftmp102 = Ftmp101 * Ftmp48;
  double Ftmp103 = 1050.0 * Ftmp12;
  double Ftmp104 = -Ftmp10 * Ftmp103 + Ftmp100 + Ftmp102;
  double Ftmp105 = (y * y * y * y);
  double Ftmp106 = Ftmp105 * Ftmp48;
  double Ftmp107 = 630.0 * Ftmp12;
  double Ftmp108 = Ftmp106 - Ftmp107 * Ftmp18 + Ftmp26;
  double Ftmp109 = (z * z * z * z);
  double Ftmp110 = Ftmp109 * Ftmp48;
  double Ftmp111 = -Ftmp107 * Ftmp21 + Ftmp110 + Ftmp26;
  double Ftmp112 = Ftmp108 * M[26];
  double Ftmp113 = Ftmp111 * M[30];
  double Ftmp114 = 1575.0 * Ftmp12;
  double Ftmp115 = Ftmp101 * Ftmp92;
  double Ftmp116 = Ftmp10 * Ftmp47;
  double Ftmp117 = Ftmp114 + Ftmp115 - 9450.0 * Ftmp116;
  double Ftmp118 = Ftmp117 * Ftmp52;
  double Ftmp119 = Ftmp105 * Ftmp92;
  double Ftmp120 = 9450.0 * Ftmp47;
  double Ftmp121 = Ftmp114 + Ftmp119 - Ftmp120 * Ftmp18;
  double Ftmp122 = Ftmp121 * M[46];
  double Ftmp123 = Ftmp109 * Ftmp92;
  double Ftmp124 = 5670.0 * Ftmp47;
  double Ftmp125 = Ftmp123 - Ftmp124 * Ftmp21 + Ftmp45;
  double Ftmp126 = Ftmp125 * M[50];
  double Ftmp127 = Ftmp117 * Ftmp58;
  double Ftmp128 = Ftmp119 - Ftmp124 * Ftmp18 + Ftmp45;
  double Ftmp129 = Ftmp128 * M[47];
  double Ftmp130 = Ftmp114 - Ftmp120 * Ftmp21 + Ftmp123;
  double Ftmp131 = Ftmp130 * M[51];
  double Ftmp132 = 3675.0 * Ftmp12;
  double Ftmp133 = Ftmp128 * M[41];
  double Ftmp134 = Ftmp125 * M[45];
  double Ftmp135 = -Ftmp18 * Ftmp45;
  double Ftmp136 = Ftmp18 * Ftmp49;
  double Ftmp137 = -Ftmp14;
  double Ftmp138 = Ftmp137 + Ftmp26;
  double Ftmp139 = Ftmp135 + Ftmp136 + Ftmp138;
  double Ftmp140 = -Ftmp21 * Ftmp45;
  double Ftmp141 = Ftmp21 * Ftmp49;
  double Ftmp142 = Ftmp138 + Ftmp140 + Ftmp141;
  double Ftmp143 = -Ftmp42;
  double Ftmp144 = Ftmp143 + Ftmp6;
  double Ftmp145 = -Ftmp30;
  double Ftmp146 = Ftmp21 * Ftmp62;
  double Ftmp147 = Ftmp145 + Ftmp146;
  double Ftmp148 = Ftmp144 + Ftmp147;
  double Ftmp149 = Ftmp18 * Ftmp93;
  double Ftmp150 = -Ftmp18 * Ftmp94;
  double Ftmp151 = Ftmp149 + Ftmp150;
  double Ftmp152 = 945.0 * Ftmp12;
  double Ftmp153 = -Ftmp10 * Ftmp94;
  double Ftmp154 = Ftmp152 + Ftmp153;
  double Ftmp155 = Ftmp52 * (Ftmp151 + Ftmp154);
  double Ftmp156 = -Ftmp21 * Ftmp94;
  double Ftmp157 = Ftmp156 + Ftmp45;
  double Ftmp158 = -Ftmp49;
  double Ftmp159 = Ftmp21 * Ftmp93;
  double Ftmp160 = Ftmp158 + Ftmp159;
  double Ftmp161 = Ftmp52 * (Ftmp157 + Ftmp160);
  double Ftmp162 = -Ftmp62;
  double Ftmp163 = Ftmp21 * Ftmp96;
  double Ftmp164 = Ftmp162 + Ftmp163;
  double Ftmp165 = Ftmp52 * (Ftmp157 + Ftmp164);
  double Ftmp166 = Ftmp58 * (Ftmp151 + Ftmp158 + Ftmp45);
  double Ftmp167 = Ftmp58 * (Ftmp154 + Ftmp156 + Ftmp159);
  double Ftmp168 = -Ftmp81;
  double Ftmp169 = Ftmp168 + Ftmp45;
  double Ftmp170 = Ftmp150 + Ftmp163;
  double Ftmp171 = Ftmp169 + Ftmp170;
  double Ftmp172 = -Ftmp18 * Ftmp89;
  double Ftmp173 = Ftmp158 + Ftmp78;
  double Ftmp174 = -Ftmp21 * Ftmp89;
  double Ftmp175 = Ftmp13 + Ftmp168;
  double Ftmp176 = x * M[2];
  double Ftmp177 = Ftmp11 + Ftmp17;
  double Ftmp178 = Ftmp15 + Ftmp19;
  double Ftmp179 = Ftmp177 * M[0];
  double Ftmp180 = Ftmp33 * x;
  double Ftmp181 = 3.0 * x;
  double Ftmp182 = Ftmp14 + Ftmp39;
  double Ftmp183 = Ftmp182 * M[20];
  double Ftmp184 = Ftmp31 * M[27];
  double Ftmp185 = 1.0 * x;
  double Ftmp186 = Ftmp65 * x;
  double Ftmp187 = 3.0 * Ftmp58;
  double Ftmp188 = Ftmp182 * M[8];
  double Ftmp189 = Ftmp72 * M[18];
  double Ftmp190 = Ftmp182 * M[7];
  double Ftmp191 = 1.0 * Ftmp18;
  double Ftmp192 = Ftmp72 * M[17];
  double Ftmp193 = Ftmp18 * z;
  double Ftmp194 = (Ftmp49 + Ftmp66) * M[20];
  double Ftmp195 = Ftmp62 + Ftmp79;
  double Ftmp196 = Ftmp38 * Ftmp52;
  double Ftmp197 = (Ftmp93 + Ftmp95) * M[35];
  double Ftmp198 = -Ftmp10 * Ftmp107 + Ftmp102 + Ftmp26;
  double Ftmp199 = Ftmp100 - Ftmp103 * Ftmp18 + Ftmp106;
  double Ftmp200 = Ftmp198 * M[16];
  double Ftmp201 = Ftmp115 - 5670.0 * Ftmp116 + Ftmp45;
  double Ftmp202 = Ftmp201 * M[33];
  double Ftmp203 = Ftmp201 * M[32];
  double Ftmp204 = Ftmp136 + Ftmp145;
  double Ftmp205 = -Ftmp10 * Ftmp45 + Ftmp26;
  double Ftmp206 = Ftmp204 + Ftmp205;
  double Ftmp207 = Ftmp137 + Ftmp141 + Ftmp144;
  double Ftmp208 = Ftmp140 + Ftmp147 + Ftmp26;
  double Ftmp209 = Ftmp149 + Ftmp162;
  double Ftmp210 = Ftmp4 * (Ftmp153 + Ftmp209 + Ftmp45);
  double Ftmp211 = Ftmp4 * (Ftmp153 + Ftmp159 + Ftmp169);
  double Ftmp212 = Ftmp4 * (Ftmp152 + Ftmp156 + Ftmp170);
  double Ftmp213 = -Ftmp10 * Ftmp89 + Ftmp78;
  double Ftmp214 = y * M[4];
  double Ftmp215 = Ftmp15 + Ftmp22;
  double Ftmp216 = Ftmp185 * M[25];
  double Ftmp217 = 1.0 * y * M[29];
  double Ftmp218 = Ftmp57 * M[42];
  double Ftmp219 = Ftmp38 * x;
  double Ftmp220 = Ftmp21 * y;
  double Ftmp221 = Ftmp21 * (Ftmp79 + Ftmp81);
  double Ftmp222 = Ftmp100 - Ftmp103 * Ftmp21 + Ftmp110;
  double Ftmp223 = Ftmp137 + Ftmp204 + Ftmp6;
  double Ftmp224 = Ftmp141 + Ftmp143 + Ftmp205;
  double Ftmp225 = Ftmp135 + Ftmp143 + Ftmp146 + Ftmp26;
#pragma omp atomic
  F[0] += Ftmp10 * Ftmp33 * Ftmp64 - Ftmp10 * Ftmp4 * (Ftmp90 + Ftmp93) * M[35] -
          Ftmp10 * Ftmp65 * Ftmp97 * y - Ftmp10 * (Ftmp14 + Ftmp74) * M[6] -
          Ftmp10 * (Ftmp115 - 13230.0 * Ftmp116 + Ftmp132) * M[31] -
          Ftmp10 * (Ftmp149 + Ftmp172 + Ftmp173) * M[34] -
          Ftmp10 * (Ftmp159 + Ftmp173 + Ftmp174) * M[36] + Ftmp104 * x * M[16] +
          Ftmp104 * M[31] + Ftmp108 * M[41] + Ftmp11 * y * M[1] + Ftmp11 * z * M[2] +
          Ftmp111 * M[45] + Ftmp112 * x + Ftmp113 * x - Ftmp118 * M[32] -
          Ftmp122 * Ftmp52 - Ftmp126 * Ftmp57 - Ftmp127 * M[33] - Ftmp129 * Ftmp58 -
          Ftmp131 * Ftmp58 - Ftmp133 * Ftmp75 - Ftmp134 * Ftmp75 + Ftmp139 * x * M[19] +
          Ftmp139 * M[34] - Ftmp14 * Ftmp4 * M[10] + Ftmp142 * x * M[21] +
          Ftmp142 * M[36] + Ftmp148 * x * M[28] + Ftmp148 * M[43] - Ftmp155 * M[37] +
          Ftmp16 * x * M[0] + Ftmp16 * M[6] - Ftmp161 * M[39] - Ftmp165 * M[48] -
          Ftmp166 * M[38] - Ftmp167 * M[40] - Ftmp171 * Ftmp58 * M[49] - Ftmp2 * y +
          Ftmp20 * M[9] + Ftmp23 * M[11] + Ftmp24 * x + Ftmp25 * x - Ftmp29 * y -
          Ftmp3 * M[2] - Ftmp32 * Ftmp33 - Ftmp34 * Ftmp36 - Ftmp37 * z -
          Ftmp38 * Ftmp41 + Ftmp4 * Ftmp51 + Ftmp4 * Ftmp7 + Ftmp4 * Ftmp72 * M[20] -
          Ftmp44 * M[25] + Ftmp52 * Ftmp83 - Ftmp53 * M[7] - Ftmp54 * x -
          Ftmp56 * Ftmp57 - Ftmp58 * Ftmp60 - Ftmp58 * Ftmp61 - Ftmp59 * M[8] +
          Ftmp64 * Ftmp65 + Ftmp71 * z + Ftmp73 * Ftmp8 - Ftmp75 * Ftmp76 -
          Ftmp75 * Ftmp77 - Ftmp75 * (Ftmp164 + Ftmp175) * M[43] + Ftmp80 * y * M[17] +
          Ftmp80 * z * M[18] + Ftmp84 * Ftmp85 - Ftmp85 * Ftmp99 * z + Ftmp86 * Ftmp87 +
          Ftmp86 * Ftmp88 + Ftmp9 * M[4];
#pragma omp atomic
  F[1] += Ftmp111 * M[50] + Ftmp113 * y - Ftmp118 * M[31] - Ftmp121 * Ftmp4 * M[47] -
          Ftmp121 * Ftmp57 * M[41] - Ftmp126 * Ftmp191 - Ftmp131 * Ftmp4 -
          Ftmp134 * Ftmp57 - Ftmp155 * M[34] - Ftmp161 * M[36] - 1.0 * Ftmp165 * M[43] +
          Ftmp176 * Ftmp4 * Ftmp6 + Ftmp177 * M[7] + Ftmp178 * y * M[3] +
          Ftmp178 * M[12] + Ftmp179 * y + Ftmp18 * Ftmp180 * Ftmp195 +
          Ftmp18 * Ftmp181 * Ftmp84 - Ftmp18 * Ftmp186 * (Ftmp90 + Ftmp96) -
          Ftmp18 * Ftmp187 * Ftmp99 - Ftmp18 * Ftmp190 + Ftmp18 * Ftmp192 -
          Ftmp18 * Ftmp197 * Ftmp58 - Ftmp18 * Ftmp203 + Ftmp18 * Ftmp83 -
          Ftmp18 * (Ftmp160 + Ftmp175) * M[39] - Ftmp18 * (Ftmp209 + Ftmp213) * M[37] -
          Ftmp18 * (Ftmp30 + Ftmp74) * M[12] -
          Ftmp18 * (Ftmp119 + Ftmp132 - 13230.0 * Ftmp18 * Ftmp47) * M[46] -
          Ftmp18 * (Ftmp164 + Ftmp174 + Ftmp78) * M[48] - Ftmp180 * Ftmp31 -
          Ftmp181 * Ftmp36 - Ftmp183 * z - Ftmp184 * z - Ftmp185 * Ftmp32 * M[9] +
          Ftmp186 * Ftmp63 + Ftmp187 * Ftmp70 - Ftmp188 * Ftmp4 + Ftmp189 * Ftmp4 +
          Ftmp19 * x * M[1] + Ftmp19 * z * M[4] - Ftmp191 * Ftmp56 + Ftmp193 * Ftmp194 +
          Ftmp193 * Ftmp195 * M[27] + Ftmp196 * Ftmp63 * M[23] + Ftmp196 * Ftmp88 +
          Ftmp198 * M[32] + Ftmp199 * y * M[26] + Ftmp199 * M[46] - Ftmp2 * x +
          Ftmp200 * y - Ftmp202 * Ftmp4 + Ftmp206 * y * M[19] + Ftmp206 * M[37] +
          Ftmp207 * y * M[21] + Ftmp207 * M[39] + Ftmp208 * y * M[28] + Ftmp208 * M[48] -
          Ftmp210 * M[38] - Ftmp211 * M[40] - Ftmp212 * M[49] + Ftmp23 * M[14] +
          Ftmp25 * y - Ftmp29 * x - Ftmp3 * M[4] - Ftmp30 * Ftmp58 * M[10] -
          Ftmp32 * z * M[13] - Ftmp4 * Ftmp61 - Ftmp44 * M[29] + Ftmp51 * Ftmp58 -
          Ftmp53 * M[6] - Ftmp57 * Ftmp77 + Ftmp58 * Ftmp7;
#pragma omp atomic
  F[2] += -Ftmp1 * Ftmp176 - Ftmp1 * Ftmp214 + Ftmp108 * M[47] + Ftmp112 * z -
          Ftmp122 * Ftmp4 - Ftmp127 * M[31] - Ftmp129 * Ftmp21 -
          Ftmp130 * Ftmp219 * M[45] - Ftmp130 * Ftmp38 * y * M[50] - Ftmp133 * Ftmp219 -
          Ftmp166 * M[34] - Ftmp167 * M[36] - Ftmp171 * Ftmp219 * M[43] +
          Ftmp176 * Ftmp22 + Ftmp177 * M[8] + Ftmp179 * z + Ftmp180 * Ftmp4 * Ftmp63 -
          Ftmp183 * y - Ftmp184 * y + Ftmp185 * Ftmp21 * Ftmp87 - Ftmp185 * Ftmp41 -
          Ftmp188 * Ftmp21 + Ftmp189 * Ftmp21 - Ftmp190 * Ftmp4 + Ftmp192 * Ftmp4 +
          Ftmp194 * Ftmp220 - Ftmp197 * Ftmp21 * Ftmp52 + Ftmp198 * M[33] +
          Ftmp20 * M[13] + Ftmp200 * z - Ftmp202 * Ftmp21 - Ftmp203 * Ftmp4 -
          Ftmp21 * Ftmp218 * Ftmp97 -
          Ftmp21 * Ftmp34 * x * (-1575.0 * Ftmp47 + Ftmp98) * M[44] - Ftmp21 * Ftmp60 -
          Ftmp21 * (Ftmp42 + Ftmp74) * M[15] -
          Ftmp21 * (Ftmp123 + Ftmp132 - 13230.0 * Ftmp67) * M[51] -
          Ftmp21 * (Ftmp13 + Ftmp158 + Ftmp209) * M[38] -
          Ftmp21 * (Ftmp159 + Ftmp168 + Ftmp213) * M[40] -
          Ftmp21 * (Ftmp163 + Ftmp168 + Ftmp172 + Ftmp78) * M[49] - Ftmp210 * M[37] -
          Ftmp211 * M[39] - Ftmp212 * M[48] + Ftmp214 * Ftmp22 + Ftmp215 * z * M[5] +
          Ftmp215 * M[15] + Ftmp216 * Ftmp221 - Ftmp216 * Ftmp43 + Ftmp217 * Ftmp221 -
          Ftmp217 * Ftmp43 + Ftmp218 * Ftmp63 - Ftmp219 * Ftmp76 + Ftmp220 * Ftmp73 +
          Ftmp222 * z * M[30] + Ftmp222 * M[51] + Ftmp223 * z * M[19] + Ftmp223 * M[38] +
          Ftmp224 * z * M[21] + Ftmp224 * M[40] + Ftmp225 * z * M[28] + Ftmp225 * M[49] +
          Ftmp24 * z + Ftmp34 * Ftmp58 * Ftmp69 * M[24] - Ftmp37 * x -
          Ftmp42 * Ftmp52 * M[10] - Ftmp44 * x * M[11] - Ftmp44 * y * M[14] +
          Ftmp51 * Ftmp52 + Ftmp52 * Ftmp7 - Ftmp54 * z - Ftmp59 * M[6] + Ftmp71 * x +
          Ftmp9 * M[1];
}

void field_m2_P2M_6(double x, double y, double z, double q, double* M) {
  double Mtmp0  = (x * x);
  double Mtmp1  = (1.0 / 2.0) * q;
  double Mtmp2  = Mtmp0 * Mtmp1;
  double Mtmp3  = q * x;
  double Mtmp4  = Mtmp3 * y;
  double Mtmp5  = Mtmp3 * z;
  double Mtmp6  = (y * y);
  double Mtmp7  = Mtmp1 * Mtmp6;
  double Mtmp8  = q * y;
  double Mtmp9  = Mtmp8 * z;
  double Mtmp10 = (z * z);
  double Mtmp11 = Mtmp1 * Mtmp10;
  double Mtmp12 = (x * x * x);
  double Mtmp13 = (1.0 / 6.0) * q;
  double Mtmp14 = Mtmp12 * Mtmp13;
  double Mtmp15 = Mtmp2 * y;
  double Mtmp16 = Mtmp7 * x;
  double Mtmp17 = Mtmp11 * x;
  double Mtmp18 = (y * y * y);
  double Mtmp19 = Mtmp13 * Mtmp18;
  double Mtmp20 = (z * z * z);
  double Mtmp21 = (x * x * x * x);
  double Mtmp22 = (1.0 / 24.0) * q;
  double Mtmp23 = Mtmp21 * Mtmp22;
  double Mtmp24 = (1.0 / 6.0) * Mtmp8;
  double Mtmp25 = Mtmp6 * q;
  double Mtmp26 = (1.0 / 4.0) * Mtmp0;
  double Mtmp27 = Mtmp25 * Mtmp26;
  double Mtmp28 = Mtmp10 * q;
  double Mtmp29 = (1.0 / 6.0) * Mtmp3;
  double Mtmp30 = (y * y * y * y);
  double Mtmp31 = Mtmp22 * Mtmp30;
  double Mtmp32 = (1.0 / 4.0) * Mtmp10;
  double Mtmp33 = (z * z * z * z);
  double Mtmp34 = (x * x * x * x * x);
  double Mtmp35 = (1.0 / 120.0) * q;
  double Mtmp36 = Mtmp34 * Mtmp35;
  double Mtmp37 = (1.0 / 24.0) * Mtmp8;
  double Mtmp38 = (1.0 / 12.0) * Mtmp12;
  double Mtmp39 = Mtmp25 * Mtmp38;
  double Mtmp40 = (1.0 / 12.0) * Mtmp18;
  double Mtmp41 = Mtmp0 * q;
  double Mtmp42 = Mtmp40 * Mtmp41;
  double Mtmp43 = Mtmp10 * Mtmp8;
  double Mtmp44 = (1.0 / 12.0) * Mtmp20;
  double Mtmp45 = (1.0 / 24.0) * Mtmp3;
  double Mtmp46 = Mtmp3 * Mtmp6;
  double Mtmp47 = (y * y * y * y * y);
  double Mtmp48 = Mtmp35 * Mtmp47;
  double Mtmp49 = (z * z * z * z * z);
  double Mtmp50 = (1.0 / 720.0) * q;
  double Mtmp51 = (1.0 / 120.0) * Mtmp8;
  double Mtmp52 = (1.0 / 48.0) * Mtmp21;
  double Mtmp53 = (1.0 / 36.0) * Mtmp12 * q;
  double Mtmp54 = (1.0 / 48.0) * Mtmp41;
  double Mtmp55 = (1.0 / 120.0) * Mtmp3;
  M[0] += Mtmp2;
  M[1] += Mtmp4;
  M[2] += Mtmp5;
  M[3] += Mtmp7;
  M[4] += Mtmp9;
  M[5] += Mtmp11;
  M[6] += -Mtmp14;
  M[7] += -Mtmp15;
  M[8] += -Mtmp2 * z;
  M[9] += -Mtmp16;
  M[10] += -Mtmp4 * z;
  M[11] += -Mtmp17;
  M[12] += -Mtmp19;
  M[13] += -Mtmp7 * z;
  M[14] += -Mtmp11 * y;
  M[15] += -Mtmp13 * Mtmp20;
  M[16] += Mtmp23;
  M[17] += Mtmp12 * Mtmp24;
  M[18] += Mtmp14 * z;
  M[19] += Mtmp27;
  M[20] += Mtmp15 * z;
  M[21] += Mtmp26 * Mtmp28;
  M[22] += Mtmp18 * Mtmp29;
  M[23] += Mtmp16 * z;
  M[24] += Mtmp17 * y;
  M[25] += Mtmp20 * Mtmp29;
  M[26] += Mtmp31;
  M[27] += Mtmp19 * z;
  M[28] += Mtmp25 * Mtmp32;
  M[29] += Mtmp20 * Mtmp24;
  M[30] += Mtmp22 * Mtmp33;
  M[31] += -Mtmp36;
  M[32] += -Mtmp21 * Mtmp37;
  M[33] += -Mtmp23 * z;
  M[34] += -Mtmp39;
  M[35] += -1.0 / 6.0 * Mtmp12 * Mtmp9;
  M[36] += -Mtmp28 * Mtmp38;
  M[37] += -Mtmp42;
  M[38] += -Mtmp27 * z;
  M[39] += -Mtmp26 * Mtmp43;
  M[40] += -Mtmp41 * Mtmp44;
  M[41] += -Mtmp30 * Mtmp45;
  M[42] += -1.0 / 6.0 * Mtmp18 * Mtmp5;
  M[43] += -Mtmp32 * Mtmp46;
  M[44] += -1.0 / 6.0 * Mtmp20 * Mtmp4;
  M[45] += -Mtmp33 * Mtmp45;
  M[46] += -Mtmp48;
  M[47] += -Mtmp31 * z;
  M[48] += -Mtmp28 * Mtmp40;
  M[49] += -Mtmp25 * Mtmp44;
  M[50] += -Mtmp33 * Mtmp37;
  M[51] += -Mtmp35 * Mtmp49;
  M[52] += Mtmp50 * (x * x * x * x * x * x);
  M[53] += Mtmp34 * Mtmp51;
  M[54] += Mtmp36 * z;
  M[55] += Mtmp25 * Mtmp52;
  M[56] += (1.0 / 24.0) * Mtmp21 * Mtmp9;
  M[57] += Mtmp28 * Mtmp52;
  M[58] += Mtmp18 * Mtmp53;
  M[59] += Mtmp39 * z;
  M[60] += Mtmp38 * Mtmp43;
  M[61] += Mtmp20 * Mtmp53;
  M[62] += Mtmp30 * Mtmp54;
  M[63] += Mtmp42 * z;
  M[64] += (1.0 / 8.0) * Mtmp0 * Mtmp10 * Mtmp25;
  M[65] += Mtmp0 * Mtmp44 * Mtmp8;
  M[66] += Mtmp33 * Mtmp54;
  M[67] += Mtmp47 * Mtmp55;
  M[68] += (1.0 / 24.0) * Mtmp30 * Mtmp5;
  M[69] += Mtmp10 * Mtmp3 * Mtmp40;
  M[70] += Mtmp44 * Mtmp46;
  M[71] += (1.0 / 24.0) * Mtmp33 * Mtmp4;
  M[72] += Mtmp49 * Mtmp55;
  M[73] += Mtmp50 * (y * y * y * y * y * y);
  M[74] += Mtmp48 * z;
  M[75] += (1.0 / 48.0) * Mtmp28 * Mtmp30;
  M[76] += (1.0 / 36.0) * Mtmp18 * Mtmp20 * q;
  M[77] += (1.0 / 48.0) * Mtmp25 * Mtmp33;
  M[78] += Mtmp49 * Mtmp51;
  M[79] += Mtmp50 * (z * z * z * z * z * z);
}
void field_m2_M2M_6(double x, double y, double z, double* M, double* Ms) {
  double Mstmp0   = x * M[0];
  double Mstmp1   = x * M[1];
  double Mstmp2   = y * M[0];
  double Mstmp3   = x * M[2];
  double Mstmp4   = z * M[0];
  double Mstmp5   = x * M[3];
  double Mstmp6   = y * M[1];
  double Mstmp7   = x * M[4];
  double Mstmp8   = y * M[2];
  double Mstmp9   = z * M[1];
  double Mstmp10  = x * M[5];
  double Mstmp11  = z * M[2];
  double Mstmp12  = y * M[3];
  double Mstmp13  = y * M[4];
  double Mstmp14  = z * M[3];
  double Mstmp15  = y * M[5];
  double Mstmp16  = z * M[4];
  double Mstmp17  = z * M[5];
  double Mstmp18  = x * M[6];
  double Mstmp19  = (x * x);
  double Mstmp20  = (1.0 / 2.0) * Mstmp19;
  double Mstmp21  = x * M[7];
  double Mstmp22  = y * M[6];
  double Mstmp23  = Mstmp0 * y;
  double Mstmp24  = x * M[8];
  double Mstmp25  = z * M[6];
  double Mstmp26  = Mstmp0 * z;
  double Mstmp27  = x * M[9];
  double Mstmp28  = y * M[7];
  double Mstmp29  = Mstmp1 * y;
  double Mstmp30  = (y * y);
  double Mstmp31  = (1.0 / 2.0) * M[0];
  double Mstmp32  = x * M[10];
  double Mstmp33  = y * M[8];
  double Mstmp34  = z * M[7];
  double Mstmp35  = Mstmp3 * y;
  double Mstmp36  = Mstmp1 * z;
  double Mstmp37  = Mstmp2 * z;
  double Mstmp38  = x * M[11];
  double Mstmp39  = z * M[8];
  double Mstmp40  = Mstmp3 * z;
  double Mstmp41  = (z * z);
  double Mstmp42  = x * M[12];
  double Mstmp43  = y * M[9];
  double Mstmp44  = Mstmp5 * y;
  double Mstmp45  = (1.0 / 2.0) * Mstmp30;
  double Mstmp46  = x * M[13];
  double Mstmp47  = y * M[10];
  double Mstmp48  = z * M[9];
  double Mstmp49  = Mstmp7 * y;
  double Mstmp50  = Mstmp5 * z;
  double Mstmp51  = Mstmp6 * z;
  double Mstmp52  = x * M[14];
  double Mstmp53  = y * M[11];
  double Mstmp54  = z * M[10];
  double Mstmp55  = Mstmp10 * y;
  double Mstmp56  = Mstmp7 * z;
  double Mstmp57  = Mstmp8 * z;
  double Mstmp58  = (1.0 / 2.0) * Mstmp41;
  double Mstmp59  = x * M[15];
  double Mstmp60  = z * M[11];
  double Mstmp61  = Mstmp10 * z;
  double Mstmp62  = y * M[12];
  double Mstmp63  = y * M[13];
  double Mstmp64  = z * M[12];
  double Mstmp65  = Mstmp12 * z;
  double Mstmp66  = y * M[14];
  double Mstmp67  = z * M[13];
  double Mstmp68  = Mstmp13 * z;
  double Mstmp69  = y * M[15];
  double Mstmp70  = z * M[14];
  double Mstmp71  = Mstmp15 * z;
  double Mstmp72  = z * M[15];
  double Mstmp73  = x * M[16];
  double Mstmp74  = (1.0 / 6.0) * (x * x * x);
  double Mstmp75  = x * M[17];
  double Mstmp76  = y * M[16];
  double Mstmp77  = Mstmp18 * y;
  double Mstmp78  = x * M[18];
  double Mstmp79  = x * M[19];
  double Mstmp80  = y * M[17];
  double Mstmp81  = Mstmp21 * y;
  double Mstmp82  = x * M[20];
  double Mstmp83  = y * M[18];
  double Mstmp84  = Mstmp24 * y;
  double Mstmp85  = x * M[21];
  double Mstmp86  = x * M[22];
  double Mstmp87  = y * M[19];
  double Mstmp88  = Mstmp27 * y;
  double Mstmp89  = (y * y * y);
  double Mstmp90  = (1.0 / 6.0) * M[0];
  double Mstmp91  = x * M[23];
  double Mstmp92  = y * M[20];
  double Mstmp93  = Mstmp32 * y;
  double Mstmp94  = x * M[24];
  double Mstmp95  = y * M[21];
  double Mstmp96  = Mstmp38 * y;
  double Mstmp97  = x * M[25];
  double Mstmp98  = (z * z * z);
  double Mstmp99  = x * M[26];
  double Mstmp100 = y * M[22];
  double Mstmp101 = Mstmp42 * y;
  double Mstmp102 = (1.0 / 6.0) * Mstmp89;
  double Mstmp103 = x * M[27];
  double Mstmp104 = y * M[23];
  double Mstmp105 = Mstmp46 * y;
  double Mstmp106 = x * M[28];
  double Mstmp107 = y * M[24];
  double Mstmp108 = Mstmp52 * y;
  double Mstmp109 = x * M[29];
  double Mstmp110 = y * M[25];
  double Mstmp111 = Mstmp59 * y;
  double Mstmp112 = (1.0 / 6.0) * Mstmp98;
  double Mstmp113 = x * M[30];
  double Mstmp114 = y * M[26];
  double Mstmp115 = y * M[27];
  double Mstmp116 = y * M[28];
  double Mstmp117 = y * M[29];
  double Mstmp118 = y * M[30];
  double Mstmp119 = (1.0 / 24.0) * (x * x * x * x);
  double Mstmp120 = (1.0 / 4.0) * Mstmp19;
  double Mstmp121 = Mstmp120 * M[0];
  double Mstmp122 = Mstmp120 * Mstmp30;
  double Mstmp123 = Mstmp120 * Mstmp41;
  double Mstmp124 = (y * y * y * y);
  double Mstmp125 = (1.0 / 24.0) * M[0];
  double Mstmp126 = (1.0 / 4.0) * Mstmp30 * Mstmp41;
  double Mstmp127 = (z * z * z * z);
  double Mstmp128 = (1.0 / 24.0) * Mstmp124;
  double Mstmp129 = (1.0 / 24.0) * Mstmp127;
#pragma omp atomic
  Ms[0] += M[0];
#pragma omp atomic
  Ms[1] += M[1];
#pragma omp atomic
  Ms[2] += M[2];
#pragma omp atomic
  Ms[3] += M[3];
#pragma omp atomic
  Ms[4] += M[4];
#pragma omp atomic
  Ms[5] += M[5];
#pragma omp atomic
  Ms[6] += Mstmp0 + M[6];
#pragma omp atomic
  Ms[7] += Mstmp1 + Mstmp2 + M[7];
#pragma omp atomic
  Ms[8] += Mstmp3 + Mstmp4 + M[8];
#pragma omp atomic
  Ms[9] += Mstmp5 + Mstmp6 + M[9];
#pragma omp atomic
  Ms[10] += Mstmp7 + Mstmp8 + Mstmp9 + M[10];
#pragma omp atomic
  Ms[11] += Mstmp10 + Mstmp11 + M[11];
#pragma omp atomic
  Ms[12] += Mstmp12 + M[12];
#pragma omp atomic
  Ms[13] += Mstmp13 + Mstmp14 + M[13];
#pragma omp atomic
  Ms[14] += Mstmp15 + Mstmp16 + M[14];
#pragma omp atomic
  Ms[15] += Mstmp17 + M[15];
#pragma omp atomic
  Ms[16] += Mstmp18 + Mstmp20 * M[0] + M[16];
#pragma omp atomic
  Ms[17] += Mstmp20 * M[1] + Mstmp21 + Mstmp22 + Mstmp23 + M[17];
#pragma omp atomic
  Ms[18] += Mstmp20 * M[2] + Mstmp24 + Mstmp25 + Mstmp26 + M[18];
#pragma omp atomic
  Ms[19] += Mstmp20 * M[3] + Mstmp27 + Mstmp28 + Mstmp29 + Mstmp30 * Mstmp31 + M[19];
#pragma omp atomic
  Ms[20] += Mstmp20 * M[4] + Mstmp32 + Mstmp33 + Mstmp34 + Mstmp35 + Mstmp36 + Mstmp37 +
            M[20];
#pragma omp atomic
  Ms[21] += Mstmp20 * M[5] + Mstmp31 * Mstmp41 + Mstmp38 + Mstmp39 + Mstmp40 + M[21];
#pragma omp atomic
  Ms[22] += Mstmp42 + Mstmp43 + Mstmp44 + Mstmp45 * M[1] + M[22];
#pragma omp atomic
  Ms[23] += Mstmp45 * M[2] + Mstmp46 + Mstmp47 + Mstmp48 + Mstmp49 + Mstmp50 + Mstmp51 +
            M[23];
#pragma omp atomic
  Ms[24] += Mstmp52 + Mstmp53 + Mstmp54 + Mstmp55 + Mstmp56 + Mstmp57 + Mstmp58 * M[1] +
            M[24];
#pragma omp atomic
  Ms[25] += Mstmp58 * M[2] + Mstmp59 + Mstmp60 + Mstmp61 + M[25];
#pragma omp atomic
  Ms[26] += Mstmp45 * M[3] + Mstmp62 + M[26];
#pragma omp atomic
  Ms[27] += Mstmp45 * M[4] + Mstmp63 + Mstmp64 + Mstmp65 + M[27];
#pragma omp atomic
  Ms[28] += Mstmp45 * M[5] + Mstmp58 * M[3] + Mstmp66 + Mstmp67 + Mstmp68 + M[28];
#pragma omp atomic
  Ms[29] += Mstmp58 * M[4] + Mstmp69 + Mstmp70 + Mstmp71 + M[29];
#pragma omp atomic
  Ms[30] += Mstmp58 * M[5] + Mstmp72 + M[30];
#pragma omp atomic
  Ms[31] += Mstmp20 * M[6] + Mstmp73 + Mstmp74 * M[0] + M[31];
#pragma omp atomic
  Ms[32] += Mstmp2 * Mstmp20 + Mstmp20 * M[7] + Mstmp74 * M[1] + Mstmp75 + Mstmp76 +
            Mstmp77 + M[32];
#pragma omp atomic
  Ms[33] += Mstmp18 * z + Mstmp20 * Mstmp4 + Mstmp20 * M[8] + Mstmp74 * M[2] + Mstmp78 +
            z * M[16] + M[33];
#pragma omp atomic
  Ms[34] += Mstmp0 * Mstmp45 + Mstmp20 * Mstmp6 + Mstmp20 * M[9] + Mstmp45 * M[6] +
            Mstmp74 * M[3] + Mstmp79 + Mstmp80 + Mstmp81 + M[34];
#pragma omp atomic
  Ms[35] += Mstmp20 * Mstmp8 + Mstmp20 * Mstmp9 + Mstmp20 * M[10] + Mstmp21 * z +
            Mstmp22 * z + Mstmp23 * z + Mstmp74 * M[4] + Mstmp82 + Mstmp83 + Mstmp84 +
            z * M[17] + M[35];
#pragma omp atomic
  Ms[36] += Mstmp0 * Mstmp58 + Mstmp11 * Mstmp20 + Mstmp20 * M[11] + Mstmp24 * z +
            Mstmp58 * M[6] + Mstmp74 * M[5] + Mstmp85 + z * M[18] + M[36];
#pragma omp atomic
  Ms[37] += Mstmp1 * Mstmp45 + Mstmp12 * Mstmp20 + Mstmp20 * M[12] + Mstmp45 * M[7] +
            Mstmp86 + Mstmp87 + Mstmp88 + Mstmp89 * Mstmp90 + M[37];
#pragma omp atomic
  Ms[38] += Mstmp13 * Mstmp20 + Mstmp14 * Mstmp20 + Mstmp20 * M[13] + Mstmp27 * z +
            Mstmp28 * z + Mstmp29 * z + Mstmp3 * Mstmp45 + Mstmp4 * Mstmp45 +
            Mstmp45 * M[8] + Mstmp91 + Mstmp92 + Mstmp93 + z * M[19] + M[38];
#pragma omp atomic
  Ms[39] += Mstmp1 * Mstmp58 + Mstmp15 * Mstmp20 + Mstmp16 * Mstmp20 + Mstmp2 * Mstmp58 +
            Mstmp20 * M[14] + Mstmp32 * z + Mstmp33 * z + Mstmp35 * z + Mstmp58 * M[7] +
            Mstmp94 + Mstmp95 + Mstmp96 + z * M[20] + M[39];
#pragma omp atomic
  Ms[40] += Mstmp17 * Mstmp20 + Mstmp20 * M[15] + Mstmp3 * Mstmp58 + Mstmp38 * z +
            Mstmp58 * M[8] + Mstmp90 * Mstmp98 + Mstmp97 + z * M[21] + M[40];
#pragma omp atomic
  Ms[41] += Mstmp100 + Mstmp101 + Mstmp102 * M[1] + Mstmp45 * Mstmp5 + Mstmp45 * M[9] +
            Mstmp99 + M[41];
#pragma omp atomic
  Ms[42] += Mstmp102 * M[2] + Mstmp103 + Mstmp104 + Mstmp105 + Mstmp42 * z + Mstmp43 * z +
            Mstmp44 * z + Mstmp45 * Mstmp7 + Mstmp45 * Mstmp9 + Mstmp45 * M[10] +
            z * M[22] + M[42];
#pragma omp atomic
  Ms[43] += Mstmp10 * Mstmp45 + Mstmp106 + Mstmp107 + Mstmp108 + Mstmp11 * Mstmp45 +
            Mstmp45 * M[11] + Mstmp46 * z + Mstmp47 * z + Mstmp49 * z + Mstmp5 * Mstmp58 +
            Mstmp58 * Mstmp6 + Mstmp58 * M[9] + z * M[23] + M[43];
#pragma omp atomic
  Ms[44] += Mstmp109 + Mstmp110 + Mstmp111 + Mstmp112 * M[1] + Mstmp52 * z + Mstmp53 * z +
            Mstmp55 * z + Mstmp58 * Mstmp7 + Mstmp58 * Mstmp8 + Mstmp58 * M[10] +
            z * M[24] + M[44];
#pragma omp atomic
  Ms[45] += Mstmp10 * Mstmp58 + Mstmp112 * M[2] + Mstmp113 + Mstmp58 * M[11] +
            Mstmp59 * z + z * M[25] + M[45];
#pragma omp atomic
  Ms[46] += Mstmp102 * M[3] + Mstmp114 + Mstmp45 * M[12] + M[46];
#pragma omp atomic
  Ms[47] += Mstmp102 * M[4] + Mstmp115 + Mstmp14 * Mstmp45 + Mstmp45 * M[13] +
            Mstmp62 * z + z * M[26] + M[47];
#pragma omp atomic
  Ms[48] += Mstmp102 * M[5] + Mstmp116 + Mstmp12 * Mstmp58 + Mstmp16 * Mstmp45 +
            Mstmp45 * M[14] + Mstmp58 * M[12] + Mstmp63 * z + z * M[27] + M[48];
#pragma omp atomic
  Ms[49] += Mstmp112 * M[3] + Mstmp117 + Mstmp13 * Mstmp58 + Mstmp17 * Mstmp45 +
            Mstmp45 * M[15] + Mstmp58 * M[13] + Mstmp66 * z + z * M[28] + M[49];
#pragma omp atomic
  Ms[50] += Mstmp112 * M[4] + Mstmp118 + Mstmp15 * Mstmp58 + Mstmp58 * M[14] +
            Mstmp69 * z + z * M[29] + M[50];
#pragma omp atomic
  Ms[51] += Mstmp112 * M[5] + Mstmp58 * M[15] + z * M[30] + M[51];
#pragma omp atomic
  Ms[52] += Mstmp119 * M[0] + Mstmp20 * M[16] + Mstmp74 * M[6] + x * M[31] + M[52];
#pragma omp atomic
  Ms[53] += Mstmp119 * M[1] + Mstmp2 * Mstmp74 + Mstmp20 * Mstmp22 + Mstmp20 * M[17] +
            Mstmp73 * y + Mstmp74 * M[7] + x * M[32] + y * M[31] + M[53];
#pragma omp atomic
  Ms[54] += Mstmp119 * M[2] + Mstmp20 * Mstmp25 + Mstmp20 * M[18] + Mstmp4 * Mstmp74 +
            Mstmp73 * z + Mstmp74 * M[8] + x * M[33] + z * M[31] + M[54];
#pragma omp atomic
  Ms[55] += Mstmp119 * M[3] + Mstmp121 * Mstmp30 + Mstmp18 * Mstmp45 + Mstmp20 * Mstmp28 +
            Mstmp20 * M[19] + Mstmp45 * M[16] + Mstmp6 * Mstmp74 + Mstmp74 * M[9] +
            Mstmp75 * y + x * M[34] + y * M[32] + M[55];
#pragma omp atomic
  Ms[56] += Mstmp119 * M[4] + Mstmp20 * Mstmp33 + Mstmp20 * Mstmp34 + Mstmp20 * Mstmp37 +
            Mstmp20 * M[20] + Mstmp74 * Mstmp8 + Mstmp74 * Mstmp9 + Mstmp74 * M[10] +
            Mstmp75 * z + Mstmp76 * z + Mstmp77 * z + Mstmp78 * y + x * M[35] +
            y * M[33] + z * M[32] + M[56];
#pragma omp atomic
  Ms[57] += Mstmp11 * Mstmp74 + Mstmp119 * M[5] + Mstmp121 * Mstmp41 + Mstmp18 * Mstmp58 +
            Mstmp20 * Mstmp39 + Mstmp20 * M[21] + Mstmp58 * M[16] + Mstmp74 * M[11] +
            Mstmp78 * z + x * M[36] + z * M[33] + M[57];
#pragma omp atomic
  Ms[58] += Mstmp0 * Mstmp102 + Mstmp102 * M[6] + Mstmp12 * Mstmp74 + Mstmp122 * M[1] +
            Mstmp20 * Mstmp43 + Mstmp20 * M[22] + Mstmp21 * Mstmp45 + Mstmp45 * M[17] +
            Mstmp74 * M[12] + Mstmp79 * y + x * M[37] + y * M[34] + M[58];
#pragma omp atomic
  Ms[59] += Mstmp122 * M[2] + Mstmp13 * Mstmp74 + Mstmp14 * Mstmp74 + Mstmp20 * Mstmp47 +
            Mstmp20 * Mstmp48 + Mstmp20 * Mstmp51 + Mstmp20 * M[23] + Mstmp24 * Mstmp45 +
            Mstmp25 * Mstmp45 + Mstmp26 * Mstmp45 + Mstmp45 * M[18] + Mstmp74 * M[13] +
            Mstmp79 * z + Mstmp80 * z + Mstmp81 * z + Mstmp82 * y + x * M[38] +
            y * M[35] + z * M[34] + M[59];
#pragma omp atomic
  Ms[60] += Mstmp123 * M[1] + Mstmp15 * Mstmp74 + Mstmp16 * Mstmp74 + Mstmp20 * Mstmp53 +
            Mstmp20 * Mstmp54 + Mstmp20 * Mstmp57 + Mstmp20 * M[24] + Mstmp21 * Mstmp58 +
            Mstmp22 * Mstmp58 + Mstmp23 * Mstmp58 + Mstmp58 * M[17] + Mstmp74 * M[14] +
            Mstmp82 * z + Mstmp83 * z + Mstmp84 * z + Mstmp85 * y + x * M[39] +
            y * M[36] + z * M[35] + M[60];
#pragma omp atomic
  Ms[61] += Mstmp0 * Mstmp112 + Mstmp112 * M[6] + Mstmp123 * M[2] + Mstmp17 * Mstmp74 +
            Mstmp20 * Mstmp60 + Mstmp20 * M[25] + Mstmp24 * Mstmp58 + Mstmp58 * M[18] +
            Mstmp74 * M[15] + Mstmp85 * z + x * M[40] + z * M[36] + M[61];
#pragma omp atomic
  Ms[62] += Mstmp1 * Mstmp102 + Mstmp102 * M[7] + Mstmp122 * M[3] + Mstmp124 * Mstmp125 +
            Mstmp20 * Mstmp62 + Mstmp20 * M[26] + Mstmp27 * Mstmp45 + Mstmp45 * M[19] +
            Mstmp86 * y + x * M[41] + y * M[37] + M[62];
#pragma omp atomic
  Ms[63] += Mstmp102 * Mstmp3 + Mstmp102 * Mstmp4 + Mstmp102 * M[8] + Mstmp122 * M[4] +
            Mstmp20 * Mstmp63 + Mstmp20 * Mstmp64 + Mstmp20 * Mstmp65 + Mstmp20 * M[27] +
            Mstmp32 * Mstmp45 + Mstmp34 * Mstmp45 + Mstmp36 * Mstmp45 + Mstmp45 * M[20] +
            Mstmp86 * z + Mstmp87 * z + Mstmp88 * z + Mstmp91 * y + x * M[42] +
            y * M[38] + z * M[37] + M[63];
#pragma omp atomic
  Ms[64] += Mstmp122 * M[5] + Mstmp123 * M[3] + Mstmp126 * M[0] + Mstmp20 * Mstmp66 +
            Mstmp20 * Mstmp67 + Mstmp20 * Mstmp68 + Mstmp20 * M[28] + Mstmp27 * Mstmp58 +
            Mstmp28 * Mstmp58 + Mstmp29 * Mstmp58 + Mstmp38 * Mstmp45 +
            Mstmp39 * Mstmp45 + Mstmp40 * Mstmp45 + Mstmp45 * M[21] + Mstmp58 * M[19] +
            Mstmp91 * z + Mstmp92 * z + Mstmp93 * z + Mstmp94 * y + x * M[43] +
            y * M[39] + z * M[38] + M[64];
#pragma omp atomic
  Ms[65] += Mstmp1 * Mstmp112 + Mstmp112 * Mstmp2 + Mstmp112 * M[7] + Mstmp123 * M[4] +
            Mstmp20 * Mstmp69 + Mstmp20 * Mstmp70 + Mstmp20 * Mstmp71 + Mstmp20 * M[29] +
            Mstmp32 * Mstmp58 + Mstmp33 * Mstmp58 + Mstmp35 * Mstmp58 + Mstmp58 * M[20] +
            Mstmp94 * z + Mstmp95 * z + Mstmp96 * z + Mstmp97 * y + x * M[44] +
            y * M[40] + z * M[39] + M[65];
#pragma omp atomic
  Ms[66] += Mstmp112 * Mstmp3 + Mstmp112 * M[8] + Mstmp123 * M[5] + Mstmp125 * Mstmp127 +
            Mstmp20 * Mstmp72 + Mstmp20 * M[30] + Mstmp38 * Mstmp58 + Mstmp58 * M[21] +
            Mstmp97 * z + x * M[45] + z * M[40] + M[66];
#pragma omp atomic
  Ms[67] += Mstmp102 * Mstmp5 + Mstmp102 * M[9] + Mstmp128 * M[1] + Mstmp42 * Mstmp45 +
            Mstmp45 * M[22] + Mstmp99 * y + x * M[46] + y * M[41] + M[67];
#pragma omp atomic
  Ms[68] += Mstmp100 * z + Mstmp101 * z + Mstmp102 * Mstmp7 + Mstmp102 * Mstmp9 +
            Mstmp102 * M[10] + Mstmp103 * y + Mstmp128 * M[2] + Mstmp45 * Mstmp46 +
            Mstmp45 * Mstmp48 + Mstmp45 * Mstmp50 + Mstmp45 * M[23] + Mstmp99 * z +
            x * M[47] + y * M[42] + z * M[41] + M[68];
#pragma omp atomic
  Ms[69] += Mstmp10 * Mstmp102 + Mstmp102 * Mstmp11 + Mstmp102 * M[11] + Mstmp103 * z +
            Mstmp104 * z + Mstmp105 * z + Mstmp106 * y + Mstmp126 * M[1] +
            Mstmp42 * Mstmp58 + Mstmp43 * Mstmp58 + Mstmp44 * Mstmp58 +
            Mstmp45 * Mstmp52 + Mstmp45 * Mstmp54 + Mstmp45 * Mstmp56 + Mstmp45 * M[24] +
            Mstmp58 * M[22] + x * M[48] + y * M[43] + z * M[42] + M[69];
#pragma omp atomic
  Ms[70] += Mstmp106 * z + Mstmp107 * z + Mstmp108 * z + Mstmp109 * y +
            Mstmp112 * Mstmp5 + Mstmp112 * Mstmp6 + Mstmp112 * M[9] + Mstmp126 * M[2] +
            Mstmp45 * Mstmp59 + Mstmp45 * Mstmp60 + Mstmp45 * Mstmp61 + Mstmp45 * M[25] +
            Mstmp46 * Mstmp58 + Mstmp47 * Mstmp58 + Mstmp49 * Mstmp58 + Mstmp58 * M[23] +
            x * M[49] + y * M[44] + z * M[43] + M[70];
#pragma omp atomic
  Ms[71] += Mstmp109 * z + Mstmp110 * z + Mstmp111 * z + Mstmp112 * Mstmp7 +
            Mstmp112 * Mstmp8 + Mstmp112 * M[10] + Mstmp113 * y + Mstmp129 * M[1] +
            Mstmp52 * Mstmp58 + Mstmp53 * Mstmp58 + Mstmp55 * Mstmp58 + Mstmp58 * M[24] +
            x * M[50] + y * M[45] + z * M[44] + M[71];
#pragma omp atomic
  Ms[72] += Mstmp10 * Mstmp112 + Mstmp112 * M[11] + Mstmp113 * z + Mstmp129 * M[2] +
            Mstmp58 * Mstmp59 + Mstmp58 * M[25] + x * M[51] + z * M[45] + M[72];
#pragma omp atomic
  Ms[73] += Mstmp102 * M[12] + Mstmp128 * M[3] + Mstmp45 * M[26] + y * M[46] + M[73];
#pragma omp atomic
  Ms[74] += Mstmp102 * Mstmp14 + Mstmp102 * M[13] + Mstmp114 * z + Mstmp128 * M[4] +
            Mstmp45 * Mstmp64 + Mstmp45 * M[27] + y * M[47] + z * M[46] + M[74];
#pragma omp atomic
  Ms[75] += Mstmp102 * Mstmp16 + Mstmp102 * M[14] + Mstmp115 * z + Mstmp126 * M[3] +
            Mstmp128 * M[5] + Mstmp45 * Mstmp67 + Mstmp45 * M[28] + Mstmp58 * Mstmp62 +
            Mstmp58 * M[26] + y * M[48] + z * M[47] + M[75];
#pragma omp atomic
  Ms[76] += Mstmp102 * Mstmp17 + Mstmp102 * M[15] + Mstmp112 * Mstmp12 +
            Mstmp112 * M[12] + Mstmp116 * z + Mstmp126 * M[4] + Mstmp45 * Mstmp70 +
            Mstmp45 * M[29] + Mstmp58 * Mstmp63 + Mstmp58 * M[27] + y * M[49] +
            z * M[48] + M[76];
#pragma omp atomic
  Ms[77] += Mstmp112 * Mstmp13 + Mstmp112 * M[13] + Mstmp117 * z + Mstmp126 * M[5] +
            Mstmp129 * M[3] + Mstmp45 * Mstmp72 + Mstmp45 * M[30] + Mstmp58 * Mstmp66 +
            Mstmp58 * M[28] + y * M[50] + z * M[49] + M[77];
#pragma omp atomic
  Ms[78] += Mstmp112 * Mstmp15 + Mstmp112 * M[14] + Mstmp118 * z + Mstmp129 * M[4] +
            Mstmp58 * Mstmp69 + Mstmp58 * M[29] + y * M[51] + z * M[50] + M[78];
#pragma omp atomic
  Ms[79] += Mstmp112 * M[15] + Mstmp129 * M[5] + Mstmp58 * M[30] + z * M[51] + M[79];
}

void field_m2_M2L_6(double x, double y, double z, double* M, double* L) {
  double R = sqrt(x * x + y * y + z * z);
  double D[80];
  double Dtmp0  = -1.0 * pow(R, -3.0);
  double Dtmp1  = (x * x);
  double Dtmp2  = pow(R, -5.0);
  double Dtmp3  = 3.0 * Dtmp2;
  double Dtmp4  = x * y;
  double Dtmp5  = x * z;
  double Dtmp6  = (y * y);
  double Dtmp7  = y * z;
  double Dtmp8  = 9.0 * Dtmp2;
  double Dtmp9  = -Dtmp8;
  double Dtmp10 = pow(R, -7.0);
  double Dtmp11 = 15.0 * Dtmp10;
  double Dtmp12 = Dtmp1 * Dtmp11;
  double Dtmp13 = -Dtmp3;
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
  double Dtmp25 = Dtmp1 * Dtmp21;
  double Dtmp26 = x * (Dtmp24 + Dtmp25);
  double Dtmp27 = -Dtmp11;
  double Dtmp28 = Dtmp21 * Dtmp6;
  double Dtmp29 = Dtmp24 + Dtmp28;
  double Dtmp30 = Dtmp17 * y;
  double Dtmp31 = Dtmp17 * z;
  double Dtmp32 = (y * y * y * y);
  double Dtmp33 = 225.0 * Dtmp10;
  double Dtmp34 = pow(R, -11.0);
  double Dtmp35 = 945.0 * Dtmp34;
  double Dtmp36 = Dtmp19 * Dtmp35;
  double Dtmp37 = Dtmp1 * Dtmp20;
  double Dtmp38 = 630.0 * Dtmp37;
  double Dtmp39 = Dtmp23 + Dtmp36 - Dtmp38;
  double Dtmp40 = -Dtmp25;
  double Dtmp41 = 315.0 * Dtmp20;
  double Dtmp42 = Dtmp41 * Dtmp6;
  double Dtmp43 = Dtmp1 * Dtmp35;
  double Dtmp44 = Dtmp43 * Dtmp6;
  double Dtmp45 = Dtmp23 + Dtmp44;
  double Dtmp46 = -Dtmp41;
  double Dtmp47 = -Dtmp28;
  double Dtmp48 = Dtmp1 * Dtmp41;
  double Dtmp49 = Dtmp32 * Dtmp35;
  double Dtmp50 = Dtmp20 * Dtmp6;
  double Dtmp51 = 630.0 * Dtmp50;
  double Dtmp52 = Dtmp23 + Dtmp49 - Dtmp51;
  double Dtmp53 = Dtmp35 * Dtmp6;
  double Dtmp54 = -Dtmp33;
  double Dtmp55 = 10395.0 * pow(R, -13.0);
  double Dtmp56 = 14175.0 * Dtmp34;
  double Dtmp57 = 1575.0 * Dtmp20;
  double Dtmp58 = Dtmp19 * Dtmp55;
  double Dtmp59 = Dtmp1 * Dtmp34;
  double Dtmp60 = x * (Dtmp57 + Dtmp58 - 9450.0 * Dtmp59);
  double Dtmp61 = 5670.0 * Dtmp59;
  double Dtmp62 = Dtmp24 - Dtmp6 * Dtmp61;
  double Dtmp63 = -2835.0 * Dtmp59;
  double Dtmp64 = Dtmp34 * Dtmp6;
  double Dtmp65 = Dtmp1 * Dtmp55 * Dtmp6;
  double Dtmp66 = -2835.0 * Dtmp64 + Dtmp65;
  double Dtmp67 = Dtmp32 * Dtmp55;
  double Dtmp68 = Dtmp57 - 9450.0 * Dtmp64 + Dtmp67;
  D[0]          = Dtmp0 + Dtmp1 * Dtmp3;
  D[1]          = Dtmp3 * Dtmp4;
  D[2]          = Dtmp3 * Dtmp5;
  D[3]          = Dtmp0 + Dtmp3 * Dtmp6;
  D[4]          = Dtmp3 * Dtmp7;
  D[5]          = -D[0] - D[3];
  D[6]          = -x * (Dtmp12 + Dtmp9);
  D[7]          = -Dtmp14 * y;
  D[8]          = -Dtmp14 * z;
  D[9]          = -Dtmp16 * Dtmp17;
  D[10]         = -Dtmp11 * Dtmp18;
  D[11]         = -D[6] - D[9];
  D[12]         = -y * (Dtmp15 + Dtmp9);
  D[13]         = -Dtmp16 * z;
  D[14]         = -D[7] - D[12];
  D[15]         = -D[8] - D[13];
  D[16]         = -Dtmp1 * Dtmp22 + Dtmp19 * Dtmp21 + Dtmp8;
  D[17]         = Dtmp26 * y;
  D[18]         = Dtmp26 * z;
  D[19]         = -Dtmp12 - Dtmp15 + Dtmp25 * Dtmp6 + Dtmp3;
  D[20]         = Dtmp7 * (Dtmp25 + Dtmp27);
  D[21]         = -D[16] - D[19];
  D[22]         = Dtmp29 * Dtmp30;
  D[23]         = Dtmp31 * (Dtmp27 + Dtmp28);
  D[24]         = -D[17] - D[22];
  D[25]         = -D[18] - D[23];
  D[26]         = Dtmp21 * Dtmp32 - Dtmp22 * Dtmp6 + Dtmp8;
  D[27]         = Dtmp29 * Dtmp7;
  D[28]         = -D[19] - D[26];
  D[29]         = -D[20] - D[27];
  D[30]         = -D[21] - D[28];
  D[31]         = -x * (Dtmp33 + Dtmp36 - 1050.0 * Dtmp37);
  D[32]         = -Dtmp39 * y;
  D[33]         = -Dtmp39 * z;
  D[34]         = -x * (Dtmp40 - Dtmp42 + Dtmp45);
  D[35]         = -Dtmp18 * (Dtmp43 + Dtmp46);
  D[36]         = -D[31] - D[34];
  D[37]         = -y * (Dtmp45 + Dtmp47 - Dtmp48);
  D[38]         = -z * (Dtmp11 + Dtmp40 + Dtmp44 + Dtmp47);
  D[39]         = -D[32] - D[37];
  D[40]         = -D[33] - D[38];
  D[41]         = -Dtmp17 * Dtmp52;
  D[42]         = -Dtmp17 * Dtmp7 * (Dtmp46 + Dtmp53);
  D[43]         = -D[34] - D[41];
  D[44]         = -D[35] - D[42];
  D[45]         = -D[36] - D[43];
  D[46]         = -y * (Dtmp33 + Dtmp49 - 1050.0 * Dtmp50);
  D[47]         = -Dtmp52 * z;
  D[48]         = -D[37] - D[46];
  D[49]         = -D[38] - D[47];
  D[50]         = -D[39] - D[48];
  D[51]         = -D[40] - D[49];
  D[52] = -Dtmp19 * Dtmp56 + 4725.0 * Dtmp37 + Dtmp54 + Dtmp55 * (x * x * x * x * x * x);
  D[53] = Dtmp60 * y;
  D[54] = Dtmp60 * z;
  D[55] = -Dtmp36 + Dtmp38 + Dtmp42 + Dtmp58 * Dtmp6 + Dtmp62;
  D[56] = Dtmp7 * (Dtmp41 + Dtmp58 - Dtmp61);
  D[57] = -D[52] - D[55];
  D[58] = Dtmp4 * (945.0 * Dtmp20 + Dtmp63 + Dtmp66);
  D[59] = Dtmp5 * (Dtmp41 - Dtmp43 + Dtmp66);
  D[60] = -D[53] - D[58];
  D[61] = -D[54] - D[59];
  D[62] = Dtmp1 * Dtmp67 + Dtmp48 - Dtmp49 + Dtmp51 + Dtmp62;
  D[63] = Dtmp7 * (Dtmp41 - Dtmp53 + Dtmp63 + Dtmp65);
  D[64] = -D[55] - D[62];
  D[65] = -D[56] - D[63];
  D[66] = -D[57] - D[64];
  D[67] = Dtmp30 * Dtmp68;
  D[68] = Dtmp31 * (Dtmp41 - 5670.0 * Dtmp64 + Dtmp67);
  D[69] = -D[58] - D[67];
  D[70] = -D[59] - D[68];
  D[71] = -D[60] - D[69];
  D[72] = -D[61] - D[70];
  D[73] = -Dtmp32 * Dtmp56 + 4725.0 * Dtmp50 + Dtmp54 + Dtmp55 * (y * y * y * y * y * y);
  D[74] = Dtmp68 * Dtmp7;
  D[75] = -D[62] - D[73];
  D[76] = -D[63] - D[74];
  D[77] = -D[64] - D[75];
  D[78] = -D[65] - D[76];
  D[79] = -D[66] - D[77];
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
          D[75] * M[75] + D[76] * M[76] + D[77] * M[77] + D[78] * M[78] + D[79] * M[79];
#pragma omp atomic
  L[1] += D[6] * M[0] + D[7] * M[1] + D[8] * M[2] + D[9] * M[3] + D[10] * M[4] +
          D[11] * M[5] + D[16] * M[6] + D[17] * M[7] + D[18] * M[8] + D[19] * M[9] +
          D[20] * M[10] + D[21] * M[11] + D[22] * M[12] + D[23] * M[13] + D[24] * M[14] +
          D[25] * M[15] + D[31] * M[16] + D[32] * M[17] + D[33] * M[18] + D[34] * M[19] +
          D[35] * M[20] + D[36] * M[21] + D[37] * M[22] + D[38] * M[23] + D[39] * M[24] +
          D[40] * M[25] + D[41] * M[26] + D[42] * M[27] + D[43] * M[28] + D[44] * M[29] +
          D[45] * M[30] + D[52] * M[31] + D[53] * M[32] + D[54] * M[33] + D[55] * M[34] +
          D[56] * M[35] + D[57] * M[36] + D[58] * M[37] + D[59] * M[38] + D[60] * M[39] +
          D[61] * M[40] + D[62] * M[41] + D[63] * M[42] + D[64] * M[43] + D[65] * M[44] +
          D[66] * M[45] + D[67] * M[46] + D[68] * M[47] + D[69] * M[48] + D[70] * M[49] +
          D[71] * M[50] + D[72] * M[51];
#pragma omp atomic
  L[2] += D[7] * M[0] + D[9] * M[1] + D[10] * M[2] + D[12] * M[3] + D[13] * M[4] +
          D[14] * M[5] + D[17] * M[6] + D[19] * M[7] + D[20] * M[8] + D[22] * M[9] +
          D[23] * M[10] + D[24] * M[11] + D[26] * M[12] + D[27] * M[13] + D[28] * M[14] +
          D[29] * M[15] + D[32] * M[16] + D[34] * M[17] + D[35] * M[18] + D[37] * M[19] +
          D[38] * M[20] + D[39] * M[21] + D[41] * M[22] + D[42] * M[23] + D[43] * M[24] +
          D[44] * M[25] + D[46] * M[26] + D[47] * M[27] + D[48] * M[28] + D[49] * M[29] +
          D[50] * M[30] + D[53] * M[31] + D[55] * M[32] + D[56] * M[33] + D[58] * M[34] +
          D[59] * M[35] + D[60] * M[36] + D[62] * M[37] + D[63] * M[38] + D[64] * M[39] +
          D[65] * M[40] + D[67] * M[41] + D[68] * M[42] + D[69] * M[43] + D[70] * M[44] +
          D[71] * M[45] + D[73] * M[46] + D[74] * M[47] + D[75] * M[48] + D[76] * M[49] +
          D[77] * M[50] + D[78] * M[51];
#pragma omp atomic
  L[3] += D[8] * M[0] + D[10] * M[1] + D[11] * M[2] + D[13] * M[3] + D[14] * M[4] +
          D[15] * M[5] + D[18] * M[6] + D[20] * M[7] + D[21] * M[8] + D[23] * M[9] +
          D[24] * M[10] + D[25] * M[11] + D[27] * M[12] + D[28] * M[13] + D[29] * M[14] +
          D[30] * M[15] + D[33] * M[16] + D[35] * M[17] + D[36] * M[18] + D[38] * M[19] +
          D[39] * M[20] + D[40] * M[21] + D[42] * M[22] + D[43] * M[23] + D[44] * M[24] +
          D[45] * M[25] + D[47] * M[26] + D[48] * M[27] + D[49] * M[28] + D[50] * M[29] +
          D[51] * M[30] + D[54] * M[31] + D[56] * M[32] + D[57] * M[33] + D[59] * M[34] +
          D[60] * M[35] + D[61] * M[36] + D[63] * M[37] + D[64] * M[38] + D[65] * M[39] +
          D[66] * M[40] + D[68] * M[41] + D[69] * M[42] + D[70] * M[43] + D[71] * M[44] +
          D[72] * M[45] + D[74] * M[46] + D[75] * M[47] + D[76] * M[48] + D[77] * M[49] +
          D[78] * M[50] + D[79] * M[51];
#pragma omp atomic
  L[4] += D[16] * M[0] + D[17] * M[1] + D[18] * M[2] + D[19] * M[3] + D[20] * M[4] +
          D[21] * M[5] + D[31] * M[6] + D[32] * M[7] + D[33] * M[8] + D[34] * M[9] +
          D[35] * M[10] + D[36] * M[11] + D[37] * M[12] + D[38] * M[13] + D[39] * M[14] +
          D[40] * M[15] + D[52] * M[16] + D[53] * M[17] + D[54] * M[18] + D[55] * M[19] +
          D[56] * M[20] + D[57] * M[21] + D[58] * M[22] + D[59] * M[23] + D[60] * M[24] +
          D[61] * M[25] + D[62] * M[26] + D[63] * M[27] + D[64] * M[28] + D[65] * M[29] +
          D[66] * M[30];
#pragma omp atomic
  L[5] += D[17] * M[0] + D[19] * M[1] + D[20] * M[2] + D[22] * M[3] + D[23] * M[4] +
          D[24] * M[5] + D[32] * M[6] + D[34] * M[7] + D[35] * M[8] + D[37] * M[9] +
          D[38] * M[10] + D[39] * M[11] + D[41] * M[12] + D[42] * M[13] + D[43] * M[14] +
          D[44] * M[15] + D[53] * M[16] + D[55] * M[17] + D[56] * M[18] + D[58] * M[19] +
          D[59] * M[20] + D[60] * M[21] + D[62] * M[22] + D[63] * M[23] + D[64] * M[24] +
          D[65] * M[25] + D[67] * M[26] + D[68] * M[27] + D[69] * M[28] + D[70] * M[29] +
          D[71] * M[30];
#pragma omp atomic
  L[6] += D[18] * M[0] + D[20] * M[1] + D[21] * M[2] + D[23] * M[3] + D[24] * M[4] +
          D[25] * M[5] + D[33] * M[6] + D[35] * M[7] + D[36] * M[8] + D[38] * M[9] +
          D[39] * M[10] + D[40] * M[11] + D[42] * M[12] + D[43] * M[13] + D[44] * M[14] +
          D[45] * M[15] + D[54] * M[16] + D[56] * M[17] + D[57] * M[18] + D[59] * M[19] +
          D[60] * M[20] + D[61] * M[21] + D[63] * M[22] + D[64] * M[23] + D[65] * M[24] +
          D[66] * M[25] + D[68] * M[26] + D[69] * M[27] + D[70] * M[28] + D[71] * M[29] +
          D[72] * M[30];
#pragma omp atomic
  L[7] += D[19] * M[0] + D[22] * M[1] + D[23] * M[2] + D[26] * M[3] + D[27] * M[4] +
          D[28] * M[5] + D[34] * M[6] + D[37] * M[7] + D[38] * M[8] + D[41] * M[9] +
          D[42] * M[10] + D[43] * M[11] + D[46] * M[12] + D[47] * M[13] + D[48] * M[14] +
          D[49] * M[15] + D[55] * M[16] + D[58] * M[17] + D[59] * M[18] + D[62] * M[19] +
          D[63] * M[20] + D[64] * M[21] + D[67] * M[22] + D[68] * M[23] + D[69] * M[24] +
          D[70] * M[25] + D[73] * M[26] + D[74] * M[27] + D[75] * M[28] + D[76] * M[29] +
          D[77] * M[30];
#pragma omp atomic
  L[8] += D[20] * M[0] + D[23] * M[1] + D[24] * M[2] + D[27] * M[3] + D[28] * M[4] +
          D[29] * M[5] + D[35] * M[6] + D[38] * M[7] + D[39] * M[8] + D[42] * M[9] +
          D[43] * M[10] + D[44] * M[11] + D[47] * M[12] + D[48] * M[13] + D[49] * M[14] +
          D[50] * M[15] + D[56] * M[16] + D[59] * M[17] + D[60] * M[18] + D[63] * M[19] +
          D[64] * M[20] + D[65] * M[21] + D[68] * M[22] + D[69] * M[23] + D[70] * M[24] +
          D[71] * M[25] + D[74] * M[26] + D[75] * M[27] + D[76] * M[28] + D[77] * M[29] +
          D[78] * M[30];
#pragma omp atomic
  L[9] += D[21] * M[0] + D[24] * M[1] + D[25] * M[2] + D[28] * M[3] + D[29] * M[4] +
          D[30] * M[5] + D[36] * M[6] + D[39] * M[7] + D[40] * M[8] + D[43] * M[9] +
          D[44] * M[10] + D[45] * M[11] + D[48] * M[12] + D[49] * M[13] + D[50] * M[14] +
          D[51] * M[15] + D[57] * M[16] + D[60] * M[17] + D[61] * M[18] + D[64] * M[19] +
          D[65] * M[20] + D[66] * M[21] + D[69] * M[22] + D[70] * M[23] + D[71] * M[24] +
          D[72] * M[25] + D[75] * M[26] + D[76] * M[27] + D[77] * M[28] + D[78] * M[29] +
          D[79] * M[30];
#pragma omp atomic
  L[10] += D[31] * M[0] + D[32] * M[1] + D[33] * M[2] + D[34] * M[3] + D[35] * M[4] +
           D[36] * M[5] + D[52] * M[6] + D[53] * M[7] + D[54] * M[8] + D[55] * M[9] +
           D[56] * M[10] + D[57] * M[11] + D[58] * M[12] + D[59] * M[13] + D[60] * M[14] +
           D[61] * M[15];
#pragma omp atomic
  L[11] += D[32] * M[0] + D[34] * M[1] + D[35] * M[2] + D[37] * M[3] + D[38] * M[4] +
           D[39] * M[5] + D[53] * M[6] + D[55] * M[7] + D[56] * M[8] + D[58] * M[9] +
           D[59] * M[10] + D[60] * M[11] + D[62] * M[12] + D[63] * M[13] + D[64] * M[14] +
           D[65] * M[15];
#pragma omp atomic
  L[12] += D[33] * M[0] + D[35] * M[1] + D[36] * M[2] + D[38] * M[3] + D[39] * M[4] +
           D[40] * M[5] + D[54] * M[6] + D[56] * M[7] + D[57] * M[8] + D[59] * M[9] +
           D[60] * M[10] + D[61] * M[11] + D[63] * M[12] + D[64] * M[13] + D[65] * M[14] +
           D[66] * M[15];
#pragma omp atomic
  L[13] += D[34] * M[0] + D[37] * M[1] + D[38] * M[2] + D[41] * M[3] + D[42] * M[4] +
           D[43] * M[5] + D[55] * M[6] + D[58] * M[7] + D[59] * M[8] + D[62] * M[9] +
           D[63] * M[10] + D[64] * M[11] + D[67] * M[12] + D[68] * M[13] + D[69] * M[14] +
           D[70] * M[15];
#pragma omp atomic
  L[14] += D[35] * M[0] + D[38] * M[1] + D[39] * M[2] + D[42] * M[3] + D[43] * M[4] +
           D[44] * M[5] + D[56] * M[6] + D[59] * M[7] + D[60] * M[8] + D[63] * M[9] +
           D[64] * M[10] + D[65] * M[11] + D[68] * M[12] + D[69] * M[13] + D[70] * M[14] +
           D[71] * M[15];
#pragma omp atomic
  L[15] += D[36] * M[0] + D[39] * M[1] + D[40] * M[2] + D[43] * M[3] + D[44] * M[4] +
           D[45] * M[5] + D[57] * M[6] + D[60] * M[7] + D[61] * M[8] + D[64] * M[9] +
           D[65] * M[10] + D[66] * M[11] + D[69] * M[12] + D[70] * M[13] + D[71] * M[14] +
           D[72] * M[15];
#pragma omp atomic
  L[16] += D[37] * M[0] + D[41] * M[1] + D[42] * M[2] + D[46] * M[3] + D[47] * M[4] +
           D[48] * M[5] + D[58] * M[6] + D[62] * M[7] + D[63] * M[8] + D[67] * M[9] +
           D[68] * M[10] + D[69] * M[11] + D[73] * M[12] + D[74] * M[13] + D[75] * M[14] +
           D[76] * M[15];
#pragma omp atomic
  L[17] += D[38] * M[0] + D[42] * M[1] + D[43] * M[2] + D[47] * M[3] + D[48] * M[4] +
           D[49] * M[5] + D[59] * M[6] + D[63] * M[7] + D[64] * M[8] + D[68] * M[9] +
           D[69] * M[10] + D[70] * M[11] + D[74] * M[12] + D[75] * M[13] + D[76] * M[14] +
           D[77] * M[15];
#pragma omp atomic
  L[18] += D[39] * M[0] + D[43] * M[1] + D[44] * M[2] + D[48] * M[3] + D[49] * M[4] +
           D[50] * M[5] + D[60] * M[6] + D[64] * M[7] + D[65] * M[8] + D[69] * M[9] +
           D[70] * M[10] + D[71] * M[11] + D[75] * M[12] + D[76] * M[13] + D[77] * M[14] +
           D[78] * M[15];
#pragma omp atomic
  L[19] += D[40] * M[0] + D[44] * M[1] + D[45] * M[2] + D[49] * M[3] + D[50] * M[4] +
           D[51] * M[5] + D[61] * M[6] + D[65] * M[7] + D[66] * M[8] + D[70] * M[9] +
           D[71] * M[10] + D[72] * M[11] + D[76] * M[12] + D[77] * M[13] + D[78] * M[14] +
           D[79] * M[15];
#pragma omp atomic
  L[20] += D[52] * M[0] + D[53] * M[1] + D[54] * M[2] + D[55] * M[3] + D[56] * M[4] +
           D[57] * M[5];
#pragma omp atomic
  L[21] += D[53] * M[0] + D[55] * M[1] + D[56] * M[2] + D[58] * M[3] + D[59] * M[4] +
           D[60] * M[5];
#pragma omp atomic
  L[22] += D[54] * M[0] + D[56] * M[1] + D[57] * M[2] + D[59] * M[3] + D[60] * M[4] +
           D[61] * M[5];
#pragma omp atomic
  L[23] += D[55] * M[0] + D[58] * M[1] + D[59] * M[2] + D[62] * M[3] + D[63] * M[4] +
           D[64] * M[5];
#pragma omp atomic
  L[24] += D[56] * M[0] + D[59] * M[1] + D[60] * M[2] + D[63] * M[3] + D[64] * M[4] +
           D[65] * M[5];
#pragma omp atomic
  L[25] += D[57] * M[0] + D[60] * M[1] + D[61] * M[2] + D[64] * M[3] + D[65] * M[4] +
           D[66] * M[5];
#pragma omp atomic
  L[26] += D[58] * M[0] + D[62] * M[1] + D[63] * M[2] + D[67] * M[3] + D[68] * M[4] +
           D[69] * M[5];
#pragma omp atomic
  L[27] += D[59] * M[0] + D[63] * M[1] + D[64] * M[2] + D[68] * M[3] + D[69] * M[4] +
           D[70] * M[5];
#pragma omp atomic
  L[28] += D[60] * M[0] + D[64] * M[1] + D[65] * M[2] + D[69] * M[3] + D[70] * M[4] +
           D[71] * M[5];
#pragma omp atomic
  L[29] += D[61] * M[0] + D[65] * M[1] + D[66] * M[2] + D[70] * M[3] + D[71] * M[4] +
           D[72] * M[5];
#pragma omp atomic
  L[30] += D[62] * M[0] + D[67] * M[1] + D[68] * M[2] + D[73] * M[3] + D[74] * M[4] +
           D[75] * M[5];
#pragma omp atomic
  L[31] += D[63] * M[0] + D[68] * M[1] + D[69] * M[2] + D[74] * M[3] + D[75] * M[4] +
           D[76] * M[5];
#pragma omp atomic
  L[32] += D[64] * M[0] + D[69] * M[1] + D[70] * M[2] + D[75] * M[3] + D[76] * M[4] +
           D[77] * M[5];
#pragma omp atomic
  L[33] += D[65] * M[0] + D[70] * M[1] + D[71] * M[2] + D[76] * M[3] + D[77] * M[4] +
           D[78] * M[5];
#pragma omp atomic
  L[34] += D[66] * M[0] + D[71] * M[1] + D[72] * M[2] + D[77] * M[3] + D[78] * M[4] +
           D[79] * M[5];
}

void field_m2_L2L_6(double x, double y, double z, double* L, double* Ls) {
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

void field_m2_L2P_6(double x, double y, double z, double* L, double* F) {
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

void field_m2_M2P_6(double x, double y, double z, double* M, double* F) {
  double R       = sqrt(x * x + y * y + z * z);
  double Ftmp0   = pow(R, -5.0);
  double Ftmp1   = 3.0 * Ftmp0;
  double Ftmp2   = Ftmp1 * M[1];
  double Ftmp3   = Ftmp1 * z;
  double Ftmp4   = y * z;
  double Ftmp5   = pow(R, -7.0);
  double Ftmp6   = 15.0 * Ftmp5;
  double Ftmp7   = Ftmp6 * M[10];
  double Ftmp8   = Ftmp4 * x;
  double Ftmp9   = Ftmp6 * Ftmp8;
  double Ftmp10  = (x * x);
  double Ftmp11  = Ftmp10 * Ftmp6;
  double Ftmp12  = pow(R, -9.0);
  double Ftmp13  = 105.0 * Ftmp12;
  double Ftmp14  = Ftmp10 * Ftmp13;
  double Ftmp15  = -9.0 * Ftmp0;
  double Ftmp16  = Ftmp11 + Ftmp15;
  double Ftmp17  = -Ftmp1;
  double Ftmp18  = (y * y);
  double Ftmp19  = Ftmp18 * Ftmp6;
  double Ftmp20  = Ftmp17 + Ftmp19;
  double Ftmp21  = (z * z);
  double Ftmp22  = Ftmp21 * Ftmp6;
  double Ftmp23  = Ftmp17 + Ftmp22;
  double Ftmp24  = Ftmp20 * M[3];
  double Ftmp25  = Ftmp23 * M[5];
  double Ftmp26  = 45.0 * Ftmp5;
  double Ftmp27  = -Ftmp26;
  double Ftmp28  = Ftmp14 + Ftmp27;
  double Ftmp29  = Ftmp28 * M[17];
  double Ftmp30  = 1.0 * y;
  double Ftmp31  = Ftmp13 * Ftmp18;
  double Ftmp32  = Ftmp27 + Ftmp31;
  double Ftmp33  = Ftmp32 * M[22];
  double Ftmp34  = 3.0 * y;
  double Ftmp35  = 35.0 * Ftmp12;
  double Ftmp36  = (Ftmp21 * Ftmp35 - 5.0 * Ftmp5) * M[24];
  double Ftmp37  = Ftmp28 * M[18];
  double Ftmp38  = 1.0 * z;
  double Ftmp39  = -Ftmp6;
  double Ftmp40  = Ftmp31 + Ftmp39;
  double Ftmp41  = Ftmp40 * M[23];
  double Ftmp42  = Ftmp13 * Ftmp21;
  double Ftmp43  = Ftmp27 + Ftmp42;
  double Ftmp44  = Ftmp38 * Ftmp43;
  double Ftmp45  = 315.0 * Ftmp12;
  double Ftmp46  = -Ftmp45;
  double Ftmp47  = pow(R, -11.0);
  double Ftmp48  = 945.0 * Ftmp47;
  double Ftmp49  = Ftmp10 * Ftmp48;
  double Ftmp50  = Ftmp46 + Ftmp49;
  double Ftmp51  = Ftmp50 * M[35];
  double Ftmp52  = x * y;
  double Ftmp53  = Ftmp28 * Ftmp52;
  double Ftmp54  = Ftmp32 * M[12];
  double Ftmp55  = Ftmp39 + Ftmp42;
  double Ftmp56  = Ftmp55 * M[14];
  double Ftmp57  = Ftmp30 * x;
  double Ftmp58  = x * z;
  double Ftmp59  = Ftmp28 * Ftmp58;
  double Ftmp60  = Ftmp40 * M[13];
  double Ftmp61  = Ftmp43 * M[15];
  double Ftmp62  = Ftmp18 * Ftmp48;
  double Ftmp63  = Ftmp46 + Ftmp62;
  double Ftmp64  = Ftmp30 * Ftmp63;
  double Ftmp65  = Ftmp64 * M[42];
  double Ftmp66  = -Ftmp13;
  double Ftmp67  = Ftmp21 * Ftmp47;
  double Ftmp68  = 315.0 * Ftmp67;
  double Ftmp69  = Ftmp66 + Ftmp68;
  double Ftmp70  = Ftmp69 * M[44];
  double Ftmp71  = Ftmp34 * Ftmp70;
  double Ftmp72  = Ftmp63 * M[27];
  double Ftmp73  = -75.0 * Ftmp5;
  double Ftmp74  = 1.0 * Ftmp10;
  double Ftmp75  = Ftmp40 * M[9];
  double Ftmp76  = Ftmp55 * M[11];
  double Ftmp77  = 525.0 * Ftmp12;
  double Ftmp78  = -Ftmp77;
  double Ftmp79  = Ftmp49 + Ftmp78;
  double Ftmp80  = Ftmp10 * y;
  double Ftmp81  = Ftmp10 * z;
  double Ftmp82  = Ftmp21 * Ftmp48;
  double Ftmp83  = Ftmp46 + Ftmp82;
  double Ftmp84  = Ftmp30 * M[29];
  double Ftmp85  = Ftmp64 * M[22];
  double Ftmp86  = Ftmp10 * Ftmp34;
  double Ftmp87  = (-Ftmp35 + Ftmp68) * M[24];
  double Ftmp88  = Ftmp10 * Ftmp38;
  double Ftmp89  = (Ftmp62 + Ftmp66) * M[23];
  double Ftmp90  = Ftmp83 * M[25];
  double Ftmp91  = 4725.0 * Ftmp47;
  double Ftmp92  = -Ftmp91;
  double Ftmp93  = pow(R, -13.0);
  double Ftmp94  = 10395.0 * Ftmp93;
  double Ftmp95  = Ftmp10 * Ftmp94;
  double Ftmp96  = 2835.0 * Ftmp47;
  double Ftmp97  = -Ftmp96;
  double Ftmp98  = Ftmp18 * Ftmp94;
  double Ftmp99  = Ftmp30 * (Ftmp97 + Ftmp98) * M[42];
  double Ftmp100 = 3465.0 * Ftmp93;
  double Ftmp101 = Ftmp100 * Ftmp21;
  double Ftmp102 = (Ftmp101 - Ftmp48) * M[44];
  double Ftmp103 = 225.0 * Ftmp5;
  double Ftmp104 = (x * x * x * x);
  double Ftmp105 = Ftmp104 * Ftmp48;
  double Ftmp106 = 1050.0 * Ftmp12;
  double Ftmp107 = -Ftmp10 * Ftmp106 + Ftmp103 + Ftmp105;
  double Ftmp108 = (y * y * y * y);
  double Ftmp109 = Ftmp108 * Ftmp48;
  double Ftmp110 = 630.0 * Ftmp12;
  double Ftmp111 = Ftmp109 - Ftmp110 * Ftmp18 + Ftmp26;
  double Ftmp112 = (z * z * z * z);
  double Ftmp113 = Ftmp112 * Ftmp48;
  double Ftmp114 = -Ftmp110 * Ftmp21 + Ftmp113 + Ftmp26;
  double Ftmp115 = Ftmp111 * M[26];
  double Ftmp116 = Ftmp114 * M[30];
  double Ftmp117 = 1575.0 * Ftmp12;
  double Ftmp118 = Ftmp104 * Ftmp93;
  double Ftmp119 = 10395.0 * Ftmp118;
  double Ftmp120 = Ftmp10 * Ftmp47;
  double Ftmp121 = 9450.0 * Ftmp120;
  double Ftmp122 = Ftmp117 + Ftmp119 - Ftmp121;
  double Ftmp123 = Ftmp122 * M[53];
  double Ftmp124 = Ftmp108 * Ftmp94;
  double Ftmp125 = 9450.0 * Ftmp47;
  double Ftmp126 = Ftmp125 * Ftmp18;
  double Ftmp127 = Ftmp117 + Ftmp124 - Ftmp126;
  double Ftmp128 = Ftmp127 * M[67];
  double Ftmp129 = (Ftmp100 * Ftmp112 + Ftmp13 - 1890.0 * Ftmp67) * M[71];
  double Ftmp130 = Ftmp122 * M[54];
  double Ftmp131 = 5670.0 * Ftmp47;
  double Ftmp132 = Ftmp131 * Ftmp18;
  double Ftmp133 = Ftmp124 - Ftmp132 + Ftmp45;
  double Ftmp134 = Ftmp133 * M[68];
  double Ftmp135 = Ftmp112 * Ftmp94;
  double Ftmp136 = Ftmp125 * Ftmp21;
  double Ftmp137 = Ftmp117 + Ftmp135 - Ftmp136;
  double Ftmp138 = Ftmp137 * Ftmp38;
  double Ftmp139 = Ftmp122 * Ftmp52;
  double Ftmp140 = Ftmp127 * M[46];
  double Ftmp141 = Ftmp131 * Ftmp21;
  double Ftmp142 = Ftmp135 - Ftmp141 + Ftmp45;
  double Ftmp143 = Ftmp142 * M[50];
  double Ftmp144 = Ftmp122 * Ftmp58;
  double Ftmp145 = Ftmp133 * M[47];
  double Ftmp146 = Ftmp137 * M[51];
  double Ftmp147 = 14175.0 * Ftmp47;
  double Ftmp148 = pow(R, -15.0);
  double Ftmp149 = 135135.0 * Ftmp148;
  double Ftmp150 = Ftmp104 * Ftmp149;
  double Ftmp151 = Ftmp10 * Ftmp93;
  double Ftmp152 = 103950.0 * Ftmp151;
  double Ftmp153 = Ftmp147 + Ftmp150 - Ftmp152;
  double Ftmp154 = Ftmp108 * Ftmp149;
  double Ftmp155 = Ftmp18 * Ftmp93;
  double Ftmp156 = 103950.0 * Ftmp155;
  double Ftmp157 = Ftmp147 + Ftmp154 - Ftmp156;
  double Ftmp158 = Ftmp157 * M[74];
  double Ftmp159 = 3675.0 * Ftmp12;
  double Ftmp160 = Ftmp133 * M[41];
  double Ftmp161 = Ftmp142 * M[45];
  double Ftmp162 = 33075.0 * Ftmp47;
  double Ftmp163 = Ftmp150 - 145530.0 * Ftmp151 + Ftmp162;
  double Ftmp164 = Ftmp112 * Ftmp149;
  double Ftmp165 = Ftmp21 * Ftmp93;
  double Ftmp166 = Ftmp147 + Ftmp164 - 103950.0 * Ftmp165;
  double Ftmp167 = Ftmp30 * M[78];
  double Ftmp168 = Ftmp10 * Ftmp30;
  double Ftmp169 = Ftmp157 * M[67];
  double Ftmp170 = 45045.0 * Ftmp112 * Ftmp148;
  double Ftmp171 = (-20790.0 * Ftmp165 + Ftmp170 + Ftmp48) * M[71];
  double Ftmp172 = 62370.0 * Ftmp155;
  double Ftmp173 = (Ftmp154 - Ftmp172 + Ftmp96) * M[68];
  double Ftmp174 = Ftmp166 * M[72];
  double Ftmp175 = -11025.0 * Ftmp12;
  double Ftmp176 = Ftmp149 * (x * x * x * x * x * x);
  double Ftmp177 = -Ftmp117;
  double Ftmp178 = Ftmp149 * (y * y * y * y * y * y);
  double Ftmp179 = 155925.0 * Ftmp93;
  double Ftmp180 = 42525.0 * Ftmp47;
  double Ftmp181 = (-Ftmp108 * Ftmp179 + Ftmp177 + Ftmp178 + Ftmp18 * Ftmp180) * M[73];
  double Ftmp182 = Ftmp149 * (z * z * z * z * z * z);
  double Ftmp183 = (-Ftmp112 * Ftmp179 + Ftmp177 + Ftmp180 * Ftmp21 + Ftmp182) * M[79];
  double Ftmp184 = -Ftmp18 * Ftmp45;
  double Ftmp185 = Ftmp18 * Ftmp49;
  double Ftmp186 = -Ftmp14;
  double Ftmp187 = Ftmp186 + Ftmp26;
  double Ftmp188 = Ftmp184 + Ftmp185 + Ftmp187;
  double Ftmp189 = -Ftmp21 * Ftmp45;
  double Ftmp190 = Ftmp21 * Ftmp49;
  double Ftmp191 = Ftmp187 + Ftmp189 + Ftmp190;
  double Ftmp192 = -Ftmp42;
  double Ftmp193 = Ftmp192 + Ftmp6;
  double Ftmp194 = -Ftmp31;
  double Ftmp195 = Ftmp21 * Ftmp62;
  double Ftmp196 = Ftmp194 + Ftmp195;
  double Ftmp197 = Ftmp193 + Ftmp196;
  double Ftmp198 = Ftmp18 * Ftmp96;
  double Ftmp199 = -Ftmp198;
  double Ftmp200 = Ftmp18 * Ftmp95;
  double Ftmp201 = Ftmp199 + Ftmp200;
  double Ftmp202 = 945.0 * Ftmp12;
  double Ftmp203 = Ftmp10 * Ftmp96;
  double Ftmp204 = -Ftmp203;
  double Ftmp205 = Ftmp202 + Ftmp204;
  double Ftmp206 = Ftmp201 + Ftmp205;
  double Ftmp207 = Ftmp206 * M[58];
  double Ftmp208 = -Ftmp49;
  double Ftmp209 = Ftmp208 + Ftmp45;
  double Ftmp210 = Ftmp21 * Ftmp96;
  double Ftmp211 = -Ftmp210;
  double Ftmp212 = Ftmp21 * Ftmp95;
  double Ftmp213 = Ftmp211 + Ftmp212;
  double Ftmp214 = Ftmp209 + Ftmp213;
  double Ftmp215 = Ftmp214 * M[60];
  double Ftmp216 = -Ftmp62;
  double Ftmp217 = Ftmp21 * Ftmp98;
  double Ftmp218 = Ftmp217 + Ftmp45;
  double Ftmp219 = Ftmp211 + Ftmp216 + Ftmp218;
  double Ftmp220 = Ftmp219 * M[69];
  double Ftmp221 = Ftmp201 + Ftmp209;
  double Ftmp222 = Ftmp221 * M[59];
  double Ftmp223 = Ftmp205 + Ftmp213;
  double Ftmp224 = Ftmp223 * M[61];
  double Ftmp225 = -Ftmp82;
  double Ftmp226 = Ftmp199 + Ftmp218 + Ftmp225;
  double Ftmp227 = Ftmp226 * M[70];
  double Ftmp228 = Ftmp206 * Ftmp52;
  double Ftmp229 = Ftmp214 * Ftmp52;
  double Ftmp230 = Ftmp221 * Ftmp58;
  double Ftmp231 = Ftmp223 * Ftmp58;
  double Ftmp232 = 31185.0 * Ftmp93;
  double Ftmp233 = Ftmp10 * Ftmp232;
  double Ftmp234 = -Ftmp233;
  double Ftmp235 = 8505.0 * Ftmp47;
  double Ftmp236 = Ftmp234 + Ftmp235;
  double Ftmp237 = Ftmp10 * Ftmp149;
  double Ftmp238 = Ftmp18 * Ftmp237;
  double Ftmp239 = -Ftmp18 * Ftmp232;
  double Ftmp240 = Ftmp238 + Ftmp239;
  double Ftmp241 = Ftmp8 * (Ftmp236 + Ftmp240);
  double Ftmp242 = Ftmp21 * Ftmp237;
  double Ftmp243 = Ftmp21 * Ftmp232;
  double Ftmp244 = -Ftmp243;
  double Ftmp245 = Ftmp242 + Ftmp244;
  double Ftmp246 = Ftmp8 * (Ftmp236 + Ftmp245);
  double Ftmp247 = Ftmp149 * Ftmp18 * Ftmp21;
  double Ftmp248 = Ftmp239 + Ftmp247;
  double Ftmp249 = Ftmp235 + Ftmp244 + Ftmp248;
  double Ftmp250 = -Ftmp18 * Ftmp91;
  double Ftmp251 = Ftmp208 + Ftmp77;
  double Ftmp252 = -Ftmp21 * Ftmp91;
  double Ftmp253 = Ftmp13 + Ftmp225;
  double Ftmp254 = Ftmp216 + Ftmp217;
  double Ftmp255 = 51975.0 * Ftmp93;
  double Ftmp256 = -Ftmp18 * Ftmp255;
  double Ftmp257 = Ftmp238 + Ftmp256;
  double Ftmp258 = Ftmp147 + Ftmp234;
  double Ftmp259 = -Ftmp95;
  double Ftmp260 = Ftmp259 + Ftmp91;
  double Ftmp261 = -Ftmp21 * Ftmp255;
  double Ftmp262 = Ftmp242 + Ftmp261;
  double Ftmp263 = -Ftmp98;
  double Ftmp264 = Ftmp247 + Ftmp96;
  double Ftmp265 = -Ftmp21 * Ftmp94;
  double Ftmp266 = -Ftmp202;
  double Ftmp267 = Ftmp203 + Ftmp266;
  double Ftmp268 = Ftmp10 * Ftmp154;
  double Ftmp269 = 62370.0 * Ftmp151;
  double Ftmp270 = -Ftmp18 * Ftmp269;
  double Ftmp271 = Ftmp268 + Ftmp270;
  double Ftmp272 = 17010.0 * Ftmp47;
  double Ftmp273 = -Ftmp108 * Ftmp232 + Ftmp18 * Ftmp272;
  double Ftmp274 = Ftmp10 * Ftmp164;
  double Ftmp275 = -Ftmp21 * Ftmp269;
  double Ftmp276 = Ftmp274 + Ftmp275;
  double Ftmp277 = -Ftmp112 * Ftmp232 + Ftmp21 * Ftmp272;
  double Ftmp278 = Ftmp147 * Ftmp18;
  double Ftmp279 = -Ftmp152 * Ftmp18 + Ftmp177;
  double Ftmp280 = -Ftmp119;
  double Ftmp281 = Ftmp150 * Ftmp18;
  double Ftmp282 = Ftmp280 + Ftmp281;
  double Ftmp283 = -Ftmp152 * Ftmp21;
  double Ftmp284 = Ftmp150 * Ftmp21;
  double Ftmp285 = Ftmp280 + Ftmp284;
  double Ftmp286 = Ftmp147 * Ftmp21 + Ftmp177;
  double Ftmp287 = -Ftmp172 * Ftmp21;
  double Ftmp288 = Ftmp287 + Ftmp46;
  double Ftmp289 = -Ftmp135;
  double Ftmp290 = Ftmp164 * Ftmp18;
  double Ftmp291 = Ftmp289 + Ftmp290;
  double Ftmp292 = -Ftmp124;
  double Ftmp293 = Ftmp154 * Ftmp21;
  double Ftmp294 = Ftmp292 + Ftmp293;
  double Ftmp295 = Ftmp18 * Ftmp242;
  double Ftmp296 = -Ftmp200 + Ftmp210 + Ftmp295;
  double Ftmp297 = Ftmp198 - Ftmp212;
  double Ftmp298 = x * M[2];
  double Ftmp299 = Ftmp11 + Ftmp17;
  double Ftmp300 = Ftmp15 + Ftmp19;
  double Ftmp301 = Ftmp299 * M[0];
  double Ftmp302 = 1.0 * x;
  double Ftmp303 = 3.0 * x;
  double Ftmp304 = Ftmp14 + Ftmp39;
  double Ftmp305 = Ftmp304 * M[20];
  double Ftmp306 = Ftmp32 * M[27];
  double Ftmp307 = Ftmp38 * x;
  double Ftmp308 = 3.0 * Ftmp58;
  double Ftmp309 = Ftmp304 * M[8];
  double Ftmp310 = Ftmp50 * M[18];
  double Ftmp311 = Ftmp304 * M[7];
  double Ftmp312 = 1.0 * Ftmp18;
  double Ftmp313 = Ftmp18 * x;
  double Ftmp314 = Ftmp50 * M[17];
  double Ftmp315 = Ftmp18 * z;
  double Ftmp316 = (Ftmp49 + Ftmp66) * M[20];
  double Ftmp317 = Ftmp62 + Ftmp78;
  double Ftmp318 = Ftmp30 * Ftmp58;
  double Ftmp319 = Ftmp18 * Ftmp302;
  double Ftmp320 = Ftmp18 * Ftmp303;
  double Ftmp321 = Ftmp18 * Ftmp38;
  double Ftmp322 = (Ftmp95 + Ftmp97) * M[35];
  double Ftmp323 = -Ftmp10 * Ftmp110 + Ftmp105 + Ftmp26;
  double Ftmp324 = Ftmp103 - Ftmp106 * Ftmp18 + Ftmp109;
  double Ftmp325 = Ftmp323 * M[16];
  double Ftmp326 = 5670.0 * Ftmp120;
  double Ftmp327 = Ftmp119 - Ftmp326 + Ftmp45;
  double Ftmp328 = Ftmp327 * M[56];
  double Ftmp329 = Ftmp127 * M[74];
  double Ftmp330 = Ftmp327 * M[33];
  double Ftmp331 = Ftmp153 * M[54];
  double Ftmp332 = Ftmp327 * M[32];
  double Ftmp333 = Ftmp18 * Ftmp47;
  double Ftmp334 = Ftmp153 * M[53];
  double Ftmp335 = (Ftmp150 - Ftmp269 + Ftmp96) * M[56];
  double Ftmp336 = Ftmp154 - 145530.0 * Ftmp155 + Ftmp162;
  double Ftmp337 = (Ftmp10 * Ftmp180 - Ftmp104 * Ftmp179 + Ftmp176 + Ftmp177) * M[52];
  double Ftmp338 = 218295.0 * Ftmp93;
  double Ftmp339 = Ftmp185 + Ftmp194;
  double Ftmp340 = -Ftmp10 * Ftmp45 + Ftmp26;
  double Ftmp341 = Ftmp339 + Ftmp340;
  double Ftmp342 = Ftmp186 + Ftmp190 + Ftmp193;
  double Ftmp343 = Ftmp189 + Ftmp196 + Ftmp26;
  double Ftmp344 = Ftmp204 + Ftmp45;
  double Ftmp345 = Ftmp200 + Ftmp216;
  double Ftmp346 = Ftmp344 + Ftmp345;
  double Ftmp347 = Ftmp346 * M[63];
  double Ftmp348 = Ftmp212 + Ftmp225;
  double Ftmp349 = Ftmp344 + Ftmp348;
  double Ftmp350 = Ftmp349 * M[65];
  double Ftmp351 = Ftmp199 + Ftmp202 + Ftmp211 + Ftmp217;
  double Ftmp352 = Ftmp351 * M[76];
  double Ftmp353 = Ftmp346 * Ftmp4;
  double Ftmp354 = Ftmp349 * Ftmp4;
  double Ftmp355 = Ftmp351 * Ftmp4;
  double Ftmp356 = -Ftmp10 * Ftmp91 + Ftmp77;
  double Ftmp357 = -51975.0 * Ftmp151;
  double Ftmp358 = Ftmp147 + Ftmp357;
  double Ftmp359 = Ftmp259 + Ftmp96;
  double Ftmp360 = Ftmp263 + Ftmp91;
  double Ftmp361 = Ftmp234 + Ftmp96;
  double Ftmp362 = Ftmp242 + Ftmp265;
  double Ftmp363 = Ftmp249 * Ftmp318;
  double Ftmp364 = Ftmp10 * Ftmp147;
  double Ftmp365 = Ftmp203 + Ftmp46;
  double Ftmp366 = Ftmp198 + Ftmp266;
  double Ftmp367 = -31185.0 * Ftmp118 + 17010.0 * Ftmp120;
  double Ftmp368 = Ftmp326 + Ftmp46;
  double Ftmp369 = Ftmp210 + Ftmp275;
  double Ftmp370 = -Ftmp156 * Ftmp21;
  double Ftmp371 = Ftmp203 - Ftmp217;
  double Ftmp372 = y * M[4];
  double Ftmp373 = Ftmp15 + Ftmp22;
  double Ftmp374 = Ftmp302 * M[25];
  double Ftmp375 = Ftmp30 * z;
  double Ftmp376 = Ftmp21 * x;
  double Ftmp377 = Ftmp21 * y;
  double Ftmp378 = Ftmp34 * Ftmp58;
  double Ftmp379 = Ftmp21 * Ftmp302;
  double Ftmp380 = Ftmp21 * (Ftmp78 + Ftmp82);
  double Ftmp381 = Ftmp103 - Ftmp106 * Ftmp21 + Ftmp113;
  double Ftmp382 = Ftmp302 * M[72];
  double Ftmp383 = Ftmp21 * (Ftmp162 + Ftmp164 - 145530.0 * Ftmp165);
  double Ftmp384 = Ftmp186 + Ftmp339 + Ftmp6;
  double Ftmp385 = Ftmp190 + Ftmp192 + Ftmp340;
  double Ftmp386 = Ftmp184 + Ftmp192 + Ftmp195 + Ftmp26;
  double Ftmp387 = Ftmp247 + Ftmp256;
  double Ftmp388 = Ftmp136 + Ftmp177;
#pragma omp atomic
  F[0] += -Ftmp10 * Ftmp4 * (Ftmp92 + Ftmp95) * M[35] + Ftmp10 * Ftmp85 -
          Ftmp10 * (Ftmp14 + Ftmp73) * M[6] -
          Ftmp10 * (Ftmp119 - 13230.0 * Ftmp120 + Ftmp159) * M[31] -
          Ftmp10 * (Ftmp200 + Ftmp250 + Ftmp251) * M[34] -
          Ftmp10 * (Ftmp212 + Ftmp251 + Ftmp252) * M[36] - Ftmp102 * Ftmp34 * Ftmp81 +
          Ftmp107 * x * M[16] + Ftmp107 * M[31] + Ftmp11 * y * M[1] + Ftmp11 * z * M[2] +
          Ftmp111 * M[41] + Ftmp114 * M[45] + Ftmp115 * x + Ftmp116 * x - Ftmp123 * y -
          Ftmp128 * Ftmp30 - Ftmp129 * Ftmp34 - Ftmp130 * z - Ftmp134 * Ftmp38 -
          Ftmp138 * M[72] - Ftmp139 * M[32] - Ftmp14 * Ftmp4 * M[10] - Ftmp140 * Ftmp52 -
          Ftmp143 * Ftmp57 - Ftmp144 * M[33] - Ftmp145 * Ftmp58 - Ftmp146 * Ftmp58 +
          Ftmp153 * Ftmp8 * M[56] + Ftmp158 * Ftmp8 + Ftmp16 * x * M[0] + Ftmp16 * M[6] -
          Ftmp160 * Ftmp74 - Ftmp161 * Ftmp74 + Ftmp163 * Ftmp80 * M[53] +
          Ftmp163 * Ftmp81 * M[54] + Ftmp166 * Ftmp167 * Ftmp58 + Ftmp168 * Ftmp169 +
          Ftmp168 * (Ftmp244 + Ftmp263 + Ftmp264) * M[69] + Ftmp171 * Ftmp86 +
          Ftmp173 * Ftmp88 + Ftmp174 * Ftmp88 + Ftmp181 * x + Ftmp183 * x +
          Ftmp188 * x * M[19] + Ftmp188 * M[34] + Ftmp191 * x * M[21] + Ftmp191 * M[36] +
          Ftmp197 * x * M[28] + Ftmp197 * M[43] - Ftmp2 * y + Ftmp20 * M[9] -
          Ftmp207 * y - Ftmp215 * y - Ftmp219 * Ftmp52 * M[48] - Ftmp220 * Ftmp30 -
          Ftmp222 * z - Ftmp224 * z - Ftmp226 * Ftmp58 * M[49] - Ftmp227 * Ftmp38 -
          Ftmp228 * M[37] - Ftmp229 * M[39] + Ftmp23 * M[11] - Ftmp230 * M[38] -
          Ftmp231 * M[40] + Ftmp24 * x + Ftmp241 * M[63] + Ftmp246 * M[65] +
          Ftmp249 * Ftmp8 * M[76] + Ftmp25 * x - Ftmp29 * y - Ftmp3 * M[2] -
          Ftmp30 * Ftmp33 - Ftmp34 * Ftmp36 - Ftmp37 * z - Ftmp38 * Ftmp41 +
          Ftmp4 * Ftmp51 + Ftmp4 * Ftmp7 - Ftmp44 * M[25] + Ftmp50 * Ftmp8 * M[20] -
          Ftmp52 * Ftmp54 - Ftmp53 * M[7] - Ftmp56 * Ftmp57 - Ftmp58 * Ftmp60 -
          Ftmp58 * Ftmp61 + Ftmp58 * Ftmp83 * Ftmp84 - Ftmp59 * M[8] + Ftmp65 * z +
          Ftmp71 * z + Ftmp72 * Ftmp8 - Ftmp74 * Ftmp75 - Ftmp74 * Ftmp76 -
          Ftmp74 * (Ftmp253 + Ftmp254) * M[43] + Ftmp79 * Ftmp80 * M[17] +
          Ftmp79 * Ftmp81 * M[18] + Ftmp80 * (Ftmp257 + Ftmp258) * M[58] +
          Ftmp80 * (Ftmp260 + Ftmp262) * M[60] - Ftmp81 * Ftmp99 +
          Ftmp81 * (Ftmp257 + Ftmp260) * M[59] + Ftmp81 * (Ftmp258 + Ftmp262) * M[61] +
          Ftmp86 * Ftmp87 + Ftmp88 * Ftmp89 + Ftmp88 * Ftmp90 +
          Ftmp88 * (Ftmp239 + Ftmp264 + Ftmp265) * M[70] + Ftmp9 * M[4] +
          x * (Ftmp267 + Ftmp271 + Ftmp273) * M[62] +
          x * (Ftmp267 + Ftmp276 + Ftmp277) * M[66] +
          x * (-218295.0 * Ftmp118 + 99225.0 * Ftmp120 + Ftmp175 + Ftmp176) * M[52] +
          x * (Ftmp121 + Ftmp278 + Ftmp279 + Ftmp282) * M[55] +
          x * (Ftmp121 + Ftmp283 + Ftmp285 + Ftmp286) * M[57] +
          x * (Ftmp132 + Ftmp210 + Ftmp288 + Ftmp294) * M[75] +
          x * (Ftmp141 + Ftmp198 + Ftmp288 + Ftmp291) * M[77] +
          x * (-Ftmp18 * Ftmp243 + Ftmp296 + Ftmp297 + Ftmp50) * M[64];
#pragma omp atomic
  F[1] += -Ftmp102 * Ftmp18 * Ftmp308 + Ftmp114 * M[50] + Ftmp116 * y - Ftmp123 * x -
          Ftmp127 * Ftmp4 * M[47] - Ftmp127 * Ftmp57 * M[41] - Ftmp128 * Ftmp302 -
          Ftmp129 * Ftmp303 - Ftmp138 * M[78] - Ftmp139 * M[31] - Ftmp143 * Ftmp312 -
          Ftmp146 * Ftmp4 + Ftmp157 * Ftmp318 * M[68] - Ftmp161 * Ftmp57 +
          Ftmp166 * Ftmp321 * M[78] + Ftmp171 * Ftmp320 + Ftmp174 * Ftmp318 -
          Ftmp18 * Ftmp311 - Ftmp18 * Ftmp322 * Ftmp58 - Ftmp18 * Ftmp332 -
          Ftmp18 * (Ftmp31 + Ftmp73) * M[12] - Ftmp18 * (Ftmp345 + Ftmp356) * M[37] -
          Ftmp18 * (Ftmp124 + Ftmp159 - 13230.0 * Ftmp333) * M[46] -
          Ftmp18 * (Ftmp208 + Ftmp212 + Ftmp253) * M[39] -
          Ftmp18 * (Ftmp252 + Ftmp254 + Ftmp77) * M[48] + Ftmp183 * y +
          Ftmp19 * x * M[1] + Ftmp19 * z * M[4] - Ftmp2 * x - Ftmp207 * x - Ftmp215 * x -
          Ftmp219 * Ftmp57 * M[43] - Ftmp220 * Ftmp302 - Ftmp228 * M[34] -
          Ftmp229 * M[36] + Ftmp23 * M[14] + Ftmp241 * M[59] + Ftmp246 * M[61] +
          Ftmp25 * y - Ftmp29 * x + Ftmp298 * Ftmp4 * Ftmp6 + Ftmp299 * M[7] -
          Ftmp3 * M[4] + Ftmp300 * y * M[3] + Ftmp300 * M[12] + Ftmp301 * y -
          Ftmp302 * Ftmp33 - Ftmp303 * Ftmp36 - Ftmp305 * z - Ftmp306 * z +
          Ftmp307 * Ftmp63 * M[42] + Ftmp308 * Ftmp70 - Ftmp309 * Ftmp4 -
          Ftmp31 * Ftmp58 * M[10] + Ftmp310 * Ftmp8 - Ftmp312 * Ftmp56 +
          Ftmp313 * Ftmp314 + Ftmp313 * Ftmp334 -
          Ftmp313 * Ftmp38 * (Ftmp92 + Ftmp98) * M[42] +
          Ftmp313 * (Ftmp240 + Ftmp358) * M[58] + Ftmp313 * (Ftmp245 + Ftmp359) * M[60] +
          Ftmp315 * Ftmp316 + Ftmp315 * Ftmp317 * M[27] + Ftmp315 * Ftmp335 +
          Ftmp315 * Ftmp336 * M[74] + Ftmp315 * (Ftmp361 + Ftmp362) * M[65] +
          Ftmp315 * (Ftmp147 + Ftmp248 + Ftmp261) * M[76] +
          Ftmp315 * (Ftmp238 + Ftmp357 + Ftmp360) * M[63] + Ftmp317 * Ftmp319 * M[22] +
          Ftmp318 * Ftmp90 + Ftmp319 * Ftmp336 * M[67] +
          Ftmp319 * (Ftmp247 + Ftmp261 + Ftmp360) * M[69] - Ftmp32 * Ftmp4 * M[13] -
          Ftmp32 * Ftmp57 * M[9] + Ftmp320 * Ftmp87 + Ftmp321 * Ftmp83 * M[29] +
          Ftmp323 * M[32] + Ftmp324 * y * M[26] + Ftmp324 * M[46] + Ftmp325 * y -
          Ftmp328 * z - Ftmp329 * z - Ftmp330 * Ftmp4 + Ftmp331 * Ftmp8 + Ftmp337 * y +
          Ftmp341 * y * M[19] + Ftmp341 * M[37] + Ftmp342 * y * M[21] + Ftmp342 * M[39] +
          Ftmp343 * y * M[28] + Ftmp343 * M[48] - Ftmp347 * z - Ftmp350 * z -
          Ftmp352 * z - Ftmp353 * M[38] - Ftmp354 * M[40] - Ftmp355 * M[49] +
          Ftmp363 * M[70] - Ftmp4 * Ftmp61 - Ftmp44 * M[29] + Ftmp51 * Ftmp58 -
          Ftmp53 * M[6] - Ftmp57 * Ftmp76 + Ftmp58 * Ftmp64 * M[23] + Ftmp58 * Ftmp7 +
          y * (Ftmp285 + Ftmp368 + Ftmp369) * M[57] +
          y * (Ftmp126 + Ftmp286 + Ftmp294 + Ftmp370) * M[75] +
          y * (Ftmp141 + Ftmp276 + Ftmp289 + Ftmp365) * M[66] +
          y * (Ftmp270 + Ftmp281 + Ftmp366 + Ftmp367) * M[55] +
          y * (Ftmp277 + Ftmp287 + Ftmp290 + Ftmp366) * M[77] +
          y * (-Ftmp108 * Ftmp338 + Ftmp175 + Ftmp178 + 99225.0 * Ftmp333) * M[73] +
          y * (-Ftmp21 * Ftmp233 + Ftmp296 + Ftmp371 + Ftmp63) * M[64] +
          y * (Ftmp126 + Ftmp268 + Ftmp279 + Ftmp292 + Ftmp364) * M[62];
#pragma omp atomic
  F[2] += -Ftmp1 * Ftmp298 - Ftmp1 * Ftmp372 + Ftmp111 * M[47] + Ftmp115 * z -
          Ftmp130 * x - Ftmp134 * Ftmp302 - Ftmp137 * Ftmp167 -
          Ftmp137 * Ftmp375 * M[50] - Ftmp137 * Ftmp382 - Ftmp138 * x * M[45] -
          Ftmp140 * Ftmp4 - Ftmp144 * M[31] - Ftmp145 * Ftmp21 + Ftmp158 * Ftmp377 -
          Ftmp160 * Ftmp307 + Ftmp167 * Ftmp383 + Ftmp169 * Ftmp318 + Ftmp173 * Ftmp379 +
          Ftmp181 * z + Ftmp20 * M[13] - Ftmp21 * Ftmp309 - Ftmp21 * Ftmp322 * Ftmp52 -
          Ftmp21 * Ftmp330 - Ftmp21 * Ftmp60 - Ftmp21 * (Ftmp348 + Ftmp356) * M[40] -
          Ftmp21 * (Ftmp42 + Ftmp73) * M[15] -
          Ftmp21 * (Ftmp13 + Ftmp208 + Ftmp345) * M[38] -
          Ftmp21 * (Ftmp135 + Ftmp159 - 13230.0 * Ftmp67) * M[51] -
          Ftmp21 * (Ftmp217 + Ftmp225 + Ftmp250 + Ftmp77) * M[49] + Ftmp22 * Ftmp298 +
          Ftmp22 * Ftmp372 - Ftmp222 * x - Ftmp224 * x - Ftmp226 * Ftmp307 * M[43] -
          Ftmp227 * Ftmp302 - Ftmp230 * M[34] - Ftmp231 * M[36] + Ftmp24 * z +
          Ftmp241 * M[58] + Ftmp246 * M[60] + Ftmp299 * M[8] + Ftmp301 * z -
          Ftmp302 * Ftmp41 - Ftmp305 * y - Ftmp306 * y - Ftmp307 * Ftmp75 +
          Ftmp310 * Ftmp376 - Ftmp311 * Ftmp4 + Ftmp314 * Ftmp8 + Ftmp316 * Ftmp377 +
          Ftmp323 * M[33] + Ftmp325 * z - Ftmp328 * y - Ftmp329 * y + Ftmp331 * Ftmp376 -
          Ftmp332 * Ftmp4 + Ftmp334 * Ftmp8 + Ftmp335 * Ftmp377 + Ftmp337 * z -
          Ftmp34 * Ftmp376 * (Ftmp101 - 1575.0 * Ftmp47) * M[44] - Ftmp347 * y -
          Ftmp350 * y - Ftmp352 * y - Ftmp353 * M[37] - Ftmp354 * M[39] -
          Ftmp355 * M[48] + Ftmp363 * M[69] - Ftmp37 * x + Ftmp373 * z * M[5] +
          Ftmp373 * M[15] + Ftmp374 * Ftmp380 - Ftmp374 * Ftmp43 -
          Ftmp375 * Ftmp43 * M[14] - Ftmp376 * Ftmp99 +
          Ftmp376 * (Ftmp240 + Ftmp359) * M[59] + Ftmp376 * (Ftmp245 + Ftmp358) * M[61] +
          Ftmp377 * Ftmp72 + Ftmp377 * (Ftmp147 + Ftmp244 + Ftmp387) * M[76] +
          Ftmp377 * (Ftmp238 + Ftmp263 + Ftmp361) * M[63] +
          Ftmp377 * (Ftmp357 + Ftmp362 + Ftmp91) * M[65] + Ftmp378 * Ftmp69 * M[24] +
          Ftmp378 * (-34650.0 * Ftmp165 + Ftmp170 + Ftmp91) * M[71] + Ftmp379 * Ftmp89 +
          Ftmp379 * (Ftmp265 + Ftmp387 + Ftmp91) * M[70] + Ftmp380 * Ftmp84 +
          Ftmp381 * z * M[30] + Ftmp381 * M[51] + Ftmp382 * Ftmp383 +
          Ftmp384 * z * M[19] + Ftmp384 * M[38] + Ftmp385 * z * M[21] + Ftmp385 * M[40] +
          Ftmp386 * z * M[28] + Ftmp386 * M[49] - Ftmp4 * Ftmp54 -
          Ftmp42 * Ftmp52 * M[10] - Ftmp43 * Ftmp84 - Ftmp44 * x * M[11] +
          Ftmp51 * Ftmp52 + Ftmp52 * Ftmp7 + Ftmp58 * Ftmp85 - Ftmp59 * M[6] +
          Ftmp65 * x + Ftmp71 * x + Ftmp9 * M[1] +
          z * (Ftmp132 + Ftmp271 + Ftmp292 + Ftmp365) * M[62] +
          z * (Ftmp198 + Ftmp270 + Ftmp282 + Ftmp368) * M[55] +
          z * (Ftmp266 + Ftmp284 + Ftmp367 + Ftmp369) * M[57] +
          z * (Ftmp278 + Ftmp291 + Ftmp370 + Ftmp388) * M[77] +
          z * (-Ftmp112 * Ftmp338 + Ftmp175 + Ftmp182 + 99225.0 * Ftmp67) * M[79] +
          z * (Ftmp210 + Ftmp266 + Ftmp273 + Ftmp287 + Ftmp293) * M[75] +
          z * (Ftmp274 + Ftmp283 + Ftmp289 + Ftmp364 + Ftmp388) * M[66] +
          z * (-Ftmp18 * Ftmp233 + Ftmp295 + Ftmp297 + Ftmp371 + Ftmp83) * M[64];
}

void field_m2_P2M_7(double x, double y, double z, double q, double* M) {
  double Mtmp0  = (x * x);
  double Mtmp1  = (1.0 / 2.0) * q;
  double Mtmp2  = Mtmp0 * Mtmp1;
  double Mtmp3  = q * x;
  double Mtmp4  = Mtmp3 * y;
  double Mtmp5  = Mtmp3 * z;
  double Mtmp6  = (y * y);
  double Mtmp7  = Mtmp1 * Mtmp6;
  double Mtmp8  = q * y;
  double Mtmp9  = Mtmp8 * z;
  double Mtmp10 = (z * z);
  double Mtmp11 = Mtmp1 * Mtmp10;
  double Mtmp12 = (x * x * x);
  double Mtmp13 = (1.0 / 6.0) * q;
  double Mtmp14 = Mtmp12 * Mtmp13;
  double Mtmp15 = Mtmp2 * y;
  double Mtmp16 = Mtmp7 * x;
  double Mtmp17 = Mtmp11 * x;
  double Mtmp18 = (y * y * y);
  double Mtmp19 = Mtmp13 * Mtmp18;
  double Mtmp20 = (z * z * z);
  double Mtmp21 = (x * x * x * x);
  double Mtmp22 = (1.0 / 24.0) * q;
  double Mtmp23 = Mtmp21 * Mtmp22;
  double Mtmp24 = (1.0 / 6.0) * Mtmp8;
  double Mtmp25 = Mtmp6 * q;
  double Mtmp26 = (1.0 / 4.0) * Mtmp0;
  double Mtmp27 = Mtmp25 * Mtmp26;
  double Mtmp28 = Mtmp10 * q;
  double Mtmp29 = (1.0 / 6.0) * Mtmp3;
  double Mtmp30 = (y * y * y * y);
  double Mtmp31 = Mtmp22 * Mtmp30;
  double Mtmp32 = (1.0 / 4.0) * Mtmp10;
  double Mtmp33 = (z * z * z * z);
  double Mtmp34 = (x * x * x * x * x);
  double Mtmp35 = (1.0 / 120.0) * q;
  double Mtmp36 = Mtmp34 * Mtmp35;
  double Mtmp37 = (1.0 / 24.0) * Mtmp8;
  double Mtmp38 = (1.0 / 12.0) * Mtmp12;
  double Mtmp39 = Mtmp25 * Mtmp38;
  double Mtmp40 = (1.0 / 12.0) * Mtmp18;
  double Mtmp41 = Mtmp0 * q;
  double Mtmp42 = Mtmp40 * Mtmp41;
  double Mtmp43 = Mtmp10 * Mtmp8;
  double Mtmp44 = (1.0 / 12.0) * Mtmp20;
  double Mtmp45 = (1.0 / 24.0) * Mtmp3;
  double Mtmp46 = Mtmp3 * Mtmp6;
  double Mtmp47 = (y * y * y * y * y);
  double Mtmp48 = Mtmp35 * Mtmp47;
  double Mtmp49 = (z * z * z * z * z);
  double Mtmp50 = (x * x * x * x * x * x);
  double Mtmp51 = (1.0 / 720.0) * q;
  double Mtmp52 = Mtmp50 * Mtmp51;
  double Mtmp53 = (1.0 / 120.0) * Mtmp8;
  double Mtmp54 = (1.0 / 48.0) * Mtmp21;
  double Mtmp55 = Mtmp25 * Mtmp54;
  double Mtmp56 = Mtmp18 * q;
  double Mtmp57 = (1.0 / 36.0) * Mtmp12;
  double Mtmp58 = Mtmp56 * Mtmp57;
  double Mtmp59 = Mtmp20 * q;
  double Mtmp60 = (1.0 / 48.0) * Mtmp41;
  double Mtmp61 = Mtmp30 * Mtmp60;
  double Mtmp62 = Mtmp0 * Mtmp10;
  double Mtmp63 = Mtmp0 * Mtmp8;
  double Mtmp64 = (1.0 / 120.0) * Mtmp3;
  double Mtmp65 = Mtmp10 * Mtmp3;
  double Mtmp66 = (y * y * y * y * y * y);
  double Mtmp67 = Mtmp51 * Mtmp66;
  double Mtmp68 = (1.0 / 48.0) * Mtmp30;
  double Mtmp69 = (1.0 / 36.0) * Mtmp20;
  double Mtmp70 = (1.0 / 48.0) * Mtmp33;
  double Mtmp71 = (z * z * z * z * z * z);
  double Mtmp72 = (1.0 / 5040.0) * q;
  double Mtmp73 = (1.0 / 720.0) * Mtmp8;
  double Mtmp74 = (1.0 / 240.0) * Mtmp34;
  double Mtmp75 = (1.0 / 144.0) * Mtmp21;
  double Mtmp76 = (1.0 / 144.0) * Mtmp30;
  double Mtmp77 = Mtmp12 * q;
  double Mtmp78 = Mtmp22 * Mtmp6;
  double Mtmp79 = (1.0 / 144.0) * Mtmp33;
  double Mtmp80 = (1.0 / 240.0) * Mtmp41;
  double Mtmp81 = (1.0 / 720.0) * Mtmp3;
  M[0] += Mtmp2;
  M[1] += Mtmp4;
  M[2] += Mtmp5;
  M[3] += Mtmp7;
  M[4] += Mtmp9;
  M[5] += Mtmp11;
  M[6] += -Mtmp14;
  M[7] += -Mtmp15;
  M[8] += -Mtmp2 * z;
  M[9] += -Mtmp16;
  M[10] += -Mtmp4 * z;
  M[11] += -Mtmp17;
  M[12] += -Mtmp19;
  M[13] += -Mtmp7 * z;
  M[14] += -Mtmp11 * y;
  M[15] += -Mtmp13 * Mtmp20;
  M[16] += Mtmp23;
  M[17] += Mtmp12 * Mtmp24;
  M[18] += Mtmp14 * z;
  M[19] += Mtmp27;
  M[20] += Mtmp15 * z;
  M[21] += Mtmp26 * Mtmp28;
  M[22] += Mtmp18 * Mtmp29;
  M[23] += Mtmp16 * z;
  M[24] += Mtmp17 * y;
  M[25] += Mtmp20 * Mtmp29;
  M[26] += Mtmp31;
  M[27] += Mtmp19 * z;
  M[28] += Mtmp25 * Mtmp32;
  M[29] += Mtmp20 * Mtmp24;
  M[30] += Mtmp22 * Mtmp33;
  M[31] += -Mtmp36;
  M[32] += -Mtmp21 * Mtmp37;
  M[33] += -Mtmp23 * z;
  M[34] += -Mtmp39;
  M[35] += -1.0 / 6.0 * Mtmp12 * Mtmp9;
  M[36] += -Mtmp28 * Mtmp38;
  M[37] += -Mtmp42;
  M[38] += -Mtmp27 * z;
  M[39] += -Mtmp26 * Mtmp43;
  M[40] += -Mtmp41 * Mtmp44;
  M[41] += -Mtmp30 * Mtmp45;
  M[42] += -1.0 / 6.0 * Mtmp18 * Mtmp5;
  M[43] += -Mtmp32 * Mtmp46;
  M[44] += -1.0 / 6.0 * Mtmp20 * Mtmp4;
  M[45] += -Mtmp33 * Mtmp45;
  M[46] += -Mtmp48;
  M[47] += -Mtmp31 * z;
  M[48] += -Mtmp28 * Mtmp40;
  M[49] += -Mtmp25 * Mtmp44;
  M[50] += -Mtmp33 * Mtmp37;
  M[51] += -Mtmp35 * Mtmp49;
  M[52] += Mtmp52;
  M[53] += Mtmp34 * Mtmp53;
  M[54] += Mtmp36 * z;
  M[55] += Mtmp55;
  M[56] += (1.0 / 24.0) * Mtmp21 * Mtmp9;
  M[57] += Mtmp28 * Mtmp54;
  M[58] += Mtmp58;
  M[59] += Mtmp39 * z;
  M[60] += Mtmp38 * Mtmp43;
  M[61] += Mtmp57 * Mtmp59;
  M[62] += Mtmp61;
  M[63] += Mtmp42 * z;
  M[64] += (1.0 / 8.0) * Mtmp25 * Mtmp62;
  M[65] += Mtmp44 * Mtmp63;
  M[66] += Mtmp33 * Mtmp60;
  M[67] += Mtmp47 * Mtmp64;
  M[68] += (1.0 / 24.0) * Mtmp30 * Mtmp5;
  M[69] += Mtmp40 * Mtmp65;
  M[70] += Mtmp44 * Mtmp46;
  M[71] += (1.0 / 24.0) * Mtmp33 * Mtmp4;
  M[72] += Mtmp49 * Mtmp64;
  M[73] += Mtmp67;
  M[74] += Mtmp48 * z;
  M[75] += Mtmp28 * Mtmp68;
  M[76] += Mtmp56 * Mtmp69;
  M[77] += Mtmp25 * Mtmp70;
  M[78] += Mtmp49 * Mtmp53;
  M[79] += Mtmp51 * Mtmp71;
  M[80] += -Mtmp72 * (x * x * x * x * x * x * x);
  M[81] += -Mtmp50 * Mtmp73;
  M[82] += -Mtmp52 * z;
  M[83] += -Mtmp25 * Mtmp74;
  M[84] += -1.0 / 120.0 * Mtmp34 * Mtmp9;
  M[85] += -Mtmp28 * Mtmp74;
  M[86] += -Mtmp56 * Mtmp75;
  M[87] += -Mtmp55 * z;
  M[88] += -Mtmp43 * Mtmp54;
  M[89] += -Mtmp59 * Mtmp75;
  M[90] += -Mtmp76 * Mtmp77;
  M[91] += -Mtmp58 * z;
  M[92] += -Mtmp10 * Mtmp12 * Mtmp78;
  M[93] += -Mtmp20 * Mtmp57 * Mtmp8;
  M[94] += -Mtmp77 * Mtmp79;
  M[95] += -Mtmp47 * Mtmp80;
  M[96] += -Mtmp61 * z;
  M[97] += -Mtmp18 * Mtmp22 * Mtmp62;
  M[98] += -Mtmp0 * Mtmp20 * Mtmp78;
  M[99] += -Mtmp63 * Mtmp70;
  M[100] += -Mtmp49 * Mtmp80;
  M[101] += -Mtmp66 * Mtmp81;
  M[102] += -1.0 / 120.0 * Mtmp47 * Mtmp5;
  M[103] += -Mtmp65 * Mtmp68;
  M[104] += -Mtmp18 * Mtmp3 * Mtmp69;
  M[105] += -Mtmp46 * Mtmp70;
  M[106] += -1.0 / 120.0 * Mtmp4 * Mtmp49;
  M[107] += -Mtmp71 * Mtmp81;
  M[108] += -Mtmp72 * (y * y * y * y * y * y * y);
  M[109] += -Mtmp67 * z;
  M[110] += -1.0 / 240.0 * Mtmp28 * Mtmp47;
  M[111] += -Mtmp59 * Mtmp76;
  M[112] += -Mtmp56 * Mtmp79;
  M[113] += -1.0 / 240.0 * Mtmp25 * Mtmp49;
  M[114] += -Mtmp71 * Mtmp73;
  M[115] += -Mtmp72 * (z * z * z * z * z * z * z);
}
void field_m2_M2M_7(double x, double y, double z, double* M, double* Ms) {
  double Mstmp0   = x * M[0];
  double Mstmp1   = x * M[1];
  double Mstmp2   = y * M[0];
  double Mstmp3   = x * M[2];
  double Mstmp4   = z * M[0];
  double Mstmp5   = x * M[3];
  double Mstmp6   = y * M[1];
  double Mstmp7   = x * M[4];
  double Mstmp8   = y * M[2];
  double Mstmp9   = z * M[1];
  double Mstmp10  = x * M[5];
  double Mstmp11  = z * M[2];
  double Mstmp12  = y * M[3];
  double Mstmp13  = y * M[4];
  double Mstmp14  = z * M[3];
  double Mstmp15  = y * M[5];
  double Mstmp16  = z * M[4];
  double Mstmp17  = z * M[5];
  double Mstmp18  = x * M[6];
  double Mstmp19  = (x * x);
  double Mstmp20  = (1.0 / 2.0) * Mstmp19;
  double Mstmp21  = x * M[7];
  double Mstmp22  = y * M[6];
  double Mstmp23  = Mstmp0 * y;
  double Mstmp24  = x * M[8];
  double Mstmp25  = z * M[6];
  double Mstmp26  = Mstmp0 * z;
  double Mstmp27  = x * M[9];
  double Mstmp28  = y * M[7];
  double Mstmp29  = Mstmp1 * y;
  double Mstmp30  = (y * y);
  double Mstmp31  = (1.0 / 2.0) * M[0];
  double Mstmp32  = x * M[10];
  double Mstmp33  = y * M[8];
  double Mstmp34  = z * M[7];
  double Mstmp35  = Mstmp3 * y;
  double Mstmp36  = Mstmp1 * z;
  double Mstmp37  = Mstmp2 * z;
  double Mstmp38  = x * M[11];
  double Mstmp39  = z * M[8];
  double Mstmp40  = Mstmp3 * z;
  double Mstmp41  = (z * z);
  double Mstmp42  = x * M[12];
  double Mstmp43  = y * M[9];
  double Mstmp44  = Mstmp5 * y;
  double Mstmp45  = (1.0 / 2.0) * Mstmp30;
  double Mstmp46  = x * M[13];
  double Mstmp47  = y * M[10];
  double Mstmp48  = z * M[9];
  double Mstmp49  = Mstmp7 * y;
  double Mstmp50  = Mstmp5 * z;
  double Mstmp51  = Mstmp6 * z;
  double Mstmp52  = x * M[14];
  double Mstmp53  = y * M[11];
  double Mstmp54  = z * M[10];
  double Mstmp55  = Mstmp10 * y;
  double Mstmp56  = Mstmp7 * z;
  double Mstmp57  = Mstmp8 * z;
  double Mstmp58  = (1.0 / 2.0) * Mstmp41;
  double Mstmp59  = x * M[15];
  double Mstmp60  = z * M[11];
  double Mstmp61  = Mstmp10 * z;
  double Mstmp62  = y * M[12];
  double Mstmp63  = y * M[13];
  double Mstmp64  = z * M[12];
  double Mstmp65  = Mstmp12 * z;
  double Mstmp66  = y * M[14];
  double Mstmp67  = z * M[13];
  double Mstmp68  = Mstmp13 * z;
  double Mstmp69  = y * M[15];
  double Mstmp70  = z * M[14];
  double Mstmp71  = Mstmp15 * z;
  double Mstmp72  = z * M[15];
  double Mstmp73  = x * M[16];
  double Mstmp74  = (x * x * x);
  double Mstmp75  = (1.0 / 6.0) * Mstmp74;
  double Mstmp76  = x * M[17];
  double Mstmp77  = y * M[16];
  double Mstmp78  = Mstmp18 * y;
  double Mstmp79  = x * M[18];
  double Mstmp80  = z * M[16];
  double Mstmp81  = Mstmp18 * z;
  double Mstmp82  = x * M[19];
  double Mstmp83  = y * M[17];
  double Mstmp84  = Mstmp21 * y;
  double Mstmp85  = x * M[20];
  double Mstmp86  = y * M[18];
  double Mstmp87  = z * M[17];
  double Mstmp88  = Mstmp24 * y;
  double Mstmp89  = Mstmp21 * z;
  double Mstmp90  = Mstmp22 * z;
  double Mstmp91  = x * M[21];
  double Mstmp92  = z * M[18];
  double Mstmp93  = Mstmp24 * z;
  double Mstmp94  = x * M[22];
  double Mstmp95  = y * M[19];
  double Mstmp96  = Mstmp27 * y;
  double Mstmp97  = (y * y * y);
  double Mstmp98  = (1.0 / 6.0) * M[0];
  double Mstmp99  = x * M[23];
  double Mstmp100 = y * M[20];
  double Mstmp101 = z * M[19];
  double Mstmp102 = Mstmp32 * y;
  double Mstmp103 = Mstmp27 * z;
  double Mstmp104 = Mstmp28 * z;
  double Mstmp105 = x * M[24];
  double Mstmp106 = y * M[21];
  double Mstmp107 = z * M[20];
  double Mstmp108 = Mstmp38 * y;
  double Mstmp109 = Mstmp32 * z;
  double Mstmp110 = Mstmp33 * z;
  double Mstmp111 = x * M[25];
  double Mstmp112 = z * M[21];
  double Mstmp113 = Mstmp38 * z;
  double Mstmp114 = (z * z * z);
  double Mstmp115 = x * M[26];
  double Mstmp116 = y * M[22];
  double Mstmp117 = Mstmp42 * y;
  double Mstmp118 = (1.0 / 6.0) * Mstmp97;
  double Mstmp119 = x * M[27];
  double Mstmp120 = y * M[23];
  double Mstmp121 = z * M[22];
  double Mstmp122 = Mstmp46 * y;
  double Mstmp123 = Mstmp42 * z;
  double Mstmp124 = Mstmp43 * z;
  double Mstmp125 = x * M[28];
  double Mstmp126 = y * M[24];
  double Mstmp127 = z * M[23];
  double Mstmp128 = Mstmp52 * y;
  double Mstmp129 = Mstmp46 * z;
  double Mstmp130 = Mstmp47 * z;
  double Mstmp131 = x * M[29];
  double Mstmp132 = y * M[25];
  double Mstmp133 = z * M[24];
  double Mstmp134 = Mstmp59 * y;
  double Mstmp135 = Mstmp52 * z;
  double Mstmp136 = Mstmp53 * z;
  double Mstmp137 = (1.0 / 6.0) * Mstmp114;
  double Mstmp138 = x * M[30];
  double Mstmp139 = z * M[25];
  double Mstmp140 = Mstmp59 * z;
  double Mstmp141 = y * M[26];
  double Mstmp142 = y * M[27];
  double Mstmp143 = z * M[26];
  double Mstmp144 = Mstmp62 * z;
  double Mstmp145 = y * M[28];
  double Mstmp146 = z * M[27];
  double Mstmp147 = Mstmp63 * z;
  double Mstmp148 = y * M[29];
  double Mstmp149 = z * M[28];
  double Mstmp150 = Mstmp66 * z;
  double Mstmp151 = y * M[30];
  double Mstmp152 = z * M[29];
  double Mstmp153 = Mstmp69 * z;
  double Mstmp154 = z * M[30];
  double Mstmp155 = x * M[31];
  double Mstmp156 = (1.0 / 24.0) * (x * x * x * x);
  double Mstmp157 = x * M[32];
  double Mstmp158 = y * M[31];
  double Mstmp159 = Mstmp73 * y;
  double Mstmp160 = x * M[33];
  double Mstmp161 = x * M[34];
  double Mstmp162 = y * M[32];
  double Mstmp163 = Mstmp76 * y;
  double Mstmp164 = (1.0 / 4.0) * Mstmp19;
  double Mstmp165 = Mstmp30 * M[0];
  double Mstmp166 = x * M[35];
  double Mstmp167 = y * M[33];
  double Mstmp168 = Mstmp79 * y;
  double Mstmp169 = x * M[36];
  double Mstmp170 = Mstmp164 * Mstmp41;
  double Mstmp171 = x * M[37];
  double Mstmp172 = y * M[34];
  double Mstmp173 = Mstmp82 * y;
  double Mstmp174 = Mstmp164 * Mstmp30;
  double Mstmp175 = x * M[38];
  double Mstmp176 = y * M[35];
  double Mstmp177 = Mstmp85 * y;
  double Mstmp178 = x * M[39];
  double Mstmp179 = y * M[36];
  double Mstmp180 = Mstmp91 * y;
  double Mstmp181 = x * M[40];
  double Mstmp182 = x * M[41];
  double Mstmp183 = y * M[37];
  double Mstmp184 = Mstmp94 * y;
  double Mstmp185 = (y * y * y * y);
  double Mstmp186 = (1.0 / 24.0) * M[0];
  double Mstmp187 = x * M[42];
  double Mstmp188 = y * M[38];
  double Mstmp189 = Mstmp99 * y;
  double Mstmp190 = x * M[43];
  double Mstmp191 = y * M[39];
  double Mstmp192 = Mstmp105 * y;
  double Mstmp193 = (1.0 / 4.0) * Mstmp41;
  double Mstmp194 = x * M[44];
  double Mstmp195 = y * M[40];
  double Mstmp196 = Mstmp111 * y;
  double Mstmp197 = x * M[45];
  double Mstmp198 = (z * z * z * z);
  double Mstmp199 = x * M[46];
  double Mstmp200 = y * M[41];
  double Mstmp201 = Mstmp115 * y;
  double Mstmp202 = (1.0 / 24.0) * Mstmp185;
  double Mstmp203 = x * M[47];
  double Mstmp204 = y * M[42];
  double Mstmp205 = Mstmp119 * y;
  double Mstmp206 = x * M[48];
  double Mstmp207 = y * M[43];
  double Mstmp208 = Mstmp125 * y;
  double Mstmp209 = Mstmp193 * Mstmp30;
  double Mstmp210 = x * M[49];
  double Mstmp211 = y * M[44];
  double Mstmp212 = Mstmp131 * y;
  double Mstmp213 = x * M[50];
  double Mstmp214 = y * M[45];
  double Mstmp215 = Mstmp138 * y;
  double Mstmp216 = (1.0 / 24.0) * Mstmp198;
  double Mstmp217 = x * M[51];
  double Mstmp218 = y * M[46];
  double Mstmp219 = y * M[47];
  double Mstmp220 = y * M[48];
  double Mstmp221 = y * M[49];
  double Mstmp222 = y * M[50];
  double Mstmp223 = y * M[51];
  double Mstmp224 = (1.0 / 120.0) * (x * x * x * x * x);
  double Mstmp225 = (1.0 / 12.0) * Mstmp74;
  double Mstmp226 = Mstmp225 * Mstmp41;
  double Mstmp227 = (1.0 / 12.0) * Mstmp19;
  double Mstmp228 = Mstmp227 * M[0];
  double Mstmp229 = Mstmp225 * Mstmp30;
  double Mstmp230 = Mstmp227 * Mstmp97;
  double Mstmp231 = Mstmp114 * Mstmp227;
  double Mstmp232 = (y * y * y * y * y);
  double Mstmp233 = (1.0 / 120.0) * M[0];
  double Mstmp234 = (1.0 / 12.0) * Mstmp41 * Mstmp97;
  double Mstmp235 = (1.0 / 12.0) * Mstmp114;
  double Mstmp236 = (z * z * z * z * z);
  double Mstmp237 = (1.0 / 120.0) * Mstmp232;
  double Mstmp238 = Mstmp235 * Mstmp30;
  double Mstmp239 = (1.0 / 120.0) * Mstmp236;
#pragma omp atomic
  Ms[0] += M[0];
#pragma omp atomic
  Ms[1] += M[1];
#pragma omp atomic
  Ms[2] += M[2];
#pragma omp atomic
  Ms[3] += M[3];
#pragma omp atomic
  Ms[4] += M[4];
#pragma omp atomic
  Ms[5] += M[5];
#pragma omp atomic
  Ms[6] += Mstmp0 + M[6];
#pragma omp atomic
  Ms[7] += Mstmp1 + Mstmp2 + M[7];
#pragma omp atomic
  Ms[8] += Mstmp3 + Mstmp4 + M[8];
#pragma omp atomic
  Ms[9] += Mstmp5 + Mstmp6 + M[9];
#pragma omp atomic
  Ms[10] += Mstmp7 + Mstmp8 + Mstmp9 + M[10];
#pragma omp atomic
  Ms[11] += Mstmp10 + Mstmp11 + M[11];
#pragma omp atomic
  Ms[12] += Mstmp12 + M[12];
#pragma omp atomic
  Ms[13] += Mstmp13 + Mstmp14 + M[13];
#pragma omp atomic
  Ms[14] += Mstmp15 + Mstmp16 + M[14];
#pragma omp atomic
  Ms[15] += Mstmp17 + M[15];
#pragma omp atomic
  Ms[16] += Mstmp18 + Mstmp20 * M[0] + M[16];
#pragma omp atomic
  Ms[17] += Mstmp20 * M[1] + Mstmp21 + Mstmp22 + Mstmp23 + M[17];
#pragma omp atomic
  Ms[18] += Mstmp20 * M[2] + Mstmp24 + Mstmp25 + Mstmp26 + M[18];
#pragma omp atomic
  Ms[19] += Mstmp20 * M[3] + Mstmp27 + Mstmp28 + Mstmp29 + Mstmp30 * Mstmp31 + M[19];
#pragma omp atomic
  Ms[20] += Mstmp20 * M[4] + Mstmp32 + Mstmp33 + Mstmp34 + Mstmp35 + Mstmp36 + Mstmp37 +
            M[20];
#pragma omp atomic
  Ms[21] += Mstmp20 * M[5] + Mstmp31 * Mstmp41 + Mstmp38 + Mstmp39 + Mstmp40 + M[21];
#pragma omp atomic
  Ms[22] += Mstmp42 + Mstmp43 + Mstmp44 + Mstmp45 * M[1] + M[22];
#pragma omp atomic
  Ms[23] += Mstmp45 * M[2] + Mstmp46 + Mstmp47 + Mstmp48 + Mstmp49 + Mstmp50 + Mstmp51 +
            M[23];
#pragma omp atomic
  Ms[24] += Mstmp52 + Mstmp53 + Mstmp54 + Mstmp55 + Mstmp56 + Mstmp57 + Mstmp58 * M[1] +
            M[24];
#pragma omp atomic
  Ms[25] += Mstmp58 * M[2] + Mstmp59 + Mstmp60 + Mstmp61 + M[25];
#pragma omp atomic
  Ms[26] += Mstmp45 * M[3] + Mstmp62 + M[26];
#pragma omp atomic
  Ms[27] += Mstmp45 * M[4] + Mstmp63 + Mstmp64 + Mstmp65 + M[27];
#pragma omp atomic
  Ms[28] += Mstmp45 * M[5] + Mstmp58 * M[3] + Mstmp66 + Mstmp67 + Mstmp68 + M[28];
#pragma omp atomic
  Ms[29] += Mstmp58 * M[4] + Mstmp69 + Mstmp70 + Mstmp71 + M[29];
#pragma omp atomic
  Ms[30] += Mstmp58 * M[5] + Mstmp72 + M[30];
#pragma omp atomic
  Ms[31] += Mstmp20 * M[6] + Mstmp73 + Mstmp75 * M[0] + M[31];
#pragma omp atomic
  Ms[32] += Mstmp2 * Mstmp20 + Mstmp20 * M[7] + Mstmp75 * M[1] + Mstmp76 + Mstmp77 +
            Mstmp78 + M[32];
#pragma omp atomic
  Ms[33] += Mstmp20 * Mstmp4 + Mstmp20 * M[8] + Mstmp75 * M[2] + Mstmp79 + Mstmp80 +
            Mstmp81 + M[33];
#pragma omp atomic
  Ms[34] += Mstmp0 * Mstmp45 + Mstmp20 * Mstmp6 + Mstmp20 * M[9] + Mstmp45 * M[6] +
            Mstmp75 * M[3] + Mstmp82 + Mstmp83 + Mstmp84 + M[34];
#pragma omp atomic
  Ms[35] += Mstmp20 * Mstmp8 + Mstmp20 * Mstmp9 + Mstmp20 * M[10] + Mstmp23 * z +
            Mstmp75 * M[4] + Mstmp85 + Mstmp86 + Mstmp87 + Mstmp88 + Mstmp89 + Mstmp90 +
            M[35];
#pragma omp atomic
  Ms[36] += Mstmp0 * Mstmp58 + Mstmp11 * Mstmp20 + Mstmp20 * M[11] + Mstmp58 * M[6] +
            Mstmp75 * M[5] + Mstmp91 + Mstmp92 + Mstmp93 + M[36];
#pragma omp atomic
  Ms[37] += Mstmp1 * Mstmp45 + Mstmp12 * Mstmp20 + Mstmp20 * M[12] + Mstmp45 * M[7] +
            Mstmp94 + Mstmp95 + Mstmp96 + Mstmp97 * Mstmp98 + M[37];
#pragma omp atomic
  Ms[38] += Mstmp100 + Mstmp101 + Mstmp102 + Mstmp103 + Mstmp104 + Mstmp13 * Mstmp20 +
            Mstmp14 * Mstmp20 + Mstmp20 * M[13] + Mstmp29 * z + Mstmp3 * Mstmp45 +
            Mstmp4 * Mstmp45 + Mstmp45 * M[8] + Mstmp99 + M[38];
#pragma omp atomic
  Ms[39] += Mstmp1 * Mstmp58 + Mstmp105 + Mstmp106 + Mstmp107 + Mstmp108 + Mstmp109 +
            Mstmp110 + Mstmp15 * Mstmp20 + Mstmp16 * Mstmp20 + Mstmp2 * Mstmp58 +
            Mstmp20 * M[14] + Mstmp35 * z + Mstmp58 * M[7] + M[39];
#pragma omp atomic
  Ms[40] += Mstmp111 + Mstmp112 + Mstmp113 + Mstmp114 * Mstmp98 + Mstmp17 * Mstmp20 +
            Mstmp20 * M[15] + Mstmp3 * Mstmp58 + Mstmp58 * M[8] + M[40];
#pragma omp atomic
  Ms[41] += Mstmp115 + Mstmp116 + Mstmp117 + Mstmp118 * M[1] + Mstmp45 * Mstmp5 +
            Mstmp45 * M[9] + M[41];
#pragma omp atomic
  Ms[42] += Mstmp118 * M[2] + Mstmp119 + Mstmp120 + Mstmp121 + Mstmp122 + Mstmp123 +
            Mstmp124 + Mstmp44 * z + Mstmp45 * Mstmp7 + Mstmp45 * Mstmp9 +
            Mstmp45 * M[10] + M[42];
#pragma omp atomic
  Ms[43] += Mstmp10 * Mstmp45 + Mstmp11 * Mstmp45 + Mstmp125 + Mstmp126 + Mstmp127 +
            Mstmp128 + Mstmp129 + Mstmp130 + Mstmp45 * M[11] + Mstmp49 * z +
            Mstmp5 * Mstmp58 + Mstmp58 * Mstmp6 + Mstmp58 * M[9] + M[43];
#pragma omp atomic
  Ms[44] += Mstmp131 + Mstmp132 + Mstmp133 + Mstmp134 + Mstmp135 + Mstmp136 +
            Mstmp137 * M[1] + Mstmp55 * z + Mstmp58 * Mstmp7 + Mstmp58 * Mstmp8 +
            Mstmp58 * M[10] + M[44];
#pragma omp atomic
  Ms[45] += Mstmp10 * Mstmp58 + Mstmp137 * M[2] + Mstmp138 + Mstmp139 + Mstmp140 +
            Mstmp58 * M[11] + M[45];
#pragma omp atomic
  Ms[46] += Mstmp118 * M[3] + Mstmp141 + Mstmp45 * M[12] + M[46];
#pragma omp atomic
  Ms[47] += Mstmp118 * M[4] + Mstmp14 * Mstmp45 + Mstmp142 + Mstmp143 + Mstmp144 +
            Mstmp45 * M[13] + M[47];
#pragma omp atomic
  Ms[48] += Mstmp118 * M[5] + Mstmp12 * Mstmp58 + Mstmp145 + Mstmp146 + Mstmp147 +
            Mstmp16 * Mstmp45 + Mstmp45 * M[14] + Mstmp58 * M[12] + M[48];
#pragma omp atomic
  Ms[49] += Mstmp13 * Mstmp58 + Mstmp137 * M[3] + Mstmp148 + Mstmp149 + Mstmp150 +
            Mstmp17 * Mstmp45 + Mstmp45 * M[15] + Mstmp58 * M[13] + M[49];
#pragma omp atomic
  Ms[50] += Mstmp137 * M[4] + Mstmp15 * Mstmp58 + Mstmp151 + Mstmp152 + Mstmp153 +
            Mstmp58 * M[14] + M[50];
#pragma omp atomic
  Ms[51] += Mstmp137 * M[5] + Mstmp154 + Mstmp58 * M[15] + M[51];
#pragma omp atomic
  Ms[52] += Mstmp155 + Mstmp156 * M[0] + Mstmp20 * M[16] + Mstmp75 * M[6] + M[52];
#pragma omp atomic
  Ms[53] += Mstmp156 * M[1] + Mstmp157 + Mstmp158 + Mstmp159 + Mstmp2 * Mstmp75 +
            Mstmp20 * Mstmp22 + Mstmp20 * M[17] + Mstmp75 * M[7] + M[53];
#pragma omp atomic
  Ms[54] += Mstmp156 * M[2] + Mstmp160 + Mstmp20 * Mstmp25 + Mstmp20 * M[18] +
            Mstmp4 * Mstmp75 + Mstmp73 * z + Mstmp75 * M[8] + z * M[31] + M[54];
#pragma omp atomic
  Ms[55] += Mstmp156 * M[3] + Mstmp161 + Mstmp162 + Mstmp163 + Mstmp164 * Mstmp165 +
            Mstmp18 * Mstmp45 + Mstmp20 * Mstmp28 + Mstmp20 * M[19] + Mstmp45 * M[16] +
            Mstmp6 * Mstmp75 + Mstmp75 * M[9] + M[55];
#pragma omp atomic
  Ms[56] += Mstmp156 * M[4] + Mstmp166 + Mstmp167 + Mstmp168 + Mstmp20 * Mstmp33 +
            Mstmp20 * Mstmp34 + Mstmp20 * Mstmp37 + Mstmp20 * M[20] + Mstmp75 * Mstmp8 +
            Mstmp75 * Mstmp9 + Mstmp75 * M[10] + Mstmp76 * z + Mstmp77 * z + Mstmp78 * z +
            z * M[32] + M[56];
#pragma omp atomic
  Ms[57] += Mstmp11 * Mstmp75 + Mstmp156 * M[5] + Mstmp169 + Mstmp170 * M[0] +
            Mstmp18 * Mstmp58 + Mstmp20 * Mstmp39 + Mstmp20 * M[21] + Mstmp58 * M[16] +
            Mstmp75 * M[11] + Mstmp79 * z + z * M[33] + M[57];
#pragma omp atomic
  Ms[58] += Mstmp0 * Mstmp118 + Mstmp118 * M[6] + Mstmp12 * Mstmp75 + Mstmp171 +
            Mstmp172 + Mstmp173 + Mstmp174 * M[1] + Mstmp20 * Mstmp43 + Mstmp20 * M[22] +
            Mstmp21 * Mstmp45 + Mstmp45 * M[17] + Mstmp75 * M[12] + M[58];
#pragma omp atomic
  Ms[59] += Mstmp13 * Mstmp75 + Mstmp14 * Mstmp75 + Mstmp174 * M[2] + Mstmp175 +
            Mstmp176 + Mstmp177 + Mstmp20 * Mstmp47 + Mstmp20 * Mstmp48 +
            Mstmp20 * Mstmp51 + Mstmp20 * M[23] + Mstmp24 * Mstmp45 + Mstmp25 * Mstmp45 +
            Mstmp26 * Mstmp45 + Mstmp45 * M[18] + Mstmp75 * M[13] + Mstmp82 * z +
            Mstmp83 * z + Mstmp84 * z + z * M[34] + M[59];
#pragma omp atomic
  Ms[60] += Mstmp15 * Mstmp75 + Mstmp16 * Mstmp75 + Mstmp170 * M[1] + Mstmp178 +
            Mstmp179 + Mstmp180 + Mstmp20 * Mstmp53 + Mstmp20 * Mstmp54 +
            Mstmp20 * Mstmp57 + Mstmp20 * M[24] + Mstmp21 * Mstmp58 + Mstmp22 * Mstmp58 +
            Mstmp23 * Mstmp58 + Mstmp58 * M[17] + Mstmp75 * M[14] + Mstmp85 * z +
            Mstmp86 * z + Mstmp88 * z + z * M[35] + M[60];
#pragma omp atomic
  Ms[61] += Mstmp0 * Mstmp137 + Mstmp137 * M[6] + Mstmp17 * Mstmp75 + Mstmp170 * M[2] +
            Mstmp181 + Mstmp20 * Mstmp60 + Mstmp20 * M[25] + Mstmp24 * Mstmp58 +
            Mstmp58 * M[18] + Mstmp75 * M[15] + Mstmp91 * z + z * M[36] + M[61];
#pragma omp atomic
  Ms[62] += Mstmp1 * Mstmp118 + Mstmp118 * M[7] + Mstmp174 * M[3] + Mstmp182 + Mstmp183 +
            Mstmp184 + Mstmp185 * Mstmp186 + Mstmp20 * Mstmp62 + Mstmp20 * M[26] +
            Mstmp27 * Mstmp45 + Mstmp45 * M[19] + M[62];
#pragma omp atomic
  Ms[63] += Mstmp118 * Mstmp3 + Mstmp118 * Mstmp4 + Mstmp118 * M[8] + Mstmp174 * M[4] +
            Mstmp187 + Mstmp188 + Mstmp189 + Mstmp20 * Mstmp63 + Mstmp20 * Mstmp64 +
            Mstmp20 * Mstmp65 + Mstmp20 * M[27] + Mstmp32 * Mstmp45 + Mstmp34 * Mstmp45 +
            Mstmp36 * Mstmp45 + Mstmp45 * M[20] + Mstmp94 * z + Mstmp95 * z +
            Mstmp96 * z + z * M[37] + M[63];
#pragma omp atomic
  Ms[64] += Mstmp100 * z + Mstmp102 * z + Mstmp165 * Mstmp193 + Mstmp170 * M[3] +
            Mstmp174 * M[5] + Mstmp190 + Mstmp191 + Mstmp192 + Mstmp20 * Mstmp66 +
            Mstmp20 * Mstmp67 + Mstmp20 * Mstmp68 + Mstmp20 * M[28] + Mstmp27 * Mstmp58 +
            Mstmp28 * Mstmp58 + Mstmp29 * Mstmp58 + Mstmp38 * Mstmp45 +
            Mstmp39 * Mstmp45 + Mstmp40 * Mstmp45 + Mstmp45 * M[21] + Mstmp58 * M[19] +
            Mstmp99 * z + z * M[38] + M[64];
#pragma omp atomic
  Ms[65] += Mstmp1 * Mstmp137 + Mstmp105 * z + Mstmp106 * z + Mstmp108 * z +
            Mstmp137 * Mstmp2 + Mstmp137 * M[7] + Mstmp170 * M[4] + Mstmp194 + Mstmp195 +
            Mstmp196 + Mstmp20 * Mstmp69 + Mstmp20 * Mstmp70 + Mstmp20 * Mstmp71 +
            Mstmp20 * M[29] + Mstmp32 * Mstmp58 + Mstmp33 * Mstmp58 + Mstmp35 * Mstmp58 +
            Mstmp58 * M[20] + z * M[39] + M[65];
#pragma omp atomic
  Ms[66] += Mstmp111 * z + Mstmp137 * Mstmp3 + Mstmp137 * M[8] + Mstmp170 * M[5] +
            Mstmp186 * Mstmp198 + Mstmp197 + Mstmp20 * Mstmp72 + Mstmp20 * M[30] +
            Mstmp38 * Mstmp58 + Mstmp58 * M[21] + z * M[40] + M[66];
#pragma omp atomic
  Ms[67] += Mstmp118 * Mstmp5 + Mstmp118 * M[9] + Mstmp199 + Mstmp200 + Mstmp201 +
            Mstmp202 * M[1] + Mstmp42 * Mstmp45 + Mstmp45 * M[22] + M[67];
#pragma omp atomic
  Ms[68] += Mstmp115 * z + Mstmp116 * z + Mstmp117 * z + Mstmp118 * Mstmp7 +
            Mstmp118 * Mstmp9 + Mstmp118 * M[10] + Mstmp202 * M[2] + Mstmp203 + Mstmp204 +
            Mstmp205 + Mstmp45 * Mstmp46 + Mstmp45 * Mstmp48 + Mstmp45 * Mstmp50 +
            Mstmp45 * M[23] + z * M[41] + M[68];
#pragma omp atomic
  Ms[69] += Mstmp10 * Mstmp118 + Mstmp11 * Mstmp118 + Mstmp118 * M[11] + Mstmp119 * z +
            Mstmp120 * z + Mstmp122 * z + Mstmp206 + Mstmp207 + Mstmp208 +
            Mstmp209 * M[1] + Mstmp42 * Mstmp58 + Mstmp43 * Mstmp58 + Mstmp44 * Mstmp58 +
            Mstmp45 * Mstmp52 + Mstmp45 * Mstmp54 + Mstmp45 * Mstmp56 + Mstmp45 * M[24] +
            Mstmp58 * M[22] + z * M[42] + M[69];
#pragma omp atomic
  Ms[70] += Mstmp125 * z + Mstmp126 * z + Mstmp128 * z + Mstmp137 * Mstmp5 +
            Mstmp137 * Mstmp6 + Mstmp137 * M[9] + Mstmp209 * M[2] + Mstmp210 + Mstmp211 +
            Mstmp212 + Mstmp45 * Mstmp59 + Mstmp45 * Mstmp60 + Mstmp45 * Mstmp61 +
            Mstmp45 * M[25] + Mstmp46 * Mstmp58 + Mstmp47 * Mstmp58 + Mstmp49 * Mstmp58 +
            Mstmp58 * M[23] + z * M[43] + M[70];
#pragma omp atomic
  Ms[71] += Mstmp131 * z + Mstmp132 * z + Mstmp134 * z + Mstmp137 * Mstmp7 +
            Mstmp137 * Mstmp8 + Mstmp137 * M[10] + Mstmp213 + Mstmp214 + Mstmp215 +
            Mstmp216 * M[1] + Mstmp52 * Mstmp58 + Mstmp53 * Mstmp58 + Mstmp55 * Mstmp58 +
            Mstmp58 * M[24] + z * M[44] + M[71];
#pragma omp atomic
  Ms[72] += Mstmp10 * Mstmp137 + Mstmp137 * M[11] + Mstmp138 * z + Mstmp216 * M[2] +
            Mstmp217 + Mstmp58 * Mstmp59 + Mstmp58 * M[25] + z * M[45] + M[72];
#pragma omp atomic
  Ms[73] += Mstmp118 * M[12] + Mstmp202 * M[3] + Mstmp218 + Mstmp45 * M[26] + M[73];
#pragma omp atomic
  Ms[74] += Mstmp118 * Mstmp14 + Mstmp118 * M[13] + Mstmp141 * z + Mstmp202 * M[4] +
            Mstmp219 + Mstmp45 * Mstmp64 + Mstmp45 * M[27] + z * M[46] + M[74];
#pragma omp atomic
  Ms[75] += Mstmp118 * Mstmp16 + Mstmp118 * M[14] + Mstmp142 * z + Mstmp202 * M[5] +
            Mstmp209 * M[3] + Mstmp220 + Mstmp45 * Mstmp67 + Mstmp45 * M[28] +
            Mstmp58 * Mstmp62 + Mstmp58 * M[26] + z * M[47] + M[75];
#pragma omp atomic
  Ms[76] += Mstmp118 * Mstmp17 + Mstmp118 * M[15] + Mstmp12 * Mstmp137 +
            Mstmp137 * M[12] + Mstmp145 * z + Mstmp209 * M[4] + Mstmp221 +
            Mstmp45 * Mstmp70 + Mstmp45 * M[29] + Mstmp58 * Mstmp63 + Mstmp58 * M[27] +
            z * M[48] + M[76];
#pragma omp atomic
  Ms[77] += Mstmp13 * Mstmp137 + Mstmp137 * M[13] + Mstmp148 * z + Mstmp209 * M[5] +
            Mstmp216 * M[3] + Mstmp222 + Mstmp45 * Mstmp72 + Mstmp45 * M[30] +
            Mstmp58 * Mstmp66 + Mstmp58 * M[28] + z * M[49] + M[77];
#pragma omp atomic
  Ms[78] += Mstmp137 * Mstmp15 + Mstmp137 * M[14] + Mstmp151 * z + Mstmp216 * M[4] +
            Mstmp223 + Mstmp58 * Mstmp69 + Mstmp58 * M[29] + z * M[50] + M[78];
#pragma omp atomic
  Ms[79] += Mstmp137 * M[15] + Mstmp216 * M[5] + Mstmp58 * M[30] + z * M[51] + M[79];
#pragma omp atomic
  Ms[80] += Mstmp156 * M[6] + Mstmp20 * M[31] + Mstmp224 * M[0] + Mstmp75 * M[16] +
            x * M[52] + M[80];
#pragma omp atomic
  Ms[81] += Mstmp155 * y + Mstmp156 * Mstmp2 + Mstmp156 * M[7] + Mstmp20 * Mstmp77 +
            Mstmp20 * M[32] + Mstmp22 * Mstmp75 + Mstmp224 * M[1] + Mstmp75 * M[17] +
            x * M[53] + y * M[52] + M[81];
#pragma omp atomic
  Ms[82] += Mstmp155 * z + Mstmp156 * Mstmp4 + Mstmp156 * M[8] + Mstmp20 * Mstmp80 +
            Mstmp20 * M[33] + Mstmp224 * M[2] + Mstmp25 * Mstmp75 + Mstmp75 * M[18] +
            x * M[54] + z * M[52] + M[82];
#pragma omp atomic
  Ms[83] += Mstmp156 * Mstmp6 + Mstmp156 * M[9] + Mstmp157 * y + Mstmp165 * Mstmp225 +
            Mstmp174 * M[6] + Mstmp20 * Mstmp83 + Mstmp20 * M[34] + Mstmp224 * M[3] +
            Mstmp28 * Mstmp75 + Mstmp45 * Mstmp73 + Mstmp45 * M[31] + Mstmp75 * M[19] +
            x * M[55] + y * M[53] + M[83];
#pragma omp atomic
  Ms[84] += Mstmp156 * Mstmp8 + Mstmp156 * Mstmp9 + Mstmp156 * M[10] + Mstmp157 * z +
            Mstmp158 * z + Mstmp159 * z + Mstmp160 * y + Mstmp20 * Mstmp86 +
            Mstmp20 * Mstmp87 + Mstmp20 * Mstmp90 + Mstmp20 * M[35] + Mstmp224 * M[4] +
            Mstmp33 * Mstmp75 + Mstmp34 * Mstmp75 + Mstmp37 * Mstmp75 + Mstmp75 * M[20] +
            x * M[56] + y * M[54] + z * M[53] + M[84];
#pragma omp atomic
  Ms[85] += Mstmp11 * Mstmp156 + Mstmp156 * M[11] + Mstmp160 * z + Mstmp170 * M[6] +
            Mstmp20 * Mstmp92 + Mstmp20 * M[36] + Mstmp224 * M[5] + Mstmp226 * M[0] +
            Mstmp39 * Mstmp75 + Mstmp58 * Mstmp73 + Mstmp58 * M[31] + Mstmp75 * M[21] +
            x * M[57] + z * M[54] + M[85];
#pragma omp atomic
  Ms[86] += Mstmp118 * Mstmp18 + Mstmp118 * M[16] + Mstmp12 * Mstmp156 +
            Mstmp156 * M[12] + Mstmp161 * y + Mstmp174 * M[7] + Mstmp20 * Mstmp95 +
            Mstmp20 * M[37] + Mstmp228 * Mstmp97 + Mstmp229 * M[1] + Mstmp43 * Mstmp75 +
            Mstmp45 * Mstmp76 + Mstmp45 * M[32] + Mstmp75 * M[22] + x * M[58] +
            y * M[55] + M[86];
#pragma omp atomic
  Ms[87] += Mstmp100 * Mstmp20 + Mstmp101 * Mstmp20 + Mstmp104 * Mstmp20 +
            Mstmp13 * Mstmp156 + Mstmp14 * Mstmp156 + Mstmp156 * M[13] + Mstmp161 * z +
            Mstmp162 * z + Mstmp163 * z + Mstmp166 * y + Mstmp174 * Mstmp4 +
            Mstmp174 * M[8] + Mstmp20 * M[38] + Mstmp229 * M[2] + Mstmp45 * Mstmp79 +
            Mstmp45 * Mstmp80 + Mstmp45 * Mstmp81 + Mstmp45 * M[33] + Mstmp47 * Mstmp75 +
            Mstmp48 * Mstmp75 + Mstmp51 * Mstmp75 + Mstmp75 * M[23] + x * M[59] +
            y * M[56] + z * M[55] + M[87];
#pragma omp atomic
  Ms[88] += Mstmp106 * Mstmp20 + Mstmp107 * Mstmp20 + Mstmp110 * Mstmp20 +
            Mstmp15 * Mstmp156 + Mstmp156 * Mstmp16 + Mstmp156 * M[14] + Mstmp166 * z +
            Mstmp167 * z + Mstmp168 * z + Mstmp169 * y + Mstmp170 * Mstmp2 +
            Mstmp170 * M[7] + Mstmp20 * M[39] + Mstmp226 * M[1] + Mstmp53 * Mstmp75 +
            Mstmp54 * Mstmp75 + Mstmp57 * Mstmp75 + Mstmp58 * Mstmp76 +
            Mstmp58 * Mstmp77 + Mstmp58 * Mstmp78 + Mstmp58 * M[32] + Mstmp75 * M[24] +
            x * M[60] + y * M[57] + z * M[56] + M[88];
#pragma omp atomic
  Ms[89] += Mstmp112 * Mstmp20 + Mstmp114 * Mstmp228 + Mstmp137 * Mstmp18 +
            Mstmp137 * M[16] + Mstmp156 * Mstmp17 + Mstmp156 * M[15] + Mstmp169 * z +
            Mstmp170 * M[8] + Mstmp20 * M[40] + Mstmp226 * M[2] + Mstmp58 * Mstmp79 +
            Mstmp58 * M[33] + Mstmp60 * Mstmp75 + Mstmp75 * M[25] + x * M[61] +
            z * M[57] + M[89];
#pragma omp atomic
  Ms[90] += Mstmp0 * Mstmp202 + Mstmp116 * Mstmp20 + Mstmp118 * Mstmp21 +
            Mstmp118 * M[17] + Mstmp171 * y + Mstmp174 * M[9] + Mstmp20 * M[41] +
            Mstmp202 * M[6] + Mstmp229 * M[3] + Mstmp230 * M[1] + Mstmp45 * Mstmp82 +
            Mstmp45 * M[34] + Mstmp62 * Mstmp75 + Mstmp75 * M[26] + x * M[62] +
            y * M[58] + M[90];
#pragma omp atomic
  Ms[91] += Mstmp118 * Mstmp24 + Mstmp118 * Mstmp25 + Mstmp118 * Mstmp26 +
            Mstmp118 * M[18] + Mstmp120 * Mstmp20 + Mstmp121 * Mstmp20 +
            Mstmp124 * Mstmp20 + Mstmp171 * z + Mstmp172 * z + Mstmp173 * z +
            Mstmp174 * Mstmp9 + Mstmp174 * M[10] + Mstmp175 * y + Mstmp20 * M[42] +
            Mstmp229 * M[4] + Mstmp230 * M[2] + Mstmp45 * Mstmp85 + Mstmp45 * Mstmp87 +
            Mstmp45 * Mstmp89 + Mstmp45 * M[35] + Mstmp63 * Mstmp75 + Mstmp64 * Mstmp75 +
            Mstmp65 * Mstmp75 + Mstmp75 * M[27] + x * M[63] + y * M[59] + z * M[58] +
            M[91];
#pragma omp atomic
  Ms[92] += Mstmp0 * Mstmp209 + Mstmp11 * Mstmp174 + Mstmp126 * Mstmp20 +
            Mstmp127 * Mstmp20 + Mstmp130 * Mstmp20 + Mstmp170 * Mstmp6 +
            Mstmp170 * M[9] + Mstmp174 * M[11] + Mstmp175 * z + Mstmp176 * z +
            Mstmp177 * z + Mstmp178 * y + Mstmp20 * M[43] + Mstmp209 * M[6] +
            Mstmp226 * M[3] + Mstmp229 * M[5] + Mstmp45 * Mstmp91 + Mstmp45 * Mstmp92 +
            Mstmp45 * Mstmp93 + Mstmp45 * M[36] + Mstmp58 * Mstmp82 + Mstmp58 * Mstmp83 +
            Mstmp58 * Mstmp84 + Mstmp58 * M[34] + Mstmp66 * Mstmp75 + Mstmp67 * Mstmp75 +
            Mstmp68 * Mstmp75 + Mstmp75 * M[28] + x * M[64] + y * M[60] + z * M[59] +
            M[92];
#pragma omp atomic
  Ms[93] += Mstmp132 * Mstmp20 + Mstmp133 * Mstmp20 + Mstmp136 * Mstmp20 +
            Mstmp137 * Mstmp21 + Mstmp137 * Mstmp22 + Mstmp137 * Mstmp23 +
            Mstmp137 * M[17] + Mstmp170 * Mstmp8 + Mstmp170 * M[10] + Mstmp178 * z +
            Mstmp179 * z + Mstmp180 * z + Mstmp181 * y + Mstmp20 * M[44] +
            Mstmp226 * M[4] + Mstmp231 * M[1] + Mstmp58 * Mstmp85 + Mstmp58 * Mstmp86 +
            Mstmp58 * Mstmp88 + Mstmp58 * M[35] + Mstmp69 * Mstmp75 + Mstmp70 * Mstmp75 +
            Mstmp71 * Mstmp75 + Mstmp75 * M[29] + x * M[65] + y * M[61] + z * M[60] +
            M[93];
#pragma omp atomic
  Ms[94] += Mstmp0 * Mstmp216 + Mstmp137 * Mstmp24 + Mstmp137 * M[18] +
            Mstmp139 * Mstmp20 + Mstmp170 * M[11] + Mstmp181 * z + Mstmp20 * M[45] +
            Mstmp216 * M[6] + Mstmp226 * M[5] + Mstmp231 * M[2] + Mstmp58 * Mstmp91 +
            Mstmp58 * M[36] + Mstmp72 * Mstmp75 + Mstmp75 * M[30] + x * M[66] +
            z * M[61] + M[94];
#pragma omp atomic
  Ms[95] += Mstmp1 * Mstmp202 + Mstmp118 * Mstmp27 + Mstmp118 * M[19] +
            Mstmp141 * Mstmp20 + Mstmp174 * M[12] + Mstmp182 * y + Mstmp20 * M[46] +
            Mstmp202 * M[7] + Mstmp230 * M[3] + Mstmp232 * Mstmp233 + Mstmp45 * Mstmp94 +
            Mstmp45 * M[37] + x * M[67] + y * M[62] + M[95];
#pragma omp atomic
  Ms[96] += Mstmp101 * Mstmp45 + Mstmp103 * Mstmp45 + Mstmp118 * Mstmp32 +
            Mstmp118 * Mstmp34 + Mstmp118 * Mstmp36 + Mstmp118 * M[20] +
            Mstmp14 * Mstmp174 + Mstmp142 * Mstmp20 + Mstmp143 * Mstmp20 +
            Mstmp144 * Mstmp20 + Mstmp174 * M[13] + Mstmp182 * z + Mstmp183 * z +
            Mstmp184 * z + Mstmp187 * y + Mstmp20 * M[47] + Mstmp202 * Mstmp3 +
            Mstmp202 * Mstmp4 + Mstmp202 * M[8] + Mstmp230 * M[4] + Mstmp45 * Mstmp99 +
            Mstmp45 * M[38] + x * M[68] + y * M[63] + z * M[62] + M[96];
#pragma omp atomic
  Ms[97] += Mstmp1 * Mstmp209 + Mstmp105 * Mstmp45 + Mstmp107 * Mstmp45 +
            Mstmp109 * Mstmp45 + Mstmp118 * Mstmp38 + Mstmp118 * Mstmp39 +
            Mstmp118 * Mstmp40 + Mstmp118 * M[21] + Mstmp12 * Mstmp170 +
            Mstmp145 * Mstmp20 + Mstmp146 * Mstmp20 + Mstmp147 * Mstmp20 +
            Mstmp16 * Mstmp174 + Mstmp170 * M[12] + Mstmp174 * M[14] + Mstmp187 * z +
            Mstmp188 * z + Mstmp189 * z + Mstmp190 * y + Mstmp20 * M[48] +
            Mstmp209 * M[7] + Mstmp230 * M[5] + Mstmp234 * M[0] + Mstmp45 * M[39] +
            Mstmp58 * Mstmp94 + Mstmp58 * Mstmp95 + Mstmp58 * Mstmp96 + Mstmp58 * M[37] +
            x * M[69] + y * M[64] + z * M[63] + M[97];
#pragma omp atomic
  Ms[98] += Mstmp100 * Mstmp58 + Mstmp102 * Mstmp58 + Mstmp111 * Mstmp45 +
            Mstmp112 * Mstmp45 + Mstmp113 * Mstmp45 + Mstmp13 * Mstmp170 +
            Mstmp137 * Mstmp27 + Mstmp137 * Mstmp28 + Mstmp137 * Mstmp29 +
            Mstmp137 * M[19] + Mstmp148 * Mstmp20 + Mstmp149 * Mstmp20 +
            Mstmp150 * Mstmp20 + Mstmp165 * Mstmp235 + Mstmp17 * Mstmp174 +
            Mstmp170 * M[13] + Mstmp174 * M[15] + Mstmp190 * z + Mstmp191 * z +
            Mstmp192 * z + Mstmp194 * y + Mstmp20 * M[49] + Mstmp209 * Mstmp3 +
            Mstmp209 * M[8] + Mstmp231 * M[3] + Mstmp45 * M[40] + Mstmp58 * Mstmp99 +
            Mstmp58 * M[38] + x * M[70] + y * M[65] + z * M[64] + M[98];
#pragma omp atomic
  Ms[99] += Mstmp1 * Mstmp216 + Mstmp105 * Mstmp58 + Mstmp106 * Mstmp58 +
            Mstmp108 * Mstmp58 + Mstmp137 * Mstmp32 + Mstmp137 * Mstmp33 +
            Mstmp137 * Mstmp35 + Mstmp137 * M[20] + Mstmp15 * Mstmp170 +
            Mstmp151 * Mstmp20 + Mstmp152 * Mstmp20 + Mstmp153 * Mstmp20 +
            Mstmp170 * M[14] + Mstmp194 * z + Mstmp195 * z + Mstmp196 * z + Mstmp197 * y +
            Mstmp2 * Mstmp216 + Mstmp20 * M[50] + Mstmp216 * M[7] + Mstmp231 * M[4] +
            Mstmp58 * M[39] + x * M[71] + y * M[66] + z * M[65] + M[99];
#pragma omp atomic
  Ms[100] += Mstmp111 * Mstmp58 + Mstmp137 * Mstmp38 + Mstmp137 * M[21] +
             Mstmp154 * Mstmp20 + Mstmp170 * M[15] + Mstmp197 * z + Mstmp20 * M[51] +
             Mstmp216 * Mstmp3 + Mstmp216 * M[8] + Mstmp231 * M[5] + Mstmp233 * Mstmp236 +
             Mstmp58 * M[40] + x * M[72] + z * M[66] + M[100];
#pragma omp atomic
  Ms[101] += Mstmp115 * Mstmp45 + Mstmp118 * Mstmp42 + Mstmp118 * M[22] + Mstmp199 * y +
             Mstmp202 * Mstmp5 + Mstmp202 * M[9] + Mstmp237 * M[1] + Mstmp45 * M[41] +
             x * M[73] + y * M[67] + M[101];
#pragma omp atomic
  Ms[102] += Mstmp118 * Mstmp46 + Mstmp118 * Mstmp48 + Mstmp118 * Mstmp50 +
             Mstmp118 * M[23] + Mstmp119 * Mstmp45 + Mstmp121 * Mstmp45 +
             Mstmp123 * Mstmp45 + Mstmp199 * z + Mstmp200 * z + Mstmp201 * z +
             Mstmp202 * Mstmp7 + Mstmp202 * Mstmp9 + Mstmp202 * M[10] + Mstmp203 * y +
             Mstmp237 * M[2] + Mstmp45 * M[42] + x * M[74] + y * M[68] + z * M[67] +
             M[102];
#pragma omp atomic
  Ms[103] += Mstmp10 * Mstmp202 + Mstmp11 * Mstmp202 + Mstmp115 * Mstmp58 +
             Mstmp116 * Mstmp58 + Mstmp117 * Mstmp58 + Mstmp118 * Mstmp52 +
             Mstmp118 * Mstmp54 + Mstmp118 * Mstmp56 + Mstmp118 * M[24] +
             Mstmp125 * Mstmp45 + Mstmp127 * Mstmp45 + Mstmp129 * Mstmp45 +
             Mstmp202 * M[11] + Mstmp203 * z + Mstmp204 * z + Mstmp205 * z +
             Mstmp206 * y + Mstmp209 * Mstmp5 + Mstmp209 * M[9] + Mstmp234 * M[1] +
             Mstmp45 * M[43] + Mstmp58 * M[41] + x * M[75] + y * M[69] + z * M[68] +
             M[103];
#pragma omp atomic
  Ms[104] += Mstmp118 * Mstmp59 + Mstmp118 * Mstmp60 + Mstmp118 * Mstmp61 +
             Mstmp118 * M[25] + Mstmp119 * Mstmp58 + Mstmp120 * Mstmp58 +
             Mstmp122 * Mstmp58 + Mstmp131 * Mstmp45 + Mstmp133 * Mstmp45 +
             Mstmp135 * Mstmp45 + Mstmp137 * Mstmp42 + Mstmp137 * Mstmp43 +
             Mstmp137 * Mstmp44 + Mstmp137 * M[22] + Mstmp206 * z + Mstmp207 * z +
             Mstmp208 * z + Mstmp209 * Mstmp7 + Mstmp209 * M[10] + Mstmp210 * y +
             Mstmp234 * M[2] + Mstmp238 * M[1] + Mstmp45 * M[44] + Mstmp58 * M[42] +
             x * M[76] + y * M[70] + z * M[69] + M[104];
#pragma omp atomic
  Ms[105] += Mstmp10 * Mstmp209 + Mstmp125 * Mstmp58 + Mstmp126 * Mstmp58 +
             Mstmp128 * Mstmp58 + Mstmp137 * Mstmp46 + Mstmp137 * Mstmp47 +
             Mstmp137 * Mstmp49 + Mstmp137 * M[23] + Mstmp138 * Mstmp45 +
             Mstmp139 * Mstmp45 + Mstmp140 * Mstmp45 + Mstmp209 * M[11] + Mstmp210 * z +
             Mstmp211 * z + Mstmp212 * z + Mstmp213 * y + Mstmp216 * Mstmp5 +
             Mstmp216 * Mstmp6 + Mstmp216 * M[9] + Mstmp238 * M[2] + Mstmp45 * M[45] +
             Mstmp58 * M[43] + x * M[77] + y * M[71] + z * M[70] + M[105];
#pragma omp atomic
  Ms[106] += Mstmp131 * Mstmp58 + Mstmp132 * Mstmp58 + Mstmp134 * Mstmp58 +
             Mstmp137 * Mstmp52 + Mstmp137 * Mstmp53 + Mstmp137 * Mstmp55 +
             Mstmp137 * M[24] + Mstmp213 * z + Mstmp214 * z + Mstmp215 * z +
             Mstmp216 * Mstmp7 + Mstmp216 * Mstmp8 + Mstmp216 * M[10] + Mstmp217 * y +
             Mstmp239 * M[1] + Mstmp58 * M[44] + x * M[78] + y * M[72] + z * M[71] +
             M[106];
#pragma omp atomic
  Ms[107] += Mstmp10 * Mstmp216 + Mstmp137 * Mstmp59 + Mstmp137 * M[25] +
             Mstmp138 * Mstmp58 + Mstmp216 * M[11] + Mstmp217 * z + Mstmp239 * M[2] +
             Mstmp58 * M[45] + x * M[79] + z * M[72] + M[107];
#pragma omp atomic
  Ms[108] += Mstmp118 * M[26] + Mstmp202 * M[12] + Mstmp237 * M[3] + Mstmp45 * M[46] +
             y * M[73] + M[108];
#pragma omp atomic
  Ms[109] += Mstmp118 * Mstmp64 + Mstmp118 * M[27] + Mstmp14 * Mstmp202 +
             Mstmp143 * Mstmp45 + Mstmp202 * M[13] + Mstmp218 * z + Mstmp237 * M[4] +
             Mstmp45 * M[47] + y * M[74] + z * M[73] + M[109];
#pragma omp atomic
  Ms[110] += Mstmp118 * Mstmp67 + Mstmp118 * M[28] + Mstmp141 * Mstmp58 +
             Mstmp146 * Mstmp45 + Mstmp16 * Mstmp202 + Mstmp202 * M[14] +
             Mstmp209 * M[12] + Mstmp219 * z + Mstmp234 * M[3] + Mstmp237 * M[5] +
             Mstmp45 * M[48] + Mstmp58 * M[46] + y * M[75] + z * M[74] + M[110];
#pragma omp atomic
  Ms[111] += Mstmp118 * Mstmp70 + Mstmp118 * M[29] + Mstmp137 * Mstmp62 +
             Mstmp137 * M[26] + Mstmp142 * Mstmp58 + Mstmp149 * Mstmp45 +
             Mstmp17 * Mstmp202 + Mstmp202 * M[15] + Mstmp209 * M[13] + Mstmp220 * z +
             Mstmp234 * M[4] + Mstmp238 * M[3] + Mstmp45 * M[49] + Mstmp58 * M[47] +
             y * M[76] + z * M[75] + M[111];
#pragma omp atomic
  Ms[112] += Mstmp118 * Mstmp72 + Mstmp118 * M[30] + Mstmp12 * Mstmp216 +
             Mstmp137 * Mstmp63 + Mstmp137 * M[27] + Mstmp145 * Mstmp58 +
             Mstmp152 * Mstmp45 + Mstmp209 * M[14] + Mstmp216 * M[12] + Mstmp221 * z +
             Mstmp234 * M[5] + Mstmp238 * M[4] + Mstmp45 * M[50] + Mstmp58 * M[48] +
             y * M[77] + z * M[76] + M[112];
#pragma omp atomic
  Ms[113] += Mstmp13 * Mstmp216 + Mstmp137 * Mstmp66 + Mstmp137 * M[28] +
             Mstmp148 * Mstmp58 + Mstmp154 * Mstmp45 + Mstmp209 * M[15] +
             Mstmp216 * M[13] + Mstmp222 * z + Mstmp238 * M[5] + Mstmp239 * M[3] +
             Mstmp45 * M[51] + Mstmp58 * M[49] + y * M[78] + z * M[77] + M[113];
#pragma omp atomic
  Ms[114] += Mstmp137 * Mstmp69 + Mstmp137 * M[29] + Mstmp15 * Mstmp216 +
             Mstmp151 * Mstmp58 + Mstmp216 * M[14] + Mstmp223 * z + Mstmp239 * M[4] +
             Mstmp58 * M[50] + y * M[79] + z * M[78] + M[114];
#pragma omp atomic
  Ms[115] += Mstmp137 * M[30] + Mstmp216 * M[15] + Mstmp239 * M[5] + Mstmp58 * M[51] +
             z * M[79] + M[115];
}

void field_m2_M2L_7(double x, double y, double z, double* M, double* L) {
  double R = sqrt(x * x + y * y + z * z);
  double D[116];
  double Dtmp0   = -1.0 * pow(R, -3.0);
  double Dtmp1   = (x * x);
  double Dtmp2   = pow(R, -5.0);
  double Dtmp3   = 3.0 * Dtmp2;
  double Dtmp4   = x * y;
  double Dtmp5   = x * z;
  double Dtmp6   = (y * y);
  double Dtmp7   = y * z;
  double Dtmp8   = 9.0 * Dtmp2;
  double Dtmp9   = -Dtmp8;
  double Dtmp10  = pow(R, -7.0);
  double Dtmp11  = 15.0 * Dtmp10;
  double Dtmp12  = Dtmp1 * Dtmp11;
  double Dtmp13  = -Dtmp3;
  double Dtmp14  = Dtmp12 + Dtmp13;
  double Dtmp15  = Dtmp11 * Dtmp6;
  double Dtmp16  = Dtmp13 + Dtmp15;
  double Dtmp17  = 1.0 * x;
  double Dtmp18  = Dtmp7 * x;
  double Dtmp19  = (x * x * x * x);
  double Dtmp20  = pow(R, -9.0);
  double Dtmp21  = 105.0 * Dtmp20;
  double Dtmp22  = 90.0 * Dtmp10;
  double Dtmp23  = 45.0 * Dtmp10;
  double Dtmp24  = -Dtmp23;
  double Dtmp25  = Dtmp1 * Dtmp21;
  double Dtmp26  = x * (Dtmp24 + Dtmp25);
  double Dtmp27  = -Dtmp11;
  double Dtmp28  = Dtmp21 * Dtmp6;
  double Dtmp29  = Dtmp24 + Dtmp28;
  double Dtmp30  = Dtmp17 * y;
  double Dtmp31  = Dtmp17 * z;
  double Dtmp32  = (y * y * y * y);
  double Dtmp33  = 225.0 * Dtmp10;
  double Dtmp34  = pow(R, -11.0);
  double Dtmp35  = 945.0 * Dtmp34;
  double Dtmp36  = Dtmp19 * Dtmp35;
  double Dtmp37  = Dtmp1 * Dtmp20;
  double Dtmp38  = 630.0 * Dtmp37;
  double Dtmp39  = Dtmp23 + Dtmp36 - Dtmp38;
  double Dtmp40  = -Dtmp25;
  double Dtmp41  = 315.0 * Dtmp20;
  double Dtmp42  = Dtmp41 * Dtmp6;
  double Dtmp43  = Dtmp1 * Dtmp35;
  double Dtmp44  = Dtmp43 * Dtmp6;
  double Dtmp45  = Dtmp23 + Dtmp44;
  double Dtmp46  = -Dtmp41;
  double Dtmp47  = -Dtmp28;
  double Dtmp48  = Dtmp1 * Dtmp41;
  double Dtmp49  = Dtmp32 * Dtmp35;
  double Dtmp50  = Dtmp20 * Dtmp6;
  double Dtmp51  = 630.0 * Dtmp50;
  double Dtmp52  = Dtmp23 + Dtmp49 - Dtmp51;
  double Dtmp53  = Dtmp35 * Dtmp6;
  double Dtmp54  = Dtmp17 * Dtmp7;
  double Dtmp55  = -Dtmp33;
  double Dtmp56  = (x * x * x * x * x * x);
  double Dtmp57  = pow(R, -13.0);
  double Dtmp58  = 10395.0 * Dtmp57;
  double Dtmp59  = 14175.0 * Dtmp34;
  double Dtmp60  = 1575.0 * Dtmp20;
  double Dtmp61  = Dtmp19 * Dtmp58;
  double Dtmp62  = Dtmp1 * Dtmp34;
  double Dtmp63  = 9450.0 * Dtmp62;
  double Dtmp64  = x * (Dtmp60 + Dtmp61 - Dtmp63);
  double Dtmp65  = 5670.0 * Dtmp62;
  double Dtmp66  = Dtmp24 - Dtmp6 * Dtmp65;
  double Dtmp67  = 945.0 * Dtmp20;
  double Dtmp68  = 2835.0 * Dtmp62;
  double Dtmp69  = -Dtmp68;
  double Dtmp70  = Dtmp34 * Dtmp6;
  double Dtmp71  = 2835.0 * Dtmp70;
  double Dtmp72  = Dtmp1 * Dtmp6;
  double Dtmp73  = Dtmp58 * Dtmp72;
  double Dtmp74  = -Dtmp71 + Dtmp73;
  double Dtmp75  = Dtmp32 * Dtmp58;
  double Dtmp76  = 9450.0 * Dtmp70;
  double Dtmp77  = Dtmp60 + Dtmp75 - Dtmp76;
  double Dtmp78  = 5670.0 * Dtmp70;
  double Dtmp79  = (y * y * y * y * y * y);
  double Dtmp80  = -11025.0 * Dtmp20;
  double Dtmp81  = 135135.0 * pow(R, -15.0);
  double Dtmp82  = Dtmp56 * Dtmp81;
  double Dtmp83  = Dtmp19 * Dtmp57;
  double Dtmp84  = -Dtmp60;
  double Dtmp85  = 42525.0 * Dtmp62 + Dtmp82 - 155925.0 * Dtmp83 + Dtmp84;
  double Dtmp86  = Dtmp19 * Dtmp81;
  double Dtmp87  = Dtmp6 * Dtmp86;
  double Dtmp88  = -Dtmp61 + Dtmp87;
  double Dtmp89  = Dtmp1 * Dtmp57;
  double Dtmp90  = 103950.0 * Dtmp89;
  double Dtmp91  = -Dtmp6 * Dtmp90 + Dtmp84;
  double Dtmp92  = -Dtmp67;
  double Dtmp93  = -62370.0 * Dtmp6 * Dtmp89;
  double Dtmp94  = Dtmp71 + Dtmp93;
  double Dtmp95  = 31185.0 * Dtmp57;
  double Dtmp96  = Dtmp32 * Dtmp81;
  double Dtmp97  = Dtmp1 * Dtmp96;
  double Dtmp98  = Dtmp68 + Dtmp93 + Dtmp97;
  double Dtmp99  = -Dtmp75;
  double Dtmp100 = Dtmp79 * Dtmp81;
  double Dtmp101 = Dtmp32 * Dtmp57;
  double Dtmp102 = Dtmp100 - 155925.0 * Dtmp101 + 42525.0 * Dtmp70 + Dtmp84;
  D[0]           = Dtmp0 + Dtmp1 * Dtmp3;
  D[1]           = Dtmp3 * Dtmp4;
  D[2]           = Dtmp3 * Dtmp5;
  D[3]           = Dtmp0 + Dtmp3 * Dtmp6;
  D[4]           = Dtmp3 * Dtmp7;
  D[5]           = -D[0] - D[3];
  D[6]           = -x * (Dtmp12 + Dtmp9);
  D[7]           = -Dtmp14 * y;
  D[8]           = -Dtmp14 * z;
  D[9]           = -Dtmp16 * Dtmp17;
  D[10]          = -Dtmp11 * Dtmp18;
  D[11]          = -D[6] - D[9];
  D[12]          = -y * (Dtmp15 + Dtmp9);
  D[13]          = -Dtmp16 * z;
  D[14]          = -D[7] - D[12];
  D[15]          = -D[8] - D[13];
  D[16]          = -Dtmp1 * Dtmp22 + Dtmp19 * Dtmp21 + Dtmp8;
  D[17]          = Dtmp26 * y;
  D[18]          = Dtmp26 * z;
  D[19]          = -Dtmp12 - Dtmp15 + Dtmp25 * Dtmp6 + Dtmp3;
  D[20]          = Dtmp7 * (Dtmp25 + Dtmp27);
  D[21]          = -D[16] - D[19];
  D[22]          = Dtmp29 * Dtmp30;
  D[23]          = Dtmp31 * (Dtmp27 + Dtmp28);
  D[24]          = -D[17] - D[22];
  D[25]          = -D[18] - D[23];
  D[26]          = Dtmp21 * Dtmp32 - Dtmp22 * Dtmp6 + Dtmp8;
  D[27]          = Dtmp29 * Dtmp7;
  D[28]          = -D[19] - D[26];
  D[29]          = -D[20] - D[27];
  D[30]          = -D[21] - D[28];
  D[31]          = -x * (Dtmp33 + Dtmp36 - 1050.0 * Dtmp37);
  D[32]          = -Dtmp39 * y;
  D[33]          = -Dtmp39 * z;
  D[34]          = -x * (Dtmp40 - Dtmp42 + Dtmp45);
  D[35]          = -Dtmp18 * (Dtmp43 + Dtmp46);
  D[36]          = -D[31] - D[34];
  D[37]          = -y * (Dtmp45 + Dtmp47 - Dtmp48);
  D[38]          = -z * (Dtmp11 + Dtmp40 + Dtmp44 + Dtmp47);
  D[39]          = -D[32] - D[37];
  D[40]          = -D[33] - D[38];
  D[41]          = -Dtmp17 * Dtmp52;
  D[42]          = -Dtmp54 * (Dtmp46 + Dtmp53);
  D[43]          = -D[34] - D[41];
  D[44]          = -D[35] - D[42];
  D[45]          = -D[36] - D[43];
  D[46]          = -y * (Dtmp33 + Dtmp49 - 1050.0 * Dtmp50);
  D[47]          = -Dtmp52 * z;
  D[48]          = -D[37] - D[46];
  D[49]          = -D[38] - D[47];
  D[50]          = -D[39] - D[48];
  D[51]          = -D[40] - D[49];
  D[52]          = -Dtmp19 * Dtmp59 + 4725.0 * Dtmp37 + Dtmp55 + Dtmp56 * Dtmp58;
  D[53]          = Dtmp64 * y;
  D[54]          = Dtmp64 * z;
  D[55]          = -Dtmp36 + Dtmp38 + Dtmp42 + Dtmp6 * Dtmp61 + Dtmp66;
  D[56]          = Dtmp7 * (Dtmp41 + Dtmp61 - Dtmp65);
  D[57]          = -D[52] - D[55];
  D[58]          = Dtmp4 * (Dtmp67 + Dtmp69 + Dtmp74);
  D[59]          = Dtmp5 * (Dtmp41 - Dtmp43 + Dtmp74);
  D[60]          = -D[53] - D[58];
  D[61]          = -D[54] - D[59];
  D[62]          = Dtmp1 * Dtmp75 + Dtmp48 - Dtmp49 + Dtmp51 + Dtmp66;
  D[63]          = Dtmp7 * (Dtmp41 - Dtmp53 + Dtmp69 + Dtmp73);
  D[64]          = -D[55] - D[62];
  D[65]          = -D[56] - D[63];
  D[66]          = -D[57] - D[64];
  D[67]          = Dtmp30 * Dtmp77;
  D[68]          = Dtmp31 * (Dtmp41 + Dtmp75 - Dtmp78);
  D[69]          = -D[58] - D[67];
  D[70]          = -D[59] - D[68];
  D[71]          = -D[60] - D[69];
  D[72]          = -D[61] - D[70];
  D[73]          = -Dtmp32 * Dtmp59 + 4725.0 * Dtmp50 + Dtmp55 + Dtmp58 * Dtmp79;
  D[74]          = Dtmp7 * Dtmp77;
  D[75]          = -D[62] - D[73];
  D[76]          = -D[63] - D[74];
  D[77]          = -D[64] - D[75];
  D[78]          = -D[65] - D[76];
  D[79]          = -D[66] - D[77];
  D[80]          = -x * (99225.0 * Dtmp62 + Dtmp80 + Dtmp82 - 218295.0 * Dtmp83);
  D[81]          = -Dtmp85 * y;
  D[82]          = -Dtmp85 * z;
  D[83]          = -x * (Dtmp59 * Dtmp6 + Dtmp63 + Dtmp88 + Dtmp91);
  D[84]          = -Dtmp18 * (Dtmp59 + Dtmp86 - Dtmp90);
  D[85]          = -D[80] - D[83];
  D[86]          = -y * (17010.0 * Dtmp62 - 31185.0 * Dtmp83 + Dtmp87 + Dtmp92 + Dtmp94);
  D[87]          = -z * (Dtmp46 + Dtmp65 + Dtmp88 + Dtmp94);
  D[88]          = -D[81] - D[86];
  D[89]          = -D[82] - D[87];
  D[90]          = -x * (-Dtmp32 * Dtmp95 + 17010.0 * Dtmp70 + Dtmp92 + Dtmp98);
  D[91] =
        -Dtmp18 * (8505.0 * Dtmp34 - Dtmp6 * Dtmp95 + Dtmp72 * Dtmp81 - 31185.0 * Dtmp89);
  D[92]  = -D[83] - D[90];
  D[93]  = -D[84] - D[91];
  D[94]  = -D[85] - D[92];
  D[95]  = -y * (Dtmp1 * Dtmp59 + Dtmp76 + Dtmp91 + Dtmp97 + Dtmp99);
  D[96]  = -z * (Dtmp46 + Dtmp78 + Dtmp98 + Dtmp99);
  D[97]  = -D[86] - D[95];
  D[98]  = -D[87] - D[96];
  D[99]  = -D[88] - D[97];
  D[100] = -D[89] - D[98];
  D[101] = -Dtmp102 * Dtmp17;
  D[102] = -Dtmp54 * (-103950.0 * Dtmp57 * Dtmp6 + Dtmp59 + Dtmp96);
  D[103] = -D[90] - D[101];
  D[104] = -D[91] - D[102];
  D[105] = -D[92] - D[103];
  D[106] = -D[93] - D[104];
  D[107] = -D[94] - D[105];
  D[108] = -y * (Dtmp100 - 218295.0 * Dtmp101 + 99225.0 * Dtmp70 + Dtmp80);
  D[109] = -Dtmp102 * z;
  D[110] = -D[95] - D[108];
  D[111] = -D[96] - D[109];
  D[112] = -D[97] - D[110];
  D[113] = -D[98] - D[111];
  D[114] = -D[99] - D[112];
  D[115] = -D[100] - D[113];
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
          D[112] * M[112] + D[113] * M[113] + D[114] * M[114] + D[115] * M[115];
#pragma omp atomic
  L[1] += D[6] * M[0] + D[7] * M[1] + D[8] * M[2] + D[9] * M[3] + D[10] * M[4] +
          D[11] * M[5] + D[16] * M[6] + D[17] * M[7] + D[18] * M[8] + D[19] * M[9] +
          D[20] * M[10] + D[21] * M[11] + D[22] * M[12] + D[23] * M[13] + D[24] * M[14] +
          D[25] * M[15] + D[31] * M[16] + D[32] * M[17] + D[33] * M[18] + D[34] * M[19] +
          D[35] * M[20] + D[36] * M[21] + D[37] * M[22] + D[38] * M[23] + D[39] * M[24] +
          D[40] * M[25] + D[41] * M[26] + D[42] * M[27] + D[43] * M[28] + D[44] * M[29] +
          D[45] * M[30] + D[52] * M[31] + D[53] * M[32] + D[54] * M[33] + D[55] * M[34] +
          D[56] * M[35] + D[57] * M[36] + D[58] * M[37] + D[59] * M[38] + D[60] * M[39] +
          D[61] * M[40] + D[62] * M[41] + D[63] * M[42] + D[64] * M[43] + D[65] * M[44] +
          D[66] * M[45] + D[67] * M[46] + D[68] * M[47] + D[69] * M[48] + D[70] * M[49] +
          D[71] * M[50] + D[72] * M[51] + D[80] * M[52] + D[81] * M[53] + D[82] * M[54] +
          D[83] * M[55] + D[84] * M[56] + D[85] * M[57] + D[86] * M[58] + D[87] * M[59] +
          D[88] * M[60] + D[89] * M[61] + D[90] * M[62] + D[91] * M[63] + D[92] * M[64] +
          D[93] * M[65] + D[94] * M[66] + D[95] * M[67] + D[96] * M[68] + D[97] * M[69] +
          D[98] * M[70] + D[99] * M[71] + D[100] * M[72] + D[101] * M[73] +
          D[102] * M[74] + D[103] * M[75] + D[104] * M[76] + D[105] * M[77] +
          D[106] * M[78] + D[107] * M[79];
#pragma omp atomic
  L[2] += D[7] * M[0] + D[9] * M[1] + D[10] * M[2] + D[12] * M[3] + D[13] * M[4] +
          D[14] * M[5] + D[17] * M[6] + D[19] * M[7] + D[20] * M[8] + D[22] * M[9] +
          D[23] * M[10] + D[24] * M[11] + D[26] * M[12] + D[27] * M[13] + D[28] * M[14] +
          D[29] * M[15] + D[32] * M[16] + D[34] * M[17] + D[35] * M[18] + D[37] * M[19] +
          D[38] * M[20] + D[39] * M[21] + D[41] * M[22] + D[42] * M[23] + D[43] * M[24] +
          D[44] * M[25] + D[46] * M[26] + D[47] * M[27] + D[48] * M[28] + D[49] * M[29] +
          D[50] * M[30] + D[53] * M[31] + D[55] * M[32] + D[56] * M[33] + D[58] * M[34] +
          D[59] * M[35] + D[60] * M[36] + D[62] * M[37] + D[63] * M[38] + D[64] * M[39] +
          D[65] * M[40] + D[67] * M[41] + D[68] * M[42] + D[69] * M[43] + D[70] * M[44] +
          D[71] * M[45] + D[73] * M[46] + D[74] * M[47] + D[75] * M[48] + D[76] * M[49] +
          D[77] * M[50] + D[78] * M[51] + D[81] * M[52] + D[83] * M[53] + D[84] * M[54] +
          D[86] * M[55] + D[87] * M[56] + D[88] * M[57] + D[90] * M[58] + D[91] * M[59] +
          D[92] * M[60] + D[93] * M[61] + D[95] * M[62] + D[96] * M[63] + D[97] * M[64] +
          D[98] * M[65] + D[99] * M[66] + D[101] * M[67] + D[102] * M[68] +
          D[103] * M[69] + D[104] * M[70] + D[105] * M[71] + D[106] * M[72] +
          D[108] * M[73] + D[109] * M[74] + D[110] * M[75] + D[111] * M[76] +
          D[112] * M[77] + D[113] * M[78] + D[114] * M[79];
#pragma omp atomic
  L[3] += D[8] * M[0] + D[10] * M[1] + D[11] * M[2] + D[13] * M[3] + D[14] * M[4] +
          D[15] * M[5] + D[18] * M[6] + D[20] * M[7] + D[21] * M[8] + D[23] * M[9] +
          D[24] * M[10] + D[25] * M[11] + D[27] * M[12] + D[28] * M[13] + D[29] * M[14] +
          D[30] * M[15] + D[33] * M[16] + D[35] * M[17] + D[36] * M[18] + D[38] * M[19] +
          D[39] * M[20] + D[40] * M[21] + D[42] * M[22] + D[43] * M[23] + D[44] * M[24] +
          D[45] * M[25] + D[47] * M[26] + D[48] * M[27] + D[49] * M[28] + D[50] * M[29] +
          D[51] * M[30] + D[54] * M[31] + D[56] * M[32] + D[57] * M[33] + D[59] * M[34] +
          D[60] * M[35] + D[61] * M[36] + D[63] * M[37] + D[64] * M[38] + D[65] * M[39] +
          D[66] * M[40] + D[68] * M[41] + D[69] * M[42] + D[70] * M[43] + D[71] * M[44] +
          D[72] * M[45] + D[74] * M[46] + D[75] * M[47] + D[76] * M[48] + D[77] * M[49] +
          D[78] * M[50] + D[79] * M[51] + D[82] * M[52] + D[84] * M[53] + D[85] * M[54] +
          D[87] * M[55] + D[88] * M[56] + D[89] * M[57] + D[91] * M[58] + D[92] * M[59] +
          D[93] * M[60] + D[94] * M[61] + D[96] * M[62] + D[97] * M[63] + D[98] * M[64] +
          D[99] * M[65] + D[100] * M[66] + D[102] * M[67] + D[103] * M[68] +
          D[104] * M[69] + D[105] * M[70] + D[106] * M[71] + D[107] * M[72] +
          D[109] * M[73] + D[110] * M[74] + D[111] * M[75] + D[112] * M[76] +
          D[113] * M[77] + D[114] * M[78] + D[115] * M[79];
#pragma omp atomic
  L[4] += D[16] * M[0] + D[17] * M[1] + D[18] * M[2] + D[19] * M[3] + D[20] * M[4] +
          D[21] * M[5] + D[31] * M[6] + D[32] * M[7] + D[33] * M[8] + D[34] * M[9] +
          D[35] * M[10] + D[36] * M[11] + D[37] * M[12] + D[38] * M[13] + D[39] * M[14] +
          D[40] * M[15] + D[52] * M[16] + D[53] * M[17] + D[54] * M[18] + D[55] * M[19] +
          D[56] * M[20] + D[57] * M[21] + D[58] * M[22] + D[59] * M[23] + D[60] * M[24] +
          D[61] * M[25] + D[62] * M[26] + D[63] * M[27] + D[64] * M[28] + D[65] * M[29] +
          D[66] * M[30] + D[80] * M[31] + D[81] * M[32] + D[82] * M[33] + D[83] * M[34] +
          D[84] * M[35] + D[85] * M[36] + D[86] * M[37] + D[87] * M[38] + D[88] * M[39] +
          D[89] * M[40] + D[90] * M[41] + D[91] * M[42] + D[92] * M[43] + D[93] * M[44] +
          D[94] * M[45] + D[95] * M[46] + D[96] * M[47] + D[97] * M[48] + D[98] * M[49] +
          D[99] * M[50] + D[100] * M[51];
#pragma omp atomic
  L[5] += D[17] * M[0] + D[19] * M[1] + D[20] * M[2] + D[22] * M[3] + D[23] * M[4] +
          D[24] * M[5] + D[32] * M[6] + D[34] * M[7] + D[35] * M[8] + D[37] * M[9] +
          D[38] * M[10] + D[39] * M[11] + D[41] * M[12] + D[42] * M[13] + D[43] * M[14] +
          D[44] * M[15] + D[53] * M[16] + D[55] * M[17] + D[56] * M[18] + D[58] * M[19] +
          D[59] * M[20] + D[60] * M[21] + D[62] * M[22] + D[63] * M[23] + D[64] * M[24] +
          D[65] * M[25] + D[67] * M[26] + D[68] * M[27] + D[69] * M[28] + D[70] * M[29] +
          D[71] * M[30] + D[81] * M[31] + D[83] * M[32] + D[84] * M[33] + D[86] * M[34] +
          D[87] * M[35] + D[88] * M[36] + D[90] * M[37] + D[91] * M[38] + D[92] * M[39] +
          D[93] * M[40] + D[95] * M[41] + D[96] * M[42] + D[97] * M[43] + D[98] * M[44] +
          D[99] * M[45] + D[101] * M[46] + D[102] * M[47] + D[103] * M[48] +
          D[104] * M[49] + D[105] * M[50] + D[106] * M[51];
#pragma omp atomic
  L[6] += D[18] * M[0] + D[20] * M[1] + D[21] * M[2] + D[23] * M[3] + D[24] * M[4] +
          D[25] * M[5] + D[33] * M[6] + D[35] * M[7] + D[36] * M[8] + D[38] * M[9] +
          D[39] * M[10] + D[40] * M[11] + D[42] * M[12] + D[43] * M[13] + D[44] * M[14] +
          D[45] * M[15] + D[54] * M[16] + D[56] * M[17] + D[57] * M[18] + D[59] * M[19] +
          D[60] * M[20] + D[61] * M[21] + D[63] * M[22] + D[64] * M[23] + D[65] * M[24] +
          D[66] * M[25] + D[68] * M[26] + D[69] * M[27] + D[70] * M[28] + D[71] * M[29] +
          D[72] * M[30] + D[82] * M[31] + D[84] * M[32] + D[85] * M[33] + D[87] * M[34] +
          D[88] * M[35] + D[89] * M[36] + D[91] * M[37] + D[92] * M[38] + D[93] * M[39] +
          D[94] * M[40] + D[96] * M[41] + D[97] * M[42] + D[98] * M[43] + D[99] * M[44] +
          D[100] * M[45] + D[102] * M[46] + D[103] * M[47] + D[104] * M[48] +
          D[105] * M[49] + D[106] * M[50] + D[107] * M[51];
#pragma omp atomic
  L[7] += D[19] * M[0] + D[22] * M[1] + D[23] * M[2] + D[26] * M[3] + D[27] * M[4] +
          D[28] * M[5] + D[34] * M[6] + D[37] * M[7] + D[38] * M[8] + D[41] * M[9] +
          D[42] * M[10] + D[43] * M[11] + D[46] * M[12] + D[47] * M[13] + D[48] * M[14] +
          D[49] * M[15] + D[55] * M[16] + D[58] * M[17] + D[59] * M[18] + D[62] * M[19] +
          D[63] * M[20] + D[64] * M[21] + D[67] * M[22] + D[68] * M[23] + D[69] * M[24] +
          D[70] * M[25] + D[73] * M[26] + D[74] * M[27] + D[75] * M[28] + D[76] * M[29] +
          D[77] * M[30] + D[83] * M[31] + D[86] * M[32] + D[87] * M[33] + D[90] * M[34] +
          D[91] * M[35] + D[92] * M[36] + D[95] * M[37] + D[96] * M[38] + D[97] * M[39] +
          D[98] * M[40] + D[101] * M[41] + D[102] * M[42] + D[103] * M[43] +
          D[104] * M[44] + D[105] * M[45] + D[108] * M[46] + D[109] * M[47] +
          D[110] * M[48] + D[111] * M[49] + D[112] * M[50] + D[113] * M[51];
#pragma omp atomic
  L[8] += D[20] * M[0] + D[23] * M[1] + D[24] * M[2] + D[27] * M[3] + D[28] * M[4] +
          D[29] * M[5] + D[35] * M[6] + D[38] * M[7] + D[39] * M[8] + D[42] * M[9] +
          D[43] * M[10] + D[44] * M[11] + D[47] * M[12] + D[48] * M[13] + D[49] * M[14] +
          D[50] * M[15] + D[56] * M[16] + D[59] * M[17] + D[60] * M[18] + D[63] * M[19] +
          D[64] * M[20] + D[65] * M[21] + D[68] * M[22] + D[69] * M[23] + D[70] * M[24] +
          D[71] * M[25] + D[74] * M[26] + D[75] * M[27] + D[76] * M[28] + D[77] * M[29] +
          D[78] * M[30] + D[84] * M[31] + D[87] * M[32] + D[88] * M[33] + D[91] * M[34] +
          D[92] * M[35] + D[93] * M[36] + D[96] * M[37] + D[97] * M[38] + D[98] * M[39] +
          D[99] * M[40] + D[102] * M[41] + D[103] * M[42] + D[104] * M[43] +
          D[105] * M[44] + D[106] * M[45] + D[109] * M[46] + D[110] * M[47] +
          D[111] * M[48] + D[112] * M[49] + D[113] * M[50] + D[114] * M[51];
#pragma omp atomic
  L[9] += D[21] * M[0] + D[24] * M[1] + D[25] * M[2] + D[28] * M[3] + D[29] * M[4] +
          D[30] * M[5] + D[36] * M[6] + D[39] * M[7] + D[40] * M[8] + D[43] * M[9] +
          D[44] * M[10] + D[45] * M[11] + D[48] * M[12] + D[49] * M[13] + D[50] * M[14] +
          D[51] * M[15] + D[57] * M[16] + D[60] * M[17] + D[61] * M[18] + D[64] * M[19] +
          D[65] * M[20] + D[66] * M[21] + D[69] * M[22] + D[70] * M[23] + D[71] * M[24] +
          D[72] * M[25] + D[75] * M[26] + D[76] * M[27] + D[77] * M[28] + D[78] * M[29] +
          D[79] * M[30] + D[85] * M[31] + D[88] * M[32] + D[89] * M[33] + D[92] * M[34] +
          D[93] * M[35] + D[94] * M[36] + D[97] * M[37] + D[98] * M[38] + D[99] * M[39] +
          D[100] * M[40] + D[103] * M[41] + D[104] * M[42] + D[105] * M[43] +
          D[106] * M[44] + D[107] * M[45] + D[110] * M[46] + D[111] * M[47] +
          D[112] * M[48] + D[113] * M[49] + D[114] * M[50] + D[115] * M[51];
#pragma omp atomic
  L[10] += D[31] * M[0] + D[32] * M[1] + D[33] * M[2] + D[34] * M[3] + D[35] * M[4] +
           D[36] * M[5] + D[52] * M[6] + D[53] * M[7] + D[54] * M[8] + D[55] * M[9] +
           D[56] * M[10] + D[57] * M[11] + D[58] * M[12] + D[59] * M[13] + D[60] * M[14] +
           D[61] * M[15] + D[80] * M[16] + D[81] * M[17] + D[82] * M[18] + D[83] * M[19] +
           D[84] * M[20] + D[85] * M[21] + D[86] * M[22] + D[87] * M[23] + D[88] * M[24] +
           D[89] * M[25] + D[90] * M[26] + D[91] * M[27] + D[92] * M[28] + D[93] * M[29] +
           D[94] * M[30];
#pragma omp atomic
  L[11] += D[32] * M[0] + D[34] * M[1] + D[35] * M[2] + D[37] * M[3] + D[38] * M[4] +
           D[39] * M[5] + D[53] * M[6] + D[55] * M[7] + D[56] * M[8] + D[58] * M[9] +
           D[59] * M[10] + D[60] * M[11] + D[62] * M[12] + D[63] * M[13] + D[64] * M[14] +
           D[65] * M[15] + D[81] * M[16] + D[83] * M[17] + D[84] * M[18] + D[86] * M[19] +
           D[87] * M[20] + D[88] * M[21] + D[90] * M[22] + D[91] * M[23] + D[92] * M[24] +
           D[93] * M[25] + D[95] * M[26] + D[96] * M[27] + D[97] * M[28] + D[98] * M[29] +
           D[99] * M[30];
#pragma omp atomic
  L[12] += D[33] * M[0] + D[35] * M[1] + D[36] * M[2] + D[38] * M[3] + D[39] * M[4] +
           D[40] * M[5] + D[54] * M[6] + D[56] * M[7] + D[57] * M[8] + D[59] * M[9] +
           D[60] * M[10] + D[61] * M[11] + D[63] * M[12] + D[64] * M[13] + D[65] * M[14] +
           D[66] * M[15] + D[82] * M[16] + D[84] * M[17] + D[85] * M[18] + D[87] * M[19] +
           D[88] * M[20] + D[89] * M[21] + D[91] * M[22] + D[92] * M[23] + D[93] * M[24] +
           D[94] * M[25] + D[96] * M[26] + D[97] * M[27] + D[98] * M[28] + D[99] * M[29] +
           D[100] * M[30];
#pragma omp atomic
  L[13] += D[34] * M[0] + D[37] * M[1] + D[38] * M[2] + D[41] * M[3] + D[42] * M[4] +
           D[43] * M[5] + D[55] * M[6] + D[58] * M[7] + D[59] * M[8] + D[62] * M[9] +
           D[63] * M[10] + D[64] * M[11] + D[67] * M[12] + D[68] * M[13] + D[69] * M[14] +
           D[70] * M[15] + D[83] * M[16] + D[86] * M[17] + D[87] * M[18] + D[90] * M[19] +
           D[91] * M[20] + D[92] * M[21] + D[95] * M[22] + D[96] * M[23] + D[97] * M[24] +
           D[98] * M[25] + D[101] * M[26] + D[102] * M[27] + D[103] * M[28] +
           D[104] * M[29] + D[105] * M[30];
#pragma omp atomic
  L[14] += D[35] * M[0] + D[38] * M[1] + D[39] * M[2] + D[42] * M[3] + D[43] * M[4] +
           D[44] * M[5] + D[56] * M[6] + D[59] * M[7] + D[60] * M[8] + D[63] * M[9] +
           D[64] * M[10] + D[65] * M[11] + D[68] * M[12] + D[69] * M[13] + D[70] * M[14] +
           D[71] * M[15] + D[84] * M[16] + D[87] * M[17] + D[88] * M[18] + D[91] * M[19] +
           D[92] * M[20] + D[93] * M[21] + D[96] * M[22] + D[97] * M[23] + D[98] * M[24] +
           D[99] * M[25] + D[102] * M[26] + D[103] * M[27] + D[104] * M[28] +
           D[105] * M[29] + D[106] * M[30];
#pragma omp atomic
  L[15] += D[36] * M[0] + D[39] * M[1] + D[40] * M[2] + D[43] * M[3] + D[44] * M[4] +
           D[45] * M[5] + D[57] * M[6] + D[60] * M[7] + D[61] * M[8] + D[64] * M[9] +
           D[65] * M[10] + D[66] * M[11] + D[69] * M[12] + D[70] * M[13] + D[71] * M[14] +
           D[72] * M[15] + D[85] * M[16] + D[88] * M[17] + D[89] * M[18] + D[92] * M[19] +
           D[93] * M[20] + D[94] * M[21] + D[97] * M[22] + D[98] * M[23] + D[99] * M[24] +
           D[100] * M[25] + D[103] * M[26] + D[104] * M[27] + D[105] * M[28] +
           D[106] * M[29] + D[107] * M[30];
#pragma omp atomic
  L[16] += D[37] * M[0] + D[41] * M[1] + D[42] * M[2] + D[46] * M[3] + D[47] * M[4] +
           D[48] * M[5] + D[58] * M[6] + D[62] * M[7] + D[63] * M[8] + D[67] * M[9] +
           D[68] * M[10] + D[69] * M[11] + D[73] * M[12] + D[74] * M[13] + D[75] * M[14] +
           D[76] * M[15] + D[86] * M[16] + D[90] * M[17] + D[91] * M[18] + D[95] * M[19] +
           D[96] * M[20] + D[97] * M[21] + D[101] * M[22] + D[102] * M[23] +
           D[103] * M[24] + D[104] * M[25] + D[108] * M[26] + D[109] * M[27] +
           D[110] * M[28] + D[111] * M[29] + D[112] * M[30];
#pragma omp atomic
  L[17] += D[38] * M[0] + D[42] * M[1] + D[43] * M[2] + D[47] * M[3] + D[48] * M[4] +
           D[49] * M[5] + D[59] * M[6] + D[63] * M[7] + D[64] * M[8] + D[68] * M[9] +
           D[69] * M[10] + D[70] * M[11] + D[74] * M[12] + D[75] * M[13] + D[76] * M[14] +
           D[77] * M[15] + D[87] * M[16] + D[91] * M[17] + D[92] * M[18] + D[96] * M[19] +
           D[97] * M[20] + D[98] * M[21] + D[102] * M[22] + D[103] * M[23] +
           D[104] * M[24] + D[105] * M[25] + D[109] * M[26] + D[110] * M[27] +
           D[111] * M[28] + D[112] * M[29] + D[113] * M[30];
#pragma omp atomic
  L[18] += D[39] * M[0] + D[43] * M[1] + D[44] * M[2] + D[48] * M[3] + D[49] * M[4] +
           D[50] * M[5] + D[60] * M[6] + D[64] * M[7] + D[65] * M[8] + D[69] * M[9] +
           D[70] * M[10] + D[71] * M[11] + D[75] * M[12] + D[76] * M[13] + D[77] * M[14] +
           D[78] * M[15] + D[88] * M[16] + D[92] * M[17] + D[93] * M[18] + D[97] * M[19] +
           D[98] * M[20] + D[99] * M[21] + D[103] * M[22] + D[104] * M[23] +
           D[105] * M[24] + D[106] * M[25] + D[110] * M[26] + D[111] * M[27] +
           D[112] * M[28] + D[113] * M[29] + D[114] * M[30];
#pragma omp atomic
  L[19] += D[40] * M[0] + D[44] * M[1] + D[45] * M[2] + D[49] * M[3] + D[50] * M[4] +
           D[51] * M[5] + D[61] * M[6] + D[65] * M[7] + D[66] * M[8] + D[70] * M[9] +
           D[71] * M[10] + D[72] * M[11] + D[76] * M[12] + D[77] * M[13] + D[78] * M[14] +
           D[79] * M[15] + D[89] * M[16] + D[93] * M[17] + D[94] * M[18] + D[98] * M[19] +
           D[99] * M[20] + D[100] * M[21] + D[104] * M[22] + D[105] * M[23] +
           D[106] * M[24] + D[107] * M[25] + D[111] * M[26] + D[112] * M[27] +
           D[113] * M[28] + D[114] * M[29] + D[115] * M[30];
#pragma omp atomic
  L[20] += D[52] * M[0] + D[53] * M[1] + D[54] * M[2] + D[55] * M[3] + D[56] * M[4] +
           D[57] * M[5] + D[80] * M[6] + D[81] * M[7] + D[82] * M[8] + D[83] * M[9] +
           D[84] * M[10] + D[85] * M[11] + D[86] * M[12] + D[87] * M[13] + D[88] * M[14] +
           D[89] * M[15];
#pragma omp atomic
  L[21] += D[53] * M[0] + D[55] * M[1] + D[56] * M[2] + D[58] * M[3] + D[59] * M[4] +
           D[60] * M[5] + D[81] * M[6] + D[83] * M[7] + D[84] * M[8] + D[86] * M[9] +
           D[87] * M[10] + D[88] * M[11] + D[90] * M[12] + D[91] * M[13] + D[92] * M[14] +
           D[93] * M[15];
#pragma omp atomic
  L[22] += D[54] * M[0] + D[56] * M[1] + D[57] * M[2] + D[59] * M[3] + D[60] * M[4] +
           D[61] * M[5] + D[82] * M[6] + D[84] * M[7] + D[85] * M[8] + D[87] * M[9] +
           D[88] * M[10] + D[89] * M[11] + D[91] * M[12] + D[92] * M[13] + D[93] * M[14] +
           D[94] * M[15];
#pragma omp atomic
  L[23] += D[55] * M[0] + D[58] * M[1] + D[59] * M[2] + D[62] * M[3] + D[63] * M[4] +
           D[64] * M[5] + D[83] * M[6] + D[86] * M[7] + D[87] * M[8] + D[90] * M[9] +
           D[91] * M[10] + D[92] * M[11] + D[95] * M[12] + D[96] * M[13] + D[97] * M[14] +
           D[98] * M[15];
#pragma omp atomic
  L[24] += D[56] * M[0] + D[59] * M[1] + D[60] * M[2] + D[63] * M[3] + D[64] * M[4] +
           D[65] * M[5] + D[84] * M[6] + D[87] * M[7] + D[88] * M[8] + D[91] * M[9] +
           D[92] * M[10] + D[93] * M[11] + D[96] * M[12] + D[97] * M[13] + D[98] * M[14] +
           D[99] * M[15];
#pragma omp atomic
  L[25] += D[57] * M[0] + D[60] * M[1] + D[61] * M[2] + D[64] * M[3] + D[65] * M[4] +
           D[66] * M[5] + D[85] * M[6] + D[88] * M[7] + D[89] * M[8] + D[92] * M[9] +
           D[93] * M[10] + D[94] * M[11] + D[97] * M[12] + D[98] * M[13] + D[99] * M[14] +
           D[100] * M[15];
#pragma omp atomic
  L[26] += D[58] * M[0] + D[62] * M[1] + D[63] * M[2] + D[67] * M[3] + D[68] * M[4] +
           D[69] * M[5] + D[86] * M[6] + D[90] * M[7] + D[91] * M[8] + D[95] * M[9] +
           D[96] * M[10] + D[97] * M[11] + D[101] * M[12] + D[102] * M[13] +
           D[103] * M[14] + D[104] * M[15];
#pragma omp atomic
  L[27] += D[59] * M[0] + D[63] * M[1] + D[64] * M[2] + D[68] * M[3] + D[69] * M[4] +
           D[70] * M[5] + D[87] * M[6] + D[91] * M[7] + D[92] * M[8] + D[96] * M[9] +
           D[97] * M[10] + D[98] * M[11] + D[102] * M[12] + D[103] * M[13] +
           D[104] * M[14] + D[105] * M[15];
#pragma omp atomic
  L[28] += D[60] * M[0] + D[64] * M[1] + D[65] * M[2] + D[69] * M[3] + D[70] * M[4] +
           D[71] * M[5] + D[88] * M[6] + D[92] * M[7] + D[93] * M[8] + D[97] * M[9] +
           D[98] * M[10] + D[99] * M[11] + D[103] * M[12] + D[104] * M[13] +
           D[105] * M[14] + D[106] * M[15];
#pragma omp atomic
  L[29] += D[61] * M[0] + D[65] * M[1] + D[66] * M[2] + D[70] * M[3] + D[71] * M[4] +
           D[72] * M[5] + D[89] * M[6] + D[93] * M[7] + D[94] * M[8] + D[98] * M[9] +
           D[99] * M[10] + D[100] * M[11] + D[104] * M[12] + D[105] * M[13] +
           D[106] * M[14] + D[107] * M[15];
#pragma omp atomic
  L[30] += D[62] * M[0] + D[67] * M[1] + D[68] * M[2] + D[73] * M[3] + D[74] * M[4] +
           D[75] * M[5] + D[90] * M[6] + D[95] * M[7] + D[96] * M[8] + D[101] * M[9] +
           D[102] * M[10] + D[103] * M[11] + D[108] * M[12] + D[109] * M[13] +
           D[110] * M[14] + D[111] * M[15];
#pragma omp atomic
  L[31] += D[63] * M[0] + D[68] * M[1] + D[69] * M[2] + D[74] * M[3] + D[75] * M[4] +
           D[76] * M[5] + D[91] * M[6] + D[96] * M[7] + D[97] * M[8] + D[102] * M[9] +
           D[103] * M[10] + D[104] * M[11] + D[109] * M[12] + D[110] * M[13] +
           D[111] * M[14] + D[112] * M[15];
#pragma omp atomic
  L[32] += D[64] * M[0] + D[69] * M[1] + D[70] * M[2] + D[75] * M[3] + D[76] * M[4] +
           D[77] * M[5] + D[92] * M[6] + D[97] * M[7] + D[98] * M[8] + D[103] * M[9] +
           D[104] * M[10] + D[105] * M[11] + D[110] * M[12] + D[111] * M[13] +
           D[112] * M[14] + D[113] * M[15];
#pragma omp atomic
  L[33] += D[65] * M[0] + D[70] * M[1] + D[71] * M[2] + D[76] * M[3] + D[77] * M[4] +
           D[78] * M[5] + D[93] * M[6] + D[98] * M[7] + D[99] * M[8] + D[104] * M[9] +
           D[105] * M[10] + D[106] * M[11] + D[111] * M[12] + D[112] * M[13] +
           D[113] * M[14] + D[114] * M[15];
#pragma omp atomic
  L[34] += D[66] * M[0] + D[71] * M[1] + D[72] * M[2] + D[77] * M[3] + D[78] * M[4] +
           D[79] * M[5] + D[94] * M[6] + D[99] * M[7] + D[100] * M[8] + D[105] * M[9] +
           D[106] * M[10] + D[107] * M[11] + D[112] * M[12] + D[113] * M[13] +
           D[114] * M[14] + D[115] * M[15];
#pragma omp atomic
  L[35] += D[80] * M[0] + D[81] * M[1] + D[82] * M[2] + D[83] * M[3] + D[84] * M[4] +
           D[85] * M[5];
#pragma omp atomic
  L[36] += D[81] * M[0] + D[83] * M[1] + D[84] * M[2] + D[86] * M[3] + D[87] * M[4] +
           D[88] * M[5];
#pragma omp atomic
  L[37] += D[82] * M[0] + D[84] * M[1] + D[85] * M[2] + D[87] * M[3] + D[88] * M[4] +
           D[89] * M[5];
#pragma omp atomic
  L[38] += D[83] * M[0] + D[86] * M[1] + D[87] * M[2] + D[90] * M[3] + D[91] * M[4] +
           D[92] * M[5];
#pragma omp atomic
  L[39] += D[84] * M[0] + D[87] * M[1] + D[88] * M[2] + D[91] * M[3] + D[92] * M[4] +
           D[93] * M[5];
#pragma omp atomic
  L[40] += D[85] * M[0] + D[88] * M[1] + D[89] * M[2] + D[92] * M[3] + D[93] * M[4] +
           D[94] * M[5];
#pragma omp atomic
  L[41] += D[86] * M[0] + D[90] * M[1] + D[91] * M[2] + D[95] * M[3] + D[96] * M[4] +
           D[97] * M[5];
#pragma omp atomic
  L[42] += D[87] * M[0] + D[91] * M[1] + D[92] * M[2] + D[96] * M[3] + D[97] * M[4] +
           D[98] * M[5];
#pragma omp atomic
  L[43] += D[88] * M[0] + D[92] * M[1] + D[93] * M[2] + D[97] * M[3] + D[98] * M[4] +
           D[99] * M[5];
#pragma omp atomic
  L[44] += D[89] * M[0] + D[93] * M[1] + D[94] * M[2] + D[98] * M[3] + D[99] * M[4] +
           D[100] * M[5];
#pragma omp atomic
  L[45] += D[90] * M[0] + D[95] * M[1] + D[96] * M[2] + D[101] * M[3] + D[102] * M[4] +
           D[103] * M[5];
#pragma omp atomic
  L[46] += D[91] * M[0] + D[96] * M[1] + D[97] * M[2] + D[102] * M[3] + D[103] * M[4] +
           D[104] * M[5];
#pragma omp atomic
  L[47] += D[92] * M[0] + D[97] * M[1] + D[98] * M[2] + D[103] * M[3] + D[104] * M[4] +
           D[105] * M[5];
#pragma omp atomic
  L[48] += D[93] * M[0] + D[98] * M[1] + D[99] * M[2] + D[104] * M[3] + D[105] * M[4] +
           D[106] * M[5];
#pragma omp atomic
  L[49] += D[94] * M[0] + D[99] * M[1] + D[100] * M[2] + D[105] * M[3] + D[106] * M[4] +
           D[107] * M[5];
#pragma omp atomic
  L[50] += D[95] * M[0] + D[101] * M[1] + D[102] * M[2] + D[108] * M[3] + D[109] * M[4] +
           D[110] * M[5];
#pragma omp atomic
  L[51] += D[96] * M[0] + D[102] * M[1] + D[103] * M[2] + D[109] * M[3] + D[110] * M[4] +
           D[111] * M[5];
#pragma omp atomic
  L[52] += D[97] * M[0] + D[103] * M[1] + D[104] * M[2] + D[110] * M[3] + D[111] * M[4] +
           D[112] * M[5];
#pragma omp atomic
  L[53] += D[98] * M[0] + D[104] * M[1] + D[105] * M[2] + D[111] * M[3] + D[112] * M[4] +
           D[113] * M[5];
#pragma omp atomic
  L[54] += D[99] * M[0] + D[105] * M[1] + D[106] * M[2] + D[112] * M[3] + D[113] * M[4] +
           D[114] * M[5];
#pragma omp atomic
  L[55] += D[100] * M[0] + D[106] * M[1] + D[107] * M[2] + D[113] * M[3] + D[114] * M[4] +
           D[115] * M[5];
}

void field_m2_L2L_7(double x, double y, double z, double* L, double* Ls) {
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

void field_m2_L2P_7(double x, double y, double z, double* L, double* F) {
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

void field_m2_M2P_7(double x, double y, double z, double* M, double* F) {
  double R       = sqrt(x * x + y * y + z * z);
  double Ftmp0   = pow(R, -5.0);
  double Ftmp1   = 3.0 * Ftmp0;
  double Ftmp2   = Ftmp1 * M[1];
  double Ftmp3   = Ftmp1 * z;
  double Ftmp4   = y * z;
  double Ftmp5   = pow(R, -7.0);
  double Ftmp6   = 15.0 * Ftmp5;
  double Ftmp7   = Ftmp6 * M[10];
  double Ftmp8   = Ftmp4 * x;
  double Ftmp9   = Ftmp6 * Ftmp8;
  double Ftmp10  = (x * x);
  double Ftmp11  = Ftmp10 * Ftmp6;
  double Ftmp12  = pow(R, -9.0);
  double Ftmp13  = 105.0 * Ftmp12;
  double Ftmp14  = Ftmp10 * Ftmp13;
  double Ftmp15  = -9.0 * Ftmp0;
  double Ftmp16  = Ftmp11 + Ftmp15;
  double Ftmp17  = -Ftmp1;
  double Ftmp18  = (y * y);
  double Ftmp19  = Ftmp18 * Ftmp6;
  double Ftmp20  = Ftmp17 + Ftmp19;
  double Ftmp21  = (z * z);
  double Ftmp22  = Ftmp21 * Ftmp6;
  double Ftmp23  = Ftmp17 + Ftmp22;
  double Ftmp24  = Ftmp20 * M[3];
  double Ftmp25  = Ftmp23 * M[5];
  double Ftmp26  = 45.0 * Ftmp5;
  double Ftmp27  = -Ftmp26;
  double Ftmp28  = Ftmp14 + Ftmp27;
  double Ftmp29  = Ftmp28 * M[17];
  double Ftmp30  = 1.0 * y;
  double Ftmp31  = Ftmp13 * Ftmp18;
  double Ftmp32  = Ftmp27 + Ftmp31;
  double Ftmp33  = Ftmp32 * M[22];
  double Ftmp34  = 3.0 * y;
  double Ftmp35  = 35.0 * Ftmp12;
  double Ftmp36  = (Ftmp21 * Ftmp35 - 5.0 * Ftmp5) * M[24];
  double Ftmp37  = Ftmp28 * M[18];
  double Ftmp38  = 1.0 * z;
  double Ftmp39  = -Ftmp6;
  double Ftmp40  = Ftmp31 + Ftmp39;
  double Ftmp41  = Ftmp40 * M[23];
  double Ftmp42  = Ftmp13 * Ftmp21;
  double Ftmp43  = Ftmp27 + Ftmp42;
  double Ftmp44  = Ftmp38 * Ftmp43;
  double Ftmp45  = 315.0 * Ftmp12;
  double Ftmp46  = -Ftmp45;
  double Ftmp47  = pow(R, -11.0);
  double Ftmp48  = 945.0 * Ftmp47;
  double Ftmp49  = Ftmp10 * Ftmp48;
  double Ftmp50  = Ftmp46 + Ftmp49;
  double Ftmp51  = Ftmp50 * M[35];
  double Ftmp52  = x * y;
  double Ftmp53  = Ftmp28 * Ftmp52;
  double Ftmp54  = Ftmp32 * M[12];
  double Ftmp55  = Ftmp39 + Ftmp42;
  double Ftmp56  = Ftmp55 * M[14];
  double Ftmp57  = Ftmp30 * x;
  double Ftmp58  = x * z;
  double Ftmp59  = Ftmp28 * Ftmp58;
  double Ftmp60  = Ftmp40 * M[13];
  double Ftmp61  = Ftmp43 * M[15];
  double Ftmp62  = Ftmp30 * z;
  double Ftmp63  = Ftmp18 * Ftmp48;
  double Ftmp64  = Ftmp46 + Ftmp63;
  double Ftmp65  = Ftmp64 * M[42];
  double Ftmp66  = Ftmp34 * z;
  double Ftmp67  = -Ftmp13;
  double Ftmp68  = Ftmp21 * Ftmp47;
  double Ftmp69  = 315.0 * Ftmp68;
  double Ftmp70  = Ftmp67 + Ftmp69;
  double Ftmp71  = Ftmp70 * M[44];
  double Ftmp72  = Ftmp64 * M[27];
  double Ftmp73  = -75.0 * Ftmp5;
  double Ftmp74  = 1.0 * Ftmp10;
  double Ftmp75  = Ftmp40 * M[9];
  double Ftmp76  = Ftmp55 * M[11];
  double Ftmp77  = 525.0 * Ftmp12;
  double Ftmp78  = -Ftmp77;
  double Ftmp79  = Ftmp49 + Ftmp78;
  double Ftmp80  = Ftmp10 * y;
  double Ftmp81  = Ftmp10 * z;
  double Ftmp82  = Ftmp21 * Ftmp48;
  double Ftmp83  = Ftmp46 + Ftmp82;
  double Ftmp84  = Ftmp30 * M[29];
  double Ftmp85  = Ftmp10 * Ftmp30;
  double Ftmp86  = Ftmp64 * M[22];
  double Ftmp87  = Ftmp10 * Ftmp34;
  double Ftmp88  = (-Ftmp35 + Ftmp69) * M[24];
  double Ftmp89  = Ftmp10 * Ftmp38;
  double Ftmp90  = (Ftmp63 + Ftmp67) * M[23];
  double Ftmp91  = Ftmp83 * M[25];
  double Ftmp92  = 4725.0 * Ftmp47;
  double Ftmp93  = -Ftmp92;
  double Ftmp94  = pow(R, -13.0);
  double Ftmp95  = 10395.0 * Ftmp94;
  double Ftmp96  = Ftmp10 * Ftmp95;
  double Ftmp97  = Ftmp93 + Ftmp96;
  double Ftmp98  = Ftmp10 * Ftmp4;
  double Ftmp99  = Ftmp30 * Ftmp81;
  double Ftmp100 = 2835.0 * Ftmp47;
  double Ftmp101 = -Ftmp100;
  double Ftmp102 = Ftmp18 * Ftmp95;
  double Ftmp103 = (Ftmp101 + Ftmp102) * M[42];
  double Ftmp104 = Ftmp34 * Ftmp81;
  double Ftmp105 = 3465.0 * Ftmp94;
  double Ftmp106 = Ftmp105 * Ftmp21;
  double Ftmp107 = (Ftmp106 - Ftmp48) * M[44];
  double Ftmp108 = 225.0 * Ftmp5;
  double Ftmp109 = (x * x * x * x);
  double Ftmp110 = Ftmp109 * Ftmp48;
  double Ftmp111 = 1050.0 * Ftmp12;
  double Ftmp112 = -Ftmp10 * Ftmp111 + Ftmp108 + Ftmp110;
  double Ftmp113 = (y * y * y * y);
  double Ftmp114 = Ftmp113 * Ftmp48;
  double Ftmp115 = 630.0 * Ftmp12;
  double Ftmp116 = Ftmp114 - Ftmp115 * Ftmp18 + Ftmp26;
  double Ftmp117 = (z * z * z * z);
  double Ftmp118 = Ftmp117 * Ftmp48;
  double Ftmp119 = -Ftmp115 * Ftmp21 + Ftmp118 + Ftmp26;
  double Ftmp120 = Ftmp116 * M[26];
  double Ftmp121 = Ftmp119 * M[30];
  double Ftmp122 = 1575.0 * Ftmp12;
  double Ftmp123 = Ftmp109 * Ftmp94;
  double Ftmp124 = 10395.0 * Ftmp123;
  double Ftmp125 = Ftmp10 * Ftmp47;
  double Ftmp126 = 9450.0 * Ftmp125;
  double Ftmp127 = Ftmp122 + Ftmp124 - Ftmp126;
  double Ftmp128 = Ftmp127 * M[53];
  double Ftmp129 = Ftmp113 * Ftmp95;
  double Ftmp130 = 9450.0 * Ftmp47;
  double Ftmp131 = Ftmp130 * Ftmp18;
  double Ftmp132 = Ftmp122 + Ftmp129 - Ftmp131;
  double Ftmp133 = Ftmp132 * M[67];
  double Ftmp134 = (Ftmp105 * Ftmp117 + Ftmp13 - 1890.0 * Ftmp68) * M[71];
  double Ftmp135 = Ftmp127 * M[54];
  double Ftmp136 = 5670.0 * Ftmp47;
  double Ftmp137 = Ftmp136 * Ftmp18;
  double Ftmp138 = Ftmp129 - Ftmp137 + Ftmp45;
  double Ftmp139 = Ftmp138 * M[68];
  double Ftmp140 = Ftmp117 * Ftmp95;
  double Ftmp141 = Ftmp130 * Ftmp21;
  double Ftmp142 = Ftmp122 + Ftmp140 - Ftmp141;
  double Ftmp143 = Ftmp142 * Ftmp38;
  double Ftmp144 = 14175.0 * Ftmp47;
  double Ftmp145 = pow(R, -15.0);
  double Ftmp146 = 135135.0 * Ftmp145;
  double Ftmp147 = Ftmp109 * Ftmp146;
  double Ftmp148 = Ftmp10 * Ftmp94;
  double Ftmp149 = 103950.0 * Ftmp148;
  double Ftmp150 = Ftmp144 + Ftmp147 - Ftmp149;
  double Ftmp151 = Ftmp150 * M[84];
  double Ftmp152 = Ftmp127 * Ftmp52;
  double Ftmp153 = Ftmp132 * M[46];
  double Ftmp154 = Ftmp136 * Ftmp21;
  double Ftmp155 = Ftmp140 - Ftmp154 + Ftmp45;
  double Ftmp156 = Ftmp155 * M[50];
  double Ftmp157 = Ftmp127 * Ftmp58;
  double Ftmp158 = Ftmp138 * M[47];
  double Ftmp159 = Ftmp142 * M[51];
  double Ftmp160 = Ftmp113 * Ftmp146;
  double Ftmp161 = Ftmp18 * Ftmp94;
  double Ftmp162 = 103950.0 * Ftmp161;
  double Ftmp163 = Ftmp144 + Ftmp160 - Ftmp162;
  double Ftmp164 = Ftmp163 * M[102];
  double Ftmp165 = Ftmp117 * Ftmp145;
  double Ftmp166 = 45045.0 * Ftmp165;
  double Ftmp167 = Ftmp21 * Ftmp94;
  double Ftmp168 = Ftmp166 - 34650.0 * Ftmp167 + Ftmp92;
  double Ftmp169 = Ftmp168 * M[106];
  double Ftmp170 = Ftmp163 * M[74];
  double Ftmp171 = 3675.0 * Ftmp12;
  double Ftmp172 = Ftmp138 * M[41];
  double Ftmp173 = Ftmp155 * M[45];
  double Ftmp174 = 33075.0 * Ftmp47;
  double Ftmp175 = 145530.0 * Ftmp148;
  double Ftmp176 = Ftmp147 + Ftmp174 - Ftmp175;
  double Ftmp177 = Ftmp117 * Ftmp146;
  double Ftmp178 = 103950.0 * Ftmp167;
  double Ftmp179 = Ftmp144 + Ftmp177 - Ftmp178;
  double Ftmp180 = Ftmp30 * M[78];
  double Ftmp181 = Ftmp163 * M[67];
  double Ftmp182 = (Ftmp166 - 20790.0 * Ftmp167 + Ftmp48) * M[71];
  double Ftmp183 = 62370.0 * Ftmp94;
  double Ftmp184 = Ftmp18 * Ftmp183;
  double Ftmp185 = (Ftmp100 + Ftmp160 - Ftmp184) * M[68];
  double Ftmp186 = Ftmp179 * M[72];
  double Ftmp187 = 363825.0 * Ftmp94;
  double Ftmp188 = pow(R, -17.0);
  double Ftmp189 = 2027025.0 * Ftmp188;
  double Ftmp190 = Ftmp109 * Ftmp189;
  double Ftmp191 = Ftmp10 * Ftmp145;
  double Ftmp192 = 1891890.0 * Ftmp191;
  double Ftmp193 = 155925.0 * Ftmp94;
  double Ftmp194 = Ftmp113 * Ftmp189;
  double Ftmp195 = Ftmp145 * Ftmp18;
  double Ftmp196 = 1351350.0 * Ftmp195;
  double Ftmp197 = (Ftmp193 + Ftmp194 - Ftmp196) * M[102];
  double Ftmp198 = 51975.0 * Ftmp94;
  double Ftmp199 = 675675.0 * Ftmp117 * Ftmp188;
  double Ftmp200 = Ftmp145 * Ftmp21;
  double Ftmp201 = (Ftmp198 + Ftmp199 - 450450.0 * Ftmp200) * M[106];
  double Ftmp202 = -11025.0 * Ftmp12;
  double Ftmp203 = (x * x * x * x * x * x);
  double Ftmp204 = Ftmp146 * Ftmp203;
  double Ftmp205 = 99225.0 * Ftmp47;
  double Ftmp206 = Ftmp10 * Ftmp205 - 218295.0 * Ftmp123 + Ftmp202 + Ftmp204;
  double Ftmp207 = -Ftmp122;
  double Ftmp208 = (y * y * y * y * y * y);
  double Ftmp209 = Ftmp146 * Ftmp208;
  double Ftmp210 = 42525.0 * Ftmp47;
  double Ftmp211 = -Ftmp113 * Ftmp193 + Ftmp18 * Ftmp210 + Ftmp207 + Ftmp209;
  double Ftmp212 = (z * z * z * z * z * z);
  double Ftmp213 = Ftmp146 * Ftmp212;
  double Ftmp214 = -Ftmp117 * Ftmp193 + Ftmp207 + Ftmp21 * Ftmp210 + Ftmp213;
  double Ftmp215 = Ftmp211 * M[73];
  double Ftmp216 = Ftmp214 * M[79];
  double Ftmp217 = -Ftmp18 * Ftmp45;
  double Ftmp218 = Ftmp18 * Ftmp49;
  double Ftmp219 = -Ftmp14;
  double Ftmp220 = Ftmp219 + Ftmp26;
  double Ftmp221 = Ftmp217 + Ftmp218 + Ftmp220;
  double Ftmp222 = -Ftmp21 * Ftmp45;
  double Ftmp223 = Ftmp21 * Ftmp49;
  double Ftmp224 = Ftmp220 + Ftmp222 + Ftmp223;
  double Ftmp225 = -Ftmp42;
  double Ftmp226 = Ftmp225 + Ftmp6;
  double Ftmp227 = -Ftmp31;
  double Ftmp228 = Ftmp21 * Ftmp63;
  double Ftmp229 = Ftmp227 + Ftmp228;
  double Ftmp230 = Ftmp226 + Ftmp229;
  double Ftmp231 = -Ftmp205;
  double Ftmp232 = Ftmp189 * Ftmp203;
  double Ftmp233 = Ftmp109 * Ftmp145;
  double Ftmp234 = 1091475.0 * Ftmp148 + Ftmp231 + Ftmp232 - 2837835.0 * Ftmp233;
  double Ftmp235 = Ftmp234 * Ftmp52;
  double Ftmp236 = Ftmp189 * Ftmp208;
  double Ftmp237 = Ftmp113 * Ftmp145;
  double Ftmp238 = 1091475.0 * Ftmp161 + Ftmp231 + Ftmp236 - 2837835.0 * Ftmp237;
  double Ftmp239 = Ftmp238 * M[108];
  double Ftmp240 = -Ftmp144;
  double Ftmp241 = Ftmp189 * Ftmp212;
  double Ftmp242 = 2027025.0 * Ftmp145;
  double Ftmp243 = 467775.0 * Ftmp94;
  double Ftmp244 = -Ftmp117 * Ftmp242 + Ftmp21 * Ftmp243 + Ftmp240 + Ftmp241;
  double Ftmp245 = Ftmp244 * M[114];
  double Ftmp246 = Ftmp234 * Ftmp58;
  double Ftmp247 = -Ftmp113 * Ftmp242 + Ftmp18 * Ftmp243 + Ftmp236 + Ftmp240;
  double Ftmp248 = Ftmp247 * M[109];
  double Ftmp249 = -2837835.0 * Ftmp165 + 1091475.0 * Ftmp167 + Ftmp231 + Ftmp241;
  double Ftmp250 = Ftmp249 * M[115];
  double Ftmp251 = -297675.0 * Ftmp47;
  double Ftmp252 = Ftmp247 * M[101];
  double Ftmp253 = Ftmp244 * M[107];
  double Ftmp254 = Ftmp100 * Ftmp18;
  double Ftmp255 = -Ftmp254;
  double Ftmp256 = Ftmp18 * Ftmp96;
  double Ftmp257 = Ftmp255 + Ftmp256;
  double Ftmp258 = 945.0 * Ftmp12;
  double Ftmp259 = Ftmp10 * Ftmp100;
  double Ftmp260 = -Ftmp259;
  double Ftmp261 = Ftmp258 + Ftmp260;
  double Ftmp262 = Ftmp257 + Ftmp261;
  double Ftmp263 = Ftmp262 * M[58];
  double Ftmp264 = -Ftmp49;
  double Ftmp265 = Ftmp264 + Ftmp45;
  double Ftmp266 = Ftmp100 * Ftmp21;
  double Ftmp267 = -Ftmp266;
  double Ftmp268 = Ftmp21 * Ftmp96;
  double Ftmp269 = Ftmp267 + Ftmp268;
  double Ftmp270 = Ftmp265 + Ftmp269;
  double Ftmp271 = Ftmp270 * M[60];
  double Ftmp272 = -Ftmp63;
  double Ftmp273 = Ftmp102 * Ftmp21;
  double Ftmp274 = Ftmp273 + Ftmp45;
  double Ftmp275 = Ftmp267 + Ftmp272 + Ftmp274;
  double Ftmp276 = Ftmp275 * M[69];
  double Ftmp277 = Ftmp257 + Ftmp265;
  double Ftmp278 = Ftmp277 * M[59];
  double Ftmp279 = Ftmp261 + Ftmp269;
  double Ftmp280 = Ftmp279 * M[61];
  double Ftmp281 = -Ftmp82;
  double Ftmp282 = Ftmp255 + Ftmp274 + Ftmp281;
  double Ftmp283 = Ftmp282 * M[70];
  double Ftmp284 = 8505.0 * Ftmp47;
  double Ftmp285 = 31185.0 * Ftmp94;
  double Ftmp286 = Ftmp10 * Ftmp285;
  double Ftmp287 = -Ftmp286;
  double Ftmp288 = Ftmp284 + Ftmp287;
  double Ftmp289 = Ftmp18 * Ftmp285;
  double Ftmp290 = -Ftmp289;
  double Ftmp291 = Ftmp10 * Ftmp146;
  double Ftmp292 = Ftmp18 * Ftmp291;
  double Ftmp293 = Ftmp290 + Ftmp292;
  double Ftmp294 = Ftmp288 + Ftmp293;
  double Ftmp295 = Ftmp294 * M[91];
  double Ftmp296 = Ftmp21 * Ftmp285;
  double Ftmp297 = -Ftmp296;
  double Ftmp298 = Ftmp21 * Ftmp291;
  double Ftmp299 = Ftmp297 + Ftmp298;
  double Ftmp300 = Ftmp288 + Ftmp299;
  double Ftmp301 = Ftmp300 * M[93];
  double Ftmp302 = Ftmp262 * Ftmp52;
  double Ftmp303 = Ftmp270 * Ftmp52;
  double Ftmp304 = Ftmp277 * Ftmp58;
  double Ftmp305 = Ftmp279 * Ftmp58;
  double Ftmp306 = Ftmp18 * Ftmp21;
  double Ftmp307 = Ftmp146 * Ftmp306;
  double Ftmp308 = Ftmp297 + Ftmp307;
  double Ftmp309 = Ftmp284 + Ftmp290 + Ftmp308;
  double Ftmp310 = Ftmp309 * M[104];
  double Ftmp311 = Ftmp294 * Ftmp8;
  double Ftmp312 = Ftmp300 * Ftmp8;
  double Ftmp313 = -Ftmp18 * Ftmp92;
  double Ftmp314 = Ftmp264 + Ftmp77;
  double Ftmp315 = -Ftmp21 * Ftmp92;
  double Ftmp316 = Ftmp13 + Ftmp281;
  double Ftmp317 = Ftmp272 + Ftmp273;
  double Ftmp318 = Ftmp18 * Ftmp198;
  double Ftmp319 = -Ftmp318;
  double Ftmp320 = Ftmp292 + Ftmp319;
  double Ftmp321 = Ftmp144 + Ftmp287;
  double Ftmp322 = -Ftmp96;
  double Ftmp323 = Ftmp322 + Ftmp92;
  double Ftmp324 = Ftmp198 * Ftmp21;
  double Ftmp325 = -Ftmp324;
  double Ftmp326 = Ftmp298 + Ftmp325;
  double Ftmp327 = -Ftmp102;
  double Ftmp328 = Ftmp100 + Ftmp327;
  double Ftmp329 = Ftmp21 * Ftmp95;
  double Ftmp330 = -Ftmp329;
  double Ftmp331 = Ftmp100 + Ftmp330;
  double Ftmp332 = Ftmp290 + Ftmp307;
  double Ftmp333 = 675675.0 * Ftmp145;
  double Ftmp334 = Ftmp18 * Ftmp333;
  double Ftmp335 = -Ftmp334;
  double Ftmp336 = Ftmp10 * Ftmp189;
  double Ftmp337 = Ftmp18 * Ftmp336;
  double Ftmp338 = 405405.0 * Ftmp191;
  double Ftmp339 = -Ftmp338;
  double Ftmp340 = Ftmp193 + Ftmp339;
  double Ftmp341 = Ftmp21 * Ftmp333;
  double Ftmp342 = -Ftmp341;
  double Ftmp343 = Ftmp21 * Ftmp336;
  double Ftmp344 = 93555.0 * Ftmp94;
  double Ftmp345 = -405405.0 * Ftmp200;
  double Ftmp346 = Ftmp344 + Ftmp345;
  double Ftmp347 = 405405.0 * Ftmp195;
  double Ftmp348 = -Ftmp347;
  double Ftmp349 = Ftmp189 * Ftmp306;
  double Ftmp350 = Ftmp348 + Ftmp349;
  double Ftmp351 = -Ftmp258;
  double Ftmp352 = Ftmp259 + Ftmp351;
  double Ftmp353 = Ftmp10 * Ftmp160;
  double Ftmp354 = 62370.0 * Ftmp148;
  double Ftmp355 = -Ftmp18 * Ftmp354;
  double Ftmp356 = Ftmp353 + Ftmp355;
  double Ftmp357 = 17010.0 * Ftmp47;
  double Ftmp358 = -Ftmp113 * Ftmp285 + Ftmp18 * Ftmp357;
  double Ftmp359 = Ftmp352 + Ftmp356 + Ftmp358;
  double Ftmp360 = Ftmp10 * Ftmp177;
  double Ftmp361 = -Ftmp21 * Ftmp354;
  double Ftmp362 = Ftmp360 + Ftmp361;
  double Ftmp363 = -Ftmp117 * Ftmp285 + Ftmp21 * Ftmp357;
  double Ftmp364 = Ftmp352 + Ftmp362 + Ftmp363;
  double Ftmp365 = Ftmp144 * Ftmp18;
  double Ftmp366 = -Ftmp149 * Ftmp18 + Ftmp207;
  double Ftmp367 = -Ftmp124;
  double Ftmp368 = Ftmp147 * Ftmp18;
  double Ftmp369 = Ftmp367 + Ftmp368;
  double Ftmp370 = Ftmp126 + Ftmp365 + Ftmp366 + Ftmp369;
  double Ftmp371 = -Ftmp149 * Ftmp21;
  double Ftmp372 = Ftmp147 * Ftmp21;
  double Ftmp373 = Ftmp367 + Ftmp372;
  double Ftmp374 = Ftmp144 * Ftmp21 + Ftmp207;
  double Ftmp375 = Ftmp126 + Ftmp371 + Ftmp373 + Ftmp374;
  double Ftmp376 = Ftmp183 * Ftmp21;
  double Ftmp377 = -Ftmp18 * Ftmp376;
  double Ftmp378 = Ftmp377 + Ftmp46;
  double Ftmp379 = -Ftmp140;
  double Ftmp380 = Ftmp177 * Ftmp18;
  double Ftmp381 = Ftmp379 + Ftmp380;
  double Ftmp382 = Ftmp154 + Ftmp254 + Ftmp378 + Ftmp381;
  double Ftmp383 = -Ftmp129;
  double Ftmp384 = Ftmp160 * Ftmp21;
  double Ftmp385 = Ftmp383 + Ftmp384;
  double Ftmp386 = Ftmp137 + Ftmp266 + Ftmp378 + Ftmp385;
  double Ftmp387 = Ftmp10 * Ftmp193;
  double Ftmp388 = -405405.0 * Ftmp237;
  double Ftmp389 = 311850.0 * Ftmp94;
  double Ftmp390 = Ftmp18 * Ftmp389;
  double Ftmp391 = Ftmp10 * Ftmp194;
  double Ftmp392 = Ftmp390 + Ftmp391;
  double Ftmp393 = -Ftmp210;
  double Ftmp394 = 1351350.0 * Ftmp191;
  double Ftmp395 = -Ftmp18 * Ftmp394;
  double Ftmp396 = Ftmp393 + Ftmp395;
  double Ftmp397 = Ftmp52 * (Ftmp387 + Ftmp388 + Ftmp392 + Ftmp396);
  double Ftmp398 = Ftmp117 * Ftmp189;
  double Ftmp399 = Ftmp10 * Ftmp398;
  double Ftmp400 = 810810.0 * Ftmp191;
  double Ftmp401 = -Ftmp21 * Ftmp400;
  double Ftmp402 = Ftmp399 + Ftmp401;
  double Ftmp403 = -Ftmp284;
  double Ftmp404 = Ftmp286 + Ftmp403;
  double Ftmp405 = -405405.0 * Ftmp165;
  double Ftmp406 = 187110.0 * Ftmp167 + Ftmp405;
  double Ftmp407 = Ftmp52 * (Ftmp402 + Ftmp404 + Ftmp406);
  double Ftmp408 = Ftmp18 * Ftmp190;
  double Ftmp409 = Ftmp18 * Ftmp193;
  double Ftmp410 = 311850.0 * Ftmp148;
  double Ftmp411 = -405405.0 * Ftmp233;
  double Ftmp412 = Ftmp410 + Ftmp411;
  double Ftmp413 = Ftmp52 * (Ftmp396 + Ftmp408 + Ftmp409 + Ftmp412);
  double Ftmp414 = -Ftmp21 * Ftmp394;
  double Ftmp415 = -Ftmp147;
  double Ftmp416 = Ftmp190 * Ftmp21;
  double Ftmp417 = Ftmp415 + Ftmp416;
  double Ftmp418 = Ftmp193 * Ftmp21;
  double Ftmp419 = Ftmp240 + Ftmp418;
  double Ftmp420 = Ftmp52 * (Ftmp149 + Ftmp414 + Ftmp417 + Ftmp419);
  double Ftmp421 = -810810.0 * Ftmp195 * Ftmp21;
  double Ftmp422 = Ftmp403 + Ftmp421;
  double Ftmp423 = Ftmp18 * Ftmp398;
  double Ftmp424 = Ftmp289 + Ftmp423;
  double Ftmp425 = Ftmp406 + Ftmp422 + Ftmp424;
  double Ftmp426 = -Ftmp196 * Ftmp21;
  double Ftmp427 = -Ftmp160;
  double Ftmp428 = Ftmp194 * Ftmp21;
  double Ftmp429 = Ftmp427 + Ftmp428;
  double Ftmp430 = Ftmp162 + Ftmp419 + Ftmp426 + Ftmp429;
  double Ftmp431 = 187110.0 * Ftmp161 + Ftmp388;
  double Ftmp432 = -Ftmp18 * Ftmp400;
  double Ftmp433 = Ftmp391 + Ftmp432;
  double Ftmp434 = Ftmp58 * (Ftmp404 + Ftmp431 + Ftmp433);
  double Ftmp435 = Ftmp393 + Ftmp414;
  double Ftmp436 = Ftmp387 + Ftmp399;
  double Ftmp437 = Ftmp21 * Ftmp389;
  double Ftmp438 = Ftmp405 + Ftmp437;
  double Ftmp439 = Ftmp58 * (Ftmp435 + Ftmp436 + Ftmp438);
  double Ftmp440 = Ftmp408 + Ftmp415;
  double Ftmp441 = Ftmp240 + Ftmp409;
  double Ftmp442 = Ftmp58 * (Ftmp149 + Ftmp395 + Ftmp440 + Ftmp441);
  double Ftmp443 = Ftmp58 * (Ftmp412 + Ftmp416 + Ftmp418 + Ftmp435);
  double Ftmp444 = -Ftmp177;
  double Ftmp445 = Ftmp423 + Ftmp444;
  double Ftmp446 = Ftmp178 + Ftmp426 + Ftmp441 + Ftmp445;
  double Ftmp447 = Ftmp296 + Ftmp428;
  double Ftmp448 = Ftmp422 + Ftmp431 + Ftmp447;
  double Ftmp449 = -Ftmp113 * Ftmp333;
  double Ftmp450 = Ftmp240 + Ftmp286;
  double Ftmp451 = -Ftmp117 * Ftmp333 + Ftmp437;
  double Ftmp452 = Ftmp18 * Ftmp187;
  double Ftmp453 = -Ftmp174;
  double Ftmp454 = -Ftmp18 * Ftmp192 + Ftmp453;
  double Ftmp455 = -Ftmp192 * Ftmp21;
  double Ftmp456 = Ftmp187 * Ftmp21 + Ftmp453;
  double Ftmp457 = Ftmp101 + Ftmp421;
  double Ftmp458 = Ftmp21 * Ftmp292;
  double Ftmp459 = -Ftmp256 + Ftmp266 + Ftmp458;
  double Ftmp460 = Ftmp254 - Ftmp268;
  double Ftmp461 = -Ftmp21 * Ftmp289 + Ftmp459 + Ftmp460 + Ftmp50;
  double Ftmp462 = Ftmp306 * Ftmp336;
  double Ftmp463 = -Ftmp292 + Ftmp462;
  double Ftmp464 = -Ftmp21 * Ftmp347 + Ftmp404;
  double Ftmp465 = -Ftmp21 * Ftmp338 + Ftmp289;
  double Ftmp466 = Ftmp52 * (Ftmp21 * Ftmp344 + Ftmp463 + Ftmp464 + Ftmp465);
  double Ftmp467 = -Ftmp298;
  double Ftmp468 = -Ftmp18 * Ftmp338 + Ftmp296 + Ftmp462;
  double Ftmp469 = Ftmp58 * (Ftmp18 * Ftmp344 + Ftmp464 + Ftmp467 + Ftmp468);
  double Ftmp470 = Ftmp324 + Ftmp463;
  double Ftmp471 = Ftmp318 + Ftmp467;
  double Ftmp472 = x * M[2];
  double Ftmp473 = Ftmp11 + Ftmp17;
  double Ftmp474 = Ftmp15 + Ftmp19;
  double Ftmp475 = Ftmp473 * M[0];
  double Ftmp476 = 1.0 * x;
  double Ftmp477 = 3.0 * x;
  double Ftmp478 = Ftmp14 + Ftmp39;
  double Ftmp479 = Ftmp478 * M[20];
  double Ftmp480 = Ftmp32 * M[27];
  double Ftmp481 = Ftmp38 * x;
  double Ftmp482 = 3.0 * Ftmp58;
  double Ftmp483 = Ftmp478 * M[8];
  double Ftmp484 = Ftmp50 * M[18];
  double Ftmp485 = Ftmp478 * M[7];
  double Ftmp486 = 1.0 * Ftmp18;
  double Ftmp487 = Ftmp18 * x;
  double Ftmp488 = Ftmp50 * M[17];
  double Ftmp489 = Ftmp18 * z;
  double Ftmp490 = (Ftmp49 + Ftmp67) * M[20];
  double Ftmp491 = Ftmp63 + Ftmp78;
  double Ftmp492 = Ftmp30 * Ftmp58;
  double Ftmp493 = Ftmp18 * Ftmp476;
  double Ftmp494 = Ftmp18 * Ftmp477;
  double Ftmp495 = Ftmp18 * Ftmp38;
  double Ftmp496 = Ftmp18 * Ftmp58;
  double Ftmp497 = (Ftmp101 + Ftmp96) * M[35];
  double Ftmp498 = Ftmp102 + Ftmp93;
  double Ftmp499 = Ftmp38 * Ftmp487;
  double Ftmp500 = Ftmp18 * Ftmp482;
  double Ftmp501 = -Ftmp10 * Ftmp115 + Ftmp110 + Ftmp26;
  double Ftmp502 = Ftmp108 - Ftmp111 * Ftmp18 + Ftmp114;
  double Ftmp503 = Ftmp501 * M[16];
  double Ftmp504 = Ftmp10 * Ftmp136;
  double Ftmp505 = Ftmp124 + Ftmp45 - Ftmp504;
  double Ftmp506 = Ftmp505 * M[56];
  double Ftmp507 = Ftmp132 * M[74];
  double Ftmp508 = Ftmp505 * M[33];
  double Ftmp509 = Ftmp150 * M[54];
  double Ftmp510 = Ftmp505 * M[32];
  double Ftmp511 = Ftmp150 * M[53];
  double Ftmp512 = (Ftmp100 + Ftmp147 - Ftmp354) * M[56];
  double Ftmp513 = 145530.0 * Ftmp161;
  double Ftmp514 = Ftmp160 + Ftmp174 - Ftmp513;
  double Ftmp515 = (Ftmp190 + Ftmp193 - Ftmp394) * M[84];
  double Ftmp516 = 1891890.0 * Ftmp195;
  double Ftmp517 = Ftmp10 * Ftmp210 - Ftmp109 * Ftmp193 + Ftmp204 + Ftmp207;
  double Ftmp518 = 218295.0 * Ftmp94;
  double Ftmp519 = -Ftmp113 * Ftmp518 + Ftmp18 * Ftmp205 + Ftmp202 + Ftmp209;
  double Ftmp520 = Ftmp517 * M[52];
  double Ftmp521 = Ftmp218 + Ftmp227;
  double Ftmp522 = -Ftmp10 * Ftmp45 + Ftmp26;
  double Ftmp523 = Ftmp521 + Ftmp522;
  double Ftmp524 = Ftmp219 + Ftmp223 + Ftmp226;
  double Ftmp525 = Ftmp222 + Ftmp229 + Ftmp26;
  double Ftmp526 = 467775.0 * Ftmp148 + Ftmp232 - 2027025.0 * Ftmp233 + Ftmp240;
  double Ftmp527 = Ftmp526 * M[82];
  double Ftmp528 = Ftmp526 * M[81];
  double Ftmp529 = Ftmp260 + Ftmp45;
  double Ftmp530 = Ftmp256 + Ftmp272;
  double Ftmp531 = Ftmp529 + Ftmp530;
  double Ftmp532 = Ftmp531 * M[63];
  double Ftmp533 = Ftmp268 + Ftmp281;
  double Ftmp534 = Ftmp529 + Ftmp533;
  double Ftmp535 = Ftmp534 * M[65];
  double Ftmp536 = Ftmp255 + Ftmp258 + Ftmp267 + Ftmp273;
  double Ftmp537 = Ftmp536 * M[76];
  double Ftmp538 = Ftmp4 * Ftmp531;
  double Ftmp539 = Ftmp4 * Ftmp534;
  double Ftmp540 = Ftmp4 * Ftmp536;
  double Ftmp541 = -Ftmp10 * Ftmp92 + Ftmp77;
  double Ftmp542 = Ftmp10 * Ftmp198;
  double Ftmp543 = -Ftmp542;
  double Ftmp544 = Ftmp144 + Ftmp543;
  double Ftmp545 = Ftmp100 + Ftmp322;
  double Ftmp546 = Ftmp327 + Ftmp92;
  double Ftmp547 = Ftmp309 * Ftmp492;
  double Ftmp548 = Ftmp337 + Ftmp348;
  double Ftmp549 = -Ftmp10 * Ftmp333 + Ftmp193;
  double Ftmp550 = Ftmp10 * Ftmp144;
  double Ftmp551 = Ftmp131 + Ftmp353 + Ftmp366 + Ftmp383 + Ftmp550;
  double Ftmp552 = Ftmp259 + Ftmp46;
  double Ftmp553 = Ftmp154 + Ftmp362 + Ftmp379 + Ftmp552;
  double Ftmp554 = Ftmp254 + Ftmp351;
  double Ftmp555 = Ftmp10 * Ftmp357 - 31185.0 * Ftmp123;
  double Ftmp556 = Ftmp355 + Ftmp368 + Ftmp554 + Ftmp555;
  double Ftmp557 = Ftmp46 + Ftmp504;
  double Ftmp558 = Ftmp266 + Ftmp361;
  double Ftmp559 = Ftmp373 + Ftmp557 + Ftmp558;
  double Ftmp560 = Ftmp363 + Ftmp377 + Ftmp380 + Ftmp554;
  double Ftmp561 = -Ftmp162 * Ftmp21;
  double Ftmp562 = Ftmp131 + Ftmp374 + Ftmp385 + Ftmp561;
  double Ftmp563 = Ftmp391 + Ftmp427;
  double Ftmp564 = Ftmp4 * (Ftmp162 + Ftmp240 + Ftmp387 + Ftmp395 + Ftmp563);
  double Ftmp565 = Ftmp4 * (Ftmp178 + Ftmp240 + Ftmp414 + Ftmp436 + Ftmp444);
  double Ftmp566 = Ftmp289 + Ftmp432;
  double Ftmp567 = Ftmp408 + Ftmp566;
  double Ftmp568 = 187110.0 * Ftmp148 + Ftmp403 + Ftmp411;
  double Ftmp569 = Ftmp4 * (Ftmp567 + Ftmp568);
  double Ftmp570 = Ftmp296 + Ftmp401;
  double Ftmp571 = Ftmp416 + Ftmp570;
  double Ftmp572 = Ftmp4 * (Ftmp568 + Ftmp571);
  double Ftmp573 = Ftmp393 + Ftmp426;
  double Ftmp574 = Ftmp4 * (Ftmp409 + Ftmp423 + Ftmp438 + Ftmp573);
  double Ftmp575 = Ftmp4 * (Ftmp388 + Ftmp390 + Ftmp418 + Ftmp428 + Ftmp573);
  double Ftmp576 = Ftmp10 * Ftmp187;
  double Ftmp577 = Ftmp101 + Ftmp286;
  double Ftmp578 = -675675.0 * Ftmp233 + Ftmp240 + Ftmp410;
  double Ftmp579 = Ftmp101 + Ftmp354;
  double Ftmp580 = Ftmp240 + Ftmp421;
  double Ftmp581 = -Ftmp21 * Ftmp516;
  double Ftmp582 = Ftmp259 - Ftmp273;
  double Ftmp583 = -Ftmp21 * Ftmp286 + Ftmp459 + Ftmp582 + Ftmp64;
  double Ftmp584 = -Ftmp307;
  double Ftmp585 = Ftmp4 * (Ftmp10 * Ftmp344 + Ftmp403 + Ftmp465 + Ftmp468 + Ftmp584);
  double Ftmp586 = Ftmp542 + Ftmp584;
  double Ftmp587 = y * M[4];
  double Ftmp588 = Ftmp15 + Ftmp22;
  double Ftmp589 = Ftmp476 * M[25];
  double Ftmp590 = Ftmp34 * x;
  double Ftmp591 = Ftmp21 * x;
  double Ftmp592 = Ftmp21 * y;
  double Ftmp593 = Ftmp34 * Ftmp58;
  double Ftmp594 = Ftmp21 * Ftmp476;
  double Ftmp595 = Ftmp21 * (Ftmp78 + Ftmp82);
  double Ftmp596 = Ftmp21 * Ftmp52;
  double Ftmp597 = Ftmp30 * Ftmp591;
  double Ftmp598 = Ftmp34 * Ftmp591;
  double Ftmp599 = Ftmp108 - Ftmp111 * Ftmp21 + Ftmp118;
  double Ftmp600 = Ftmp476 * M[72];
  double Ftmp601 = 145530.0 * Ftmp167;
  double Ftmp602 = Ftmp21 * (Ftmp174 + Ftmp177 - Ftmp601);
  double Ftmp603 = -Ftmp117 * Ftmp518 + Ftmp202 + Ftmp205 * Ftmp21 + Ftmp213;
  double Ftmp604 = Ftmp219 + Ftmp521 + Ftmp6;
  double Ftmp605 = Ftmp223 + Ftmp225 + Ftmp522;
  double Ftmp606 = Ftmp217 + Ftmp225 + Ftmp228 + Ftmp26;
  double Ftmp607 = Ftmp330 + Ftmp92;
  double Ftmp608 = Ftmp137 + Ftmp356 + Ftmp383 + Ftmp552;
  double Ftmp609 = Ftmp141 + Ftmp207;
  double Ftmp610 = Ftmp360 + Ftmp371 + Ftmp379 + Ftmp550 + Ftmp609;
  double Ftmp611 = Ftmp254 + Ftmp355 + Ftmp369 + Ftmp557;
  double Ftmp612 = Ftmp351 + Ftmp372 + Ftmp555 + Ftmp558;
  double Ftmp613 = Ftmp365 + Ftmp381 + Ftmp561 + Ftmp609;
  double Ftmp614 = Ftmp266 + Ftmp351 + Ftmp358 + Ftmp377 + Ftmp384;
  double Ftmp615 = Ftmp453 + Ftmp601;
  double Ftmp616 = -Ftmp18 * Ftmp286 + Ftmp458 + Ftmp460 + Ftmp582 + Ftmp83;
#pragma omp atomic
  F[0] +=
        -Ftmp10 * (Ftmp14 + Ftmp73) * M[6] -
        Ftmp10 * (Ftmp124 - 13230.0 * Ftmp125 + Ftmp171) * M[31] -
        Ftmp10 * (Ftmp256 + Ftmp313 + Ftmp314) * M[34] -
        Ftmp10 * (Ftmp268 + Ftmp314 + Ftmp315) * M[36] -
        Ftmp10 * (Ftmp402 + Ftmp450 + Ftmp451) * M[94] -
        Ftmp10 * (1964655.0 * Ftmp148 + Ftmp232 - 3648645.0 * Ftmp233 + Ftmp251) * M[80] -
        Ftmp10 * (Ftmp175 + Ftmp417 + Ftmp455 + Ftmp456) * M[85] -
        Ftmp10 * (Ftmp175 + Ftmp440 + Ftmp452 + Ftmp454) * M[83] -
        Ftmp10 * (Ftmp392 + Ftmp432 + Ftmp449 + Ftmp450) * M[90] -
        Ftmp10 * (-Ftmp21 * Ftmp334 + Ftmp470 + Ftmp471 + Ftmp97) * M[92] -
        Ftmp103 * Ftmp99 - Ftmp104 * Ftmp107 - Ftmp104 * Ftmp201 + Ftmp11 * y * M[1] +
        Ftmp11 * z * M[2] + Ftmp112 * x * M[16] + Ftmp112 * M[31] + Ftmp116 * M[41] +
        Ftmp119 * M[45] + Ftmp120 * x + Ftmp121 * x - Ftmp128 * y - Ftmp133 * Ftmp30 -
        Ftmp134 * Ftmp34 - Ftmp135 * z - Ftmp139 * Ftmp38 - Ftmp14 * Ftmp4 * M[10] -
        Ftmp143 * M[72] + Ftmp150 * Ftmp8 * M[56] + Ftmp151 * Ftmp4 - Ftmp152 * M[32] -
        Ftmp153 * Ftmp52 - Ftmp156 * Ftmp57 - Ftmp157 * M[33] - Ftmp158 * Ftmp58 -
        Ftmp159 * Ftmp58 + Ftmp16 * x * M[0] + Ftmp16 * M[6] + Ftmp164 * Ftmp62 +
        Ftmp169 * Ftmp66 + Ftmp170 * Ftmp8 - Ftmp172 * Ftmp74 - Ftmp173 * Ftmp74 +
        Ftmp176 * Ftmp80 * M[53] + Ftmp176 * Ftmp81 * M[54] + Ftmp179 * Ftmp180 * Ftmp58 +
        Ftmp181 * Ftmp85 + Ftmp182 * Ftmp87 + Ftmp185 * Ftmp89 + Ftmp186 * Ftmp89 -
        Ftmp197 * Ftmp99 - Ftmp2 * y + Ftmp20 * M[9] + Ftmp206 * x * M[52] +
        Ftmp206 * M[80] + Ftmp211 * M[101] + Ftmp214 * M[107] + Ftmp215 * x +
        Ftmp216 * x + Ftmp221 * x * M[19] + Ftmp221 * M[34] + Ftmp224 * x * M[21] +
        Ftmp224 * M[36] + Ftmp23 * M[11] + Ftmp230 * x * M[28] + Ftmp230 * M[43] -
        Ftmp235 * M[81] - Ftmp239 * Ftmp52 + Ftmp24 * x - Ftmp245 * Ftmp57 -
        Ftmp246 * M[82] - Ftmp248 * Ftmp58 + Ftmp25 * x - Ftmp250 * Ftmp58 -
        Ftmp252 * Ftmp74 - Ftmp253 * Ftmp74 - Ftmp263 * y - Ftmp271 * y -
        Ftmp275 * Ftmp52 * M[48] - Ftmp276 * Ftmp30 - Ftmp278 * z - Ftmp280 * z -
        Ftmp282 * Ftmp58 * M[49] - Ftmp283 * Ftmp38 - Ftmp29 * y + Ftmp295 * Ftmp4 -
        Ftmp3 * M[2] - Ftmp30 * Ftmp33 + Ftmp301 * Ftmp4 - Ftmp302 * M[37] -
        Ftmp303 * M[39] - Ftmp304 * M[38] - Ftmp305 * M[40] + Ftmp309 * Ftmp8 * M[76] +
        Ftmp310 * Ftmp62 + Ftmp311 * M[63] + Ftmp312 * M[65] - Ftmp34 * Ftmp36 +
        Ftmp359 * x * M[62] + Ftmp359 * M[90] + Ftmp364 * x * M[66] + Ftmp364 * M[94] -
        Ftmp37 * z + Ftmp370 * x * M[55] + Ftmp370 * M[83] + Ftmp375 * x * M[57] +
        Ftmp375 * M[85] - Ftmp38 * Ftmp41 + Ftmp382 * x * M[77] + Ftmp382 * M[105] +
        Ftmp386 * x * M[75] + Ftmp386 * M[103] - Ftmp397 * M[95] + Ftmp4 * Ftmp51 +
        Ftmp4 * Ftmp7 - Ftmp407 * M[99] - Ftmp413 * M[86] - Ftmp420 * M[88] -
        Ftmp425 * Ftmp52 * M[112] - Ftmp430 * Ftmp52 * M[110] - Ftmp434 * M[96] -
        Ftmp439 * M[100] - Ftmp44 * M[25] - Ftmp442 * M[87] - Ftmp443 * M[89] -
        Ftmp446 * Ftmp58 * M[113] - Ftmp448 * Ftmp58 * M[111] + Ftmp461 * x * M[64] +
        Ftmp461 * M[92] - Ftmp466 * M[97] - Ftmp469 * M[98] + Ftmp50 * Ftmp8 * M[20] -
        Ftmp52 * Ftmp54 - Ftmp53 * M[7] - Ftmp56 * Ftmp57 - Ftmp58 * Ftmp60 -
        Ftmp58 * Ftmp61 + Ftmp58 * Ftmp83 * Ftmp84 - Ftmp59 * M[8] + Ftmp62 * Ftmp65 +
        Ftmp66 * Ftmp71 + Ftmp72 * Ftmp8 - Ftmp74 * Ftmp75 - Ftmp74 * Ftmp76 -
        Ftmp74 * (Ftmp316 + Ftmp317) * M[43] -
        Ftmp74 * (Ftmp184 + Ftmp296 + Ftmp429 + Ftmp457) * M[103] -
        Ftmp74 * (Ftmp289 + Ftmp376 + Ftmp445 + Ftmp457) * M[105] +
        Ftmp79 * Ftmp80 * M[17] + Ftmp79 * Ftmp81 * M[18] +
        Ftmp80 * (Ftmp320 + Ftmp321) * M[58] + Ftmp80 * (Ftmp323 + Ftmp326) * M[60] +
        Ftmp81 * (Ftmp320 + Ftmp323) * M[59] + Ftmp81 * (Ftmp321 + Ftmp326) * M[61] +
        Ftmp85 * Ftmp86 + Ftmp85 * (Ftmp308 + Ftmp328) * M[69] + Ftmp87 * Ftmp88 +
        Ftmp89 * Ftmp90 + Ftmp89 * Ftmp91 + Ftmp89 * (Ftmp331 + Ftmp332) * M[70] +
        Ftmp9 * M[4] - Ftmp97 * Ftmp98 * M[35] -
        Ftmp98 * (Ftmp187 + Ftmp190 - Ftmp192) * M[84] -
        Ftmp98 * (Ftmp335 + Ftmp337 + Ftmp340) * M[91] -
        Ftmp98 * (Ftmp340 + Ftmp342 + Ftmp343) * M[93] -
        Ftmp99 * (Ftmp346 + Ftmp350) * M[104];
#pragma omp atomic
  F[1] += -Ftmp107 * Ftmp500 + Ftmp119 * M[50] + Ftmp121 * y - Ftmp128 * x -
          Ftmp132 * Ftmp4 * M[47] - Ftmp132 * Ftmp57 * M[41] - Ftmp133 * Ftmp476 -
          Ftmp134 * Ftmp477 - Ftmp143 * M[78] + Ftmp151 * Ftmp58 - Ftmp152 * M[31] -
          Ftmp156 * Ftmp486 - Ftmp159 * Ftmp4 + Ftmp163 * Ftmp492 * M[68] +
          Ftmp164 * Ftmp481 + Ftmp169 * Ftmp482 - Ftmp173 * Ftmp57 +
          Ftmp179 * Ftmp495 * M[78] - Ftmp18 * Ftmp485 - Ftmp18 * Ftmp510 -
          Ftmp18 * Ftmp528 - Ftmp18 * (Ftmp31 + Ftmp73) * M[12] -
          Ftmp18 * (Ftmp530 + Ftmp541) * M[37] - Ftmp18 * (Ftmp567 + Ftmp578) * M[86] -
          Ftmp18 * (Ftmp129 + Ftmp171 - 13230.0 * Ftmp18 * Ftmp47) * M[46] -
          Ftmp18 * (Ftmp264 + Ftmp268 + Ftmp316) * M[39] -
          Ftmp18 * (Ftmp315 + Ftmp317 + Ftmp77) * M[48] -
          Ftmp18 * (Ftmp417 + Ftmp570 + Ftmp579) * M[88] -
          Ftmp18 * (Ftmp424 + Ftmp451 + Ftmp580) * M[112] -
          Ftmp18 * (1964655.0 * Ftmp161 + Ftmp236 - 3648645.0 * Ftmp237 + Ftmp251) *
                M[108] -
          Ftmp18 * (Ftmp376 + Ftmp402 + Ftmp444 + Ftmp577) * M[99] -
          Ftmp18 * (Ftmp429 + Ftmp456 + Ftmp513 + Ftmp581) * M[110] -
          Ftmp18 * (Ftmp454 + Ftmp513 + Ftmp563 + Ftmp576) * M[95] -
          Ftmp18 * (-Ftmp10 * Ftmp341 + Ftmp470 + Ftmp498 + Ftmp586) * M[97] +
          Ftmp182 * Ftmp494 + Ftmp186 * Ftmp492 + Ftmp19 * x * M[1] + Ftmp19 * z * M[4] -
          Ftmp2 * x - Ftmp201 * Ftmp500 + Ftmp214 * M[114] + Ftmp216 * y +
          Ftmp23 * M[14] - Ftmp235 * M[80] - Ftmp238 * Ftmp4 * M[109] -
          Ftmp238 * Ftmp57 * M[101] - Ftmp245 * Ftmp486 + Ftmp25 * y - Ftmp250 * Ftmp4 -
          Ftmp253 * Ftmp57 - Ftmp263 * x - Ftmp271 * x - Ftmp275 * Ftmp57 * M[43] -
          Ftmp276 * Ftmp476 - Ftmp29 * x + Ftmp295 * Ftmp58 - Ftmp3 * M[4] +
          Ftmp301 * Ftmp58 - Ftmp302 * M[34] - Ftmp303 * M[36] - Ftmp31 * Ftmp58 * M[10] +
          Ftmp310 * Ftmp481 + Ftmp311 * M[59] + Ftmp312 * M[61] - Ftmp32 * Ftmp4 * M[13] -
          Ftmp32 * Ftmp57 * M[9] - Ftmp33 * Ftmp476 - Ftmp36 * Ftmp477 - Ftmp397 * M[90] +
          Ftmp4 * Ftmp472 * Ftmp6 - Ftmp4 * Ftmp483 - Ftmp4 * Ftmp508 - Ftmp4 * Ftmp527 -
          Ftmp4 * Ftmp61 - Ftmp407 * M[94] - Ftmp413 * M[83] - Ftmp420 * M[85] -
          Ftmp425 * Ftmp57 * M[105] - Ftmp430 * Ftmp57 * M[103] - Ftmp44 * M[29] -
          Ftmp466 * M[92] + Ftmp473 * M[7] + Ftmp474 * y * M[3] + Ftmp474 * M[12] +
          Ftmp475 * y - Ftmp479 * z - Ftmp480 * z + Ftmp481 * Ftmp65 + Ftmp482 * Ftmp71 +
          Ftmp484 * Ftmp8 - Ftmp486 * Ftmp56 + Ftmp487 * Ftmp488 + Ftmp487 * Ftmp511 +
          Ftmp487 * (Ftmp293 + Ftmp544) * M[58] + Ftmp487 * (Ftmp299 + Ftmp545) * M[60] +
          Ftmp489 * Ftmp490 + Ftmp489 * Ftmp491 * M[27] + Ftmp489 * Ftmp512 +
          Ftmp489 * Ftmp514 * M[74] + Ftmp489 * (Ftmp144 + Ftmp325 + Ftmp332) * M[76] +
          Ftmp489 * (Ftmp287 + Ftmp298 + Ftmp331) * M[65] +
          Ftmp489 * (Ftmp292 + Ftmp543 + Ftmp546) * M[63] + Ftmp491 * Ftmp493 * M[22] +
          Ftmp492 * Ftmp64 * M[23] + Ftmp492 * Ftmp91 + Ftmp493 * Ftmp514 * M[67] +
          Ftmp493 * (Ftmp307 + Ftmp325 + Ftmp546) * M[69] + Ftmp494 * Ftmp88 +
          Ftmp495 * Ftmp83 * M[29] - Ftmp496 * Ftmp497 - Ftmp496 * Ftmp515 -
          Ftmp496 * (Ftmp548 + Ftmp549) * M[91] -
          Ftmp496 * (Ftmp339 + Ftmp343 + Ftmp346) * M[93] - Ftmp498 * Ftmp499 * M[42] -
          Ftmp499 * (Ftmp187 + Ftmp194 - Ftmp516) * M[102] -
          Ftmp499 * (Ftmp193 + Ftmp342 + Ftmp350) * M[104] + Ftmp501 * M[32] +
          Ftmp502 * y * M[26] + Ftmp502 * M[46] + Ftmp503 * y - Ftmp506 * z -
          Ftmp507 * z + Ftmp509 * Ftmp8 + Ftmp51 * Ftmp58 + Ftmp517 * M[81] +
          Ftmp519 * y * M[73] + Ftmp519 * M[108] + Ftmp520 * y + Ftmp523 * y * M[19] +
          Ftmp523 * M[37] + Ftmp524 * y * M[21] + Ftmp524 * M[39] + Ftmp525 * y * M[28] +
          Ftmp525 * M[48] - Ftmp53 * M[6] - Ftmp532 * z - Ftmp535 * z - Ftmp537 * z -
          Ftmp538 * M[38] - Ftmp539 * M[40] - Ftmp540 * M[49] + Ftmp547 * M[70] +
          Ftmp551 * y * M[62] + Ftmp551 * M[95] + Ftmp553 * y * M[66] + Ftmp553 * M[99] +
          Ftmp556 * y * M[55] + Ftmp556 * M[86] + Ftmp559 * y * M[57] + Ftmp559 * M[88] +
          Ftmp560 * y * M[77] + Ftmp560 * M[112] + Ftmp562 * y * M[75] +
          Ftmp562 * M[110] - Ftmp564 * M[96] - Ftmp565 * M[100] - Ftmp569 * M[87] -
          Ftmp57 * Ftmp76 - Ftmp572 * M[89] - Ftmp574 * M[113] - Ftmp575 * M[111] +
          Ftmp58 * Ftmp7 + Ftmp583 * y * M[64] + Ftmp583 * M[97] - Ftmp585 * M[98];
#pragma omp atomic
  F[2] += -Ftmp1 * Ftmp472 - Ftmp1 * Ftmp587 - Ftmp103 * Ftmp597 + Ftmp116 * M[47] +
          Ftmp120 * z - Ftmp135 * x - Ftmp139 * Ftmp476 - Ftmp142 * Ftmp180 -
          Ftmp142 * Ftmp600 - Ftmp142 * Ftmp62 * M[50] - Ftmp143 * x * M[45] +
          Ftmp151 * Ftmp52 - Ftmp153 * Ftmp4 - Ftmp157 * M[31] - Ftmp158 * Ftmp21 +
          Ftmp164 * Ftmp57 + Ftmp168 * Ftmp593 * M[71] + Ftmp169 * Ftmp590 +
          Ftmp170 * Ftmp592 - Ftmp172 * Ftmp481 + Ftmp180 * Ftmp602 + Ftmp181 * Ftmp492 +
          Ftmp185 * Ftmp594 - Ftmp197 * Ftmp597 + Ftmp20 * M[13] - Ftmp21 * Ftmp248 -
          Ftmp21 * Ftmp483 - Ftmp21 * Ftmp508 - Ftmp21 * Ftmp527 - Ftmp21 * Ftmp60 -
          Ftmp21 * (Ftmp42 + Ftmp73) * M[15] - Ftmp21 * (Ftmp533 + Ftmp541) * M[40] -
          Ftmp21 * (Ftmp571 + Ftmp578) * M[89] -
          Ftmp21 * (Ftmp13 + Ftmp264 + Ftmp530) * M[38] -
          Ftmp21 * (Ftmp140 + Ftmp171 - 13230.0 * Ftmp68) * M[51] -
          Ftmp21 * (Ftmp440 + Ftmp566 + Ftmp579) * M[87] -
          Ftmp21 * (-3648645.0 * Ftmp165 + 1964655.0 * Ftmp167 + Ftmp241 + Ftmp251) *
                M[115] -
          Ftmp21 * (Ftmp184 + Ftmp427 + Ftmp433 + Ftmp577) * M[96] -
          Ftmp21 * (Ftmp273 + Ftmp281 + Ftmp313 + Ftmp77) * M[49] -
          Ftmp21 * (Ftmp390 + Ftmp447 + Ftmp449 + Ftmp580) * M[111] -
          Ftmp21 * (Ftmp445 + Ftmp452 + Ftmp581 + Ftmp615) * M[113] -
          Ftmp21 * (Ftmp399 + Ftmp444 + Ftmp455 + Ftmp576 + Ftmp615) * M[100] -
          Ftmp21 * (-Ftmp10 * Ftmp334 + Ftmp329 + Ftmp462 + Ftmp471 + Ftmp586 + Ftmp93) *
                M[98] +
          Ftmp211 * M[109] + Ftmp215 * z + Ftmp22 * Ftmp472 + Ftmp22 * Ftmp587 -
          Ftmp239 * Ftmp4 + Ftmp24 * z - Ftmp246 * M[80] - Ftmp249 * Ftmp481 * M[107] -
          Ftmp249 * Ftmp62 * M[114] - Ftmp252 * Ftmp481 - Ftmp278 * x - Ftmp280 * x -
          Ftmp282 * Ftmp481 * M[43] - Ftmp283 * Ftmp476 + Ftmp295 * Ftmp52 +
          Ftmp301 * Ftmp52 - Ftmp304 * M[34] - Ftmp305 * M[36] + Ftmp310 * Ftmp57 +
          Ftmp311 * M[58] + Ftmp312 * M[60] - Ftmp37 * x - Ftmp4 * Ftmp485 -
          Ftmp4 * Ftmp510 - Ftmp4 * Ftmp528 - Ftmp4 * Ftmp54 - Ftmp41 * Ftmp476 -
          Ftmp42 * Ftmp52 * M[10] - Ftmp43 * Ftmp589 - Ftmp43 * Ftmp62 * M[14] -
          Ftmp43 * Ftmp84 - Ftmp434 * M[90] - Ftmp439 * M[94] - Ftmp44 * x * M[11] -
          Ftmp442 * M[83] - Ftmp443 * M[85] - Ftmp446 * Ftmp481 * M[105] -
          Ftmp448 * Ftmp481 * M[103] - Ftmp469 * M[92] + Ftmp473 * M[8] + Ftmp475 * z -
          Ftmp479 * y - Ftmp480 * y - Ftmp481 * Ftmp75 + Ftmp484 * Ftmp591 +
          Ftmp488 * Ftmp8 + Ftmp490 * Ftmp592 + Ftmp492 * Ftmp86 - Ftmp497 * Ftmp596 +
          Ftmp501 * M[33] + Ftmp503 * z - Ftmp506 * y - Ftmp507 * y + Ftmp509 * Ftmp591 +
          Ftmp51 * Ftmp52 + Ftmp511 * Ftmp8 + Ftmp512 * Ftmp592 - Ftmp515 * Ftmp596 +
          Ftmp517 * M[82] + Ftmp52 * Ftmp7 + Ftmp520 * z - Ftmp532 * y - Ftmp535 * y -
          Ftmp537 * y - Ftmp538 * M[37] - Ftmp539 * M[39] - Ftmp540 * M[48] +
          Ftmp547 * M[69] - Ftmp564 * M[95] - Ftmp565 * M[99] - Ftmp569 * M[86] +
          Ftmp57 * Ftmp65 - Ftmp572 * M[88] - Ftmp574 * M[112] - Ftmp575 * M[110] -
          Ftmp585 * M[97] + Ftmp588 * z * M[5] + Ftmp588 * M[15] + Ftmp589 * Ftmp595 -
          Ftmp59 * M[6] + Ftmp590 * Ftmp71 + Ftmp591 * (Ftmp293 + Ftmp545) * M[59] +
          Ftmp591 * (Ftmp299 + Ftmp544) * M[61] + Ftmp592 * Ftmp72 +
          Ftmp592 * (Ftmp144 + Ftmp308 + Ftmp319) * M[76] +
          Ftmp592 * (Ftmp287 + Ftmp292 + Ftmp328) * M[63] +
          Ftmp592 * (Ftmp298 + Ftmp543 + Ftmp607) * M[65] + Ftmp593 * Ftmp70 * M[24] +
          Ftmp594 * Ftmp90 + Ftmp594 * (Ftmp307 + Ftmp319 + Ftmp607) * M[70] +
          Ftmp595 * Ftmp84 - Ftmp596 * (Ftmp339 + Ftmp344 + Ftmp548) * M[91] -
          Ftmp596 * (Ftmp343 + Ftmp345 + Ftmp549) * M[93] -
          Ftmp597 * (Ftmp193 + Ftmp335 + Ftmp345 + Ftmp349) * M[104] -
          Ftmp598 * (Ftmp106 - 1575.0 * Ftmp47) * M[44] -
          Ftmp598 * (Ftmp199 - 630630.0 * Ftmp200 + 121275.0 * Ftmp94) * M[106] +
          Ftmp599 * z * M[30] + Ftmp599 * M[51] + Ftmp600 * Ftmp602 +
          Ftmp603 * z * M[79] + Ftmp603 * M[115] + Ftmp604 * z * M[19] + Ftmp604 * M[38] +
          Ftmp605 * z * M[21] + Ftmp605 * M[40] + Ftmp606 * z * M[28] + Ftmp606 * M[49] +
          Ftmp608 * z * M[62] + Ftmp608 * M[96] + Ftmp610 * z * M[66] + Ftmp610 * M[100] +
          Ftmp611 * z * M[55] + Ftmp611 * M[87] + Ftmp612 * z * M[57] + Ftmp612 * M[89] +
          Ftmp613 * z * M[77] + Ftmp613 * M[113] + Ftmp614 * z * M[75] +
          Ftmp614 * M[111] + Ftmp616 * z * M[64] + Ftmp616 * M[98] + Ftmp9 * M[1];
}

template <>
void P2M<2, 3>(double x, double y, double z, double q, double* M, int order) {
  switch (order) {
    case 3:
      field_m2_P2M_3(x, y, z, q, M);
      break;
    case 4:
      field_m2_P2M_4(x, y, z, q, M);
      break;
    case 5:
      field_m2_P2M_5(x, y, z, q, M);
      break;
    case 6:
      field_m2_P2M_6(x, y, z, q, M);
      break;
    case 7:
      field_m2_P2M_7(x, y, z, q, M);
      break;
  }
}
template <>
void M2M<2, 3>(double x, double y, double z, double* M, double* Ms, int order) {
  switch (order) {
    case 3:
      field_m2_M2M_3(x, y, z, M, Ms);
      break;
    case 4:
      field_m2_M2M_4(x, y, z, M, Ms);
      break;
    case 5:
      field_m2_M2M_5(x, y, z, M, Ms);
      break;
    case 6:
      field_m2_M2M_6(x, y, z, M, Ms);
      break;
    case 7:
      field_m2_M2M_7(x, y, z, M, Ms);
      break;
  }
}
template <>
void M2L<2, 3>(double x, double y, double z, double* M, double* L, int order) {
  switch (order) {
    case 3:
      field_m2_M2L_3(x, y, z, M, L);
      break;
    case 4:
      field_m2_M2L_4(x, y, z, M, L);
      break;
    case 5:
      field_m2_M2L_5(x, y, z, M, L);
      break;
    case 6:
      field_m2_M2L_6(x, y, z, M, L);
      break;
    case 7:
      field_m2_M2L_7(x, y, z, M, L);
      break;
  }
}
template <>
void L2L<2, 3>(double x, double y, double z, double* L, double* Ls, int order) {
  switch (order) {
    case 3:
      field_m2_L2L_3(x, y, z, L, Ls);
      break;
    case 4:
      field_m2_L2L_4(x, y, z, L, Ls);
      break;
    case 5:
      field_m2_L2L_5(x, y, z, L, Ls);
      break;
    case 6:
      field_m2_L2L_6(x, y, z, L, Ls);
      break;
    case 7:
      field_m2_L2L_7(x, y, z, L, Ls);
      break;
  }
}
template <>
void L2P<2, 3>(double x, double y, double z, double* L, double* F, int order) {
  switch (order) {
    case 3:
      field_m2_L2P_3(x, y, z, L, F);
      break;
    case 4:
      field_m2_L2P_4(x, y, z, L, F);
      break;
    case 5:
      field_m2_L2P_5(x, y, z, L, F);
      break;
    case 6:
      field_m2_L2P_6(x, y, z, L, F);
      break;
    case 7:
      field_m2_L2P_7(x, y, z, L, F);
      break;
  }
}
template <>
void M2P<2, 3>(double x, double y, double z, double* M, double* F, int order) {
  switch (order) {
    case 3:
      field_m2_M2P_3(x, y, z, M, F);
      break;
    case 4:
      field_m2_M2P_4(x, y, z, M, F);
      break;
    case 5:
      field_m2_M2P_5(x, y, z, M, F);
      break;
    case 6:
      field_m2_M2P_6(x, y, z, M, F);
      break;
    case 7:
      field_m2_M2P_7(x, y, z, M, F);
      break;
  }
}

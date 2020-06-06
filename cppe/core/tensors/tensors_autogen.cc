// maximum order, maximum order for Thole-damped tensors, maximum expansion order of pow
// maxorder = 6, maxorder_damped = 6, maxpow = 5
#include "tensors_autogen.hh"
namespace libcppe {
namespace tensors {
Eigen::VectorXd T0(const Eigen::Vector3d& rij) {
  double x = rij(0);
  double y = rij(1);
  double z = rij(2);
  Eigen::VectorXd result(1);
  result[0] = pow(x * x + y * y + z * z, -0.5);  //
  return result;
}
Eigen::VectorXd T1(const Eigen::Vector3d& rij) {
  double x = rij(0);
  double y = rij(1);
  double z = rij(2);
  Eigen::VectorXd result(3);
  double x0 = pow(x * x + y * y + z * z, -1.5);
  result[0] = -x * x0;  // x
  result[1] = -x0 * y;  // y
  result[2] = -x0 * z;  // z
  return result;
}
Eigen::VectorXd T2(const Eigen::Vector3d& rij) {
  double x = rij(0);
  double y = rij(1);
  double z = rij(2);
  Eigen::VectorXd result(6);
  double x0 = x * x;
  double x1 = y * y;
  double x2 = z * z;
  double x3 = x0 + x1 + x2;
  double x4 = pow(x3, -1.5);
  double x5 = 3.0 * 1.0 / x3;
  double x6 = pow(x3, -2.5);
  double x7 = 3.0 * x * x6;
  result[0] = x4 * (x0 * x5 - 1.0);  // xx
  result[1] = x7 * y;                // xy
  result[2] = x7 * z;                // xz
  result[3] = x4 * (x1 * x5 - 1.0);  // yy
  result[4] = 3.0 * x6 * y * z;      // yz
  result[5] = x4 * (x2 * x5 - 1.0);  // zz
  return result;
}
Eigen::VectorXd T3(const Eigen::Vector3d& rij) {
  double x = rij(0);
  double y = rij(1);
  double z = rij(2);
  Eigen::VectorXd result(10);
  double x0  = x * x;
  double x1  = y * y;
  double x2  = z * z;
  double x3  = x0 + x1 + x2;
  double x4  = 5.0 * 1.0 / x3;
  double x5  = x0 * x4;
  double x6  = 3.0 * pow(x3, -2.5);
  double x7  = x * x6;
  double x8  = x6 * (x5 - 1.0);
  double x9  = x1 * x4;
  double x10 = x9 - 1.0;
  double x11 = x2 * x4;
  double x12 = x11 - 1.0;
  double x13 = x6 * y;
  double x14 = x6 * z;
  result[0]  = -x7 * (x5 - 3.0);                   // xxx
  result[1]  = -x8 * y;                            // xxy
  result[2]  = -x8 * z;                            // xxz
  result[3]  = -x10 * x7;                          // xyy
  result[4]  = -15.0 * x * pow(x3, -3.5) * y * z;  // xyz
  result[5]  = -x12 * x7;                          // xzz
  result[6]  = -x13 * (x9 - 3.0);                  // yyy
  result[7]  = -x10 * x14;                         // yyz
  result[8]  = -x12 * x13;                         // yzz
  result[9]  = -x14 * (x11 - 3.0);                 // zzz
  return result;
}
Eigen::VectorXd T4(const Eigen::Vector3d& rij) {
  double x = rij(0);
  double y = rij(1);
  double z = rij(2);
  Eigen::VectorXd result(15);
  double x0  = x * x;
  double x1  = y * y;
  double x2  = z * z;
  double x3  = x0 + x1 + x2;
  double x4  = 1.0 / x3;
  double x5  = x0 * x4;
  double x6  = 35.0 * 1.0 / (x3 * x3);
  double x7  = 3.0 * pow(x3, -2.5);
  double x8  = pow(x3, -3.5);
  double x9  = 15.0 * x8 * y;
  double x10 = 7.0 * x5;
  double x11 = x * (x10 - 3.0);
  double x12 = 15.0 * x8 * z;
  double x13 = 5.0 * x4;
  double x14 = -x1 * x13;
  double x15 = x0 * x6;
  double x16 = 1.0 - 5.0 * x5;
  double x17 = x9 * z;
  double x18 = -x13 * x2;
  double x19 = 7.0 * x4;
  double x20 = x1 * x19;
  double x21 = x20 - 3.0;
  double x22 = x * x9;
  double x23 = x * x12;
  double x24 = x19 * x2;
  double x25 = x24 - 3.0;
  double x26 = 30.0 * x4;
  result[0]  = x7 * (-30.0 * x5 + x6 * (x * x * x * x) + 3.0);  // xxxx
  result[1]  = x11 * x9;                                        // xxxy
  result[2]  = x11 * x12;                                       // xxxz
  result[3]  = x7 * (x1 * x15 + x14 + x16);                     // xxyy
  result[4]  = x17 * (x10 - 1.0);                               // xxyz
  result[5]  = x7 * (x15 * x2 + x16 + x18);                     // xxzz
  result[6]  = x21 * x22;                                       // xyyy
  result[7]  = x23 * (x20 - 1.0);                               // xyyz
  result[8]  = x22 * (x24 - 1.0);                               // xyzz
  result[9]  = x23 * x25;                                       // xzzz
  result[10] = x7 * (-x1 * x26 + x6 * (y * y * y * y) + 3.0);   // yyyy
  result[11] = x17 * x21;                                       // yyyz
  result[12] = x7 * (x1 * x2 * x6 + x14 + x18 + 1.0);           // yyzz
  result[13] = x17 * x25;                                       // yzzz
  result[14] = x7 * (-x2 * x26 + x6 * (z * z * z * z) + 3.0);   // zzzz
  return result;
}
Eigen::VectorXd T5(const Eigen::Vector3d& rij) {
  double x = rij(0);
  double y = rij(1);
  double z = rij(2);
  Eigen::VectorXd result(21);
  double x0  = x * x;
  double x1  = y * y;
  double x2  = z * z;
  double x3  = x0 + x1 + x2;
  double x4  = 1.0 / x3;
  double x5  = x0 * x4;
  double x6  = 1.0 / (x3 * x3);
  double x7  = x6 * (x * x * x * x);
  double x8  = pow(x3, -3.5);
  double x9  = 15.0 * x8;
  double x10 = x * x9;
  double x11 = 45.0 * x8;
  double x12 = x11 * (-14.0 * x5 + 21.0 * x7 + 1.0);
  double x13 = 21.0 * x4;
  double x14 = -x1 * x13;
  double x15 = 63.0 * x6;
  double x16 = x0 * x15;
  double x17 = x1 * x16;
  double x18 = -7.0 * x5;
  double x19 = x18 + 3.0;
  double x20 = 315.0 * x * pow(x3, -4.5) * y * z;
  double x21 = -x13 * x2;
  double x22 = x16 * x2;
  double x23 = 7.0 * x4;
  double x24 = -x1 * x23;
  double x25 = x17 + x24;
  double x26 = 3.0 - 21.0 * x5;
  double x27 = x9 * y;
  double x28 = x18 + 1.0;
  double x29 = x9 * z;
  double x30 = x2 * x23;
  double x31 = -x30;
  double x32 = x22 + x31;
  double x33 = x1 * x4;
  double x34 = y * y * y * y;
  double x35 = 21.0 * x6;
  double x36 = -14.0 * x33 + x34 * x35 + 1.0;
  double x37 = x * x11;
  double x38 = x1 * x2;
  double x39 = x2 * x4;
  double x40 = z * z * z * z;
  double x41 = x35 * x40 - 14.0 * x39 + 1.0;
  double x42 = x15 * x38 + 3.0;
  result[0]  = -x10 * (-70.0 * x5 + 63.0 * x7 + 15.0);  // xxxxx
  result[1]  = -x12 * y;                                // xxxxy
  result[2]  = -x12 * z;                                // xxxxz
  result[3]  = -x10 * (x14 + x17 + x19);                // xxxyy
  result[4]  = -x20 * (3.0 * x5 - 1.0);                 // xxxyz
  result[5]  = -x10 * (x19 + x21 + x22);                // xxxzz
  result[6]  = -x27 * (x25 + x26);                      // xxyyy
  result[7]  = -x29 * (x25 + x28);                      // xxyyz
  result[8]  = -x27 * (x28 + x32);                      // xxyzz
  result[9]  = -x29 * (x26 + x32);                      // xxzzz
  result[10] = -x36 * x37;                              // xyyyy
  result[11] = -x20 * (3.0 * x33 - 1.0);                // xyyyz
  result[12] = -x10 * (2.0 * x33 * (4.0 * x39 - 1.0) + 20.0 * x38 * x6 +
                       (x30 - 1.0) * (5.0 * x33 - 1.0));  // xyyzz
  result[13] = -x20 * (3.0 * x39 - 1.0);                  // xyzzz
  result[14] = -x37 * x41;                                // xzzzz
  result[15] = -x27 * (x15 * x34 - 70.0 * x33 + 15.0);    // yyyyy
  result[16] = -x11 * x36 * z;                            // yyyyz
  result[17] = -x27 * (x21 + x24 + x42);                  // yyyzz
  result[18] = -x29 * (x14 + x31 + x42);                  // yyzzz
  result[19] = -x11 * x41 * y;                            // yzzzz
  result[20] = -x29 * (x15 * x40 - 70.0 * x39 + 15.0);    // zzzzz
  return result;
}
Eigen::VectorXd T6(const Eigen::Vector3d& rij) {
  double x = rij(0);
  double y = rij(1);
  double z = rij(2);
  Eigen::VectorXd result(28);
  double x0  = x * x;
  double x1  = y * y;
  double x2  = z * z;
  double x3  = x0 + x1 + x2;
  double x4  = 1.0 / x3;
  double x5  = x0 * x4;
  double x6  = x * x * x * x;
  double x7  = 1.0 / (x3 * x3);
  double x8  = x6 * x7;
  double x9  = 1.0 / (x3 * x3 * x3);
  double x10 = 231.0 * x9;
  double x11 = pow(x3, -3.5);
  double x12 = 45.0 * x11;
  double x13 = 315.0 * y;
  double x14 = 33.0 * x8;
  double x15 = pow(x3, -4.5);
  double x16 = x * x15;
  double x17 = x16 * (x14 - 30.0 * x5 + 5.0);
  double x18 = 315.0 * z;
  double x19 = 7.0 * x4;
  double x20 = x1 * x19;
  double x21 = x10 * x6;
  double x22 = 14.0 * x5 - 21.0 * x8;
  double x23 = x0 * x1;
  double x24 = 126.0 * x7;
  double x25 = -x23 * x24 - 1.0;
  double x26 = x15 * z;
  double x27 = x13 * x26;
  double x28 = x0 * x2;
  double x29 = -x24 * x28;
  double x30 = x19 * x2 - 1.0;
  double x31 = x1 * x4;
  double x32 = 3.0 * x31;
  double x33 = -x32;
  double x34 = 11.0 * x7;
  double x35 = 1.0 - 3.0 * x5;
  double x36 = 945.0 * x16;
  double x37 = 33.0 * x7;
  double x38 = x23 * x37;
  double x39 = x16 * x18;
  double x40 = x2 * x4;
  double x41 = 9.0 * x40;
  double x42 = x28 * x37;
  double x43 = x13 * x16;
  double x44 = 3.0 * x40;
  double x45 = -x44;
  double x46 = 7.0 * x5;
  double x47 = y * y * y * y;
  double x48 = x0 * x10;
  double x49 = 21.0 * x7;
  double x50 = x47 * x49;
  double x51 = 14.0 * x31;
  double x52 = -x50 + x51;
  double x53 = 1.0 - 9.0 * x5;
  double x54 = 63.0 * x7;
  double x55 = x1 * x2;
  double x56 = z * z * z * z;
  double x57 = 14.0 * x40 - x49 * x56 - 1.0;
  double x58 = -30.0 * x31 + x37 * x47 + 5.0;
  double x59 = x31 * (4.0 * x40 - 1.0);
  double x60 = x37 * x56;
  double x61 = -30.0 * x40 + x60 + 5.0;
  double x62 = 315.0 * x7;
  double x63 = -x24 * x55;
  result[0]  = x12 * (pow(x, 6) * x10 + 105.0 * x5 - 315.0 * x8 - 5.0);  // xxxxxx
  result[1]  = x13 * x17;                                                // xxxxxy
  result[2]  = x17 * x18;                                                // xxxxxz
  result[3]  = x12 * (x1 * x21 + x20 + x22 + x25);                       // xxxxyy
  result[4]  = x27 * (x14 - 18.0 * x5 + 1.0);                            // xxxxyz
  result[5]  = x12 * (x2 * x21 + x22 + x29 + x30);                       // xxxxzz
  result[6]  = x36 * y * (x23 * x34 + x33 + x35);                        // xxxyyy
  result[7]  = x39 * (-9.0 * x31 + x35 + x38);                           // xxxyyz
  result[8]  = x43 * (x35 - x41 + x42);                                  // xxxyzz
  result[9]  = x36 * z * (x28 * x34 + x35 + x45);                        // xxxzzz
  result[10] = x12 * (x25 + x46 + x47 * x48 + x52);                      // xxyyyy
  result[11] = x27 * (x33 + x38 + x53);                                  // xxyyyz
  result[12] = 15.0 * x11 *
               (693.0 * x2 * x23 * x9 + x20 - x23 * x54 - x28 * x54 + x30 + x46 -
                x54 * x55);                                        // xxyyzz
  result[13] = x27 * (x42 + x45 + x53);                            // xxyzzz
  result[14] = x12 * (x29 + x46 + x48 * x56 + x57);                // xxzzzz
  result[15] = x43 * x58;                                          // xyyyyy
  result[16] = x39 * (4.0 * x31 * (x32 - 1.0) + x50 - x51 + 1.0);  // xyyyyz
  result[17] = 105.0 * x16 * y *
               (28.0 * x55 * x7 + 2.0 * x59 + (x20 - 3.0) * (x41 - 1.0));  // xyyyzz
  result[18] = 45.0 * x16 * z *
               (10.0 * x30 * x31 + 8.0 * x31 * (2.0 * x40 - 1.0) + 10.0 * x59 +
                7.0 * (5.0 * x31 - 1.0) * (x44 - 1.0));                  // xyyzzz
  result[19] = x43 * (-18.0 * x40 + x60 + 1.0);                          // xyzzzz
  result[20] = x39 * x61;                                                // xzzzzz
  result[21] = x12 * (x10 * pow(y, 6) + 105.0 * x31 - x47 * x62 - 5.0);  // yyyyyy
  result[22] = x27 * x58;                                                // yyyyyz
  result[23] = x12 * (x10 * x2 * x47 + x30 + x52 + x63);                 // yyyyzz
  result[24] = 945.0 * x26 * y * (x33 + x34 * x55 + x45 + 1.0);          // yyyzzz
  result[25] = x12 * (x1 * x10 * x56 + x20 + x57 + x63);                 // yyzzzz
  result[26] = x27 * x61;                                                // yzzzzz
  result[27] = x12 * (x10 * pow(z, 6) + 105.0 * x40 - x56 * x62 - 5.0);  // zzzzzz
  return result;
}
Eigen::VectorXd T0_damp_thole(const Eigen::Vector3d& rij, double a) {
  double x = rij(0);
  double y = rij(1);
  double z = rij(2);
  Eigen::VectorXd result(1);
  double x0 = sqrt(x * x + y * y + z * z);
  double x1 = a * x0;
  result[0] = -(0.5 * (x1 + 2.0) * exp(-x1) - 1.0) * 1.0 / x0;  //
  return result;
}
Eigen::VectorXd T1_damp_thole(const Eigen::Vector3d& rij, double a) {
  double x = rij(0);
  double y = rij(1);
  double z = rij(2);
  Eigen::VectorXd result(3);
  double x0 = x * x + y * y + z * z;
  double x1 = a * sqrt(x0);
  double x2 = -x1;
  double x3 = exp(x2);
  double x4 = -0.5 * a * x3 * (x2 - 1.0) * 1.0 / x0 +
              0.5 * pow(x0, -1.5) * (x3 * (x1 + 2.0) - 2.0);
  result[0] = x * x4;  // x
  result[1] = x4 * y;  // y
  result[2] = x4 * z;  // z
  return result;
}
Eigen::VectorXd T2_damp_thole(const Eigen::Vector3d& rij, double a) {
  double x = rij(0);
  double y = rij(1);
  double z = rij(2);
  Eigen::VectorXd result(6);
  double x0  = x * x;
  double x1  = y * y;
  double x2  = z * z;
  double x3  = x0 + x1 + x2;
  double x4  = sqrt(x3);
  double x5  = a * x4;
  double x6  = exp(-x5);
  double x7  = a * x6;
  double x8  = x7 * (x5 + 1.0) * 1.0 / (x3 * x3);
  double x9  = 1.0 / x3;
  double x10 = 3.0 * x9;
  double x11 = pow(x3, -1.5);
  double x12 = x5 + 2.0;
  double x13 = x12 * x6 - 2.0;
  double x14 = 0.5 * x11 * x13;
  double x15 = x0 * x11;
  double x16 = a * x9;
  double x17 = 2.0 * x16;
  double x18 = x12 * x16;
  double x19 = 1.0 / x4;
  double x20 = -x12 * x19 + x19;
  double x21 = 0.5 * x19 * x7;
  double x22 = x11 * x12;
  double x23 = 1.5 * x13 * pow(x3, -2.5) + x21 * (-x11 - x17 + x18 + x22) + x8;
  double x24 = x * x23;
  result[0]  = -x0 * x8 - x14 * (x0 * x10 - 1.0) -
              x21 * (-x0 * x17 + x0 * x18 + x12 * x15 - x15 + x20);  // xx
  result[1] = -x24 * y;                                              // xy
  result[2] = -x24 * z;                                              // xz
  result[3] = -x1 * x8 - x14 * (x1 * x10 - 1.0) -
              x21 * (-x1 * x11 - x1 * x17 + x1 * x18 + x1 * x22 + x20);  // yy
  result[4] = -x23 * y * z;                                              // yz
  result[5] = -x14 * (x10 * x2 - 1.0) - x2 * x8 -
              x21 * (-x11 * x2 - x17 * x2 + x18 * x2 + x2 * x22 + x20);  // zz
  return result;
}
Eigen::VectorXd T3_damp_thole(const Eigen::Vector3d& rij, double a) {
  double x = rij(0);
  double y = rij(1);
  double z = rij(2);
  Eigen::VectorXd result(10);
  double x0  = x * x;
  double x1  = y * y;
  double x2  = z * z;
  double x3  = x0 + x1 + x2;
  double x4  = 1.0 / x3;
  double x5  = x0 * x4;
  double x6  = 3.0 * x5 - 1.0;
  double x7  = sqrt(x3);
  double x8  = a * x7;
  double x9  = exp(-x8);
  double x10 = x8 + 1.0;
  double x11 = a * 1.0 / (x3 * x3);
  double x12 = x10 * x11;
  double x13 = 3.0 * x12 * x9;
  double x14 = x8 + 2.0;
  double x15 = x14 * x9 - 2.0;
  double x16 = pow(x3, -2.5);
  double x17 = 3.0 * x16;
  double x18 = x15 * x17;
  double x19 = pow(x3, -1.5);
  double x20 = 3.0 * x19;
  double x21 = x0 * x19;
  double x22 = x14 * x21;
  double x23 = a * x4;
  double x24 = 2.0 * x23;
  double x25 = x14 * x23;
  double x26 = 1.0 / x7;
  double x27 = -x14 * x26 + x26;
  double x28 = -x0 * x24 + x0 * x25 - x21 + x22 + x27;
  double x29 = a * x9;
  double x30 = x28 * x29;
  double x31 = x0 * x17;
  double x32 = 6.0 * x11;
  double x33 = a * a;
  double x34 = x20 * x33;
  double x35 = x14 * x31;
  double x36 = x11 * x14;
  double x37 = 3.0 * x36;
  double x38 = -x14 * x20 + x20 + 6.0 * x23 - 3.0 * x25;
  double x39 = x26 * x29;
  double x40 = 0.5 * x;
  double x41 = x10 * x29 * 1.0 / (x3 * x3 * x3);
  double x42 = 3.0 * x41;
  double x43 = x15 * pow(x3, -3.5);
  double x44 = 3.0 * x43;
  double x45 = 0.5 * x9;
  double x46 = x12 * x45;
  double x47 = 1.5 * x15 * x16;
  double x48 = x14 * x19;
  double x49 = -x19 + x48;
  double x50 = x29 * (-x24 + x25 + x49);
  double x51 = 0.5 * x19;
  double x52 = x33 * x4 * x45;
  double x53 = 5.0 * x11;
  double x54 = 2.0 * x36;
  double x55 = -x23 + x49;
  double x56 = 0.5 * x39;
  double x57 = x0 * x42 + x0 * x44 + x21 * x50 + x28 * x52 + x30 * x51 + x46 * x6 +
               x47 * x6 - x56 * (x0 * x53 - x0 * x54 + x21 * x33 + x31 - x35 + x55);
  double x58 = 7.5 * x43;
  double x59 = x19 + x24 - x25 - x48;
  double x60 = x1 * x17;
  double x61 = x14 * x60;
  double x62 = x1 * x48;
  double x63 = -x1 * x32 - x1 * x34 + x1 * x37 + x33 * x62 - x60 + x61;
  double x64 = 4.5 * x41;
  double x65 = -x46 - x47;
  double x66 = x1 * x19;
  double x67 = -x1 * x24 + x1 * x25 + x27 + x62 - x66;
  double x68 = x29 * x51;
  double x69 = x50 * x66 + x67 * x68;
  double x70 = x14 * x17;
  double x71 = x33 * x48;
  double x72 = x19 * x2;
  double x73 = -x2 * x24 + x2 * x25 + x2 * x48 + x27 - x72;
  double x74 = -x17 * x2 - x2 * x32 - x2 * x34 + x2 * x37 + x2 * x70 + x2 * x71;
  double x75 = x2 * x58 + x2 * x64 + x50 * x72 + x56 * (x59 + x74) + x65 + x68 * x73;
  double x76 = x1 * x4;
  double x77 = 3.0 * x76 - 1.0;
  double x78 = x20 * x29;
  double x79 = x2 * x4;
  result[0] =
        x40 *
        (x13 * x6 + x18 * (5.0 * x5 - 3.0) + x20 * x30 +
         x39 * (-x0 * x32 - x0 * x34 + x0 * x37 + x22 * x33 - x31 + x35 + x38));  // xxx
  result[1] = x57 * y;                                                            // xxy
  result[2] = x57 * z;                                                            // xxz
  result[3] = x * (x1 * x58 + x1 * x64 + x56 * (x59 + x63) + x65 + x69);          // xyy
  result[4] = x40 * y * z *
              (x20 * x50 + x39 * (-x17 - x32 - x34 + x37 + x70 + x71) + 9.0 * x41 +
               15.0 * x43);  // xyz
  result[5] = x * x75;       // xzz
  result[6] =
        0.5 * y *
        (x13 * x77 + x18 * (5.0 * x76 - 3.0) + x39 * (x38 + x63) + x67 * x78);  // yyy
  result[7] =
        z * (x1 * x42 + x1 * x44 + x46 * x77 + x47 * x77 + x52 * x67 -
             x56 * (x1 * x53 - x1 * x54 + x33 * x66 + x55 + x60 - x61) + x69);  // yyz
  result[8] = x75 * y;                                                          // yzz
  result[9] = 0.5 * z *
              (x13 * (3.0 * x79 - 1.0) + x18 * (5.0 * x79 - 3.0) + x39 * (x38 + x74) +
               x73 * x78);  // zzz
  return result;
}
Eigen::VectorXd T4_damp_thole(const Eigen::Vector3d& rij, double a) {
  double x = rij(0);
  double y = rij(1);
  double z = rij(2);
  Eigen::VectorXd result(15);
  double x0   = x * x;
  double x1   = y * y;
  double x2   = z * z;
  double x3   = x0 + x1 + x2;
  double x4   = 1.0 / x3;
  double x5   = x0 * x4;
  double x6   = 5.0 * x5 - 3.0;
  double x7   = a * 1.0 / (x3 * x3 * x3);
  double x8   = x0 * x7;
  double x9   = sqrt(x3);
  double x10  = a * x9;
  double x11  = exp(-x10);
  double x12  = x10 + 1.0;
  double x13  = x11 * x12;
  double x14  = 6.0 * x13;
  double x15  = x * x * x * x;
  double x16  = 1.0 / (x3 * x3);
  double x17  = 35.0 * x16;
  double x18  = x10 + 2.0;
  double x19  = x11 * x18 - 2.0;
  double x20  = pow(x3, -2.5);
  double x21  = 1.5 * x20;
  double x22  = x19 * x21;
  double x23  = 3.0 * x5 - 1.0;
  double x24  = pow(x3, -1.5);
  double x25  = 3.0 * x24;
  double x26  = x0 * x24;
  double x27  = x18 * x26;
  double x28  = a * x4;
  double x29  = 2.0 * x28;
  double x30  = x18 * x28;
  double x31  = 1.0 / x9;
  double x32  = -x18 * x31 + x31;
  double x33  = -x0 * x29 + x0 * x30 - x26 + x27 + x32;
  double x34  = x11 * x33;
  double x35  = a * x34;
  double x36  = a * x16;
  double x37  = 6.0 * x36;
  double x38  = a * a;
  double x39  = x25 * x38;
  double x40  = x27 * x38;
  double x41  = 3.0 * x36;
  double x42  = x18 * x41;
  double x43  = 3.0 * x20;
  double x44  = x0 * x43;
  double x45  = x18 * x44;
  double x46  = -x44 + x45;
  double x47  = x18 * x25;
  double x48  = 6.0 * x28;
  double x49  = 3.0 * x30;
  double x50  = x25 - x47 + x48 - x49;
  double x51  = -x0 * x37 - x0 * x39 + x0 * x42 + x40 + x46 + x50;
  double x52  = a * x11;
  double x53  = x26 * x52;
  double x54  = pow(x3, -3.5);
  double x55  = 15.0 * x54;
  double x56  = x15 * x55;
  double x57  = 18.0 * x20;
  double x58  = x0 * x57;
  double x59  = 30.0 * x7;
  double x60  = x38 * x57;
  double x61  = a * a * a;
  double x62  = x16 * x61;
  double x63  = x15 * x62;
  double x64  = x26 * x38;
  double x65  = x0 * x36;
  double x66  = 18.0 * x36;
  double x67  = x18 * x66;
  double x68  = x15 * x18;
  double x69  = x20 * x38;
  double x70  = 6.0 * x69;
  double x71  = 15.0 * x7;
  double x72  = -x25 + x47 - x48 + x49;
  double x73  = 0.5 * x11;
  double x74  = a * x31 * x73;
  double x75  = 9.0 * x12;
  double x76  = 1.0 / (x3 * x3 * x3 * x3);
  double x77  = x0 * x52;
  double x78  = x76 * x77;
  double x79  = x19 * pow(x3, -4.5);
  double x80  = x0 * x79;
  double x81  = x11 * x7;
  double x82  = x12 * x81;
  double x83  = 4.5 * x82;
  double x84  = 1.5 * x82;
  double x85  = x19 * x54;
  double x86  = 7.5 * x85;
  double x87  = x18 * x24;
  double x88  = -x24 + x87;
  double x89  = -x29 + x30 + x88;
  double x90  = a * x24;
  double x91  = x11 * x90;
  double x92  = 1.5 * x91;
  double x93  = x89 * x92;
  double x94  = 4.5 * x20;
  double x95  = 1.5 * x16 * x38;
  double x96  = x34 * x95;
  double x97  = 5.0 * x65;
  double x98  = 2.0 * x18 * x65;
  double x99  = -x28 + x88;
  double x100 = x11 * (x44 - x45 + x64 + x97 - x98 + x99);
  double x101 = x100 * x90;
  double x102 = 0.5 * x91;
  double x103 = x38 * x4;
  double x104 = x103 * x73;
  double x105 = 12.0 * x0;
  double x106 = x18 * x7;
  double x107 = x0 * x55;
  double x108 = x107 * x18;
  double x109 = x107 - x108;
  double x110 = -x39;
  double x111 = 9.0 * x20;
  double x112 = x111 * x18;
  double x113 = x110 - x111 + x112 + x18 * x37 - 15.0 * x36;
  double x114 = x * (-1.5 * x101 + x102 * x51 + x104 * x51 + x23 * x83 + x23 * x93 +
                     x35 * x94 + x6 * x84 + x6 * x86 -
                     x74 * (x0 * x62 - x105 * x106 + x105 * x69 + x109 + x113 -
                            x38 * x45 + 27.0 * x8) +
                     x75 * x78 + 15.0 * x80 + x96);
  double x115 = x1 * x41;
  double x116 = 5.0 * x69;
  double x117 = x0 * x116;
  double x118 = 23.0 * x8;
  double x119 = 8.0 * x18;
  double x120 = x119 * x8;
  double x121 = x24 - x87;
  double x122 = x1 * x43;
  double x123 = x122 * x18;
  double x124 = -x122 + x123;
  double x125 = x121 + x124 + x28;
  double x126 = x46 - x64 - x97 + x98;
  double x127 = 30.0 * x80;
  double x128 = x23 * x86;
  double x129 = x1 * x24;
  double x130 = a * x100;
  double x131 = x1 * x4;
  double x132 = x100 * x38;
  double x133 = x121 + x29 - x30;
  double x134 = x1 * x87;
  double x135 = x134 * x38;
  double x136 = -x1 * x37 - x1 * x39 + x115 * x18 + x124 + x135;
  double x137 = x133 + x136;
  double x138 = x21 * x35;
  double x139 = -x1 * x29 + x1 * x30 - x129 + x134 + x32;
  double x140 = x102 * x139;
  double x141 = x61 * x73;
  double x142 = x141 * x33;
  double x143 = 21.0 * x12;
  double x144 = x143 * x78;
  double x145 = x20 * x89;
  double x146 = 6.0 * x145;
  double x147 = x146 * x77;
  double x148 = x1 * x7;
  double x149 = 3.0 * x13;
  double x150 = x148 * x149;
  double x151 = 3.0 * x85;
  double x152 = x0 * x151 + x102 * x33 + x104 * x33 + x149 * x8 + x22 * x23;
  double x153 = 3.0 * x82;
  double x154 = x153 * x23;
  double x155 = x102 * x89;
  double x156 = x38 * x87;
  double x157 = x18 * x43;
  double x158 = x157 - x43;
  double x159 = x110 + x156 + x158 - x37 + x42;
  double x160 = y * z;
  double x161 = x157 * x2 - x2 * x43;
  double x162 = x161 - x2 * x41;
  double x163 = x2 * x24;
  double x164 = x2 * x4;
  double x165 = x156 * x2;
  double x166 = x161 + x165 - x2 * x37 - x2 * x39 + x2 * x42;
  double x167 = x133 + x166;
  double x168 = -x163 - x2 * x29 + x2 * x30 + x2 * x87 + x32;
  double x169 = x102 * x168;
  double x170 = x18 * x36;
  double x171 = x111 - x112 - 9.0 * x170 + 9.0 * x24 * x38 - x38 * x47 + x66;
  double x172 = x1 * x55;
  double x173 = x1 * x57;
  double x174 = x1 * x62;
  double x175 = x172 * x18;
  double x176 = x1 * x69;
  double x177 = 15.0 * x1;
  double x178 = -x1 * x59 + x106 * x177 - x172 - x173 * x38 + x174 * x18 - 4.0 * x174 +
                x175 + 6.0 * x176 * x18;
  double x179 = x137 * x91;
  double x180 = x1 * x52;
  double x181 = x89 * x94;
  double x182 = 52.5 * x79;
  double x183 = x180 * x76;
  double x184 = 30.0 * x12;
  double x185 = x1 * x182 + x183 * x184;
  double x186 = x136 + x50;
  double x187 = x139 * x52;
  double x188 = x102 * x186 + x187 * x94;
  double x189 = -x75 * x81 - 22.5 * x85 - x93;
  double x190 = x * y;
  double x191 = -x156 - x157 + x37 + x39 - x42 + x43;
  double x192 = x187 * x21;
  double x193 = x159 * x52;
  double x194 = 7.5 * x145;
  double x195 = -x153 - x155 - x86;
  double x196 = x * z;
  double x197 = x167 * x91;
  double x198 = x2 * x55;
  double x199 = x2 * x57;
  double x200 = x2 * x62;
  double x201 = x18 * x2;
  double x202 = x18 * x198 + x18 * x200 - x198 - x199 * x38 - x2 * x59 - 4.0 * x200 +
                x201 * x70 + x201 * x71;
  double x203 = x168 * x52;
  double x204 = x2 * x52;
  double x205 = x182 * x2 + x184 * x204 * x76;
  double x206 = x166 + x50;
  double x207 = x102 * x206 + x181 * x204 + x189 + 1.5 * x197 + x203 * x94 + x205 +
                x74 * (x171 + x202);
  double x208 = 5.0 * x131 - 3.0;
  double x209 = y * y * y * y;
  double x210 = 3.0 * x131 - 1.0;
  double x211 = x129 * x52;
  double x212 = x209 * x55;
  double x213 = x209 * x62;
  double x214 = x129 * x38;
  double x215 = x1 * x36;
  double x216 = x18 * x209;
  double x217 = 5.0 * x215;
  double x218 = 2.0 * x1 * x170;
  double x219 = x122 - x123 + x214 + x217 - x218 + x99;
  double x220 = x11 * x139 * x95;
  double x221 = x1 * x2;
  double x222 = x221 * x52;
  double x223 = x2 * x210;
  double x224 = x163 * x52;
  double x225 = x148 * x2;
  double x226 = z * z * z * z;
  double x227 = x226 * x55;
  double x228 = x226 * x62;
  double x229 = x18 * x226;
  result[0]   = -x14 * x6 * x8 - x22 * (x15 * x17 - 30.0 * x5 + 3.0) - x23 * x25 * x35 -
              2.0 * x51 * x53 -
              x74 * (-x0 * x67 - x15 * x59 - x15 * x60 + x18 * x56 - x18 * x58 +
                     x18 * x63 - 6.0 * x40 - x56 + x58 - 4.0 * x63 + 18.0 * x64 +
                     36.0 * x65 + x68 * x70 + x68 * x71 + x72);  // xxxx
  result[1] = -x114 * y;                                         // xxxy
  result[2] = -x114 * z;                                         // xxxz
  result[3] = -x1 * x127 - x1 * x128 - x1 * x138 - x1 * x144 - x1 * x147 - x1 * x96 +
              x129 * x130 - x129 * x142 + x131 * x132 - x137 * x53 - x140 * x23 -
              x150 * x23 + x152 +
              x74 * (x1 * x107 - x1 * x108 + x1 * x117 + x1 * x118 - x1 * x120 - x115 +
                     x125 + x126);  // xxyy
  result[4] = -x160 * (-x100 * x103 - x101 + x127 + x128 + x138 + x142 * x24 + x144 +
                       x147 + x154 + x155 * x23 + x159 * x53 -
                       x74 * (x109 + x117 + x118 - x120 + x158 - x41) + x96);  // xxyz
  result[5] = -x127 * x2 - x128 * x2 + x130 * x163 + x132 * x164 - x138 * x2 -
              x142 * x163 - x144 * x2 - x147 * x2 + x152 - x154 * x2 - x167 * x53 -
              x169 * x23 - x2 * x96 +
              x74 * (x107 * x2 - x108 * x2 + x117 * x2 + x118 * x2 - x120 * x2 + x121 +
                     x126 + x162 + x28);  // xxzz
  result[6]  = -x190 * (1.5 * x179 + x180 * x181 + x185 + x188 + x189 +
                       x74 * (x171 + x178));  // xyyy
  result[7]  = -x196 * (x129 * x193 + x179 + x180 * x194 + x185 + x192 + x195 +
                       x74 * (x178 + x191));  // xyyz
  result[8]  = -x190 * (x163 * x193 + x194 * x204 + x195 + x197 + x203 * x21 + x205 +
                       x74 * (x191 + x202));  // xyzz
  result[9]  = -x196 * x207;                   // xzzz
  result[10] = -x14 * x148 * x208 - 2.0 * x186 * x211 - x187 * x210 * x25 -
               x22 * (-30.0 * x131 + x17 * x209 + 3.0) -
               x74 * (-x1 * x67 - 6.0 * x135 - x173 * x18 + x173 + x18 * x212 +
                      x18 * x213 - x209 * x59 - x209 * x60 - x212 - 4.0 * x213 +
                      18.0 * x214 + 36.0 * x215 + x216 * x70 + x216 * x71 + x72);  // yyyy
  result[11] = -x160 * (x104 * x186 + x177 * x79 + x183 * x75 + x188 + x208 * x84 +
                        x208 * x86 + x210 * x83 + x210 * x93 - x219 * x92 + x220 -
                        x74 * (x113 - x123 * x38 - 12.0 * x148 * x18 + 27.0 * x148 +
                               x172 + x174 - x175 + 12.0 * x176));  // yyyz
  result[12] = x1 * x151 + x104 * x139 + x11 * x164 * x219 * x38 - x139 * x141 * x163 +
               x140 - x143 * x222 * x76 - x146 * x222 + x150 - x153 * x223 - x167 * x211 -
               x169 * x210 - x192 * x2 - x2 * x220 + x210 * x22 + x219 * x224 -
               30.0 * x221 * x79 - x223 * x86 +
               x74 * (x116 * x221 - x119 * x225 + x125 + x162 + x172 * x2 - x175 * x2 -
                      x214 - x217 + x218 + 23.0 * x225);  // yyzz
  result[13] = -x160 * x207;                              // yzzz
  result[14] =
        -6.0 * x2 * x82 * (5.0 * x164 - 3.0) - x203 * x25 * (3.0 * x164 - 1.0) -
        2.0 * x206 * x224 - x22 * (-30.0 * x164 + x17 * x226 + 3.0) -
        x74 * (18.0 * x163 * x38 - 6.0 * x165 - x18 * x199 + x18 * x227 + x18 * x228 +
               x199 + 36.0 * x2 * x36 - x2 * x67 - x226 * x59 - x226 * x60 - x227 -
               4.0 * x228 + x229 * x70 + x229 * x71 + x72);  // zzzz
  return result;
}
Eigen::VectorXd T5_damp_thole(const Eigen::Vector3d& rij, double a) {
  double x = rij(0);
  double y = rij(1);
  double z = rij(2);
  Eigen::VectorXd result(21);
  double x0  = x * x;
  double x1  = y * y;
  double x2  = z * z;
  double x3  = x0 + x1 + x2;
  double x4  = 1.0 / x3;
  double x5  = x0 * x4;
  double x6  = x * x * x * x;
  double x7  = 1.0 / (x3 * x3);
  double x8  = x6 * x7;
  double x9  = -30.0 * x5 + 35.0 * x8 + 3.0;
  double x10 = 1.0 / (x3 * x3 * x3);
  double x11 = a * x10;
  double x12 = sqrt(x3);
  double x13 = a * x12;
  double x14 = exp(-x13);
  double x15 = x13 + 1.0;
  double x16 = x14 * x15;
  double x17 = x11 * x16;
  double x18 = 7.5 * x17;
  double x19 = pow(x3, -3.5);
  double x20 = x13 + 2.0;
  double x21 = x14 * x20 - 2.0;
  double x22 = x19 * x21;
  double x23 = 7.5 * x22;
  double x24 = pow(x3, -1.5);
  double x25 = x0 * x24;
  double x26 = x20 * x25;
  double x27 = a * x4;
  double x28 = 2.0 * x27;
  double x29 = x20 * x27;
  double x30 = 1.0 / x12;
  double x31 = -x20 * x30 + x30;
  double x32 = -x0 * x28 + x0 * x29 - x25 + x26 + x31;
  double x33 = x14 * x32;
  double x34 = pow(x3, -2.5);
  double x35 = a * x34;
  double x36 = 5.0 * x5 - 3.0;
  double x37 = 15.0 * x36;
  double x38 = 3.0 * x5 - 1.0;
  double x39 = a * x7;
  double x40 = 6.0 * x39;
  double x41 = a * a;
  double x42 = 3.0 * x24;
  double x43 = x41 * x42;
  double x44 = x26 * x41;
  double x45 = 3.0 * x39;
  double x46 = x20 * x45;
  double x47 = 3.0 * x34;
  double x48 = x0 * x47;
  double x49 = x20 * x48;
  double x50 = -x48 + x49;
  double x51 = x20 * x42;
  double x52 = 6.0 * x27;
  double x53 = 3.0 * x29;
  double x54 = x42 - x51 + x52 - x53;
  double x55 = -x0 * x40 - x0 * x43 + x0 * x46 + x44 + x50 + x54;
  double x56 = a * x24;
  double x57 = x14 * x56;
  double x58 = 5.0 * x57;
  double x59 = 15.0 * x19;
  double x60 = x59 * x6;
  double x61 = 18.0 * x34;
  double x62 = x0 * x61;
  double x63 = 30.0 * x11;
  double x64 = x41 * x6;
  double x65 = x20 * x62;
  double x66 = a * a * a;
  double x67 = x66 * x8;
  double x68 = x25 * x41;
  double x69 = x0 * x39;
  double x70 = 18.0 * x39;
  double x71 = x20 * x70;
  double x72 = x34 * x41;
  double x73 = 6.0 * x20;
  double x74 = x72 * x73;
  double x75 = x11 * x20;
  double x76 = 15.0 * x75;
  double x77 = -x42 + x51 - x52 + x53;
  double x78 = -x0 * x71 + x20 * x60 + x20 * x67 - 6.0 * x44 - x6 * x63 + x6 * x74 +
               x6 * x76 - x60 - x61 * x64 + x62 - x65 - 4.0 * x67 + 18.0 * x68 +
               36.0 * x69 + x77;
  double x79  = x57 * x78;
  double x80  = pow(x3, -4.5);
  double x81  = 105.0 * x80;
  double x82  = x6 * x81;
  double x83  = x0 * x19;
  double x84  = 150.0 * x83;
  double x85  = a * 1.0 / (x3 * x3 * x3 * x3);
  double x86  = x6 * x85;
  double x87  = x19 * x64;
  double x88  = x10 * x66;
  double x89  = x6 * x88;
  double x90  = a * a * a * a;
  double x91  = x34 * x90;
  double x92  = x6 * x91;
  double x93  = x66 * x7;
  double x94  = x0 * x93;
  double x95  = x20 * x82;
  double x96  = x0 * x72;
  double x97  = x0 * x11;
  double x98  = x20 * x97;
  double x99  = 60.0 * x20;
  double x100 = 10.0 * x20;
  double x101 = 45.0 * x19;
  double x102 = x101 * x20;
  double x103 = x20 * x86;
  double x104 = 45.0 * x34;
  double x105 = x24 * x41;
  double x106 = x20 * x24;
  double x107 = x106 * x41;
  double x108 = x20 * x39;
  double x109 = x104 * x20 - x104 - 45.0 * x105 + 15.0 * x107 + 45.0 * x108 - 90.0 * x39;
  double x110 = 0.5 * x14;
  double x111 = a * x30;
  double x112 = x110 * x111;
  double x113 = a * x14;
  double x114 = 60.0 * x113;
  double x115 = x15 * 1.0 / (x3 * x3 * x3 * x3 * x3);
  double x116 = x114 * x115;
  double x117 = x0 * x85;
  double x118 = x117 * x16;
  double x119 = x21 * x80;
  double x120 = x0 * x119;
  double x121 = 30.0 * x120;
  double x122 = 1.5 * x14;
  double x123 = x11 * x122 * x15;
  double x124 = x106 - x24;
  double x125 = x124 - x28 + x29;
  double x126 = x14 * x35;
  double x127 = x125 * x126;
  double x128 = 6.0 * x127;
  double x129 = a * x33;
  double x130 = 18.0 * x83;
  double x131 = 9.0 * x34;
  double x132 = a * x131;
  double x133 = x33 * x38;
  double x134 = x41 * x7;
  double x135 = 3.0 * x134;
  double x136 = 5.0 * x69;
  double x137 = 2.0 * x108;
  double x138 = x0 * x137;
  double x139 = x124 - x27;
  double x140 = x136 - x138 + x139 + x48 - x49 + x68;
  double x141 = x14 * x140;
  double x142 = x0 * x55;
  double x143 = 6.0 * x126;
  double x144 = x134 * x14;
  double x145 = 2.0 * x144;
  double x146 = 12.0 * x96;
  double x147 = 27.0 * x97;
  double x148 = 12.0 * x98;
  double x149 = x41 * x49;
  double x150 = 15.0 * x83;
  double x151 = x150 * x20;
  double x152 = x150 - x151;
  double x153 = -x43;
  double x154 = 15.0 * x39;
  double x155 = x20 * x40;
  double x156 = x131 * x20;
  double x157 = -x131 + x156;
  double x158 = x153 - x154 + x155 + x157;
  double x159 = x14 * (x146 + x147 - x148 - x149 + x152 + x158 + x94);
  double x160 = a * x159;
  double x161 = x4 * x41;
  double x162 = x110 * x161;
  double x163 = 90.0 * x83;
  double x164 = 30.0 * x20;
  double x165 = 4.0 * x20;
  double x166 = x131 - x156;
  double x167 = x154 - x155 + x43;
  double x168 = x166 + x167;
  double x169 = -a * x141 * x38 * x42 + x0 * x128 * x36 -
                x112 * (-90.0 * x103 + x163 * x20 - x163 - x164 * x87 - x165 * x89 +
                        x168 + x41 * x65 + x82 + 195.0 * x86 + 105.0 * x87 + 22.0 * x89 +
                        x92 - 6.0 * x94 - x95 - 72.0 * x96 - 162.0 * x97 + 72.0 * x98) +
                x116 * x6 + 30.0 * x118 * x36 + x121 * (7.0 * x5 - 3.0) + x123 * x9 +
                x129 * x130 + x132 * x133 + x133 * x135 + x142 * x143 + x142 * x145 -
                2.0 * x160 * x25 + x162 * x78 + x23 * x9 + 0.5 * x79;
  double x170 = x21 * pow(x3, -5.5);
  double x171 = x0 * x170;
  double x172 = 210.0 * x171;
  double x173 = x1 * x101;
  double x174 = 69.0 * x11;
  double x175 = 15.0 * x72;
  double x176 = x173 * x20;
  double x177 = 7.0 * x88;
  double x178 = x0 * x177;
  double x179 = 24.0 * x75;
  double x180 = x41 * x83;
  double x181 = 72.0 * x180;
  double x182 = 177.0 * x117;
  double x183 = x117 * x20;
  double x184 = 72.0 * x183;
  double x185 = x1 * x151;
  double x186 = x0 * x81;
  double x187 = x1 * x186;
  double x188 = -x187 * x20 + x187;
  double x189 = -x150 + x151;
  double x190 = -x146 - x147 + x148 + x149 + x168 + x189 - x94;
  double x191 = 52.5 * x119;
  double x192 = x191 * x36;
  double x193 = x1 * x24;
  double x194 = x193 * x66;
  double x195 = x110 * x55;
  double x196 = x1 * x4;
  double x197 = x159 * x41;
  double x198 = x132 * x141;
  double x199 = x135 * x141;
  double x200 = x122 * x55;
  double x201 = x200 * x35;
  double x202 = x1 * x106;
  double x203 = -x1 * x28 + x1 * x29 - x193 + x202 + x31;
  double x204 = x122 * x35;
  double x205 = x203 * x204;
  double x206 = -x106 + x24;
  double x207 = x206 + x28 - x29;
  double x208 = x1 * x40;
  double x209 = x1 * x43;
  double x210 = x202 * x41;
  double x211 = x1 * x45;
  double x212 = x20 * x211;
  double x213 = x1 * x47;
  double x214 = x20 * x213;
  double x215 = -x213 + x214;
  double x216 = -x208 - x209 + x210 + x212 + x215;
  double x217 = x207 + x216;
  double x218 = 1.5 * x57;
  double x219 = x218 * x38;
  double x220 = x134 * x200;
  double x221 = 1.5 * x33;
  double x222 = x34 * x66;
  double x223 = x221 * x222;
  double x224 = x10 * x41;
  double x225 = x224 * x33;
  double x226 = 10.5 * x225;
  double x227 = 22.5 * x19;
  double x228 = x129 * x227;
  double x229 = x132 * x14;
  double x230 = x125 * x229;
  double x231 = x230 * x38;
  double x232 = x16 * x85;
  double x233 = x232 * x37;
  double x234 = x113 * x125;
  double x235 = x130 * x234;
  double x236 = x0 * x1;
  double x237 = x113 * x115;
  double x238 = 120.0 * x237;
  double x239 = -x134 * x221;
  double x240 = 22.5 * x232;
  double x241 = x240 * x38;
  double x242 = x1 * x241;
  double x243 = x239 + x242;
  double x244 = 5.0 * x96;
  double x245 = 23.0 * x97;
  double x246 = 8.0 * x98;
  double x247 = x206 + x215 + x27;
  double x248 = -x136 + x138 + x50 - x68;
  double x249 = x1 * x150 + x1 * x244 + x1 * x245 - x1 * x246 - x185 - x211 + x247 + x248;
  double x250 = 4.5 * x33;
  double x251 = 4.5 * x17;
  double x252 = -x250 * x35 - x251 * x38;
  double x253 = -x218 * x249 + x252;
  double x254 = 0.5 * x57;
  double x255 = -9.0 * x118 - 15.0 * x120 - x162 * x55 - x23 * x36 - x254 * x55;
  double x256 = x0 * x238;
  double x257 = x20 * x47;
  double x258 = x257 - x47;
  double x259 = x107 + x153 + x258 - x40 + x46;
  double x260 = x218 * x259;
  double x261 = x152 + x244 + x245 - x246 + x258 - x45;
  double x262 = x24 * x66;
  double x263 = x110 * x262;
  double x264 = x186 * x20;
  double x265 = x151 * x41;
  double x266 = x * y * z;
  double x267 = x186 * x2 - x2 * x264;
  double x268 = x101 * x2;
  double x269 = x102 * x2;
  double x270 = -x174 * x2 - x175 * x2 + x179 * x2 - x268 + x269;
  double x271 = x2 * x24;
  double x272 = x271 * x66;
  double x273 = x2 * x4;
  double x274 = x106 * x2 - x2 * x28 + x2 * x29 - x271 + x31;
  double x275 = x204 * x274;
  double x276 = x2 * x40;
  double x277 = x2 * x43;
  double x278 = x107 * x2;
  double x279 = x2 * x46;
  double x280 = x2 * x47;
  double x281 = x2 * x257;
  double x282 = -x280 + x281;
  double x283 = -x276 - x277 + x278 + x279 + x282;
  double x284 = x207 + x283;
  double x285 = x2 * x241;
  double x286 = x239 + x285;
  double x287 = -x2 * x45 + x282;
  double x288 = x150 * x2 - x151 * x2 + x2 * x244 + x2 * x245 - x2 * x246 + x206 + x248 +
                x27 + x287;
  double x289 = -x218 * x288 + x252;
  double x290 = 5.0 * x19;
  double x291 = x1 * x290;
  double x292 = x1 * x11;
  double x293 = 35.0 * x80;
  double x294 = x236 * x293;
  double x295 = x1 * x180;
  double x296 = x1 * x117;
  double x297 = x1 * x183;
  double x298 = -x257 + x47;
  double x299 = x189 - x244 - x245 + x246 + x45;
  double x300 = x298 + x299;
  double x301 = x111 * x122;
  double x302 = x122 * x161;
  double x303 = x1 * x59;
  double x304 = x20 * x303;
  double x305 = -x303 + x304;
  double x306 = x166 + x305;
  double x307 = 9.0 * x105;
  double x308 = 9.0 * x108;
  double x309 = x41 * x51;
  double x310 = x307 - x308 - x309 + x70;
  double x311 = x1 * x63;
  double x312 = x1 * x61;
  double x313 = x312 * x41;
  double x314 = x1 * x93;
  double x315 = 4.0 * x314;
  double x316 = x20 * x314;
  double x317 = x1 * x72;
  double x318 = x317 * x73;
  double x319 = x1 * x76;
  double x320 = -x311 - x313 - x315 + x316 + x318 + x319;
  double x321 = x306 + x310 + x320;
  double x322 = x113 * x25;
  double x323 = x216 + x54;
  double x324 = x254 * x323;
  double x325 = x113 * x83;
  double x326 = x203 * x325;
  double x327 = x217 * x229;
  double x328 = 4.5 * x203;
  double x329 = x126 * x328;
  double x330 = x0 * x234;
  double x331 = 22.5 * x22;
  double x332 = 1.5 * x141;
  double x333 = -x0 * x230 - 54.0 * x118 - 90.0 * x120 - x134 * x250 + x140 * x218 +
                x161 * x332 - x221 * x262 - x331 * x38;
  double x334 = 315.0 * x171;
  double x335 = x1 * x191;
  double x336 = x110 * x7 * x90;
  double x337 = x32 * x336;
  double x338 = x33 * x66;
  double x339 = 4.5 * x141;
  double x340 = x1 * x339;
  double x341 = 7.5 * x1;
  double x342 = x129 * x19;
  double x343 = 195.0 * x237;
  double x344 = x1 * x334 + x1 * x337 - x134 * x340 - x194 * x332 + x213 * x338 +
                x225 * x341 + x236 * x343 + x335 * x38 - x340 * x35 + x341 * x342;
  double x345 = x298 + x305;
  double x346 = -x107 + x40 + x43 - x46;
  double x347 = x320 + x345 + x346;
  double x348 = x217 * x254;
  double x349 = x113 * x261;
  double x350 = x14 * x41;
  double x351 = x261 * x350;
  double x352 = x113 * x48;
  double x353 = x213 * x234;
  double x354 = x143 * x259;
  double x355 = 51.0 * x234 * x83;
  double x356 = -18.0 * x118 - x121 - x123 * x38 + x140 * x162 + x140 * x254 -
                x221 * x35 - x23 * x38 - x234 * x48 - x263 * x32;
  double x357 = x11 * x2;
  double x358 = x180 * x2;
  double x359 = x117 * x2;
  double x360 = x183 * x2;
  double x361 = x2 * x59;
  double x362 = -x361;
  double x363 = x20 * x361;
  double x364 = x298 + x362 + x363;
  double x365 = x2 * x63;
  double x366 = x2 * x61;
  double x367 = x366 * x41;
  double x368 = x2 * x93;
  double x369 = 4.0 * x368;
  double x370 = x20 * x368;
  double x371 = x2 * x72;
  double x372 = x371 * x73;
  double x373 = x2 * x76;
  double x374 = -x365 - x367 - x369 + x370 + x372 + x373;
  double x375 = x346 + x364 + x374;
  double x376 = x254 * x284;
  double x377 = x274 * x325;
  double x378 = x234 * x280;
  double x379 = x143 * x2;
  double x380 = x191 * x2;
  double x381 = x2 * x339;
  double x382 = 7.5 * x2;
  double x383 = x0 * x2 * x343 - x134 * x381 + x2 * x334 + x2 * x337 + x225 * x382 -
                x272 * x332 + x280 * x338 + x342 * x382 - x35 * x381 + x38 * x380;
  double x384 = x2 * x293;
  double x385 = x0 * x384;
  double x386 = x2 * x290;
  double x387 = x20 * x386 - 5.0 * x357 - x386;
  double x388 = x166 + x310 + x362 + x363 + x374;
  double x389 = x283 + x54;
  double x390 = x254 * x389;
  double x391 = x229 * x284;
  double x392 = 4.5 * x126;
  double x393 = x274 * x392;
  double x394 = 315.0 * x119;
  double x395 = y * y * y * y;
  double x396 = 472.5 * x170;
  double x397 = 90.0 * x19;
  double x398 = x1 * x397;
  double x399 = x20 * x398;
  double x400 = 108.0 * x72;
  double x401 = 90.0 * x11;
  double x402 = x1 * x20;
  double x403 = 36.0 * x20;
  double x404 = x157 - x307 + x308 + x309 - x70;
  double x405 = x395 * x81;
  double x406 = 210.0 * x85;
  double x407 = x395 * x41;
  double x408 = 135.0 * x19;
  double x409 = 40.0 * x88;
  double x410 = x395 * x91;
  double x411 = x20 * x405;
  double x412 = x395 * x88;
  double x413 = x395 * x85;
  double x414 = 105.0 * x20;
  double x415 = x100 * x412 + x102 * x407 + x20 * x410 - x395 * x406 - x395 * x409 -
                x405 - x407 * x408 - 5.0 * x410 + x411 + x413 * x414;
  double x416 = x203 * x229;
  double x417 = x113 * x42;
  double x418 = 135.0 * x232;
  double x419 = x113 * x193;
  double x420 = 2.0 * x419;
  double x421 = x19 * x234;
  double x422 = 30.0 * x421;
  double x423 = x113 * x203;
  double x424 = 262.5 * x237;
  double x425 = x251 + x331;
  double x426 = x395 * x59;
  double x427 = x395 * x93;
  double x428 = x193 * x41;
  double x429 = 36.0 * x39;
  double x430 = x1 * x429 - x1 * x71 - x20 * x312 + x20 * x426 + x20 * x427 - 6.0 * x210 +
                x312 - x395 * x63 + x395 * x74 + x395 * x76 - x407 * x61 - x426 -
                4.0 * x427 + 18.0 * x428 + x77;
  double x431 = x1 * x323;
  double x432 = x143 * x431 + x254 * x430;
  double x433 = x1 * x396;
  double x434 = x1 * x81;
  double x435 = x1 * x406;
  double x436 = x408 * x41;
  double x437 = x1 * x436;
  double x438 = x1 * x409;
  double x439 = x1 * x91;
  double x440 = 5.0 * x439;
  double x441 = x20 * x434;
  double x442 = x20 * x439;
  double x443 = x100 * x88;
  double x444 = x1 * x443;
  double x445 = x176 * x41;
  double x446 = x402 * x85;
  double x447 = 105.0 * x446;
  double x448 = x41 * x61;
  double x449 = x20 * x93;
  double x450 = x101 - x102 - x20 * x448 + x401 - 3.0 * x449 + 54.0 * x72 - 45.0 * x75 +
                12.0 * x93;
  double x451 = x204 * x323;
  double x452 = x113 * x227;
  double x453 = x203 * x452;
  double x454 = x1 * x392;
  double x455 = 52.5 * x421;
  double x456 = -157.5 * x119 - 13.5 * x127 - 67.5 * x232 - x260;
  double x457 = x2 * x240;
  double x458 = x1 * x2;
  double x459 = x19 * x458;
  double x460 = x113 * x274;
  double x461 = x19 * x460;
  double x462 = x19 * x423;
  double x463 = x382 * x462;
  double x464 = x2 * x392;
  double x465 = x113 * x271;
  double x466 = x2 * x303;
  double x467 = x2 * x304;
  double x468 = x213 - x214;
  double x469 = x2 * x434;
  double x470 = x2 * x441;
  double x471 = x303 - x304;
  double x472 = x2 * x81;
  double x473 = x2 * x91;
  double x474 = x414 * x85;
  double x475 = z * z * z * z;
  double x476 = x475 * x59;
  double x477 = x2 * x429 - x2 * x71 - x20 * x366 + x20 * x476 + 18.0 * x271 * x41 -
                6.0 * x278 + x366 - x448 * x475 + x449 * x475 - x475 * x63 + x475 * x74 +
                x475 * x76 - 4.0 * x475 * x93 - x476 + x77;
  double x478 = x2 * x397;
  double x479 = x475 * x81;
  double x480 = x475 * x91;
  double x481 = x102 * x41 * x475 + x20 * x479 + x20 * x480 - x406 * x475 - x409 * x475 -
                x436 * x475 + x443 * x475 + x474 * x475 - x479 - 5.0 * x480;
  double x482 = x112 * (-x2 * x20 * x401 + x2 * x400 - x20 * x478 + 180.0 * x357 +
                        24.0 * x368 - 6.0 * x370 - x371 * x403 + x404 + x478 + x481) +
                x2 * x391 - x2 * x394 - x2 * x418 - x229 * x274 - x234 * x366 +
                x254 * x477 + x268 * x460 - x284 * x417 + x379 * x389 +
                2.0 * x388 * x465 + x396 * x475 + x422 * x475 + x424 * x475 + x425;
  double x483 = x395 * x7;
  double x484 = -30.0 * x196 + 35.0 * x483 + 3.0;
  double x485 = 5.0 * x196 - 3.0;
  double x486 = 15.0 * x126;
  double x487 = 3.0 * x196 - 1.0;
  double x488 = 2.5 * x57;
  double x489 = 150.0 * x19;
  double x490 = x1 * x489;
  double x491 = x1 * x75;
  double x492 = x1 * x119;
  double x493 = x19 * x407;
  double x494 = 5.0 * x1 * x39;
  double x495 = x1 * x137;
  double x496 = x139 + x428 + x468 + x494 - x495;
  double x497 = 12.0 * x317;
  double x498 = 27.0 * x292;
  double x499 = 12.0 * x491;
  double x500 = x214 * x41;
  double x501 = x158 + x314 + x471 + x497 + x498 - x499 - x500;
  double x502 = x135 * x14;
  double x503 = x1 * x232;
  double x504 = x170 * x458;
  double x505 = x122 * x203;
  double x506 = x41 * x459;
  double x507 = x458 * x85;
  double x508 = x2 * x446;
  double x509 = x2 * x496;
  double x510 = x2 * x203;
  double x511 = x14 * x224;
  double x512 = 5.0 * x317;
  double x513 = 23.0 * x292;
  double x514 = 8.0 * x491;
  double x515 = x2 * x512 + x2 * x513 - x2 * x514 + x247 + x287 - x428 + x466 - x467 -
                x494 + x495;
  double x516 = -x218 * x515 - x251 * x487 - x329 + x457 * x487;
  double x517 = x1 * x384;
  double x518 = x475 * x7;
  double x519 = x2 * x489;
  result[0]   = x * (x112 * (x100 * x89 - x100 * x94 + x102 * x64 + 105.0 * x103 + x109 -
                           x20 * x84 + x20 * x92 - x82 + x84 - 210.0 * x86 - 135.0 * x87 -
                           40.0 * x89 - 5.0 * x92 + 40.0 * x94 + x95 - x96 * x99 +
                           180.0 * x96 + 300.0 * x97 - 150.0 * x98) +
                   x18 * x9 + x23 * (-70.0 * x5 + 63.0 * x8 + 15.0) + x33 * x35 * x37 +
                   x38 * x55 * x58 + 2.5 * x79);  // xxxxx
  result[1]   = x169 * y;                           // xxxxy
  result[2]   = x169 * z;                           // xxxxz
  result[3] =
        x * (x1 * x172 + x1 * x192 - x1 * x198 - x1 * x199 + x1 * x201 + x1 * x220 +
             x1 * x223 + x1 * x226 + x1 * x228 + x1 * x231 + x1 * x233 + x1 * x235 -
             x112 * (-x1 * x174 - x1 * x175 + x1 * x178 + x1 * x179 + x1 * x181 +
                     x1 * x182 - x1 * x184 - x173 + x176 - x185 * x41 + x188 + x190) -
             x160 * x193 + x194 * x195 - x196 * x197 + x205 * x36 + x217 * x219 +
             x236 * x238 + x243 + x253 + x255);  // xxxyy
  result[4] = x266 * (-x112 * (-x101 + x102 - x174 - x175 + x178 + x179 + x181 + x182 -
                               x184 + x186 - x264 - x265) +
                      x125 * x204 * x36 - x159 * x161 - x159 * x56 + x172 + x192 - x198 -
                      x199 + x201 - x218 * x261 + x220 + x223 + x226 + x228 + x231 +
                      x233 + x235 + x241 + x256 + x260 * x38 + x263 * x55);  // xxxyz
  result[5] = x * (-x112 * (x178 * x2 + x181 * x2 + x182 * x2 - x184 * x2 + x190 -
                            x2 * x265 + x267 + x270) -
                   x160 * x271 + x172 * x2 + x192 * x2 + x195 * x272 - x197 * x273 -
                   x198 * x2 - x199 * x2 + x2 * x201 + x2 * x220 + x2 * x223 + x2 * x226 +
                   x2 * x228 + x2 * x231 + x2 * x233 + x2 * x235 + x2 * x256 +
                   x219 * x284 + x255 + x275 * x36 + x286 + x289);  // xxxzz
  result[6] =
        y * (x0 * x327 + x173 * x330 + x242 - x249 * x302 + x253 -
             x301 * (x20 * x291 - x20 * x294 - x291 - 5.0 * x292 + x294 + 11.0 * x295 +
                     51.0 * x296 - 16.0 * x297 + x300) +
             x321 * x322 + x324 * x38 + 9.0 * x326 + x329 * x38 + x333 + x344);  // xxyyy
  result[7] = z * (x1 * x355 -
                   x112 * (x188 - 15.0 * x292 + 33.0 * x295 + 153.0 * x296 - 48.0 * x297 +
                           x299 + x345) -
                   x162 * x249 - x193 * x349 - x196 * x351 + x205 * x38 + x217 * x352 +
                   x236 * x354 + x243 - x249 * x254 + x322 * x347 + 3.0 * x326 + x344 +
                   x348 * x38 + x353 * x38 + x356);  // xxyyz
  result[8] = y * (x0 * x259 * x379 -
                   x112 * (x267 + x299 - 15.0 * x357 + 33.0 * x358 + 153.0 * x359 -
                           48.0 * x360 + x364) -
                   x162 * x288 + x2 * x355 - x254 * x288 - x271 * x349 - x273 * x351 +
                   x275 * x38 + x284 * x352 + x286 + x322 * x375 + x356 + x376 * x38 +
                   3.0 * x377 + x378 * x38 + x383);  // xxyzz
  result[9] =
        z * (x0 * x391 + x268 * x330 + x285 - x288 * x302 + x289 -
             x301 * (-x20 * x385 + x300 + 11.0 * x358 + 51.0 * x359 - 16.0 * x360 + x385 +
                     x387) +
             x322 * x388 + x333 + 9.0 * x377 + x38 * x390 + x38 * x393 + x383);  // xxzzz
  result[10] = x * (x1 * x327 - x1 * x394 - x1 * x418 +
                    x112 * (x1 * x400 + 180.0 * x292 + 24.0 * x314 - 6.0 * x316 -
                            x317 * x403 + x398 - x399 - x401 * x402 + x404 + x415) +
                    x173 * x423 - x217 * x417 - x234 * x312 + x321 * x420 + x395 * x396 +
                    x395 * x422 + x395 * x424 - x416 + x425 + x432);  // xyyyy
  result[11] = x266 * (x1 * x424 + x1 * x455 +
                       x112 * (-x434 - x435 - x437 - x438 - x440 + x441 + x442 + x444 +
                               x445 + x447 + x450) +
                       x218 * x347 + x259 * x454 + x321 * x57 + x327 + x433 + x451 +
                       x453 + x456);  // xyyyz
  result[12] = x * (-x1 * x240 +
                    x112 * (-x2 * x435 - x2 * x437 - x2 * x438 - x2 * x440 + x2 * x442 +
                            x2 * x444 + x2 * x445 + x2 * x447 + x259 + x311 + x313 +
                            x315 - x316 - x318 - x319 + x361 - x363 + x365 + x367 + x369 -
                            x370 - x372 - x373 - x469 + x470 + x471) +
                    x114 * x125 * x459 + x123 + x2 * x433 - x205 + x217 * x464 + x23 +
                    x254 * (x125 - x2 * x311 - x2 * x313 - x2 * x315 + x2 * x316 +
                            x2 * x318 + x2 * x319 + x208 + x209 - x210 - x212 + x276 +
                            x277 - x278 - x279 + x280 - x281 - x466 + x467 + x468) -
                    x275 + x284 * x454 - x335 + x341 * x461 + x347 * x465 - x348 - x353 +
                    x354 * x458 + x375 * x419 - x376 - x378 - x380 + x424 * x458 - x457 +
                    x463);  // xyyzz
  result[13] =
        x266 * (x112 * (-x2 * x406 - x2 * x409 - x2 * x436 + x2 * x443 + x2 * x474 +
                        x20 * x472 + x20 * x473 + x269 * x41 + x450 - x472 - 5.0 * x473) +
                x2 * x396 + x2 * x424 + x2 * x455 + x204 * x389 + x218 * x375 +
                x259 * x464 + x274 * x452 + x388 * x57 + x391 + x456);  // xyzzz
  result[14] = x * x482;                                                // xzzzz
  result[15] =
        y * (x112 * (x109 - x20 * x490 + 300.0 * x292 + 40.0 * x314 - 10.0 * x316 -
                     x317 * x99 + 180.0 * x317 + x415 + x490 - 150.0 * x491) +
             x18 * x484 + x203 * x485 * x486 + x23 * (-70.0 * x196 + 63.0 * x483 + 15.0) +
             x323 * x487 * x58 + x430 * x488);  // yyyyy
  result[16] =
        z * (x1 * x128 * x485 + 18.0 * x1 * x462 -
             x112 * (-x164 * x493 - x165 * x412 + x168 + x20 * x313 - 90.0 * x20 * x413 -
                     162.0 * x292 - 6.0 * x314 - 72.0 * x317 - x398 + x399 + x405 + x410 -
                     x411 + 22.0 * x412 + 195.0 * x413 + 72.0 * x491 + 105.0 * x493) +
             x116 * x395 + x123 * x484 + x145 * x431 + x162 * x430 + x203 * x487 * x502 +
             x23 * x484 + x416 * x487 - x417 * x487 * x496 - x420 * x501 + x432 +
             30.0 * x485 * x503 + 30.0 * x492 * (7.0 * x196 - 3.0));  // yyyyz
  result[17] =
        y *
        (x110 * x272 * x323 -
         x112 * (x167 + x177 * x458 + x270 + x306 - x314 - x41 * x467 + x469 - x470 -
                 x497 - x498 + x499 + x500 + 72.0 * x506 + 177.0 * x507 - 72.0 * x508) +
         x122 * x134 * x2 * x323 - x134 * x505 - x162 * x323 + x2 * x222 * x505 +
         x2 * x230 * x487 + 15.0 * x2 * x232 * x485 + x2 * x451 + x2 * x453 +
         x218 * x284 * x487 - x229 * x509 - x23 * x485 + 18.0 * x234 * x459 +
         x238 * x458 - x273 * x350 * x501 + x275 * x485 - x324 + x380 * x485 -
         x465 * x501 - 15.0 * x492 - x502 * x509 - 9.0 * x503 + 210.0 * x504 +
         10.5 * x510 * x511 + x516);  // yyyzz
  result[18] =
        z * (-x1 * x230 + x1 * x391 + 9.0 * x1 * x461 - x122 * x272 * x496 +
             x14 * x203 * x280 * x66 - x144 * x328 - 4.5 * x144 * x509 +
             x173 * x2 * x234 + x203 * x382 * x511 + x218 * x496 - x262 * x505 -
             x301 * (-x20 * x517 + x345 + x387 + x45 + 11.0 * x506 + 51.0 * x507 -
                     16.0 * x508 - x512 - x513 + x514 + x517) +
             x302 * x496 - x302 * x515 - x331 * x487 + x336 * x510 + x343 * x458 +
             x380 * x487 + x388 * x419 + x390 * x487 + x393 * x487 + x463 - x464 * x496 -
             90.0 * x492 - 54.0 * x503 + 315.0 * x504 + x516);  // yyzzz
  result[19] = x482 * y;                                        // yzzzz
  result[20] =
        z *
        (x112 * (x109 - 150.0 * x2 * x75 - x20 * x519 + 300.0 * x357 + 40.0 * x368 -
                 10.0 * x370 - x371 * x99 + 180.0 * x371 + x481 + x519) +
         x18 * (-30.0 * x273 + 35.0 * x518 + 3.0) +
         x23 * (-70.0 * x273 + 63.0 * x518 + 15.0) + x274 * x486 * (5.0 * x273 - 3.0) +
         x389 * x58 * (3.0 * x273 - 1.0) + x477 * x488);  // zzzzz
  return result;
}
Eigen::VectorXd T6_damp_thole(const Eigen::Vector3d& rij, double a) {
  double x = rij(0);
  double y = rij(1);
  double z = rij(2);
  Eigen::VectorXd result(28);
  double x0  = x * x;
  double x1  = y * y;
  double x2  = z * z;
  double x3  = x0 + x1 + x2;
  double x4  = 1.0 / x3;
  double x5  = x0 * x4;
  double x6  = x * x * x * x;
  double x7  = 1.0 / (x3 * x3);
  double x8  = x6 * x7;
  double x9  = -70.0 * x5 + 63.0 * x8 + 15.0;
  double x10 = 1.0 / (x3 * x3 * x3 * x3);
  double x11 = a * x10;
  double x12 = x0 * x11;
  double x13 = sqrt(x3);
  double x14 = a * x13;
  double x15 = exp(-x14);
  double x16 = x14 + 1.0;
  double x17 = x15 * x16;
  double x18 = x12 * x17;
  double x19 = 45.0 * x18;
  double x20 = pow(x, 6);
  double x21 = 1.0 / (x3 * x3 * x3);
  double x22 = 231.0 * x21;
  double x23 = pow(x3, -3.5);
  double x24 = x14 + 2.0;
  double x25 = x15 * x24 - 2.0;
  double x26 = x23 * x25;
  double x27 = 22.5 * x26;
  double x28 = -30.0 * x5 + 35.0 * x8 + 3.0;
  double x29 = pow(x3, -2.5);
  double x30 = a * x29;
  double x31 = pow(x3, -1.5);
  double x32 = x0 * x31;
  double x33 = x24 * x32;
  double x34 = a * x4;
  double x35 = 2.0 * x34;
  double x36 = x24 * x34;
  double x37 = 1.0 / x13;
  double x38 = -x24 * x37 + x37;
  double x39 = -x0 * x35 + x0 * x36 - x32 + x33 + x38;
  double x40 = x15 * x39;
  double x41 = 22.5 * x40;
  double x42 = x15 * x30;
  double x43 = a * x7;
  double x44 = 6.0 * x43;
  double x45 = a * a;
  double x46 = 3.0 * x31;
  double x47 = x45 * x46;
  double x48 = x33 * x45;
  double x49 = 3.0 * x43;
  double x50 = x24 * x49;
  double x51 = 3.0 * x29;
  double x52 = x0 * x51;
  double x53 = x24 * x52;
  double x54 = -x52 + x53;
  double x55 = x24 * x46;
  double x56 = 6.0 * x34;
  double x57 = 3.0 * x36;
  double x58 = x46 - x55 + x56 - x57;
  double x59 = -x0 * x44 - x0 * x47 + x0 * x50 + x48 + x54 + x58;
  double x60 = 30.0 * x59;
  double x61 = 5.0 * x5 - 3.0;
  double x62 = x0 * x61;
  double x63 = 3.0 * x5 - 1.0;
  double x64 = 15.0 * x23;
  double x65 = x6 * x64;
  double x66 = 18.0 * x29;
  double x67 = x0 * x66;
  double x68 = a * x21;
  double x69 = 30.0 * x68;
  double x70 = x45 * x6;
  double x71 = x24 * x67;
  double x72 = a * a * a;
  double x73 = x72 * x8;
  double x74 = x32 * x45;
  double x75 = x0 * x43;
  double x76 = 18.0 * x43;
  double x77 = x0 * x24;
  double x78 = x29 * x45;
  double x79 = x24 * x78;
  double x80 = 6.0 * x79;
  double x81 = 15.0 * x68;
  double x82 = x24 * x6;
  double x83 = -x46 + x55 - x56 + x57;
  double x84 = x24 * x65 + x24 * x73 - 6.0 * x48 - x6 * x69 + x6 * x80 - x65 - x66 * x70 +
               x67 - x71 - 4.0 * x73 + 18.0 * x74 + 36.0 * x75 - x76 * x77 + x81 * x82 +
               x83;
  double x85  = a * x31;
  double x86  = x15 * x85;
  double x87  = 7.5 * x86;
  double x88  = x0 * x23;
  double x89  = 150.0 * x88;
  double x90  = x11 * x6;
  double x91  = x23 * x70;
  double x92  = x21 * x72;
  double x93  = x6 * x92;
  double x94  = a * a * a * a;
  double x95  = x29 * x94;
  double x96  = x6 * x95;
  double x97  = x7 * x72;
  double x98  = x0 * x97;
  double x99  = 180.0 * x78;
  double x100 = x0 * x68;
  double x101 = x24 * x96;
  double x102 = 150.0 * x24;
  double x103 = x0 * x78;
  double x104 = 60.0 * x24;
  double x105 = x24 * x98;
  double x106 = x24 * x93;
  double x107 = 45.0 * x23;
  double x108 = x107 * x24;
  double x109 = x24 * x90;
  double x110 = pow(x3, -4.5);
  double x111 = x110 * x6;
  double x112 = 105.0 * x111;
  double x113 = x112 * x24;
  double x114 = -x112 + x113;
  double x115 = 45.0 * x29;
  double x116 = 90.0 * x43;
  double x117 = x31 * x45;
  double x118 = 45.0 * x117;
  double x119 = x115 * x24;
  double x120 = x24 * x31;
  double x121 = x120 * x45;
  double x122 = 15.0 * x121;
  double x123 = x24 * x43;
  double x124 = 45.0 * x123;
  double x125 = -x115 - x116 - x118 + x119 + x122 + x124;
  double x126 = x0 * x99 - x100 * x102 + 300.0 * x100 + x101 - x103 * x104 - 10.0 * x105 +
                10.0 * x106 + x108 * x70 + 105.0 * x109 + x114 + x125 - x24 * x89 + x89 -
                210.0 * x90 - 135.0 * x91 - 40.0 * x93 - 5.0 * x96 + 40.0 * x98;
  double x127 = a * x46;
  double x128 = x127 * x15;
  double x129 = pow(x3, -5.5);
  double x130 = 945.0 * x129;
  double x131 = x130 * x20;
  double x132 = 675.0 * x88;
  double x133 = 1575.0 * x111;
  double x134 = a * 1.0 / (x3 * x3 * x3 * x3 * x3);
  double x135 = 1890.0 * x134;
  double x136 = 1260.0 * x45;
  double x137 = x110 * x136;
  double x138 = x10 * x72;
  double x139 = 420.0 * x20;
  double x140 = x23 * x94;
  double x141 = 75.0 * x140;
  double x142 = a * a * a * a * a;
  double x143 = x142 * x21;
  double x144 = x143 * x20;
  double x145 = 675.0 * x24;
  double x146 = x24 * x64;
  double x147 = x146 * x94;
  double x148 = x20 * x24;
  double x149 = 105.0 * x138;
  double x150 = 270.0 * x78;
  double x151 = x24 * x45;
  double x152 = x110 * x151;
  double x153 = 945.0 * x134;
  double x154 = x115 + x116 + x118 - x119 - x122 - x124;
  double x155 = 0.5 * x15;
  double x156 = a * x37;
  double x157 = x155 * x156;
  double x158 = 7.0 * x5;
  double x159 = x158 - 3.0;
  double x160 = x134 * x17;
  double x161 = x0 * x160;
  double x162 = x0 * x25;
  double x163 = x129 * x162;
  double x164 = 210.0 * x163;
  double x165 = x11 * x17;
  double x166 = x165 * x28;
  double x167 = 7.5 * x165;
  double x168 = x110 * x25;
  double x169 = 52.5 * x168;
  double x170 = 7.5 * x42;
  double x171 = x120 - x31;
  double x172 = x171 - x35 + x36;
  double x173 = x172 * x28;
  double x174 = x0 * x110;
  double x175 = a * x40;
  double x176 = x174 * x175;
  double x177 = a * x23;
  double x178 = x177 * x40;
  double x179 = x21 * x45;
  double x180 = x179 * x40;
  double x181 = 5.0 * x75;
  double x182 = 2.0 * x123;
  double x183 = x0 * x182;
  double x184 = x171 - x34;
  double x185 = x181 - x183 + x184 + x52 - x53 + x74;
  double x186 = 15.0 * x42;
  double x187 = a * x15;
  double x188 = x187 * x88;
  double x189 = x188 * x60;
  double x190 = x59 * x63;
  double x191 = x45 * x7;
  double x192 = x15 * x191;
  double x193 = 5.0 * x192;
  double x194 = 12.0 * x103;
  double x195 = 27.0 * x100;
  double x196 = x100 * x24;
  double x197 = 12.0 * x196;
  double x198 = x45 * x53;
  double x199 = x0 * x64;
  double x200 = x199 * x24;
  double x201 = x199 - x200;
  double x202 = -x47;
  double x203 = 15.0 * x43;
  double x204 = x24 * x44;
  double x205 = 9.0 * x29;
  double x206 = x205 * x24;
  double x207 = -x205 + x206;
  double x208 = x202 - x203 + x204 + x207;
  double x209 = x194 + x195 - x197 - x198 + x201 + x208 + x98;
  double x210 = 5.0 * x86;
  double x211 = x15 * x84;
  double x212 = x211 * x30;
  double x213 = x191 * x211;
  double x214 = 90.0 * x88;
  double x215 = 162.0 * x100;
  double x216 = 72.0 * x103;
  double x217 = 6.0 * x98;
  double x218 = 22.0 * x93;
  double x219 = x214 * x24;
  double x220 = 105.0 * x91;
  double x221 = 195.0 * x90;
  double x222 = 90.0 * x109;
  double x223 = 30.0 * x24;
  double x224 = x223 * x91;
  double x225 = 4.0 * x106;
  double x226 = x45 * x71;
  double x227 = 72.0 * x196;
  double x228 = x205 - x206;
  double x229 = x203 - x204 + x47;
  double x230 = x228 + x229;
  double x231 = x15 * (x112 - x113 - x214 - x215 - x216 - x217 + x218 + x219 + x220 +
                       x221 - x222 - x224 - x225 + x226 + x227 + x230 + x96);
  double x232 = x231 * x85;
  double x233 = 0.5 * x86;
  double x234 = x4 * x45;
  double x235 = x155 * x234;
  double x236 = 1050.0 * x174;
  double x237 = x130 * x6;
  double x238 = 1050.0 * x45;
  double x239 = x237 * x24;
  double x240 = x0 * x92;
  double x241 = 10.0 * x95;
  double x242 = x140 * x6;
  double x243 = x138 * x6;
  double x244 = x134 * x6;
  double x245 = x24 * x244;
  double x246 = x111 * x45;
  double x247 = x24 * x246;
  double x248 = 5.0 * x23;
  double x249 = x248 * x94;
  double x250 = x45 * x88;
  double x251 = 300.0 * x24;
  double x252 = x12 * x24;
  double x253 = 225.0 * x23;
  double x254 = x24 * x253;
  double x255 = 15.0 * x97;
  double x256 = x24 * x68;
  double x257 = -x119 * x45 + x253 - x254 + x255 - 180.0 * x256 + 405.0 * x68 + x99;
  double x258 =
        x * (x126 * x233 + x126 * x235 -
             x157 * (-x0 * x241 - x104 * x243 + x111 * x238 - 1950.0 * x12 + x143 * x6 +
                     x236 * x24 - x236 + x237 - x238 * x88 - x239 + 40.0 * x24 * x240 -
                     220.0 * x240 + 35.0 * x242 + 285.0 * x243 + 1785.0 * x244 -
                     840.0 * x245 - 315.0 * x247 - x249 * x82 + x250 * x251 +
                     900.0 * x252 + x257) +
             150.0 * x159 * x161 + x164 * (9.0 * x5 - 5.0) + 37.5 * x166 + x167 * x9 +
             x169 * x9 + x170 * x173 + 150.0 * x176 + 75.0 * x178 * x61 +
             15.0 * x180 * x61 - x185 * x186 * x61 + x186 * x190 + x189 + x190 * x193 -
             x209 * x210 * x63 + 7.5 * x212 + 2.5 * x213 - 2.5 * x232);
  double x259 = 630.0 * x174;
  double x260 = x1 * x259;
  double x261 = 15.0 * x78;
  double x262 = x1 * x261;
  double x263 = 69.0 * x68;
  double x264 = x1 * x263;
  double x265 = 1062.0 * x12;
  double x266 = 432.0 * x250;
  double x267 = 42.0 * x240;
  double x268 = 24.0 * x256;
  double x269 = x1 * x268;
  double x270 = x1 * x23;
  double x271 = x270 * x94;
  double x272 = 9.0 * x6;
  double x273 = 162.0 * x243;
  double x274 = 825.0 * x246;
  double x275 = 1665.0 * x244;
  double x276 = 720.0 * x245;
  double x277 = 210.0 * x247;
  double x278 = x24 * x243;
  double x279 = 24.0 * x278;
  double x280 = x219 * x45;
  double x281 = 432.0 * x252;
  double x282 = x1 * x107;
  double x283 = x24 * x282;
  double x284 = x282 - x283;
  double x285 = x114 + x208 + x214 + x215 + x216 + x217 - x218 - x219 - x220 - x221 +
                x222 + x224 + x225 - x226 - x227 - x96;
  double x286 = 12.0 * x4;
  double x287 = x0 * x1;
  double x288 = 42.0 * x7;
  double x289 = 3.0 - x158;
  double x290 = x174 * x25;
  double x291 = 30.0 * x290;
  double x292 = x1 * x110;
  double x293 = x25 * x292;
  double x294 = 52.5 * x28;
  double x295 = x1 * x31;
  double x296 = a * x231;
  double x297 = x1 * x4;
  double x298 = x231 * x45;
  double x299 = x129 * x25;
  double x300 = 300.0 * x159;
  double x301 = x299 * x300;
  double x302 = 7.0 * x240;
  double x303 = 72.0 * x250;
  double x304 = 177.0 * x12;
  double x305 = 72.0 * x252;
  double x306 = x1 * x200;
  double x307 = 105.0 * x174;
  double x308 = x1 * x307;
  double x309 = x24 * x308;
  double x310 = x308 - x309;
  double x311 = -x199 + x200;
  double x312 = -x194 - x195 + x197 + x198 + x230 + x311 - x98;
  double x313 = x1 * x302 + x1 * x303 + x1 * x304 - x1 * x305 - x262 - x264 + x269 -
                x282 + x283 - x306 * x45 + x310 + x312;
  double x314 = x187 * x32;
  double x315 = 2.0 * x314;
  double x316 = x1 * x49;
  double x317 = 5.0 * x78;
  double x318 = 23.0 * x68;
  double x319 = 8.0 * x196;
  double x320 = -x120 + x31;
  double x321 = x1 * x51;
  double x322 = x24 * x321;
  double x323 = -x321 + x322;
  double x324 = x320 + x323 + x34;
  double x325 = -x181 + x183 + x54 - x74;
  double x326 =
        x1 * x199 - x1 * x319 + x287 * x317 + x287 * x318 - x306 - x316 + x324 + x325;
  double x327 = x15 * x326;
  double x328 = x127 * x327;
  double x329 = 1.5 * x212;
  double x330 = x1 * x120;
  double x331 = -x1 * x35 + x1 * x36 - x295 + x330 + x38;
  double x332 = 1.5 * x15;
  double x333 = x30 * x332;
  double x334 = x28 * x333;
  double x335 = 1.5 * x213;
  double x336 = x295 * x72;
  double x337 = x155 * x84;
  double x338 = x16 * x187 / pow(x3, 6);
  double x339 = 840.0 * x338;
  double x340 = x339 * x6;
  double x341 = 180.0 * x176;
  double x342 = x172 * x187;
  double x343 = 120.0 * x342;
  double x344 = x111 * x343;
  double x345 = x175 * x63;
  double x346 = x10 * x40;
  double x347 = x346 * x45;
  double x348 = 36.0 * x347;
  double x349 = 21.0 * x180 * x63;
  double x350 = x1 * x11;
  double x351 = x17 * x350;
  double x352 = 15.0 * x28;
  double x353 = 14.0 * x179;
  double x354 = x15 * x59;
  double x355 = x287 * x354;
  double x356 = x320 + x35 - x36;
  double x357 = x1 * x44;
  double x358 = x1 * x47;
  double x359 = x330 * x45;
  double x360 = x24 * x316;
  double x361 = x323 - x357 - x358 + x359 + x360;
  double x362 = x356 + x361;
  double x363 = 6.0 * x42;
  double x364 = x363 * x62;
  double x365 = x40 * x63;
  double x366 = x365 * x72;
  double x367 = x29 * x72;
  double x368 = 2.0 * x367;
  double x369 = 4.0 * x192;
  double x370 = x209 * x287;
  double x371 = x185 * x192;
  double x372 = 6.0 * x371 * x63;
  double x373 = 12.0 * x42;
  double x374 = x1 * x66;
  double x375 = x185 * x187;
  double x376 = x375 * x63;
  double x377 = 36.0 * x185 * x188;
  double x378 = 210.0 * x61;
  double x379 = x134 * x287;
  double x380 = x17 * x379;
  double x381 = 60.0 * x17;
  double x382 = x159 * x381;
  double x383 = x172 * x188;
  double x384 = 60.0 * x383 * x61;
  double x385 = 7.5 * x26;
  double x386 = x0 * x59;
  double x387 = 2.0 * x192;
  double x388 = x205 * x40;
  double x389 = 18.0 * x88;
  double x390 = a * x388 * x63 + x175 * x389 + 30.0 * x18 * x61 + 3.0 * x191 * x365 +
                x233 * x84 + x235 * x84 + x244 * x381 + x28 * x385 + x363 * x386 +
                x386 * x387;
  double x391 = x24 * x51;
  double x392 = x391 - x51;
  double x393 = x121 + x202 + x392 - x44 + x50;
  double x394 = x0 * x317;
  double x395 = x0 * x318;
  double x396 = x392 - x49;
  double x397 = x201 - x319 + x394 + x395 + x396;
  double x398 = x15 * x386;
  double x399 = x0 * x209;
  double x400 = x24 * x307;
  double x401 = x200 * x45;
  double x402 = -x107 + x108;
  double x403 = x15 * (-x261 - x263 + x268 + x302 + x303 + x304 - x305 + x307 - x400 -
                       x401 + x402);
  double x404 = a * x403;
  double x405 = x31 * x72;
  double x406 = x155 * x405;
  double x407 = 210.0 * x174;
  double x408 = 315.0 * x129;
  double x409 = x408 * x6;
  double x410 = 8.0 * x256;
  double x411 = x24 * x407;
  double x412 = -x146 + x64;
  double x413 = x156 * x332;
  double x414 = y * z;
  double x415 = x2 * x259;
  double x416 = x2 * x23;
  double x417 = x416 * x94;
  double x418 = x2 * x261;
  double x419 = x2 * x263;
  double x420 = x2 * x268;
  double x421 = x107 * x2;
  double x422 = x108 * x2;
  double x423 = x421 - x422;
  double x424 = x418 + x419 - x420 + x423;
  double x425 = -x2 * x286;
  double x426 = x0 * x2;
  double x427 = x110 * x2;
  double x428 = x25 * x427;
  double x429 = x2 * x31;
  double x430 = x2 * x4;
  double x431 = x2 * x307;
  double x432 = x2 * x400;
  double x433 = x431 - x432;
  double x434 = -x418 - x419 + x420 - x421 + x422;
  double x435 =
        x2 * x302 + x2 * x303 + x2 * x304 - x2 * x305 - x2 * x401 + x312 + x433 + x434;
  double x436 = x2 * x51;
  double x437 = x2 * x391;
  double x438 = -x436 + x437;
  double x439 = -x2 * x49 + x438;
  double x440 = x199 * x2 - x2 * x200 - x2 * x319 + x2 * x394 + x2 * x395 + x320 + x325 +
                x34 + x439;
  double x441 = x15 * x440;
  double x442 = x127 * x441;
  double x443 = x120 * x2 - x2 * x35 + x2 * x36 + x38 - x429;
  double x444 = x429 * x72;
  double x445 = x11 * x2;
  double x446 = x17 * x445;
  double x447 = x354 * x426;
  double x448 = x2 * x44;
  double x449 = x2 * x47;
  double x450 = x121 * x2;
  double x451 = x2 * x50;
  double x452 = x438 - x448 - x449 + x450 + x451;
  double x453 = x356 + x452;
  double x454 = x209 * x426;
  double x455 = x2 * x66;
  double x456 = x134 * x426;
  double x457 = x17 * x456;
  double x458 = pow(x3, -6.5);
  double x459 = x25 * x458;
  double x460 = 2835.0 * x459;
  double x461 = x287 * x460;
  double x462 = 13.5 * x42;
  double x463 = x1 * x248;
  double x464 = 5.0 * x68;
  double x465 = 35.0 * x174;
  double x466 = x1 * x465;
  double x467 = 11.0 * x250;
  double x468 = 51.0 * x12;
  double x469 = 16.0 * x252;
  double x470 = -x391 + x51;
  double x471 = x311 + x319 - x394 - x395 + x49;
  double x472 = x470 + x471;
  double x473 = -x1 * x464 + x1 * x467 + x1 * x468 - x1 * x469 + x24 * x463 - x24 * x466 -
                x463 + x466 + x472;
  double x474 = 4.5 * x86;
  double x475 = 4.5 * x192;
  double x476 = -x326 * x475;
  double x477 = 1.5 * x86;
  double x478 = x270 * x45;
  double x479 = x138 * x287;
  double x480 = x24 * x350;
  double x481 = x174 * x45;
  double x482 = x1 * x481;
  double x483 = x24 * x379;
  double x484 = x287 * x408;
  double x485 = -x307 + x400;
  double x486 = -x24 * x484 + x484 + x485;
  double x487 = -x268;
  double x488 = 105.0 * x292;
  double x489 = x24 * x488;
  double x490 = x107 - x108;
  double x491 = -x488 + x489 + x490;
  double x492 = x261 + x263 + x487 + x491;
  double x493 = -x302 - x303 - x304 + x305 + x401;
  double x494 = x234 * x332;
  double x495 = x1 * x299;
  double x496 = 472.5 * x495;
  double x497 = x496 * x61;
  double x498 = x1 * x59;
  double x499 = x7 * x94;
  double x500 = x155 * x499;
  double x501 = x498 * x500;
  double x502 = x15 * x72;
  double x503 = x502 * x59;
  double x504 = x321 * x503;
  double x505 = x40 * x72;
  double x506 = 18.0 * x270;
  double x507 = x505 * x506;
  double x508 = x187 * x362;
  double x509 = x508 * x88;
  double x510 = x187 * x331;
  double x511 = x174 * x510;
  double x512 = 67.5 * x375;
  double x513 = x270 * x512;
  double x514 = 31.5 * x179;
  double x515 = x15 * x185;
  double x516 = x1 * x515;
  double x517 = x514 * x516;
  double x518 = 4.5 * x42;
  double x519 = x1 * x209;
  double x520 = x518 * x519;
  double x521 = x475 * x519;
  double x522 = 4.5 * x367;
  double x523 = x516 * x522;
  double x524 = x209 * x332;
  double x525 = x336 * x524;
  double x526 = x361 + x58;
  double x527 = x333 * x526;
  double x528 = x1 * x64;
  double x529 = x1 * x146;
  double x530 = -x528 + x529;
  double x531 = x228 + x530;
  double x532 = 9.0 * x117;
  double x533 = 9.0 * x123;
  double x534 = x45 * x55;
  double x535 = x532 - x533 - x534 + x76;
  double x536 = x1 * x69;
  double x537 = x374 * x45;
  double x538 = x1 * x97;
  double x539 = 4.0 * x538;
  double x540 = x24 * x538;
  double x541 = x1 * x80;
  double x542 = x1 * x81;
  double x543 = x24 * x542;
  double x544 = -x536 - x537 - x539 + x540 + x541 + x543;
  double x545 = x531 + x535 + x544;
  double x546 = x477 * x63;
  double x547 = x332 * x39;
  double x548 = x21 * x94;
  double x549 = x547 * x548;
  double x550 = x1 * x549;
  double x551 = 7.5 * x59;
  double x552 = x187 * x270;
  double x553 = x551 * x552;
  double x554 = x15 * x179;
  double x555 = 7.5 * x554;
  double x556 = x498 * x555;
  double x557 = x462 * x63;
  double x558 = x15 * x331;
  double x559 = 22.5 * x177;
  double x560 = x558 * x559;
  double x561 = 85.5 * x347;
  double x562 = x1 * x561;
  double x563 = 157.5 * x175;
  double x564 = x292 * x563;
  double x565 = x174 * x342;
  double x566 = x1 * x565;
  double x567 = 1575.0 * x338;
  double x568 = x287 * x567;
  double x569 = x270 * x342;
  double x570 = 67.5 * x569 * x63;
  double x571 = x1 * x160;
  double x572 = 157.5 * x571;
  double x573 = x572 * x63;
  double x574 = x572 * x61;
  double x575 = 157.5 * x168;
  double x576 = 4.5 * x59;
  double x577 = 4.5 * x40;
  double x578 = x332 * x59;
  double x579 = x165 * x63;
  double x580 = 22.5 * x165;
  double x581 = -315.0 * x161 - 630.0 * x163 - x172 * x557 - 67.5 * x178 - 31.5 * x180 +
                x185 * x462 + x185 * x475 - x192 * x576 + x209 * x477 + x209 * x494 -
                x367 * x577 - 27.0 * x383 - x405 * x578 - x42 * x576 - x575 * x61 -
                67.5 * x579 - x580 * x61;
  double x582 = x * y;
  double x583 = 315.0 * x292;
  double x584 = x24 * x583;
  double x585 = x130 * x287;
  double x586 = x24 * x585;
  double x587 = x261 + x263 + x487 + x493;
  double x588 = x485 + x490 + x587;
  double x589 = 33.0 * x250;
  double x590 = x1 * x589;
  double x591 = 153.0 * x12;
  double x592 = x1 * x591;
  double x593 = 48.0 * x252;
  double x594 = x1 * x593;
  double x595 = x470 + x530;
  double x596 = x310 + x471 - x542 + x590 + x592 - x594 + x595;
  double x597 = x477 * x596;
  double x598 = x191 * x332;
  double x599 = x326 * x598;
  double x600 = x326 * x518;
  double x601 = x403 * x45;
  double x602 = x1 * x397;
  double x603 = 3.0 * x192;
  double x604 = x187 * x205;
  double x605 = x177 * x558;
  double x606 = 7.5 * x61;
  double x607 = x362 * x518;
  double x608 = x607 * x63;
  double x609 = x333 * x61;
  double x610 = -x121 + x44 + x47 - x50;
  double x611 = x544 + x595 + x610;
  double x612 = x477 * x611;
  double x613 = x188 * x393;
  double x614 = 18.0 * x613;
  double x615 = x342 * x61;
  double x616 = x393 * x604;
  double x617 = x1 * x616;
  double x618 = 22.5 * x178;
  double x619 = x172 * x518;
  double x620 = x619 * x63;
  double x621 = 105.0 * x161 + x164 + x167 * x61 + x169 * x61 + 10.5 * x180 -
                x185 * x518 - x185 * x598 + x191 * x578 - x209 * x233 - x209 * x235 +
                x30 * x578 + x367 * x547 + 9.0 * x383 + x406 * x59 + 22.5 * x579 + x618 +
                x620;
  double x622 = x * z;
  double x623 = 315.0 * x427;
  double x624 = x416 * x45;
  double x625 = x24 * x623;
  double x626 = x130 * x426;
  double x627 = x138 * x426;
  double x628 = x24 * x445;
  double x629 = x2 * x481;
  double x630 = x24 * x456;
  double x631 = x426 * x460;
  double x632 = x2 * x299;
  double x633 = 472.5 * x632;
  double x634 = x61 * x633;
  double x635 = x2 * x81;
  double x636 = x2 * x589;
  double x637 = x2 * x591;
  double x638 = x2 * x593;
  double x639 = x2 * x64;
  double x640 = -x639;
  double x641 = x146 * x2;
  double x642 = x470 + x640 + x641;
  double x643 = x433 + x471 - x635 + x636 + x637 - x638 + x642;
  double x644 = x477 * x643;
  double x645 = x440 * x598;
  double x646 = x440 * x518;
  double x647 = 18.0 * x416;
  double x648 = x505 * x647;
  double x649 = x187 * x443;
  double x650 = x174 * x649;
  double x651 = x187 * x453;
  double x652 = x651 * x88;
  double x653 = x436 * x503;
  double x654 = x2 * x397;
  double x655 = x427 * x563;
  double x656 = x2 * x561;
  double x657 = x187 * x416;
  double x658 = x551 * x657;
  double x659 = x15 * x443;
  double x660 = x2 * x551 * x554;
  double x661 = x453 * x518;
  double x662 = x63 * x661;
  double x663 = x2 * x69;
  double x664 = x45 * x455;
  double x665 = x2 * x97;
  double x666 = 4.0 * x665;
  double x667 = x24 * x665;
  double x668 = x2 * x80;
  double x669 = x24 * x635;
  double x670 = -x663 - x664 - x666 + x667 + x668 + x669;
  double x671 = x610 + x642 + x670;
  double x672 = x477 * x671;
  double x673 = x2 * x549;
  double x674 = x2 * x500;
  double x675 = x59 * x674;
  double x676 = x444 * x524;
  double x677 = x2 * x209;
  double x678 = x518 * x677;
  double x679 = x475 * x677;
  double x680 = x2 * x515;
  double x681 = x522 * x680;
  double x682 = x514 * x680;
  double x683 = x416 * x512;
  double x684 = x426 * x567;
  double x685 = x2 * x565;
  double x686 = x2 * x616;
  double x687 = x160 * x2;
  double x688 = 157.5 * x687;
  double x689 = x63 * x688;
  double x690 = x61 * x688;
  double x691 = 67.5 * x416;
  double x692 = x342 * x691;
  double x693 = x63 * x692;
  double x694 = x2 * x465;
  double x695 = x2 * x248;
  double x696 = -x2 * x464 + x24 * x695 - x695;
  double x697 = x2 * x467 + x2 * x468 - x2 * x469 - x24 * x694 + x472 + x694 + x696;
  double x698 = -x440 * x475;
  double x699 = x408 * x426;
  double x700 = -x24 * x699 + x485 + x699;
  double x701 = 105.0 * x427;
  double x702 = x24 * x701;
  double x703 = -x701 + x702;
  double x704 = x490 + x703;
  double x705 = -153.0 * x445 - 33.0 * x624 + 48.0 * x628;
  double x706 = x452 + x58;
  double x707 = x333 * x706;
  double x708 = x228 + x535 + x640 + x641 + x670;
  double x709 = x559 * x659;
  double x710 = y * y * y * y;
  double x711 = 3780.0 * x162 * x458;
  double x712 = 3.0 * x234;
  double x713 = x25 * x63;
  double x714 = 1890.0 * x299;
  double x715 = 472.5 * x299;
  double x716 = x63 * x715;
  double x717 = x110 * x710;
  double x718 = 35.0 * x717;
  double x719 = 30.0 * x270;
  double x720 = x11 * x710;
  double x721 = x0 * x408;
  double x722 = x710 * x721;
  double x723 = 306.0 * x12;
  double x724 = 66.0 * x250;
  double x725 = x45 * x710;
  double x726 = 93.0 * x174;
  double x727 = x134 * x710;
  double x728 = 443.0 * x0;
  double x729 = 128.0 * x77;
  double x730 = 96.0 * x252;
  double x731 = 108.0 * x78;
  double x732 = 180.0 * x68;
  double x733 = 90.0 * x68;
  double x734 = x1 * x733;
  double x735 = 36.0 * x79;
  double x736 = 105.0 * x717;
  double x737 = -x736;
  double x738 = 90.0 * x270;
  double x739 = x24 * x738;
  double x740 = x24 * x736;
  double x741 = x737 + x738 - x739 + x740;
  double x742 = x207 - x532 + x533 + x534 - x76;
  double x743 = 135.0 * x23;
  double x744 = 40.0 * x92;
  double x745 = x710 * x95;
  double x746 = x24 * x745;
  double x747 = x710 * x92;
  double x748 = 10.0 * x24;
  double x749 = 105.0 * x24;
  double x750 = x108 * x725 - x710 * x744 + x720 * x749 - 210.0 * x720 - x725 * x743 -
                5.0 * x745 + x746 + x747 * x748;
  double x751 = x1 * x731 + x1 * x732 - x1 * x735 - x24 * x734 + 24.0 * x538 -
                6.0 * x540 + x741 + x742 + x750;
  double x752 = 18.0 * x45;
  double x753 = x7 * x752;
  double x754 = x185 * x502;
  double x755 = x295 * x754;
  double x756 = x40 * x548;
  double x757 = 5.0 * x756;
  double x758 = x7 * x710;
  double x759 = 2.0 * x94;
  double x760 = x515 * x759;
  double x761 = x1 * x327;
  double x762 = 3.0 * x1;
  double x763 = x40 * x499;
  double x764 = x187 * x295;
  double x765 = 6.0 * x473;
  double x766 = x15 * x45;
  double x767 = x297 * x766;
  double x768 = x331 * x63;
  double x769 = 9.0 * x192;
  double x770 = x515 * x710;
  double x771 = 12.0 * x367;
  double x772 = 30.0 * x177;
  double x773 = 30.0 * x21;
  double x774 = x515 * x773;
  double x775 = x1 * x180;
  double x776 = 52.5 * x175;
  double x777 = 52.5 * x346;
  double x778 = x41 * x72;
  double x779 = x23 * x778;
  double x780 = x710 * x97;
  double x781 = x146 * x710;
  double x782 = x295 * x45;
  double x783 = 36.0 * x43;
  double x784 = x24 * x76;
  double x785 = x24 * x81;
  double x786 = x1 * x783 - x1 * x784 - x24 * x374 + x24 * x780 - 6.0 * x359 + x374 -
                x64 * x710 - x66 * x725 - x69 * x710 + x710 * x785 + x710 * x80 -
                4.0 * x780 + x781 + 18.0 * x782 + x83;
  double x787 = x233 * x786;
  double x788 = x142 * x155 * x29;
  double x789 = x39 * x788;
  double x790 = 2205.0 * x338;
  double x791 = x0 * x790;
  double x792 = 420.0 * x565;
  double x793 = 210.0 * x17;
  double x794 = x63 * x793;
  double x795 = 180.0 * x1;
  double x796 = x510 * x63;
  double x797 = x1 * x526;
  double x798 = 12.0 * x188;
  double x799 = x363 * x797;
  double x800 = x351 * x63;
  double x801 = -x19 - x191 * x577 - x27 * x63 - 90.0 * x290 - x30 * x577 - x405 * x547;
  double x802 = 3780.0 * x459;
  double x803 = x287 * x802;
  double x804 = 35.0 * x292;
  double x805 = x24 * x804;
  double x806 = x412 - x589 - x591 + x593 + x81;
  double x807 = x332 * x405;
  double x808 = x496 * x63;
  double x809 = 210.0 * x350;
  double x810 = 135.0 * x478;
  double x811 = x1 * x744;
  double x812 = x1 * x95;
  double x813 = 5.0 * x812;
  double x814 = x24 * x812;
  double x815 = x1 * x92;
  double x816 = x748 * x815;
  double x817 = x283 * x45;
  double x818 = 105.0 * x480;
  double x819 = 12.0 * x97;
  double x820 = 54.0 * x78;
  double x821 = 45.0 * x256;
  double x822 = x45 * x66;
  double x823 = x24 * x822;
  double x824 = x24 * x97;
  double x825 = 3.0 * x824;
  double x826 = x733 + x819 + x820 - x821 - x823 - x825;
  double x827 = x491 - x809 - x810 - x811 - x813 + x814 + x816 + x817 + x818 + x826;
  double x828 = x233 * x545;
  double x829 = x375 * x719;
  double x830 = 30.0 * x179;
  double x831 = 2.0 * x499;
  double x832 = 3.0 * x188;
  double x833 = x187 * x52;
  double x834 = 5.0 * x1;
  double x835 = x0 * x604;
  double x836 = x332 * x397;
  double x837 = x270 * x778;
  double x838 = 52.5 * x347;
  double x839 = x187 * x393;
  double x840 = x0 * x839;
  double x841 = 210.0 * x571;
  double x842 = x287 * x790;
  double x843 = 22.5 * x172;
  double x844 = x552 * x843;
  double x845 = 45.0 * x165;
  double x846 = -x130 * x162 - 495.0 * x161 - 22.5 * x180 + x185 * x604 + 9.0 * x371 -
                99.0 * x383 - x388 * x72 - x393 * x835 + x397 * x477 + x397 * x494 +
                x46 * x754 - x499 * x547 - x575 * x63 - x618 - x620 - x63 * x845;
  double x847 = x446 * x63;
  double x848 = 52.5 * x63;
  double x849 = x1 * x2;
  double x850 = x134 * x849;
  double x851 = 30.0 * x1;
  double x852 = 30.0 * x2;
  double x853 = 3.0 * x88;
  double x854 = 7.5 * x175;
  double x855 = x180 * x2;
  double x856 = x2 * x292;
  double x857 = x443 * x63;
  double x858 = x63 * x649;
  double x859 = 3.0 * x2;
  double x860 = x429 * x754;
  double x861 = x397 * x849;
  double x862 = x515 * x849;
  double x863 = x187 * x429;
  double x864 = x155 * x444;
  double x865 = x430 * x766;
  double x866 = x2 * x528;
  double x867 = x2 * x529;
  double x868 = x321 - x322;
  double x869 = x172 - x2 * x536 - x2 * x537 - x2 * x539 + x2 * x540 + x2 * x543 + x357 +
                x358 - x359 - x360 + x436 - x437 + x448 + x449 - x450 - x451 +
                x80 * x849 - x866 + x867 + x868;
  double x870 = x233 * x63;
  double x871 = x2 * x350;
  double x872 = x2 * x24;
  double x873 = x2 * x488;
  double x874 = -x873;
  double x875 = x2 * x489;
  double x876 = x528 - x529;
  double x877 = x874 + x875 + x876;
  double x878 = x639 - x641;
  double x879 = x877 + x878;
  double x880 = -x2 * x809;
  double x881 = -x2 * x810;
  double x882 = -x2 * x811;
  double x883 = -x2 * x813;
  double x884 = x2 * x814;
  double x885 = x2 * x816;
  double x886 = x2 * x817;
  double x887 = x2 * x818;
  double x888 = x536 + x537 + x539 - x540 - x541 - x543 + x880 + x881 + x882 + x883 +
                x884 + x885 + x886 + x887;
  double x889 = x663 + x664 + x666 - x667 - x668 - x669;
  double x890 = x393 + x879 + x888 + x889;
  double x891 = 35.0 * x427;
  double x892 = 210.0 * x445;
  double x893 = 135.0 * x624;
  double x894 = x2 * x744;
  double x895 = x2 * x95;
  double x896 = 5.0 * x895;
  double x897 = x24 * x895;
  double x898 = x2 * x92;
  double x899 = x748 * x898;
  double x900 = x422 * x45;
  double x901 = 105.0 * x628;
  double x902 = x704 + x826 - x892 - x893 - x894 - x896 + x897 + x899 + x900 + x901;
  double x903 = x233 * x708;
  double x904 = 30.0 * x416;
  double x905 = 210.0 * x687;
  double x906 = x657 * x843;
  double x907 = z * z * z * z;
  double x908 = x721 * x907;
  double x909 = x45 * x907;
  double x910 = x134 * x907;
  double x911 = x110 * x907;
  double x912 = 35.0 * x911;
  double x913 = x11 * x907;
  double x914 = -x24 * x904 + x24 * x912 + x663 + x904 - x912 - 35.0 * x913;
  double x915 = 90.0 * x416;
  double x916 = x2 * x733;
  double x917 = 105.0 * x911;
  double x918 = x907 * x95;
  double x919 = x24 * x918;
  double x920 = x907 * x92;
  double x921 = x108 * x909 + x24 * x917 - x743 * x909 - x744 * x907 + x748 * x920 +
                x749 * x913 - 210.0 * x913 - x917 - 5.0 * x918 + x919;
  double x922 = x2 * x731 + x2 * x732 - x2 * x735 - x24 * x915 - x24 * x916 +
                24.0 * x665 - 6.0 * x667 + x742 + x915 + x921;
  double x923 = x7 * x907;
  double x924 = x2 * x440;
  double x925 = x46 * x502;
  double x926 = 6.0 * x697;
  double x927 = x515 * x907;
  double x928 = x146 * x907;
  double x929 = x2 * x783 - x2 * x784 - x24 * x455 + x429 * x752 - 6.0 * x450 + x455 -
                x64 * x907 - x69 * x907 + x785 * x907 + x80 * x907 - x822 * x907 +
                x824 * x907 + x83 - 4.0 * x907 * x97 + x928;
  double x930 = 180.0 * x2;
  double x931 = x2 * x651;
  double x932 = x2 * x706;
  double x933 = x373 * x708;
  double x934 = x363 * x932;
  double x935 = 1050.0 * x292;
  double x936 = x24 * x935;
  double x937 = x24 * x478;
  double x938 = 100.0 * x24;
  double x939 = x1 * x241;
  double x940 = -x150 + x24 * x255 - x253 + x254 + 225.0 * x256 - 450.0 * x68 +
                90.0 * x79 - 60.0 * x97;
  double x941 = x130 * x710;
  double x942 = 420.0 * x138;
  double x943 = x143 * x710;
  double x944 = x24 * x941;
  double x945 = x138 * x710;
  double x946 = 420.0 * x151;
  double x947 = x24 * x727;
  double x948 = -x135 * x710 - x136 * x717 - x141 * x710 + x24 * x943 - x710 * x942 +
                x717 * x946 + x749 * x945 + x781 * x94 - x941 - 6.0 * x943 + x944 +
                945.0 * x947;
  double x949 = x186 * x526;
  double x950 = x751 * x86;
  double x951 = x1 * x545;
  double x952 = x292 * x510;
  double x953 = x342 * x717;
  double x954 = 150.0 * x270;
  double x955 = 300.0 * x68;
  double x956 = 150.0 * x256;
  double x957 = 60.0 * x79;
  double x958 = x1 * x955 - x1 * x956 - x1 * x957 + x1 * x99 + x125 - x24 * x954 +
                40.0 * x538 - 10.0 * x540 + x737 + x740 + x750 + x954;
  double x959 = x170 * x786 + x233 * x958;
  double x960 = 225.0 * x165 + 787.5 * x168 + x42 * x843;
  double x961 = 5197.5 * x459;
  double x962 = x253 * x342;
  double x963 = 2835.0 * x338;
  double x964 = -x1 * x962 + x710 * x961 + x710 * x963;
  double x965 = x187 * x526 * x719;
  double x966 = 630.0 * x292;
  double x967 = x24 * x966;
  double x968 = 630.0 * x350;
  double x969 = x402 - x733 - x819 - x820 + x821 + x823 + x825;
  double x970 = x333 * x786;
  double x971 = 2.0 * x764;
  double x972 = x604 * x611;
  double x973 = x15 * x772;
  double x974 = x393 * x973;
  double x975 = x575 + x619 + x845;
  double x976 = 630.0 * x445;
  double x977 = x130 * x849;
  double x978 = x143 * x849;
  double x979 = x24 * x849;
  double x980 = -x135 * x849 - x136 * x856 + x149 * x979 + x153 * x979 -
                75.0 * x2 * x271 + x24 * x977 + x24 * x978 - x849 * x942 + x856 * x946 +
                x867 * x94 + x969 - x977 - 6.0 * x978;
  double x981 = 112.5 * x342;
  double x982 = x1 * x671;
  double x983 = x2 * x518;
  double x984 = x292 * x649;
  double x985 = x427 * x510;
  double x986 = 7.5 * x526 * x657 + 157.5 * x985;
  double x987 = x342 * x856;
  double x988 = x2 * x282 * x839 + x477 * x890 + x518 * x869 - x560 - x709 + x849 * x961 +
                x849 * x963 + x975 + 577.5 * x987;
  double x989 = x24 * x537;
  double x990 = -x989;
  double x991 = x552 * x706;
  double x992 = x270 * x651;
  double x993 = x86 * x922;
  double x994 = 630.0 * x427;
  double x995 = x24 * x624;
  double x996 = x130 * x907;
  double x997 = x143 * x907;
  double x998 = x24 * x907;
  double x999 = -x135 * x907 - x136 * x911 - x141 * x907 + x149 * x998 + x153 * x998 +
                x24 * x996 + x24 * x997 - x907 * x942 + x911 * x946 + x928 * x94 - x996 -
                6.0 * x997;
  double x1000 = x2 * x604;
  double x1001 = x2 * x708;
  double x1002 = x342 * x911;
  double x1003 = -x2 * x962 + x907 * x961 + x907 * x963;
  double x1004 = 150.0 * x416;
  double x1005 = -x1004 * x24 + x1004 + x125 + x2 * x955 - x2 * x956 - x2 * x957 +
                 x2 * x99 + 40.0 * x665 - 10.0 * x667 + x921;
  double x1006 = 1050.0 * x427;
  double x1007 = x1001 * x186 + 262.5 * x1002 + x1003 + x1005 * x233 - x115 * x651 +
                 x157 * (-x1006 * x24 + x1006 - x241 * x872 + 2100.0 * x445 +
                         1350.0 * x624 - 1050.0 * x628 + 50.0 * x895 - x898 * x938 +
                         400.0 * x898 + x940 - 450.0 * x995 + x999) +
                 x170 * x929 - x186 * x706 - x210 * x708 - x253 * x649 +
                 75.0 * x416 * x651 + 525.0 * x427 * x649 - 4725.0 * x632 +
                 75.0 * x657 * x706 - 2100.0 * x687 + x960 + 2.5 * x993;
  double x1008 = -70.0 * x297 + 63.0 * x758 + 15.0;
  double x1009 = 45.0 * x351;
  double x1010 = pow(y, 6);
  double x1011 = -30.0 * x297 + 35.0 * x758 + 3.0;
  double x1012 = 22.5 * x331;
  double x1013 = 5.0 * x297 - 3.0;
  double x1014 = 30.0 * x1013;
  double x1015 = 3.0 * x297 - 1.0;
  double x1016 = x1010 * x130;
  double x1017 = 675.0 * x270;
  double x1018 = 1575.0 * x717;
  double x1019 = x1 * x68;
  double x1020 = x1 * x78;
  double x1021 = x1010 * x143;
  double x1022 = x23 * x725;
  double x1023 = x24 * x720;
  double x1024 = x1010 * x24;
  double x1025 = 420.0 * x152;
  double x1026 = x1 * x256;
  double x1027 = 162.0 * x1019;
  double x1028 = 72.0 * x1020;
  double x1029 = 6.0 * x538;
  double x1030 = 22.0 * x747;
  double x1031 = 105.0 * x1022;
  double x1032 = 195.0 * x720;
  double x1033 = 90.0 * x1023;
  double x1034 = x1022 * x223;
  double x1035 = 4.0 * x24 * x747;
  double x1036 = 72.0 * x1026;
  double x1037 = -x1027 - x1028 - x1029 + x1030 + x1031 + x1032 - x1033 - x1034 - x1035 +
                 x1036 + x230 + x736 - x738 + x739 - x740 + x745 + x989;
  double x1038 = x43 * x834;
  double x1039 = x1 * x182;
  double x1040 = x1038 - x1039 + x184 + x782 + x868;
  double x1041 = 12.0 * x1020;
  double x1042 = 27.0 * x1019;
  double x1043 = 12.0 * x1026;
  double x1044 = x322 * x45;
  double x1045 = x1041 + x1042 - x1043 - x1044 + x208 + x538 + x876;
  double x1046 = x179 * x558;
  double x1047 = 7.0 * x297;
  double x1048 = x1047 - 3.0;
  double x1049 = x2 * x495;
  double x1050 = 180.0 * x856;
  double x1051 = x10 * x766;
  double x1052 = x1015 * x331;
  double x1053 = x421 * x510;
  double x1054 = x1046 * x2;
  double x1055 = x436 * x502;
  double x1056 = x1040 * x187 * x455;
  double x1057 = 6.0 * x1040;
  double x1058 = x15 * x526 * x849;
  double x1059 = x1045 * x849;
  double x1060 = -x1038 + x1039 + x317 * x849 + x318 * x849 + x324 - x410 * x849 + x439 -
                 x782 + x866 - x867;
  double x1061 = x1060 * x128;
  double x1062 = 7.0 * x815;
  double x1063 = 72.0 * x478;
  double x1064 = 177.0 * x350;
  double x1065 = 72.0 * x480;
  double x1066 = -x1041 - x1042 + x1043 + x1044 + x1062 * x2 + x1063 * x2 + x1064 * x2 -
                 x1065 * x2 + x229 + x434 - x45 * x867 + x531 - x538 + x873 - x875;
  double x1067 = x2 * x478;
  double x1068 = x427 * x725;
  double x1069 = x2 * x480;
  double x1070 = x1 * x318;
  double x1071 = x1 * x317;
  double x1072 = x2 * x805;
  double x1073 = x1 * x410;
  double x1074 = 11.0 * x1067 - 16.0 * x1069 - x1070 - x1071 - x1072 + x1073 + x2 * x804 +
                 x49 + x595 + x696 + 51.0 * x871;
  double x1075 = x408 * x849;
  double x1076 = x331 * x502;
  double x1077 = x1040 * x15;
  double x1078 = x1077 * x2;
  double x1079 = x1015 * x462;
  double x1080 = x2 * x331;
  double x1081 = x1 * x907;
  double x1082 = x1060 * x2;
  double x1083 = 6.0 * x1074;
  double x1084 = 210.0 * x856;
  double x1085 = x1081 * x408;
  double x1086 = x1081 * x134;
  double x1087 = pow(z, 6);
  double x1088 = x1087 * x130;
  double x1089 = 675.0 * x416;
  double x1090 = 1575.0 * x911;
  double x1091 = x1087 * x143;
  double x1092 = x23 * x909;
  double x1093 = x1087 * x24;
  result[0] =
        -x0 * x126 * x128 -
        x157 * (x100 * x145 - 1350.0 * x100 - 15.0 * x101 - x102 * x93 - 810.0 * x103 +
                45.0 * x105 - 1575.0 * x109 + x131 * x24 - x131 + x132 * x24 - x132 -
                x133 * x24 + x133 - x135 * x20 - x137 * x20 - x138 * x139 + x139 * x152 -
                x141 * x20 + x144 * x24 - 6.0 * x144 - x145 * x91 + x147 * x20 +
                x148 * x149 + x148 * x153 + x150 * x77 + x154 + 3150.0 * x90 +
                2025.0 * x91 + 600.0 * x93 + 75.0 * x96 - 180.0 * x98) -
        x19 * x9 - x27 * (x20 * x22 + 105.0 * x5 - 315.0 * x8 - 5.0) - x28 * x30 * x41 -
        x42 * x60 * x62 - x63 * x84 * x87;  // xxxxxx
  result[1] = -x258 * y;                    // xxxxxy
  result[2] = -x258 * z;                    // xxxxxz
  result[3] = -x1 * x189 - x1 * x329 - x1 * x335 - x1 * x340 - x1 * x341 - x1 * x344 -
              x1 * x349 + x1 * x372 + x1 * x377 - x1 * x384 +
              x157 * (x1 * x237 - x1 * x239 - x1 * x265 - x1 * x266 - x1 * x267 +
                      x1 * x273 + x1 * x274 + x1 * x275 - x1 * x276 - x1 * x277 -
                      x1 * x279 + x1 * x280 + x1 * x281 + x24 * x260 - x260 + x262 +
                      x264 - x269 + x271 * x272 + x284 + x285) -
              x282 * x345 - x287 * x301 - x287 * x348 -
              x291 * (-x1 * x286 + x287 * x288 + x289) - x293 * x294 + x295 * x296 +
              x297 * x298 + x313 * x315 - x321 * x366 + x328 * x63 - x331 * x334 -
              x336 * x337 - x351 * x352 - x353 * x355 - x355 * x368 - x362 * x364 +
              x369 * x370 + x370 * x373 + x374 * x376 - x378 * x380 - x379 * x382 +
              x390;  // xxxxyy
  result[4] =
        -x414 *
        (x0 * x134 * x382 + x0 * x348 + x107 * x345 - x128 * x397 * x63 + x161 * x378 +
         x163 * x300 + 180.0 * x163 * (x158 - 2.0) + 15.0 * x166 + x169 * x28 +
         x173 * x333 + x189 - x231 * x234 - x232 - 2.0 * x32 * x404 + x329 + x335 + x340 +
         x341 + x344 + x349 + x353 * x398 + x364 * x393 + x366 * x51 + x368 * x398 -
         x369 * x399 - x372 - x373 * x399 - x376 * x66 - x377 + x384 + x406 * x84 -
         x413 * (-354.0 * x12 + x223 * x250 - x24 * x409 - 14.0 * x240 + 3.0 * x242 +
                 54.0 * x243 + 555.0 * x244 - 240.0 * x245 + 275.0 * x246 - 70.0 * x247 -
                 144.0 * x250 + 144.0 * x252 - 8.0 * x278 + x317 + x318 - x407 + x409 -
                 x410 + x411 + x412));  // xxxxyz
  result[5] =
        x157 * (x2 * x237 - x2 * x239 - x2 * x265 - x2 * x266 - x2 * x267 + x2 * x273 +
                x2 * x274 + x2 * x275 - x2 * x276 - x2 * x277 - x2 * x279 + x2 * x280 +
                x2 * x281 + x24 * x415 + x272 * x417 + x285 - x415 + x424) -
        x189 * x2 - x2 * x329 - x2 * x335 - x2 * x340 - x2 * x341 - x2 * x344 -
        x2 * x349 + x2 * x372 + x2 * x377 - x2 * x384 -
        x291 * (x288 * x426 + x289 + x425) - x294 * x428 + x296 * x429 + x298 * x430 -
        x301 * x426 + x315 * x435 - x334 * x443 - x337 * x444 - x345 * x421 -
        x348 * x426 - x352 * x446 - x353 * x447 - x364 * x453 - x366 * x436 -
        x368 * x447 + x369 * x454 + x373 * x454 + x376 * x455 - x378 * x457 -
        x382 * x456 + x390 + x442 * x63;  // xxxxzz
  result[6] =
        -x582 *
        (-x313 * x477 - x313 * x494 - x326 * x462 + x362 * x557 -
         x413 * (-x151 * x466 - 153.0 * x350 + 507.0 * x379 - 33.0 * x478 + 19.0 * x479 +
                 48.0 * x480 + 192.0 * x482 - 192.0 * x483 + x486 + x492 + x493) +
         x461 - x473 * x474 + x476 + x497 + x501 + x504 + x507 + 27.0 * x509 +
         45.0 * x511 - x513 - x517 - x520 - x521 - x523 - x525 + x527 * x61 +
         x545 * x546 + x550 + x553 + x556 + x560 * x61 + x562 + x564 + 270.0 * x566 +
         x568 + x570 + x573 + x574 + x581);  // xxxyyy
  result[7] =
        x622 *
        (-x1 * x614 +
         x157 * (-x151 * x308 - 459.0 * x350 + 1521.0 * x379 - 99.0 * x478 + 57.0 * x479 +
                 144.0 * x480 + 576.0 * x482 - 576.0 * x483 - x583 + x584 + x585 - x586 +
                 x588) +
         x233 * x313 + x235 * x313 + x295 * x404 + x297 * x601 - x362 * x609 - x461 -
         x497 - x501 - x504 - x507 - 9.0 * x509 - 15.0 * x511 + x513 + x517 + x520 +
         x521 + x523 + x525 - x528 * x615 - x550 - x553 - x556 - x562 - x564 -
         300.0 * x566 - x568 - x570 - x573 - x574 + x597 + x599 + x600 + x602 * x603 +
         x602 * x604 - x605 * x606 - x608 - x612 * x63 - x617 * x63 + x621);  // xxxyyz
  result[8] =
        x582 * (x157 * (-x24 * x626 - x432 * x45 - 459.0 * x445 + 1521.0 * x456 + x588 -
                        x623 - 99.0 * x624 + x625 + x626 + 57.0 * x627 + 144.0 * x628 +
                        576.0 * x629 - 576.0 * x630) -
                x177 * x606 * x659 - x2 * x614 + x233 * x435 + x235 * x435 + x404 * x429 +
                x430 * x601 - x453 * x609 + x603 * x654 + x604 * x654 - x615 * x639 +
                x621 - x63 * x672 - x63 * x686 - x631 - x634 + x644 + x645 + x646 - x648 -
                15.0 * x650 - 9.0 * x652 - x653 - x655 - x656 - x658 - x660 - x662 -
                x673 - x675 + x676 + x678 + x679 + x681 + x682 + x683 - x684 -
                300.0 * x685 - x689 - x690 - x693);  // xxxyzz
  result[9]  = -x622 * (-x413 * (-x151 * x694 + 507.0 * x456 + x587 + 19.0 * x627 +
                                192.0 * x629 - 192.0 * x630 + x700 + x704 + x705) -
                       x435 * x477 - x435 * x494 - x440 * x462 + x453 * x557 -
                       x474 * x697 + x546 * x708 + x581 + x61 * x707 + x61 * x709 + x631 +
                       x634 + x648 + 45.0 * x650 + 27.0 * x652 + x653 + x655 + x656 +
                       x658 + x660 + x673 + x675 - x676 - x678 - x679 - x681 - x682 -
                       x683 + x684 + 270.0 * x685 + x689 + x690 + x693 + x698);  // xxxzzz
  result[10] = a * x205 * x761 - x1 * x214 * x508 + x1 * x326 * x769 + x175 * x282 -
               x282 * x796 - x287 * x373 * x545 + x287 * x714 - x314 * x751 -
               x327 * x712 - x328 - x374 * x375 + x374 * x505 + 990.0 * x380 +
               x383 * x795 + x389 * x510 +
               x413 * (-x1 * x407 + x1 * x411 - x1 * x723 - x1 * x724 + x1 * x730 +
                       x24 * x718 - x24 * x719 - x24 * x722 + x397 + x536 - x718 + x719 -
                       35.0 * x720 + x722 + x725 * x726 + x727 * x728 - x727 * x729) +
               x46 * x72 * x761 + x508 * x67 - x511 * x795 - x516 * x753 + x583 * x713 +
               x604 * x768 - x63 * x787 - x63 * x799 - x710 * x711 - x710 * x716 -
               x710 * x757 - x710 * x779 - x710 * x789 - x710 * x791 - x710 * x792 -
               x717 * x776 + x725 * x774 - x725 * x777 - x727 * x794 - 6.0 * x755 +
               x758 * x760 + x762 * x763 + x764 * x765 + x765 * x767 + x770 * x771 +
               x770 * x772 + 45.0 * x775 - x797 * x798 + 90.0 * x800 + x801;  // xxyyyy
  result[11] = -x414 * (x1 * x789 + x1 * x838 + x282 * x840 + x292 * x776 + x314 * x827 -
                        x326 * x807 - x336 * x836 -
                        x413 * (-35.0 * x350 + 443.0 * x379 + 93.0 * x482 - 128.0 * x483 +
                                x486 - x804 + x805 + x806) -
                        x473 * x477 - x473 * x494 - x475 * x602 + x476 - x494 * x596 +
                        54.0 * x509 + 90.0 * x511 - x516 * x771 - x516 * x830 -
                        x516 * x831 - x518 * x602 + x526 * x832 + x527 * x63 +
                        x545 * x833 + x560 * x63 + 510.0 * x566 - x597 - x600 + x608 +
                        x611 * x835 + x63 * x828 + x63 * x841 + x63 * x844 + x756 * x834 +
                        x803 + x808 - x829 + x837 + x842 + x846);  // xxyyyz
  result[12] =
        x1 * x333 * x440 + x1 * x39 * x500 + x1 * x645 - 21.0 * x1 * x652 +
        x155 * x336 * x440 +
        x157 * (1329.0 * x2 * x379 + x2 * x585 - x2 * x586 - x308 + x309 -
                384.0 * x379 * x872 + x397 - x431 + x432 + 279.0 * x481 * x849 + x542 -
                x590 - x592 + x594 + x635 - x636 - x637 + x638 - 105.0 * x871 + x879) -
        15.0 * x18 + 6.0 * x192 * x861 + x2 * x326 * x333 - x2 * x342 * x63 * x719 -
        21.0 * x2 * x509 + x2 * x599 - x2 * x803 - x2 * x808 + x2 * x829 - x2 * x837 -
        x2 * x842 - x233 * x326 - x233 * x440 - x235 * x326 - x235 * x440 + x25 * x484 +
        x25 * x699 + x270 * x854 - 7.5 * x270 * x858 - x287 * x363 * x671 - x291 +
        x293 * x848 + 2.0 * x295 * x502 * x654 - x30 * x547 - x314 * x890 - x321 * x375 +
        x321 * x505 - x321 * x63 * x651 + x326 * x864 + x333 * x768 + x333 * x857 -
        x363 * x426 * x611 + x363 * x861 - x371 * x762 - x371 * x859 - x375 * x436 +
        165.0 * x380 + x383 * x851 + x383 * x852 - x385 * x63 - x39 * x406 - x39 * x598 +
        x39 * x674 - 7.5 * x416 * x796 + x416 * x854 + x428 * x848 + x436 * x505 -
        x436 * x508 * x63 + 165.0 * x457 + x508 * x52 + x510 * x853 - x511 * x852 +
        x52 * x651 - 540.0 * x565 * x849 + x596 * x863 + x596 * x865 -
        60.0 * x613 * x849 + x643 * x764 + x643 * x767 + x649 * x853 - x650 * x851 -
        x755 - x757 * x849 + x771 * x862 + 7.5 * x775 - x776 * x856 - x789 * x849 -
        x794 * x850 + 15.0 * x800 + x830 * x862 + x831 * x862 - x838 * x849 +
        15.0 * x847 + 7.5 * x855 - x860 - x869 * x870;  // xxyyzz
  result[13] =
        -x414 *
        (x2 * x757 + x2 * x789 + x2 * x838 + x314 * x902 - x375 * x904 -
         x413 * (x24 * x891 - 35.0 * x445 + 443.0 * x456 + 93.0 * x629 - 128.0 * x630 +
                 x700 + x806 - x891) +
         x416 * x778 + x421 * x840 + x426 * x790 + x426 * x802 + x427 * x776 -
         x440 * x807 - x444 * x836 - x475 * x654 - x477 * x697 - x494 * x643 -
         x494 * x697 - x518 * x654 + x63 * x633 + x63 * x707 + x63 * x709 + x63 * x903 +
         x63 * x905 + x63 * x906 - x644 - x646 + 90.0 * x650 + 54.0 * x652 + x662 +
         x671 * x835 - x680 * x771 - x680 * x830 - x680 * x831 + 510.0 * x685 + x698 +
         x706 * x832 + x708 * x833 + x846);  // xxyzzz
  result[14] =
        x175 * x421 - x214 * x931 - x314 * x922 - x375 * x455 + x383 * x930 +
        x389 * x649 +
        x413 * (-x2 * x407 + x2 * x411 - x2 * x723 - x2 * x724 + x2 * x730 - x24 * x908 +
                x397 + x726 * x909 + x728 * x910 - x729 * x910 + x908 + x914) -
        x421 * x858 + x426 * x714 - x426 * x933 - x441 * x712 - x442 + x455 * x505 +
        990.0 * x457 + x604 * x857 + x604 * x924 + x623 * x713 - x63 * x934 -
        x650 * x930 + x651 * x67 - x680 * x753 - x711 * x907 - x716 * x907 - x757 * x907 +
        x760 * x923 + x763 * x859 + x769 * x924 + x771 * x927 + x772 * x927 +
        x774 * x909 - x776 * x911 - x777 * x909 - x779 * x907 - x789 * x907 -
        x791 * x907 - x792 * x907 - x794 * x910 - x798 * x932 + x801 + 90.0 * x847 +
        45.0 * x855 - 6.0 * x860 + x863 * x926 + x865 * x926 - x870 * x929 +
        x924 * x925;  // xxzzzz
  result[15] =
        -x582 * (-x115 * x508 +
                 x157 * (-x24 * x939 + 2100.0 * x350 + 1350.0 * x478 - 1050.0 * x480 +
                         50.0 * x812 - x815 * x938 + 400.0 * x815 + x935 - x936 -
                         450.0 * x937 + x940 + x948) +
                 x186 * x951 - x210 * x545 - x253 * x510 + 75.0 * x270 * x508 -
                 4725.0 * x495 + 75.0 * x526 * x552 - 2100.0 * x571 - x949 + 2.5 * x950 +
                 525.0 * x952 + 262.5 * x953 + x959 + x960 + x964);  // xyyyyy
  result[16] =
        -x622 *
        (x1 * x972 - x107 * x510 - x128 * x611 +
         x157 * (-x104 * x815 - x24 * x968 + 1260.0 * x350 + 810.0 * x478 + 30.0 * x812 -
                 6.0 * x814 + 240.0 * x815 - 270.0 * x937 + x948 + x966 - x967 + x969) +
         x373 * x951 - x374 * x839 - 2835.0 * x495 - x508 * x66 + x508 * x738 +
         x510 * x583 - 1260.0 * x571 + x710 * x974 + x827 * x971 + x950 + 472.5 * x953 +
         x964 + x965 + x970 + x975);  // xyyyyz
  result[17] =
        -x582 *
        (x157 * (-x223 * x898 - x24 * x893 - x437 * x94 + x488 - x489 + x623 +
                 405.0 * x624 - x625 - 315.0 * x628 + x809 + x810 + x811 + x813 - x814 -
                 x816 - x817 - x818 + 15.0 * x895 + 120.0 * x898 + x976 + x980) +
         x2 * x972 +
         x233 * (x2 * x819 + x2 * x820 - x2 * x821 - x24 * x664 + x423 - 3.0 * x667 +
                 x742 + x877 + x888 + x916) -
         x416 * x981 - x453 * x604 - x496 + x508 * x691 + x518 * x982 - x527 +
         x545 * x983 - x607 - 1417.5 * x632 + x651 * x719 - x672 - x686 - 630.0 * x687 +
         x827 * x863 - x828 - x841 - x844 + 52.5 * x984 + x986 + x988);  // xyyyzz
  result[18] = -x622 * (x1 * x518 * x708 +
                        x157 * (-x223 * x815 - x24 * x810 - x322 * x94 + 405.0 * x478 -
                                315.0 * x480 + x583 - x584 + x701 - x702 + 15.0 * x812 +
                                120.0 * x815 + x892 + x893 + x894 + x896 - x897 - x899 -
                                x900 - x901 + x968 + x980) +
                        x233 * (x1 * x819 + x1 * x820 - x1 * x821 + x284 - 3.0 * x540 +
                                x734 + x742 + x874 + x875 + x878 + x880 + x881 + x882 +
                                x883 + x884 + x885 + x886 + x887 + x889 + x990) -
                        x270 * x981 - x362 * x604 - 1417.5 * x495 + x508 * x904 -
                        630.0 * x571 + x604 * x982 + x611 * x983 - x612 - x617 - x633 -
                        x661 - x707 + x764 * x902 - x903 - x905 - x906 + 157.5 * x984 +
                        52.5 * x985 + x988 + 7.5 * x991 + 67.5 * x992);  // xyyzzz
  result[19] = -x582 * (x1000 * x671 + x1001 * x373 + 472.5 * x1002 + x1003 -
                        x107 * x649 - x128 * x671 +
                        x157 * (-x104 * x898 - x24 * x976 - x24 * x994 + 1260.0 * x445 +
                                810.0 * x624 + 30.0 * x895 - 6.0 * x897 + 240.0 * x898 +
                                x969 + x994 - 270.0 * x995 + x999) +
                        x187 * x706 * x904 + x333 * x929 - x455 * x839 + x623 * x649 -
                        2835.0 * x632 - x651 * x66 + x651 * x915 - 1260.0 * x687 +
                        2.0 * x863 * x902 + x907 * x974 + x975 + x993);  // xyzzzz
  result[20] = -x1007 * x622;                                            // xzzzzz
  result[21] = -x1 * x128 * x958 - x1008 * x1009 - x1011 * x1012 * x42 -
               x1014 * x42 * x797 - x1015 * x786 * x87 -
               x157 * (x1 * x150 * x24 + x1010 * x1025 - x1010 * x135 - x1010 * x137 -
                       x1010 * x141 + x1010 * x147 - x1010 * x942 + x1016 * x24 - x1016 +
                       x1017 * x24 - x1017 - x1018 * x24 + x1018 - 1350.0 * x1019 -
                       x102 * x747 - 810.0 * x1020 + x1021 * x24 - 6.0 * x1021 -
                       x1022 * x145 + 2025.0 * x1022 - 1575.0 * x1023 + x1024 * x149 +
                       x1024 * x153 + 675.0 * x1026 + x154 - 180.0 * x538 + 45.0 * x540 +
                       3150.0 * x720 + 75.0 * x745 - 15.0 * x746 + 600.0 * x747) -
               x27 * (x1010 * x22 + 105.0 * x297 - 315.0 * x758 - 5.0);  // yyyyyy
  result[22] =
        -x414 *
        (x1008 * x167 + x1008 * x169 + 37.5 * x1011 * x165 + x1011 * x170 * x172 -
         x1013 * x1040 * x186 + 15.0 * x1013 * x1046 + 75.0 * x1013 * x605 -
         x1015 * x1045 * x210 + x1015 * x193 * x526 + x1015 * x949 - 2.5 * x1037 * x86 +
         150.0 * x1048 * x571 -
         x157 * (-x104 * x945 + 35.0 * x140 * x710 - 315.0 * x151 * x717 - x238 * x270 +
                 x238 * x717 - x24 * x249 * x710 + x24 * x811 + x251 * x478 + x257 -
                 1950.0 * x350 + 900.0 * x480 + 1785.0 * x727 - 220.0 * x815 - x935 +
                 x936 - x939 + x941 + x943 - x944 + 285.0 * x945 - 840.0 * x947) +
         2.5 * x192 * x786 + x235 * x958 + 210.0 * x495 * (9.0 * x297 - 5.0) +
         150.0 * x952 + x959 + x965);  // yyyyyz
  result[23] =
        -x1 * x1013 * x363 * x453 - x1011 * x333 * x443 + x1011 * x385 -
        52.5 * x1011 * x428 - 15.0 * x1011 * x446 - 60.0 * x1013 * x2 * x569 -
        x1013 * x793 * x850 + x1014 * x351 - x1015 * x1053 - 21.0 * x1015 * x1054 +
        x1015 * x1056 + x1015 * x1057 * x192 * x2 + x1015 * x1061 + x1037 * x863 +
        x1037 * x865 + 36.0 * x1040 * x2 * x552 - 300.0 * x1048 * x1049 -
        x1048 * x381 * x850 - x1050 * x510 - 36.0 * x1051 * x331 * x849 - x1052 * x1055 +
        x1052 * x603 + x1052 * x604 - x1058 * x353 - x1058 * x368 + x1059 * x369 +
        x1059 * x373 + x1066 * x971 +
        x157 * (x1027 + x1028 + x1029 - x1030 - x1031 - x1032 + x1033 + x1034 + x1035 -
                x1036 - 432.0 * x1067 - 210.0 * x1068 * x24 + 825.0 * x1068 +
                432.0 * x1069 + x2 * x45 * x739 + 1665.0 * x2 * x727 - 42.0 * x2 * x815 +
                x2 * x941 - x2 * x944 + 162.0 * x2 * x945 - 720.0 * x2 * x947 -
                x2 * x966 + x2 * x967 + x208 + 9.0 * x417 * x710 + x424 + x741 - x745 -
                1062.0 * x871 - 24.0 * x872 * x945 + x990) -
        x2 * x339 * x710 - x2 * x598 * x786 - x2 * x965 - x2 * x970 + x235 * x786 -
        30.0 * x293 * (-x1047 + x288 * x849 + x425 + 3.0) - x343 * x427 * x710 +
        x381 * x727 + x387 * x797 + x506 * x510 - x786 * x864 + x787 + x799;  // yyyyzz
  result[24] = -x414 *
               (-x1013 * x575 - x1013 * x580 + x1013 * x633 + x1013 * x688 +
                x1013 * x707 + x1013 * x709 - 67.5 * x1015 * x165 + x1015 * x477 * x708 +
                x1015 * x688 + x1015 * x692 + x1040 * x462 + x1040 * x475 -
                67.5 * x1040 * x657 - x1045 * x2 * x475 - x1045 * x332 * x444 +
                x1045 * x477 + x1045 * x494 - x1045 * x983 + 85.5 * x1051 * x1080 +
                x1055 * x526 - x1060 * x462 - x1060 * x475 - x1066 * x477 - x1066 * x494 -
                x1074 * x474 + x1076 * x647 - x1078 * x514 - x1078 * x522 - x1079 * x172 +
                x1079 * x453 + x1080 * x332 * x548 + x2 * x526 * x555 -
                x413 * (-x1062 - x1063 - x1064 + x1065 - x1072 * x45 - x1075 * x24 +
                        x1075 + 19.0 * x138 * x849 - 192.0 * x24 * x850 + x45 * x529 +
                        192.0 * x45 * x856 + x492 + x703 + x705 + 507.0 * x850) +
                x460 * x849 - x475 * x526 - 630.0 * x495 - x514 * x558 - x518 * x526 -
                x522 * x558 + x526 * x674 - x526 * x807 + x567 * x849 - 27.0 * x569 -
                315.0 * x571 - 67.5 * x605 + 45.0 * x984 + x986 + 270.0 * x987 +
                27.0 * x992);  // yyyzzz
  result[25] =
        -52.5 * x10 * x558 * x909 + x1000 * x1060 - x1009 - x1012 * x23 * x502 * x907 -
        x1015 * x233 * x929 + x1015 * x25 * x623 - x1015 * x27 - x1015 * x421 * x649 +
        x1015 * x443 * x604 + 90.0 * x1015 * x446 - x1015 * x715 * x907 -
        x1015 * x793 * x910 - x1015 * x934 + x1040 * x907 * x973 + 1890.0 * x1049 -
        x1050 * x649 + x1053 + 45.0 * x1054 - x1056 - x1057 * x429 * x502 -
        x1060 * x15 * x712 - x1061 + x1076 * x455 + x1077 * x759 * x923 +
        x1077 * x771 * x907 + x1077 * x773 * x909 - x1078 * x753 - x1081 * x790 -
        x1081 * x802 + x1082 * x769 + x1082 * x925 + x1083 * x863 + x1083 * x865 +
        990.0 * x17 * x850 - 12.0 * x2 * x991 - 420.0 * x292 * x342 * x907 - 90.0 * x293 -
        x331 * x475 - x331 * x518 - x331 * x788 * x907 - x331 * x807 + x374 * x651 +
        x413 * (-66.0 * x1067 + 96.0 * x1069 + x1070 + x1071 - x1073 + x1084 * x24 -
                x1084 - x1085 * x24 + x1085 - 128.0 * x1086 * x24 + 443.0 * x1086 +
                93.0 * x292 * x909 + x396 - 306.0 * x871 + x876 + x914) +
        x499 * x558 * x859 + x506 * x649 - 52.5 * x510 * x911 - 5.0 * x548 * x558 * x907 +
        x569 * x930 - x738 * x931 - x764 * x922 - x849 * x933;  // yyzzzz
  result[26] = -x1007 * x414;                                   // yzzzzz
  result[27] =
        -x1005 * x128 * x2 -
        x157 * (-x102 * x920 + x1025 * x1087 - x1087 * x135 - x1087 * x137 -
                x1087 * x141 + x1087 * x147 - x1087 * x942 + x1088 * x24 - x1088 +
                x1089 * x24 - x1089 - x1090 * x24 + x1090 + x1091 * x24 - 6.0 * x1091 -
                x1092 * x145 + 2025.0 * x1092 + x1093 * x149 + x1093 * x153 +
                x150 * x872 + x154 + 675.0 * x2 * x256 - 1350.0 * x2 * x68 -
                810.0 * x2 * x78 - 1575.0 * x24 * x913 - 180.0 * x665 + 45.0 * x667 +
                3150.0 * x913 + 75.0 * x918 - 15.0 * x919 + 600.0 * x920) -
        x27 * (x1087 * x22 + 105.0 * x430 - 315.0 * x923 - 5.0) -
        22.5 * x42 * x443 * (-30.0 * x430 + 35.0 * x923 + 3.0) -
        x42 * x706 * x852 * (5.0 * x430 - 3.0) -
        45.0 * x446 * (-70.0 * x430 + 63.0 * x923 + 15.0) -
        x87 * x929 * (3.0 * x430 - 1.0);  // zzzzzz
  return result;
}
}  // namespace tensors
}  // namespace libcppe

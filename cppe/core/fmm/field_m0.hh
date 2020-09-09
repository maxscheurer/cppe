#pragma once
void field_m0_P2M_1(double x, double y, double z, double q, double* M);
void field_m0_M2M_1(double x, double y, double z, double* M, double* Ms);
void field_m0_M2L_1(double x, double y, double z, double* M, double* L);
void field_m0_L2L_1(double x, double y, double z, double* L, double* Ls);
void field_m0_L2P_1(double x, double y, double z, double* L, double* F);
void field_m0_M2P_1(double x, double y, double z, double* M, double* F);
template <int m_order, int osize>
void P2P(double x, double y, double z, double* S, double* F);
void field_m0_P2M_2(double x, double y, double z, double q, double* M);
void field_m0_M2M_2(double x, double y, double z, double* M, double* Ms);
void field_m0_M2L_2(double x, double y, double z, double* M, double* L);
void field_m0_L2L_2(double x, double y, double z, double* L, double* Ls);
void field_m0_L2P_2(double x, double y, double z, double* L, double* F);
void field_m0_M2P_2(double x, double y, double z, double* M, double* F);
void field_m0_P2M_3(double x, double y, double z, double q, double* M);
void field_m0_M2M_3(double x, double y, double z, double* M, double* Ms);
void field_m0_M2L_3(double x, double y, double z, double* M, double* L);
void field_m0_L2L_3(double x, double y, double z, double* L, double* Ls);
void field_m0_L2P_3(double x, double y, double z, double* L, double* F);
void field_m0_M2P_3(double x, double y, double z, double* M, double* F);
void field_m0_P2M_4(double x, double y, double z, double q, double* M);
void field_m0_M2M_4(double x, double y, double z, double* M, double* Ms);
void field_m0_M2L_4(double x, double y, double z, double* M, double* L);
void field_m0_L2L_4(double x, double y, double z, double* L, double* Ls);
void field_m0_L2P_4(double x, double y, double z, double* L, double* F);
void field_m0_M2P_4(double x, double y, double z, double* M, double* F);
void field_m0_P2M_5(double x, double y, double z, double q, double* M);
void field_m0_M2M_5(double x, double y, double z, double* M, double* Ms);
void field_m0_M2L_5(double x, double y, double z, double* M, double* L);
void field_m0_L2L_5(double x, double y, double z, double* L, double* Ls);
void field_m0_L2P_5(double x, double y, double z, double* L, double* F);
void field_m0_M2P_5(double x, double y, double z, double* M, double* F);
void field_m0_P2M_6(double x, double y, double z, double q, double* M);
void field_m0_M2M_6(double x, double y, double z, double* M, double* Ms);
void field_m0_M2L_6(double x, double y, double z, double* M, double* L);
void field_m0_L2L_6(double x, double y, double z, double* L, double* Ls);
void field_m0_L2P_6(double x, double y, double z, double* L, double* F);
void field_m0_M2P_6(double x, double y, double z, double* M, double* F);
void field_m0_P2M_7(double x, double y, double z, double q, double* M);
void field_m0_M2M_7(double x, double y, double z, double* M, double* Ms);
void field_m0_M2L_7(double x, double y, double z, double* M, double* L);
void field_m0_L2L_7(double x, double y, double z, double* L, double* Ls);
void field_m0_L2P_7(double x, double y, double z, double* L, double* F);
void field_m0_M2P_7(double x, double y, double z, double* M, double* F);
template <int m_order, int osize>
void P2M(double x, double y, double z, double q, double* M, int order);
template <int m_order, int osize>
void M2M(double x, double y, double z, double* M, double* Ms, int order);
template <int m_order, int osize>
void M2L(double x, double y, double z, double* M, double* L, int order);
template <int m_order, int osize>
void L2L(double x, double y, double z, double* L, double* Ls, int order);
template <int m_order, int osize>
void L2P(double x, double y, double z, double* L, double* F, int order);
template <int m_order, int osize>
void M2P(double x, double y, double z, double* M, double* F, int order);

#ifndef INCLUDED_kuiper_H  //include-guard
#define INCLUDED_kuiper_H  //include-guard


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/stat.h>

#define DIRECTORY ../data/Interaction_20Myr_1Hill_group/N100_1/  //ファイル保存用のディレクトリ
#define _STR(str) #str
#define STR(str) _STR(str)

#define Np 4  //惑星の数
#define Nt 100  //テスト粒子の数

#define RAND_SEED 1

#define TAU_MOVE (2.0*M_PI*5.0E5)  //0.5Myr  惑星移動の時間
#define TAU_MASS (2.0*M_PI*5.0E6)  //5Myr  質量を増加させる時間

#define T_MAX (2.0*M_PI*2.0E7)  //20Myr  全計算時間

#define TIME_INTERVAL 2.0*M_PI*1.0E3,2.0*M_PI*1.0E4,2.0*M_PI*1.0E5,2.0*M_PI*1.0E6,2.0*M_PI*1.0E7,T_MAX  //t_ene配列の中身
#define TIME_INTERVAL_MAX 6  //t_ene配列の要素数

#define G 1
#define M_0 1
#define EPSILON 0.0
#define ETA 0.05
#define ITE_MAX 3

#define INNER_AXIS 28.0  //テスト粒子の長軸半径　下限
#define OUTER_AXIS 52.0  //テスト粒子の長軸半径　上限
#define UPPER_ECC 0.05  //テスト粒子の離心率　上限
#define UPPER_ICC 0.05  //テスト粒子の軌道傾斜角　上限

/*惑星の番号*/
#define JUP 1
#define SAT 2
#define URA 3
#define NEP 4

#define R_MIN 1  //0以外のときr_min計算
#define INTERACTION 1  //0以外のとき惑星同士の相互作用を入れる
#define SWAP 1  //0以外のとき1R_H以内に入った粒子を一番後ろの粒子と交換。そしてNをN-1としloopの回数を減らす。R_MINとセット。

/*0以外のときファイル作成*/
#define ENERGY_FILE 0
#define ORBITALELEMENTS_FILE 1



struct orbital_elements{
  char name[30];
  double m;
  double e;
  double axis;
  double u;
  double I;
  double OMEGA;
  double omega;
  double M;
  double R_H;
  double axis_0;
  double axis_f;
  double m_f;
  double delta_axis;
  int judge;
  int orinum;
};

double MeanLongitude_to_EccentricAnomaly(int i,struct orbital_elements ele[]);

double AngularMomentum(int i,int k,struct orbital_elements ele[],double x_0[][4],double v_0[][4],double abs_L_0,int N);

double InnerProduct(int i,int k,double x_0[][4],double v_0[][4],double r_dot_v[]);

double RadiusFromCenter(int i,int k,double x_0[][4],double r_0[]);

double SquareOfVelocity(int i,int k,double v_0[][4],double v2_0[]);

double Timestep_i_0(int i,int k,double a_0[][4],double adot[][4],double abs_a[],double abs_adot[],double dt_[]);

double Timestep_i_sys(int i_sys,int k,double a[][4],double adot[][4],double adot2_dt2[][4],double adot3_dt3[][4],double abs_a[],double abs_adot[],double abs_adot2[],double abs_adot3[],double dt_[]);

double SquareOfRaletiveDistance(int i,int j,double x_0[][4]);

double SquareOfRaletiveVelocity(int i,int j,double v_0[][4]);

double RaletiveInnerProduct(int i,int j,double x_0[][4],double v_0[][4]);

double Acceleration_ij(int i,int j,int k,struct orbital_elements ele[],double x_0[][4],double abs_r2[]);

double dAcceleration_ij(int i,int j,int k,struct orbital_elements ele[],double x_0[][4],double v_0[][4],double r_dot_v_ij[],double abs_r2[]);

double External_Acceleration(int i,int k,double x_0[][4],
double r_0[]);

double External_dAcceleration(int i,int k,double x_0[][4],double v_0[][4],double r_0[],double r_dot_v[]);

double All_Acceleration(int i,int j,int k,struct orbital_elements ele[],double x_0[][4],double r_0[],double abs_r2[],double a_0[][4]);

double All_dAcceleration(int i,int j,int k,struct orbital_elements ele[],double x_0[][4],double v_0[][4],double r_dot_v[],double r_dot_v_ij[],double r_0[],double abs_r2[],double adot_0[][4]);

double RotatingCoordinate_X(struct orbital_elements ele[],double x_0[][4],double t_sys);

double RotatingCoordinate_Y(struct orbital_elements ele[],double x_0[][4],double t_sys);

int Calculate_OrbitalElements(int i,int k,double x_c[][4],double v_c[][4],struct orbital_elements ele[],double P[][4],double Q[][4],double r_c[],double v2_c[],double r_dot_v[],double step,double t_sys,int N);

double Calculate_P(int i,int k,struct orbital_elements ele[]);

double Calculate_Q(int i,int k,struct orbital_elements ele[]);

double Calculate_Energy(int i,int j,struct orbital_elements ele[],double x_c[][4],double v_c[][4],double v2_c[],double r_c[],double abs_r[],double abs_r2[],double abs_v[],double abs_v2[],double r_dot_v_ij[],double E[],double E_tot,int N);

void InitialCondition(int i,int k,double P[][4],double Q[][4],double x_0[][4],double v_0[][4],double r_0[],struct orbital_elements ele[]);

void Predictor(int i,int k,double x_0[][4],double v_0[][4],double a_0[][4],double adot_0[][4],double x_p[][4],double v_p[][4],double r_p[],double v2_p[],double Dt[]);

void Corrector_sys(int i_sys,int j,int k,struct orbital_elements ele[],double x_p[][4],double v_p[][4],double r_p[],double v2_p[],double x_c[][4],double v_c[][4],double r_c[],double v2_c[],double a_0[][4],double adot_0[][4],double a[][4],double adot[][4],double adot2_dt2[][4],double adot3_dt3[][4],double abs_r[],double abs_r2[],double abs_v[],double abs_v2[],double r_dot_v_ij[],double r_dot_v[],double t_sys,double dt_[]);

void Iteration_sys(int i_sys,int j,int k,struct orbital_elements ele[],double x_p[][4],double v_p[][4],double x_c[][4],double v_c[][4],double r_c[],double v2_c[],double a_0[][4],double adot_0[][4],double a[][4],double adot[][4],double adot2_dt2[][4],double adot3_dt3[][4],double abs_r[],double abs_r2[],double abs_v[],double abs_v2[],double r_dot_v_ij[],double r_dot_v[],double t_sys,double dt_[]);



#endif //include-guard

#include "kuiper.h"

double MeanLongitude_to_EccentricAnomaly(int i,struct orbital_elements ele[]){
  ele[i].u = ele[i].M;
  ele[i].u += (ele[i].e - ele[i].e*ele[i].e*ele[i].e/8.0 + ele[i].e*ele[i].e*ele[i].e*ele[i].e*ele[i].e/192.0)*sin(ele[i].M);
  ele[i].u += (ele[i].e*ele[i].e/2.0 - ele[i].e*ele[i].e*ele[i].e*ele[i].e/6.0 + ele[i].e*ele[i].e*ele[i].e*ele[i].e*ele[i].e*ele[i].e/48.0)*sin(2.0*ele[i].M);
  ele[i].u += (ele[i].e*ele[i].e*ele[i].e*3.0/8.0 - ele[i].e*ele[i].e*ele[i].e*ele[i].e*ele[i].e*27.0/128.0)*sin(3.0*ele[i].M);
  ele[i].u += (ele[i].e*ele[i].e*ele[i].e*ele[i].e/3.0 - ele[i].e*ele[i].e*ele[i].e*ele[i].e*ele[i].e*ele[i].e*4.0/15.0)*sin(4.0*ele[i].M);
  ele[i].u += ele[i].e*ele[i].e*ele[i].e*ele[i].e*ele[i].e*125.0/384.0*sin(5.0*ele[i].M);
  ele[i].u += ele[i].e*ele[i].e*ele[i].e*ele[i].e*ele[i].e*ele[i].e*27.0/80.0*sin(6.0*ele[i].M);
  return ele[i].u;
}



/*r_i,v_iの内積*/
double InnerProduct(int i,int k,double x_0[][4],double v_0[][4],double r_dot_v[]){ 
  r_dot_v[i] = 0.0;
  for(k=1;k<=3;k++){
    r_dot_v[i] += x_0[i][k]*v_0[i][k]; 
  }
  return r_dot_v[i];
}

/*中心星からの距離の2乗*/
double RadiusFromCenter(int i,int k,double x_0[][4],double r_0[]){
  r_0[i] = 0.0;
  for(k=1;k<=3;k++){
    r_0[i] += x_0[i][k]*x_0[i][k];  
  }
  r_0[i] = sqrt(r_0[i]);
  //printf("r_0[%d]=%f\n",i,r_0[i]);
  return r_0[i];
}

/*速度の2乗*/
double SquareOfVelocity(int i,int k,double v_0[][4],double v2_0[]){
  v2_0[i] = 0.0;
  for(k=1;k<=3;k++){
    v2_0[i] += v_0[i][k]*v_0[i][k]; 
  }
  return v2_0[i];
}

/*相対距離の2乗*/
double SquareOfRaletiveDistance(int i,int j,double x_0[][4]){
  return (x_0[j][1] - x_0[i][1])*(x_0[j][1] - x_0[i][1]) + (x_0[j][2] - x_0[i][2])*(x_0[j][2] - x_0[i][2]) + (x_0[j][3] - x_0[i][3])*(x_0[j][3] - x_0[i][3]);
}

/*相対速度の2乗*/
double SquareOfRaletiveVelocity(int i,int j,double v_0[][4]){
  return (v_0[j][1] - v_0[i][1])*(v_0[j][1] - v_0[i][1]) + (v_0[j][2] - v_0[i][2])*(v_0[j][2] - v_0[i][2]) + (v_0[j][3] - v_0[i][3])*(v_0[j][3] - v_0[i][3]);
}

/*r_ij, v_ijの内積*/
double RaletiveInnerProduct(int i,int j,double x_0[][4],double v_0[][4]){
  return (x_0[j][1] - x_0[i][1])*(v_0[j][1] - v_0[i][1]) + (x_0[j][2] - x_0[i][2])*(v_0[j][2] - v_0[i][2]) + (x_0[j][3] - x_0[i][3])*(v_0[j][3] - v_0[i][3]);
}



/*回転座標系への変換 x*/
double RotatingCoordinate_X(struct orbital_elements ele[],double x_0[][4],double t_sys){
  return x_0[2][1]*cos(sqrt(G*M_0/ele[1].axis/ele[1].axis/ele[1].axis)*t_sys) + x_0[2][2]*sin(sqrt(G*M_0/ele[1].axis/ele[1].axis/ele[1].axis)*t_sys);
}

/*回転座標系への変換 y*/
double RotatingCoordinate_Y(struct orbital_elements ele[],double x_0[][4],double t_sys){
  return -x_0[2][1]*sin(sqrt(G*M_0/ele[1].axis/ele[1].axis/ele[1].axis)*t_sys) + x_0[2][2]*cos(sqrt(G*M_0/ele[1].axis/ele[1].axis/ele[1].axis)*t_sys);
}








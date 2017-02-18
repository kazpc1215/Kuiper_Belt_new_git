#include "kuiper.h"


/*エネルギー計算*/
double Calculate_Energy(int i,int j,struct orbital_elements ele[],double x_c[][4],double v_c[][4],double v2_c[],double r_c[],double abs_r[],double abs_r2[],double abs_v[],double abs_v2[],double r_dot_v_ij[],double E[],double E_tot,int N){
  double rij1;
  double r1;
  for(i=1;i<=N;++i){
    E[i] = 0.5*ele[i].m*v2_c[i];
    for(j=1;j<=N;++j){
      if(i!=j){
	abs_r2[j] = SquareOfRaletiveDistance(i,j,x_c); //絶対値2乗
	abs_v2[j] = SquareOfRaletiveVelocity(i,j,v_c);
	r_dot_v_ij[j] = RaletiveInnerProduct(i,j,x_c,v_c);  //r_ij,v_ijの内積 
	    
	abs_r[j] = sqrt(abs_r2[j]); //絶対値
	abs_v[j] = sqrt(abs_v2[j]);
	
	rij1 = 1.0/abs_r[j];
	
	E[i] += -0.5*G*ele[i].m*ele[j].m*rij1;  //エネルギー計算	
      }
    }  //j loop
    r1 = 1.0/r_c[i];
    E_tot += - G*M_0*ele[i].m*r1 + E[i];
  }  //i loop
  return E_tot;
}

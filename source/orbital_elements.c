#include "kuiper.h"

/*軌道要素計算*/
int Calculate_OrbitalElements(int i,int k,double x_c[][4],double v_c[][4],struct orbital_elements ele[],double P[][4],double Q[][4],double r_c[],double v2_c[],double r_dot_v[],double step,double t_sys,int N){
  
  double esin_u;
  double ecos_u;
  double sin_I;
  double cos_I;
  double sin_omega;
  double cos_omega;
  double sin_OMEGA;
  double cos_OMEGA;
  double radian;
  
  ele[i].axis = 1.0/(2.0/r_c[i] - v2_c[i]/(G*M_0));

#if SWAP
  if(ele[i].axis<=0.0){  //axisがマイナスになったとき=楕円軌道ではなくなった　swapする
    ele[i].judge = 0;
    FILE *fporbit;
    char orbit[100];
    sprintf(orbit,"%s%s.dat",STR(DIRECTORY),ele[i].name);
    fporbit = fopen(orbit,"a");
    if(fporbit==NULL){
      printf("orbit error\n");
      return -1;
    }
    fprintf(fporbit,"step=%e\tt=%e[yr]\t%s is dead\tswap with %s\n",step,t_sys,ele[i].name,ele[N].name);
    fclose(fporbit);
    printf("step=%e\tt=%e[yr]\t%s is dead\tswap with %s\n",step,t_sys,ele[i].name,ele[N].name);
    fclose(fporbit);
  }
  
#endif
     
  ele[i].e = sqrt((1.0-r_c[i]/ele[i].axis)*(1.0-r_c[i]/ele[i].axis) + r_dot_v[i]*r_dot_v[i]/(G*M_0*ele[i].axis));

  if(ele[i].e==0.0){
    ele[i].u = 0.0;
  }else{
    esin_u = r_dot_v[i]/sqrt(G*M_0*ele[i].axis);
    ecos_u = 1.0-r_c[i]/ele[i].axis;
    radian = atan2(esin_u,ecos_u);
    if(radian<0.0){
      ele[i].u = radian + 2.0*M_PI;
    }else{
      ele[i].u = radian;
    }
  }
  
  for(k=1;k<=3;++k){
    P[i][k] = x_c[i][k]*cos(ele[i].u)/r_c[i] - sqrt(ele[i].axis/(G*M_0))*v_c[i][k]*sin(ele[i].u);
    Q[i][k] = (x_c[i][k]*sin(ele[i].u)/r_c[i] + sqrt(ele[i].axis/(G*M_0))*v_c[i][k]*(cos(ele[i].u)-ele[i].e))/sqrt(1.0-ele[i].e);
  }
  
  sin_I = sqrt(P[i][3]*P[i][3] + Q[i][3]*Q[i][3]);
  cos_I = P[i][1]*Q[i][2] - P[i][2]*Q[i][1];
  radian = atan2(sin_I,cos_I);
  if(radian<0.0){
    ele[i].I = radian + 2.0*M_PI;
  }else{
    ele[i].I = radian;
  }
  
  sin_omega = P[i][3]/sin_I;
  cos_omega = Q[i][3]/sin_I;
  radian = atan2(sin_omega,cos_omega);
  if(radian<0.0){
    ele[i].omega = radian + 2.0*M_PI;
  }else{
    ele[i].omega = radian;
  }

  sin_OMEGA = (P[i][2]*Q[i][3] - Q[i][2]*P[i][3])/sin_I;
  cos_OMEGA = (P[i][1]*Q[i][3] - Q[i][1]*P[i][3])/sin_I;
  radian = atan2(sin_OMEGA,cos_OMEGA);
  if(radian<0.0){
    ele[i].OMEGA = radian + 2.0*M_PI;
  }else{
    ele[i].OMEGA = radian;
  }	  
  
  if(i<=Np){
    ele[i].R_H = ele[i].axis*cbrt(ele[i].m/((double)M_0)/3.0);
  }else{
    ele[i].R_H = 0.0;
  }
  

  return 0;
}

/*P計算*/
double Calculate_P(int i,int k,struct orbital_elements ele[]){
  if(k==1){
    return cos(ele[i].omega)*cos(ele[i].OMEGA) - sin(ele[i].omega)*sin(ele[i].OMEGA)*cos(ele[i].I);
  }else if(k==2){
    return cos(ele[i].omega)*sin(ele[i].OMEGA) + sin(ele[i].omega)*cos(ele[i].OMEGA)*cos(ele[i].I);
  }else{
    return sin(ele[i].omega)*sin(ele[i].I);
  }
}

/*Q計算*/
double Calculate_Q(int i,int k,struct orbital_elements ele[]){
  if(k==1){
    return -sin(ele[i].omega)*cos(ele[i].OMEGA) - cos(ele[i].omega)*sin(ele[i].OMEGA)*cos(ele[i].I);
  }else if(k==2){
    return -sin(ele[i].omega)*sin(ele[i].OMEGA) + cos(ele[i].omega)*cos(ele[i].OMEGA)*cos(ele[i].I);
  }else{
    return cos(ele[i].omega)*sin(ele[i].I);
  }
}



/*初期位置、速度計算*/
void InitialCondition(int i,int k,double P[][4],double Q[][4],double x_0[][4],double v_0[][4],double r_0[],struct orbital_elements ele[]){
  for(k=1;k<=3;k++){
    P[i][k] = Calculate_P(i,k,ele);
    Q[i][k] = Calculate_Q(i,k,ele);
	
    x_0[i][k] = ele[i].axis*P[i][k]*(cos(ele[i].u)-ele[i].e) + ele[i].axis*sqrt(1.0-ele[i].e*ele[i].e)*Q[i][k]*sin(ele[i].u);
  }
  //printf("x=%f\ty=%f\tz=%f\n",x_0[i][1],x_0[i][2],x_0[i][3]);
           
  r_0[i] = RadiusFromCenter(i,k,x_0,r_0);  //中心星からの距離


      
  for(k=1;k<=3;++k){
    v_0[i][k] = sqrt(G*M_0/ele[i].axis)/r_0[i]*(-ele[i].axis*P[i][k]*sin(ele[i].u) + ele[i].axis*sqrt(1.0-ele[i].e*ele[i].e)*Q[i][k]*cos(ele[i].u));
  }
  //printf("vx=%f\tvy=%f\tvz=%f\n",v_0[i][1],v_0[i][2],v_0[i][3]);
}

#include "kuiper.h"
 


int main(void){

  int N = Np + Nt;
  int i,i_sys,j,k,ite,interval,Ntemp;
  double t_sys;
  double t_ene[TIME_INTERVAL_MAX]={TIME_INTERVAL};

  double t_[N+1],dt_[N+1],Dt[N+1];
  double step=0.0;
  
  double x_0[N+1][4],r_0[N+1],v_0[N+1][4],v2_0[N+1];
  double x_p[N+1][4],r_p[N+1],v_p[N+1][4],v2_p[N+1];
  double x_c[N+1][4],r_c[N+1],v_c[N+1][4],v2_c[N+1];
  
  double abs_r[N+1],abs_r2[N+1],abs_v[N+1],abs_v2[N+1],r_dot_v_ij[N+1],r_dot_v[N+1];
  
  double a_0[N+1][4],adot_0[N+1][4];
  double a[N+1][4],adot[N+1][4],adot2_dt2[N+1][4],adot3_dt3[N+1][4];

  double abs_a[N+1],abs_adot[N+1],abs_adot2[N+1],abs_adot3[N+1];
  
  double E[N+1];
  double E_tot_0,E_tot;
  
  double abs_L_0,abs_L;
  
  double P[N+1][4],Q[N+1][4];
  
  double r_min_RH[Np+1],hill;

  struct orbital_elements ele[N+1];
  
  mkdir(STR(DIRECTORY), 0755);  //ディレクトリを作成  755 = rwxr-xr-x

  srand(RAND_SEED);

  //初期値

  /*
    Jupiter Mean Orbital Elements (J2000)
    
    Semimajor axis (AU)                  5.20336301  
    Orbital eccentricity                 0.04839266   
    Orbital inclination (deg)            1.30530  
    Longitude of ascending node (deg)  100.55615   
    Longitude of perihelion (deg)       14.75385   
    Mean Longitude (deg)                34.40438  
    Mass (g)                          1898.19e27
  */

  sprintf(ele[JUP].name,"Jupiter"); 
  ele[JUP].m_f = 1898.19/1988500.0;
  ele[JUP].e = 0.04839266;  //離心率  
  ele[JUP].axis_f = 5.20336301;  //長半径 
  ele[JUP].I = 1.30530*M_PI/180.0;  //軌道傾斜角
  ele[JUP].OMEGA = 100.55615*M_PI/180.0;  //昇交点経度
  ele[JUP].omega = 14.75385*M_PI/180.0;  //近日点引数
  ele[JUP].M = 34.40438*M_PI/180.0;  //平均近点離角
  ele[JUP].u = MeanLongitude_to_EccentricAnomaly(JUP,ele);  //離心近点離角
  ele[JUP].delta_axis = 0.2;
  

  /*
    Saturn Mean Orbital Elements (J2000)
    
    Semimajor axis (AU)                  9.53707032  
    Orbital eccentricity                 0.05415060   
    Orbital inclination (deg)            2.48446   
    Longitude of ascending node (deg)  113.71504   
    Longitude of perihelion (deg)       92.43194   
    Mean Longitude (deg)                49.94432
    Mass (g)                           568.34e27
  */

  sprintf(ele[SAT].name,"Saturn"); 
  ele[SAT].m_f = 568.34/1988500.0;
  ele[SAT].e = 0.05415060;  //離心率  
  ele[SAT].axis_f = 9.53707032;  //長半径 
  ele[SAT].I = 2.48446*M_PI/180.0;  //軌道傾斜角
  ele[SAT].OMEGA = 113.71504*M_PI/180.0;  //昇交点経度
  ele[SAT].omega = 92.43194*M_PI/180.0;  //近日点引数
  ele[SAT].M = 49.94432*M_PI/180.0;  //平均近点離角
  ele[SAT].u = MeanLongitude_to_EccentricAnomaly(SAT,ele);  //離心近点離角
  ele[SAT].delta_axis = -0.8;


  /*
    Uranus Mean Orbital Elements (J2000)
    
    Semimajor axis (AU)                 19.19126393  
    Orbital eccentricity                 0.04716771   
    Orbital inclination (deg)            0.76986   
    Longitude of ascending node (deg)   74.22988  
    Longitude of perihelion (deg)      170.96424  
    Mean Longitude (deg)               313.23218
    Mass (g)                            86.813e27
  */

  sprintf(ele[URA].name,"Uranus"); 
  ele[URA].m_f = 86.813/1988500.0;
  ele[URA].e = 0.04716771;  //離心率  
  ele[URA].axis_f = 19.19126393;  //長半径 
  ele[URA].I = 0.76986*M_PI/180.0;  //軌道傾斜角
  ele[URA].OMEGA = 74.22988*M_PI/180.0;  //昇交点経度
  ele[URA].omega = 170.96424*M_PI/180.0;  //近日点引数
  ele[URA].M = 313.23218*M_PI/180.0;  //平均近点離角
  ele[URA].u = MeanLongitude_to_EccentricAnomaly(URA,ele);  //離心近点離角
  ele[URA].delta_axis = -3.0;


  /*
    Neptune Mean Orbital Elements (J2000)
    
    Semimajor axis (AU)                 30.06896348  
    Orbital eccentricity                 0.00858587   
    Orbital inclination (deg)            1.76917  
    Longitude of ascending node (deg)  131.72169   
    Longitude of perihelion (deg)       44.97135  
    Mean Longitude (deg)               304.88003
    Mass (g)                           102.413e27
  */

  sprintf(ele[NEP].name,"Neptune"); 
  ele[NEP].m_f = 102.413/1988500.0;
  ele[NEP].e = 0.0858587;  //離心率  
  ele[NEP].axis_f = 30.06896348;  //長半径 
  ele[NEP].I = 1.76917*M_PI/180.0;  //軌道傾斜角
  ele[NEP].OMEGA = 131.72169*M_PI/180.0;  //昇交点経度
  ele[NEP].omega = 44.97135*M_PI/180.0;  //近日点引数
  ele[NEP].M = 304.88003*M_PI/180.0;  //平均近点離角
  ele[NEP].u = MeanLongitude_to_EccentricAnomaly(NEP,ele);  //離心近点離角
  ele[NEP].delta_axis = -7.0;


  for(i=1;i<=Np;++i){
    //ele[i].u = ((double)rand())/((double)RAND_MAX+1.0)*2*M_PI;
    ele[i].m = 0.0;
    ele[i].axis_0 = ele[i].axis_f + ele[i].delta_axis; 
    ele[i].axis = ele[i].axis_0;
    ele[i].R_H = ele[i].axis*cbrt(ele[i].m/((double)M_0)/3.0);
    ele[i].judge = 1;
    ele[i].orinum = i;
  }


  for(i=Np+1;i<=N;++i){   
 
    sprintf(ele[i].name,"test_particle%03d",i-Np);
    ele[i].m = 0.0;
    ele[i].e = ((double)rand())/((double)RAND_MAX+1.0)*UPPER_ECC;  //離心率
    ele[i].axis = ((double)rand())/((double)RAND_MAX+1.0)*(OUTER_AXIS-INNER_AXIS)+INNER_AXIS ;  //長半径
    ele[i].u = ((double)rand())/((double)RAND_MAX+1.0)*2.0*M_PI;  //離心近点離角
    ele[i].I = ((double)rand())/((double)RAND_MAX+1.0)*UPPER_ICC;  //軌道傾斜角
    ele[i].OMEGA = ((double)rand())/((double)RAND_MAX+1.0)*2.0*M_PI;  //昇交点経度
    ele[i].omega = ((double)rand())/((double)RAND_MAX+1.0)*2.0*M_PI;  //近日点引数
    ele[i].R_H = 0.0;
    ele[i].judge = 1;
    ele[i].orinum = i;
  }


  for(i=1;i<=N;++i){  
    InitialCondition(i,k,P,Q,x_0,v_0,r_0,ele);  //初期位置、速度
         
    printf("%s\tx_0[%d][1]=%f\tx_0[%d][2]=%f\tx_0[%d][3]=%f\n",ele[i].name,i,x_0[i][1],i,x_0[i][2],i,x_0[i][3]);
    printf("%s\tv_0[%d][1]=%f\tv_0[%d][2]=%f\tv_0[%d][3]=%f\n",ele[i].name,i,v_0[i][1],i,v_0[i][2],i,v_0[i][3]);
  }  //i loop
    
    
    
    
#if ORBITALELEMENTS_FILE
  FILE *fporbit;   //初期の軌道要素をファイルへ書き出し
  char orbit[100];
  for(i=1;i<=N;++i){
    sprintf(orbit,"%s%s.dat",STR(DIRECTORY),ele[i].name);
    fporbit = fopen(orbit,"w");
    if(fporbit==NULL){
      printf("orbit 0 error\n");
      return -1;
    }

    fprintf(fporbit,"#t\te\ta\tu\tI\tOMEGA\tomega\tR_H\n");
    fprintf(fporbit,"%e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",0.0,ele[i].e,ele[i].axis,ele[i].u,ele[i].I,ele[i].OMEGA,ele[i].omega,ele[i].R_H);
  
    fclose(fporbit);
  }
#endif
    


  for(i=1;i<=N;++i){
    v2_0[i] = SquareOfVelocity(i,k,v_0,v2_0);  //速度の2乗
  }


  E_tot_0 = 0.0;  
  E_tot_0 = Calculate_Energy(i,j,ele,x_0,v_0,v2_0,r_0,abs_r,abs_r2,abs_v,abs_v2,r_dot_v_ij,E,E_tot_0,N);  //初期エネルギー,その他  コメントアウトしちゃだめなやつ
  //printf("%e\t%.15e\n",0.0,E_tot_0);

#if ENERGY_FILE
  FILE *fpEne;   //初期エネルギーをファイルへ書き出し
  char Ene[30];
  sprintf(Ene,"%sENERGY.dat",STR(DIRECTORY));
  fpEne = fopen(Ene,"w");
  if(fpEne==NULL){
    printf("Ene 0 error\n");
    return -1;
  }
  fprintf(fpEne,"#t\tE_tot\tE error\n");
  fprintf(fpEne,"%e\t%.15e\t%e\n",0.0,E_tot_0,0.0);
  fclose(fpEne);

  abs_L_0 = AngularMomentum(i,k,ele,x_0,v_0,abs_L_0,N);  //角運動量の大きさ
  //printf("abs_L_0=%.15e\n",abs_L_0);
#endif
      
 
  for(i=1;i<=N;++i){     
    r_dot_v[i] = InnerProduct(i,k,x_0,v_0,r_dot_v);  //r_i,v_iの内積
  }

 
  for(i=1;i<=Np;++i){  //惑星の加速度
    for(k=1;k<=3;++k){
  
#if INTERACTION
    a_0[i][k] = All_Acceleration(i,j,k,ele,x_0,r_0,abs_r2,a_0);  //初期の加速度
    adot_0[i][k] = All_dAcceleration(i,j,k,ele,x_0,v_0,r_dot_v,r_dot_v_ij,r_0,abs_r2,adot_0);
#else
    a_0[i][k] = External_Acceleration(i,k,x_0,r_0);
    adot_0[i][k] = External_dAcceleration(i,k,x_0,v_0,r_0,r_dot_v);
#endif
    a_0[i][k] += v_0[i][k]/sqrt(v2_0[i])/TAU_MOVE*(sqrt(G*M_0/ele[i].axis_0) - sqrt(G*M_0/ele[i].axis_f));
  
    //printf("a_0[%d][%d]=%f\tadot_0[%d][%d]=%f\n",i,k,a_0[i][k],i,k,adot_0[i][k]);
    }  //k loop
  }  //i loop
  
  
  for(i=Np+1;i<=N;++i){  //test particle の加速度
    for(k=1;k<=3;++k){
    a_0[i][k] = All_Acceleration(i,j,k,ele,x_0,r_0,abs_r2,a_0);  //初期の加速度
    adot_0[i][k] = All_dAcceleration(i,j,k,ele,x_0,v_0,r_dot_v,r_dot_v_ij,r_0,abs_r2,adot_0);
    //printf("a_0[%d][%d]=%f\tadot_0[%d][%d]=%f\n",i,k,a_0[i][k],i,k,adot_0[i][k]);
    }
  }
  
#if R_MIN
  for(i=1;i<=Np;i++){
    r_min_RH[i] = 10000.0;  //適当な距離[RH]
  }
#endif


  for(i=1;i<=N;++i){ 
    dt_[i] = Timestep_i_0(i,k,a_0,adot_0,abs_a,abs_adot,dt_);  //初期のタイムステップ計算
    //printf("initial dt_[%d]=%e\n",i,dt_[i]);
  }
    


  for(i=1;i<=N;++i){  
    if(i==1){
      t_sys = dt_[1];
      i_sys = 1;
    }else if(dt_[i] < t_sys){
      t_sys = dt_[i];  //dt_i が最小のものを選ぶ
      i_sys = i;  //i_sysを選ぶ
    } 
  }


  for(i=1;i<=N;++i){
    t_[i] = 0.0;
  }



  //printf("-----\n");
    
  ////////////////////////////////////////////////////////////////////////////////////

  interval = 0;

  while(TIME_INTERVAL_MAX>interval){
    
    if(t_sys<t_ene[interval]){  
      //individual timestep

      for(i=1;i<=Np;++i){  //惑星の質量を線形的に増加 10TAU = 5My
	if(t_sys<=TAU_MASS){
	  ele[i].m = ele[i].m_f*t_sys/TAU_MASS;
	}else{
	  ele[i].m = ele[i].m_f;
	}
      }

      for(i=1;i<=N;++i){ 
	Dt[i] = t_sys - t_[i]; 
	Predictor(i,k,x_0,v_0,a_0,adot_0,x_p,v_p,r_p,v2_p,Dt);  //予測子 t_sysにおけるすべての粒子を計算
      }      
	
      for(i=1;i<=N;++i){
	if(i==i_sys){  
	  Corrector_sys(i_sys,j,k,ele,x_p,v_p,r_p,v2_p,x_c,v_c,r_c,v2_c,a_0,adot_0,a,adot,adot2_dt2,adot3_dt3,abs_r,abs_r2,abs_v,abs_v2,r_dot_v_ij,r_dot_v,t_sys,dt_);  //修正子 i_sys のみ
	}else{
	  for(k=1;k<=3;++k){  //i_sys 以外の粒子は予測子を使う
	    x_c[i][k] = x_p[i][k];
	    v_c[i][k] = v_p[i][k];
	  }
	  r_c[i] = RadiusFromCenter(i,k,x_c,r_c);  //中心からの距離
	  v2_c[i] = SquareOfVelocity(i,k,v_c,v2_c);  //速度の2乗
	}
      }  //i loop
	
	     
      for(ite=1;ite<=ITE_MAX;++ite){  //iteration 3回 
	Iteration_sys(i_sys,j,k,ele,x_p,v_p,x_c,v_c,r_c,v2_c,a_0,adot_0,a,adot,adot2_dt2,adot3_dt3,abs_r,abs_r2,abs_v,abs_v2,r_dot_v_ij,r_dot_v,t_sys,dt_);  // i_sys のみ	
      }
	

#if R_MIN
      for(i=1;i<=Np;i++){
	ele[i].R_H = ele[i].axis*cbrt(ele[i].m/((double)M_0)/3.0);
	hill = ele[i].R_H;
	for(j=Np+1;j<=N;++j){
	  if(i!=j){
	    abs_r2[j] = SquareOfRaletiveDistance(i,j,x_c); //絶対値2乗
	    abs_r[j] = sqrt(abs_r2[j]); //絶対値
#endif
	    
#if R_MIN * SWAP
	    if(abs_r[j]<=hill){
	      ele[j].judge = 0;
	      abs_r[j] = abs_r[Np+1];  //r_minをリセットする用　Np+1の粒子の値を使う
	      
	      sprintf(orbit,"%s%s.dat",STR(DIRECTORY),ele[j].name);
	      fporbit = fopen(orbit,"a");
	      if(fporbit==NULL){
		printf("orbit error\n");
		return -1;
	      }
	      fprintf(fporbit,"step=%e\tt=%e[yr]\t%s is dead\tswap with %s\n",step,t_sys,ele[j].name,ele[N].name);
	      fclose(fporbit);
	      
	      printf("step=%e\tt=%e[yr]\t%s is dead\tswap with %s\n",step,t_sys,ele[j].name,ele[N].name);
	      continue;  //以下の処理をスキップ
	    }
#endif
	  
#if R_MIN
	    if(abs_r[j] < r_min_RH[i]*hill){  //最近接距離計算 
	      r_min_RH[i] = abs_r[j]/hill;  //Hill半径で規格化
	    } 
	  }
	}  
      }       
#endif

	
    }else{  
      //t_ene[interval] ですべての粒子をそろえ、エネルギー、軌道要素等計算
      
      t_sys = t_ene[interval];
      
      for(i=1;i<=Np;++i){  //惑星の質量を線形的に増加 10TAU = 5Myr
	if(t_sys<=TAU_MASS){
	  ele[i].m = ele[i].m_f*t_sys/TAU_MASS;
	}else{
	  ele[i].m = ele[i].m_f;
	}
      }
      
      
      for(i=1;i<=N;++i){
	Dt[i] = t_ene[interval] - t_[i];
	dt_[i] = Dt[i];
      }

      for(i=1;i<=N;++i){ 	  
	Predictor(i,k,x_0,v_0,a_0,adot_0,x_p,v_p,r_p,v2_p,Dt);  //予測子 t_sysにおけるすべての粒子を計算
      }
	
      for(i=1;i<=N;++i){
	Corrector_sys(i,j,k,ele,x_p,v_p,r_p,v2_p,x_c,v_c,r_c,v2_c,a_0,adot_0,a,adot,adot2_dt2,adot3_dt3,abs_r,abs_r2,abs_v,abs_v2,r_dot_v_ij,r_dot_v,t_sys,dt_);
      }

      for(ite=1;ite<=ITE_MAX;++ite){  //iteration 3回 
	for(i=1;i<=N;++i){	  
	  Iteration_sys(i,j,k,ele,x_p,v_p,x_c,v_c,r_c,v2_c,a_0,adot_0,a,adot,adot2_dt2,adot3_dt3,abs_r,abs_r2,abs_v,abs_v2,r_dot_v_ij,r_dot_v,t_sys,dt_);
	}
      }


#if ENERGY_FILE

      E_tot = 0.0;
      E_tot = Calculate_Energy(i,j,ele,x_c,v_c,v2_c,r_c,abs_r,abs_r2,abs_v,abs_v2,r_dot_v_ij,E,E_tot,N);  //エネルギー計算
	
      sprintf(Ene,"%sENERGY.dat",STR(DIRECTORY));
      fpEne = fopen(Ene,"a");
      if(fpEne==NULL){
	printf("Ene error\n");
	return -1;
      }
      fprintf(fpEne,"%e\t%.15e\t%.15e\n",t_sys,E_tot,(E_tot-E_tot_0)/fabs(E_tot_0));
      fclose(fpEne);

      abs_L = AngularMomentum(i,k,ele,x_c,v_c,abs_L,N);
#endif
	


#if R_MIN
      for(i=1;i<=Np;i++){
	ele[i].R_H = ele[i].axis*cbrt(ele[i].m/((double)M_0)/3.0);
	hill = ele[i].R_H;
	for(j=Np+1;j<=N;++j){
	  if(i!=j){
	    abs_r2[j] = SquareOfRaletiveDistance(i,j,x_c); //絶対値2乗
	    abs_r[j] = sqrt(abs_r2[j]); //絶対値
#endif
	    
#if R_MIN * SWAP
	    if(abs_r[j]<=hill){
	      ele[j].judge = 0;
	      abs_r[j] = abs_r[Np+1];  //r_minをリセットする用　Np+1の粒子の値を使う
	      
	      sprintf(orbit,"%s%s.dat",STR(DIRECTORY),ele[j].name);
	      fporbit = fopen(orbit,"a");
	      if(fporbit==NULL){
		printf("orbit error\n");
		return -1;
	      }
	      fprintf(fporbit,"step=%e\tt=%e[yr]\t%s is dead\tswap with %s\n",step,t_sys,ele[j].name,ele[N].name);
	      fclose(fporbit);
	      
	      printf("step=%e\tt=%e[yr]\t%s is dead\tswap with %s\n",step,t_sys,ele[j].name,ele[N].name);
	      continue;  //以下の処理をスキップ
	    }
#endif
	  
#if R_MIN
	    if(abs_r[j] < r_min_RH[i]*hill){  //最近接距離計算 
	      r_min_RH[i] = abs_r[j]/hill;  //Hill半径で規格化
	    } 
	  }
	}  
      }       
#endif

      
      
      
#if ORBITALELEMENTS_FILE
      
      for(i=1;i<=N;++i){
	Calculate_OrbitalElements(i,k,x_c,v_c,ele,P,Q,r_c,v2_c,r_dot_v,step,t_sys,N);  //軌道要素計算  ファイルへ書き出し
	if(ele[i].judge == 0){
	  continue;  //以下の処理をスキップ
	}
	sprintf(orbit,"%s%s.dat",STR(DIRECTORY),ele[i].name);
	fporbit = fopen(orbit,"a");
	if(fporbit==NULL){
	  printf("orbit error\n");
	  return -1;
	}
	
	fprintf(fporbit,"%e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",t_sys,ele[i].e,ele[i].axis,ele[i].u,ele[i].I,ele[i].OMEGA,ele[i].omega,ele[i].R_H);
	
	fclose(fporbit);
      }
#endif
      
      interval++;
      
    }
    
    /*for(i=1;i<=N;++i){
      printf("%s\tx_c[%d][1]=%f\tx_c[%d][2]=%f\tx_c[%d][3]=%f\n",ele[i].name,i,x_c[i][1],i,x_c[i][2],i,x_c[i][3]);
      printf("%s\tv_c[%d][1]=%f\tv_c[%d][2]=%f\tv_c[%d][3]=%f\n",ele[i].name,i,v_c[i][1],i,v_c[i][2],i,v_c[i][3]);
      }*/
    
    
#if SWAP
    //消したい粒子を交換
    //0番目の要素はコピーに使うだけ
    Ntemp = N;
    for(i=Ntemp;i>=Np+1;--i){
      if(ele[i].judge==0){	
	ele[0] = ele[i];
	ele[i] = ele[N];
	ele[N] = ele[0];
	
	t_[0] = t_[i];
	t_[i] = t_[N];
	t_[N] = t_[0];
	
	dt_[0] = dt_[i];
	dt_[i] = dt_[N];
	dt_[N] = dt_[0];
	
	for(k=1;k<=3;++k){
	  x_0[0][k] = x_0[i][k];
	  x_0[i][k] = x_0[N][k];
	  x_0[N][k] = x_0[0][k];
	  
	  x_c[0][k] = x_c[i][k];
	  x_c[i][k] = x_c[N][k];
	  x_c[N][k] = x_c[0][k];
	  
	  v_0[0][k] = v_0[i][k];
	  v_0[i][k] = v_0[N][k];
	  v_0[N][k] = v_0[0][k];
	  
	  v_c[0][k] = v_c[i][k];
	  v_c[i][k] = v_c[N][k];
	  v_c[N][k] = v_c[0][k];
	  
	  a_0[0][k] = a_0[i][k];
	  a_0[i][k] = a_0[N][k];
	  a_0[N][k] = a_0[0][k];
	  
	  a[0][k] = a[i][k];
	  a[i][k] = a[N][k];
	  a[N][k] = a[0][k];
	  
	  adot_0[0][k] = adot_0[i][k];
	  adot_0[i][k] = adot_0[N][k];
	  adot_0[N][k] = adot_0[0][k];
	  
	  adot[0][k] = adot[i][k];
	  adot[i][k] = adot[N][k];
	  adot[N][k] = adot[0][k];
	  
	  adot2_dt2[0][k] = adot2_dt2[i][k];
	  adot2_dt2[i][k] = adot2_dt2[N][k];
	  adot2_dt2[N][k] = adot2_dt2[0][k];
	  
	  adot3_dt3[0][k] = adot3_dt3[i][k];
	  adot3_dt3[i][k] = adot3_dt3[N][k];
	  adot3_dt3[N][k] = adot3_dt3[0][k];	
	}
	N--;
      }  //if
    }  //i loop
    
#endif
    
    
    
    if(fmod(step,1.0E5)==0.0){
      //printf("i_sys=%03d\tt=%.15e\tE=%.15e\tL=%.15e\tr_min=%.15e\n",i_sys,t_sys,E_tot,abs_L,r_min);  //全エネルギー,全角運動量
      printf("step=%e\tN=%d\ti_sys=%03d\tt=%.2e[yr]",step,N,i_sys,t_sys/2.0/M_PI);
      for(i=1;i<=Np;i++){
        printf("\tr_min[%d]=%.6f[RH]",i,r_min_RH[i]);
      }
      printf("\n");
    }
    
    
    t_[i_sys] += dt_[i_sys];  //i_sys のみ時間更新
    
    
    

    dt_[i_sys] = Timestep_i_sys(i_sys,k,a,adot,adot2_dt2,adot3_dt3,abs_a,abs_adot,abs_adot2,abs_adot3,dt_);  //i_sys のみタイムステップ計算
    
    
    
    for(k=1;k<=3;++k){  //i_sys のみ更新
      x_0[i_sys][k] = x_c[i_sys][k];
      v_0[i_sys][k] = v_c[i_sys][k];
      a_0[i_sys][k] = a[i_sys][k];
      adot_0[i_sys][k] = adot[i_sys][k];
    }
    
    
    t_sys = t_[1] + dt_[1];
    i_sys = 1;
    for(i=2;i<=N;++i){  
      if((t_[i] + dt_[i]) < t_sys){
	t_sys = t_[i] + dt_[i];  //dt_i が最小のものを選ぶ
	i_sys = i;  //i_sysを選ぶ
      }
    }
    
    
    
    step+=1.0;
  }  //t loop
  
#if SWAP
  printf("dead\n");
  for(i=1;i<=Np+Nt;i++){
    if(ele[i].judge==0){
      printf("[%d] ",ele[i].orinum);
    }
  }
  printf("\n");
#endif
  
#if ENERGY_FILE
  printf("dt_[i_sys]=%e\tE_error=%.15e\tL_error=%.15e\n",dt_[i_sys],(E_tot-E_tot_0)/fabs(E_tot_0),(abs_L-abs_L_0)/abs_L_0);
#endif
  
  printf("step=%e\n",step);
  
  for(i=0;i<TIME_INTERVAL_MAX;i++){
    printf("t_ene[%d]=%e\n",i,t_ene[i]);
  } 
}

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#define GYOU 2213  //行数
#define N 1000
#define BIN_MIN 25
#define BIN_MAX 70
#define Max 3


struct orbital_elements{
  char number[10];
  char name[30];
  double t;
  double omega;
  double Omega;
  double inc;
  double ecc;
  double axis;
};


int main(void){
  
  int i,sum,sumtotal,j;
  double bin,delta=0.3; //0.3AU刻み
  struct orbital_elements ele[N+1];

  char buf[200];
  FILE *fp;
  char filename[50];
    
  
    
  for(i=1;i<=100;i++){
    sprintf(filename,"./N100_1/test_particle%03d.dat",i);
    if((fp = fopen(filename,"r")) == NULL){
      fprintf(stderr,"Can't open %s.\n",filename);
      exit(-1);
    }
      
    fgets(buf,sizeof(buf),fp);  //読み飛ばし
    for(j=1;j<=Max;j++){
      fgets(buf,sizeof(buf),fp);
    }
      
    if(fgets(buf,sizeof(buf),fp) != NULL){
      sscanf(buf,"%lf\t%*lf\t%lf\t%*lf\t%*lf\t%*lf\t%*lf\t%*lf",&ele[i].t,&ele[i].axis);
    }/*else{
       continue;
      }*/
    
    //printf("i=%d\taxis=%lf\n",i,ele[i].axis);
    
    
    fclose(fp);
  }
  
  for(i=101;i<=200;i++){
    sprintf(filename,"./N100_2/test_particle%03d.dat",i-100);
    if((fp = fopen(filename,"r")) == NULL){
      fprintf(stderr,"Can't open %s.\n",filename);
      exit(-1);
    }
    
    fgets(buf,sizeof(buf),fp);  //読み飛ばし
    for(j=1;j<=Max;j++){
      fgets(buf,sizeof(buf),fp);
    }
    
    if(fgets(buf,sizeof(buf),fp) != NULL){
      sscanf(buf,"%*lf\t%*lf\t%lf\t%*lf\t%*lf\t%*lf\t%*lf\t%*lf",&ele[i].axis);
    }/*else{
      continue;
      }*/
    
    //printf("i=%d\taxis=%lf\n",i,ele[i].axis);
    
    
    fclose(fp);
  }

  for(i=201;i<=300;i++){
    sprintf(filename,"./N100_3/test_particle%03d.dat",i-200);
    if((fp = fopen(filename,"r")) == NULL){
      fprintf(stderr,"Can't open %s.\n",filename);
      exit(-1);
    }
    
    fgets(buf,sizeof(buf),fp);  //読み飛ばし
    for(j=1;j<=Max;j++){
      fgets(buf,sizeof(buf),fp);
    }
    
    if(fgets(buf,sizeof(buf),fp) != NULL){
      sscanf(buf,"%*lf\t%*lf\t%lf\t%*lf\t%*lf\t%*lf\t%*lf\t%*lf",&ele[i].axis);
    }/*else{
      continue;
      }*/
    
    //printf("i=%d\taxis=%lf\n",i,ele[i].axis);
    
    fclose(fp);
  }

  for(i=301;i<=400;i++){
    sprintf(filename,"./N100_4/test_particle%03d.dat",i-300);
    if((fp = fopen(filename,"r")) == NULL){
      fprintf(stderr,"Can't open %s.\n",filename);
      exit(-1);
    }
    
    fgets(buf,sizeof(buf),fp);  //読み飛ばし
    for(j=1;j<=Max;j++){
      fgets(buf,sizeof(buf),fp);
    }
    
    if(fgets(buf,sizeof(buf),fp) != NULL){
      sscanf(buf,"%*lf\t%*lf\t%lf\t%*lf\t%*lf\t%*lf\t%*lf\t%*lf",&ele[i].axis);
    }/*else{
      continue;
      }*/
    
    //printf("i=%d\taxis=%lf\n",i,ele[i].axis);
    
    fclose(fp);
  }
  
  for(i=401;i<=500;i++){
    sprintf(filename,"./N100_5/test_particle%03d.dat",i-400);
    if((fp = fopen(filename,"r")) == NULL){
      fprintf(stderr,"Can't open %s.\n",filename);
      exit(-1);
    }
    
    fgets(buf,sizeof(buf),fp);  //読み飛ばし
    for(j=1;j<=Max;j++){
      fgets(buf,sizeof(buf),fp);
    }
    
    if(fgets(buf,sizeof(buf),fp) != NULL){
      sscanf(buf,"%*lf\t%*lf\t%lf\t%*lf\t%*lf\t%*lf\t%*lf\t%*lf",&ele[i].axis);
    }/*else{
       continue;
       }*/
    
    //printf("i=%d\taxis=%lf\n",i,ele[i].axis);
    
    fclose(fp);
  }
  
  for(i=501;i<=600;i++){
    sprintf(filename,"./N100_6/test_particle%03d.dat",i-500);
    if((fp = fopen(filename,"r")) == NULL){
      fprintf(stderr,"Can't open %s.\n",filename);
      exit(-1);
    }
    
    fgets(buf,sizeof(buf),fp);  //読み飛ばし
    for(j=1;j<=Max;j++){
      fgets(buf,sizeof(buf),fp);
    }
    
    if(fgets(buf,sizeof(buf),fp) != NULL){
      sscanf(buf,"%*lf\t%*lf\t%lf\t%*lf\t%*lf\t%*lf\t%*lf\t%*lf",&ele[i].axis);
    }/*else{
      continue;
      }*/
    
    //printf("i=%d\taxis=%lf\n",i,ele[i].axis);
    
    fclose(fp);
  }

  for(i=601;i<=700;i++){
    sprintf(filename,"./N100_7/test_particle%03d.dat",i-600);
    if((fp = fopen(filename,"r")) == NULL){
      fprintf(stderr,"Can't open %s.\n",filename);
      exit(-1);
    }
    
    fgets(buf,sizeof(buf),fp);  //読み飛ばし
    for(j=1;j<=Max;j++){
      fgets(buf,sizeof(buf),fp);
    }
    
    if(fgets(buf,sizeof(buf),fp) != NULL){
      sscanf(buf,"%*lf\t%*lf\t%lf\t%*lf\t%*lf\t%*lf\t%*lf\t%*lf",&ele[i].axis);
    }/*else{
      continue;
      }*/
    
    //printf("i=%d\taxis=%lf\n",i,ele[i].axis);
    
    fclose(fp);
  }

  for(i=701;i<=800;i++){
    sprintf(filename,"./N100_8/test_particle%03d.dat",i-700);
    if((fp = fopen(filename,"r")) == NULL){
      fprintf(stderr,"Can't open %s.\n",filename);
      exit(-1);
    }
    
    fgets(buf,sizeof(buf),fp);  //読み飛ばし
    for(j=1;j<=Max;j++){
      fgets(buf,sizeof(buf),fp);
    }
    
    if(fgets(buf,sizeof(buf),fp) != NULL){
      sscanf(buf,"%*lf\t%*lf\t%lf\t%*lf\t%*lf\t%*lf\t%*lf\t%*lf",&ele[i].axis);
    }/*else{
      continue;
      }*/
    
    //printf("i=%d\taxis=%lf\n",i,ele[i].axis);
    
    fclose(fp);
  }

  for(i=801;i<=900;i++){
    sprintf(filename,"./N100_9/test_particle%03d.dat",i-800);
    if((fp = fopen(filename,"r")) == NULL){
      fprintf(stderr,"Can't open %s.\n",filename);
      exit(-1);
    }
    
    fgets(buf,sizeof(buf),fp);  //読み飛ばし
    for(j=1;j<=Max;j++){
      fgets(buf,sizeof(buf),fp);
    }
    
    if(fgets(buf,sizeof(buf),fp) != NULL){
      sscanf(buf,"%*lf\t%*lf\t%lf\t%*lf\t%*lf\t%*lf\t%*lf\t%*lf",&ele[i].axis);
    }/*else{
      continue;
      }*/
    
    //printf("i=%d\taxis=%lf\n",i,ele[i].axis);
    
    fclose(fp);
  }

  for(i=901;i<=1000;i++){
    sprintf(filename,"./N100_10/test_particle%03d.dat",i-900);
    if((fp = fopen(filename,"r")) == NULL){
      fprintf(stderr,"Can't open %s.\n",filename);
      exit(-1);
    }
    
    fgets(buf,sizeof(buf),fp);  //読み飛ばし
    for(j=1;j<=Max;j++){
      fgets(buf,sizeof(buf),fp);
    }
    
    if(fgets(buf,sizeof(buf),fp) != NULL){
      sscanf(buf,"%*lf\t%*lf\t%lf\t%*lf\t%*lf\t%*lf\t%*lf\t%*lf",&ele[i].axis);
    }/*else{
      continue;
      }*/
    
    //printf("i=%d\taxis=%lf\n",i,ele[i].axis);
    
    fclose(fp);
  }




  ////////////////////////////////////////////////////////////////////////
    
  bin=BIN_MIN;
  sumtotal = 0;
  while(bin<=BIN_MAX){
    sum=0;
    for(i=1;i<=N;i++){
      if((ele[i].axis>=bin) && (ele[i].axis<bin+delta)){
	sum++;
	sumtotal++;
      }
    }
    printf("%f\t%d\n",bin+delta/2.0,sum);
    //printf("%f\t%d\n",bin,sum);
    bin+=delta;
  }
    
  printf("time=%e\tsumtotal=%d\n",ele[1].t,sumtotal);
  

}

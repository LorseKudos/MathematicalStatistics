#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define T_MIN 0.0 // 積分範囲の最小値
#define T_MAX 1.0 // 積分範囲の最大値

double function (double ,double ,double); // 被積分関数 1つ目の変数はx,2つ目の変数はパラメータa,3つ目の変数はパラメータb
double simpson (double,double,double); //シンプソン則 1つ目の変数は積分区間の上限,2つ目の変数はパラメータa,3つ目の変数はパラメータb
double beta_func(double,double,double); //ベータ関数 1つ目の変数はx,2つ目の変数はパラメータa,3つ目の変数はパラメータb
double Pr(double ,double); //累積分布関数Pr[T<=t] 一つ目の変数がt,二つ目の変数は自由度

int main(void){
  FILE *fp;
  char *fname = "Energy_GDP.csv";
  int NUMBER = 136; //データ数
  double x[NUMBER],y[NUMBER]; //xがGDPのlog10,yがEnergyのlog10
  double ave_x = 0,ave_y = 0; //x,yの平均
  double a_hat,b_hat; //回帰係数alpha,beta
  double denomi_b_hat = 0,nume_b_hat = 0; // //回帰係数betaの分母,分子
  int k = 0; //ループ変数
  double variance = 0; //推定誤差の不偏分散
  double sum_pow_x = 0; //(x_i-ave_x)の二乗和
  double t; //t値
  double p; //p値
  double alpha = 0.05; //有意水準

  double energy,gdp;

  fp = fopen(fname,"r");
  while(!feof(fp)){
    //データ読み込み
    fscanf(fp,"%lf,%lf\n",&energy,&gdp);
    x[k] = log10(gdp);
    y[k] = log10(energy);
    ave_x += x[k];
    ave_y += y[k];
    k++;
  }
  fclose(fp);

  ave_x = ave_x / NUMBER;
  ave_y = ave_y / NUMBER;

  for(k=0;k<NUMBER;k++){
    nume_b_hat += (y[k] - ave_y)*(x[k] - ave_x);
    denomi_b_hat += pow(x[k] - ave_x,2);
  }

  b_hat = nume_b_hat / denomi_b_hat;
  a_hat = ave_y - b_hat * ave_x;

  printf("alpha = %.20lf,beta = %lf\n",a_hat,b_hat);

  for(k=0;k<NUMBER;k++){
    variance += pow(y[k] - (a_hat + b_hat * x[k]),2);
    sum_pow_x += pow(x[k] - ave_x,2);
  }

  variance /= NUMBER - 2;
  printf("variance = %lf\n",variance);

  t = b_hat / sqrt(variance/sum_pow_x);

  p = 1 - Pr(t,NUMBER-2);

  printf("t = %lf,p = %.20lf\n",t,p);

  if(p <= alpha){
    printf("pが有意水準%lf以下なのでH0は棄却される\n",alpha);
  }else{
    printf("pが有意水準%lfより大きいのでH0は棄却されない\n",alpha);
  }

  return 0;
}

double function (double t,double a,double b){
  double result;
  result = pow(t,a-1) * pow(1-t,b-1);
  return result;
}

double simpson(double x,double a,double b){
  double integral; //積分結果
  int div_Num; // 分割数
  double h = 0.000004; //刻み幅
  double t;
  int i;
  div_Num = (x - T_MIN) / h;
  t = T_MIN;
  integral = function(0,a,b);
  t += h;

  if(x == 0) return 0;

  for (i=1; i<div_Num; i++) {
    integral += (3 + pow(-1,i+1))*function(t,a,b);
    t += h;
  }

  integral += function(x,a,b);
  integral *= h/3.0;

  return integral;
}

double beta_func(double x,double a,double b){
  return simpson(x,a,b)/simpson(1.0,a,b);
}

double Pr(double t,double nu){
  return beta_func((t+sqrt(t*t+nu))/(2*sqrt(t*t+nu)),nu/2.0,nu/2.0);
}
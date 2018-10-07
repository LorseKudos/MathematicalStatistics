#include <stdio.h>
#include <math.h>

#define T_MIN 0.0 // 積分範囲の最小値
#define T_MAX 1.0 // 積分範囲の最大値

double function (double ,double ,double); // 被積分関数 1つ目の変数はx,2つ目の変数はパラメータa,3つ目の変数はパラメータb
double simpson (double,double,double); //シンプソン則 1つ目の変数は積分区間の上限,2つ目の変数はパラメータa,3つ目の変数はパラメータb
double beta_func(double,double,double); //ベータ関数 1つ目の変数はx,2つ目の変数はパラメータa,3つ目の変数はパラメータb
double Pr(double ,double); //累積分布関数Pr[T<=t] 一つ目の変数がt,二つ目の変数は自由度

int main(void){
  unsigned int I = 2312; //乗算法の初期値
  unsigned int a = 69621; //乗算法のパラメータ
  double U_1,U_2; //一様乱数
  double N_1,N_2; //標準正規乱数
  double mu = 4; //平均
  double sigma = 3; //標準偏差
  double sum = 0,pow_sum = 0; //和、二乗和
  double ave,pow_ave,variance; //平均,二乗平均,不偏分散
  int k; //ループ変数
  int n = 100; //乱数生成個数
  double alpha = 0.05; //有意水準
  double t; //t値
  double p; //p値

  for(k=0;k<n;k++){
    I *= a;
    U_1 =  ((double) I) / (pow(2,32) - 1) ; //一様乱数
    I *= a;
    U_2 =  ((double) I) / (pow(2,32) - 1) ; //一様乱数
    N_1 = sigma * pow(-2.0 * log(U_1),0.5) * cos(2.0 * M_PI * U_2) + mu; //平均mu,分散sigma^2の正規乱数
    N_2 = sigma * pow(-2.0 * log(U_1),0.5) * sin(2.0 * M_PI * U_2) + mu;

    sum += N_1;
    pow_sum += N_1*N_1;
  }

  ave = sum / n;
  pow_ave = pow_sum / n;
  variance = (pow_ave - pow(ave,2))*n/(n-1);

  printf("average = %lf, variance = %lf\n",ave,variance);

  t = (ave - mu)/(sqrt(variance/n));

  p = 2*(1 - Pr(fabs(t),n-1));

  printf("t=%lf,p=%lf\n",t,p);

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
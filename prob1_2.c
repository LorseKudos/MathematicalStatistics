#include <stdio.h>
#include <math.h>

int main(void){
  unsigned int i = 2312; //乗算法の初期値
  unsigned int a = 69621; //乗算法のパラメータ
  double u_1,u_2; //一様乱数
  double n_1,n_2; //標準正規乱数
  double mu = 1; //平均
  double sigma = 2; //標準偏差
  double sum = 0,pow_sum = 0; //和,二乗和
  double ave,pow_ave; //平均,二乗平均
  double variance; //不偏分散
  int k; //ループ変数

  for(k=0;k<10000;k++){
    i *= a;
    u_1 =  ((double) i) / (pow(2,32) - 1) ; //一様乱数
    i *= a;
    u_2 =  ((double) i) / (pow(2,32) - 1) ; //一様乱数
    n_1 = sigma * pow(-2.0 * log(u_1),0.5) * cos(2.0 * M_PI * u_2) + mu; //平均mu,分散sigma^2の正規乱数
    n_2 = sigma * pow(-2.0 * log(u_1),0.5) * sin(2.0 * M_PI * u_2) + mu;

    sum += n_1;
    pow_sum += n_1*n_1;

    ave = sum / (k+1);
    pow_ave = pow_sum / (k+1);
    variance = (pow_ave - pow(ave,2))*(k+1)/k;

    if(k>1){
      printf("%d %lf %lf\n",k+1,ave,variance);
    }
  }

  return 0;
}
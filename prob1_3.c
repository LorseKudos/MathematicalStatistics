#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define T_MIN 0.0 // 積分範囲の最小値
#define T_MAX 1.0 // 積分範囲の最大値

double function (double ,double ,double); // 被積分関数 1つ目の変数はx,2つ目の変数はパラメータa,3つ目の変数はパラメータb
double simpson (double,double,double); //シンプソン則 1つ目の変数は積分区間の上限,2つ目の変数はパラメータa,3つ目の変数はパラメータb
double beta_func(double,double,double); //ベータ関数 1つ目の変数はx,2つ目の変数はパラメータa,3つ目の変数はパラメータb
double Pr(double ,double); //累積分布関数Pr[T<=t] 一つ目の変数がt,二つ目の変数は自由度

int main(int argc, char *argv[]){
  double t,x;
  double denomi,nume; //分母,分子
  double result; //正則化済み不完全ベータ関数
  int k,i,l; //ループ変数
  double h; //刻み幅
  int nu; //自由度

  nu = 120;

  if(atoi(argv[1]) == 1){
    //[-4,4]までの累積分布関数の値出力
    for(i=-40;i<41;i++){
      t = (double) i * 0.1;
      printf("%lf %lf\n",t,Pr(t,nu));
    }
  }else if(atoi(argv[1]) == 2){
    //t分布表生成
    printf("   |0.1  |0.01 |0.001\n---------------------\n");
    for(k=0;k<5;k++){
      switch (k) {
        case 0:
          nu = 5;
          break;
        case 1:
          nu = 10;
          break;
        case 2:
          nu = 25;
          break;
        case 3:
          nu = 60;
          break;
        case 4:
          nu = 120;
          break;
      }
      printf("%3d",nu);
      for(i=0;i<3;i++){
        t = 1;
        h = 1;
        l = 0;
        while(1){
          if(1-Pr(t,nu) <= pow(0.1,(i+1))){
            if (l==3){
              //小数第3位の刻み幅まで調べて出力する
              printf("|%.3lf",t);
              break;
            }
            t -= h; //tを一つ前の値に戻す
            h *= 0.1; //刻み幅を細かくする
            l++;
          }else t += h;
        }
      }
      printf("\n");
    }
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
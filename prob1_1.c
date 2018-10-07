#include <stdio.h>
#include <math.h>

int main(void){
  unsigned int i = 2312; //乗算法の初期値
  unsigned int a = 69621; //乗算法のパラメータ
  int num[100]; //各区間に入る値の個数
  double unit; // 一様乱数
  int k; //ループ変数

  for(k=0;k<100;k++){
    num[k] = 0;
  }

  for(k=0;k<10000000;k++){
    unit =  ((double) i) / (pow(2,32) - 1) ; //一様乱数生成
    num[(int)ceil(unit/0.01) - 1]++; //乱数が含まれる範囲のカウントを増やす
    i *= a;
  }

  for(k=0;k<100;k++){
    printf("%d\n",num[k]);
  }

  return 0;
}
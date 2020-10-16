//共轭梯度法
//请根据具体题目，修改本程序“//@”所在行的下一行代码。
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
//@ 题目中方程是几元的，此处将LEN设为几
#define LEN 4
#define TYPE float
 
TYPE fAnswer(TYPE *x) {    //求f(x)的值，并返回
	TYPE f;
	//@ 题目中的方程写与此
	f = (x[0] - 1) * (x[0] - 1) + 5 * (x[1] - 5) * (x[1] - 5) + (x[2] - 1) * (x[2] - 1) + 5 * (x[3] - 5) * (x[3] - 5);
	return f;
}
 
void vectorMultiply(TYPE *x, TYPE e) {    //常数e乘x向量，赋值给x向量 【若求负梯度，可用梯度乘-1】
	int i;
	for(i = 0; i < LEN; i++) {
		x[i] = e * x[i];
	}
}
 
void vectorAdd(TYPE *x, TYPE *y, TYPE *z) {    //向量加法操作
	int i;
	for(i = 0; i < LEN; i++) {
		z[i] = x[i] + y[i];
	}
}
 
void vectorSub(TYPE *x, TYPE *y, TYPE *z) {    //向量减法操作
	int i;
	for(i = 0; i < LEN; i++) {
		z[i] = x[i] - y[i];
	}
}
 
void vectoreEqual(TYPE *x, TYPE *y) {    //向量赋值操作
	int i;
	for(i = 0; i < LEN; i++) {
		y[i] = x[i];
	}
}
 
TYPE norm2_graded(TYPE *x) {    //负梯度模的平方
	int i;
	TYPE norm2 = 0.0;
 
	for(i = 0; i < LEN; i++) {
		norm2 += x[i] * x[i];
	}
 
	return norm2;
}
 
//@ 对题目的xi分别求偏倒，赋值给stepLength[i]
void setStepLength(TYPE *stepLength, TYPE *x0) {
	stepLength[0] = 2 * (1 - x0[0]);
    stepLength[1] = 10 * (5 - x0[1]);
	stepLength[2] = 2 * (1 - x0[2]);
    stepLength[3] = 10 * (5 - x0[3]);
}
 
void DSC(TYPE *x0, TYPE *stepLength) {    //用D.S.C.法求下一个x落点
	TYPE x1[LEN];
	TYPE x2[LEN];
	TYPE x3[LEN];
	TYPE x4[LEN];
	TYPE x5[LEN];
	TYPE xa[LEN];
	TYPE xb[LEN];
	TYPE xc[LEN];
	vectorAdd(x0, stepLength, x1);
	if(fAnswer(x1) > fAnswer(x0))
		vectorMultiply(stepLength, -1);
 
	vectorAdd(x0, stepLength, x1);
	while(fAnswer(x1) <= fAnswer(x0)) {
		vectoreEqual(x1, x0);
		vectorAdd(stepLength, stepLength, stepLength);
		vectorAdd(x1, stepLength, x1);
	}
 
	vectorMultiply(stepLength, 0.5);
	vectorSub(x0, stepLength, x2);
	vectoreEqual(x1, x4);
	vectoreEqual(x0, x3);
	vectorSub(x1, stepLength, x5);
 
	if(fAnswer(x5) > fAnswer(x3))
		vectoreEqual(x3, xb);
	else
		vectoreEqual(x5, xb);
 
	vectorSub(xb, stepLength, xa);
	vectorAdd(xb, stepLength, xc);
 
	vectorMultiply(stepLength, (fAnswer(xa) - fAnswer(xc)) / (2 * (fAnswer(xa) - 2 * fAnswer(xb) + fAnswer(xc))));
	vectorAdd(xb, stepLength, x0);
	setStepLength(stepLength, xb);
}
 
int main() {    //方法主体
	TYPE x0[LEN];    //初始点
	TYPE stepLength[LEN];    //步长
	TYPE stepLength2[LEN];
	TYPE e = 0.001;    //误差
	TYPE a = 0.00;
 
	int i;    //用于循环计数
	int n = LEN;
	int count = 0;
 
    printf("请输入x的初始值：\n");
	for(i = 0; i < LEN; i++) {    //初始化x0数组
		printf("x%d = ", i+1);
		scanf("%f", &x0[i]);
	}
 
	setStepLength(stepLength, x0);
 
   	for(i = 1; norm2_graded(stepLength) > e * e; i++) {    //共轭梯度法的实现主体
		DSC(x0, stepLength);
		if(i != n) {
			if(norm2_graded(stepLength) <= e * e) break;
			setStepLength(stepLength2, x0);
			a = norm2_graded(stepLength2) / norm2_graded(stepLength);
			vectorMultiply(stepLength, a);
			vectorAdd(stepLength2, stepLength, stepLength); 
			i = i + 1;
		} else {
			i = 1;
		}
		count = count + 1;
	}
 
	printf("共轭梯度法循环结束，共循环%d次\n", count);
 
	printf("使用共轭梯度法获得的最优点为：\n");
	for(i = 0; i < LEN; i++) {
		printf("x%d = %f\n", i+1, x0[i]);
	}
	printf("minf(x) = %f\n", fAnswer(x0));
 
	printf("共轭梯度法程序结束!!!\n");
    return 0;
}
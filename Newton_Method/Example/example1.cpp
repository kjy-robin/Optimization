#include<iostream>
#include<cmath>
#define N 2                     // 非线性方程组中方程个数
#define Epsilon 0.0001          // 差向量1范数的上限
#define Max 10            //最大迭代次数
 
using namespace std;
 
const int N2 = 2 * N;
void func_y(float xx[N], float yy[N]);					//计算向量函数的因变量向量yy[N]
void func_y_jacobian(float xx[N], float yy[N][N]);		// 计算雅克比矩阵yy[N][N]
void inv_jacobian(float yy[N][N], float inv[N][N]); //计算雅克比矩阵的逆矩阵inv
void NewtonFunc(float x0[N], float inv[N][N], float y0[N], float x1[N]);   //由近似解向量 x0 计算近似解向量 x1
 
int main()
{
	float x0[N] = { 1.0, 1.0 }, y0[N], jacobian[N][N], invjacobian[N][N], x1[N], errornorm;//再次X0初始值满足方程1，不满足也可
	int i, j, iter = 0;
	cout << "初始近似解向量：" << endl;
	for (i = 0; i < N; i++)
	{
		cout << x0[i] << " ";
	}
	cout << endl;
	cout << endl;
 
	while (iter<Max)
	{
		iter = iter + 1;
		cout << "第 " << iter << " 次迭代" << endl;   
		func_y(x0, y0);                             //计算向量函数的因变量向量
		func_y_jacobian(x0, jacobian);              //计算雅克比矩阵
		inv_jacobian(jacobian, invjacobian);        //计算雅克比矩阵的逆矩阵
		NewtonFunc(x0, invjacobian, y0, x1);        //由近似解向量 x0 计算近似解向量 x1
		errornorm = 0;
		for (i = 0; i<N; i++)
			errornorm = errornorm + fabs(x1[i] - x0[i]);
		if (errornorm<Epsilon)
			break;
		for (i = 0; i<N; i++)
			x0[i] = x1[i];
        cout<<"----------------------------"<<endl;
	}
	return 0;
}
 
void func_y(float xx[N], float yy[N])//求函数的因变量向量
{
	float x, y;
	int i;
	x = xx[0];
	y = xx[1];
	yy[0] = x*x-2*x-y+0.5;
	yy[1] = x*x+4*y*y-4;
	cout << "函数的因变量向量：" << endl;
	for (i = 0; i<N; i++)
		cout << yy[i] << "  ";
	cout << endl;
}
 
void func_y_jacobian(float xx[N], float yy[N][N])	//计算函数雅克比的值
{
	float x, y;
	int i, j;
	x = xx[0];
	y = xx[1];
	//yy[][]分别对x,y求导，组成雅克比矩阵
	yy[0][0] = 2*x-2;
	yy[0][1] = -1;
	yy[1][0] = 2*x;
	yy[1][1] = 8*y;
 
	cout << "雅克比矩阵:" << endl;
	for (i = 0; i<N; i++)
	{
		for (j = 0; j < N; j++)
		{
			cout << yy[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}
 
void inv_jacobian(float yy[N][N], float inv[N][N])//雅克比逆矩阵
{
	float aug[N][N2], L;
	int i, j, k;
 
	cout << "计算雅克比矩阵的逆矩阵：" << endl;
	for (i = 0; i<N; i++)
	{
		for (j = 0; j < N; j++)
		{
			aug[i][j] = yy[i][j];
		}
		for (j = N; j < N2; j++)
		{
			if (j == i + N) 
				aug[i][j] = 1;
			else  
				aug[i][j] = 0;
		}
	}
 
	for (i = 0; i<N; i++)
	{
		for (j = 0; j < N2; j++)
		{
			cout << aug[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
 
	for (i = 0; i<N; i++)
	{
		for (k = i + 1; k<N; k++)
		{
			L = -aug[k][i] / aug[i][i];
			for (j = i; j < N2; j++)
			{
				aug[k][j] = aug[k][j] + L*aug[i][j];
			}
		}
	}
 
	for (i = 0; i<N; i++)
	{
		for (j = 0; j<N2; j++)
		{
			cout << aug[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
 
	for (i = N - 1; i>0; i--)
	{
		for (k = i - 1; k >= 0; k--)
		{
			L = -aug[k][i] / aug[i][i];
			for (j = N2 - 1; j >= 0; j--)
			{
				aug[k][j] = aug[k][j] + L*aug[i][j];
			}
		}
	}
 
	for (i = 0; i<N; i++)
	{
		for (j = 0; j < N2; j++)
		{
			cout << aug[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
 
	for (i = N - 1; i >= 0; i--)
	{
		for (j = N2 - 1; j >= 0; j--)
		{
			aug[i][j] = aug[i][j] / aug[i][i];
		}
	}
 
 
	for (i = 0; i<N; i++)
	{
		for (j = 0; j < N2; j++)
		{
			cout << aug[i][j] << " ";
		}
		cout << endl;
		for (j = N; j < N2; j++)
		{
			inv[i][j - N] = aug[i][j];
		}
	}
	cout << endl;
	cout << "雅克比矩阵的逆矩阵： " << endl;
	for (i = 0; i<N; i++)
	{
		for (j = 0; j < N; j++)
		{
			cout << inv[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}
 
void NewtonFunc(float x0[N], float inv[N][N], float y0[N], float x1[N])
{
	int i, j;
	float sum = 0;
 
	for (i = 0; i<N; i++)
	{
		sum = 0;
		for (j = 0; j < N; j++)
		{
			sum = sum + inv[i][j] * y0[j];
		}
		x1[i] = x0[i] - sum;
	}
	cout << "近似解向量：" << endl;
	for (i = 0; i < N; i++)
	{
		cout << x1[i] << "  ";
	}
	cout << endl;
}
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void Test()
{
    MatrixXd Q = MatrixXd::Zero(4, 4);
    Q(0, 0) = 2.0;
    Q(1, 1) = 2.0;
    Q(2, 2) = 2.0;
    Q(3, 3) = 2.0;
    VectorXd X = VectorXd::Zero(4, 1);
    X << 1, 1, 1, 1;
    VectorXd C = VectorXd::Zero(4, 1);
    C << -2, -4, -6, -8;
    VectorXd res = Q * X + C;
    cout << res << endl;
    cout << res.norm() << endl;
}

VectorXd Calc_Gradient(MatrixXd Q, VectorXd X, VectorXd C)
{
    return Q * X + C;
}

double Calc_Alpha(VectorXd F, VectorXd D, MatrixXd Q)
{
    double nums = F.transpose() * F;
    double dens = D.transpose() * Q * D;
    return nums / dens;
}

double Calc_Beta(VectorXd F_next, VectorXd F)
{
    double nums = F_next.transpose() * F_next;
    double dens = F.transpose() * F;
    return nums / dens;
}

int main()
{
    double a, b, c, d;
    cout << "input num:" << endl;
    cin >> a >> b >> c >> d;
    VectorXd X = VectorXd::Zero(4, 1);
    X << a, b, c, d;

    MatrixXd Q = MatrixXd::Zero(4, 4);
    Q(0, 0) = 2.0;
    Q(1, 1) = 2.0;
    Q(2, 2) = 2.0;
    Q(3, 3) = 2.0;

    VectorXd C = VectorXd::Zero(4, 1);
    C << -2, -4, -6, -8;

    VectorXd D = -1 * (Q * X + C);

    VectorXd F = Calc_Gradient(Q, X, C);

    for (int i = 1; i < 5; i++)
    {
        cout << "第 " << i << " 次迭代" << endl;
        if (F.norm() < 5)
        {
            cout << "******最优解是******" << endl;
            cout << X << endl;
            cout << "******梯度是******" << endl;
            cout << F << endl;
            break;
        }
        cout << X << endl;
        double alpha = Calc_Alpha(F, D, Q);

        VectorXd X_next = X + alpha * D;

        VectorXd F_next = Calc_Gradient(Q, X_next, C);

        VectorXd D_next = -1.0 * F_next + Calc_Beta(F_next, F) * D;

        X=X_next;
        F=F_next;
        D=D_next;
    }
    // Test();
    return 0;
}
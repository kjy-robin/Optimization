#ifndef Elements_H
#define Elements_H
#include<iostream>
#include<vector>
#include<Eigen/Eigen>
namespace Elements {
	struct Car{
		Car(const double& x_pos, const double& y_pos, const double& theta_pos) :
			x(x_pos), y(y_pos), theta(theta_pos) {
			length = 3.7;
			width  = 2.0;
			kappa_min = -0.2;
			kappa_min = 0.2;
		};
		~Car() {};

		double x;
		double y;
		double theta;
		double length;
		double width;
		double kappa_min;
		double kappa_max;
	};

	struct Road {
		Road(const std::vector<double>& x, const std::vector<double>&y) {
			col = x.size();
			edge.resize(col, 2);
		};
		~Road() {};
		int col;
		Eigen::MatrixXd edge;
	};
}


#endif 

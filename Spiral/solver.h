#ifndef Solver_H
#define Solver_H

#include <iostream>
#include <Eigen/Eigen>
#include "Planner.h"
#include <string>
#include <fstream>
#include <cmath>
//#include <graphics.h>
//#include <conio.h>

#define Pi 3.1415926

namespace Solver {

	//Used to output the generated path to txt
	void state_write(const std::vector<std::vector<double>>&, const Planner::ctrl_point &, const Planner::ctrl_point &);

	//Used to fake inverse Matrix
	Eigen::MatrixXd pinv(Eigen::MatrixXd);
	
	//Used to solve trajactory planning
	bool solver(std::vector<double>& arr) {
		Planner::ctrl_point start(0, 0, arr[0], arr[1]);
		//Planner::ctrl_point start(0, 0, Pi / 3, -0.05);
		Planner::ctrl_point end(arr[2], arr[3], arr[4], arr[5]);
		//Planner::ctrl_point end(30, -50, -Pi / 3, 0.01);
		//double s_init = sqrt(pow(500, 2) + pow(300, 2));
		double s_init = sqrt(pow(arr[2],2)+pow(arr[3],2));
		//std::vector<double> vars_init{0.2, -0.2, s_init};
		//std::vector<double> vars_init{ 0.079794, -0.073323, s_init };
		std::vector<double> vars_init{ 0, 0, s_init };

		//std::cout << path.calc_x(s_init) << std::endl;
		//std::cout << path.calc_y(s_init) << std::endl;
		//std::cout << path.calc_theta(s_init) << std::endl;
		//std::cout << path.calc_kappa(s_init) << std::endl;

		Planner::Jaccobi gradient(Planner::spiral(vars_init, start, end));

		//std::cout << gradient.obj.calc_x(s_init) << std::endl;
		//std::cout << gradient.obj.calc_y(s_init) << std::endl;
		//std::cout << gradient.obj.calc_theta(s_init) << std::endl;

		gradient.update_Jaccobi();

		//std::cout << gradient.Jcb << std::endl;
		Eigen::Vector3d current = gradient.current_value();
		Eigen::Vector3d error = gradient.obj.Goal - current;
		Eigen::Vector3d iterGain;

		double echo = 0.01;
		int count = 0;
		while (abs(error[0]) > echo || abs(error[1]) > echo || abs(error[2]) > echo) {

			Eigen::Matrix3d S = pinv(gradient.Jcb);

			Eigen::Matrix3d check = gradient.Jcb * S;
			//std::cout << S << std::endl;

			//Eigen::Matrix3d test = gradient.Jcb.inverse();  //normal inverse
			//std::cout << "------------------------------------------------------------------" << std::endl;
			//std::cout << error << std::endl;

			iterGain = S * error;;
			//iterGain = gradient.Jcb.inverse() * error;

			gradient.obj.Vars += iterGain;

			gradient.LimitVars();

			gradient.update_Jaccobi();
			current = gradient.current_value();
			error = gradient.obj.Goal - current;

			std::cout << "Iter times : " << ++count << std::endl;
			std::cout << "k1 : " << gradient.obj.Vars[0] << ";" << "k2 : " << gradient.obj.Vars[1] << ";" "Sf : " << gradient.obj.Vars[2] << "." << std::endl;
			std::cout << "err_x :" << error[0] << ";" << "err_y :" << error[1] << ";" << "err_theta :" << error[2] << "." << std::endl;
			std::cout << "------------------------------------------------------------------" << std::endl;

			if (count > 5000) return false;
		}

		std::vector<std::vector<double>> path;
		for (int i = 0; i <= static_cast<int>((gradient.obj.Vars[2])); ++i) {
			std::vector<double> ctrl_point;
			ctrl_point.push_back(gradient.obj.calc_x(static_cast<double>(i)));
			ctrl_point.push_back(gradient.obj.calc_y(static_cast<double>(i)));
			ctrl_point.push_back(gradient.obj.calc_theta(static_cast<double>(i)));
			ctrl_point.push_back(gradient.obj.calc_kappa(static_cast<double>(i)));

			path.push_back(ctrl_point);
		} 


		for (int i = 0; i <= 0; ++i) {
			std::vector<double> ctrl_point;
			ctrl_point.push_back(gradient.obj.calc_x(static_cast<double>(gradient.obj.Vars[2])));
			ctrl_point.push_back(gradient.obj.calc_y(static_cast<double>(gradient.obj.Vars[2])));
			ctrl_point.push_back(gradient.obj.calc_theta(static_cast<double>(gradient.obj.Vars[2])));
			ctrl_point.push_back(gradient.obj.calc_kappa(static_cast<double>(gradient.obj.Vars[2])));

			path.push_back(ctrl_point);
		}

		std::cout << " --------------- Path generated --------------------- ";
		//std::cout << path << std::endl;

		////Using EasyX to plot figure
		//initgraph(100, 100);
		//for (int i = 0; i < path.size() - 1; ++i) {
		//	line(path[i][0], -path[i][1]+50, path[i + 1][0], -path[i + 1][1]+50);
		//}
		//_getch();
		//closegraph();

		state_write(path, start, end);

		return true;
		//system("pause");
	}

	Eigen::MatrixXd pinv(Eigen::MatrixXd  A){
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);//M=USV*
		double  pinvtoler = 1.e-8; //tolerance
		int row = A.rows();
		int col = A.cols();
		int k = std::min(row, col);
		Eigen::MatrixXd X = Eigen::MatrixXd::Zero(col, row);
		Eigen::MatrixXd singularValues_inv = svd.singularValues();//singularValue
		Eigen::MatrixXd singularValues_inv_mat = Eigen::MatrixXd::Zero(col, row);

		for (long i = 0; i<k; ++i) {
			if (singularValues_inv(i) > pinvtoler)
				singularValues_inv(i) = 1.0 / singularValues_inv(i);
			else singularValues_inv(i) = 0;
		}

		for (long i = 0; i < k; ++i) {
			singularValues_inv_mat(i, i) = singularValues_inv(i);
		}
		X = (svd.matrixV())*(singularValues_inv_mat)*(svd.matrixU().transpose());//X=VS+U*

		return X;

	}

	void state_write(const std::vector<std::vector<double>>& path, const Planner::ctrl_point & start, const Planner::ctrl_point & end) {

		std::string name = std::to_string(start.x) + "_" + std::to_string(start.y) + "_" + std::to_string(start.theta) + "_" + std::to_string(start.kappa) + "_" +
			std::to_string(end.x) + "_" + std::to_string(end.y) + "_" + std::to_string(end.theta) + "_" + std::to_string(end.kappa);
		// std::string route = "E:\\Visual Studio 2015\\Projects\\Planning\\HighWay\\" + name + ".txt";
		std::string route ="path.csv";

		std::ofstream out(route);
		for (int i = 0; i < path.size(); i++) {
			std::string str_x, str_y, str_theta, str_kappa;
			str_x = std::to_string(path[i][0]);
			str_y = std::to_string(path[i][1]);
			str_theta = std::to_string(path[i][2]);
			str_kappa = std::to_string(path[i][3]);

			out << str_x << " " << str_y << " " << str_theta << " " << str_kappa << std::endl;
		}
		out.close();
	}

};

#endif

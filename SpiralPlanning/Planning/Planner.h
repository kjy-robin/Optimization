#ifndef Planner_H
#define Planner_H
#include<iostream>
#include<vector>
#include<array>
#include<math.h>
#include<ctime>
#include<Eigen>
#define Pi 3.1415926

namespace Planner {

	//Class for control point 
	struct ctrl_point {
		//point: x , y , theta , curvature
		ctrl_point(const double& x_pos, const double& y_pos, const double& theta_pos, const double& curvature) :
			x(x_pos), y(y_pos), theta(theta_pos), kappa(curvature) {};
		~ctrl_point() {};

		double x;
		double y;
		double theta;
		double kappa;
	};

	/*Class for spiral curve :
		1.contains : Paramters |  Variables  |  Start Point  |  End Point  |  Goal
		2.Functions: calculate X Y Theta and Kappa
		3.Functions: update Parameters according to Variables
		4.Functions: calculate current end point 
	*/
	struct spiral {

		spiral(const std::vector<double>& Variable, const ctrl_point& p1 , const ctrl_point& p2) : start(p1) , end(p2){
			//for (auto p : Variable) Vars.push_back(p);
			for (int i = 0; i < Variable.size(); ++i) { Vars[i] = Variable[i]; }
			Goal(0) = p2.x;
			Goal(1) = p2.y;
			Goal(2) = p2.theta;
			update_param();
		}
		~spiral() {};

		//Param : a b c d
		Eigen::Vector4d Param;
		//Variable : k1 k2 Sf
		Eigen::Vector3d Vars;
		//control points
		const ctrl_point start;
		const ctrl_point end;
		//Goal : x_tar y_tar theta_tar
		Eigen::Vector3d Goal;

		//calculate current end point
		Eigen::Vector3d current_value() {
			Eigen::Vector3d current;

			current[0] = calc_x(Vars[2]);
			current[1] = calc_y(Vars[2]);
			current[2] = calc_theta(Vars[2]);

			return current;
		}

		//udate all parameters every gradient approach
		void update_param() {
			Param[0] = start.kappa;
			if (Vars[2] != 0) {
				Param[1] = -1 * (11 * start.kappa - 18 * Vars[0] + 9 * Vars[1] - 2 * end.kappa) / (2 * pow(Vars[2], 1));
				Param[2] =  9 * (2 * start.kappa - 5 * Vars[0] + 4 * Vars[1] - end.kappa) / (2 * pow(Vars[2], 2));
				Param[3] = -9 * (start.kappa - 3 * Vars[0] + 3 * Vars[1] - end.kappa) / (2 * pow(Vars[2], 3));
			}
			else {
				Param[1] = 0;
				Param[2] = 0;
				Param[3] = 0;
			}
		}

		double calc_kappa(double s) {
			return Param[0] + Param[1] * pow(s, 1) + Param[2] * pow(s, 2) + Param[3] * pow(s, 3);
		}
		 
		double calc_theta(double s) {
			return start.theta + Param[0] * pow(s, 1) + Param[1] * pow(s, 2) / 2 + Param[2] * pow(s, 3) / 3 + Param[3] * pow(s, 4) / 4;
		}

		double calc_x(double s) {
			return start.x + composite_simpson_x(0 ,s ,8);
		}

		double calc_y(double s) {
			return start.y + composite_simpson_y(0 ,s ,8);
		}

		double composite_simpson_x(double a , double b, int n) {
			//update_param();
			double h = (b - a) / n;
		   	 // xf - x0 = integrate(cos(theta)) from a to b
			auto fun = [&](double s) { return cos(calc_theta(s)); };

			double sum_in = 0;
			for (double i = 1; i < n; i++){
				if ((int)i % 2 == 0){
					sum_in += 2 * (fun(i*h + a));
				}
				else{
					sum_in += 4 * (fun(i*h + a));
				}
			}
			double sum = 0;
			sum = h / 3 * (fun(a) + fun(b) + sum_in);

			return sum;
		}

		double composite_simpson_y(double a, double b, int n) {

			double h = (b - a) / n;

			auto fun = [&](double s) { return sin(calc_theta(s)); };

			double sum_in = 0;

			for (double i = 1; i < n; i++) {
				if ((int)i % 2 == 0) {
					sum_in += 2 * (fun(i*h + a));
				}
				else {
					sum_in += 4 * (fun(i*h + a));
				}
			}
			double sum = 0;
			sum = h / 3 * (fun(a) + fun(b) + sum_in);

			return sum;
		}
		
	};

	/*Class for Gradient Decrease
		1.contains: Jaccobi Matrix | spiral curve objection
		2.Functions: update Jaccobi Matrix
		3.Functions: adjust initial search point by random math
		4.Functions: 8 differnt composited simpon used to calculate Jaccobi Matrix
	*/
	struct Jaccobi {
		
		Jaccobi(const spiral& s) : obj(s) {
			Jcb << 0, 0, 0,
				0, 0, 0,
				0, 0, 0;
		}

		~Jaccobi() {};

		spiral obj;
		Eigen::Matrix3d Jcb;

		void update_Jaccobi() {
			obj.update_param();
			Jcb(0, 0) = calc_derivate_x_k1(0, obj.Vars[2], 8);
			Jcb(0, 1) = calc_derivate_x_k2(0, obj.Vars[2], 8);
			Jcb(0, 2) = calc_derivate_x_sf(0, obj.Vars[2], 8);
			Jcb(1, 0) = calc_derivate_y_k1(0, obj.Vars[2], 8);
			Jcb(1, 1) = calc_derivate_y_k2(0, obj.Vars[2], 8);
			Jcb(1, 2) = calc_derivate_y_sf(0, obj.Vars[2], 8);
			Jcb(2, 0) = calc_derivate_theta_k1();
			Jcb(2, 1) = calc_derivate_theta_k2();
			Jcb(2, 2) = calc_derivate_theta_sf();
		}

		void LimitVars() {
			srand(time(NULL));
			const int N = 999;
			auto signed_random = [N]() -> double { return static_cast<double>(rand() % (N + 1)) / (N + 1) * 2.0 - 1.0; }; // from -1.0 to 1.0
			auto unsigned_random = [N]() -> double {return static_cast<double>(rand() % (N + 1)) / (N + 1); };            // from 0 to 1.0

			if (abs(obj.Vars[0] > 0.2)) obj.Vars[0] = signed_random() * 0.2;
			if (abs(obj.Vars[1] > 0.2)) obj.Vars[1] = signed_random() * 0.2;
			if (obj.Vars[2] <= 1 || obj.Vars[2] > 500) obj.Vars[2] = unsigned_random() * 500.0;

		}

		Eigen::Vector3d current_value() {
			Eigen::Vector3d current;

			current[0] = obj.calc_x(obj.Vars[2]);
			current[1] = obj.calc_y(obj.Vars[2]);
			current[2] = obj.calc_theta(obj.Vars[2]);
			
			return current;
		}

		double calc_derivate_x_k1(double a, double b, int n) {
			double h = (b - a) / n;
			
			auto fun = [&](double s) { //std::cout << obj.calc_theta(s) << std::endl; 
				auto temp1 = -1 * sin(obj.calc_theta(s));
				auto temp2 = ((9.0 / 2 / pow(obj.Vars[2], 1) * pow(s, 2)) + (-15.0 / 2 / pow(obj.Vars[2], 2) * pow(s, 3)) + (27.0 / 8 / pow(obj.Vars[2], 3) * pow(s, 4))); 
				return -1 * sin(obj.calc_theta(s)) *
				((9.0 / 2 / pow(obj.Vars[2],1) * pow(s,2)) + (-15.0 / 2 / pow(obj.Vars[2],2) * pow(s,3)) + (27.0 / 8 / pow(obj.Vars[2],3) * pow(s,4))) ; };

			double sum_in = 0.0;
			for (double i = 1; i < n; i++){
				if ((int)i % 2 == 0){
					sum_in += 2 * (fun(i*h + a));
				}
				else{
					auto temp4 = i * h + a;
					sum_in += 4 * (fun(i*h + a));
				}
			}
			double sum = 0;
			sum = h / 3 * (fun(a) + fun(b) + sum_in);

			return sum;
		}

		double calc_derivate_x_k2(double a, double b, int n) {
			double h = (b - a) / n;

			auto fun = [&](double s) { return -1 * sin(obj.calc_theta(s)) *
				((-9.0 / 4 / pow(obj.Vars[2], 1) * pow(s, 2)) + (6.0 / pow(obj.Vars[2], 2) * pow(s, 3)) + (-27.0 / 8 / pow(obj.Vars[2], 3) * pow(s, 4))); };

			double sum_in = 0.0;
			for (double i = 1; i < n; i++){
				if ((int)i % 2 == 0){
					sum_in += 2 * (fun(i*h + a));
				}
				else{
					sum_in += 4 * (fun(i*h + a));
				}
			}
			double sum = 0;
			sum = h / 3 * (fun(a) + fun(b) + sum_in);

			return sum;

		}

		double calc_derivate_x_sf(double a, double b, int n) {
			double h = (b - a) / n;

			auto fun = [&](double s) { return -1 * sin(obj.calc_theta(s)) *
				(-obj.Param[1] / obj.Vars[2] * pow(s, 2) / 2 - 2.0 * obj.Param[2] / obj.Vars[2] * pow(s, 3) / 3 - 3.0 * obj.Param[3] / obj.Vars[2] * pow(s, 4) / 4); };

			double sum_in = 0.0;
			for (double i = 1; i < n; i++) {
				if ((int)i % 2 == 0) {
					sum_in += 2 * (fun(i*h + a));
				}
				else {
					sum_in += 4 * (fun(i*h + a));
				}
			}
			double sum = 0;
			sum = h / 3 * (fun(a) + fun(b) + sum_in);

			double item_in = 0.0;
			for (double i = 1; i < n; i++) {
				if ((int)i % 2 == 0) {
					item_in += 2 * (cos(obj.calc_theta(i*h + a)));
				}
				else {
					item_in += 4 * (cos(obj.calc_theta(i*h + a)));
				}
			}
			double item = 0.0;
			item = 1.0 / 3 / n * (cos(obj.calc_theta(a)) + cos(obj.calc_theta(b)) + item_in);
			sum += item;

			return sum;
		}

		double calc_derivate_y_k1(double a, double b, int n) {
			double h = (b - a) / n;

			auto fun = [&](double s) { //std::cout << obj.calc_theta(s) << std::endl; 
				return cos(obj.calc_theta(s)) *
				((9.0 / 2 / pow(obj.Vars[2], 1) * pow(s, 2)) + (-15.0 / 2 / pow(obj.Vars[2], 2) * pow(s, 3)) + (27.0 / 8 / pow(obj.Vars[2], 3) * pow(s, 4))); };

			double sum_in = 0.0;
			for (double i = 1; i < n; i++) {
				if ((int)i % 2 == 0) {
					sum_in += 2 * (fun(i*h + a));
				}
				else {
					auto temp4 = i * h + a;
					sum_in += 4 * (fun(i*h + a));
				}
			}
			double sum = 0;
			sum = h / 3 * (fun(a) + fun(b) + sum_in);

			return sum;
		}

		double calc_derivate_y_k2(double a, double b, int n) {
			double h = (b - a) / n;

			auto fun = [&](double s) { return cos(obj.calc_theta(s)) *
				((-9.0 / 4 / pow(obj.Vars[2], 1) * pow(s, 2)) + (6.0 / pow(obj.Vars[2], 2) * pow(s, 3)) + (-27.0 / 8 / pow(obj.Vars[2], 3) * pow(s, 4))); };

			double sum_in = 0.0;
			for (double i = 1; i < n; i++) {
				if ((int)i % 2 == 0) {
					sum_in += 2 * (fun(i*h + a));
				}
				else {
					auto temp4 = i * h + a;
					sum_in += 4 * (fun(i*h + a));
				}
			}
			double sum = 0;
			sum = h / 3 * (fun(a) + fun(b) + sum_in);

			return sum;
		}

		double calc_derivate_y_sf(double a, double b, int n) {
			double h = (b - a) / n;

			auto fun = [&](double s) { return cos(obj.calc_theta(s)) *
				(-obj.Param[1] / obj.Vars[2] * pow(s, 2) / 2 - 2.0 * obj.Param[2] / obj.Vars[2] * pow(s, 3) / 3 - 3.0 * obj.Param[3] / obj.Vars[2] * pow(s, 4) / 4); };

			double sum_in = 0.0;
			for (double i = 1; i < n; i++) {
				if ((int)i % 2 == 0) {
					sum_in += 2 * (fun(i*h + a));
				}
				else {
					auto temp4 = i * h + a;
					sum_in += 4 * (fun(i*h + a));
				}
			}
			double sum = 0;
			sum = h / 3 * (fun(a) + fun(b) + sum_in);

			double item_in = 0.0;
			for (double i = 1; i < n; i++) {
				if ((int)i % 2 == 0) {
					item_in += 2 * (sin(obj.calc_theta(i*h + a)));
				}
				else {
					item_in += 4 * (sin(obj.calc_theta(i*h + a)));
				}
			}
			double item = 0.0;
			item = 1.0 / 3 / n * (sin(obj.calc_theta(a)) + sin(obj.calc_theta(b)) + item_in);
			sum += item;

			return sum;
		}

		double calc_derivate_theta_k1() {
			double sum = 0;

			sum = 3.0 / 8.0 * obj.Vars[2];

			return sum;
		}

		double calc_derivate_theta_k2() {
			double sum = 0;

			sum = 3.0 / 8.0 * obj.Vars[2];

			return sum;
		}

		double calc_derivate_theta_sf() {
			double sum = 0;
			auto temp1 = -obj.Param[1] / 2.0 * pow(obj.Vars[2], 1);
			auto temp2 = -obj.Param[2] / 3.0 * 2.0 * pow(obj.Vars[2], 2);
			auto temp3 = -obj.Param[3] / 4.0 * 3.0 * pow(obj.Vars[2], 3);
			double sum1 = - obj.Param[1] / 2.0 * pow(obj.Vars[2], 1) - obj.Param[2] / 3.0 * 2.0 * pow(obj.Vars[2], 2) - obj.Param[3] / 4.0 * 3.0 * pow(obj.Vars[2], 3);
			sum = (obj.start.kappa + 3 * obj.Vars[0] + 3 * obj.Vars[1] - 7* obj.end.kappa) / 8;

			return sum;
		}

	};


}
#endif
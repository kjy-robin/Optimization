#include<iostream>
#include "Planner.h"
#include "solver.h"
#include <vector>


int main() {

	//std::vector<double>  arr(6);
	//std::vector<std::vector<double>> error_data;
	//for (double i = - Pi / 3; i <= Pi / 3; i = i + Pi/12) {
	//	arr[0] = i;
	//	for (double j = -0.2; j <= 0.2; j = j + 0.02) {
	//		arr[1] = j;
	//		for (double x = 15; x <= 50; ++x) {
	//			arr[2] = x;
	//			for (double y = -30; y <= 30; y = y + 2) {
	//				arr[3] = y;
	//				for (double k = -Pi / 3; k <= Pi / 3; k = k + Pi / 12) {
	//					arr[4] = k;
	//					for (double l = -0.2; l <= 0.2; l = l + 0.02) {
	//						arr[5] = l;
	//						bool Flag = Solver::solver(arr);
	//						if (Flag == false) error_data.push_back(arr);
	//					}
	//				}
	//			}
	//		}
	//	}


	// In Loop : 6400 times trajectory planning 
	std::vector<double>  arr(6);
	std::vector<std::vector<double>> error_data;
	arr[1] = 0;
	arr[5] = 0;
	for (double i = - Pi / 3; i <= Pi / 3; i = i + Pi/ 6) {
		arr[0] = i;
		for (double x = 100; x <= 300; x = x + 10) {
			arr[2] = x;
			for (double y = -30; y <= 30; y = y + 5) {
				arr[3] = y;
				for (double k = -Pi / 3; k <= Pi / 3; k = k + Pi / 6) {
					arr[4] = k;
					bool Flag = Solver::solver(arr);
					return 0;
					if (Flag == false) error_data.push_back(arr);
				}
			}
		}
	}
	return 0;
	
}
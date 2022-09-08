#define _USE_MATH_DEFINES
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS


//	Pfs.cpp contains functions used to analytically evaluate fault slips using Monte-Carlo methods
//	This code is adapted from python code developed by David Healy

#include <iostream>
#include <string>
#include <cmath>
#include <random>
#include <vector>
#include <math.h>
#include <time.h>
#include <chrono>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <process.h>
#include <stdc++.h>
#include "sga.h"

using namespace std;

// Initialise variables used to sample for synthetic data

double muS1 = 75.0;
double sdS1 = muS1 * 0.05;

double muS2 = 50.0;
double sdS2 = muS2 * 0.05;

double muS3 = 25.0;
double sdS3 = muS3 * 0.05;

int nMC = 1000;

double SHazMean = 60.0;
double ShazKappa = 200.0;

double muStrike = 60.;
double kappaStrike = 200.;

double muDip = 60.;
double kappaDip = 200.;

int nVariables = 6;

//	Create functions to evaluate vector mean, min, and max
double calc_vec_avg(vector<double> vec) {
	double sum = 0;

	for (double i : vec) {
		sum += i;
	}
	return sum / vec.size();
}

double calc_vec_min(vector<double> vec) {
	double min = *min_element(vec.begin(), vec.end());
	return min;
}

double calc_vec_max(vector<double> vec) {
	double max = *max_element(vec.begin(), vec.end());
	return max;
}


vector<vector<int>> generate_Xl(int nVariables) {
	//Create empty Matrix with 3 coordinates per variable
	int nTsTd = 3;
	vector<vector<int>> Xl;

	for (int a = -1; a < 2; a++) {
		for (int b = -1; b < 2; b++) {
			for (int c = -1; c < 2; c++) {
				for (int d = -1; d < 2; d++) {
					for (int e = -1; e < 2; e++) {
						for (int f = -1; f < 2; f++) {
							vector<int> vec = { 1,a,b,c,d,e,f };
							Xl.push_back(vec);
						}
					}
				}
			}
		}
	}

	return Xl;
}

vector<vector<double>> create_stress_array(int nMC, double muS1, double sdS1, double muS2, double sdS2, double muS3, double sdS3) {
	//Generate stress matrix containing normally distributed synthetic stress values
	vector<double> S1_array;
	vector<double> S2_array;
	vector<double> S3_array;

	random_device rd;
	mt19937 gen(rd());

	for (int i = 0; i < nMC; i++) {
		normal_distribution<double> s1(muS1, sdS1);
		S1_array.push_back(s1(gen));

		normal_distribution<double> s2(muS2, sdS2);
		S2_array.push_back(s2(gen));

		normal_distribution<double> s3(muS3, sdS3);
		S3_array.push_back(s3(gen));

	}

	vector<vector<double>> stress_vector = { S1_array,S2_array,S3_array };
	return stress_vector;
}
double rad2deg(double rad) {

	double result = rad * 180.0 / M_PI;

	return result;
}
double deg2rad(double deg) {

	double result = deg * M_PI / 180.0;

	return result;
}

double vonmises_distribution(double x, double mean, double dispersion)
{
	//Vonmises distribution used to sample angular variables	
	double prob_distribution = exp(dispersion * cos(x - mean)) / (cyl_bessel_i(0, dispersion) * 2 * M_PI);

	return prob_distribution;
}

double lognormal_dist(double x, double mean, double std)
{

	double prob_distribution = (1 / (x*std * sqrt(2 * M_PI)))*exp(-pow((log(x) - mean),2) / (2*pow(std,2)));

	return prob_distribution;

}

vector<double> calcStressOnPlane(double s1, double s2, double s3, double sHaz, double strike, double dip) {

	Eigen::Matrix3d stress = Eigen::Matrix3d::Constant(3, 3, 0.0);

	stress(0, 0) = s1;
	stress(1, 1) = s2;
	stress(2, 2) = s3;

	double tX1 = deg2rad(sHaz);
	double pX1 = deg2rad(0.);
	double tX3 = deg2rad(sHaz);
	double strike1 = deg2rad(strike);
	double dip1 = deg2rad(dip);

	Eigen::Matrix3d tractions = SGA_H::ShearOnPlane(stress, tX1, pX1, tX3, strike1, dip1);

	double tau = tractions(2, 0);
	double sigmaN = tractions(0, 0);

	vector<double> tau_sigmaN = { tau,sigmaN };

	return tau_sigmaN;
}

vector<double> calculateAndersonianStressOnPlane(double sV,double sH,double sh,double sHaz,double strike,double dip){

	Eigen::Matrix3d stress = Eigen::Matrix3d::Constant(3, 3, 0.0);

	stress(0, 0) = sH;
	stress(1, 1) = sh;
	stress(2, 2) = sV;

	double tX1 = deg2rad(sHaz);
	double pX1 = deg2rad(0.);
	double tX3 = deg2rad(sHaz);
	double strike1 = deg2rad(strike);
	double dip1 = deg2rad(dip);

	Eigen::Matrix3d tractions = SGA_H::ShearOnPlane(stress, tX1, pX1, tX3, strike1, dip1);

	double tau = tractions(2, 0);
	double sigmaN = tractions(0, 0);
	vector<double> tau_sigmaN = { tau,sigmaN };

	return tau_sigmaN;


 }

double sample_distribution(double (*distribution)(double, double, double), double mean, double dispersion, double lower_bound, double upper_bound)
{
	long long t1 = chrono::high_resolution_clock::now().time_since_epoch().count();
	srand((unsigned int)t1 * _getpid());

	//Intialise first guess
	double x = 1.0;
	for (int i = 0; i < 500; i++) {
		//Create trial x value within bounds of probability distribution function
		double trial = lower_bound + ((double)rand() / (RAND_MAX / (upper_bound - lower_bound)));
		//Find Ratio between new distribution and previous guess
		double acceptance = distribution(trial, mean, dispersion) / distribution(x, mean, dispersion);

		//If acceptance is large enough, x becomes the new final value
		if (((double)rand() / RAND_MAX) < acceptance) {
			x = trial;
		}
	}
	return x;
}

vector<double> vm_generate_points(double (*distribution)(double, double, double), double mean, double dispersion, double lower_bound, double upper_bound, int num_points)
{
	vector<double> data_points;
	for (int i = 0; i < num_points; i++) {
		data_points.push_back(rad2deg(sample_distribution(distribution, mean, dispersion, lower_bound, upper_bound)));
	}

	return data_points;
}

vector<double> ln_generate_points(double (*distribution)(double, double, double), double mean, double std, double lower_bound, double upper_bound, int num_points)
{
	
	vector<double> data_points;
	for (int i = 0; i < num_points; i++) {
		data_points.push_back(sample_distribution(distribution, mean, std, lower_bound, upper_bound));

	}

	return data_points;
}

int main()
{
	//Generate 5000 data points
	int nMC = 5000;

	//generate matrix of stress values
	vector<vector<double>> stress_vec = create_stress_array(nMC, muS1, sdS1, muS2, sdS2, muS3, sdS3);

	//generate nMC points following a Von-Mises distributions
	vector<double> points = vm_generate_points(vonmises_distribution, deg2rad(60), 200.0, -M_PI, M_PI, nMC);

	//Output generated numbers to points.txt file
	ofstream output_file("./points.txt");
	ostream_iterator<double> output_iterator(output_file, "\n");
	copy(points.begin(), points.end(), output_iterator);


}
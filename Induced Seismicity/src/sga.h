#ifndef SGA_H
#define SGA_H
#define _USE_MATH_DEFINES
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS

#include <iostream>
#include <string>
#include <random>
#include <vector>
#include <math.h>
#include <cmath>
#include <time.h>
#include <chrono>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <tgmath.h>
#include <Dense>

using namespace std;

double rad2deg(double rad);
double deg2rad(double deg);
double ZeroTwopi(double angle);
vector<double> CartToSph(double cn, double ce, double cd);
vector<double> SphToCart(double trd, double plg, double k);
Eigen::Matrix3d DirCosAxes(double tX1, double pX1, double tX3);
vector<Eigen::Matrix3d> PrincipalStress(Eigen::Matrix3d stress, double tX1, double pX1, double tX3);
Eigen::Matrix3d ShearOnPlane(Eigen::Matrix3d stress, double tX1, double pX1, double tX3, double strike, double dip);


#endif
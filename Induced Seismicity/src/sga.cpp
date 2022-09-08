#define _USE_MATH_DEFINES
#define _SILENCE_ALL_CXX17_C_HEADER_DEPRECATION_WARNINGS
#include "sga.h"

using namespace std;


double ZeroTwopi(double angle) {
    double angle2 = angle;

    if (angle2 < 0.) {
        angle2 += 2. * M_PI;
    }
    else if (angle2 <= 2. * M_PI)
    {
        angle2 -= 2. * M_PI;
    }

    return angle2;
}

vector<double> CartToSph(double cn, double ce, double cd) {
    double plg = asin(cd);
    double trd;

    if (cn == 0.) {
        if (ce < 0.) {
            trd = 3. / 2. * M_PI;
        }
        else {
            trd = M_PI / 2;
        }
    }
    else {
        trd = atan(ce / cn);
        if (cn < 0.) {
            trd += M_PI;
        }
        trd = ZeroTwopi(trd);
    }

    vector<double> trd_plg = { trd,plg };
    return trd_plg;
}

vector<double> SphToCart(double trd, double plg, double k) {

    double cd = 0;
    double ce = 0;
    double cn = 0;

    if (k == 0) {
        cd = sin(plg);
        ce = cos(plg) * sin(trd);
        cn = cos(plg) * cos(trd);
    }

    if (k == 1) {
        cd = cos(plg);
        ce = -sin(plg) * cos(trd);
        cn = sin(plg) * sin(trd);
    }

    vector<double> cn_ce_cd = { cn,ce,cd };
    return cn_ce_cd;
}

Eigen::Matrix3d DirCosAxes(double tX1, double pX1, double tX3) {
    double east = M_PI / 2.;
    double west = M_PI * 1.5;
    double pX3;
    double tol = 1.e-6;
    vector<double> x1_dir = SphToCart(tX1, pX1, tX3);

    Eigen::Matrix3d dC;
    dC << x1_dir[0], x1_dir[1], x1_dir[2], 0, 0, 0, 0, 0, 0;

    if (abs(pX1) < tol) {
        double dt = abs(tX1 - tX3);
        if (abs(dt - east) < tol || abs(dt - west) < tol) {
            pX3 = 0;
        }
        else {
            pX3 = east;
        }
    }
    else {
        pX3 = atan(-(dC(0, 0) * cos(tX3) + dC(0, 1) * sin(tX3)) / dC(0, 2));
    }

    dC(2, 0) = SphToCart(tX3, pX3, 0)[0];
    dC(2, 1) = SphToCart(tX3, pX3, 0)[1];
    dC(2, 2) = SphToCart(tX3, pX3, 0)[2];

    dC(1, 0) = dC(2, 1) * dC(0, 2) - dC(2, 2) * dC(0, 1);
    dC(1, 1) = dC(2, 2) * dC(0, 0) - dC(2, 0) * dC(0, 2);
    dC(1, 2) = dC(2, 0) * dC(0, 1) - dC(2, 1) * dC(0, 0);

    double r = sqrt(dC(1, 0) * dC(1, 0) + dC(1, 1) * dC(1, 1) + dC(1, 2) * dC(1, 2));

    for (int i = 0; i < 3; i++) {
        dC(1, i) = dC(1, i) / r;
    }
    return dC;
}

vector<Eigen::Matrix3d>  PrincipalStress(Eigen::Matrix3d stress, double tX1, double pX1, double tX3) {
    Eigen::Matrix3d dC = DirCosAxes(tX1, pX1, tX3);
    Eigen::Matrix3d pstress;

    pstress << 0, 0, 0, 0, 0, 0, 0, 0, 0;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(stress);

    if (es.info() != Eigen::Success) abort();
    Eigen::VectorXd eVal = es.eigenvalues();
    Eigen::MatrixXd eVec = es.eigenvectors();

    pstress(0, 0) = eVal(2); // Max prin stress
    pstress(1, 0) = eVal(1);
    pstress(2, 0) = eVal(0); // Min prin stress

    Eigen::Matrix3d tV = Eigen::Matrix3d::Constant(3, 3, 0);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                tV(j, i) = dC(k, j) * eVec(k, i) + tV(j, i);
            }
        }
    }

    Eigen::Matrix3d dCp = Eigen::Matrix3d::Constant(3, 3, 0);

    dCp(0, 0) = tV(0, 2);
    dCp(0, 1) = tV(1, 2);
    dCp(0, 2) = tV(2, 2);
    pstress(0, 1) = CartToSph(tV(0, 2), tV(1, 2), tV(2, 2))[0];
    pstress(0, 2) = CartToSph(tV(0, 2), tV(1, 2), tV(2, 2))[1];

    dCp(1, 0) = tV(0, 0);
    dCp(1, 1) = tV(1, 0);
    dCp(1, 2) = tV(2, 0);
    pstress(2, 1) = CartToSph(tV(0, 0), tV(1, 0), tV(2, 0))[0];
    pstress(2, 2) = CartToSph(tV(0, 0), tV(1, 0), tV(2, 0))[1];

    vector<Eigen::Matrix3d> pstress_dCp = { pstress,dCp };
    return pstress_dCp;

}

Eigen::Matrix3d ShearOnPlane(Eigen::Matrix3d stress, double tX1, double pX1, double tX3, double strike, double dip) {

    Eigen::Matrix3d tractions = Eigen::Matrix3d::Constant(3, 3, 0);
    Eigen::Matrix3d dCTT = Eigen::Matrix3d::Constant(3, 3, 0);

    Eigen::Matrix3d pstress = PrincipalStress(stress, tX1, pX1, tX3)[0];
    Eigen::Matrix3d dCp = PrincipalStress(stress, tX1, pX1, tX3)[1];

    for (int i = 0; i < 3; i++) {
        stress(i, i) = pstress(i, 0);
    }

    Eigen::Vector3d p = { 0,0,0 };

    p(0) = SphToCart(strike, dip, 1)[0];
    p(1) = SphToCart(strike, dip, 1)[1];
    p(2) = SphToCart(strike, dip, 1)[2];

    Eigen::Vector3d pT = { 0,0,0 };

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            pT(i) = dCp(i, j) * p(j) + pT(i);
        }
    }

    Eigen::Vector3d T = { 0,0,0 };

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            T(i) = stress(i, j) * pT(j) + T(i);
        }
    }

    Eigen::Vector3d B = { 0,0,0 };
    B(0) = T(1) * pT(2) - T(2) * pT(1);
    B(1) = T(2) * pT(0) - T(0) * pT(2);
    B(2) = T(0) * pT(1) - T(1) * pT(0);

    Eigen::Vector3d S = { 0,0,0 };
    S(0) = pT(1) * B(2) - pT(2) * B(1);
    S(1) = pT(2) * B(0) - pT(0) * B(2);
    S(2) = pT(0) * B(1) - pT(1) * B(0);

    double rB = sqrt(B(0) * B(0) + B(1) * B(1) + B(2) * B(2));
    double rS = sqrt(B(0) * B(0) + B(1) * B(1) + B(2) * B(2));

    for (int i = 0; i < 3; i++) {
        B(i) = B(i) / rB;
        S(i) = S(i) / rS;
    }

    Eigen::Matrix3d a = Eigen::Matrix3d::Constant(3, 3, 0);
    a(0, 0) = pT(0);
    a(0, 1) = pT(1);
    a(0, 2) = pT(2);

    a(1, 0) = B(0);
    a(1, 1) = B(1);
    a(1, 2) = B(2);

    a(2, 0) = S(0);
    a(2, 1) = S(1);
    a(2, 2) = S(2);

    for (int i = 0; i < 3; i++) {
        tractions(i, 0) = stress(0, 0) * a(0, 0) * a(i, 0) + stress(1, 1) * a(0, 1) * a(i, 1) + stress(2, 2) * a(i, 2) * a(i, 2);
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            dCTT(0, i) = dCp(j, i) * pT(j) + dCTT(0, i);
            dCTT(1, i) = dCp(j, i) * B(j) + dCTT(1, i);
            dCTT(2, i) = dCp(j, i) * S(j) + dCTT(2, i);
        }
    }

    tractions(0, 1) = CartToSph(dCTT(0, 0), dCTT(0, 1), dCTT(0, 2))[0];
    tractions(0, 2) = CartToSph(dCTT(0, 0), dCTT(0, 1), dCTT(0, 2))[1];

    tractions(1, 1) = CartToSph(dCTT(1, 0), dCTT(1, 1), dCTT(1, 2))[0];
    tractions(1, 2) = CartToSph(dCTT(1, 0), dCTT(1, 1), dCTT(1, 2))[1];

    tractions(2, 1) = CartToSph(dCTT(2, 0), dCTT(2, 1), dCTT(2, 2))[0];
    tractions(2, 2) = CartToSph(dCTT(2, 0), dCTT(2, 1), dCTT(2, 2))[1];

    return tractions;

}
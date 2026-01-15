#pragma once
#include <vector>
#include <cmath>
#include <algorithm>

inline double s_curve(double r){
    return 10*std::pow(r,3) - 15*std::pow(r,4) + 6*std::pow(r,5);
}

inline double s_curve_d1(double r){
    // ds/dr
    return 30*std::pow(r,2) - 60*std::pow(r,3) + 30*std::pow(r,4);
}

inline double s_curve_d2(double r){
    // d2s/dr2
    return 60*r - 180*std::pow(r,2) + 120*std::pow(r,3);
}

inline double s_curve_d3(double r){
    // d3s/dr3
    return 60 - 360*r + 360*std::pow(r,2);
}

inline std::vector<std::vector<double>> plan_minjerk(
    const std::vector<double>& q0,
    const std::vector<double>& q1,
    double T, double dt)
{
    size_t dof = q0.size();
    int N = std::max(2, (int)std::round(T/dt));
    std::vector<std::vector<double>> out;
    out.reserve(N+1);

    for(int k=0;k<=N;k++){
        double t = k*dt;
        double r = std::clamp(t / std::max(T, 1e-9), 0.0, 1.0);
        double s = s_curve(r);
        std::vector<double> row(1 + (int)dof);
        row[0] = t;
        for(size_t i=0;i<dof;i++){
            row[1 + (int)i] = q0[i] + (q1[i]-q0[i]) * s;
        }
        out.push_back(std::move(row));
    }
    return out;
}


struct PMPPoint {
    double t;
    std::vector<double> q;
    std::vector<double> dq;
    std::vector<double> ddq;
    std::vector<double> u;   
};

inline std::vector<PMPPoint> plan_pmp_minimum_jerk(
    const std::vector<double>& q0,
    const std::vector<double>& q1,
    double T, double dt)
{
    size_t dof = q0.size();
    int N = std::max(2, (int)std::round(T/dt));
    std::vector<PMPPoint> out;
    out.reserve(N+1);

    for(int k=0;k<=N;k++){
        double t = k*dt;
        double r = std::clamp(t / std::max(T, 1e-9), 0.0, 1.0);

        double s   = s_curve(r);
        double s1  = s_curve_d1(r);
        double s2  = s_curve_d2(r);
        double s3  = s_curve_d3(r);

        PMPPoint p;
        p.t = t;
        p.q.resize(dof);
        p.dq.resize(dof);
        p.ddq.resize(dof);
        p.u.resize(dof);

        for(size_t i=0;i<dof;i++){
            double dq = q1[i] - q0[i];

            p.q[i]   = q0[i] + dq * s;
            p.dq[i]  = dq * s1 / T;
            p.ddq[i] = dq * s2 / (T*T);
            p.u[i]   = dq * s3 / (T*T*T);  
        }
        out.push_back(std::move(p));
    }
    return out;
}



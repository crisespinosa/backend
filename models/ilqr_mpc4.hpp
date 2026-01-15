#pragma once
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>

struct ILQRWeights4 {
    double Q_pos = 50.0;
    double Q_vel = 3.0;
    double R_u   = 0.2;
    double Qf_pos = 200.0;
    double Qf_vel = 20.0;
};

struct ILQROptions4 {
    int N = 40;
    int iters = 15;
    double dt = 0.02;
    double umax = 25.0;
    double tol = 1e-6;
};

struct ILQRResult4 {
    double dt{};
    std::vector<double> t;                       // N+1
    std::vector<std::array<double,4>> q;         // N+1
    std::vector<std::array<double,4>> dq;        // N+1
    std::vector<std::array<double,4>> u;         // N
    double final_cost{};
};

namespace ilqr4 {

using Vec8 = std::array<double,8>;   // [q0..q3, v0..v3]
using Vec4 = std::array<double,4>;

inline double clamp(double v, double lo, double hi) {
    return std::max(lo, std::min(hi, v));
}
inline Vec4 clamp_u(const Vec4& u, double umax) {
    return { clamp(u[0], -umax, umax),
             clamp(u[1], -umax, umax),
             clamp(u[2], -umax, umax),
             clamp(u[3], -umax, umax) };
}

inline Vec8 f_step(const Vec8& x, const Vec4& u, double dt) {
    Vec8 xn = x;
    for(int i=0;i<4;i++){
        xn[i]   = x[i]   + dt * x[4+i];   // q_next
        xn[4+i] = x[4+i] + dt * u[i];     // v_next
    }
    return xn;
}

inline double stage_cost(const Vec8& x, const Vec4& u, const Vec4& q_ref, const ILQRWeights4& w, double dt) {
    double c = 0.0;
    for(int i=0;i<4;i++){
        double e = x[i] - q_ref[i];
        c += w.Q_pos*e*e;
        double v = x[4+i];
        c += w.Q_vel*v*v;
        c += w.R_u*u[i]*u[i];
    }
    return c * dt;
}

inline double terminal_cost(const Vec8& x, const Vec4& q_ref, const ILQRWeights4& w) {
    double c = 0.0;
    for(int i=0;i<4;i++){
        double e = x[i] - q_ref[i];
        c += w.Qf_pos*e*e;
        double v = x[4+i];
        c += w.Qf_vel*v*v;
    }
    return c;
}

inline double rollout_cost(
    const Vec8& x0,
    const std::vector<Vec4>& u_seq,
    const Vec4& q_ref,
    const ILQROptions4& opt,
    const ILQRWeights4& w,
    std::vector<Vec8>* x_out=nullptr
){
    Vec8 x=x0;
    double J=0.0;
    if(x_out){ x_out->clear(); x_out->reserve(opt.N+1); x_out->push_back(x); }
    for(int k=0;k<opt.N;k++){
        Vec4 uk = clamp_u(u_seq[k], opt.umax);
        J += stage_cost(x, uk, q_ref, w, opt.dt);
        x = f_step(x, uk, opt.dt);
        if(x_out) x_out->push_back(x);
    }
    J += terminal_cost(x, q_ref, w);
    return J;
}

inline ILQRResult4 solve_mpc_ilqr(
    const Vec4& q0,
    const Vec4& dq0,
    const Vec4& q_ref,
    const ILQROptions4& opt,
    const ILQRWeights4& w
){
    // init u=0
    std::vector<Vec4> u(opt.N, Vec4{0,0,0,0});
    Vec8 x0{ q0[0],q0[1],q0[2],q0[3], dq0[0],dq0[1],dq0[2],dq0[3] };

    std::vector<Vec8> x_traj;
    double J = rollout_cost(x0, u, q_ref, opt, w, &x_traj);

    for(int it=0; it<opt.iters; ++it){
        
        std::array<double,8> Pdiag = { w.Qf_pos,w.Qf_pos,w.Qf_pos,w.Qf_pos, w.Qf_vel,w.Qf_vel,w.Qf_vel,w.Qf_vel };
        Vec8 p{};
        {
            const Vec8& xT = x_traj.back();
            for(int i=0;i<4;i++){
                p[i]   = 2.0*w.Qf_pos*(xT[i]-q_ref[i]);
                p[4+i] = 2.0*w.Qf_vel*(xT[4+i]);
            }
        }

        std::vector<std::array<double,32>> K(opt.N); // 4x8 = 32
        std::vector<Vec4> kff(opt.N, Vec4{0,0,0,0});

        for(int kidx=opt.N-1; kidx>=0; --kidx){
            const Vec8& xk = x_traj[kidx];
            Vec4 uk = clamp_u(u[kidx], opt.umax);

            Vec8 lx{};
            for(int i=0;i<4;i++){
                lx[i]   = 2.0*w.Q_pos*(xk[i]-q_ref[i]) * opt.dt;
                lx[4+i] = 2.0*w.Q_vel*(xk[4+i]) * opt.dt;
            }
            Vec4 lu{ 2.0*w.R_u*uk[0]*opt.dt, 2.0*w.R_u*uk[1]*opt.dt, 2.0*w.R_u*uk[2]*opt.dt, 2.0*w.R_u*uk[3]*opt.dt };

            Vec4 Qu{
                lu[0] + opt.dt * p[4],
                lu[1] + opt.dt * p[5],
                lu[2] + opt.dt * p[6],
                lu[3] + opt.dt * p[7]
            };

            Vec8 Qx{};
            for(int i=0;i<4;i++){
                Qx[i]   = lx[i] + p[i];
                Qx[4+i] = lx[4+i] + opt.dt*p[i] + p[4+i];
            }

            Vec4 Quu{
                2.0*w.R_u*opt.dt + (opt.dt*opt.dt)*Pdiag[4],
                2.0*w.R_u*opt.dt + (opt.dt*opt.dt)*Pdiag[5],
                2.0*w.R_u*opt.dt + (opt.dt*opt.dt)*Pdiag[6],
                2.0*w.R_u*opt.dt + (opt.dt*opt.dt)*Pdiag[7]
            };

            std::array<double,32> Qux{}; // 4x8
            Qux[0*8 + 4] = opt.dt * Pdiag[4];
            Qux[1*8 + 5] = opt.dt * Pdiag[5];
            Qux[2*8 + 6] = opt.dt * Pdiag[6];
            Qux[3*8 + 7] = opt.dt * Pdiag[7];

            Vec4 k_local{
                -Qu[0]/Quu[0],
                -Qu[1]/Quu[1],
                -Qu[2]/Quu[2],
                -Qu[3]/Quu[3]
            };

            std::array<double,32> K_local{};
            for(int r=0;r<4;r++){
                double inv = -1.0/Quu[r];
                for(int c=0;c<8;c++){
                    K_local[r*8 + c] = inv * Qux[r*8 + c];
                }
            }

            kff[kidx]=k_local;
            K[kidx]=K_local;

            p = Qx;
            for(int i=0;i<4;i++){
                Pdiag[i]   = (2.0*w.Q_pos*opt.dt) + Pdiag[i];
                Pdiag[4+i] = (2.0*w.Q_vel*opt.dt) + Pdiag[4+i] + (opt.dt*opt.dt)*Pdiag[i];
            }
        }

        double bestJ = J;
        std::vector<Vec4> bestU = u;

        const double alphas[] = {1.0,0.5,0.25,0.1};
        for(double a: alphas){
            Vec8 x = x0;
            std::vector<Vec4> u_new(opt.N);

            for(int kidx=0;kidx<opt.N;kidx++){
                Vec8 dx{};
                for(int i=0;i<8;i++) dx[i]= x[i]-x_traj[kidx][i];

                const auto& Kk = K[kidx];
                Vec4 Ku{};
                for(int r=0;r<4;r++){
                    double s=0.0;
                    for(int c=0;c<8;c++) s += Kk[r*8+c]*dx[c];
                    Ku[r]=s;
                }

                Vec4 u_try{
                    u[kidx][0] + a*kff[kidx][0] + Ku[0],
                    u[kidx][1] + a*kff[kidx][1] + Ku[1],
                    u[kidx][2] + a*kff[kidx][2] + Ku[2],
                    u[kidx][3] + a*kff[kidx][3] + Ku[3]
                };
                u_try = clamp_u(u_try, opt.umax);
                u_new[kidx]=u_try;
                x = f_step(x, u_try, opt.dt);
            }

            double Jnew = rollout_cost(x0, u_new, q_ref, opt, w, nullptr);
            if(Jnew < bestJ){ bestJ=Jnew; bestU=std::move(u_new); }
        }

        double improvement = J - bestJ;
        u = std::move(bestU);
        J = bestJ;
        rollout_cost(x0, u, q_ref, opt, w, &x_traj);
        if(improvement < opt.tol) break;
    }

    ILQRResult4 res;
    res.dt = opt.dt;
    res.final_cost = J;
    res.t.resize(opt.N+1);
    res.q.resize(opt.N+1);
    res.dq.resize(opt.N+1);
    res.u.resize(opt.N);

    for(int k=0;k<=opt.N;k++){
        res.t[k]=k*opt.dt;
        res.q[k]={ x_traj[k][0],x_traj[k][1],x_traj[k][2],x_traj[k][3] };
        res.dq[k]={ x_traj[k][4],x_traj[k][5],x_traj[k][6],x_traj[k][7] };
        if(k<opt.N) res.u[k]=clamp_u(u[k], opt.umax);
    }
    return res;
}

} // namespace ilqr4

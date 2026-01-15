#include "ArmController.h"
#include <drogon/HttpAppFramework.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <json/json.h>
#include <iostream>

#include "trajectory.hpp"         // plan_minjerk, plan_pmp_minimum_jerk
#include "models/ilqr_mpc4.hpp"    // iLQR/MPC 4DOF

using namespace drogon;


static void ik_2link_3angles(double x, double y,
                            double &a1, double &a2, double &a3)
{
    double L = 1.0;
    double D = std::sqrt(x*x + y*y);

  
    if (D > 2 * L) {
        x *= (2 * L / D);
        y *= (2 * L / D);
        D = 2 * L;
    }

    double cos_elbow = (x*x + y*y - 2*L*L) / (2*L*L);
    cos_elbow = std::max(-1.0, std::min(1.0, cos_elbow));
    double elbow = std::acos(cos_elbow);

    double k1 = L + L * cos_elbow;
    double k2 = L * std::sin(elbow);
    double shoulder = std::atan2(y, x) - std::atan2(k2, k1);

    double wrist = -shoulder - elbow;

    a1 = shoulder;
    a2 = elbow;
    a3 = wrist;
}

// ----- CONSTRUCTOR -----
ArmController::ArmController()
    : dyn_(4)  
{
    dyn_.setState({0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0});
}


void ArmController::handleArmInput(const HttpRequestPtr &req,
                                   std::function<void (const HttpResponsePtr &)> &&callback)
{
    auto json = req->getJsonObject();
    if (!json || !json->isMember("x") || !json->isMember("y")) {
        auto resp = HttpResponse::newHttpJsonResponse(Json::Value("Not enough parameters x or y"));
        resp->setStatusCode(k400BadRequest);
        callback(resp);
        return;
    }

    double x = (*json)["x"].asDouble();
    double y = (*json)["y"].asDouble();

    double a1, a2, a3;
    ik_2link_3angles(x, y, a1, a2, a3);

    std::cout << "Angles IK:\n";
    std::cout << "  a1: " << a1 << " rad\n";
    std::cout << "  a2: " << a2 << " rad\n";
    std::cout << "  a3: " << a3 << " rad\n";

    Json::Value result;
    result["x"] = x;
    result["y"] = y;
    result["shoulder_angle"] = a1;
    result["elbow_angle"]    = a2;
    result["wrist_angle"]    = a3;

    auto resp = HttpResponse::newHttpJsonResponse(result);
    callback(resp);
}


void ArmController::handlePlan(const HttpRequestPtr &req,
                               std::function<void (const HttpResponsePtr &)> &&callback)
{
    auto json = req->getJsonObject();
    if (!json || !json->isMember("x") || !json->isMember("y")) {
        auto resp = HttpResponse::newHttpJsonResponse(Json::Value("Not enough parameters x or y"));
        resp->setStatusCode(k400BadRequest);
        callback(resp);
        return;
    }

    double x  = (*json)["x"].asDouble();
    double y  = (*json)["y"].asDouble();
    double T  = json->isMember("T")  ? (*json)["T"].asDouble()  : 1.0;
    double dt = json->isMember("dt") ? (*json)["dt"].asDouble() : 0.02;

    double a1, a2, a3;
    ik_2link_3angles(x, y, a1, a2, a3);

    std::vector<double> q_target = { a1, a2, a3 };

    auto q0 = dyn_.state().q;
    if (q0.size() < 3) q0 = {0.0, 0.0, 0.0, 0.0};

    std::vector<double> q0_3 = { q0[0], q0[1], q0[2] };

    auto table = plan_minjerk(q0_3, q_target, T, dt);

    auto st = dyn_.state();
    std::vector<double> q4 = st.q;
    std::vector<double> dq4 = st.dq;
    if (q4.size() < 4) q4 = {0,0,0,0};
    if (dq4.size() < 4) dq4 = {0,0,0,0};

    q4[0]=a1; q4[1]=a2; q4[2]=a3;
    dq4[0]=0; dq4[1]=0; dq4[2]=0;
    dyn_.setState(q4, dq4);

    Json::Value out(Json::objectValue);
    out["x"]  = x;
    out["y"]  = y;
    out["dt"] = dt;
    out["trajectory"] = Json::arrayValue;

    for (const auto &row : table) {
        Json::Value item(Json::objectValue);
        item["t"] = row[0];

        Json::Value q(Json::arrayValue);
        for (size_t i = 1; i < row.size(); ++i) q.append(row[i]);

        item["q"] = q;
        out["trajectory"].append(item);
    }

    auto resp = HttpResponse::newHttpJsonResponse(out);
    callback(resp);
}


void ArmController::handlePlanPMP(const HttpRequestPtr &req,
                                  std::function<void (const HttpResponsePtr &)> &&callback)
{
    auto json = req->getJsonObject();
    if (!json || !json->isMember("x") || !json->isMember("y")) {
        auto resp = HttpResponse::newHttpJsonResponse(Json::Value("Not enough parameters x or y"));
        resp->setStatusCode(k400BadRequest);
        callback(resp);
        return;
    }

    double x  = (*json)["x"].asDouble();
    double y  = (*json)["y"].asDouble();
    double T  = json->isMember("T")  ? (*json)["T"].asDouble()  : 1.0;
    double dt = json->isMember("dt") ? (*json)["dt"].asDouble() : 0.02;

    double a1, a2, a3;
    ik_2link_3angles(x, y, a1, a2, a3);

    std::vector<double> q_target = { a1, a2, a3 };

    auto st = dyn_.state();
    if (st.q.size() < 4) dyn_.setState({0,0,0,0}, {0,0,0,0});

    auto q0 = dyn_.state().q;
    std::vector<double> q0_3 = { q0[0], q0[1], q0[2] };

    auto pmp_traj = plan_pmp_minimum_jerk(q0_3, q_target, T, dt);


    st = dyn_.state();
    std::vector<double> q4 = st.q;
    std::vector<double> dq4 = st.dq;
    if (q4.size() < 4) q4 = {0,0,0,0};
    if (dq4.size() < 4) dq4 = {0,0,0,0};

    q4[0]=a1; q4[1]=a2; q4[2]=a3;
    dq4[0]=0; dq4[1]=0; dq4[2]=0;
    dyn_.setState(q4, dq4);

    Json::Value out(Json::objectValue);
    out["x"]  = x;
    out["y"]  = y;
    out["dt"] = dt;
    out["trajectory"] = Json::arrayValue;

    for (const auto &p : pmp_traj) {
        Json::Value item(Json::objectValue);
        item["t"] = p.t;

        Json::Value q_json(Json::arrayValue);
        Json::Value dq_json(Json::arrayValue);
        Json::Value ddq_json(Json::arrayValue);
        Json::Value u_json(Json::arrayValue);

        for (size_t i = 0; i < p.q.size(); ++i) {
            q_json.append(p.q[i]);
            dq_json.append(p.dq[i]);
            ddq_json.append(p.ddq[i]);
            u_json.append(p.u[i]);
        }

        item["q"]   = q_json;
        item["dq"]  = dq_json;
        item["ddq"] = ddq_json;
        item["u"]   = u_json;

        out["trajectory"].append(item);
    }

    auto resp = HttpResponse::newHttpJsonResponse(out);
    callback(resp);
}


void ArmController::handleMPCiLQR4(const HttpRequestPtr &req,
                                  std::function<void (const HttpResponsePtr &)> &&callback)
{
    auto json = req->getJsonObject();
    if (!json || !json->isMember("x") || !json->isMember("y")) {
        auto resp = HttpResponse::newHttpJsonResponse(Json::Value("Not enough parameters x or y"));
        resp->setStatusCode(k400BadRequest);
        callback(resp);
        return;
    }

    double x  = (*json)["x"].asDouble();
    double y  = (*json)["y"].asDouble();
    int N     = json->isMember("N")     ? (*json)["N"].asInt()    : 40;
    int iters = json->isMember("iters") ? (*json)["iters"].asInt(): 15;
    double dt = json->isMember("dt")    ? (*json)["dt"].asDouble(): 0.02;

    // IK -> 3 targets
    double a1, a2, a3;
    ik_2link_3angles(x, y, a1, a2, a3);

   
    auto st = dyn_.state();
    if (st.q.size() < 4 || st.dq.size() < 4) {
        dyn_.setState({0,0,0,0},{0,0,0,0});
        st = dyn_.state();
    }

    ilqr4::Vec4 q0  { st.q[0],  st.q[1],  st.q[2],  st.q[3]  };
    ilqr4::Vec4 dq0 { st.dq[0], st.dq[1], st.dq[2], st.dq[3] };

    
    double roll_keep = q0[3];
    ilqr4::Vec4 q_ref { a1, a2, a3, roll_keep };

    ILQROptions4 opt;
    opt.N = N;
    opt.iters = iters;
    opt.dt = dt;
    opt.umax = 25.0;

    ILQRWeights4 w; // default

    auto sol = ilqr4::solve_mpc_ilqr(q0, dq0, q_ref, opt, w);

    if (!sol.u.empty()) {
        ilqr4::Vec8 x0s { q0[0],q0[1],q0[2],q0[3], dq0[0],dq0[1],dq0[2],dq0[3] };
        auto x1 = ilqr4::f_step(x0s, sol.u[0], dt);
        dyn_.setState({x1[0],x1[1],x1[2],x1[3]}, {x1[4],x1[5],x1[6],x1[7]});
    }

    Json::Value out(Json::objectValue);
    out["x"] = x;
    out["y"] = y;
    out["dt"] = dt;
    out["cost"] = sol.final_cost;
    out["trajectory"] = Json::arrayValue;

    for (int k=0; k<=opt.N; k++) {
        Json::Value item(Json::objectValue);
        item["t"] = sol.t[k];

        Json::Value qj(Json::arrayValue), dqj(Json::arrayValue);
        for (int i=0;i<4;i++){ qj.append(sol.q[k][i]); dqj.append(sol.dq[k][i]); }
        item["q"] = qj;
        item["dq"] = dqj;

        if (k < opt.N) {
            Json::Value uj(Json::arrayValue);
            for (int i=0;i<4;i++) uj.append(sol.u[k][i]);
            item["u"] = uj;
        }

        out["trajectory"].append(item);
    }

    auto resp = HttpResponse::newHttpJsonResponse(out);
    callback(resp);
}


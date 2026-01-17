#include "ArmController.h"
#include <drogon/HttpAppFramework.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <json/json.h>
#include <iostream>

#include "trajectory.hpp"         // plan_minjerk, plan_pmp_minimum_jerk
#include "models/ilqr_mpc4.hpp"   // iLQR/MPC 4DOF

using namespace drogon;

// ------------------------------------------------------------
// Helper: vector q (tam 3/4/6) -> JSON q[6]
// ------------------------------------------------------------
static Json::Value to_q6_json(const std::vector<double> &q_in)
{
    Json::Value q(Json::arrayValue);
    for (int i = 0; i < 6; ++i) {
        double v = (i < (int)q_in.size()) ? q_in[i] : 0.0;
        q.append(v);
    }
    return q;
}

// ------------------------------------------------------------
// IK simple 2D (3 ángulos) - (tu ejemplo actual)
// ------------------------------------------------------------
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

// ------------------------------------------------------------
// POST /arm   -> devuelve q[6] directo (IK)
// ------------------------------------------------------------
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

    std::cout << "Angles IK:\n"
              << "  a1: " << a1 << " rad\n"
              << "  a2: " << a2 << " rad\n"
              << "  a3: " << a3 << " rad\n";

    // q1..q3 de IK; q4..q6 por ahora 0
    std::vector<double> q_vec = {a1, a2, a3, 0.0, 0.0, 0.0};

    Json::Value result(Json::objectValue);
    result["x"] = x;
    result["y"] = y;
    result["q"] = to_q6_json(q_vec);
    result["unit"] = "rad";

    auto resp = HttpResponse::newHttpJsonResponse(result);
    callback(resp);
}

// ------------------------------------------------------------
// POST /arm/plan  -> trajectory con q[6] (min-jerk)
// body: {x,y,T,dt}
// ------------------------------------------------------------
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

    // actualiza estado interno (solo 3 dof en este ejemplo)
    auto st = dyn_.state();
    std::vector<double> q4  = st.q;
    std::vector<double> dq4 = st.dq;
    if (q4.size()  < 4) q4  = {0,0,0,0};
    if (dq4.size() < 4) dq4 = {0,0,0,0};

    q4[0]=a1; q4[1]=a2; q4[2]=a3;
    dq4[0]=0; dq4[1]=0; dq4[2]=0;
    dyn_.setState(q4, dq4);

    Json::Value out(Json::objectValue);
    out["x"]  = x;
    out["y"]  = y;
    out["dt"] = dt;
    out["unit"] = "rad";
    out["trajectory"] = Json::arrayValue;

    // row: [t, q1, q2, q3]
    for (const auto &row : table) {
        Json::Value item(Json::objectValue);
        item["t"] = row[0];

        std::vector<double> q3;
        for (size_t i = 1; i < row.size(); ++i) q3.push_back(row[i]);

        item["q"] = to_q6_json(q3); // -> q[6]
        out["trajectory"].append(item);
    }

    auto resp = HttpResponse::newHttpJsonResponse(out);
    callback(resp);
}

// ------------------------------------------------------------
// POST /arm/plan_pmp -> trajectory con q,dq,ddq,u (todos q[6])
// body: {x,y,T,dt}
// ------------------------------------------------------------
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

    // actualiza estado interno
    st = dyn_.state();
    std::vector<double> q4  = st.q;
    std::vector<double> dq4 = st.dq;
    if (q4.size()  < 4) q4  = {0,0,0,0};
    if (dq4.size() < 4) dq4 = {0,0,0,0};

    q4[0]=a1; q4[1]=a2; q4[2]=a3;
    dq4[0]=0; dq4[1]=0; dq4[2]=0;
    dyn_.setState(q4, dq4);

    Json::Value out(Json::objectValue);
    out["x"]  = x;
    out["y"]  = y;
    out["dt"] = dt;
    out["unit"] = "rad";
    out["trajectory"] = Json::arrayValue;

    for (const auto &p : pmp_traj) {
        Json::Value item(Json::objectValue);
        item["t"] = p.t;

        item["q"]   = to_q6_json(p.q);
        item["dq"]  = to_q6_json(p.dq);
        item["ddq"] = to_q6_json(p.ddq);
        item["u"]   = to_q6_json(p.u);

        out["trajectory"].append(item);
    }

    auto resp = HttpResponse::newHttpJsonResponse(out);
    callback(resp);
}
void ArmController::handlePlanPMPQ(const HttpRequestPtr &req,
                                   std::function<void (const HttpResponsePtr &)> &&callback)
{
    auto json = req->getJsonObject();
    if (!json || !json->isMember("q_target") || !(*json)["q_target"].isArray()) {
        auto resp = HttpResponse::newHttpJsonResponse(
            Json::Value("Not enough parameters: q_target (array)"));
        resp->setStatusCode(k400BadRequest);
        callback(resp);
        return;
    }

    // Trajectory params
    double T  = json->isMember("T")  ? (*json)["T"].asDouble()  : 1.0;
    double dt = json->isMember("dt") ? (*json)["dt"].asDouble() : 0.02;

    // Read q_target[6] (rad)
    std::vector<double> q_target6;
    for (const auto &v : (*json)["q_target"]) q_target6.push_back(v.asDouble());
    if (q_target6.size() < 6) q_target6.resize(6, 0.0);

    // Current state from dyn_ (we store at least 4 in your project)
    auto st = dyn_.state();
    if (st.q.size() < 4) {
        dyn_.setState({0,0,0,0}, {0,0,0,0});
        st = dyn_.state();
    }

    // We plan PMP for first 3 joints (because your PMP planner is 3-DOF now)
    std::vector<double> q0_3 = { st.q[0], st.q[1], st.q[2] };
    std::vector<double> qT_3 = { q_target6[0], q_target6[1], q_target6[2] };

    auto pmp_traj = plan_pmp_minimum_jerk(q0_3, qT_3, T, dt);

    // Update internal state to last target (store in dyn_)
    std::vector<double> q4  = st.q;   if (q4.size()  < 4) q4  = {0,0,0,0};
    std::vector<double> dq4 = st.dq;  if (dq4.size() < 4) dq4 = {0,0,0,0};

    q4[0] = qT_3[0];
    q4[1] = qT_3[1];
    q4[2] = qT_3[2];
    dq4[0]=0; dq4[1]=0; dq4[2]=0;
    dyn_.setState(q4, dq4);

    // Build response
    Json::Value out(Json::objectValue);
    out["dt"] = dt;
    out["unit"] = "rad";
    out["trajectory"] = Json::arrayValue;

    for (const auto &p : pmp_traj) {
        Json::Value item(Json::objectValue);
        item["t"] = p.t;

        // Compose full q[6]:
        // First 3 from PMP, last 3 from q_target6 (so wrists follow target constant)
        std::vector<double> q6 = { p.q[0], p.q[1], p.q[2], q_target6[3], q_target6[4], q_target6[5] };
        item["q"] = to_q6_json(q6);

        out["trajectory"].append(item);
    }

    auto resp = HttpResponse::newHttpJsonResponse(out);
    callback(resp);
}

// ------------------------------------------------------------
// POST /arm/mpc_ilqr4 (opcional) -> trajectory q[6], dq[6], u[6]
// ------------------------------------------------------------
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
    out["unit"] = "rad";
    out["cost"] = sol.final_cost;
    out["trajectory"] = Json::arrayValue;

    for (int k=0; k<=opt.N; k++) {
        Json::Value item(Json::objectValue);
        item["t"] = sol.t[k];

        std::vector<double> q4v, dq4v;
        for (int i=0;i<4;i++){ q4v.push_back(sol.q[k][i]); dq4v.push_back(sol.dq[k][i]); }

        item["q"]  = to_q6_json(q4v);
        item["dq"] = to_q6_json(dq4v);

        if (k < opt.N) {
            std::vector<double> u4v;
            for (int i=0;i<4;i++) u4v.push_back(sol.u[k][i]);
            item["u"] = to_q6_json(u4v);
        }

        out["trajectory"].append(item);
    }

    auto resp = HttpResponse::newHttpJsonResponse(out);
    callback(resp);
}


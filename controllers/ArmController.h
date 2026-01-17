#pragma once

#include <drogon/HttpController.h>
#include <functional>
#include "dynamics.hpp"   // SimpleDynamics

class ArmController : public drogon::HttpController<ArmController> {
public:
    ArmController();

    METHOD_LIST_BEGIN
        ADD_METHOD_TO(ArmController::handleArmInput,   "/arm",       drogon::Post);
        ADD_METHOD_TO(ArmController::handlePlan,       "/arm/plan",      drogon::Post);
        ADD_METHOD_TO(ArmController::handlePlanPMP,    "/arm/plan_pmp",  drogon::Post);
        ADD_METHOD_TO(ArmController::handlePlanPMPQ,   "/arm/plan_pmp_q",drogon::Post);
        ADD_METHOD_TO(ArmController::handleMPCiLQR4,   "/arm/mpc_ilqr4", drogon::Post);
    METHOD_LIST_END

    void handleArmInput(const drogon::HttpRequestPtr &,
                        std::function<void (const drogon::HttpResponsePtr &)> &&);

    void handlePlan(const drogon::HttpRequestPtr &,
                    std::function<void (const drogon::HttpResponsePtr &)> &&);

    void handlePlanPMP(const drogon::HttpRequestPtr &,
                       std::function<void (const drogon::HttpResponsePtr &)> &&);
    void handlePlanPMPQ(const drogon::HttpRequestPtr &,
                    std::function<void (const drogon::HttpResponsePtr &)> &&);
                    
    void handleMPCiLQR4(const drogon::HttpRequestPtr &,
                        std::function<void (const drogon::HttpResponsePtr &)> &&);

private:
    SimpleDynamics dyn_;  
};



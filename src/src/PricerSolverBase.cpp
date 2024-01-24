// MIT License

// Copyright (c) 2021 Daniel Kowalczyk

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "PricerSolverBase.hpp"
#include <algorithm>                    // for min, __fill_fn
#include <cmath>                        // for fabs
#include <limits>                       // for numeric_limits
#include <memory>                       // for __shared_ptr_a...
#include <range/v3/view/enumerate.hpp>  // for enumerate_fn
#include <range/v3/view/zip.hpp>        // for zip
#include <span>                         // for span
#include <vector>                       // for vector
#include "Column.h"                     // for Column
#include "Instance.h"                   // for Instance
#include "gurobi_c++.h"                 // for GRBModel, GRBEnv
#include "gurobi_c.h"                   // for GRB_INFEASIBLE

PricerSolverBase::PricerSolverBase(const Instance& instance)
    : jobs(instance.jobs),
      convex_constr_id(instance.nb_jobs),
      convex_rhs(instance.nb_machines),
      env(genv),
      model(*env),
      reformulation_model(instance.nb_jobs, instance.nb_machines),
      is_integer_solution(false),
      constLB(0.0),
      UB(std::numeric_limits<int>::max()),
      x_bar(std::vector<std::vector<BddCoeff>>(instance.nb_jobs,
                                               std::vector<BddCoeff>())) {
    try {
        model.set(GRB_IntParam_Method, GRB_METHOD_AUTO);
        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        model.set(GRB_IntParam_Presolve, GRB_PRESOLVE_AGGRESSIVE);
    } catch (const GRBException& e) {
        fmt::print("Error code = {}\n", e.getErrorCode());
        fmt::print("{}", e.getMessage());
    } catch (...) {
        fmt::print("Exception during optimization\n");
    }
}

PricerSolverBase::PricerSolverBase(const PricerSolverBase& other)
    : jobs(other.jobs),
      convex_constr_id(other.convex_constr_id),
      convex_rhs(other.convex_rhs),
      problem_name(other.problem_name),
      env(other.env),
      model(*genv),
      reformulation_model(other.reformulation_model),
      is_integer_solution(other.is_integer_solution),
      constLB(other.constLB),
      UB(other.UB),
      x_bar(other.x_bar) {
    model.set(GRB_IntParam_Method, GRB_METHOD_AUTO);
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model.set(GRB_IntParam_Presolve, GRB_PRESOLVE_AGGRESSIVE);
}

PricerSolverBase::~PricerSolverBase() = default;

auto PricerSolverBase::add_constraints() -> int {
    int val = 0;

    return val;
}

auto PricerSolverBase::evaluate_mip_model() -> bool {
    int opt_status = model.get(GRB_IntAttr_Status);

    double obj_val{};
    switch (opt_status) {
        case GRB_OPTIMAL:
            obj_val = model.get(GRB_DoubleAttr_ObjVal);
            update_UB(obj_val);
            return obj_val < UB;
        case GRB_INFEASIBLE:
        case GRB_INF_OR_UNBD:
        case GRB_UNBOUNDED:
            return false;
        default:
            return true;
    }
}


auto PricerSolverBase::compute_sub_optimal_duals(
    const std::span<const double>&              lambda,
    const std::vector<std::shared_ptr<Column>>& columns) -> bool {
    GRBModel sub_optimal(*genv);
    // sub_optimal.set(GRB_IntParam_OutputFlag, 1);
    sub_optimal.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    std::vector<GRBVar> beta{convex_constr_id};
    std::vector<GRBVar> eta;
    double              LB{};
    auto                removed = false;

    for (auto& it : beta) {
        it = sub_optimal.addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
    }

    auto last = sub_optimal.addVar(
        0.0, GRB_INFINITY, -static_cast<double>(convex_rhs), GRB_CONTINUOUS);

    // std::span<const double> aux_cols{lambda, columns.size()};

    for (auto&& [set, x] : ranges::views::zip(columns, lambda)) {
        if (x > EPS_SOLVER) {
            // auto* tmp = columns[i].get();
            eta.emplace_back(sub_optimal.addVar(0.0, GRB_INFINITY, 1.0, 'C'));
            GRBLinExpr expr = -last;
            for (auto& it : set->job_list) {
                expr += beta[it->job];
            }
            expr += eta.back();
            sub_optimal.addConstr(expr, '=',
                                  set->total_weighted_completion_time);
            LB += set->total_weighted_completion_time * x;
        } else {
            // auto*      tmp = columns[i].get();
            GRBLinExpr expr = -last;
            for (auto& it : set->job_list) {
                expr += beta[it->job];
            }
            sub_optimal.addConstr(expr, '<',
                                  set->total_weighted_completion_time);
        }
    }

    GRBLinExpr expr = -static_cast<double>(convex_rhs) * last;
    for (auto& it : beta) {
        expr += it;
    }
    sub_optimal.addConstr(expr, '>', LB - RC_FIXING);

    auto cont = false;
    do {
        sub_optimal.update();
        sub_optimal.optimize();

        auto aux_pi = beta |
                  ranges::views::transform([&](const auto& tmp) -> double {
                      return tmp.get(GRB_DoubleAttr_X);
                  }) |
                  ranges::to_vector;
        aux_pi.push_back(last.get(GRB_DoubleAttr_X));
        auto pi = std::span<const double>{aux_pi};

        auto sol = pricing_algorithm(pi);
        auto rc = sol.cost + pi.back();
        for (auto& it : sol.jobs) {
            rc -= pi[it->job];
        }

        if (rc < -RC_FIXING) {
            GRBLinExpr expr_pricing = -last;
            for (auto& it : sol.jobs) {
                expr_pricing += beta[it->job];
            }
            sub_optimal.addConstr(expr_pricing, '<', sol.cost);
            cont = true;
        } else {
            cont = false;
            calculate_constLB(pi);
            removed = evaluate_nodes(pi);
        }

    } while (cont);

    return removed;
}

void PricerSolverBase::remove_constraints(int first, int nb_del) {
    reformulation_model.delete_constraints(first, nb_del);
}

void PricerSolverBase::update_rows_coeff([[maybe_unused]] size_t first) {}

auto PricerSolverBase::print_num_paths() -> boost::multiprecision::cpp_int {
    return 0UL;
}

[[maybe_unused]] auto PricerSolverBase::get_UB() const -> double {
    return UB;
}

void PricerSolverBase::update_UB(double _ub) {
    if (_ub < UB) {
        UB = _ub;
    }
}

[[maybe_unused]] auto PricerSolverBase::get_int_attr_model(enum MIP_Attr c)
    -> int {
    int val = -1;
    switch (c) {
        case MIP_Attr_Nb_Vars:
            val = model.get(GRB_IntAttr_NumVars);
            break;
        case MIP_Attr_Nb_Constr:
            val = model.get(GRB_IntAttr_NumConstrs);
            break;
        case MIP_Attr_Status:
            val = model.get(GRB_IntAttr_Status);
            break;
        default:
            break;
    }

    return val;
}

auto PricerSolverBase::get_dbl_attr_model(enum MIP_Attr c) -> double {
    double val{};
    int    status = model.get(GRB_IntAttr_Status);
    if (status != GRB_INF_OR_UNBD && status != GRB_INFEASIBLE &&
        status != GRB_UNBOUNDED) {
        switch (c) {
            case MIP_Attr_Obj_Bound:
                val = model.get(GRB_DoubleAttr_ObjBound);
                break;
            case MIP_Attr_Obj_Bound_LP:
                val = model.get(GRB_DoubleAttr_ObjBoundC);
                break;
            case MIP_Attr_Mip_Gap:
                val = model.get(GRB_DoubleAttr_MIPGap);
                break;
            case MIP_Attr_Run_Time:
                val = model.get(GRB_DoubleAttr_Runtime);
                break;
            case MIP_Attr_Nb_Simplex_Iter:
                val = model.get(GRB_DoubleAttr_IterCount);
                break;
            case MIP_Attr_Nb_Nodes:
                val = model.get(GRB_DoubleAttr_NodeCount);
                break;
            default:
                val = std::numeric_limits<double>::max();
                break;
        }
    } else {
        switch (c) {
            case MIP_Attr_Run_Time:
                val = model.get(GRB_DoubleAttr_Runtime);
                break;
            case MIP_Attr_Nb_Simplex_Iter:
                val = model.get(GRB_DoubleAttr_IterCount);
                break;
            case MIP_Attr_Nb_Nodes:
                val = model.get(GRB_DoubleAttr_NodeCount);
                break;
            default:
                val = std::numeric_limits<double>::max();
                break;
        }
    }
    return val;
}


auto PricerSolverBase::compute_reduced_cost(const PricingSolution&   sol,
                                            std::span<const double>& pi,
                                            double* lhs) -> double {
    double result = sol.cost;
    // auto      nb_constraints = reformulation_model.get_nb_constraints();
    std::span aux_lhs{lhs, reformulation_model.size()};
    std::ranges::fill(aux_lhs, 0.0);

    for (auto it : sol.jobs | ranges::views::transform(
                                  [](const auto& tmp) { return tmp->job; })) {
        VariableKeyBase k(it, 0);
        for (const auto&& [constr, aux_pi, aux_lhs_] :
             ranges::views::zip(reformulation_model, pi, aux_lhs)) {
            if (constr == reformulation_model[convex_constr_id]) {
                continue;
            }
            auto coeff = (*constr)(k);

            if (std::fabs(coeff) > EPS_SOLVER) {
                result -= coeff * aux_pi;
                aux_lhs_ += coeff;
            }
        }
    }

    double          dual = pi[convex_constr_id];
    auto*           constr = reformulation_model[convex_constr_id].get();
    VariableKeyBase k(0, 0, true);
    double          coeff = (*constr)(k);
    result -= coeff * dual;
    aux_lhs[convex_constr_id] += coeff;

    return result;
}

void PricerSolverBase::compute_lhs(const PricingSolution& sol, double* lhs) {
    std::span aux_lhs{lhs, reformulation_model.size()};
    std::ranges::fill(aux_lhs, 0.0);

    for (auto it : sol.jobs | ranges::views::transform(
                                  [](const auto& tmp) { return tmp->job; })) {
        VariableKeyBase k(it, 0);
        for (const auto&& [c, constr] :
             reformulation_model | ranges::views::enumerate) {
            if (c == convex_constr_id) {
                continue;
            }
            auto coeff = (*constr)(k);

            if (fabs(coeff) > EPS_SOLVER) {
                aux_lhs[c] += coeff;
            }
        }
    }

    auto*           constr = reformulation_model[convex_constr_id].get();
    VariableKeyBase k(0, 0, true);
    double          coeff = (*constr)(k);
    aux_lhs[convex_constr_id] += coeff;
}

void PricerSolverBase::compute_lhs(const Column& sol, double* lhs) {
    std::span aux_lhs{lhs, reformulation_model.size()};
    std::ranges::fill(aux_lhs, 0.0);

    for (auto it : sol.job_list | ranges::views::transform([](const auto& tmp) {
                       return tmp->job;
                   })) {
        VariableKeyBase k(it, 0);
        for (const auto&& [c, constr] :
             reformulation_model | ranges::views::enumerate) {
            if (c == convex_constr_id) {
                continue;
            }
            auto coeff = (*constr)(k);

            if (fabs(coeff) > EPS_SOLVER) {
                aux_lhs[c] += coeff;
            }
        }
    }

    auto*           constr = reformulation_model[convex_constr_id].get();
    VariableKeyBase k(0, 0, true);
    double          coeff = (*constr)(k);
    aux_lhs[convex_constr_id] += coeff;
}

auto PricerSolverBase::compute_reduced_cost_simple(const PricingSolution&   sol,
                                                   std::span<const double>& pi)
    -> double {
    double result = sol.cost;

    for (const auto& it : sol.jobs) {
        VariableKeyBase k(it->job, 0);
        for (const auto&& [constr, pi_aux] :
             ranges::views::zip(reformulation_model, pi)) {
            if (constr == reformulation_model[convex_constr_id]) {
                continue;
            }
            auto coeff = (*constr)(k);

            if (std::fabs(coeff) > EPS_SOLVER) {
                result -= coeff * pi_aux;
            }
        }
    }

    double          dual = pi[convex_constr_id];
    auto*           constr = reformulation_model[convex_constr_id].get();
    VariableKeyBase k(0, 0, true);
    double          coeff = (*constr)(k);
    result -= coeff * dual;

    return result;
}

auto PricerSolverBase::compute_lagrange(const PricingSolution&     sol,
                                        const std::vector<double>& pi)
    -> double {
    double result = sol.cost;
    double dual_bound = 0.0;

    for (const auto& it : sol.jobs) {
        VariableKeyBase k(it->job, 0);
        auto            dual = pi[it->job];
        auto*           constr = reformulation_model[it->job].get();
        auto            coeff = (*constr)(k);

        if (std::fabs(coeff) > EPS_SOLVER) {
            result -= coeff * dual;
        }

        for (auto c = convex_constr_id + 1; c < reformulation_model.size();
             c++) {
            double dual_ = pi[c];
            double coeff_ = (*reformulation_model[c])(k);

            if (std::fabs(coeff_) > EPS_SOLVER) {
                result -= coeff_ * dual_;
            }
        }
    }

    result = std::min(0.0, result);

    for (const auto&& [c, constr] :
         reformulation_model | ranges::views::enumerate) {
        if (c == convex_constr_id) {
            continue;
        }

        dual_bound += constr->get_rhs() * pi[c];
    }

    result = -reformulation_model[convex_constr_id]->get_rhs() * result;
    result = dual_bound + result;

    return result;
}

auto PricerSolverBase::compute_lagrange(const PricingSolution&         sol,
                                        const std::span<const double>& pi)
    -> double {
    double result = sol.cost;
    double dual_bound = 0.0;

    for (const auto& it : sol.jobs) {
        VariableKeyBase k(it->job, 0);
        auto            dual = pi[it->job];
        auto*           constr = reformulation_model[it->job].get();
        auto            coeff = (*constr)(k);

        if (fabs(coeff) > EPS_SOLVER) {
            result -= coeff * dual;
        }

        for (auto c = convex_constr_id + 1; c < reformulation_model.size();
             c++) {
            double dual_ = pi[c];
            double coeff_ = (*reformulation_model[c])(k);

            if (fabs(coeff_) > EPS_SOLVER) {
                result -= coeff_ * dual_;
            }
        }
    }

    result = std::min(0.0, result);

    for (const auto&& [c, constr] :
         reformulation_model | ranges::views::enumerate) {
        if (c == convex_constr_id) {
            continue;
        }

        dual_bound += constr->get_rhs() * pi[c];
    }

    result = -reformulation_model[convex_constr_id]->get_rhs() * result;
    result = dual_bound + result;

    return result;
}

auto PricerSolverBase::compute_subgradient(const PricingSolution& sol,
                                           double* subgradient) -> double {
    std::span aux_subgradient{subgradient, reformulation_model.size()};
    auto      rhs = -reformulation_model[convex_constr_id]->get_rhs();

    for (size_t i = 0; const auto& constr : reformulation_model) {
        aux_subgradient[i] = constr->get_rhs();
        ++i;
    }

    for (const auto& it : sol.jobs) {
        VariableKeyBase k(it->job, 0);
        auto*           constr = reformulation_model[it->job].get();
        auto            coeff = (*constr)(k);

        if (fabs(coeff) > EPS_SOLVER) {
            aux_subgradient[k.get_j()] -= coeff * rhs;
        }

        for (auto c = convex_constr_id + 1; c < reformulation_model.size();
             c++) {
            auto coeff_ = (*reformulation_model[c])(k);

            if (fabs(coeff_) > EPS_SOLVER) {
                aux_subgradient[c] -= coeff_ * rhs;
            }
        }
    }

    aux_subgradient[convex_constr_id] += rhs;

    return 0.0;
}


void PricerSolverBase::calculate_constLB(std::span<const double>& pi) {
    constLB = 0.0;
    for (const auto&& [constr, aux_pi] :
         ranges::views::zip(reformulation_model, pi)) {
        if (constr == reformulation_model[convex_constr_id]) {
            continue;
        }
        constLB += constr->get_rhs() * aux_pi;
    }
}

const std::shared_ptr<GRBEnv> PricerSolverBase::genv =
    std::make_shared<GRBEnv>();

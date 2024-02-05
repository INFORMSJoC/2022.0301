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

#include "PricingStabilization.hpp"
#include <fmt/core.h>              // for print
#include <algorithm>               // for max, min
#include <cmath>                   // for sqrt
#include <memory>                  // for make_unique, unique_ptr
#include <range/v3/view/take.hpp>  // for take
#include <range/v3/view/zip.hpp>   // for zip
#include <utility>                 // for move
#include <vector>                  // for vector
#include "DebugLvl.hpp"            // for debug_lvl
#include "ModelInterface.hpp"      // for ReformulationModel
#include "PricerSolverBase.hpp"    // for PricerSolverBase
#include "or-utils/util.h"         // for dbg_lvl, SQR

/**
 * @brief Construct a new Pricing Stabilization Base:: Pricing Stabilization
 * Base object
 *
 * @param _solver
 */
PricingStabilizationBase::PricingStabilizationBase(
    PricerSolverBase*        _solver,
    std::span<const double>& _pi_out)
    : solver(_solver),
      pi_out(_pi_out),
      pi_in(_solver->convex_constr_id + 1, 0.0),
      pi_sep(_solver->convex_constr_id + 1, 0.0) {}

/**
@brief solve the pricing problem without stabilization
 *
 * @param _eta_out
 * @param _pi_out
 * @param _lhs
 */
void PricingStabilizationBase::solve(double _eta_out,
                                     //  const std::vector<double>& _pi_out,
                                     double* _lhs) {
    solver->calculate_constLB(pi_out);
    sol = solver->pricing_algorithm(pi_out);
    reduced_cost = solver->compute_reduced_cost(sol, pi_out, _lhs);
    eta_sep = solver->compute_lagrange(sol, pi_out);
    update_stab_center = false;
    eta_out = _eta_out;
    continueLP = (eta_in < eta_out);
    if (eta_sep > eta_in) {
        eta_in = eta_sep;
        pi_in.assign(pi_out.begin(), pi_out.end());
        update_stab_center = true;
    }

    if (debug_lvl(0)) {
        fmt::print(
            R"(alpha = {1}, result of primal bound and Lagragian bound: out = {2:.{0}f}, in = {3:.{0}f} UB = {4}
)",
            4, solver->get_nb_vertices(), eta_out, eta_in, solver->UB);
    }
    iterations++;
}

auto PricingStabilizationBase::get_sol() -> PricingSolution& {
    return sol;
}

auto PricingStabilizationBase::get_reduced_cost() const -> double {
    return reduced_cost;
}

auto PricingStabilizationBase::get_eps_stab_solver() -> double {
    return ETA_DIFF;
}

auto PricingStabilizationBase::get_update_stab_center() const -> bool {
    return update_stab_center;
}

auto PricingStabilizationBase::get_eta_in() -> double {
    return eta_in;
}

auto PricingStabilizationBase::get_eta_sep() -> double {
    return eta_in;
}

auto PricingStabilizationBase::stopping_criteria() -> bool {
    return eta_out - eta_in > ETA_DIFF;
}

void PricingStabilizationBase::update_duals() {
    if (solver->reformulation_model.size() != pi_sep.size()) {
        pi_sep.resize(solver->reformulation_model.size(), 0.0);
        pi_in.resize(solver->reformulation_model.size(), 0.0);
    }
}

auto PricingStabilizationBase::reduced_cost_fixing() -> bool {
    previous_fix = iterations;
    auto aux_sep = std::span<const double>{pi_sep};
    solver->calculate_constLB(aux_sep);
    return solver->evaluate_nodes(aux_sep);
}

void PricingStabilizationBase::remove_constraints(int first, int nb_del) {
    auto it_sep = pi_sep.begin() + first;
    pi_sep.erase(it_sep, it_sep + nb_del);
    auto it_in = pi_in.begin() + first;
    pi_in.erase(it_in, it_in + nb_del);
}

auto PricingStabilizationBase::do_reduced_cost_fixing() -> bool {
    return ((iterations == 1 || !continueLP ||
             previous_fix <= iterations - RC_FIXING_RATE));
}

void PricingStabilizationBase::update_continueLP(double obj) {
    eta_out = obj;
    continueLP = (ETA_DIFF < eta_out - eta_in);
}

auto PricingStabilizationBase::clone(PricerSolverBase*        _solver,
                                     std::span<const double>& _pi_out)
    -> std::unique_ptr<PricingStabilizationBase> {
    return std::make_unique<PricingStabilizationBase>(_solver, _pi_out);
}
/**
 * @brief Wentgnes stabilization technique
 *
 */

/**
 * @brief Construct a new Pricing Stabilization Stat:: Pricing Stabilization
 * Stat object
 *
 * @param _solver
 */
PricingStabilizationStat::PricingStabilizationStat(
    PricerSolverBase*        _solver,
    std::span<const double>& _pi_out)
    : PricingStabilizationBase(_solver, _pi_out) {}

/**
@brief Solve pricing problem with Wentgnes stabilization
 *
 * @param _eta_out
 * @param _pi_out
 * @param _lhs_coeff
 */
void PricingStabilizationStat::solve(double _eta_out, double* _lhs_coeff) {
    k = 0;
    bool miss_pricing = true;
    update = 0;

    eta_out = _eta_out;
    iterations++;
    update_stab_center = false;

    auto aux_sep = std::span<const double>{ pi_sep };
    do {
        k += 1;
        alphabar =
            (hasstabcenter) ? std::max(0.0, 1.0 - k * (1.0 - alpha)) : 0.0;
        compute_pi_eta_sep(alphabar);
        solver->calculate_constLB(aux_sep);
        PricingSolution aux_sol = solver->pricing_algorithm(aux_sep);

        eta_sep = solver->compute_lagrange(aux_sol, aux_sep);
        reduced_cost =
            solver->compute_reduced_cost(aux_sol, pi_out, _lhs_coeff);

        if (reduced_cost < EPS_RC) {
            sol = std::move(aux_sol);
            update = 1;
            miss_pricing = false;
        }
    } while (miss_pricing && alphabar > 0); /** miss_pricing check */

    continueLP = (ETA_DIFF < eta_out - eta_in);

    if (eta_sep > eta_in) {
        hasstabcenter = true;
        eta_in = eta_sep;
        pi_in = pi_sep;
        update_stab_center = true;
    }

    if (iterations % (solver->convex_constr_id) == 0 && debug_lvl(0)) {
        fmt::print(
            R"(alpha = {1:.{0}f}, result of primal bound and Lagragian bound: out = {2:.{0}f}, in = {3:.{0}f}
)",
            4, alpha, eta_out, eta_in);
    }
}

auto PricingStabilizationStat::get_eta_in() -> double {
    return eta_in;
}

auto PricingStabilizationStat::stopping_criteria() -> bool {
    return eta_in <= (eta_out - ETA_DIFF);
}

void PricingStabilizationStat::update_duals() {
    if (solver->reformulation_model.size() != pi_sep.size()) {
        pi_sep.resize(solver->reformulation_model.size(), 0.0);
        // pi_out.resize(solver->reformulation_model.size(), 0.0);
        pi_in.resize(solver->reformulation_model.size(), 0.0);
    }
}

void PricingStabilizationStat::remove_constraints(int first, int nb_del) {
    auto it_in = pi_in.begin() + first;
    pi_in.erase(it_in, it_in + nb_del);
    // auto it_out = pi_out.begin() + first;
    // pi_out.erase(it_out, it_out + nb_del);
    auto it_sep = pi_sep.begin() + first;
    pi_sep.erase(it_sep, it_sep + nb_del);
}

void PricingStabilizationStat::set_alpha(double _alpha) {
    alpha = _alpha;
}

void PricingStabilizationStat::compute_pi_eta_sep(double _alpha_bar) {
    double beta_bar = 1.0 - _alpha_bar;

    for (auto&& [sep, in, out] :
         ranges::views::zip(pi_sep, pi_in, pi_out) |
             ranges::views::take(solver->reformulation_model.size())) {
        sep = _alpha_bar * in + beta_bar * out;
    }

    eta_sep = _alpha_bar * (eta_in) + beta_bar * (eta_out);
}

auto PricingStabilizationStat::clone(PricerSolverBase*        _solver,
                                     std::span<const double>& _pi_out)
    -> std::unique_ptr<PricingStabilizationBase> {
    return std::make_unique<PricingStabilizationStat>(_solver, _pi_out);
}

/**
 * Dynamic stabilization technique
 */

PricingStabilizationDynamic::PricingStabilizationDynamic(
    PricerSolverBase*        _solver,
    std::span<const double>& _pi_out)
    : PricingStabilizationStat(_solver, _pi_out),
      subgradient(_solver->convex_constr_id + 1, 0.0)

{}

void PricingStabilizationDynamic::solve(double _eta_out, double* _lhs) {
    k = 0;
    double result_sep{};
    bool   miss_pricing = true;
    update = 0;

    eta_out = _eta_out;
    iterations++;
    update_stab_center = false;

    auto aux_sep = std::span<const double>(pi_sep.data(), pi_sep.size());
    do {
        k += 1;
        alphabar = hasstabcenter ? std::max(0.0, 1.0 - k * (1 - alpha)) : 0.0;
        compute_pi_eta_sep(alphabar);
        auto aux_sol = solver->pricing_algorithm(aux_sep);
        result_sep = solver->compute_lagrange(aux_sol, aux_sep);
        reduced_cost = solver->compute_reduced_cost(aux_sol, pi_out, _lhs);

        if (reduced_cost < EPS_RC) {
            solver->compute_subgradient(aux_sol, subgradient.data());
            adjust_alpha();
            sol = std::move(aux_sol);
            alpha = alphabar;
            update = 1;
            miss_pricing = false;
        }
    } while (miss_pricing && alphabar > 0.0);

    if (result_sep > eta_in) {
        hasstabcenter = true;
        eta_in = result_sep;
        pi_in = pi_sep;
        update_stab_center = true;
        if (iterations % RC_FIXING_RATE == 0) {
            solver->calculate_constLB(aux_sep);
            solver->evaluate_nodes(aux_sep);
        }
    }

    if (iterations % solver->convex_constr_id == 0 && debug_lvl(0)) {
        fmt::print(
            R"(alpha = {1:.{0}f}, result of primal bound and Lagragian bound: out = {2:.{0}f}, in = {3:.{0}f}
)",
            4, alpha, eta_out, eta_in);
    }
}
void PricingStabilizationDynamic::compute_subgradient(
    const PricingSolution& _sol) {
    solver->compute_subgradient(_sol, subgradient.data());

    // subgradientnorm = 0.0;

    // for (int i = 0; i < nb_rows; ++i) {
    //     double sqr = SQR(subgradient_in[i]);

    //     if (sqr > 0.00001) {
    //         subgradientnorm += sqr;
    //     }
    // }

    // subgradientnorm = SQRT(subgradientnorm);
}

void PricingStabilizationDynamic::adjust_alpha() {
    auto sum = 0.0;

    for (auto i = 0UL; i < solver->reformulation_model.size(); ++i) {
        sum += subgradient[i] * (pi_out[i] - pi_in[i]);
    }

    if (sum > 0) {
        alphabar = std::max<double>(0, alphabar - ALPHA_CHG);
    } else {
        alphabar =
            std::min<double>(ALPHA_MAX, alphabar + (1 - alphabar) * ALPHA_CHG);
    }
}

auto PricingStabilizationDynamic::clone(PricerSolverBase*        _solver,
                                        std::span<const double>& _pi_out)
    -> std::unique_ptr<PricingStabilizationBase> {
    return std::make_unique<PricingStabilizationDynamic>(_solver, _pi_out);
}

void PricingStabilizationDynamic::update_duals() {
    if (solver->reformulation_model.size() != pi_sep.size()) {
        pi_sep.resize(solver->reformulation_model.size(), 0.0);
        // pi_out.resize(solver->reformulation_model.size(), 0.0);
        pi_in.resize(solver->reformulation_model.size(), 0.0);
        subgradient.resize(solver->reformulation_model.size(), 0.0);
    }
}

void PricingStabilizationDynamic::remove_constraints(int first, int nb_del) {
    auto it_in = pi_in.begin() + first;
    pi_in.erase(it_in, it_in + nb_del);
    // auto it_out = pi_out.begin() + first;
    // pi_out.erase(it_out, it_out + nb_del);
    auto it_sep = pi_sep.begin() + first;
    pi_sep.erase(it_sep, it_sep + nb_del);
    auto it_sub = subgradient.begin() + first;
    subgradient.erase(it_sub, it_sub + nb_del);
}

/**
@brief Construct a new Pricing Stabilization Hybrid:: Pricing Stabilization
Hybrid object
 *
 * @param pricer_solver
 */
PricingStabilizationHybrid::PricingStabilizationHybrid(
    PricerSolverBase*        pricer_solver,
    std::span<const double>& _pi_out)
    : PricingStabilizationDynamic(pricer_solver, _pi_out),
      subgradient_in(pricer_solver->convex_constr_id + 1) {}

void PricingStabilizationHybrid::update_duals() {
    if (solver->reformulation_model.size() != pi_sep.size()) {
        pi_sep.resize(solver->reformulation_model.size(), 0.0);
        // pi_out.resize(solver->reformulation_model.size(), 0.0);
        pi_in.resize(solver->reformulation_model.size(), 0.0);
        subgradient_in.resize(solver->reformulation_model.size(), 0.0);
        subgradient.resize(solver->reformulation_model.size(), 0.0);
    }
}

auto PricingStabilizationHybrid::clone(PricerSolverBase*        _solver,
                                       std::span<const double>& _pi_out)
    -> std::unique_ptr<PricingStabilizationBase> {
    return std::make_unique<PricingStabilizationBase>(_solver, _pi_out);
}

void PricingStabilizationHybrid::remove_constraints(int first, int nb_del) {
    auto it_in = pi_in.begin() + first;
    pi_in.erase(it_in, it_in + nb_del);
    // auto it_out = pi_out.begin() + first;
    // pi_out.erase(it_out, it_out + nb_del);
    auto it_sep = pi_sep.begin() + first;
    pi_sep.erase(it_sep, it_sep + nb_del);
    auto it_sub = subgradient.begin() + first;
    subgradient.erase(it_sub, it_sub + nb_del);
    auto it_sub_in = subgradient_in.begin() + first;
    subgradient_in.erase(it_sub_in, it_sub_in + nb_del);
}

void PricingStabilizationHybrid::solve(double _eta_out,
                                       //    const std::vector<double>& _pi_out,
                                       double* _lhs) {
    update = 0;
    bool stabilized = false;
    eta_out = _eta_out;
    iterations++;

    do {
        update_hybrid();

        stabilized = is_stabilized();

        for (auto i = 0UL; i < solver->reformulation_model.size(); ++i) {
            pi_sep[i] = compute_dual(i);
        }

        auto            aux_sep = std::span<const double>{pi_sep};
        PricingSolution aux_sol = solver->pricing_algorithm(aux_sep);

        eta_sep = solver->compute_lagrange(aux_sol, aux_sep);
        reduced_cost = solver->compute_reduced_cost(aux_sol, pi_out, _lhs);

        update_stabcenter(aux_sol);

        if (reduced_cost < EPS_RC) {
            if (in_mispricing_schedule) {
                in_mispricing_schedule = false;
            }
            sol = std::move(aux_sol);
            update_subgradientproduct();
            update_alpha();
        } else {
            if (stabilized) {
                in_mispricing_schedule = true;
                update_alpha_misprice();
            } else {
                in_mispricing_schedule = false;
            }
        }
    } while (in_mispricing_schedule && stabilized);

    if (iterations % solver->convex_constr_id == 0 && debug_lvl(0)) {
        fmt::print(
            R"(alpha = {1:.{0}f}, result of primal bound and Lagragian bound: out = {2:.{0}f}, in = {3:.{0}f}
)",
            4, alpha, eta_out, eta_in);
    }
}
void PricingStabilizationHybrid::calculate_dualdiffnorm() {
    dualdiffnorm = 0.0;

    for (auto i = 0UL; i < solver->reformulation_model.size(); ++i) {
        double dualdiff = SQR(pi_in[i] - pi_out[i]);
        if (dualdiff > EPS_STAB) {
            dualdiffnorm += dualdiff;
        }
    }

    dualdiffnorm = std::sqrt(dualdiffnorm);
}

void PricingStabilizationHybrid::calculate_beta() {
    beta = 0.0;
    for (auto i = 0UL; i < solver->convex_constr_id; ++i) {
        double dualdiff = std::abs(pi_out[i] - pi_in[i]);
        double product = dualdiff * std::abs(subgradient_in[i]);

        if (product > EPS_STAB) {
            beta += product;
        }
    }

    if (subgradientnorm > EPS_STAB) {
        beta = beta / (subgradientnorm * dualdiffnorm);
    }
}

void PricingStabilizationHybrid::calculate_hybridfactor() {
    double aux_norm = 0.0;
    for (auto i = 0UL; i < solver->convex_constr_id; ++i) {
        double aux_double =
            SQR((beta - 1.0) * (pi_out[i] - pi_in[i]) +
                beta * (subgradient_in[i] * dualdiffnorm / subgradientnorm));
        if (aux_double > EPS_STAB) {
            aux_norm += aux_double;
        }
    }
    aux_norm = std::sqrt(aux_norm);

    hybridfactor = ((1 - alpha) * dualdiffnorm) / aux_norm;
}

auto PricingStabilizationHybrid::is_stabilized() -> bool {
    if (in_mispricing_schedule) {
        return alphabar > 0.0;
    }
    return alpha > 0.0;
}

auto PricingStabilizationHybrid::compute_dual(auto i) -> double {
    double usedalpha = alpha;
    double usedbeta = beta;

    if (in_mispricing_schedule) {
        usedalpha = alphabar;
        usedbeta = 0.0;
    }

    if (hasstabcenter && (usedbeta == 0.0 || usedalpha == 0.0)) {
        return usedalpha * pi_in[i] + (1.0 - usedalpha) * pi_out[i];
    }

    if (hasstabcenter && usedbeta > 0.0) {
        double dual = pi_in[i] +
                      hybridfactor *
                          (beta * (pi_in[i] + subgradient_in[i] * dualdiffnorm /
                                                  subgradientnorm) +
                           (1.0 - beta) * pi_out[i] - pi_in[i]);
        return std::max(dual, 0.0);
    }

    return pi_out[i];
}

void PricingStabilizationHybrid::update_stabcenter(
    const PricingSolution& _sol) {
    if (eta_sep > eta_in) {
        pi_in = pi_sep;
        compute_subgradient_norm(_sol);
        eta_in = eta_sep;
        hasstabcenter = true;
        update_stab_center = true;
        auto aux_sep = std::span<const double>{pi_sep};
        solver->calculate_constLB(aux_sep);
        solver->evaluate_nodes(aux_sep);
    }
}

void PricingStabilizationHybrid::compute_subgradient_norm(
    const PricingSolution& _sol) {
    // double* subgradient_in = &g_array_index(pd->subgradient_in, double,
    // 0); double* rhs = &g_array_index(pd->rhs, double, 0);
    // fill_dbl(subgradient_in, pd->nb_rows, 1.0);
    // subgradient_in[pd->nb_jobs] = 0.0;
    solver->compute_subgradient(_sol, subgradient_in.data());

    // for (guint i = 0; i < sol.jobs->len; i++) {
    //     Job* tmp_j = reinterpret_cast<Job*>(g_ptr_array_index(sol.jobs,
    //     i)); subgradient_in[tmp_j->job] += rhs[pd->nb_jobs] * 1.0;
    // }

    subgradientnorm = 0.0;

    for (auto i = 0UL; i < solver->reformulation_model.size(); ++i) {
        double sqr = SQR(subgradient_in[i]);

        if (sqr > EPS_STAB) {
            subgradientnorm += sqr;
        }
    }

    subgradientnorm = std::sqrt(subgradientnorm);
}
void PricingStabilizationHybrid::update_subgradientproduct() {
    subgradientproduct = 0.0;
    for (auto i = 0UL; i < solver->reformulation_model.size(); ++i) {
        subgradientproduct += (pi_out[i] - pi_in[i]) * subgradient_in[i];
    }
}

void PricingStabilizationHybrid::update_alpha_misprice() {
    k++;
    alphabar = std::max(0.0, 1 - k * (1.0 - alpha));
}

void PricingStabilizationHybrid::update_alpha() {
    if (subgradientproduct > 0.0) {
        alpha = std::max(0.0, alpha - ALPHA_CHG);
    } else {
        alpha = std::min(ALPHA_MAX, alpha + (1.0 - alpha) * ALPHA_CHG);
    }
}

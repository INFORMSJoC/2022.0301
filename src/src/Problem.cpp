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

#include "Problem.h"
#include <fmt/core.h>                   // for print
#include <boost/timer/timer.hpp>        // for auto_cpu_timer
#include <cstddef>                      // for size_t
#include <memory>                       // for make_unique, unique_ptr, shar...
#include <utility>                      // for move
#include "BranchBoundTree.hpp"          // for BranchBoundTree
#include "Instance.h"                   // for Instance
#include "LocalSearch.hpp"              // for LocalSearchData, PerturbOperator
#include "NodeData.h"                   // for NodeData
#include "PricerSolverArcTimeDP.hpp"    // for PricerSolverArcTimeDp
#include "PricerSolverBase.hpp"         // for PricerSolverBase
#include "PricerSolverBddBackward.hpp"  // for PricerSolverBddBackwardCycle
#include "PricerSolverBddForward.hpp"   // for PricerSolverBddCycle, PricerS...
#include "PricerSolverSimpleDP.hpp"     // for PricerSolverSimpleDp
#include "PricingStabilization.hpp"     // for PricingStabilizationBase, Pri...
#include "Solution.hpp"                 // for Sol
#include "Statistics.h"                 // for Statistics, Statistics::bb_timer

Problem::~Problem() = default;

Problem::Problem(int argc, const char** argv)
    : parms(argc, argv),
      stat(parms),
      instance(parms),
      root_pd(std::make_unique<NodeData>(this)),
      status(no_sol) {
    /**
     * @brief
     * Finding heuristic solutions to the problem or start without feasible
     * solutions
     */
    stat.start_resume_timer(Statistics::heuristic_timer);
    if (parms.use_heuristic.value()) {
        heuristic();
    } else {
        Sol best_sol(instance);
        best_sol.construct_edd(instance.jobs);
        fmt::print("Solution Constructed with EDD heuristic:\n");
        best_sol.print_solution();
        best_sol.canonical_order(instance.intervals);
        fmt::print("Solution in canonical order: \n");
        best_sol.print_solution();
        opt_sol = best_sol;
    }
    stat.suspend_timer(Statistics::heuristic_timer);

    /**
     * @brief
     * Create pricing solver
     */
    stat.start_resume_timer(Statistics::build_dd_timer);
    switch (parms.pricing_solver.value()) {
        case bdd_solver_simple:
            root_pd->solver = std::make_unique<PricerSolverBddSimple>(instance);
            break;
        case bdd_solver_cycle:
            root_pd->solver = std::make_unique<PricerSolverBddCycle>(instance);
            break;
        case bdd_solver_backward_simple:
            root_pd->solver =
                std::make_unique<PricerSolverBddBackwardSimple>(instance);
            break;
        case bdd_solver_backward_cycle:
            root_pd->solver =
                std::make_unique<PricerSolverBddBackwardCycle>(instance);
            break;
        case dp_solver:
            root_pd->solver = std::make_unique<PricerSolverSimpleDp>(instance);
            break;
        case ati_solver:
            root_pd->solver = std::make_unique<PricerSolverArcTimeDp>(instance);
            break;
        case dp_bdd_solver:
            root_pd->solver = std::make_unique<PricerSolverSimpleDp>(instance);
            break;
        default:
            root_pd->solver =
                std::make_unique<PricerSolverBddBackwardCycle>(instance);
            break;
    }
    stat.suspend_timer(Statistics::build_dd_timer);
    stat.first_size_graph = root_pd->solver->get_nb_edges();

    /**
     * @brief
     * Initial stabilization method
     */
    auto* tmp_solver = root_pd->solver.get();
    switch (parms.stab_technique.value()) {
        case stab_wentgnes:
            root_pd->solver_stab = std::make_unique<PricingStabilizationStat>(
                tmp_solver, root_pd->pi);
            break;
        case stab_dynamic:
            root_pd->solver_stab =
                std::make_unique<PricingStabilizationDynamic>(tmp_solver,
                                                              root_pd->pi);
            break;
        case stab_hybrid:
            root_pd->solver_stab = std::make_unique<PricingStabilizationHybrid>(
                tmp_solver, root_pd->pi);
            break;
        case no_stab:
            root_pd->solver_stab = std::make_unique<PricingStabilizationBase>(
                tmp_solver, root_pd->pi);
            break;
        default:
            root_pd->solver_stab = std::make_unique<PricingStabilizationStat>(
                tmp_solver, root_pd->pi);
            break;
    }

    if (instance.best_known_sol.has_value()) {
        opt_sol.tw = *instance.best_known_sol - instance.off;
        root_pd->solver->update_UB(*instance.best_known_sol - instance.off);
        root_pd->upper_bound = *instance.best_known_sol - instance.off;
        stat.root_upper_bound = *instance.best_known_sol;
        stat.global_upper_bound = *instance.best_known_sol;

    }
    else
    {
    root_pd->solver->update_UB(opt_sol.tw);
    root_pd->solver_stab->set_alpha(parms.alpha.value());
    root_pd->upper_bound = opt_sol.tw;
    stat.root_upper_bound = opt_sol.tw + instance.off;
    stat.global_upper_bound = opt_sol.tw + instance.off;
    }

    /**
     * @brief
     * Initialization of the B&B tree
     */
    if (!parms.use_mip_solver.value()) {
        stat.start_resume_timer(Statistics::bb_timer);
        tree = std::make_unique<BranchBoundTree>(std::move(root_pd), 0, 1);
        tree->explore();
        stat.global_upper_bound =
            static_cast<int>(tree->get_UB()) + instance.off;
        stat.global_lower_bound =
            static_cast<int>(tree->get_LB()) + instance.off;
        stat.rel_error = (stat.global_upper_bound - stat.global_lower_bound) /
                         (EPS + stat.global_lower_bound);
        stat.suspend_timer(Statistics::bb_timer);
    } else {
        root_pd->build_rmp();
        root_pd->solve_relaxation();
        root_pd->stat.start_resume_timer(Statistics::lb_root_timer);
        root_pd->compute_lower_bound();
        stat.suspend_timer(Statistics::lb_root_timer);
        tmp_solver->project_sol_on_original_variables(opt_sol);
        tmp_solver->build_mip();
        // set_lb(pd->lower_bound);
        // set_obj_value(pd->LP_lower_bound);
    }
    if (parms.print_csv.value()) {
        to_csv();
    }
}

void Problem::solve() {
    tree->explore();
    if (parms.print_csv.value()) {
        to_csv();
    }
}

void Problem::heuristic() {
    // auto                         ILS = instance.nb_jobs / 2;
    auto IR = static_cast<size_t>(parms.nb_iterations_rvnd.value());
    boost::timer::auto_cpu_timer timer_heuristic;

    Sol best_sol{instance};
    best_sol.construct_edd(instance.jobs);
    fmt::print("Solution Constructed with EDD heuristic:\n");
    best_sol.print_solution();
    best_sol.canonical_order(instance.intervals);
    fmt::print("Solution in canonical order: \n");
    best_sol.print_solution();

    /** Local Search */
    auto local = LocalSearchData(instance.nb_jobs, instance.nb_machines);
    local.RVND(best_sol);
    /** Perturbation operator */
    PerturbOperator perturb{};

    best_sol.canonical_order(instance.intervals);
    fmt::print("Solution after local search:\n");
    best_sol.print_solution();
    root_pd->add_solution_to_colpool(best_sol);

    Sol sol{instance};
    sol = best_sol;
    for (auto i = 0UL; i < IR && best_sol.tw != 0; ++i) {
        Sol sol1{instance};
        sol1 = sol;
        perturb(sol1);
        local.RVND(sol1);

        if (sol1.tw < sol.tw) {
            sol = sol1;
            if (sol.tw < best_sol.tw) {
                best_sol = sol;
                root_pd->add_solution_to_colpool(best_sol);
                i /= 2;
            }
        }
    }

    fmt::print("Best new heuristics\n");
    best_sol.print_solution();
    opt_sol = best_sol;
}

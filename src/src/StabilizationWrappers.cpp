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

#include <memory>                    // for shared_ptr, unique_ptr, make_shared
#include <utility>                   // for move
#include <vector>                    // for vector
#include "Column.h"                  // for ScheduleSet
#include "NodeData.h"                // for NodeData
#include "Parms.h"                   // for Parms, yes_reduced_cost
#include "PricerSolverBase.hpp"      // for PricerSolverBase
#include "PricingSolution.hpp"       // for PricingSolution
#include "PricingStabilization.hpp"  // for PricingStabilizationBase
#include "Statistics.h"              // for Statistics, Statistics::reduced_...

auto NodeData::solve_pricing() -> int {
    int val = 0;

    solver_stab->solve(LP_lower_bound_BB, lhs_coeff.data());

    if (solver_stab->get_update_stab_center()) {
        if (solver_stab->do_reduced_cost_fixing() &&
            parms.reduce_cost_fixing.value()) {
            stat.start_resume_timer(Statistics::reduced_cost_fixing_timer);
            if (solver_stab->reduced_cost_fixing()) {
                delete_infeasible_columns();
            }
            stat.suspend_timer(Statistics::reduced_cost_fixing_timer);
        }
    }

    if (!solver_stab->continueLP) {
        stat.start_resume_timer(Statistics::reduced_cost_fixing_timer);
        if (solver_stab->reduced_cost_fixing()) {
            delete_infeasible_columns();
        }
        stat.suspend_timer(Statistics::reduced_cost_fixing_timer);
    }

    if (solver_stab->get_reduced_cost() < -EPS_BOUND &&
        solver_stab->continueLP &&
        (solver_stab->get_eta_in() < upper_bound - 1.0 + EPS_BOUND)) {
        localColPool.emplace_back(
            std::make_shared<Column>(std::move(solver_stab->get_sol())));
        // val = add_lhs_column_to_rmp(
        //     localColPool.back().get()->total_weighted_completion_time);

        val = add_lhs_column_to_rmp(
            localColPool.back()->total_weighted_completion_time, lhs_coeff);
    } else {
        stat.start_resume_timer(Statistics::reduced_cost_fixing_timer);
        if (solver_stab->reduced_cost_fixing()) {
            delete_infeasible_columns();
        }
        stat.suspend_timer(Statistics::reduced_cost_fixing_timer);
    }
    solve_relaxation();
    solver_stab->update_continueLP(LP_lower_bound);

    return val;
}

void NodeData::solve_farkas_dbl() {
    auto s = solver->farkas_pricing(pi);

    if (s.obj < EPS) {
        localColPool.emplace_back(std::make_shared<Column>(std::move(s)));
        add_lhs_column_to_rmp(
            localColPool.back()->total_weighted_completion_time);
    }

    solve_relaxation();
    solver_stab->update_continueLP(LP_lower_bound);
}

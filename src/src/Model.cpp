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

#include <OsiGrbSolverInterface.hpp>
#include <array>
#include <range/v3/algorithm/for_each.hpp>  // for for_each, for_...
#include <range/v3/numeric/iota.hpp>        // for iota, iota_fn
#include <range/v3/view/enumerate.hpp>      // for enumerate_fn
#include <range/v3/view/take.hpp>           // for take_view, take
#include <span>                             // for span
#include <vector>                           // for vector
#include "Column.h"                         // for ScheduleSet
#include "NodeData.h"                       // for NodeData
#include "PricerSolverBase.hpp"             // for PricerSolverBase
#include "gurobi_c.h"                       // for GRB_INFINITY

auto NodeData::add_lhs_column_to_rmp(double cost) -> int {
    id_row.clear();
    coeff_row.clear();

    // for (auto const&& [j, it] : lhs_coeff | ranges::views::enumerate) {
    //     if (std::abs(it) > EPS_BOUND) {
    //         id_row.emplace_back(static_cast<int>(j));
    //         coeff_row.emplace_back(it);
    //     }
    // }

    for (int j = 0; const auto& it : lhs_coeff) {
        if (std::abs(it) > EPS_BOUND) {
            id_row.emplace_back(j);
            coeff_row.emplace_back(it);
        }
        j++;
    }

    osi_rmp->addCol(static_cast<int>(id_row.size()), id_row.data(),
                    coeff_row.data(), 0.0, osi_rmp->getInfinity(), cost);

    return 0;
}

auto NodeData::add_lhs_column_to_rmp(double                     cost,
                                     const std::vector<double>& _lhs) -> int {
    id_row.clear();
    coeff_row.clear();

    for (auto const&& [j, it] : _lhs | ranges::views::enumerate) {
        if (std::abs(it) > EPS_BOUND) {
            id_row.emplace_back(static_cast<int>(j));
            coeff_row.emplace_back(it);
        }
    }

    osi_rmp->addCol(static_cast<int>(id_row.size()), id_row.data(),
                    coeff_row.data(), 0.0, osi_rmp->getInfinity(), cost);

    return 0;
}

void NodeData::create_assignment_constraints() {
    id_assignment_constraint = osi_rmp->getNumRows();
    std::vector<int>    start(nb_jobs + 1, 0);
    std::vector<double> rhs_tmp(nb_jobs, 1.0);
    std::vector<double> ub_rhs(nb_jobs, osi_rmp->getInfinity());

    osi_rmp->addRows(static_cast<int>(nb_jobs), start.data(), nullptr, nullptr,
                     rhs_tmp.data(), ub_rhs.data());
}

void NodeData::create_convex_contraint() {
    id_convex_constraint = osi_rmp->getNumRows();
    auto               rhs_lb = -static_cast<double>(nb_machines);
    auto               rhs_ub = osi_rmp->getInfinity();
    std::array<int, 2> start = {0, 0};

    osi_rmp->addRows(1, start.data(), nullptr, nullptr, &rhs_lb, &rhs_ub);

    id_valid_cuts = osi_rmp->getNumRows();
}

void NodeData::create_artificial_cols() {
    id_art_var_assignment = osi_rmp->getNumCols();
    id_art_var_convex = nb_jobs;
    id_art_var_cuts = static_cast<int>(nb_jobs) + 1;
    id_next_var_cuts = id_art_var_cuts;
    auto nb_vars = nb_jobs + 1 + max_nb_cuts;
    id_pseudo_schedules = static_cast<int>(nb_vars);

    std::vector<double> lb(nb_vars, 0.0);
    std::vector<double> ub(nb_vars, GRB_INFINITY);
    std::vector<double> obj(nb_vars, 100.0 * (opt_sol.tw + 1));
    std::vector<int> start_vars(nb_vars + 1, static_cast<int>(nb_jobs + 1UL));

    auto nz = static_cast<int>(nb_jobs + 1UL);

    std::vector<int> rows_ind(nz);
    ranges::iota(rows_ind, 0);
    std::vector<double> coeff_vals(nz, 1.0);
    coeff_vals[nb_jobs] = -1.0;
    ranges::iota(start_vars | ranges::views::take(nz), 0);

    osi_rmp->addCols(static_cast<int>(nb_vars), start_vars.data(),
                     rows_ind.data(), coeff_vals.data(), lb.data(), ub.data(),
                     obj.data());

    id_pseudo_schedules = osi_rmp->getNumCols();
}

void NodeData::add_cols_local_pool() {
    std::vector<double> _lhs(osi_rmp->getNumRows());
    ranges::for_each(localColPool, [&](auto& it) {
        solver->compute_lhs(*it.get(), _lhs.data());
        add_lhs_column_to_rmp(it->total_weighted_completion_time, _lhs);
    });
}

auto NodeData::build_rmp() -> int {
    /**
     * Set up messegahandler osi solver
     */
    osi_rmp->messageHandler()->setLogLevel(0);
    osi_rmp->messageHandler()->setPrefix(false);
    osi_rmp->setHintParam(OsiDoReducePrint, true, OsiHintTry);

    create_assignment_constraints();
    create_convex_contraint();
    create_artificial_cols();

    prune_duplicated_sets();
    add_cols_local_pool();
    osi_rmp->initialSolve();

    rhs = std::span<const double>(osi_rmp->getRightHandSide(),
                                  static_cast<size_t>(osi_rmp->getNumRows()));
    lhs_coeff.resize(nb_jobs + 1, 0.0);
    id_row.reserve(nb_jobs + 1);
    coeff_row.reserve(nb_jobs + 1);

    return 0;
}

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

#include <fmt/core.h>
#include <algorithm>                           // for min, __for_eac...
#include <cassert>                             // for assert
#include <memory>                              // for unique_ptr
#include <range/v3/action/action.hpp>          // for operator|=
#include <range/v3/action/remove_if.hpp>       // for remove_if, rem...
#include <range/v3/algorithm/sort.hpp>         // for sort, sort_fn
#include <range/v3/numeric/inner_product.hpp>  // for inner_product
#include <range/v3/range/conversion.hpp>       // for to_container::fn
#include <range/v3/view/drop.hpp>              // for drop, drop_fn
#include <range/v3/view/filter.hpp>            // for filter_view
#include <range/v3/view/transform.hpp>         // for transform_view
#include <range/v3/view/zip.hpp>               // for zip_view, zip
#include <utility>                             // for move, pair
#include <vector>                              // for vector
#include "Column.h"                            // for Column
#include "DebugLvl.hpp"                        // for debug_lvl
#include "NodeData.h"                          // for NodeData
#include "Parms.h"                             // for Parms
#include "PricerSolverBase.hpp"                // for PricerSolverBase
#include "PricingStabilization.hpp"            // for PricingStabili...
#include "Statistics.h"                        // for Statistics

/** Help function for column generation */
void NodeData::print_ages() {
    fmt::print("AGES:");

    std::ranges::for_each(localColPool,
                          [](auto const& it) { fmt::print(" {}", it->age); });

    fmt::print("\n");
}

auto NodeData::grow_ages() -> int {
    int val = 0;
    // lp_interface_get_nb_cols(RMP.get(), &nb_cols);
    auto nb_cols = osi_rmp->getNumCols();
    assert((nb_cols - id_pseudo_schedules) ==
           static_cast<int>(localColPool.size()));
    if (!localColPool.empty()) {
        column_status.resize(nb_cols);
        row_status.resize(osi_rmp->getNumRows());
        osi_rmp->getBasisStatus(column_status.data(), row_status.data());

        zero_count = 0;

        auto zipped = ranges::views::zip(column_status, localColPool);

        std::ranges::for_each(zipped, [&](auto&& it) {
            if (it.first == lp_interface_LOWER ||
                it.first == lp_interface_FREE) {
                it.second->age++;

                if (it.second->age > retirementage) {
                    zero_count++;
                }
            } else {
                it.second->age = 0;
            }
        });
    }

    return val;
}

auto NodeData::delete_unused_rows() -> int {
    int              val = 0;
    std::vector<int> del_indices{};

    auto osi_slack = std::span<const double>(osi_rmp->getRowActivity(),
                                             osi_rmp->getNumRows());

    int it = id_valid_cuts;
    int first_del = -1;
    int last_del = -1;
    for (auto s : osi_slack | ranges::views::drop(id_valid_cuts)) {
        if (std::fabs(s) < EPS) {
            if (first_del != -1) {
                val = delete_unused_rows_range(first_del, last_del);
                it = it - (last_del - first_del);
                first_del = last_del = -1;
            } else {
                it++;
            }
        } else {
            if (first_del == -1) {
                first_del = it;
                last_del = first_del;
            } else {
                last_del++;
            }
            it++;
        }
    }

    if (first_del != -1) {
        delete_unused_rows_range(first_del, last_del);
    }

    call_update_rows_coeff();

    return val;
}

auto NodeData::delete_old_columns() -> int {
    int val = 0;
    int min_numdel = static_cast<int>(
        floor(static_cast<int>(nb_jobs) * min_nb_del_row_ratio));
    /** pd->zero_count can be deprecated! */
    zero_count = 0;

    std::ranges::for_each(localColPool, [&](auto const& it) {
        if (it->age > 0) {
            zero_count++;
        }
    });

    if (zero_count > min_numdel) {
        int              iter = 0;
        std::vector<int> dellist{};
        // lp_interface_get_nb_cols(RMP.get(), &nb_cols);
        auto nb_cols = osi_rmp->getNumCols();
        assert((nb_cols - id_pseudo_schedules) ==
               static_cast<int>(localColPool.size()));

        localColPool |= ranges::actions::remove_if([&](const auto& it) {
            auto aux = false;
            if (it->age > retirementage) {
                dellist.emplace_back(iter + id_pseudo_schedules);
                aux = true;
            }
            iter++;
            return aux;
        });

        // lp_interface_delete_cols_array(RMP.get(), dellist.data(),
        //                                static_cast<int>(dellist.size()));
        osi_rmp->deleteCols(static_cast<int>(dellist.size()), dellist.data());

        if (debug_lvl(1)) {
            fmt::print("Deleted {} out of {} columns with age > {}.\n",
                       zero_count, localColPool.size(), retirementage);
        }
        // lp_interface_get_nb_cols(RMP.get(), &nb_cols);
        nb_cols = osi_rmp->getNumCols();
        assert(static_cast<int>(localColPool.size()) ==
               (nb_cols - id_pseudo_schedules));
        zero_count = 0;
        if (!dellist.empty()) {
            solve_relaxation();
        }
    }

    return val;
}

auto NodeData::delete_infeasible_columns() -> int {
    /** pd->zero_count can be deprecated! */
    zero_count = 0;

    int iter = 0;
    // nb_cols = osi_rmp->getNumCols();
    std::vector<int> dellist{};

    localColPool |= ranges::actions::remove_if([&](auto& it) {
        auto val = false;
        if (!solver->check_column(it.get())) {
            dellist.emplace_back(iter + id_pseudo_schedules);
            val = true;
        }
        iter++;
        return val;
    });

    // lp_interface_delete_cols_array(RMP.get(), dellist.data(),
    //                                static_cast<int>(dellist.size()));
    // lp_interface_get_nb_cols(RMP.get(), &nb_cols);
    osi_rmp->deleteCols(static_cast<int>(dellist.size()), dellist.data());
    auto nb_cols = osi_rmp->getNumCols();
    assert((nb_cols - id_pseudo_schedules) ==
           static_cast<int>(localColPool.size()));

    if (debug_lvl(1)) {
        fmt::print(
            "Deleted {} out of {} columns(infeasible columns after "
            "reduce "
            "cost "
            "fixing).\n",
            dellist.size(), localColPool.size());
    }

    // lp_interface_get_nb_cols(RMP.get(), &nb_cols);
    nb_cols = osi_rmp->getNumCols();
    assert(static_cast<int>(localColPool.size()) ==
           (nb_cols - id_pseudo_schedules));
    if (debug_lvl(2)) {
        fmt::print("number of cols = {}\n", nb_cols - id_pseudo_schedules);
    }

    if (!dellist.empty()) {
        solve_relaxation();
    }

    return 0;
}

// void g_make_pi_feasible(gpointer data, gpointer user_data) {
//     ScheduleSet* x = static_cast<ScheduleSet*>(data);
//     NodeData*    pd = static_cast<NodeData*>(user_data);
//     Job*         tmp_j = nullptr;

//     double colsum = .0;

//     for (guint i = 0; i < x->job_list->len; ++i) {
//         tmp_j = static_cast<Job*>(g_ptr_array_index(x->job_list, i));
//         if (signbit(pd->pi[tmp_j->job])) {
//             pd->pi[tmp_j->job] = 0.0;
//         }

//         colsum += pd->pi[tmp_j->job];
//         colsum = nextafter(colsum, DBL_MAX);
//     }

//     if (!signbit(pd->pi[pd->nb_jobs])) {
//         pd->pi[pd->nb_jobs] = 0.0;
//     }

//     colsum += pd->pi[pd->nb_jobs];
//     colsum = nextafter(colsum, DBL_MAX);

//     if (colsum > x->total_weighted_completion_time) {
//         double newcolsum = .0;
//         for (guint i = 0; i < x->job_list->len; ++i) {
//             tmp_j = static_cast<Job*>(g_ptr_array_index(x->job_list, i));
//             pd->pi[tmp_j->job] /= colsum;
//             pd->pi[tmp_j->job] *= x->total_weighted_completion_time;
//             newcolsum += pd->pi[tmp_j->job];
//         }

//         pd->pi[pd->nb_jobs] /= colsum;
//         pd->pi[pd->nb_jobs] *= x->total_weighted_completion_time;
//         newcolsum += pd->pi[pd->nb_jobs];

//         if (dbg_lvl() > 1) {
//             fmt::print(
//                 R"(Decreased column sum of {} from  {:30.20f} to
//                 {:30.20f}
// )",
//                 x->id, colsum, newcolsum);
//         }
//     }
// }

// MAYBE_UNUSED
// void make_pi_feasible(NodeData* pd) {
//     g_ptr_array_foreach(pd->localColPool, g_make_pi_feasible, pd);

//     std::for_each()
// }

// void g_make_pi_feasible_farkas(gpointer data, gpointer user_data) {
//     ScheduleSet* x = static_cast<ScheduleSet*>(data);
//     NodeData*    pd = static_cast<NodeData*>(user_data);

//     double colsum = .0;

//     for (guint i = 0; i < x->job_list->len; ++i) {
//         double tmp = pd->pi[i];
//         if (signbit(tmp)) {
//             tmp = 0.0;
//         }

//         colsum += tmp;
//         colsum = nextafter(colsum, DBL_MAX);
//     }

//     colsum += pd->pi[pd->nb_jobs];

//     if (colsum > x->total_weighted_completion_time) {
//         double newcolsum = .0;
//         for (guint i = 0; i < x->job_list->len; ++i) {
//             double tmp = pd->pi[i];
//             tmp /= colsum;
//             tmp *= x->total_weighted_completion_time;
//             newcolsum += tmp;
//         }

//         double tmp = pd->pi[pd->nb_jobs];
//         tmp /= colsum;
//         tmp *= x->total_weighted_completion_time;
//         newcolsum += tmp;

//         if (dbg_lvl() > 1) {
//             fmt::print(
//                 R"(Decreased column sum of {} from  {:30.20f} to
//                 {:30.20f}
// )",
//                 x->id, colsum, newcolsum);
//         }
//     }
// }

// MAYBE_UNUSED
// void NodeData::make_pi_feasible_farkas_pricing() {
//     g_ptr_array_foreach(localColPool, g_make_pi_feasible_farkas, this);
// }

auto NodeData::compute_objective() -> int {
    LP_lower_bound_dual = 0.0;

    /** compute lower bound with the dual variables */
    assert(osi_rmp->getNumRows() == static_cast<int>(pi.size()));

    LP_lower_bound_dual = ranges::inner_product(pi, rhs, 0.0);
    LP_lower_bound_dual -= EPS_BOUND;

    /** Get the LP lower bound and compute the lower bound of WCT */
    // lp_interface_objval(RMP.get(), &(LP_lower_bound));
    LP_lower_bound = osi_rmp->getObjValue();
    LP_lower_bound -= EPS_BOUND;
    // CCcheck_val_2(val, "lp_interface_objval failed");
    lower_bound = (ceil(LP_lower_bound_dual) < ceil(LP_lower_bound))
                      ? static_cast<int>(ceil(LP_lower_bound_dual))
                      : static_cast<int>(ceil(LP_lower_bound));
    LP_lower_bound_BB = std::min(LP_lower_bound, LP_lower_bound_dual);
    LP_lower_min = std::min(LP_lower_min, LP_lower_bound_BB);

    if (iterations % (nb_jobs) == 0 && debug_lvl(0)) {
        fmt::print(
            "Current primal LP objective: {:19.16f}  (LP_dual-bound "
            "{:19.16f}, "
            "lowerbound = {}, eta_in = {}, eta_out = {}).\n",
            LP_lower_bound + instance.off, LP_lower_bound_dual + instance.off,
            lower_bound + instance.off,
            solver_stab->get_eta_in() + instance.off,
            LP_lower_bound + instance.off);
    }

    return 0;
}

auto NodeData::solve_relaxation() -> int {
    auto val = 0;
    auto status_RMP = false;
    auto real_time_solve_lp = 0.0;

    /** Compute LP relaxation */
    real_time_solve_lp = getRealTime();
    stat.start_resume_timer(Statistics::solve_lp_timer);
    // val = lp_interface_optimize(RMP.get(), &status_RMP);
    osi_rmp->resolve();
    // status_RMP = osi_rmp->get
    // CCcheck_val_2(val, "lp_interface_optimize failed");
    stat.suspend_timer(Statistics::solve_lp_timer);

    if (debug_lvl(1)) {
        fmt::print("Simplex took {} seconds.\n", real_time_solve_lp);
    }

    if (debug_lvl(1)) {
        print_ages();
    }

    status_RMP = osi_rmp->isProvenOptimal();
    if (status_RMP) {
        grow_ages();
        pi = std::span<const double>{
            osi_rmp->getRowPrice(), static_cast<size_t>(osi_rmp->getNumRows())};
        compute_objective();
    } else if (osi_rmp->isProvenPrimalInfeasible()) {
        /**
         * @brief get dual variables
         *
         */
    }

    // switch (status_RMP) {
    //     case LP_INTERFACE_OPTIMAL:
    //         /** grow ages of the different columns */
    //         grow_ages();
    //         /** get the dual variables and make them feasible */
    //         lp_interface_pi(RMP.get(), pi.data());
    //         /** Compute the objective function */
    //         compute_objective();
    //         break;

    //     case LP_INTERFACE_INFEASIBLE:
    //         /** get the dual variables and make them feasible */
    //         val = lp_interface_pi_inf(RMP.get(), pi.data());
    //         break;
    // }

    return val;
}

auto NodeData::estimate_lower_bound(size_t _iter) -> int {
    auto val = 0;
    auto has_cols = true;
    auto status_RMP = false;

    if (debug_lvl(1)) {
        fmt::print(
            R"(Starting compute_lower_bound with lb {} and ub %d at depth {}
)",
            lower_bound, upper_bound, depth);
    }

    stat.start_resume_timer(Statistics::lb_timer);

    /**
     * Construction of new solution if localPoolColPool is empty
     */
    if (localColPool.empty()) {
        add_solution_to_colpool(opt_sol);
    }

    has_cols = true;
    // has_cuts = 0;
    while ((iterations < _iter) && has_cols &&
           stat.total_time_nano_sec(Statistics::cputime_timer) <=
               parms.branching_cpu_limit.value()) {
        /**
         * Solve the pricing problem
         */

        
            stat.start_resume_timer(Statistics::pricing_timer);
        status_RMP = osi_rmp->isProvenOptimal();

        if (status_RMP) {
            status = LP_bound_estimated;
            val = solve_pricing();
        } else if (osi_rmp->isProvenPrimalInfeasible()) {
            status = infeasible;
            solve_farkas_dbl();
        }

        ++iterations;

        stat.suspend_timer(Statistics::pricing_timer);

        if (status_RMP) {
            has_cols =
                (solver_stab->stopping_criteria() &&
                 solver_stab->get_eta_in() < upper_bound - 1.0 + EPS_BOUND);
        } else {
            has_cols = false;
            status = infeasible;
        }
    }

    stat.suspend_timer(Statistics::lb_timer);

    return val;
}

auto NodeData::compute_lower_bound() -> int {
    auto val = 0;
    auto has_cols{true};
    // auto has_cuts = 0;
    auto status_RMP = false;
    auto old_LP_bound = 0.0;

    if (debug_lvl(1)) {
        fmt::print(
            R"(Starting compute_lower_bound with lb {} and ub %d at depth {}
)",
            lower_bound, upper_bound, depth);
    }

    stat.start_resume_timer(Statistics::lb_timer);

    /**
     * Construction of new solution if localPoolColPool is empty
     */
    if (localColPool.empty()) {
        add_solution_to_colpool(opt_sol);
    }

    // delete_infeasible_columns();
    auto refined{false};

    solve_relaxation();
    do {
        has_cols = true;
        refined = false;
        while ((iterations < NB_CG_ITERATIONS) && has_cols &&
               solver->structure_feasible()) {
            /**
             * Delete old columns
             */
            status_RMP = osi_rmp->isProvenOptimal();
            if (zero_count > nb_jobs * min_nb_del_row_ratio && status_RMP) {
                delete_old_columns();
            }

            /**
             * Solve the pricing problem
             */
            stat.start_resume_timer(Statistics::pricing_timer);
            // val = lp_interface_status(RMP.get(), &status_RMP);

            if (status_RMP) {
                ++iterations;
                status = LP_bound_estimated;
                solve_pricing();
            } else {
                status = infeasible;
                solve_farkas_dbl();
            }

            // switch (status_RMP) {
            //     case GRB_OPTIMAL:
            //         ++iterations;
            //         status = infeasible;
            //         val = solve_pricing();
            //         break;

            //     case GRB_INFEASIBLE:
            //         solve_farkas_dbl();
            //         break;
            // }

            stat.suspend_timer(Statistics::pricing_timer);

            if (status_RMP) {
                // case GRB_OPTIMAL:
                has_cols =
                    (solver_stab->stopping_criteria() &&
                     solver_stab->get_eta_in() < upper_bound - 1.0 + EPS_BOUND);
                // nb_new_sets = 0;
                // || nb_non_improvements > 5;  // ||
                // (ceil(eta_in - 0.00001) >= eta_out);

                //     break;

                // case GRB_INFEASIBLE:
                //     has_cols = (nb_new_sets == 0);
                //     // nb_new_sets = 0;
                //     break;
            }
        }

        if (status_RMP) {
            // case GRB_OPTIMAL:
            if (debug_lvl(1)) {
                fmt::print(
                    R"(Found lb = {} ({}) upper_bound = {} (iterations = {}).
)",
                    lower_bound, LP_lower_bound, upper_bound, iterations);
            }

            /**
             * Compute the objective function
             */
            solve_relaxation();

            if (!localColPool.empty() && solver->structure_feasible()) {
                status = LP_bound_computed;
                construct_lp_sol_from_rmp();
                if (parms.suboptimal_duals.value()) {
                    refined =
                        solver->compute_sub_optimal_duals(lambda, localColPool);
                    delete_infeasible_columns();
                    if (refined) {
                        status = LP_bound_estimated;
                    }
                }
                if (parms.refine_bdd.value() &&
                    nb_non_improvements < NB_NON_IMPROVEMENTS && !refined) {
                    if (std::abs(old_LP_bound - LP_lower_bound) < EPS) {
                        nb_non_improvements++;
                    } else {
                        nb_non_improvements = 0;
                    }
                    old_LP_bound = LP_lower_bound;
                    refined = refinement();
                    if (refined) {
                        status = LP_bound_estimated;
                    }
                }
                generate_cuts();
            } else {
                status = infeasible;
                LP_lower_bound_dual = LP_lower_bound = LP_lower_bound_BB =
                    upper_bound;
                lower_bound = upper_bound;
            }
            //     break;

            // case GRB_INFEASIBLE:
            //     status = infeasible;
            //     // lp_interface_write(RMP.get(), "infeasible_RMP.lp");
            //     osi_rmp->writeLp("infeasible_RMP");
            //     // lp_interface_compute_IIS(RMP.get());
        }
    } while (refined);

    if (solver->print_num_paths() < NB_PATHS && parms.enumerate.value() &&
        status == LP_bound_computed) {
        fmt::print("computing elementary paths\n");
        solver->enumerate_columns();
        solver->enumerate_columns(pi);
    }

    // if (iterations < NB_CG_ITERATIONS &&
    //     stat.total_timer(Statistics::cputime_timer) <=
    //         parms.branching_cpu_limit) {
    // } else {
    //     switch (status_RMP) {
    //         case GRB_OPTIMAL:
    //             status = LP_bound_computed;
    //             break;

    //         case GRB_INFEASIBLE:
    //             status = infeasible;
    //             break;
    //     }
    // }
    // } while (depth == 1);

    if (depth == 0UL) {
        stat.global_lower_bound =
            std::max(lower_bound + instance.off, stat.global_lower_bound);
        stat.root_lower_bound = stat.global_lower_bound;
        stat.root_upper_bound = stat.global_upper_bound;
        stat.root_rel_error = static_cast<double>(stat.global_upper_bound -
                                                  stat.global_lower_bound) /
                              (stat.global_lower_bound + EPS);
        stat.size_graph_after_reduced_cost_fixing = solver->get_nb_edges();
        stat.nb_generated_col_root = iterations;
    }

    stat.nb_generated_col += iterations;
    stat.suspend_timer(Statistics::lb_timer);

    return val;
}

auto NodeData::refinement() -> bool {
    // auto status_RMP = 0;
    auto refined{false};

    // lp_interface_status(RMP.get(), &status_RMP);
    std::vector<std::pair<std::shared_ptr<Column>, double>> paths;
    if (osi_rmp->isProvenOptimal()) {
        // switch (status_RMP) {
        // case GRB_OPTIMAL:
        // lp_interface_get_nb_cols(RMP.get(), &nb_cols);
        auto nb_cols = osi_rmp->getNumCols();
        assert(static_cast<int>(localColPool.size()) ==
               (nb_cols - id_pseudo_schedules));
        // lambda.resize(localColPool.size(), 0.0);
        // lp_interface_x(RMP.get(), lambda.data(), id_pseudo_schedules);
        lambda = std::span<const double>{osi_rmp->getColSolution(),
                                         static_cast<size_t>(nb_cols)};

        for (auto&& [it, x] :
             ranges::views::zip(
                 localColPool,
                 lambda | ranges::views::drop(id_pseudo_schedules)) |
                 //  ranges::views::zip(localColPool, lambda) |
                 ranges::views::filter(
                     [](const auto x) { return x > EPS; },
                     [](const auto& tmp) { return tmp.second; })) {
            paths.emplace_back(it, x);
        }

        ranges::sort(paths, std::greater<>{},
                     [](const auto& tmp) { return tmp.second; });

        refined = solver->refinement_structure(
            paths | ranges::views::transform([](const auto& tmp) {
                return tmp.first;
            }) |
            ranges::to_vector);

        if (refined) {
            delete_infeasible_columns();
            solve_relaxation();
        }
        // break;
    }

    return refined;
}

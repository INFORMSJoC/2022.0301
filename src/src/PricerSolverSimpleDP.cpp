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

#include "PricerSolverSimpleDP.hpp"
#include <fmt/core.h>
#include <gurobi_c++.h>                            // for GRBLinExpr, GRBModel
#include <algorithm>                               // for fill, find, max
#include <boost/graph/adjacency_list.hpp>          // for adjacency_list
#include <boost/graph/detail/adjacency_list.hpp>   // for edges, get, num_edges
#include <boost/graph/detail/edge.hpp>             // for operator!=, operat...
#include <boost/graph/graph_selectors.hpp>         // for bidirectionalS
#include <boost/graph/graphviz.hpp>                // for write_graphviz
#include <boost/iterator/iterator_facade.hpp>      // for operator!=, operat...
#include <boost/multiprecision/cpp_int.hpp>        // for cpp_int
#include <boost/pending/property.hpp>              // for no_property
#include <cstddef>                                 // for size_t
#include <fstream>                                 // for operator<<, basic_...
#include <iostream>                                // for cerr
#include <limits>                                  // for numeric_limits
#include <list>                                    // for operator==
#include <range/v3/iterator/basic_iterator.hpp>    // for basic_iterator
#include <range/v3/iterator/reverse_iterator.hpp>  // for reverse_cursor
#include <range/v3/view/iota.hpp>                  // for iota_view, ints
#include <range/v3/view/reverse.hpp>               // for reverse_fn, revers...
#include <range/v3/view/view.hpp>                  // for operator|, view_cl...
#include <span>                                    // for span
#include <string>                                  // for char_traits, opera...
#include <utility>                                 // for move
#include <vector>                                  // for vector, vector<>::...
#include "Column.h"                                // for ScheduleSet
#include "Instance.h"                              // for Instance
#include "Job.h"                                   // for Job
#include "ModelInterface.hpp"                      // for ReformulationModel
#include "PricerSolverBase.hpp"                    // for PricerSolverBase
/**
 * Pricersolver for the TI index formulation
 */
PricerSolverSimpleDp::PricerSolverSimpleDp(const Instance& instance)
    : PricerSolverBase(instance),
      Hmax(instance.H_max),
      size_graph{},
      A(Hmax + 1),
      F(Hmax + 1),
      backward_F(Hmax + 1),
      TI_x(convex_constr_id * (Hmax + 1), GRBVar()),
      take(convex_constr_id * (Hmax + 1), false),
      lp_x(convex_constr_id * (Hmax + 1), 0.0),
      solution_x(convex_constr_id * (Hmax + 1), 0.0) {
    init_table();
}

void PricerSolverSimpleDp::init_table() {
    backward_graph = std::vector<std::vector<Job*>>(Hmax + 1);
    forward_graph = std::vector<std::vector<Job*>>(Hmax + 1);

    for (auto t = 0UL; t < Hmax + 1; t++) {
        for (auto i = 1UL; i < convex_constr_id + 1; i++) {
            auto  j = i - 1;
            auto* job = jobs[j].get();

            if (t >= static_cast<size_t>(job->processing_time)) {
                forward_graph[t].push_back(job);
                size_graph++;
            }

            if (t + job->processing_time <= Hmax) {
                backward_graph[t].push_back(job);
            }
        }
    }

    fmt::print("Number of arcs in TI formulation = {}\n", size_graph);
}

auto PricerSolverSimpleDp::evaluate_nodes(
    [[maybe_unused]] std::span<const double>& pi) -> bool {
    forward_evaluator(pi);
    backward_evaluator(pi);

    auto counter = size_t{};
    auto x = size_t{};
    auto num_removed = 0;

    for (auto t = size_t{}; t < Hmax + 1; ++t) {
        auto it = forward_graph[t].begin();
        while (it != forward_graph[t].end()) {
            double result = F[t - (*it)->processing_time] +
                            (*it)->weighted_tardiness(t) - pi[(*it)->job] +
                            backward_F[t];

            if (constLB + result +
                    static_cast<double>(convex_rhs - 1) * F[Hmax] >
                UB - 1.0 + RC_FIXING) {
                --size_graph;
                it = forward_graph[t].erase(it);
                ++num_removed;
            } else {
                it++;
            }
        }

        x += forward_graph[t].size();

        auto iter = backward_graph[t].begin();
        while (iter != backward_graph[t].end()) {
            double result = F[t] + (*iter)->weighted_tardiness_start(t) -
                            pi[(*iter)->job] +
                            backward_F[t + (*iter)->processing_time];
            if (constLB + result +
                    static_cast<double>(convex_rhs - 1) * F[Hmax] >
                UB - 1.0 + RC_FIXING) {
                take[(*iter)->job * (Hmax + 1) + t] = false;
                iter = backward_graph[t].erase(iter);
            } else {
                take[(*iter)->job * (Hmax + 1) + t] = true;
                iter++;
            }
        }

        counter += backward_graph[t].size();
    }

    return (num_removed > 0);
}


void PricerSolverSimpleDp::build_mip() {
    try {
        fmt::print("Building Mip model for the TI formulation\n");

        /** Constructing variables */
        for (auto t = size_t{}; t < Hmax + 1; t++) {
            for (auto& it : backward_graph[t]) {
                double cost = it->weighted_tardiness_start(t);
                double ub = take[(it->job) * (Hmax + 1) + t] ? 1.0 : 0.0;
                TI_x[it->job * (Hmax + 1) + t] =
                    model.addVar(0.0, ub, cost, 'B');
            }
        }

        model.update();

        /** Assignment variables */
        std::vector<GRBLinExpr> assignment(convex_constr_id, GRBLinExpr());
        std::vector<char>       sense(convex_constr_id, '=');
        std::vector<double>     rhs(convex_constr_id, 1.0);

        for (auto t = 0UL; t <= Hmax; t++) {
            for (auto& it : backward_graph[t]) {
                assignment[it->job] += TI_x[it->job * (Hmax + 1) + t];
            }
        }

        std::unique_ptr<GRBConstr> assignment_constrs(
            model.addConstrs(assignment.data(), sense.data(), rhs.data(),
                             nullptr, static_cast<int>(convex_constr_id)));

        model.update();

        std::vector<GRBLinExpr> interval_constr(Hmax + 1, GRBLinExpr());
        std::vector<char>       interval_sense(Hmax + 1);
        std::vector<double>     interval_rhs(Hmax + 1);

        for (auto t = 0UL; t <= Hmax; t++) {
            auto add_constraint = false;
            for (auto& it : backward_graph[t]) {
                for (auto s = std::max(0UL, t - it->processing_time); s <= t;
                     s++) {
                    if (std::find(backward_graph[s].begin(),
                                  backward_graph[s].end(),
                                  it) != backward_graph[s].end()) {
                        interval_constr[t] += TI_x[(it->job) * (Hmax + 1) + s];
                        add_constraint = true;
                    }
                }
            }

            if (add_constraint) {
                interval_sense[t] = '<';
                interval_rhs[t] = static_cast<double>(convex_rhs);

                model.addConstr(interval_constr[t], interval_sense[t],
                                interval_rhs[t]);
            }
        }

        model.update();

    } catch (GRBException& e) {
        std::cerr << e.getMessage() << '\n';
    }

    for (auto t = 0UL; t <= Hmax; t++) {
        for (auto& it : backward_graph[t]) {
            TI_x[it->job * (Hmax + 1) + t].set(
                GRB_DoubleAttr_Start, solution_x[(it->job) * (Hmax + 1) + t]);
        }
    }

    model.write("ti_" + problem_name + "_" + std::to_string(convex_rhs) +
                "correct.lp");
    model.optimize();
}

void PricerSolverSimpleDp::forward_evaluator(std::span<const double>& _pi) {
    /** Initialisation */
    F[0] = 0.0;
    A[0] = nullptr;

    for (auto t = 1UL; t < Hmax + 1; t++) {
        F[t] = std::numeric_limits<double>::max() / 2;
        A[t] = nullptr;
    }

    /** Recursion */
    for (auto t = size_t{1}; t < Hmax + 1; t++) {
        for (auto& it : forward_graph[t]) {
            if (F[t - it->processing_time] +
                    static_cast<double>(it->weighted_tardiness(t)) -
                    _pi[it->job] <=
                F[t]) {
                F[t] = F[t - it->processing_time] + it->weighted_tardiness(t) -
                       _pi[it->job];
                A[t] = it;
            }
        }

        if (F[t - 1] <= F[t]) {
            F[t] = F[t - 1];
        }
    }
}

void PricerSolverSimpleDp::forward_evaluator(double* _pi) {
    /** Initialisation */
    std::span aux_pi{_pi, reformulation_model.size()};
    F[0] = 0.0;
    A[0] = nullptr;

    for (auto t = 1UL; t < Hmax + 1; t++) {
        F[t] = std::numeric_limits<double>::max() / 2;
        A[t] = nullptr;
    }

    /** Recursion */
    for (auto t = size_t{1}; t < Hmax + 1; t++) {
        for (auto& it : forward_graph[t]) {
            if (F[t - it->processing_time] +
                    static_cast<double>(it->weighted_tardiness(t)) -
                    aux_pi[it->job] <=
                F[t]) {
                F[t] = F[t - it->processing_time] + it->weighted_tardiness(t) -
                       aux_pi[it->job];
                A[t] = it;
            }
        }

        if (F[t - 1] <= F[t]) {
            F[t] = F[t - 1];
        }
    }
}

void PricerSolverSimpleDp::backward_evaluator(std::span<const double>& _pi) {
    backward_F[Hmax] = 0.0;

    for (auto t = 0UL; t < Hmax; t++) {
        backward_F[t] = std::numeric_limits<double>::max() / 2;
    }

    for (auto t :
         ranges::views::ints(size_t{}, Hmax) | ranges::views::reverse) {
        for (auto& it : backward_graph[t]) {
            auto tt = t + it->processing_time;
            if (backward_F[tt] +
                    static_cast<double>(it->weighted_tardiness(tt)) -
                    _pi[it->job] <=
                backward_F[t]) {
                backward_F[t] =
                    backward_F[tt] +
                    static_cast<double>(it->weighted_tardiness(tt)) -
                    _pi[it->job];
            }
        }

        if (backward_F[t + 1] <= backward_F[t]) {
            backward_F[t] = backward_F[t + 1];
        }
    }
}

void PricerSolverSimpleDp::backward_evaluator(double* _pi) {
    backward_F[Hmax] = 0.0;
    std::span aux_pi{_pi, reformulation_model.size()};

    for (auto t = 0UL; t < Hmax; t++) {
        backward_F[t] = std::numeric_limits<double>::max() / 2;
    }

    for (auto t :
         ranges::views::ints(size_t{}, Hmax) | ranges::views::reverse) {
        for (auto& it : backward_graph[t]) {
            auto tt = t + it->processing_time;
            if (backward_F[tt] +
                    static_cast<double>(it->weighted_tardiness(tt)) -
                    aux_pi[it->job] <=
                backward_F[t]) {
                backward_F[t] =
                    backward_F[tt] +
                    static_cast<double>(it->weighted_tardiness(tt)) -
                    aux_pi[it->job];
            }
        }

        if (backward_F[t + 1] <= backward_F[t]) {
            backward_F[t] = backward_F[t + 1];
        }
    }
}

auto PricerSolverSimpleDp::pricing_algorithm(std::span<const double>& _pi)
    -> PricingSolution {
    PricingSolution opt_sol;
    opt_sol.cost = 0;
    std::vector<Job*> v;

    forward_evaluator(_pi);

    /** Find optimal solution */
    opt_sol.obj = std::numeric_limits<double>::max();

    for (auto i = 0UL; i < Hmax + 1; i++) {
        if (F[i] < opt_sol.obj) {
            opt_sol.C_max = i;
            opt_sol.obj = F[i];
        }
    }

    auto t_min = opt_sol.C_max;

    /** Construct the solution */
    while (A[t_min] != nullptr) {
        Job* job = A[t_min];
        v.push_back(A[t_min]);
        opt_sol.cost += A[t_min]->weighted_tardiness(t_min);
        t_min -= job->processing_time;
    }

    auto it = v.rbegin();

    for (; it != v.rend(); ++it) {
        opt_sol.jobs.push_back(*it);
    }

    /** Free the memory */
    return opt_sol;
}


auto PricerSolverSimpleDp::farkas_pricing(
    [[maybe_unused]] std::span<const double>& _pi) -> PricingSolution {
    PricingSolution opt_sol;

    return opt_sol;
}

void PricerSolverSimpleDp::construct_lp_sol_from_rmp(
    const std::span<const double>&              lambda,
    const std::vector<std::shared_ptr<Column>>& columns) {
    // std::span aux_columns{columns->pdata, columns->len};
    std::fill(lp_x.begin(), lp_x.end(), 0.0);
    for (auto k = 0UL; k < columns.size(); k++) {
        if (lambda[k] > EPS_SOLVER) {
            auto* tmp = columns[k].get();
            int   t = 0;
            // std::span aux_jobs{tmp->job_list->pdata, tmp->job_list->len};
            for (auto& it : tmp->job_list) {
                lp_x[(it->job) * (Hmax + 1) + t] += lambda[k];
                t += it->processing_time;
            }
        }
    }
}

// void PricerSolverSimpleDp::create_dot_zdd([[maybe_unused]] const char* name)
// {
//     boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS>
//         graph;
//     for (auto t = 0UL; t <= Hmax; t++) {
//         boost::add_vertex(graph);
//     }

//     for (auto t = 0UL; t < Hmax; t++) {
//         for (auto& it : backward_graph[t]) {
//             boost::add_edge(t, t + it->processing_time, graph);
//         }
//     }
//     auto file_name = "TI_representation_" + problem_name + "_" +
//                      std::to_string(convex_rhs) + ".gv";
//     auto otf = std::ofstream(file_name);
//     boost::write_graphviz(otf, graph);
//     otf.close();
// }

auto PricerSolverSimpleDp::get_nb_edges() -> size_t {
    size_t nb_edges = 0U;
    for (auto& it : forward_graph) {
        nb_edges += it.size();
    }
    return nb_edges;
}

auto PricerSolverSimpleDp::get_nb_vertices() -> size_t {
    size_t nb_vertices = 0U;
    for (auto& it : forward_graph) {
        if (!it.empty()) {
            nb_vertices++;
        }
    }
    return nb_vertices;
}

auto PricerSolverSimpleDp::print_num_paths() -> boost::multiprecision::cpp_int {
    return 0;
}

auto PricerSolverSimpleDp::check_column([[maybe_unused]] Column const* set)
    -> bool {
    // int t = 0;
    // for(unsigned int j = 0; j < set->len; j++) {
    //     Job* tmp_j = static_cast<Job*>() g_ptr_array_index()atic_cast<set,
    //     j>())); t += tmp_j->processing_time;

    //     if (t > Hmax) {
    //         return false;
    //     }

    //     auto it = std::find(forward_graph[t].begin(),
    //     forward_graph[t].end(),tmp_j); if (it == forward_graph[t].end()) {
    //         return false;
    //     }

    // }

    return true;
}

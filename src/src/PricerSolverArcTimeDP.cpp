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

#include "PricerSolverArcTimeDP.hpp"
#include <limits>                       // for numeric_limits
#include <range/v3/algorithm/find.hpp>  // for find
#include <range/v3/view/filter.hpp>     // for filter
#include <range/v3/view/reverse.hpp>    // for reverse_fn, revers...
#include <range/v3/view/take.hpp>       // for take
#include <range/v3/view/zip.hpp>        // for zip
#include <span>                         // for span
#include <string>                       // for char_traits, opera...
#include <vector>                       // for vector, vector<>::...
#include "Column.h"                     // for Column
#include "Instance.h"                   // for Instance
#include "PricerSolverBase.hpp"         // for PricerSolverBase
#include "gurobi_c++.h"                 // for GRBLinExpr, GRBModel
#include "gurobi_c.h"                   // for GRB_EQUAL, GRB_BINARY

using ranges::to_vector;

PricerSolverArcTimeDp::PricerSolverArcTimeDp(const Instance& instance)
    : PricerSolverBase(instance),
      Hmax(instance.H_max),
      n(instance.nb_jobs),
      size_graph(0U),
      vector_jobs(instance.jobs | ranges::views::transform([](const auto& tmp) {
                      return tmp.get();
                  }) |
                  to_vector),
      nb_edges_removed{},
      lp_x((n + 1) * (n + 1) * (Hmax + 1), 0.0),
      solution_x((n + 1) * (n + 1) * (Hmax + 1), 0.0) {
    // for (auto i : ranges::views::ints(0UL, n)) {V
    //     vector_jobs.push_back((jobs)[i].get());
    // }
    j0.job = n;
    vector_jobs.push_back(&j0);

    init_table();
}

void PricerSolverArcTimeDp::init_table() {
    graph = vector3d_jobs(n + 1);
    reversed_graph = vector3d_jobs(n + 1);
    for (auto j = 0UL; j < n + 1; j++) {
        reversed_graph[j] = vector2d_jobs(Hmax + 1);
    }

    forward_F = vector2d_dbl(convex_constr_id + 1);
    backward_F = vector2d_dbl(convex_constr_id + 1);

    for (unsigned i = 0; i < convex_constr_id + 1; ++i) {
        forward_F[i] = vector1d_dbl(Hmax + 1, 0.0);
        backward_F[i] = vector1d_dbl(Hmax + 1, 0.0);
    }

    A = vector2d_jobs(convex_constr_id + 1);
    for (unsigned i = 0; i < convex_constr_id + 1; ++i) {
        A[i] = vector1d_jobs(Hmax + 1, nullptr);
    }

    B = vector2d_int(convex_constr_id + 1);
    for (unsigned i = 0; i < convex_constr_id + 1; ++i) {
        B[i] = vector1d_int(Hmax + 1);
    }

    arctime_x = vector_grb_var(convex_constr_id + 1);
    for (auto i = 0UL; i < n + 1; i++) {
        arctime_x[i] = vector2d_grb_var(convex_constr_id + 1);
        for (auto j = 0UL; j < n + 1; j++) {
            arctime_x[i][j] = vector1d_grb_var(Hmax + 1);
        }
    }

    for (auto j = 0UL; j < n; ++j) {
        graph[j] = vector2d_jobs(Hmax + 1);
        Job* tmp = jobs[j].get();
        for (auto t : ranges::views::ints(size_t{}, Hmax + 1)) {
            for (auto& it : vector_jobs) {
                int p = (it->job != j) ? it->processing_time : 1;
                if (it != tmp && static_cast<int>(t) >= p &&
                    t + tmp->processing_time <= Hmax) {
                    graph[j][t].push_back(it);
                    ++size_graph;
                }
            }
        }
    }

    graph[n] = vector2d_jobs(Hmax + 1);
    for (auto t : ranges::views::ints(size_t{}, Hmax + 1)) {
        for (auto& it : vector_jobs) {
            if (static_cast<int>(t) >= it->processing_time) {
                graph[n][t].push_back(it);
                ++size_graph;
            }
        }
    }

    /**
     * Remove all not needed arcs from the sets
     */
    for (auto i = 0UL; i < n - 1; ++i) {
        Job* tmp_i = vector_jobs[i];
        for (auto j = i + 1; j < n; ++j) {
            Job* tmp_j = vector_jobs[j];
            for (size_t t = tmp_i->processing_time;
                 t <= Hmax - tmp_j->processing_time; ++t) {
                if (delta1(i, j, static_cast<int>(t)) >= 0) {
                    remove_arc(i, j, t);
                    size_graph--;
                } else {
                    remove_arc(
                        j, i,
                        t - tmp_i->processing_time + tmp_j->processing_time);
                    size_graph--;
                }
            }
        }
    }

    for (auto j = 0UL; j < n; ++j) {
        Job* tmp_j = vector_jobs[j];
        for (auto t : ranges::views::ints(
                 static_cast<size_t>(tmp_j->processing_time), Hmax)) {
            if (delta2(j, t) <= 0) {
                remove_arc(n, j, t - tmp_j->processing_time + 1);
                --size_graph;
            } else {
                remove_arc(j, n, t);
                --size_graph;
            }
        }
    }

    for (auto j = 0UL; j < n + 1; ++j) {
        Job* tmp = vector_jobs[j];
        for (auto t : ranges::views::ints(size_t{}, Hmax + 1)) {
            if (graph[j][t].empty()) {
                forward_F[j][t] = std::numeric_limits<double>::max() / 2;
            } else {
                for (auto& it : graph[j][t]) {
                    reversed_graph[it->job][t].push_back(tmp);
                }
            }
        }
    }

    for (auto j = 0UL; j < n + 1; j++) {
        for (auto t : ranges::views::ints(size_t{}, Hmax + 1)) {
            if (reversed_graph[j][t].empty()) {
                backward_F[j][t] = std::numeric_limits<double>::max() / 2;
            }
        }
    }

    fmt::print("Number of arcs in ATI formulation = {}\n", size_graph);
}


auto PricerSolverArcTimeDp::evaluate_nodes(
    [[maybe_unused]] std::span<const double>& pi) -> bool {
    forward_evaluator(pi);
    backward_evaluator(pi);
    auto num_edges_removed = 0;

    for (auto* tmp : vector_jobs | ranges::views::take(n)) {
        for (auto t = size_t{}; t <= Hmax - tmp->processing_time; t++) {
            auto it = graph[tmp->job][t].begin();
            while (it != graph[tmp->job][t].end()) {
                double result =
                    forward_F[(*it)->job][t - (*it)->processing_time] +
                    tmp->weighted_tardiness_start(t) - pi[tmp->job] +
                    backward_F[tmp->job][t + tmp->processing_time];
                if (result +
                        static_cast<double>(convex_rhs - 1) *
                            (forward_F[n][Hmax]) +
                        constLB >
                    UB - 1.0 + RC_FIXING) {
                    size_graph--;
                    num_edges_removed++;
                    auto pend =
                        ranges::find(reversed_graph[(*it)->job][t], tmp);
                    reversed_graph[(*it)->job][t].erase(pend);
                    it = graph[tmp->job][t].erase(it);
                } else {
                    it++;
                }
            }
        }
    }

    return (num_edges_removed > 0);
}

void PricerSolverArcTimeDp::build_mip() {
    fmt::print("Building Mip model for the arcTI formulation\n");

    /** Constructing variables */
    for (auto j = size_t{}; j < n + 1; j++) {
        for (auto t = size_t{}; t + vector_jobs[j]->processing_time <= Hmax;
             t++) {
            for (auto& it : graph[j][t]) {
                double cost = vector_jobs[j]->weighted_tardiness_start(
                    static_cast<int>(t));
                double tmp = (it->job == vector_jobs[j]->job)
                                 ? static_cast<double>(convex_rhs)
                                 : 1.0;
                auto   s =
                    (it->job == vector_jobs[j]->job) ? GRB_INTEGER : GRB_BINARY;
                arctime_x[it->job][j][t] = model.addVar(0.0, tmp, cost, s);
            }
        }
    }

    model.update();

    /** Assignment variables */
    std::vector<GRBLinExpr> assignment(convex_constr_id, GRBLinExpr());
    std::vector<char>       sense(convex_constr_id, GRB_GREATER_EQUAL);
    std::vector<double>     rhs(convex_constr_id, 1.0);

    for (auto j = size_t{}; j < n; j++) {
        for (auto t = size_t{}; t <= Hmax - vector_jobs[j]->processing_time;
             t++) {
            for (auto& it : graph[j][t]) {
                assignment[j] += arctime_x[it->job][j][t];
            }
        }
    }

    std::unique_ptr<GRBConstr> assignment_constrs(
        model.addConstrs(assignment.data(), sense.data(), rhs.data(), nullptr,
                         static_cast<int>(convex_constr_id)));

    for (auto i = 0UL; i < n; i++) {
        for (auto t = 0UL; t <= Hmax - vector_jobs[i]->processing_time; t++) {
            GRBLinExpr expr{};
            for (auto& it : graph[i][t]) {
                expr += arctime_x[it->job][i][t];
            }

            for (auto& it :
                 reversed_graph[i][t + vector_jobs[i]->processing_time]) {
                expr -=
                    arctime_x[i][it->job][t + vector_jobs[i]->processing_time];
            }
            model.addConstr(expr, GRB_EQUAL, 0);
        }
    }

    for (auto t : ranges::views::ints(size_t{}, Hmax)) {
        GRBLinExpr expr{};
        for (auto& it : graph[n][t]) {
            expr += arctime_x[it->job][n][t];
        }

        for (auto& it : reversed_graph[n][t + 1]) {
            expr -= arctime_x[n][it->job][t + 1];
        }

        model.addConstr(expr, GRB_EQUAL, 0);
    }

    GRBLinExpr expr{};
    for (auto& it : reversed_graph[n][0]) {
        expr += arctime_x[n][it->job][0];
    }
    model.addConstr(expr, GRB_EQUAL, static_cast<double>(convex_rhs));

    for (auto j = 0UL; j < n + 1; j++) {
        for (auto t : ranges::views::ints(
                 size_t{}, Hmax - vector_jobs[j]->processing_time + 1)) {
            for (auto& it : graph[j][t]) {
                arctime_x[it->job][j][t].set(
                    GRB_DoubleAttr_Start,
                    solution_x[(it->job) * (Hmax + 1) * (n + 1) +
                               j * (Hmax + 1) + t]);
                arctime_x[it->job][j][t].set(
                    GRB_DoubleAttr_PStart,
                    lp_x[(it->job) * (Hmax + 1) * (n + 1) + j * (Hmax + 1) +
                         t]);
            }
        }
    }

    model.write("ati_" + problem_name + "_" + std::to_string(convex_rhs) +
                ".lp");
    model.optimize();

    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        for (auto j = 0UL; j < n; j++) {
            for (auto t = 0UL; t <= Hmax - vector_jobs[j]->processing_time;
                 t++) {
                for (auto& it : graph[j][t]) {
                    auto a = arctime_x[it->job][j][t].get(GRB_DoubleAttr_X);
                    if (a > 0) {
                        fmt::print("{} {} {}\n", j, t,
                                   jobs[j]->processing_time);
                    }
                }
            }
        }
    }
}

PricerSolverArcTimeDp::~PricerSolverArcTimeDp() = default;

void PricerSolverArcTimeDp::backward_evaluator(double* _pi) {
    // backward_F[n][T] = 0;
    backward_F[n][Hmax] = 0.0;
    std::span aux_pi{_pi, reformulation_model.size()};

    for (auto t :
         ranges::views::ints(size_t{}, Hmax) | ranges::views::reverse) {
        for (auto i = 0UL; i <= n; ++i) {
            // Job* tmp = vector_jobs[i];
            backward_F[i][t] = ((i == n) && (t == Hmax))
                                   ? 0.0
                                   : std::numeric_limits<double>::max() / 2;
            auto it = reversed_graph[i][t].begin();

            if (!reversed_graph[i][t].empty() && t <= Hmax) {
                double reduced_cost = ((*it)->job == n)
                                          ? (*it)->weighted_tardiness_start(t)
                                          : (*it)->weighted_tardiness_start(t) -
                                                aux_pi[(*it)->job];
                int    tt = ((*it)->job != n) ? (*it)->processing_time
                            : (*it)->job == i ? 1
                                              : 0;
                backward_F[i][t] =
                    backward_F[(*it)->job][t + tt] + reduced_cost;
                it++;
                while (it != reversed_graph[i][t].end()) {
                    reduced_cost = ((*it)->job == n)
                                       ? (*it)->weighted_tardiness_start(t)
                                       : (*it)->weighted_tardiness_start(t) -
                                             aux_pi[(*it)->job];

                    tt = ((*it)->job != n) ? (*it)->processing_time
                         : (*it)->job == i ? 1
                                           : 0;
                    double result =
                        backward_F[(*it)->job][t + tt] + reduced_cost;

                    if (backward_F[i][t] >= result) {
                        backward_F[i][t] = result;
                    }
                    it++;
                }

            } else {
                backward_F[i][t] = std::numeric_limits<double>::max() / 2;
            }
        }
    }
}

void PricerSolverArcTimeDp::backward_evaluator(std::span<const double>& _pi) {
    // backward_F[n][T] = 0;
    backward_F[n][Hmax] = 0.0;

    for (auto t :
         ranges::views::ints(size_t{}, Hmax) | ranges::views::reverse) {
        for (auto i = 0UL; i <= n; ++i) {
            // Job* tmp = vector_jobs[i];
            backward_F[i][t] = ((i == n) && (t == Hmax))
                                   ? 0.0
                                   : std::numeric_limits<double>::max() / 2;
            auto it = reversed_graph[i][t].begin();

            if (!reversed_graph[i][t].empty() && t <= Hmax) {
                double reduced_cost =
                    ((*it)->job == n)
                        ? (*it)->weighted_tardiness_start(t)
                        : (*it)->weighted_tardiness_start(t) - _pi[(*it)->job];
                int tt = ((*it)->job != n) ? (*it)->processing_time
                         : (*it)->job == i ? 1
                                           : 0;
                backward_F[i][t] =
                    backward_F[(*it)->job][t + tt] + reduced_cost;
                it++;
                while (it != reversed_graph[i][t].end()) {
                    reduced_cost = ((*it)->job == n)
                                       ? (*it)->weighted_tardiness_start(t)
                                       : (*it)->weighted_tardiness_start(t) -
                                             _pi[(*it)->job];

                    tt = ((*it)->job != n) ? (*it)->processing_time
                         : (*it)->job == i ? 1
                                           : 0;
                    double result =
                        backward_F[(*it)->job][t + tt] + reduced_cost;

                    if (backward_F[i][t] >= result) {
                        backward_F[i][t] = result;
                    }
                    it++;
                }

            } else {
                backward_F[i][t] = std::numeric_limits<double>::max() / 2;
            }
        }
    }
}

void PricerSolverArcTimeDp::forward_evaluator(double* _pi) {
    forward_F[n][0] = 0;
    std::span aux_pi{_pi, reformulation_model.size()};

    for (auto t : ranges::views::ints(int{}, static_cast<int>(Hmax) + 1)) {
        for (auto j = 0UL; j <= n; ++j) {
            Job* tmp = vector_jobs[j];
            A[j][t] = nullptr;
            B[j][t] = -1;
            forward_F[j][t] = ((j == n) && (t == 0))
                                  ? 0
                                  : std::numeric_limits<double>::max() / 2;
            auto it = graph[j][t].begin();
            if (!graph[j][t].empty() && t <= Hmax - tmp->processing_time) {
                double reduced_cost =
                    (j == n) ? tmp->weighted_tardiness_start(t)
                             : tmp->weighted_tardiness_start(t) - aux_pi[j];
                forward_F[j][t] =
                    forward_F[(*it)->job][t - (*it)->processing_time] +
                    reduced_cost;

                A[j][t] = (*it);
                B[j][t] = t - (*it)->processing_time;
                it++;
                while (it != graph[j][t].end()) {
                    reduced_cost =
                        (j == n) ? tmp->weighted_tardiness_start(t)
                                 : tmp->weighted_tardiness_start(t) - aux_pi[j];
                    double result =
                        ((*it)->job != vector_jobs[j]->job)
                            ? forward_F[(*it)->job][t - (*it)->processing_time]
                            : forward_F[(*it)->job][t - 1];
                    result += reduced_cost;
                    if (forward_F[j][t] >= result) {
                        forward_F[j][t] = result;
                        A[j][t] = (*it);
                        B[j][t] = ((*it)->job != vector_jobs[j]->job)
                                      ? t - (*it)->processing_time
                                      : t - 1;
                    }
                    it++;
                }
            }
        }
    }
}

void PricerSolverArcTimeDp::forward_evaluator(std::span<const double>& _pi) {
    forward_F[n][0] = 0;

    for (auto t : ranges::views::ints(int{}, static_cast<int>(Hmax) + 1)) {
        for (auto j = 0UL; j <= n; ++j) {
            Job* tmp = vector_jobs[j];
            A[j][t] = nullptr;
            B[j][t] = -1;
            forward_F[j][t] = ((j == n) && (t == 0))
                                  ? 0
                                  : std::numeric_limits<double>::max() / 2;
            auto it = graph[j][t].begin();
            if (!graph[j][t].empty() && t <= Hmax - tmp->processing_time) {
                double reduced_cost =
                    (j == n) ? tmp->weighted_tardiness_start(t)
                             : tmp->weighted_tardiness_start(t) - _pi[j];
                forward_F[j][t] =
                    forward_F[(*it)->job][t - (*it)->processing_time] +
                    reduced_cost;

                A[j][t] = (*it);
                B[j][t] = t - (*it)->processing_time;
                it++;
                while (it != graph[j][t].end()) {
                    reduced_cost =
                        (j == n) ? tmp->weighted_tardiness_start(t)
                                 : tmp->weighted_tardiness_start(t) - _pi[j];
                    double result =
                        ((*it)->job != vector_jobs[j]->job)
                            ? forward_F[(*it)->job][t - (*it)->processing_time]
                            : forward_F[(*it)->job][t - 1];
                    result += reduced_cost;
                    if (forward_F[j][t] >= result) {
                        forward_F[j][t] = result;
                        A[j][t] = (*it);
                        B[j][t] = ((*it)->job != vector_jobs[j]->job)
                                      ? t - (*it)->processing_time
                                      : t - 1;
                    }
                    it++;
                }
            }
        }
    }
}


auto PricerSolverArcTimeDp::pricing_algorithm(std::span<const double>& _pi)
    -> PricingSolution {
    PricingSolution   sol(_pi[n]);
    std::vector<Job*> v;

    forward_evaluator(_pi);

    auto job = n;
    auto T = Hmax;

    while (T > 0) {
        auto aux_job = A[job][T]->job;
        int  aux_T = B[job][T];
        if (aux_job != n) {
            v.push_back(vector_jobs[aux_job]);
            sol.C_max += vector_jobs[aux_job]->processing_time;
            sol.cost += vector_jobs[aux_job]->weighted_tardiness_start(aux_T);
            sol.obj += vector_jobs[aux_job]->weighted_tardiness_start(aux_T) -
                       _pi[aux_job];
        }
        job = aux_job;
        T = aux_T;
    }

    sol.C_max = 0;
    for (auto& it : v | ranges::views::reverse) {
        sol.jobs.push_back(it);
    }

    return sol;
}


auto PricerSolverArcTimeDp::farkas_pricing(
    [[maybe_unused]] std::span<const double>& _pi) -> PricingSolution {
    PricingSolution opt_sol;

    return opt_sol;
}


void PricerSolverArcTimeDp::construct_lp_sol_from_rmp(
    const std::span<const double>&              lambda,
    const std::vector<std::shared_ptr<Column>>& columns) {
    std::fill(lp_x.begin(), lp_x.end(), 0.0);
    auto positive = [](const auto& tmp) { return (tmp.first > 0.0); };

    for (auto&& [x, col] : ranges::views::zip(lambda, columns) |
                               ranges::views::filter(positive)) {
        auto counter = 0UL;
        auto i = n;
        auto t = 0UL;
        while (t < Hmax + 1) {
            Job* tmp_j = nullptr;
            auto j = n;

            if (counter < col->job_list.size()) {
                tmp_j = col->job_list[counter];
                j = tmp_j->job;
            }

            lp_x[i * (n + 1) * (Hmax + 1) + j * (Hmax + 1) + t] += x;

            if (tmp_j == nullptr) {
                i = n;
                t += 1;
            } else {
                i = j;
                t += tmp_j->processing_time;
                counter++;
            }
        }
    }
}

auto PricerSolverArcTimeDp::get_nb_edges() -> size_t {
    auto nb_edges = size_t{};
    for (auto& it : graph) {
        for (auto& it_in : it) {
            nb_edges += it_in.size();
        }
    }
    return nb_edges;
}

auto PricerSolverArcTimeDp::get_nb_vertices() -> size_t {
    size_t nb_vertices = 0U;
    for (auto& it : graph) {
        for (auto& it_in : it) {
            if (!it_in.empty()) {
                nb_vertices++;
            }
        }
    }
    return nb_vertices;
}

auto PricerSolverArcTimeDp::print_num_paths()
    -> boost::multiprecision::cpp_int {
    return 0;
}

auto PricerSolverArcTimeDp::check_column(Column const* col) -> bool {
    // std::span aux_set{set->pdata, set->len};
    const auto& set = col->job_list;
    auto        counter = set.size() - 1;
    auto        i = n;
    auto        t = 0UL;

    while (t < Hmax + 1) {
        Job* tmp_j = nullptr;
        auto j = n;

        if (counter < set.size()) {
            j = set[counter]->job;
        }

        if (std::find(graph[j][t].begin(), graph[j][t].end(), vector_jobs[i]) ==
            graph[j][t].end()) {
            return true;
        }

        if (tmp_j == nullptr) {
            i = n;
            t += 1;
        } else {
            i = j;
            t += tmp_j->processing_time;
            counter--;
        }
    }

    return false;
}

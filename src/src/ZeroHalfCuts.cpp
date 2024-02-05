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

#include "ZeroHalfCuts.hpp"
#include <fmt/format.h>               // for format, print
#include <gurobi_c++.h>               // for GRBModel, GRBVar
#include <algorithm>                  // for __for_each_fn, for_each
#include <array>                      // for array, array<>::valu...
#include <cmath>                      // for floor
#include <functional>                 // for identity
#include <memory>                     // for shared_ptr, unique_ptr
#include <range/v3/view/drop.hpp>     // for drop, drop_fn
#include <range/v3/view/filter.hpp>   // for filter
#include <range/v3/view/iota.hpp>     // for iota_view, iota_view...
#include <range/v3/view/join.hpp>     // for join_view<>::cursor
#include <range/v3/view/reverse.hpp>  // for reverse_view, revers...
#include <range/v3/view/take.hpp>     // for take_view, take, tak...
#include <range/v3/view/view.hpp>     // for operator|, view_closure
#include <range/v3/view/zip.hpp>      // for zip
#include <utility>                    // for move
#include <vector>                     // for vector
#include "ModelInterface.hpp"         // for ConstraintGeneric
#include "ModernDD/NodeBddTable.hpp"  // for NodeTableEntity
#include "ModernDD/NodeId.hpp"        // for NodeId
#include "gurobi_c.h"                 // for GRB_INFINITY, GRB_PR...

ZeroHalfCuts::ZeroHalfCuts(size_t                    _nb_jobs,
                           size_t                    _nb_machines,
                           ReformulationModel*       _rmp_model,
                           NodeId const&             _root,
                           NodeTableEntity<NodeBdd>* _table)
    : env(std::make_unique<GRBEnv>()),
      model(std::make_unique<GRBModel>(*env)),
      nb_jobs(_nb_jobs),
      nb_machines(_nb_machines),
      rmp_model(_rmp_model),
      root(_root),
      table(_table),
      jobs_var(nb_jobs) {
    // Limit how many solutions to collect
    // model->set(GRB_IntParam_PoolSolutions, 124);

    // Limit the search space by setting a gap for the worst possible solution
    // that will be accepted
    // model->set(GRB_DoubleParam_PoolGap, 0.10);

    // do a systematic search for the k-best solutions
    model->set(GRB_IntParam_PoolSearchMode, GRB_PRESOLVE_AGGRESSIVE);
    model->set(GRB_IntParam_Method, GRB_METHOD_AUTO);
    model->set(GRB_IntParam_Threads, 1);
    model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model->set(GRB_IntParam_Presolve, GRB_PRESOLVE_AGGRESSIVE);
    model->set(GRB_DoubleParam_Heuristics, 1 / HALF);
    model->set(GRB_IntParam_MIPFocus, 1);
    model->set(GRB_DoubleParam_TimeLimit, TIMELIMIT);

    generate_model();
}

auto ZeroHalfCuts::add_cuts() -> bool {
    return true;
}

void ZeroHalfCuts::generate_model() {
    try {
        for (auto i : ranges::views::ints(size_t{}, nb_jobs)) {
            jobs_var[i] =
                model->addVar(0.0, 1.0, 0.0, 'B', fmt::format("jobs_{}", i));
        }

        auto& root_node = table->node(root);
        auto& terminal_node = table->node(1);
        terminal_node.update_lp_visited(false);
        root_node.set_sigma(
            model->addVar(0.0, 1.0, 0.0, 'B', fmt::format("sigma_root")));
        node_ids.push_back(root);

        dfs(root);

        GRBLinExpr          expr = 0;
        std::vector<double> coeffs(nb_jobs, 1.0);
        auto                m = nb_machines % 2 == 0 ? 0.0 : 1.0;
        q = model->addVar(0.0, GRB_INFINITY, 0.0, 'I');

        expr.addTerms(coeffs.data(), jobs_var.data(),
                      static_cast<int>(coeffs.size()));
        expr += m * root_node.get_sigma() - m * terminal_node.get_sigma() -
                HALF * q;
        model->addConstr(expr, '=', 1.0);
        model->update();
        init_table();

    } catch (GRBException& e) {
        fmt::print("Error code = {0} {1:-^{2}} {3}", e.getErrorCode(), "",
                   ALIGN, e.getMessage());
    }
}

void ZeroHalfCuts::init_table() {
    for (auto& it : (*table) | ranges::views::join) {
        // it.in_degree = {};
        it.reset_in_degree();
        // it.in_degree_1 = 0;
        it.reset_in_edges();
        it.reset_coeff_cut();
    }

    for (auto& it : *table | ranges::views::take(root.row() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        auto& n0 = table->node(it[0]);
        auto& n1 = table->node(it[1]);

        // n0.in_edges[0].push_back(it.ptr_node_id);
        n0.update_in_degree(false);
        // n0.in_degree[0]++;
        // n1.in_edges[1].push_back(it.ptr_node_id);
        // n1.in_degree[1]++;
        n1.update_in_degree(true);
    }
}

void ZeroHalfCuts::init_coeff_cut() {
    for (auto& it : node_ids) {
        auto& node = table->node(it);
        init_coeff_node(&node);
    }

    for (auto& iter : node_ids_lift) {
        auto& node = table->node(iter);
        init_coeff_node(&node);
    }
    node_ids_lift.clear();
}

void ZeroHalfCuts::init_coeff_node(NodeBdd* node) {
    for (auto k : {false, true}) {
        node->reset_coeff_cut(k);
        // for (auto& it : node->get_in_edges(k)) {
        // auto aux = it.lock();
        // if (aux) {
        //     auto& aux_node = table->node(*aux);
        //     aux_node.reset_coeff_cut(k);
        // }
        // }
    }
}

void ZeroHalfCuts::construct_cut() {
    std::unique_ptr<GenericData> data = std::make_unique<GenericData>();

    auto add_coeff_constr = [&](const auto& it) {
        auto& node = table->node(it);

        for (auto k : {false, true}) {
            auto coeff = floor(node.get_coeff_cut(k) / HALF);
            if (coeff > EPS_CUT) {
                data->add_coeff_hash_table(node.get_nb_job(), node.get_weight(),
                                           k, -coeff);
            }

            for (auto& iter : node.get_in_edges(k)) {
                auto aux = iter.lock();
                if (aux) {
                    auto& aux_node = table->node(*aux);
                    auto  coeff_in = floor(aux_node.get_coeff_cut(k) / 2);
                    if (coeff_in < -EPS_CUT) {
                        data->add_coeff_hash_table(aux_node.get_nb_job(),
                                                   aux_node.get_weight(), k,
                                                   -coeff_in);
                    }
                }
            }
        }
    };

    // auto print_node_ids = [&](const auto& it) {
    //     auto& node = table->node(it);
    //     if (it > 1 && node.sigma.get(GRB_DoubleAttr_Xn) > 0.0) {
    //         fmt::print("node {} {} {} {} {} | ", node.get_nb_job(),
    //                    node.get_weight(), node.sigma.get(GRB_DoubleAttr_Xn),
    //                    node.lp_x[0], node.lp_x[1]);
    //     } else if (node.sigma.get(GRB_DoubleAttr_Xn) > 0.0) {
    //         fmt::print("Terminal Node | ");
    //     }
    // };

    // std::ranges::for_each(node_ids, add_coeff_constr);
    // std::for_each(node_ids_lift.begin(), node_ids_lift.end(),
    // add_coeff_constr); std::for_each(node_ids.begin(), node_ids.end(),
    // print_node_ids); fmt::print("\n");

    auto& root_node = table->node(root);
    auto& terminal_node = table->node(1);
    auto  rhs = 0.0;
    rhs += (root_node.get_sigma().get(GRB_DoubleAttr_Xn) > EPS_CUT)
               ? static_cast<double>(nb_machines)
               : 0.0;
    rhs += (terminal_node.get_sigma().get(GRB_DoubleAttr_Xn) > EPS_CUT)
               ? -static_cast<double>(nb_machines)
               : 0.0;

    for (auto& it : jobs_var) {
        auto x = it.get(GRB_DoubleAttr_Xn);
        if (x > EPS_CUT) {
            rhs += 1.0;
        }
    }
    // fmt::print("test rhs = {}\n", rhs);
    std::shared_ptr<ConstraintGeneric> constr{
        std::make_shared<ConstraintGeneric>(data.release(),
                                            -floor(rhs / HALF))};
    // data->list_coeff();
    // fmt::print("RHS = {}\n",
    //            -floor((2.0 * q.get(GRB_DoubleAttr_Xn) + 1 + rhs) / 2.0));
    for (auto& it : cut_list) {
        if (*constr == *it) {
            return;
        }
    }
    cut_list.push_back(std::move(constr));
}

auto ZeroHalfCuts::get_cut_list()
    -> std::vector<std::shared_ptr<ConstraintGeneric>> {
    return cut_list;
}

void ZeroHalfCuts::generate_cuts() {
    try {
        model->optimize();
        auto status = model->get(GRB_IntAttr_Status);

        if (status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE ||
            status == GRB_UNBOUNDED) {
            fmt::print(
                "The model cannot be solved because it is infeasible or "
                "unbounded\n");
            return;
        }

        if (status != GRB_OPTIMAL) {
            fmt::print("Optimization was stopped with status {}\n", status);
        }

        // Print number of solutions stored
        auto nb_solutions = model->get(GRB_IntAttr_SolCount);
        fmt::print("Number of solutions found: {}\n", nb_solutions);

        for (auto i = 0; i < nb_solutions; i++) {
            model->set(GRB_IntParam_SolutionNumber, i);
            init_coeff_cut();

            bool add = false;
            auto calc_coeff_cut = [&](const auto& it) {
                auto& node = table->node(it);
                auto  x = node.get_sigma().get(GRB_DoubleAttr_Xn);
                if (x > EPS_CUT) {
                    if (it > 1) {
                        // node.coeff_cut[0] += 1.0;
                        // node.coeff_cut[1] += 1.0;
                        for (auto k : {false, true}) {
                            node.update_coeff_cut(k);
                        }
                        // fmt::print("Node {} {}\n", node.get_nb_job(),
                        //    node.get_weight());
                        add = true;
                    }

                    for (auto j : {false, true}) {
                        for (auto& it_aux : node.get_in_edges(j)) {
                            auto aux = it_aux.lock();
                            if (aux) {
                                auto& aux_node = table->node(*aux);
                                aux_node.reduce_coeff_cut(j);
                            }
                        }
                    }
                }
            };

            for (auto k = root.row(); k > 0; k--) {
                for (auto& it : (*table)[k]) {
                    auto j = it.get_nb_job();
                    auto x = jobs_var[j].get(GRB_DoubleAttr_Xn);
                    if (x > EPS_CUT) {
                        it.update_coeff_cut(true);
                    }
                }
            }

            // for (auto i = 0; i < nb_jobs; i++) {
            //     auto x = jobs_var[i].get(GRB_DoubleAttr_Xn);
            //     if (x > 1e-4) {
            //         fmt::print("{} ", i);
            //     }
            // }
            // fmt::print("\n");

            // std::ranges::for_each(node_ids, calc_coeff_cut);

            // for (auto& iter : node_ids) {
            //     auto& node = table->node(iter);
            //     if (node.branch[0] <= 1) {
            //         continue;
            //     }
            //     auto x = node.sigma.get(GRB_DoubleAttr_Xn);
            //     if ((x > 1e-4 && node.lp_x[0] < 1e-6)) {
            //         auto& child = table->node(node.branch[0]);
            //         for (int j = 0; j < 2; j++) {
            //             child.coeff_cut[j] += 1.0;
            //             // if (j) {
            //             //     auto y = jobs_var[child.get_nb_job()].get(
            //             //         GRB_DoubleAttr_Xn);
            //             //     if (y > 1e-4) {
            //             //         child.coeff_cut[j] += 1.0;
            //             //     }
            //             // }
            //             for (auto& it : child.in_edges[j]) {
            //                 auto aux = it.lock();
            //                 if (aux) {
            //                     auto& aux_node = table->node(*aux);
            //                     aux_node.coeff_cut[j] -= 1.0;
            //                     // auto y =
            //                     jobs_var[aux_node.get_nb_job()].get(
            //                     //     GRB_DoubleAttr_Xn);
            //                     // if (y == 1.0 && !aux_node.lp_visited) {
            //                     //     aux_node.coeff_cut[j] += 1.0;
            //                     // }
            //                 }
            //             }
            //         }
            //         node_ids_lift.push_back(node.branch[0]);
            //         dfs_lift(node.branch[0]);
            //     }
            // }
            if (add) {
                construct_cut();
            }
            // fmt::print("---------------------------------------------\n");
        }
        // getchar();

    } catch (GRBException& e) {
        fmt::print("Error code = {0} {1:-^{2}} {3}", e.getErrorCode(), "",
                   ALIGN, e.getMessage());
    }
}

void ZeroHalfCuts::dfs(const NodeId& v) {
    auto& node = table->node(v);
    node.update_lp_visited(true);
    std::array<bool, 2> edges_type{false, true};

    for (auto [child_id, x, y, r, high] :
         ranges::views::zip(node, node.get_lp_x(), node.get_y(), node.get_r(),
                            edges_type) |
             ranges::views::filter(
                 [](const auto& tmp) { return tmp > EPS_CUT; },
                 [](const auto& tmp) { return std::get<1>(tmp); })) {
        // if (x > EPS_CUT) {
        auto& child = table->node(child_id);
        if (!child.get_lp_visited()) {
            if (child_id == 1) {
                child.set_sigma(model->addVar(0.0, 1.0, 0.0, 'B',
                                              fmt::format("sigma_terminal")));
            } else {
                child.set_sigma(
                    model->addVar(0.0, 1.0, 0.0, 'B',
                                  fmt::format("sigma_{}_{}", child.get_nb_job(),
                                              child.get_weight())));
            }

            node_ids.push_back(child_id);
            auto& s_source = node.get_sigma();
            auto& s_head = child.get_sigma();
            auto  str_y =
                fmt::format("y_{}_{}", node.get_nb_job(), node.get_weight());
            auto str_r =
                fmt::format("r_{}_{}", node.get_nb_job(), node.get_weight());
            y = model->addVar(0.0, GRB_INFINITY, x, 'B', str_y);
            r = model->addVar(0.0, GRB_INFINITY, 0.0, 'I', str_r);
            GRBLinExpr expr = -s_head + s_source - y - HALF * r;
            if (high) {
                expr += jobs_var[node.get_nb_job()];
            }
            model->addConstr(expr, '=', 0.0);
            dfs(child_id);
        } else {
            auto& s_source = node.get_sigma();
            auto& s_head = child.get_sigma();
            auto  str_y =
                fmt::format("y_{}_{}", node.get_nb_job(), node.get_weight());
            auto str_r =
                fmt::format("r_{}_{}", node.get_nb_job(), node.get_weight());
            y = model->addVar(0.0, GRB_INFINITY, x, 'B', str_y);
            r = model->addVar(0.0, GRB_INFINITY, 0.0, 'I', str_r);
            GRBLinExpr expr = -s_head + s_source - y - HALF * r;
            if (high) {
                expr += jobs_var[node.get_nb_job()];
            }

            model->addConstr(expr, '=', 0.0);
        }
        // }
    }
}

void ZeroHalfCuts::dfs_lift(const NodeId& v) {
    auto& node = table->node(v);

    for (auto [child, x] : ranges::views::zip(node, node.get_lp_x())) {
        if (child <= 1) {
            continue;
        }
        auto& child_node = table->node(child);

        if (x == 0.0 && !child_node.get_lp_visited()) {
            // child_node.coeff_cut[0] += 1.0;
            // child_node.coeff_cut[1] += 1.0;

            for (auto j : {false, true}) {
                child_node.update_coeff_cut(j);
            }

            for (auto j : {false, true}) {
                for (auto& it : child_node.get_in_edges(j)) {
                    auto aux = it;
                    if (aux) {
                        // auto& aux_node = table->node(*aux);
                        // aux_node.reduce_coeff_cut(j);
                    }
                }
            }
            node_ids_lift.push_back(child);
            dfs_lift(child);
        }
    }
}

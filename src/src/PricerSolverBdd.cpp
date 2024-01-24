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

#include "PricerSolverBdd.hpp"
#include <fmt/core.h>                               // for print
#include <gurobi_c++.h>                             // for GRBLinExpr
#include <algorithm>                                // for find, min, all_of
#include <array>                                    // for array, array<>...
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <boost/graph/adjacency_list.hpp>           // for source, target
#include <boost/graph/detail/adjacency_list.hpp>    // for undirected_edg...
#include <boost/graph/detail/edge.hpp>              // for edge_desc_impl
#include <boost/graph/graphviz.hpp>                 // for write_graphviz
#include <boost/iterator/iterator_facade.hpp>       // for operator!=
#include <boost/multiprecision/cpp_int.hpp>         // for big_cpp
#include <cassert>                                  // for assert
#include <cmath>                                    // for fabs
#include <cstddef>                                  // for size_t
#include <limits>                                   // for numeric_limits
#include <list>                                     // for operator==, list
#include <memory>                                   // for allocator, mak...
#include <ostream>                                  // for operator<<
#include <range/v3/all.hpp>                         // all range-v3 include
#include <set>                                      // for operator==
#include <span>                                     // for span
#include <string>                                   // for char_traits
#include <tuple>                                    // for tuple, get
#include <vector>                                   // for vector, _Bit_r...
#include "CardinalityPaths.hpp"                     // for CardinalityPaths
#include "CglGomory.hpp"                            // for CglGomory
#include "CglZeroHalf.hpp"                          // for CglZeroHalf
#include "Column.h"                                 // for ScheduleSet
#include "DebugLvl.hpp"                             // for DebugLvl
#include "Instance.h"                               // for Instance
#include "Job.h"                                    // for Job, value_dif...
#include "MipGraph.hpp"                             // for MipGraph, Colo...
#include "ModelInterface.hpp"                       // for BddCoeff, Refo...
#include "ModernDD/NodeId.hpp"                      // for NodeId
#include "NodeBdd.hpp"                              // for NodeBdd
#include "NodeData.h"                               // for NodeData
#include "OsiGrbSolverInterface.hpp"                // for OsiGrbSolverInterface
#include "PricerConstruct.hpp"                      // for PricerConstruct
#include "PricerSolverBase.hpp"                     // for PricerSolverBa...
#include "PricingSolution.hpp"                      // for PricingSolution
#include "Solution.hpp"                             // for Machine, VecJo...
#include "gurobi_c.h"                               // for GRB_EQUAL, GRB...
#include "or-utils/lp.h"                            // for lp_interface_g...
#include "or-utils/util.h"

PricerSolverBdd::PricerSolverBdd(const Instance& instance)
    : PricerSolverBase(instance),
      decision_diagram(PricerConstruct(instance)),
      size_graph{decision_diagram.size()},
      ordered_jobs_new(instance.vector_pair),
      original_model(reformulation_model),
      H_min{static_cast<size_t>(instance.H_min)},
      H_max(instance.H_max) {
    auto& table = *(decision_diagram.getDiagram());

    auto i = decision_diagram.root().row();
    ordered_jobs_new |=
        ranges::actions::remove_if([&]([[maybe_unused]] const auto& tmp) {
            bool remove = ranges::all_of(
                table[i], [&](const auto& n) { return n[1] == 0; });
            --i;
            return remove;
        });

    if (debug_lvl(0)) {
        fmt::print("{0: <{2}}{1}\n", "The new number of layers",
                   ordered_jobs_new.size(), ALIGN);
    }

    decision_diagram.compressBdd();
    init_table();
    cleanup_arcs();
    bottom_up_filtering();
    topdown_filtering();
    construct_mipgraph();
    init_coeff_constraints();
    if (debug_lvl(0)) {
        fmt::print("ENDING CONSTRUCTION\n\n");
    }
}

PricerSolverBdd::PricerSolverBdd(const PricerSolverBdd& src)
    : PricerSolverBase(src),
      decision_diagram(src.decision_diagram),
      size_graph(src.size_graph),
      nb_removed_edges(src.nb_removed_edges),
      nb_removed_nodes(src.nb_removed_nodes),
      ordered_jobs_new(src.ordered_jobs_new),
      mip_graph(src.mip_graph),
      original_model(src.original_model),
      H_min(src.H_min),
      H_max(src.H_max) {
    cleanup_arcs();
    bottom_up_filtering();
    topdown_filtering();
    construct_mipgraph();
}

PricerSolverBdd::~PricerSolverBdd() = default;

void PricerSolverBdd::construct_mipgraph() {
    mip_graph.clear();
    auto& table = *(decision_diagram.getDiagram());

    auto index{0U};
    for (auto&& [i, row] :
         table | ranges::views::take(decision_diagram.root().row() + 1) |
             ranges::views::enumerate | ranges::views::reverse) {
        for (auto&& [j, it] : row | ranges::views::enumerate) {
            if (NodeId(i, j) != 0 && (it.any_of_calc())) {
                it.set_key(
                    boost::add_vertex({index++, NodeId(i, j)}, mip_graph));
            }
        }
    }

    auto edge_index{0U};

    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        if (it[0] != 0 && it.get_calc(false)) {
            auto& n0 = table.node(it[0]);
            add_edge(it.get_key(), n0.get_key(),
                     {edge_index++, false, GRBVar()}, mip_graph);
        }

        if (it[1] != 0 && it.get_calc(true)) {
            auto& n1 = table.node(it[1]);
            add_edge(it.get_key(), n1.get_key(), {edge_index++, true, GRBVar()},
                     mip_graph);
        }
    }
}

void PricerSolverBdd::init_coeff_constraints() {
    auto& table = *(decision_diagram.getDiagram());
    original_model.clear_all_coeff();

    for (auto i = decision_diagram.root().row(); i > 0UL; i--) {
        for (auto& it : table[i]) {
            it.add_coeff_list_clear();
            for (const auto&& [c, constr] :
                 original_model | ranges::views::take(convex_constr_id) |
                     ranges::views::enumerate) {
                VariableKeyBase key_aux(it.get_nb_job(), it.get_weight());
                auto            coeff = (*constr.get_constr())(key_aux);
                if (fabs(coeff) > EPS_SOLVER) {
                    auto ptr_coeff{std::make_shared<BddCoeff>(
                        it.get_nb_job(), it.get_weight(), coeff, 0.0, c)};
                    constr.add_coeff_to_list(ptr_coeff);
                    it.add_coeff_list(ptr_coeff, true);
                }
            }
        }
    }

    auto&           root = decision_diagram.root();
    auto&           root_node = table.node(root);
    VariableKeyBase key_aux(root_node.get_nb_job(), root_node.get_weight(),
                            true);
    auto*           constr = reformulation_model[convex_constr_id].get();
    auto            coeff = (*constr)(key_aux);
    if (fabs(coeff) > EPS_SOLVER) {
        auto ptr_coeff_high{std::make_shared<BddCoeff>(
            root_node.get_nb_job(), root_node.get_weight(), coeff, true, true)};
        original_model[convex_constr_id].add_coeff_to_list(ptr_coeff_high);
        auto ptr_coeff_low{std::make_shared<BddCoeff>(root_node.get_nb_job(),
                                                      root_node.get_weight(),
                                                      coeff, false, true)};
        original_model[convex_constr_id].add_coeff_to_list(ptr_coeff_low);
    }
}

void PricerSolverBdd::update_coeff_constraints() {
    auto nb_constr = original_model.get_nb_constraints();

    for (auto& j : reformulation_model | ranges::views::drop(nb_constr)) {
        original_model.add_constraint(j);
    }

    auto& table = *(decision_diagram.getDiagram());
    for (auto i = decision_diagram.root().row(); i > 0; i--) {
        for (auto& it : table[i]) {
            for (auto&& [c, constr] : original_model |
                                          ranges::views::drop(nb_constr) |
                                          ranges::views::enumerate) {
                VariableKeyBase key_high{it.get_nb_job(), it.get_weight(),
                                         true};

                auto coeff_high = (*constr.get_constr())(key_high);
                if (fabs(coeff_high) > EPS_SOLVER) {
                    auto ptr_coeff{std::make_shared<BddCoeff>(
                        it.get_nb_job(), it.get_weight(), coeff_high, 0.0,
                        nb_constr + c)};
                    constr.add_coeff_to_list(ptr_coeff);
                    it.add_coeff_list(ptr_coeff, true);
                }

                BddCoeff key_low(it.get_nb_job(), it.get_weight(), 0.0, 0.0,
                                 nb_constr + c, false);
                auto     coeff_low = (*constr.get_constr())(key_low);
                if (fabs(coeff_low) > EPS_SOLVER) {
                    auto ptr_coeff{std::make_shared<BddCoeff>(
                        it.get_nb_job(), it.get_weight(), coeff_low, 0.0,
                        nb_constr + c, false)};
                    constr.add_coeff_to_list(ptr_coeff);
                    it.add_coeff_list(ptr_coeff, false);
                }
            }
        }
    }
}

void PricerSolverBdd::init_table() {
    auto& table = *(decision_diagram.getDiagram());
    /** init table */
    auto& root = table.node(decision_diagram.root());
    root.init_node(0UL, true);
    root.set_node_id_label(decision_diagram.root());
    root.reset_all(convex_constr_id);

    for (auto i :
         ranges::views::ints(size_t{}, decision_diagram.root().row() + 1) |
             ranges::views::reverse) {
        for (auto it = 0UL; it < table[i].size(); it++) {
            if (i != 0UL) {
                auto  layer = ordered_jobs_new.size() - i;
                auto& tmp_pair = ordered_jobs_new[layer];
                auto  node_id = NodeId(i, it);
                auto& node = table.node(node_id);
                auto* aux_job = tmp_pair.first;
                auto  w = node.get_weight();
                auto  p = static_cast<size_t>(aux_job->processing_time);

                auto& n0 = table.node(node[0]);
                auto& n1 = table.node(node[1]);
                n0.set_node_id_label(node[0]);
                n1.set_node_id_label(node[1]);
                node.set_job_label(aux_job);

                node.set_ptr_node_id(i, it);
                node.set_job(aux_job);
                n0.init_node(w);
                n1.init_node(w + p);
                node.set_cost(
                    static_cast<double>(aux_job->weighted_tardiness_start(w)));

                n0.update_in_degree(false);
                n1.update_in_degree(true);

            } else {
                auto& node = table.node(NodeId(i, it));
                node.set_job(nullptr);
            }
        }
    }
}

void PricerSolverBdd::insert_constraints_lp(NodeData* pd) {
    auto nb_rows = pd->osi_rmp->getNumRows();
    auto nb_new_constraints = reformulation_model.size() -
                              static_cast<size_t>(pd->osi_rmp->getNumRows());

    fmt::print("nb rows initial {} {} {}\n", nb_rows,
               reformulation_model.size(), nb_new_constraints);

    assert((nb_new_constraints <=
            (pd->id_pseudo_schedules - pd->id_next_var_cuts)));
    std::vector<int>    starts(nb_new_constraints + 1);
    std::vector<char>   sense(nb_new_constraints);
    std::vector<double> rhs(nb_new_constraints);
    std::vector<double> rhs_ub(nb_new_constraints, pd->osi_rmp->getInfinity());
    std::vector<int>    column_ind;
    std::vector<double> coeff;

    int pos = 0;
    for (auto&& [c, constr] : reformulation_model |
                                  ranges::views::drop(nb_rows) |
                                  ranges::views::enumerate) {
        sense[c] = constr->get_sense();
        starts[c] = pos;
        rhs[c] = constr->get_rhs();

        if (rhs[c] != 0.0) {
            pos++;
            column_ind.push_back(pd->id_next_var_cuts++);
            if (sense[c] == '>') {
                coeff.push_back(1.0);
            } else {
                coeff.push_back(-1.0);
            }
        }

        for (auto&& [i, it] : pd->localColPool | ranges::views::enumerate) {
            auto& table = *(decision_diagram.getDiagram());
            auto  tmp_nodeid(decision_diagram.root());

            auto coeff_val = 0.0;
            auto job_it = it->job_list.begin();
            while (tmp_nodeid > 1) {
                auto& tmp_node = table.node(tmp_nodeid);
                auto* tmp_j =
                    (job_it != it->job_list.end()) ? *job_it : nullptr;

                VariableKeyBase key(tmp_node.get_nb_job(),
                                    tmp_node.get_weight(),
                                    tmp_j == tmp_node.get_job());
                coeff_val += (*constr)(key);
                if (key.get_high()) {
                    tmp_nodeid = tmp_node[1];
                    ++job_it;
                } else {
                    tmp_nodeid = tmp_node[0];
                }
            }

            assert((tmp_nodeid == 1));

            if (fabs(coeff_val) > EPS_SOLVER) {
                column_ind.push_back(pd->id_pseudo_schedules +
                                     static_cast<int>(i));
                coeff.push_back(coeff_val);
                pos++;
            }
        }
    }

    starts[nb_new_constraints] = pos;

    pd->osi_rmp->addRows(nb_new_constraints, starts.data(), column_ind.data(),
                         coeff.data(), nullptr, rhs.data());
    pd->id_row.resize(reformulation_model.size(), 0);
    pd->coeff_row.resize(reformulation_model.size(), 0.0);
}


auto PricerSolverBdd::compute_reduced_cost(const PricingSolution&   sol,
                                           std::span<const double>& pi,
                                           double* lhs) -> double {
    double    result = sol.cost;
    auto&     table = *decision_diagram.getDiagram();
    auto      tmp_nodeid(decision_diagram.root());
    auto      it = sol.jobs.begin();
    std::span aux_lhs{lhs, reformulation_model.size()};
    ranges::fill(aux_lhs, 0.0);

    while (tmp_nodeid > 1) {
        auto& tmp_node = table.node(tmp_nodeid);
        auto* tmp_j = (it != sol.jobs.end()) ? *it : nullptr;

        VariableKeyBase key(tmp_node.get_nb_job(), tmp_node.get_weight(),
                            tmp_j == tmp_node.get_job());
        if (key.get_high()) {
            tmp_nodeid = tmp_node[1];
            ++it;
            auto* constr = reformulation_model[key.get_j()].get();
            auto  dual = pi[key.get_j()];
            auto  coeff = (*constr)(key);

            if (fabs(coeff) > EPS_SOLVER) {
                result -= coeff * dual;
                aux_lhs[key.get_j()] += coeff;
            }
        } else {
            tmp_nodeid = tmp_node[0];
        }

        for (auto c = convex_constr_id + 1; c < reformulation_model.size();
             ++c) {
            auto* constr = reformulation_model[c].get();
            auto  dual = pi[c];
            auto  coeff = (*constr)(key);

            if (fabs(coeff) > EPS_SOLVER) {
                result -= coeff * dual;
                aux_lhs[c] += coeff;
            }
        }
    }

    auto*           constr = reformulation_model[convex_constr_id].get();
    auto            dual = pi[convex_constr_id];
    VariableKeyBase k(0, 0, true);
    auto            coeff = (*constr)(k);
    result -= coeff * dual;
    aux_lhs[convex_constr_id] += coeff;

    return result;
}

void PricerSolverBdd::compute_lhs(const PricingSolution& sol, double* lhs) {
    auto&     table = *decision_diagram.getDiagram();
    auto      tmp_nodeid(decision_diagram.root());
    auto      it = sol.jobs.begin();
    std::span aux_lhs{lhs, reformulation_model.size()};
    ranges::fill(aux_lhs, 0.0);

    while (tmp_nodeid > 1) {
        auto& tmp_node = table.node(tmp_nodeid);
        auto* tmp_j = (it != sol.jobs.end()) ? *it : nullptr;

        VariableKeyBase key(tmp_node.get_nb_job(), tmp_node.get_weight(),
                            tmp_j == tmp_node.get_job());
        if (key.get_high()) {
            tmp_nodeid = tmp_node[1];
            ++it;
            auto* constr = reformulation_model[key.get_j()].get();
            auto  coeff = (*constr)(key);

            if (fabs(coeff) > EPS_SOLVER) {
                aux_lhs[key.get_j()] += coeff;
            }
        } else {
            tmp_nodeid = tmp_node[0];
        }

        for (auto c = convex_constr_id + 1; c < reformulation_model.size();
             ++c) {
            auto* constr = reformulation_model[c].get();
            auto  coeff = (*constr)(key);

            if (fabs(coeff) > EPS_SOLVER) {
                aux_lhs[c] += coeff;
            }
        }
    }

    auto*           constr = reformulation_model[convex_constr_id].get();
    VariableKeyBase k(0, 0, true);
    auto            coeff = (*constr)(k);
    aux_lhs[convex_constr_id] += coeff;
}

void PricerSolverBdd::compute_lhs(const Column& sol, double* lhs) {
    auto&     table = *decision_diagram.getDiagram();
    auto      tmp_nodeid(decision_diagram.root());
    auto      it = sol.job_list.begin();
    std::span aux_lhs{lhs, reformulation_model.size()};
    ranges::fill(aux_lhs, 0.0);

    while (tmp_nodeid > 1) {
        auto& tmp_node = table.node(tmp_nodeid);
        auto* tmp_j = (it != sol.job_list.end()) ? *it : nullptr;

        VariableKeyBase key(tmp_node.get_nb_job(), tmp_node.get_weight(),
                            tmp_j == tmp_node.get_job());
        if (key.get_high()) {
            tmp_nodeid = tmp_node[1];
            ++it;
            auto* constr = reformulation_model[key.get_j()].get();
            auto  coeff = (*constr)(key);

            if (fabs(coeff) > EPS_SOLVER) {
                aux_lhs[key.get_j()] += coeff;
            }
        } else {
            tmp_nodeid = tmp_node[0];
        }

        for (auto c = convex_constr_id + 1; c < reformulation_model.size();
             ++c) {
            auto* constr = reformulation_model[c].get();
            auto  coeff = (*constr)(key);

            if (fabs(coeff) > EPS_SOLVER) {
                aux_lhs[c] += coeff;
            }
        }
    }

    auto*           constr = reformulation_model[convex_constr_id].get();
    VariableKeyBase k(0, 0, true);
    auto            coeff = (*constr)(k);
    aux_lhs[convex_constr_id] += coeff;
}

auto PricerSolverBdd::compute_subgradient(const PricingSolution& sol,
                                          double* sub_gradient) -> double {
    double    result = sol.cost;
    auto&     table = *decision_diagram.getDiagram();
    auto      tmp_nodeid(decision_diagram.root());
    auto      it = sol.jobs.begin();
    auto      nb_constraints = reformulation_model.size();
    auto      rhs = -reformulation_model[convex_constr_id]->get_rhs();
    std::span aux_subgradient{sub_gradient, nb_constraints};

    for (auto&& [constr, subgradient] :
         ranges::views::zip(reformulation_model, aux_subgradient)) {
        subgradient = constr->get_rhs();
    }

    while (tmp_nodeid > 1) {
        auto& tmp_node = table.node(tmp_nodeid);
        auto* tmp_j = it != sol.jobs.end() ? *it : nullptr;

        VariableKeyBase key(tmp_node.get_nb_job(), tmp_node.get_weight(),
                            tmp_j == tmp_node.get_job());
        if (key.get_high()) {
            tmp_nodeid = tmp_node[1];
            ++it;
            auto* constr = reformulation_model[key.get_j()].get();
            auto  coeff = (*constr)(key);

            if (fabs(coeff) > EPS_SOLVER) {
                aux_subgradient[key.get_j()] -= coeff * rhs;
            }
        } else {
            tmp_nodeid = tmp_node[0];
        }

        for (auto c = convex_constr_id + 1; c < reformulation_model.size();
             c++) {
            auto* constr = reformulation_model[c].get();
            auto  coeff = (*constr)(key);

            if (fabs(coeff) > EPS_SOLVER) {
                aux_subgradient[c] -= coeff * rhs;
            }
        }
    }

    aux_subgradient[convex_constr_id] += rhs;
    assert(aux_subgradient[convex_constr_id] == 0.0);
    assert(tmp_nodeid == 1);

    return result;
}

auto PricerSolverBdd::refinement_structure(
    const std::vector<std::shared_ptr<Column>>& paths) -> bool {
    bool  refined_structure = false;
    auto& table = *decision_diagram.getDiagram();
    for (const auto& path : paths) {
        std::vector<NodeId> P;
        std::vector<bool>   L;
        std::vector<NodeId> S;
        std::vector<size_t> index_P;

        S.resize(convex_constr_id, NodeId());
        index_P.resize(convex_constr_id, 0UL);

        P.push_back(decision_diagram.root());
        auto job_it = path->job_list.begin();
        Job* conflict_job{nullptr};

        while (P.back() > 1) {
            auto& tmp_node{table.node(P.back())};
            auto* tmp_job{job_it != path->job_list.end() ? *job_it : nullptr};

            if (tmp_job == tmp_node.get_job() && tmp_job != nullptr) {
                if (S[tmp_job->job] != 0) {
                    conflict_job = tmp_job;
                    break;
                }
                S[tmp_job->job] = P.back();
                index_P[tmp_job->job] = P.size() - 1;
                P.push_back(tmp_node[1]);
                L.push_back(true);
                ++job_it;
            } else {
                P.push_back(tmp_node[0]);
                L.push_back(false);
            }
        }

        if (conflict_job != nullptr) {
            auto nodeid_new = P[index_P[conflict_job->job]];
            for (auto&& label :
                 L | ranges::views::drop(index_P[conflict_job->job])) {
                auto& p = table.node(nodeid_new);
                auto& c = label ? table.node(p[1]) : table.node(p[0]);
                auto  w = c;
                nodeid_new = label
                                 ? NodeId(p[1].row(), table[p[1].row()].size())
                                 : NodeId(p[0].row(), table[p[0].row()].size());

                if (conflict_job == c.get_job()) {
                    w[0] = c[0];
                    w[1] = 0;
                }

                w.set_node_id_label(nodeid_new);
                w.reset_all(convex_constr_id);
                table[label ? p[1].row() : p[0].row()].emplace_back(w);

                if(label) {
                    p[1] = nodeid_new;
                } else {
                    p[0] = nodeid_new;
                }
            }
            refined_structure = true;
        }
    }

    if (refined_structure) {
        decision_diagram.compressBdd();
        nb_removed_nodes -= size_graph;
        size_graph = decision_diagram.size();
        bottom_up_filtering();
        topdown_filtering();
        cleanup_arcs();
        construct_mipgraph();
    }

    return refined_structure;
}

void PricerSolverBdd::enumerate_columns() {
    auto& table = *decision_diagram.getDiagram();
    auto  cursor = 0;
    auto  begin = true;
    std::vector<std::tuple<NodeId, bool, boost::dynamic_bitset<>>> path{};
    auto iterations = boost::multiprecision::cpp_int{};
    auto empty = boost::dynamic_bitset<>{convex_constr_id, 0};

    do {
        auto f = begin ? decision_diagram.root() : NodeId(0, 0);
        for (;;) {
            begin = false;
            while (f > 1) { /* down */
                auto const& s = table[f.row()][f.col()];
                auto& set = !path.empty() ? std::get<2>(path.back()) : empty;

                if (s[0] != 0) {
                    cursor = static_cast<int>(path.size());
                    path.emplace_back(f, false, set);
                    f = s[0];
                } else if (!set[s.get_nb_job()]) {
                    path.emplace_back(f, true, set);
                    std::get<2>(path.back())[s.get_nb_job()] = true;
                    f = s[1];
                } else {
                    f = 0;
                }
            }

            if (f == 1) {
                ++iterations;
                break; /* found */
            }

            for (; cursor >= 0; --cursor) { /* up */
                auto&       sel = path[cursor];
                auto const& ss =
                    table[std::get<0>(sel).row()][std::get<0>(sel).col()];
                if (!std::get<1>(sel) && ss[1] != 0 &&
                    !std::get<2>(sel)[ss.get_nb_job()]) {
                    f = std::get<0>(sel);
                    std::get<1>(sel) = true;
                    path.resize(cursor + 1);
                    std::get<2>(path.back()) = std::get<2>(sel);
                    std::get<2>(path.back())[ss.get_nb_job()] = true;
                    f = decision_diagram.child(f, 1);
                    break;
                }
            }

            if (cursor < 0) { /* end() state */
                fmt::print("number of elementary paths {}\n", iterations.str());
                return;
            }
        }
    } while (true);
}


void PricerSolverBdd::enumerate_columns(std::span<const double>& _pi) {
    auto& table = *get_decision_diagram().getDiagram();
    auto  root_id = get_decision_diagram().root();
    compute_labels(_pi);
    auto reduced_cost = table.node(root_id).backward_label[0].get_f();

    auto cursor = 0;
    auto begin = true;

    std::vector<std::tuple<NodeId, bool, boost::dynamic_bitset<>>> path{};
    auto iterations = boost::multiprecision::cpp_int{};
    auto empty = boost::dynamic_bitset<>{convex_constr_id, 0};
    auto aux_nb_machines = static_cast<double>(convex_rhs - 1);

    do {
        auto f = begin ? root_id : NodeId(0, 0);
        for (;;) {
            begin = false;
            while (f > 1) { /* down */
                auto& s = table[f.row()][f.col()];
                auto& set = !path.empty() ? std::get<2>(path.back()) : empty;
                auto  rc = (constLB + aux_nb_machines * reduced_cost +
                           evaluate_rc_arc(s)) < UB;

                if (s[0] != 0) {
                    cursor = static_cast<int>(path.size());
                    path.emplace_back(f, false, set);
                    f = s[0];
                } else if (!set[s.get_nb_job()] && rc) {
                    path.emplace_back(f, true, set);
                    std::get<2>(path.back())[s.get_nb_job()] = true;
                    f = s[1];
                } else {
                    f = 0;
                }
            }

            if (f == 1) {
                ++iterations;
                break; /* found */
            }

            for (; cursor >= 0; --cursor) { /* up */
                auto& sel = path[cursor];
                auto& ss =
                    table[std::get<0>(sel).row()][std::get<0>(sel).col()];
                auto rc = (constLB + aux_nb_machines * reduced_cost +
                           evaluate_rc_arc(ss)) < UB;
                if (!std::get<1>(sel) && ss[1] != 0 &&
                    !std::get<2>(sel)[ss.get_nb_job()] && rc) {
                    f = std::get<0>(sel);
                    std::get<1>(sel) = true;
                    path.resize(cursor + 1);
                    std::get<2>(path.back()) = std::get<2>(sel);
                    std::get<2>(path.back())[ss.get_nb_job()] = true;
                    f = get_decision_diagram().child(f, 1);
                    break;
                }
            }

            if (cursor < 0) { /* end() state */
                fmt::print("number of elementary paths with rc {}\n",
                           iterations.str());
                return;
            }
        }
    } while (true);
}

auto PricerSolverBdd::compute_lagrange(const PricingSolution&     sol,
                                       const std::vector<double>& pi)
    -> double {
    double result = sol.cost;
    auto   dual_bound = 0.0;

    auto& table = *decision_diagram.getDiagram();
    auto  tmp_nodeid(decision_diagram.root());

    auto it = sol.jobs.begin();
    while (tmp_nodeid > 1) {
        auto&           tmp_node = table.node(tmp_nodeid);
        auto*           tmp_j = it != sol.jobs.end() ? *it : nullptr;
        VariableKeyBase key(tmp_node.get_nb_job(), tmp_node.get_weight(),
                            tmp_j == tmp_node.get_job());
        if (key.get_high()) {
            auto dual = pi[key.get_j()];
            auto coeff = (*reformulation_model[key.get_j()])(key);

            if (fabs(coeff) > EPS_SOLVER) {
                result -= coeff * dual;
            }

            ++it;
            tmp_nodeid = tmp_node[1];
        } else {
            tmp_nodeid = tmp_node[0];
        }

        for (auto&& [constr, pi_tmp] :
             ranges::views::zip(reformulation_model, pi) |
                 ranges::views::drop(convex_constr_id + 1)) {
            auto coeff = (*constr)(key);

            if (fabs(coeff) > EPS_SOLVER) {
                result -= coeff * pi_tmp;
            }
        }
    }

    result = std::min(0.0, result);

    for (const auto&& [constr, pi_aux] :
         ranges::views::zip(reformulation_model, pi)) {
        if (constr == reformulation_model[convex_constr_id]) {
            continue;
        }

        dual_bound += constr->get_rhs() * pi_aux;
    }

    result = -reformulation_model[convex_constr_id]->get_rhs() * result;
    result = dual_bound + result;

    return result;
}

auto PricerSolverBdd::compute_lagrange(const PricingSolution&         sol,
                                       const std::span<const double>& pi)
    -> double {
    double result = sol.cost;
    auto   dual_bound = 0.0;

    auto& table = *decision_diagram.getDiagram();
    auto  tmp_nodeid(decision_diagram.root());

    auto it = sol.jobs.begin();
    while (tmp_nodeid > 1) {
        auto&           tmp_node = table.node(tmp_nodeid);
        auto*           tmp_j = it != sol.jobs.end() ? *it : nullptr;
        VariableKeyBase key(tmp_node.get_nb_job(), tmp_node.get_weight(),
                            tmp_j == tmp_node.get_job());
        if (key.get_high()) {
            auto dual = pi[key.get_j()];
            auto coeff = (*reformulation_model[key.get_j()])(key);

            if (fabs(coeff) > EPS_SOLVER) {
                result -= coeff * dual;
            }

            ++it;
            tmp_nodeid = tmp_node[1];
        } else {
            tmp_nodeid = tmp_node[0];
        }

        for (auto&& [constr, pi_tmp] :
             ranges::views::zip(reformulation_model, pi) |
                 ranges::views::drop(convex_constr_id + 1)) {
            auto coeff = (*constr)(key);

            if (fabs(coeff) > EPS_SOLVER) {
                result -= coeff * pi_tmp;
            }
        }
    }

    result = std::min(0.0, result);

    for (const auto&& [constr, pi_aux] :
         ranges::views::zip(reformulation_model, pi)) {
        if (constr == reformulation_model[convex_constr_id]) {
            continue;
        }

        dual_bound += constr->get_rhs() * pi_aux;
    }

    result = -reformulation_model[convex_constr_id]->get_rhs() * result;
    result = dual_bound + result;

    return result;
}

[[maybe_unused]] void PricerSolverBdd::remove_layers_init() {
    auto& table = *(decision_diagram.getDiagram());

    auto i = decision_diagram.root().row();
    ordered_jobs_new |=
        ranges::actions::remove_if([&]([[maybe_unused]] const auto& tmp) {
            bool remove = std::ranges::all_of(
                table[i], [&](const auto& n) { return n[1] == 0; });
            --i;
            return remove;
        });

    if (debug_lvl(0)) {
        fmt::print("{0: <{2}}{1}\n", "The new number of layers",
                   ordered_jobs_new.size(), ALIGN);
    }
}

void PricerSolverBdd::remove_layers() {
    // auto& table = *(decision_diagram.getDiagram());

    // auto i = decision_diagram.topLevel();

    // ordered_jobs_new |= ranges::actions::remove_if([&](const auto& tmp) {
    //     auto remove = true;

    //     for (auto& iter : table[i]) {
    //         if (iter.calc[1]) {
    //             remove = false;
    //         } else {
    //             auto& cur_node_1 = iter[1];
    //             cur_node_1 = 0;
    //         }
    //     }
    //     --i;
    //     return remove;
    // });

    // if (debug_lvl(0)) {
    //     fmt::print("{0: <{2}}{1}\n", "The new number of layers",
    //                ordered_jobs_new.size(), ALIGN);
    // }
}

void PricerSolverBdd::remove_edges() {
    auto& table = *(decision_diagram.getDiagram());

    /** remove the unnecessary nodes of the bdd */
    for (auto& iter : table | ranges::views::drop(1) | ranges::views::join) {
        if (!iter.get_calc(true)) {
            auto& cur_node_1 = iter[1];
            iter.reset_ptr_node_id();
            cur_node_1 = 0;
        }

        if (!iter.get_calc(false)) {
            auto& cur_node_0 = iter[0];
            cur_node_0 = 0;
        }
    }

    decision_diagram.compressBdd();
    nb_removed_nodes -= size_graph;
    size_graph = decision_diagram.size();
}

[[maybe_unused]] void PricerSolverBdd::print_representation_file() {
    auto& table = *(decision_diagram.getDiagram());
    auto  index_edge{std::vector<std::vector<size_t>>(convex_constr_id,
                                                     std::vector<size_t>())};

    auto outfile_file_mip_str =
        problem_name + "_" + std::to_string(convex_rhs) + ".txt";
    std::ofstream out_file_mip(outfile_file_mip_str);

    out_file_mip << boost::num_vertices(mip_graph) << " "
                 << boost::num_edges(mip_graph) << " " << convex_constr_id
                 << " " << convex_rhs << "\n\n";

    for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
        auto& head =
            table.node(mip_graph[source(*it.first, mip_graph)].node_id);
        auto& n = table.node(mip_graph[target(*it.first, mip_graph)].node_id);
        if (mip_graph[*it.first].high) {
            double cost =
                head.get_job()->weighted_tardiness_start(head.get_weight());
            out_file_mip << head.get_key() << " " << n.get_key() << " " << cost
                         << '\n';
            index_edge[head.get_nb_job()].push_back(mip_graph[*it.first].id);
        } else {
            out_file_mip << head.get_key() << " " << n.get_key() << " " << 0.0
                         << '\n';
        }
    }

    out_file_mip << '\n';

    for (auto i = 0UL; i < convex_constr_id; i++) {
        out_file_mip << index_edge[i].size() << " ";
        for (auto& it : index_edge[i]) {
            out_file_mip << it << " ";
        }
        out_file_mip << '\n';
    }

    out_file_mip << '\n';

    for (auto it = vertices(mip_graph); it.first != it.second; it.first++) {
        auto const& node_id = mip_graph[*it.first].node_id;
        auto&       n = table.node(node_id);
        if (node_id > 1) {
            out_file_mip << n.get_nb_job() << " " << n.get_weight() << '\n';
        }
    }

    out_file_mip << "99 99\n";

    out_file_mip.close();
}

auto PricerSolverBdd::evaluate_nodes(std::span<const double>& pi) -> bool {
    auto& table = *get_decision_diagram().getDiagram();
    compute_labels(pi);
    auto reduced_cost =
        table.node(get_decision_diagram().root()).backward_label[0].get_f();

    auto removed_edges = false;
    auto nb_removed_edges_evaluate = 0;
    auto aux_nb_machines = static_cast<double>(convex_rhs - 1);

    /** check for each node the Lagrangian dual */
    for (auto& it :
         table | ranges::views::take(get_decision_diagram().topLevel() + 1) |
             ranges::views::drop(1) | ranges::views::reverse |
             ranges::views::join) {
        if (((constLB + aux_nb_machines * reduced_cost + evaluate_rc_arc(it) >
              UB - 1.0 + RC_FIXING) ||
             (it.get_job()->weighted_tardiness_start(it.get_weight()) >
              UB - 1.0 + RC_FIXING)) &&
            (it.get_calc(true))) {
            it.update_calc(true);
            removed_edges = true;
            add_nb_removed_edges();
            nb_removed_edges_evaluate++;
        }
        if (((constLB + aux_nb_machines * reduced_cost +
                  evaluate_rc_low_arc(it) >
              UB - 1.0 + RC_FIXING)) &&
            (it.get_calc(false))) {
            it.update_calc(false);
            removed_edges = true;
            add_nb_removed_edges();
            nb_removed_edges_evaluate++;
        }
    }

    if (removed_edges) {
        if (debug_lvl(0)) {
            fmt::print("{0: <{2}}{1}\n",
                       "Number of edges removed by evaluate "
                       "nodes",
                       nb_removed_edges_evaluate, ALIGN);
            fmt::print("{0: <{2}}{1}\n", "Total number of edges removed",
                       get_nb_removed_edges(), ALIGN);
            fmt::print("{0: <{2}}{1}\n", "Number of edges", get_nb_edges(),
                       ALIGN);
        }
        remove_layers();
        remove_edges();
        bottom_up_filtering();
        topdown_filtering();
        cleanup_arcs();
        construct_mipgraph();
    }

    return removed_edges;
}


void PricerSolverBdd::build_mip() {
    try {
        fmt::print("Building Mip model for the extended formulation:\n");
        auto& table = *(decision_diagram.getDiagram());

        /** Constructing variables */
        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            auto& high = mip_graph[*it.first].high;
            auto& x = mip_graph[*it.first].x;
            auto& n =
                table.node(mip_graph[source(*it.first, mip_graph)].node_id);
            if (high) {
                double cost =
                    n.get_job()->weighted_tardiness_start(n.get_weight());
                x = model.addVar(0.0, 1.0, cost, GRB_BINARY);
            } else {
                x = model.addVar(0.0, static_cast<double>(convex_rhs), 0.0,
                                 GRB_INTEGER);
            }
        }

        model.update();

        /** Assignment constraints */
        auto assignment{
            std::vector<GRBLinExpr>(convex_constr_id, GRBLinExpr())};
        auto sense{std::vector<char>(convex_constr_id, GRB_EQUAL)};
        auto rhs{std::vector<double>(convex_constr_id, 1.0)};

        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            auto& high = mip_graph[*it.first].high;
            if (high) {
                auto& x = mip_graph[*it.first].x;
                auto& n =
                    table.node(mip_graph[source(*it.first, mip_graph)].node_id);
                assignment[n.get_nb_job()] += x;
            }
        }

        std::unique_ptr<GRBConstr> assignment_constrs(
            model.addConstrs(assignment.data(), sense.data(), rhs.data(),
                             nullptr, static_cast<int>(convex_constr_id)));
        model.update();

        /** Flow constraints */
        auto num_vertices = boost::num_vertices(mip_graph);
        auto flow_conservation_constr{
            std::vector<GRBLinExpr>(num_vertices, GRBLinExpr())};
        auto sense_flow{std::vector<char>(num_vertices, GRB_EQUAL)};
        auto rhs_flow(std::vector<double>(num_vertices, 0));

        for (auto it = vertices(mip_graph); it.first != it.second; ++it.first) {
            const auto node_id = mip_graph[*it.first].node_id;
            const auto vertex_key = mip_graph[*it.first].index;
            sense_flow[vertex_key] = GRB_EQUAL;
            auto out_edges_it = boost::out_edges(*it.first, mip_graph);

            for (; out_edges_it.first != out_edges_it.second;
                 ++out_edges_it.first) {
                flow_conservation_constr[vertex_key] -=
                    mip_graph[*out_edges_it.first].x;
            }

            auto in_edges_it = boost::in_edges(*it.first, mip_graph);

            for (; in_edges_it.first != in_edges_it.second;
                 ++in_edges_it.first) {
                flow_conservation_constr[vertex_key] +=
                    mip_graph[*in_edges_it.first].x;
            }

            if (node_id == decision_diagram.root()) {
                rhs_flow[vertex_key] = -static_cast<double>(convex_rhs);
            } else if (node_id == 1) {
                rhs_flow[vertex_key] = static_cast<double>(convex_rhs);
            } else {
                rhs_flow[vertex_key] = 0.0;
            }
        }

        std::unique_ptr<GRBConstr> flow_constrs(model.addConstrs(
            flow_conservation_constr.data(), sense_flow.data(), rhs_flow.data(),
            nullptr, static_cast<int>(num_vertices)));
        model.update();
        for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
            auto& high = mip_graph[*it.first].high;
            auto& x = mip_graph[*it.first].x;
            auto& n =
                table.node(mip_graph[source(*it.first, mip_graph)].node_id);
            x.set(GRB_DoubleAttr_Start, n.get_best_sol_x(high));
            // if (high) {
            //     x.set(GRB_DoubleAttr_Start, n.best_sol_x[1]);

            // } else {
            //     x.set(GRB_DoubleAttr_Start, n.best_sol_x[0]);
            // }
        }
        model.optimize();
    } catch (GRBException& e) {
        fmt::print("Error code = {}\n", e.getErrorCode());
        fmt::print("{}", e.getMessage());
    } catch (...) {
        fmt::print("Exception during optimization\n");
    }
}

void PricerSolverBdd::cleanup_arcs() {
    auto& table = *(decision_diagram.getDiagram());

    table.node(0).reset_backward_distance(std::numeric_limits<int>::min());
    table.node(1).reset_backward_distance(0);
    auto removed_edges = false;
    auto nb_edges_removed_tmp = 0;

    BackwardDistance backward_distance_eval;

    decision_diagram.compute_labels_backward(backward_distance_eval);

    /** remove the unnecessary nodes of the bdd */
    for (auto& iter : table |
                          ranges::views::take(decision_diagram.topLevel() + 1) |
                          ranges::views::drop(1) | ranges::views::join) {
        // if (iter.get_weight() + iter.backward_distance[0] < H_min &&
        //     iter.branch[0] != 0) {
        //     iter.calc[0] = false;
        //     removed_edges = true;
        //     nb_edges_removed_tmp++;
        //     nb_removed_edges++;
        // }

        if (iter.get_weight() + iter.get_backward_distance()[1] < H_min) {
            iter.update_calc(true);
            removed_edges = true;
            nb_edges_removed_tmp++;
            nb_removed_edges++;
        }
    }

    if (removed_edges) {
        if (debug_lvl(0)) {
            fmt::print("{0: <{2}}{1}\n", "Number of edges removed by clean up",
                       nb_edges_removed_tmp, ALIGN);
            fmt::print("{0: <{2}}{1}\n", "Total number of edges removed",
                       get_nb_removed_edges(), ALIGN);
        }
        remove_layers();
        remove_edges();
    }
}

void PricerSolverBdd::topdown_filtering() {
    auto  removed_edges = false;
    auto  nb_edges_removed_tmp = 0;
    auto& table = *(decision_diagram.getDiagram());
    auto& root = table.node(decision_diagram.root());
    root.init_node(0, true);
    for (auto& it : table |
                        ranges::views::take(decision_diagram.root().row() + 1) |
                        ranges::views::reverse | ranges::views::join) {
        it.reset_visited();
        it.reset_all(convex_constr_id);
    }

    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        auto& n0 = table.node(it[0]);

        if (n0.get_visited()) {
            n0.intersect_all(it.get_all());
        } else {
            n0.update_visited();
            n0.reset_all(convex_constr_id);
            n0.union_all(it.get_all());
        }

        auto& n1 = table.node(it[1]);

        if (n1.get_visited()) {
            if (n1.is_element(it.get_nb_job())) {
                n1.intersect_all(it.get_all());
                n1.add_element(it.get_nb_job());
            } else {
                n1.intersect_all(it.get_all());
            }
        } else {
            n1.reset_all(convex_constr_id);
            n1.union_all(it.get_all());
            n1.add_element(it.get_nb_job());
            n1.update_visited();
        }
    }

    for (auto i :
         ranges::views::ints(size_t{1}, decision_diagram.root().row() + 1) |
             ranges::views::reverse) {
        for (auto& it : table[i]) {
            if (it.is_element(it.get_nb_job())) {
                removed_edges = true;
                it.update_calc(true);
                nb_removed_edges++;
                nb_edges_removed_tmp++;
            }
        }
    }

    if (removed_edges) {
        remove_edges();
        cleanup_arcs();
    }
}

void PricerSolverBdd::bottom_up_filtering() {
    auto  removed_edges = false;
    auto  nb_edges_removed_tmp = 0;
    auto& table = *(decision_diagram.getDiagram());
    for (auto i :
         ranges::views::ints(size_t{1}, decision_diagram.root().row() + 1) |
             ranges::views::reverse) {
        for (auto& it : table[i]) {
            it.reset_visited();
            it.reset_all(convex_constr_id);
        }
    }

    table.node(0).reset_all(convex_constr_id);
    table.node(1).reset_all(convex_constr_id);
    table.node(0).get_all().flip();

    for (auto i :
         ranges::views::ints(size_t{1}, decision_diagram.root().row() + 1)) {
        for (auto& it : table[i]) {
            it.add_element(it.get_nb_job());
            it.union_all(table.node(it[1]).get_all());
            it.intersect_all(table.node(it[0]).get_all());
        }
    }

    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        if (table.node(it[1]).is_element(it.get_nb_job())) {
            removed_edges = true;
            it.update_calc(true);
            nb_removed_edges++;
            nb_edges_removed_tmp++;
        }
    }

    if (removed_edges) {
        remove_edges();
        cleanup_arcs();
    }
}

[[maybe_unused]] void PricerSolverBdd::check_infeasible_arcs() {
    /** init table */
    auto  removed_edges = false;
    auto  nb_edges_removed_tmp = 0;
    auto& table = *(decision_diagram.getDiagram());

    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::reverse | ranges::views::join) {
        it.reset_visited();
        it.reset_all(convex_constr_id);
        it.reset_calc();
    }

    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        auto& n0 = table.node(it[0]);
        n0.union_all(it.get_all());
        auto& n1 = table.node(it[1]);
        n1.add_element(it.get_nb_job());
    }

    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        if (!it.all_is_empty() &&
            it.get_first_all() != boost::dynamic_bitset<>::npos) {
            auto index = it.get_first_all();

            auto max = value_diff_Fij(it.get_weight(), it.get_job(),
                                      jobs[index].get());
            while (index != boost::dynamic_bitset<>::npos && max < 0) {
                index = it.get_next_all(index);
                if (index != boost::dynamic_bitset<>::npos) {
                    int a = value_diff_Fij(it.get_weight(), it.get_job(),
                                           jobs[index].get());
                    if (a > max) {
                        max = a;
                    }
                }
            }

            if (max < 0) {
                removed_edges = true;
                it.update_calc(true);
                nb_removed_edges++;
                nb_edges_removed_tmp++;
            }
        }
    }

    if (removed_edges) {
        fmt::print("removing edges based on order\n");
        fmt::print("Number edges removed order = {}\n", nb_edges_removed_tmp);
        fmt::print("Number edges removed total = {}\n", nb_removed_edges);
        remove_layers();
        remove_edges();
        cleanup_arcs();
    }
}

void PricerSolverBdd::equivalent_paths_filtering() {
    /** init table */
    auto  removed_edges = false;
    auto  nb_edges_removed_tmp = 0;
    auto& table = *(decision_diagram.getDiagram());
    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        it.reset_visited();
        it.reset_all(convex_constr_id);
        auto& n0 = table.node(it[0]);
        n0.update_in_degree(false);
        auto& n1 = table.node(it[1]);
        n1.update_in_degree(true);
    }

    std::vector<size_t> vertices;

    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        if (it.alternative_paths()) {
            vertices.push_back(it.get_key());
        }
    }

    for (auto& it : vertices) {
        auto              start_v = mip_graph[it].node_id;
        auto              num_vertices = boost::num_vertices(mip_graph);
        std::list<NodeId> queue;

        queue.push_back(start_v);
        std::vector<bool>                    visited(num_vertices, false);
        std::vector<bool>                    edge_visited(num_vertices, false);
        std::vector<boost::dynamic_bitset<>> all(
            num_vertices, boost::dynamic_bitset(convex_constr_id, 0));
        std::vector<int> C(num_vertices, 0);

        // auto& tmp_n = table.node(start_v);
        // visited[tmp_n.get_key()];
        auto stop = false;

        while (!queue.empty()) {
            auto currVertex = queue.front();
            queue.pop_front();
            auto iter =
                boost::in_edges(table.node(currVertex).get_key(), mip_graph);

            for (; iter.first != iter.second; iter.first++) {
                auto adjVertex =
                    mip_graph[source(*iter.first, mip_graph)].node_id;
                auto& n = table.node(adjVertex);
                auto* data = static_cast<EdgeData*>(iter.first->get_property());
                auto  high = data->high;

                if (!visited[n.get_key()]) {
                    visited[n.get_key()] = true;
                    queue.push_back(adjVertex);
                    if (high) {
                        auto& tmp_node = table.node(n[1]);
                        all[n.get_key()] |= all[tmp_node.get_key()];
                        all[n.get_key()][n.get_nb_job()] = true;
                        C[n.get_key()] = C[tmp_node.get_key()] +
                                         n.get_job()->weighted_tardiness(
                                             tmp_node.get_weight());
                        edge_visited[n.get_key()] = true;
                    } else {
                        auto& tmp_node = table.node(n[0]);
                        all[n.get_key()] |= all[tmp_node.get_key()];
                        C[n.get_key()] = C[tmp_node.get_key()];
                    }
                } else {
                    boost::dynamic_bitset<> tmp;
                    int                     tmp_C{};
                    if (high) {
                        auto& tmp_node = table.node(n[1]);
                        tmp = all[tmp_node.get_key()];
                        tmp[n.get_nb_job()] = true;
                        tmp_C = C[tmp_node.get_key()] +
                                n.get_job()->weighted_tardiness(
                                    tmp_node.get_weight());
                    } else {
                        auto& tmp_node = table.node(n[0]);
                        tmp = all[tmp_node.get_key()];
                        tmp_C = C[tmp_node.get_key()];
                    }

                    if (all[n.get_key()] == tmp) {
                        NodeId cur{};
                        NodeId prev = adjVertex;
                        if (high) {
                            if (tmp_C > C[n.get_key()]) {
                                cur = n[1];
                            } else {
                                cur = n[0];
                            }
                        } else {
                            if (tmp_C > C[n.get_key()]) {
                                cur = n[0];
                            } else {
                                cur = n[1];
                            }
                        }

                        while (cur != start_v) {
                            auto& node = table.node(cur);
                            if (node.alternative_paths()) {
                                break;
                            }
                            prev = cur;
                            assert(cur != NodeId(0, 0));
                            assert(cur != NodeId(0, 1));
                            if (edge_visited[node.get_key()]) {
                                cur = node[1];
                            } else {
                                cur = node[0];
                            }
                        }

                        auto& node_delete = table.node(prev);
                        if (edge_visited[node_delete.get_key()] &&
                            cur == start_v) {
                            node_delete.update_calc(true);
                            removed_edges = true;
                            nb_edges_removed_tmp++;
                        } else if (cur == start_v) {
                            node_delete.update_calc(false);
                            removed_edges = true;
                            nb_edges_removed_tmp++;
                        }
                    }
                    stop = true;
                    break;
                }
            }

            if (stop) {
                break;
            }
        }
    }

    if (removed_edges) {
        fmt::print(
            "Number of edges removed by equivalent_path_filtering = "
            "{}\nNumber "
            "of edges removed in total = {}\n",
            nb_edges_removed_tmp, nb_removed_edges);

        remove_layers();
        remove_edges();
        cleanup_arcs();
        construct_mipgraph();
    }
}

void PricerSolverBdd::construct_lp_sol_from_rmp(
    const std::span<const double>&              lambda,
    const std::vector<std::shared_ptr<Column>>& columns) {
    auto& table = *(decision_diagram.getDiagram());
    for (auto& it : table |
                        ranges::views::take(decision_diagram.root().row() + 1) |
                        ranges::views::join) {
        it.reset_lp_x();
    }

    set_is_integer_solution(true);
    for (auto&& [x, tmp] : ranges::views::zip(lambda, columns)) {
        if (x > EPS_SOLVER) {
            // auto* tmp = columns[i].get();
            auto it = tmp->job_list.begin();
            auto tmp_nodeid(decision_diagram.root());

            while (tmp_nodeid > 1) {
                Job*  tmp_j = (it != tmp->job_list.end()) ? *it : nullptr;
                auto& tmp_node = table.node(tmp_nodeid);

                auto high = (tmp_j == tmp_node.get_job());
                tmp_node.update_lp_x(x, high);
                tmp_nodeid = high ? tmp_node[1] : tmp_node[0];
                if (high) {
                    ++it;
                }
            }

            assert(tmp_nodeid == 1);
        }
    }

    ranges::for_each(x_bar, [](auto& tmp) { tmp.clear(); });
    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        it.update_lp_visited(false);
        auto value = it.get_lp_x()[1];
        if (value > EPS_SOLVER) {
            x_bar[it.get_nb_job()].emplace_back(it.get_nb_job(),
                                                it.get_weight(), 0.0, value);
            if (safe_frac(value) > EPS_SOLVER) {
                set_is_integer_solution(false);
            }
        }
    }

    if (is_integer_solution && debug_lvl(1)) {
        fmt::print("FOUND INTEGER SOLUTION\n\n");
    }

    // ColorWriterEdgeX  edge_writer(mip_graph, &table);
    // ColorWriterVertex vertex_writer(mip_graph, table);
    // auto              file_name = "lp_solution_" + problem_name + "_" +
    //                  std::to_string(convex_rhs) + ".gv";
    // std::ofstream outf(file_name);
    // boost::write_graphviz(outf, mip_graph, vertex_writer, edge_writer);
    // outf.close();
}


void PricerSolverBdd::project_sol_on_original_variables(const Sol& _sol) {
    auto& table = *(decision_diagram.getDiagram());
    for (const auto& m : _sol.machines) {
        auto tmp_nodeid(decision_diagram.root());
        auto it = m.job_list.begin();

        while (tmp_nodeid > 1) {
            auto* tmp_j = (it != m.job_list.end()) ? *it : nullptr;
            auto& tmp_node = table.node(tmp_nodeid);

            auto high = (tmp_j == tmp_node.get_job());
            tmp_node.update_best_sol_x(high);
            tmp_nodeid = high ? tmp_node[1] : tmp_node[0];
            if (high) {
                ++it;
            }
        }
        assert(tmp_nodeid == 1);
    }
}

auto PricerSolverBdd::calculate_job_time()
    -> std::vector<std::vector<BddCoeff>>& {
    for (auto&& it : x_bar) {
        std::ranges::sort(it, std::less<>{},
                          [](const auto& aux) { return aux.get_t(); });
        auto value = 0.0;
        for (auto& x : it) {
            value += x.get_value();
            x.update_cum_value(value);
        }
    }

    return x_bar;
}

void PricerSolverBdd::split_job_time(size_t _job, int _time, bool _left) {
    auto& table = *(decision_diagram.getDiagram());
    auto  removed_edges = false;

    for (auto& it : table |
                        ranges::views::take(decision_diagram.topLevel() + 1) |
                        ranges::views::drop(1) | ranges::views::reverse |
                        ranges::views::join) {
        if (_left) {
            if (it.get_weight() <= static_cast<size_t>(_time) &&
                it.get_nb_job() == _job) {
                it.update_calc(true);
                removed_edges = true;
            }
        } else {
            if (it.get_weight() > static_cast<size_t>(_time) &&
                it.get_nb_job() == _job) {
                it.update_calc(true);
                removed_edges = true;
            }
        }
    }

    if (removed_edges) {
        decision_diagram.compressBdd();
        remove_layers();
        remove_edges();
        bottom_up_filtering();
        topdown_filtering();
        cleanup_arcs();
        construct_mipgraph();

        if (debug_lvl(0)) {
            auto&             table_bis = *(decision_diagram.getDiagram());
            ColorWriterVertex vertex_writer(mip_graph, table_bis);
            auto file_name = fmt::format("split_solution_{}_{}_{}_{}.gv",
                                         problem_name, _job, _time, _left);
            std::ofstream output_stream(file_name);
            boost::write_graphviz(output_stream, mip_graph, vertex_writer);
            output_stream.close();
        }
    }
}

auto PricerSolverBdd::build_model() -> std::unique_ptr<OsiSolverInterface> {
    auto& table = *(decision_diagram.getDiagram());
    auto  aux_model = std::make_unique<OsiGrbSolverInterface>();
    nb_vertices = 0UL;
    nb_edges = 0UL;
    auto nz = 0UL;
    auto view_table = ranges::views::take(decision_diagram.topLevel() + 1) |
                      ranges::views::drop(1) | ranges::views::reverse |
                      ranges::views::join;

    for (auto& it : table | ranges::views::join) {
        it.reset_in_edges();
    }

    for (auto& it : table | view_table) {
        it.set_key_model(nb_vertices);
        ++nb_vertices;
        if (it[0] != 0 && it.get_calc(false)) {
            auto& n0 = table.node(it[0]);
            n0.add_in_edge(false, nb_edges);
            it.set_key_edge(false, nb_edges);
            ++nb_edges;
        }

        if (it[1] != 0 && it.get_calc(true)) {
            auto& n1 = table.node(it[1]);
            n1.add_in_edge(true, nb_edges);
            it.set_key_edge(true, nb_edges);
            ++nb_edges;
        }
    }

    auto terminal_node = table.node(1);
    terminal_node.set_key_model(nb_vertices);
    ++nb_vertices;

    auto cost = std::vector<double>(nb_edges);
    auto high_edge = std::vector<bool>(nb_edges);
    auto rhs_lb = std::vector<double>(nb_edges);
    auto rhs_ub = std::vector<double>(nb_edges);
    auto nodes = std::vector<std::pair<size_t, size_t>>(nb_edges);
    auto edge_job = std::vector<int>(nb_edges);

    nb_edges = 0;
    for (auto& it : table | view_table) {
        for (auto& tmp : it.get_in_edges(false)) {
            nodes[tmp].second = it.get_key_model();
        }

        for (auto& tmp : it.get_in_edges(true)) {
            nodes[tmp].second = it.get_key_model();
        }

        if (it[0] != 0 && it.get_calc(false)) {
            cost[nb_edges] = 0;
            rhs_lb[nb_edges] = 0;
            rhs_ub[nb_edges] =
                aux_model->getInfinity();  // static_cast<double>(convex_rhs);
            high_edge[nb_edges] = false;
            edge_job[nb_edges] = static_cast<int>(it.get_nb_job());
            nodes[nb_edges].first = it.get_key_model();
            ++nb_edges;
            nz += 2;
        }

        if (it[1] != 0 && it.get_calc(true)) {
            cost[nb_edges] =
                it.get_job()->weighted_tardiness_start(it.get_weight());
            rhs_lb[nb_edges] = 0;
            rhs_ub[nb_edges] = aux_model->getInfinity();
            high_edge[nb_edges] = true;
            edge_job[nb_edges] = static_cast<int>(it.get_nb_job());
            nodes[nb_edges].first = it.get_key_model();
            ++nb_edges;
            nz += 3;
        }
    }

    for (auto& tmp : terminal_node.get_in_edges(false)) {
        nodes[tmp].second = terminal_node.get_key_model();
    }

    for (auto& tmp : terminal_node.get_in_edges(true)) {
        nodes[tmp].second = terminal_node.get_key_model();
    }

    std::vector<int>    start_vars(nb_edges + 1);
    std::vector<int>    rows(nz);
    std::vector<double> coeff(nz);

    auto i = 0UL;
    nz = 0;
    for (auto&& [j, h, n] : ranges::views::zip(edge_job, high_edge, nodes)) {
        start_vars[i] = nz;
        if (h) {
            rows[nz] = j;
            coeff[nz] = 1.0;
            nz += 1;
        }
        rows[nz] = static_cast<int>(convex_constr_id + n.first);
        coeff[nz] = 1.0;
        nz += 1;
        rows[nz] = static_cast<int>(convex_constr_id + n.second);
        coeff[nz] = -1.0;
        nz += 1;
        ++i;
    }
    start_vars[i] = nz;

    std::vector<int>    start(nb_vertices + 1 + convex_constr_id, 0);
    std::vector<double> rhs(nb_vertices + convex_constr_id, 0.0);
    for (auto& it : rhs | ranges::views::take(convex_constr_id)) {
        it = 1.0;
    }
    rhs[convex_constr_id] = static_cast<double>(convex_rhs);
    rhs[convex_constr_id + terminal_node.get_key_model()] =
        -static_cast<double>(convex_rhs);

    aux_model->addRows(static_cast<int>(nb_vertices + convex_constr_id),
                       start.data(), nullptr, nullptr, rhs.data(), rhs.data());
    aux_model->addCols(static_cast<int>(nb_edges), start_vars.data(),
                       rows.data(), coeff.data(), rhs_lb.data(), rhs_ub.data(),
                       cost.data());

    return aux_model;
}

auto PricerSolverBdd::add_constraints() -> int {
    // auto& table = *(decision_diagram.getDiagram());
    // auto  aux_model = build_model();
    // fmt::print("number of nodes = {}\n", table.size());
    // auto sol_x = std::vector<double>(aux_model->getNumCols(), 0.0);

    // for (auto& it : table |
    //                     ranges::views::take(decision_diagram.topLevel() + 1)
    //                     | ranges::views::drop(1) | ranges::views::reverse |
    //                     ranges::views::join) {
    //     sol_x[it.get_key_edge(false)] = it.get_lp_x(false);
    //     sol_x[it.get_key_edge(true)] = it.get_lp_x(true);
    // }

    // // aux_model->setColSolution(sol_x.data());
    // aux_model->messageHandler()->setLogLevel(0);
    // for (int i = 0; i < aux_model->getNumCols(); ++i) {
    //     aux_model->setInteger(i);
    // }
    // aux_model->initialSolve();
    // OsiCuts     cuts;
    // CglZeroHalf cg;
    // cg.refreshSolver(aux_model.get());
    // cg.generateCuts(*aux_model, cuts);
    //     fmt::print("new LB {}\n", aux_model->getObjValue());
    // if (cuts.sizeRowCuts() > 0) {
    //     aux_model->applyCuts(cuts);
    //     aux_model->resolve();
    //     fmt::print("new LB {}\n", aux_model->getObjValue());
    // }
    // fmt::print("Number of Cuts {} {}\n", cuts.sizeRowCuts(),
    //            aux_model->getObjValue());

    // fmt::print("------------------------------\n");
    // aux_model->writeLp("test");
    // construct_mipgraph();
    // ColorWriterEdgeX  edge_writer(mip_graph, &table);
    // ColorWriterVertex vertex_writer(mip_graph, table);
    // auto              file_name = "lp_solution_" + problem_name + "_" +
    //                  std::to_string(convex_rhs) + ".gv";
    // std::ofstream outf(file_name);
    // boost::write_graphviz(outf, mip_graph, vertex_writer, edge_writer);
    // outf.close();
    // getchar();

    return 0;
}

void PricerSolverBdd::remove_constraints(int first, int nb_del) {
    original_model.delete_constraints(first, nb_del);
    reformulation_model.delete_constraints(first, nb_del);
}

void PricerSolverBdd::update_rows_coeff(size_t first) {
    for (auto k :
         ranges::views::ints(first, original_model.get_nb_constraints())) {
        auto& aux = original_model.get_coeff_list(k);
        for (auto& it : aux) {
            it->set_row(k);
        }
    }
}

auto PricerSolverBdd::check_column(Column const* col) -> bool {
    const auto& set = col->job_list;
    auto&       table = *(decision_diagram.getDiagram());
    auto        tmp_nodeid(decision_diagram.root());
    auto        it = set.begin();
    auto        counter = col->job_list.size();

    while (tmp_nodeid > 1) {
        auto& tmp_node = table.node(tmp_nodeid);
        auto* tmp_j = (it != set.end()) ? *it : nullptr;

        if (tmp_j == tmp_node.get_job()) {
            tmp_nodeid = tmp_node[1];
            --counter;
            ++it;
        } else {
            tmp_nodeid = tmp_node[0];
        }
    }

    return (tmp_nodeid == 1 && it == set.end() && counter == 0);
}

auto PricerSolverBdd::evaluate_rc_low_arc(NodeBdd& n) -> double {
    auto& table = *(get_decision_diagram().getDiagram());
    auto& child = table.node(n[0]);
    return n.forward_label[0].get_f() + child.backward_label[0].get_f() +
           n.get_reduced_cost()[0];
}

auto PricerSolverBdd::get_nb_edges() -> size_t {
    return num_edges(mip_graph);
}

auto PricerSolverBdd::get_nb_vertices() -> size_t {
    return num_vertices(mip_graph);
}

auto PricerSolverBdd::structure_feasible() -> bool {
    return ((decision_diagram.size() != 0) &&
            decision_diagram.root().row() >= convex_constr_id);
}

auto PricerSolverBdd::print_num_paths() -> boost::multiprecision::cpp_int {
    auto evaluator = CardinalityPaths();
    return decision_diagram.evaluate_backward(evaluator);
}

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

#include "PricerSolverZdd.hpp"
#include <fmt/core.h>                              // for print
#include <gurobi_c++.h>                            // for GRBException
#include <algorithm>                               // for remove_if, fill
#include <array>                                   // for array
#include <boost/graph/detail/adjacency_list.hpp>   // for num_edges
#include <boost/multiprecision/cpp_int.hpp>        // for cpp_int
#include <range/v3/iterator/reverse_iterator.hpp>  // for reverse_cursor
#include <range/v3/view/iota.hpp>                  // for iota_view, ints
#include <range/v3/view/reverse.hpp>               // for reverse_fn, revers...
#include <range/v3/view/zip.hpp>                   // for zip_view
#include <set>                                     // for operator==, _Rb_tr...
#include <span>                                    // for span
#include <string>                                  // for char_traits, opera...
#include <vector>                                  // for vector<>::iterator
#include "Column.h"                                // for ScheduleSet
#include "DebugLvl.hpp"                            // for debug_lvl
#include "Instance.h"                              // for Instance
#include "Job.h"                                   // for Job
#include "ModernDD/NodeBddStructure.hpp"           // for DdStructure, DdStr...
#include "ModernDD/NodeBddTable.hpp"               // for NodeTableEntity
#include "ModernDD/NodeId.hpp"                     // for NodeId
#include "PricerConstruct.hpp"                     // for PricerConstruct
#include "PricingSolution.hpp"                     // for PricingSolution
#include "ZddNode.hpp"                             // for NodeZdd, SubNodeZdd

PricerSolverZdd::PricerSolverZdd(const Instance& instance)
    : PricerSolverBase(instance),
      decision_diagram(
          std::make_unique<DdStructure<NodeZdd<>>>(PricerConstruct(instance))),
      size_graph{decision_diagram->size()} {
    remove_layers_init();
    // decision_diagram->compressBdd();
    // decision_diagram->reduceZdd();
    init_table();
    construct_mipgraph();
    lp_x = std::vector<double>(num_edges(mip_graph), 0.0);
    solution_x = std::vector<double>(num_edges(mip_graph), 0.0);
}

void PricerSolverZdd::construct_mipgraph() {
    //     mip_graph.clear();
    //     auto&             table = *(decision_diagram->getDiagram());
    //     NodeZddIdAccessor vertex_node_zdd_id_list(
    //         get(boost::vertex_color_t(), mip_graph));
    //     NodeIdAccessor vertex_nodeid_list(get(boost::vertex_name_t(),
    //     mip_graph)); NodeMipIdAccessor vertex_mip_id_list(
    //         get(boost::vertex_degree_t(), mip_graph));
    //     EdgeTypeAccessor edge_type_list(get(boost::edge_weight_t(),
    //     mip_graph));

    //     for (int i = decision_diagram->topLevel(); i >= 0; i--) {
    //         for (size_t j = 0; j < table[i].size(); j++) {
    //             auto n{NodeId(i, j)};
    //             if (n.row() != 0) {
    //                 for (auto& it : table[i][j].list) {
    //                     auto key = add_vertex(mip_graph);
    //                     it->key = key;
    //                     vertex_mip_id_list[it->key] = key;
    //                     vertex_nodeid_list[it->key] = it->node_id;
    //                     vertex_node_zdd_id_list[it->key] = it;
    //                 }
    //             } else {
    //                 if (n != 0) {
    //                     auto terminal_node = add_vertex(mip_graph);
    //                     for (auto& it : table[i][j].list) {
    //                         it->key = terminal_node;
    //                         vertex_mip_id_list[terminal_node] =
    //                         terminal_node; vertex_nodeid_list[terminal_node]
    //                         = it->node_id;
    //                         vertex_node_zdd_id_list[terminal_node] = it;
    //                     }
    //                 }
    //             }
    //         }
    //     }

    //     int count = 0;

    //     for (int i = decision_diagram->topLevel(); i > 0; i--) {
    // #ifndef NDEBUG
    //         auto job = (static_cast<job_interval_pair*>(
    //                         g_ptr_array_index(ordered_jobs, ordered_jobs->len
    //                         - i)))
    //                        ->j;
    // #endif  // NDEBUG

    //         for (auto& it : table[i]) {
    //             if (it.branch[0] != 0) {
    //                 for (auto& iter : it.list) {
    //                     auto n = iter->n;
    //                     assert(iter->weight == n->weight);
    //                     auto a = add_edge(iter->key, n->key, mip_graph);
    //                     put(edge_type_list, a.first, false);
    //                     iter->low_edge_key = count;
    //                     put(boost::edge_index_t(), mip_graph, a.first,
    //                     count++);
    //                 }
    //             }

    //             if (it.branch[1] != 0) {
    //                 for (auto& iter : it.list) {
    //                     auto y = iter->y;
    //                     assert(iter->weight + job->processing_time ==
    //                     y->weight); auto a = add_edge(iter->key, y->key,
    //                     mip_graph); put(edge_type_list, a.first, true);
    //                     iter->high_edge_key = count;
    //                     put(boost::edge_index_t(), mip_graph, a.first,
    //                     count++);
    //                 }
    //             }
    //         }
    //     }

    //     std::cout << "Number of vertices = " << num_vertices(mip_graph) <<
    //     '\n'; std::cout << "Number of edges = " << num_edges(mip_graph) <<
    //     '\n';
}

void PricerSolverZdd::init_table() {
    auto& table = *(decision_diagram->getDiagram());
    /** init table */
    auto& n = table.node(decision_diagram->root());
    // std::span span_ordered_job{ordered_jobs->pdata, ordered_jobs->len};
    n.add_sub_node(0, decision_diagram->root(), true, false);
    n.set_node_id(decision_diagram->root());
    if (table.node(decision_diagram->root()).list.size() > 1) {
        table.node(decision_diagram->root()).list.pop_back();
    }

    for (auto i :
         ranges::views::ints(size_t{}, decision_diagram->root().row() + 1) |
             ranges::views::reverse) {
        size_t layer = ordered_jobs_new.size() - i;
        auto&  tmp_pair = ordered_jobs_new[layer];
        // static_cast<job_interval_pair*>(span_ordered_job[layer]);

        for (auto& it : table[i]) {
            if (i != 0UL) {
                it.set_job(tmp_pair.first);
                auto& n0 = table.node(it[0]);
                auto& n1 = table.node(it[1]);
                int   p = it.get_job()->processing_time;
                it.child[0] = table.node_ptr(it[0]);
                it.child[1] = table.node_ptr(it[1]);
                for (auto& iter : it.list) {
                    int w = iter->weight;
                    iter->n = n0.add_weight(w, it[0]);
                    iter->y = n1.add_weight(w + p, it[1]);
                }
            } else {
                it.set_job(nullptr);
            }
        }
    }
}

auto PricerSolverZdd::farkas_pricing(
    [[maybe_unused]] std::span<const double>& pi) -> PricingSolution {
    PricingSolution sol;

    return sol;
}

void PricerSolverZdd::remove_layers_init() {
    auto& table = *(decision_diagram->getDiagram());

    auto i = decision_diagram->root().row();
    ordered_jobs_new.erase(
        std::remove_if(ordered_jobs_new.begin(), ordered_jobs_new.end(),
                       [&]([[maybe_unused]] const auto& tmp) {
                           bool remove = std::ranges::all_of(
                               table[i],
                               [&](const auto& n) { return n[1] == 0; });
                           --i;
                           return remove;
                       }),
        ordered_jobs_new.end());

    if (debug_lvl(0)) {
        fmt::print("{0: <{2}}{1}\n", "The new number of layers",
                   ordered_jobs_new.size(), ALIGN);
    }

    fmt::print("The new number of layers = {}\n", ordered_jobs_new.size());
}

void PricerSolverZdd::remove_layers() {
    // int   first_del = -1;
    // int   last_del = -1;
    // int   it = 0;
    auto& table = *(decision_diagram->getDiagram());

    // /** remove the unnecessary layers of the bdd */
    // for (int i = decision_diagram->topLevel(); i > 0; i--) {
    //     bool remove_layer = true;

    //     for (auto& iter : table[i]) {
    //         auto end =
    //             std::remove_if(iter.list.begin(), iter.list.end(),
    //                            [](std::shared_ptr<SubNodeZdd<>> const& n) {
    //                                return !(n->calc_yes);
    //                            });
    //         iter.list.erase(end, iter.list.end());

    //         if (iter.list.empty()) {
    //             NodeId& cur_node_1 = iter.branch[1];
    //             cur_node_1 = 0;
    //         } else {
    //             remove_layer = false;
    //         }
    //     }

    //     if (!remove_layer) {
    //         if (first_del != -1) {
    //             g_ptr_array_remove_range(ordered_jobs, first_del,
    //                                      last_del - first_del + 1);
    //             it = it - (last_del - first_del);
    //             first_del = last_del = -1;
    //         } else {
    //             it++;
    //         }
    //     } else {
    //         if (first_del == -1) {
    //             first_del = it;
    //             last_del = first_del;
    //         } else {
    //             last_del++;
    //         }

    //         it++;
    //     }
    // }

    // if (first_del != -1) {
    //     g_ptr_array_remove_range(ordered_jobs, first_del,
    //                              last_del - first_del + 1);
    // }
    auto i = decision_diagram->root().row();
    ordered_jobs_new.erase(
        std::remove_if(ordered_jobs_new.begin(), ordered_jobs_new.end(),
                       [&]([[maybe_unused]] const auto& tmp) {
                           bool remove_layer = true;

                           for (auto& iter : table[i]) {
                               auto end = std::remove_if(
                                   iter.list.begin(), iter.list.end(),
                                   [](std::shared_ptr<SubNodeZdd<>> const& n) {
                                       return !(n->calc_yes);
                                   });
                               iter.list.erase(end, iter.list.end());

                               if (iter.list.empty()) {
                                   NodeId& cur_node_1 = iter[1];
                                   cur_node_1 = 0;
                               } else {
                                   remove_layer = false;
                               }
                           }

                           --i;
                           return remove_layer;
                       }),
        ordered_jobs_new.end());

    fmt::print("The new number of layers = {}\n", ordered_jobs_new.size());
}

void PricerSolverZdd::remove_edges() {
    // decision_diagram->compressBdd();
    // decision_diagram->reduceZdd();
    nb_removed_nodes -= size_graph;
    size_graph = decision_diagram->size();
    fmt::print("The new size of ZDD = {}\n", size_graph);
}

void PricerSolverZdd::build_mip() {
    try {
        // fmt::print("Building Mip model for the extended formulation:\n");
        // // NodeIdAccessor vertex_nodeid_list(
        // //     get(boost::vertex_name_t(), mip_graph));
        // // NodeMipIdAccessor vertex_mip_id_list(
        // //     get(boost::vertex_degree_t(), mip_graph));
        // // EdgeTypeAccessor edge_type_list(get(boost::edge_weight_t(),
        // // mip_graph)); EdgeVarAccessor
        // // edge_var_list(get(boost::edge_weight2_t(), mip_graph));

        // /** Constructing variables */
        // for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
        //     if (mip_graph[*it.first].high) {
        //         auto& n = decision_diagram->getDiagram()->node(
        //             mip_graph[source(*it.first, mip_graph)].node_id);
        //         Job* job = n.get_job();

        //         double cost = value_Fj(n.+ job->processing_time, job);
        //         edge_var_list[*it.first].x =
        //             model.addVar(0.0, 1.0, cost, GRB_CONTINUOUS);
        //     } else {
        //         edge_var_list[*it.first].x = model.addVar(
        //             0.0, static_cast<double>(convex_rhs), 0.0,
        //             GRB_CONTINUOUS);
        //     }
        // }

        // model.update();
        // /** Assignment constraints */
        // std::vector<GRBLinExpr> assignment(convex_constr_id, GRBLinExpr());
        // std::vector<char>       sense(convex_constr_id, GRB_GREATER_EQUAL);
        // std::vector<double>     rhs(convex_constr_id, 1.0);

        // for (auto it = edges(mip_graph); it.first != it.second; it.first++) {
        //     auto high = edge_type_list[*it.first];

        //     if (high) {
        //         auto& n = get(boost::vertex_color_t(), mip_graph,
        //                       source(*it.first, mip_graph));
        //         assignment[n->get_job()->job] += edge_var_list[*it.first].x;
        //     }
        // }

        // std::unique_ptr<GRBConstr> assignment_constrs(
        //     model.addConstrs(assignment.data(), sense.data(), rhs.data(),
        //                      nullptr, convex_constr_id));
        // model.update();
        // /** Flow constraints */
        // size_t num_vertices = boost::num_vertices(mip_graph);
        // // boost::num_vertices(mip_graph) - table[0][1].list.size() + 1;
        // std::vector<GRBLinExpr> flow_conservation_constr(num_vertices,
        //                                                  GRBLinExpr());
        // std::vector<char>       sense_flow(num_vertices, GRB_EQUAL);
        // std::vector<double>     rhs_flow(num_vertices, 0.0);

        // for (auto it = vertices(mip_graph); it.first != it.second;
        // ++it.first) {
        //     const auto node_id = vertex_nodeid_list[*it.first];
        //     const auto vertex_key = vertex_mip_id_list[*it.first];
        //     auto       out_edges_it = boost::out_edges(*it.first, mip_graph);

        //     for (; out_edges_it.first != out_edges_it.second;
        //          ++out_edges_it.first) {
        //         flow_conservation_constr[vertex_key] -=
        //             edge_var_list[*out_edges_it.first].x;
        //     }

        //     auto in_edges_it = boost::in_edges(*it.first, mip_graph);

        //     for (; in_edges_it.first != in_edges_it.second;
        //          ++in_edges_it.first) {
        //         flow_conservation_constr[vertex_key] +=
        //             edge_var_list[*in_edges_it.first].x;
        //     }

        //     if (node_id == decision_diagram->root()) {
        //         rhs_flow[vertex_key] = static_cast<double>(-convex_rhs);
        //     } else if (node_id.row() == 0) {
        //         rhs_flow[vertex_key] = static_cast<double>(convex_rhs);
        //     }
        // }

        // std::unique_ptr<GRBConstr> flow_constrs(
        //     model.addConstrs(flow_conservation_constr.data(),
        //     sense_flow.data(),
        //                      rhs_flow.data(), nullptr, num_vertices));
        // model.update();
        // model.write("zdd_" + problem_name + "_" + std::to_string(convex_rhs)
        // +
        //             ".lp");
        // model.optimize();
    } catch (GRBException& e) {
        fmt::print("Error code = {}\n", e.getErrorCode());
        fmt::print("{}\n", e.getMessage());
    }
    catch (...) {
        fmt::print("Exception during optimization\n");
    }
}

void PricerSolverZdd::construct_lp_sol_from_rmp(
    const std::span<const double>&              lambda,
    const std::vector<std::shared_ptr<Column>>& columns) {
    auto& table = *(decision_diagram->getDiagram());
    // std::span aux_sets{columns->pdata, columns->len};
    std::fill(lp_x.begin(), lp_x.end(), 0.0);
    for (auto i = 0UL; i < columns.size(); ++i) {
        if (lambda[i] > EPS_SOLVER) {
            auto  counter = 0UL;
            auto* tmp = columns[i].get();
            // std::span aux_jobs{tmp->job_list->pdata, tmp->job_list->pdata};
            NodeId                        tmp_nodeid(decision_diagram->root());
            std::shared_ptr<SubNodeZdd<>> tmp_sub_node =
                table.node(tmp_nodeid).list[0];

            while (tmp_nodeid > 1) {
                Job* tmp_j{};

                if (counter < tmp->job_list.size()) {
                    tmp_j = tmp->job_list[counter];
                } else {
                    tmp_j = nullptr;
                }

                NodeZdd<>& tmp_node = table.node(tmp_nodeid);

                if (tmp_j == tmp_node.get_job()) {
                    lp_x[tmp_sub_node->high_edge_key] += lambda[i];
                    tmp_nodeid = tmp_node[1];
                    tmp_sub_node = tmp_sub_node->y;
                    counter++;
                } else {
                    lp_x[tmp_sub_node->low_edge_key] += lambda[i];
                    tmp_nodeid = tmp_node[0];
                    tmp_sub_node = tmp_sub_node->n;
                }
            }
        }
    }

    // ColorWriterEdge edge_writer(g, x);
    // ColorWriterVertex vertex_writer(g, table);
    // std::ofstream outf("min.gv");
    // boost::write_graphviz(outf, g, vertex_writer, edge_writer);
    // outf.close();
}

auto PricerSolverZdd::check_column(Column const* col) -> bool {
    auto        weight = 0UL;
    const auto& set = col->job_list;
    // std::span aux_jobs{set->pdata, set->len};
    auto&  table = *(decision_diagram->getDiagram());
    NodeId tmp_nodeid(decision_diagram->root());

    for (unsigned j = 0; j < set.size(); ++j) {
        Job* tmp_j = set[j];
        while (tmp_nodeid > 1) {
            NodeZdd<>& tmp_node = table.node(tmp_nodeid);

            if (tmp_j == tmp_node.get_job()) {
                tmp_nodeid = tmp_node[1];
                weight += 1;

                if (j + 1 != weight) {
                    return false;
                }
            } else {
                tmp_nodeid = tmp_node[0];
            }
        }
    }

    return (weight == set.size());
}

auto PricerSolverZdd::get_nb_edges() -> size_t {
    return num_edges(mip_graph);
}

auto PricerSolverZdd::get_nb_vertices() -> size_t {
    return decision_diagram->size();
}

auto PricerSolverZdd::print_num_paths() -> boost::multiprecision::cpp_int {
    return 0;
}

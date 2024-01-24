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

#include "PricerSolverBddForward.hpp"
#include <fmt/core.h>                            // for print
#include <array>                                 // for array
#include <range/v3/iterator/basic_iterator.hpp>  // for operator!=
#include <range/v3/view/join.hpp>                // for join_view
#include <range/v3/view/view.hpp>                // for operator|
#include "Instance.h"                            // for Instance
#include "Label.hpp"                             // for Label
#include "ModernDD/NodeBddStructure.hpp"         // for DdStructure
#include "ModernDD/NodeBddTable.hpp"             // for NodeTableEntity, Tab...
#include "PricerSolverBase.hpp"                  // for PricerSolverBase::ALIGN
#include "PricerSolverBdd.hpp"                   // for PricerSolverBdd

/**
 *  bdd solver pricersolver for the flow formulation
 */
PricerSolverBddSimple::PricerSolverBddSimple(const Instance& instance)
    : PricerSolverBdd(instance) {
    fmt::print("{0: <{1}}{2}\n", "Constructing BDD with evaluator:", ALIGN,
               "Forward Simple Evaluator");
    fmt::print("{0: <{1}}{2}\n", "Number of vertices BDD", ALIGN,
               get_nb_vertices());
    fmt::print("{0: <{1}}{2}\n", "Number of edges BDD", ALIGN, get_nb_edges());
}

auto PricerSolverBddSimple::pricing_algorithm(std::span<const double>& _pi)
    -> PricingSolution {
    evaluator.set_aux_pi(_pi);
    return get_decision_diagram().evaluate_forward(evaluator);
}

auto PricerSolverBddSimple::farkas_pricing(std::span<const double>& _pi)
    -> PricingSolution {
    farkas_evaluator.set_aux_pi(_pi);
    return get_decision_diagram().evaluate_backward(farkas_evaluator);
}

void PricerSolverBddSimple::compute_labels(std::span<const double>& _pi) {
    evaluator.set_aux_pi(_pi);
    reversed_evaluator.set_aux_pi(_pi);
    get_decision_diagram().compute_labels_forward(evaluator);
    get_decision_diagram().compute_labels_backward(reversed_evaluator);
}

auto PricerSolverBddSimple::evaluate_rc_arc(NodeBdd& n) -> double {
    auto& table = *(get_decision_diagram().getDiagram());

    auto& child = table.node(n[1]);
    return n.forward_label[0].get_f() + child.backward_label[0].get_f() +
           n.get_reduced_cost()[1];
}

// bool PricerSolverBddSimple::evaluate_nodes(double* pi) {
//     auto& table = *(get_decision_diagram().getDiagram());
//     compute_labels(pi);
//     auto reduced_cost = table.node(1).forward_label[0].get_f();
//     auto removed_edges = false;
//     auto nb_removed_edges_evaluate = 0;

//     /** check for each node the Lagrangian dual */
//     for (auto& it :
//          table | ranges::views::take(get_decision_diagram().topLevel() + 1) |
//              ranges::views ::drop(1) | ranges::views::reverse |
//              ranges::views::join) {
//         auto&  child = table.node(it[1]);
//         double result = it.forward_label[0].get_f() +
//                         child.backward_label[0].get_f() + it.reduced_cost[1];
//         auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
//         if (constLB + aux_nb_machines * reduced_cost + result >
//                 UB + RC_FIXING &&
//             (it.calc[1])) {
//             it.calc[1] = false;
//             add_nb_removed_edges();
//             removed_edges = true;
//             nb_removed_edges_evaluate++;
//         }
//     }

//     if (removed_edges) {
//         fmt::print("Number of edges removed by evaluate nodes {0:<{1}}\n",
//                    nb_removed_edges_evaluate, ALIGN_HALF);
//         fmt::print("Total number of edges removed {0:<{1}}\n",
//                    get_nb_removed_edges(), ALIGN_HALF);
//         fmt::print("Number of edges {0:<{1}}\n", get_nb_edges(), ALIGN_HALF);
//         remove_layers();
//         remove_edges();
//         bottom_up_filtering();
//         topdown_filtering();
//         cleanup_arcs();
//         construct_mipgraph();
//     }

//     return removed_edges;
// }

/**
 * bdd solver pricersolver for the flow formulation that takes care of the
 * consecutive jobs
 */
PricerSolverBddCycle::PricerSolverBddCycle(const Instance& instance)
    : PricerSolverBdd(instance) {
    fmt::print("{0: <{1}}{2}\n", "Constructing BDD with evaluator:", ALIGN,
               "Forward Cycle Evaluator");
    fmt::print("{0: <{1}}{2}\n", "Number of vertices BDD", ALIGN,
               get_nb_vertices());
    fmt::print("{0: <{1}}{2}\n", "Number of edges BDD", ALIGN, get_nb_edges());
}

auto PricerSolverBddCycle::pricing_algorithm(std::span<const double>& _pi)
    -> PricingSolution {
    evaluator.set_aux_pi(_pi);
    return get_decision_diagram().evaluate_forward(evaluator);
}

auto PricerSolverBddCycle::farkas_pricing(std::span<const double>& _pi)
    -> PricingSolution {
    farkas_evaluator.set_aux_pi(_pi);
    return get_decision_diagram().evaluate_backward(farkas_evaluator);
}

void PricerSolverBddCycle::compute_labels(std::span<const double>& _pi) {
    evaluator.set_aux_pi(_pi);
    reversed_evaluator.set_aux_pi(_pi);
    get_decision_diagram().compute_labels_forward(evaluator);
    get_decision_diagram().compute_labels_backward(reversed_evaluator);
}

auto PricerSolverBddCycle::evaluate_rc_arc(NodeBdd& n) -> double {
    auto& table = *(get_decision_diagram().getDiagram());
    auto* job = n.get_job();
    auto& child = table.node(n[1]);

    if (n.forward_label[0].prev_job_forward() != job &&
        child.backward_label[0].prev_job_backward() != job) {
        return n.forward_label[0].get_f() + child.backward_label[0].get_f() +
               n.get_reduced_cost()[1];
    }

    if (n.forward_label[0].prev_job_forward() == job &&
        child.backward_label[0].prev_job_backward() != job) {
        return n.forward_label[1].get_f() + child.backward_label[0].get_f() +
               n.get_reduced_cost()[1];
    }

    if (n.forward_label[0].prev_job_forward() != job &&
        child.backward_label[0].prev_job_backward() == job) {
        return n.forward_label[0].get_f() + child.backward_label[1].get_f() +
               n.get_reduced_cost()[1];
    }

    return n.forward_label[1].get_f() + child.backward_label[1].get_f() +
           n.get_reduced_cost()[1];
}

// bool PricerSolverBddCycle::evaluate_nodes(double* pi) {
//     auto& table = *(get_decision_diagram().getDiagram());
//     compute_labels(pi);
//     auto reduced_cost = table.node(1).forward_label[0].get_f();
//     auto removed_edges = false;
//     auto nb_removed_edges_evaluate = 0;

//     /** check for each node the Lagrangian dual */
//     for (auto& it :
//          table | ranges::views::take(get_decision_diagram().topLevel() + 1) |
//              ranges::views ::drop(1) | ranges::views::reverse |
//              ranges::views::join) {
//         auto* job = it.get_job();
//         auto  result{0.0};
//         auto& child = table.node(it[1]);

//         if (it.forward_label[0].prev_job_forward() != job &&
//             child.backward_label[0].prev_job_backward() != job) {
//             result = it.forward_label[0].get_f() +
//                      child.backward_label[0].get_f() + it.reduced_cost[1];
//         } else if (it.forward_label[0].prev_job_forward() == job &&
//                    child.backward_label[0].prev_job_backward() != job) {
//             result = it.forward_label[1].get_f() +
//                      child.backward_label[0].get_f() + it.reduced_cost[1];
//         } else if (it.forward_label[0].prev_job_forward() != job &&
//                    child.backward_label[0].prev_job_backward() == job) {
//             result = it.forward_label[0].get_f() +
//                      child.backward_label[1].get_f() + it.reduced_cost[1];
//         } else {
//             result = it.forward_label[1].get_f() +
//                      child.backward_label[1].get_f() + it.reduced_cost[1];
//         }

//         auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
//         if (constLB + aux_nb_machines * reduced_cost + result >
//                 UB + RC_FIXING &&
//             (it.calc[1])) {
//             it.calc[1] = false;
//             removed_edges = true;
//             add_nb_removed_edges();
//             nb_removed_edges_evaluate++;
//         }
//     }

//     if (removed_edges) {
//         fmt::print("Number of edges removed by evaluate nodes {0: <{1}}\n",
//                    nb_removed_edges_evaluate, ALIGN_HALF);
//         fmt::print("Total number of edges removed {0: <{1}}\n",
//                    get_nb_removed_edges(), ALIGN_HALF);
//         fmt::print("Number of edges {0: <{1}}\n", get_nb_edges(),
//         ALIGN_HALF); remove_layers(); remove_edges(); bottom_up_filtering();
//         topdown_filtering();
//         cleanup_arcs();
//         construct_mipgraph();
//     }

//     return removed_edges;
// }

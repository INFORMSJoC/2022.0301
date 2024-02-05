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

#include "PricerSolverZddBackward.hpp"
#include <fmt/core.h>                            // for print
#include <array>                                 // for array
#include <iostream>                              // for operator<<, basic_os...
#include <memory>                                // for __shared_ptr_access
#include <range/v3/iterator/basic_iterator.hpp>  // for operator!=, basic_it...
#include <range/v3/view/drop.hpp>                // for drop, drop_fn
#include <range/v3/view/join.hpp>                // for join_view, join_view...
#include <range/v3/view/subrange.hpp>            // for subrange
#include <range/v3/view/take.hpp>                // for take_view, take, tak...
#include <range/v3/view/view.hpp>                // for operator|, view_closure
#include <vector>                                // for vector
#include "Instance.h"                            // for Instance
#include "Job.h"                                 // for Job
#include "Label.hpp"                             // for Label
#include "ModernDD/NodeBddStructure.hpp"         // for DdStructure
#include "ModernDD/NodeBddTable.hpp"             // for NodeTableEntity, Tab...
#include "PricerSolverBase.hpp"                  // for PricerSolverBase::RC...
#include "ZddNode.hpp"                           // for SubNodeZdd, NodeZdd

/**
 *  bdd solver pricersolver for the flow formulation
 */

PricerSolverZddBackwardSimple::PricerSolverZddBackwardSimple(
    const Instance& instance)
    : PricerSolverZdd(instance) {
    fmt::print("Constructing ZDD with Backward Simple evaluator\n");
    fmt::print("number vertices ZDD = {}\n", get_nb_vertices());
    fmt::print("number edges ZDD = {}\n", get_nb_edges());
    evaluator = BackwardZddSimpleDouble();
    reversed_evaluator = ForwardZddSimpleDouble();
}

auto PricerSolverZddBackwardSimple::pricing_algorithm(
    std::span<const double>& _pi) -> PricingSolution {
    evaluator.initialize_pi(_pi);
    return decision_diagram->evaluate_backward(evaluator);
}

void PricerSolverZddBackwardSimple::compute_labels(
    std::span<const double>& _pi) {
    evaluator.initialize_pi(_pi);
    reversed_evaluator.initialize_pi(_pi);

    decision_diagram->compute_labels_backward(evaluator);
    decision_diagram->compute_labels_forward(reversed_evaluator);
}

auto PricerSolverZddBackwardSimple::evaluate_nodes(std::span<const double>& pi)
    -> bool {
    auto& table = *(decision_diagram->getDiagram());
    compute_labels(pi);
    double reduced_cost =
        table.node(decision_diagram->root()).list[0]->backward_label[0].get_f();

    nb_removed_edges = 0;

    // /** check for each node the Lagrangian dual */
    for (auto& it : table |
                        ranges::views::take(decision_diagram->topLevel() + 1) |
                        ranges::views ::drop(1) | ranges::views::join) {
        for (auto& iter : it.list) {
            int    w = iter->get_weight();
            Job*   job = it.get_job();
            double result = iter->forward_label[0].get_f() +
                            iter->y->backward_label[0].get_f() -
                            job->weighted_tardiness_start(w) + pi[job->job];
            auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
            if (constLB + aux_nb_machines * reduced_cost + result >
                    UB - 1 + RC_FIXING &&
                (iter->calc_yes)) {
                iter->calc_yes = false;
                nb_removed_edges++;
            }
        }
    }

    fmt::print("removed edges = {}\n", nb_removed_edges);

    return (nb_removed_edges > 0);
}
PricerSolverZddBackwardCycle::PricerSolverZddBackwardCycle(
    const Instance& instance)
    : PricerSolverZdd(instance) {
    fmt::print("Constructing ZDD with Backward ZddCycle evaluator\n");
    fmt::print("number vertices ZDD = {}\n", get_nb_vertices());
    fmt::print("number edges ZDD = {}\n", get_nb_edges());
    evaluator = BackwardZddCycleDouble();
    reversed_evaluator = ForwardZddCycleDouble();
}

auto PricerSolverZddBackwardCycle::pricing_algorithm(
    std::span<const double>& _pi) -> PricingSolution {
    evaluator.initialize_pi(_pi);
    return decision_diagram->evaluate_backward(evaluator);
}

void PricerSolverZddBackwardCycle::compute_labels(
    std::span<const double>& _pi) {
    evaluator.initialize_pi(_pi);
    reversed_evaluator.initialize_pi(_pi);

    decision_diagram->compute_labels_backward(evaluator);
    decision_diagram->compute_labels_forward(reversed_evaluator);
}

auto PricerSolverZddBackwardCycle::evaluate_nodes(std::span<const double>& pi)
    -> bool {
    auto& table = *(decision_diagram->getDiagram());
    compute_labels(pi);
    double reduced_cost =
        table.node(decision_diagram->root()).list[0]->backward_label[0].get_f();
    nb_removed_edges = 0;

    /** check for each node the Lagrangian dual */
    for (auto& it : table |
                        ranges::views::take(decision_diagram->topLevel() + 1) |
                        ranges::views ::drop(1) | ranges::views::join) {
        auto* job = it.get_job();
        if (job == nullptr) {
            continue;
        }

        for (auto& iter : it.list) {
            auto w = iter->get_weight();

            auto aux_nb_machines = static_cast<double>(convex_rhs - 1);
            if (iter->forward_label[0].prev_job_forward() != job &&
                iter->y->backward_label[0].prev_job_backward() != job) {
                double result = iter->forward_label[0].get_f() +
                                iter->y->backward_label[0].get_f() -
                                job->weighted_tardiness_start(w) + pi[job->job];
                if (constLB + aux_nb_machines * reduced_cost + result >
                        UB + RC_FIXING &&
                    (iter->calc_yes)) {
                    iter->calc_yes = false;
                    nb_removed_edges++;
                }
            } else if (iter->forward_label[0].prev_job_forward() == job &&
                       iter->y->backward_label[0].prev_job_backward() != job) {
                double result = iter->forward_label[1].get_f() +
                                iter->y->backward_label[0].get_f() -
                                job->weighted_tardiness_start(w) + pi[job->job];
                if (constLB + aux_nb_machines * reduced_cost + result >
                        UB + RC_FIXING &&
                    (iter->calc_yes)) {
                    iter->calc_yes = false;
                    nb_removed_edges++;
                }
            } else if (iter->forward_label[0].prev_job_forward() != job &&
                       iter->y->backward_label[0].prev_job_backward() == job) {
                double result = iter->forward_label[0].get_f() +
                                iter->y->backward_label[1].get_f() -
                                job->weighted_tardiness_start(w) + pi[job->job];
                if (constLB + aux_nb_machines * reduced_cost + result >
                        UB + RC_FIXING &&
                    (iter->calc_yes)) {
                    iter->calc_yes = false;
                    nb_removed_edges++;
                }
            } else {
                double result = iter->forward_label[1].get_f() +
                                iter->y->backward_label[1].get_f() -
                                job->weighted_tardiness_start(w) + pi[job->job];
                if (constLB + aux_nb_machines * reduced_cost + result >
                        UB + RC_FIXING &&
                    (iter->calc_yes)) {
                    iter->calc_yes = false;
                    nb_removed_edges++;
                }
            }
        }
    }

    fmt::print("removed edges = {}\n", nb_removed_edges);

    return (nb_removed_edges > 0);
}

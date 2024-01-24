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

#include "PricerSolverBddBackward.hpp"
#include <fmt/core.h>            // for print
#include <span>                  // for span
#include "DebugLvl.hpp"          // for debug_lvl
#include "Instance.h"            // for Instance
#include "NodeBdd.hpp"           // for NodeBdd
#include "PricerSolverBase.hpp"  // for PricerSolverBase::ALIGN
#include "PricerSolverBdd.hpp"   // for PricerSolverBdd
// #include "orutils/util.h"        // for dbg_lvl

/**
 * backward bdd pricersolver for the flow formulation that takes care of the
 * consecutive jobs
 */

PricerSolverBddBackwardSimple::PricerSolverBddBackwardSimple(
    const Instance& instance)
    : PricerSolverBdd(instance) {
    if (debug_lvl(0)) {
        fmt::print("{0: <{1}}{2}\n", "Constructing BDD with evaluator:", ALIGN,
                   "Backward Simple Evaluator");
        fmt::print(R"({0: <{1}}{2}
)",
                   "Number of vertices BDD", ALIGN, get_nb_vertices());
        fmt::print("{0: <{1}}{2}\n", "Number of edges BDD", ALIGN,
                   get_nb_edges());
    }
    evaluator.set_table(&(*(get_decision_diagram().getDiagram())));
}


auto PricerSolverBddBackwardSimple::pricing_algorithm(
    std::span<const double>& pi) -> PricingSolution {
    evaluator.set_aux_pi(pi);
    return get_decision_diagram().evaluate_backward(evaluator);
}


auto PricerSolverBddBackwardSimple::farkas_pricing(std::span<const double>& _pi)
    -> PricingSolution {
    farkas_evaluator.set_aux_pi(_pi);
    return get_decision_diagram().evaluate_backward(farkas_evaluator);
}

void PricerSolverBddBackwardSimple::compute_labels(
    std::span<const double>& _pi) {
    evaluator.set_aux_pi(_pi);
    reversed_evaluator.set_aux_pi(_pi);
    get_decision_diagram().compute_labels_backward(evaluator);
    get_decision_diagram().compute_labels_forward(reversed_evaluator);
}

auto PricerSolverBddBackwardSimple::evaluate_rc_arc(NodeBdd& node) -> double {
    auto& table = *(get_decision_diagram().getDiagram());
    auto& child = table.node(node[1]);
    return node.forward_label[0].get_f() + child.backward_label[0].get_f() +
           node.get_reduced_cost()[1];
}

/**
 * Simple backward bdd pricersolver for the flow formulation
 */
PricerSolverBddBackwardCycle::PricerSolverBddBackwardCycle(
    const Instance& instance)
    : PricerSolverBdd(instance) {
    if (debug_lvl(0)) {
        fmt::print("{0: <{1}}{2}\n", "Constructing BDD with evaluator:", ALIGN,
                   "Backward Cycle Evaluator");
        fmt::print("{0: <{1}}{2}\n", "Number of vertices BDD", ALIGN,
                   get_nb_vertices());
        fmt::print("{0: <{1}}{2}\n", "Number of edges BDD", ALIGN,
                   get_nb_edges());
    }
    evaluator.set_table(&(*(get_decision_diagram().getDiagram())));
}


auto PricerSolverBddBackwardCycle::pricing_algorithm(
    std::span<const double>& _pi) -> PricingSolution {
    evaluator.set_aux_pi(_pi);
    return get_decision_diagram().evaluate_backward(evaluator);
}


auto PricerSolverBddBackwardCycle::farkas_pricing(std::span<const double>& _pi)
    -> PricingSolution {
    farkas_evaluator.set_aux_pi(_pi);
    return get_decision_diagram().evaluate_backward(farkas_evaluator);
}

void PricerSolverBddBackwardCycle::compute_labels(
    std::span<const double>& _pi) {
    evaluator.set_aux_pi(_pi);
    reversed_evaluator.set_aux_pi(_pi);
    get_decision_diagram().compute_labels_backward(evaluator);
    get_decision_diagram().compute_labels_forward(reversed_evaluator);
}

auto PricerSolverBddBackwardCycle::evaluate_rc_arc(NodeBdd& n) -> double {
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

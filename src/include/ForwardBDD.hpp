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

#ifndef FORWARD_BDD_HPP
#define FORWARD_BDD_HPP
#include <limits>                         // for numeric_limits
#include <range/v3/action/remove_if.hpp>  // for remove_if
#include <span>                           // for span
#include "ModernDD/NodeBddEval.hpp"       // for Eval
#include "NodeBdd.hpp"                    // for NodeBdd
#include "PricingSolution.hpp"            // for PricingSolution

class ForwardBddBase : public Eval<NodeBdd, PricingSolution> {
    std::span<const double> aux_pi;

   public:
    ForwardBddBase() = default;
    ForwardBddBase(const ForwardBddBase& src) = default;
    ForwardBddBase(ForwardBddBase&&) noexcept = default;
    ~ForwardBddBase() override = default;

    auto operator=(const ForwardBddBase&) -> ForwardBddBase& = default;
    auto operator=(ForwardBddBase&&) noexcept -> ForwardBddBase& = default;

    void set_aux_pi(std::span<const double>& _pi) { aux_pi = _pi; }

    [[nodiscard]] auto get_aux_pi() const -> std::span<const double> {
        return aux_pi;
    }

    void initialize_node(NodeBdd& n) const override = 0;
    void initialize_root_node(NodeBdd& n) const override = 0;
    void evalNode(NodeBdd& n) const override = 0;

    auto get_objective(NodeBdd& n) const -> PricingSolution override {
        PricingSolution sol(0.0);
        auto*           ptr_node = &(n.forward_label[0]);
        auto*           table_tmp = Eval<NodeBdd, PricingSolution>::get_table();

        while (ptr_node->get_previous() != nullptr) {
            auto* aux_prev_node = ptr_node->get_previous();
            auto& node = table_tmp->node(aux_prev_node->get_node_id());
            auto* aux_job = node.get_job();
            sol.C_max += aux_job->processing_time;
            sol.push_job_back(aux_job, node.get_weight(),
                              node.get_reduced_cost()[1]);
            ptr_node = aux_prev_node;
        }

        sol.reverse_jobs();

        return sol;
    }
};

class ForwardBddCycle : public ForwardBddBase {
   public:
    ForwardBddCycle() = default;
    ForwardBddCycle(const ForwardBddCycle& src) = default;
    ForwardBddCycle(ForwardBddCycle&&) noexcept = default;
    ~ForwardBddCycle() override = default;

    auto operator=(const ForwardBddCycle&) -> ForwardBddCycle& = default;
    auto operator=(ForwardBddCycle&&) noexcept -> ForwardBddCycle& = default;

    void initialize_node(NodeBdd& n) const override {
        if (n.get_weight() == 0) {
            n.forward_label[0].forward_update(0, nullptr, false);
            n.forward_label[1].reset();
        } else {
            for (auto& it : n.forward_label) {
                it.reset();
            }
        }
    }

    void initialize_root_node(NodeBdd& n) const override {
        n.forward_label[0].get_f() = 0;
        n.forward_label[1].set_f(std::numeric_limits<double>::max());
    }

    void evalNode(NodeBdd& n) const override {
        auto*       tmp_j = n.get_job();
        auto*       table_tmp = Eval<NodeBdd, PricingSolution>::get_table();
        auto&       p0 = table_tmp->node(n[0]);
        auto&       p1 = table_tmp->node(n[1]);
        auto dual = ForwardBddBase::get_aux_pi();

        n.reset_reduced_costs();
        n.adjust_reduced_costs(dual[tmp_j->job], true);
        n.adjust_reduced_costs(0.0, false);

        // for (auto& list : n.coeff_list) {
        //     list |= ranges::actions::remove_if([&](auto& it) {
        //         auto aux = it.lock();
        //         if (aux) {
        //             n.adjust_reduced_costs(
        //                 aux->get_coeff() * dual[aux->get_row()],
        //                 aux->get_high());
        //             return false;
        //         }
        //         return true;
        //     });
        // }

        /**
         * High edge calculation
         */
        auto* prev = n.forward_label[0].prev_job_forward();
        auto* aux1 = p1.forward_label[0].prev_job_forward();

        if (prev != tmp_j) {
            auto g = n.forward_label[0].get_f() + n.get_reduced_cost()[1];
            if (g < p1.forward_label[0].get_f()) {
                if (aux1 != tmp_j) {
                    p1.forward_label[1].forward_update(p1.forward_label[0]);
                }
                p1.forward_label[0].forward_update(g, &(n.forward_label[0]),
                                                   true);
            } else if ((g < p1.forward_label[1].get_f()) && (aux1 != tmp_j)) {
                p1.forward_label[1].forward_update(g, &(n.forward_label[0]),
                                                   true);
            }
        } else {
            auto g = n.forward_label[1].get_f() + n.get_reduced_cost()[1];
            prev = n.forward_label[1].prev_job_forward();

            if (g < p1.forward_label[0].get_f()) {
                if (aux1 != tmp_j) {
                    p1.forward_label[1].forward_update(p1.forward_label[0]);
                }
                p1.forward_label[0].forward_update(g, &(n.forward_label[1]),
                                                   true);
            } else if ((g < p1.forward_label[1].get_f()) && (aux1 != tmp_j)) {
                p1.forward_label[1].forward_update(g, &(n.forward_label[1]),
                                                   true);
            }
        }

        /**
         * Low edge calculation
         */
        aux1 = p0.forward_label[0].prev_job_forward();
        auto g = n.forward_label[0].get_f() + n.get_reduced_cost()[0];
        auto g1 = n.forward_label[1].get_f() + n.get_reduced_cost()[0];
        if (g < p0.forward_label[0].get_f()) {
            if (prev != aux1) {
                p0.forward_label[1].forward_update(p0.forward_label[0]);
            }
            p0.forward_label[0].forward_update(g, n.forward_label[0]);
            if (g1 < p0.forward_label[1].get_f()) {
                p0.forward_label[1].forward_update(g1, n.forward_label[1]);
            }
        } else if ((g < p0.forward_label[1].get_f()) && (aux1 != prev)) {
            p0.forward_label[1].forward_update(g, n.forward_label[0]);
        } else if ((g1 < p0.forward_label[1].get_f())) {
            p0.forward_label[1].forward_update(g1, n.forward_label[1]);
        }
    }
};

class ForwardBddSimple : public ForwardBddBase {
   public:
    ForwardBddSimple() = default;

    ForwardBddSimple(ForwardBddSimple&&) noexcept = default;
    ForwardBddSimple(const ForwardBddSimple&) = default;
    ~ForwardBddSimple() override = default;

    auto operator=(const ForwardBddSimple&) -> ForwardBddSimple& = default;
    auto operator=(ForwardBddSimple&&) noexcept -> ForwardBddSimple& = default;

    void initialize_node(NodeBdd& n) const override {
        if (n.get_weight() == 0) {
            n.forward_label[0].forward_update(0, nullptr, false);
        } else {
            n.forward_label[0].reset();
        }
    }

    void initialize_root_node(NodeBdd& n) const override {
        n.forward_label[0].get_f() = 0;
    }

    void evalNode(NodeBdd& n) const override {
        auto* table_tmp = Eval<NodeBdd, PricingSolution>::get_table();
        auto& p0 = table_tmp->node(n[0]);
        auto& p1 = table_tmp->node(n[1]);
        n.reset_reduced_costs();
        auto dual = ForwardBddBase::get_aux_pi();

        for (auto& list : n.get_coeff_list()) {
            list |= ranges::actions::remove_if([&](auto& it) {
                auto aux = it.lock();
                if (aux) {
                    n.adjust_reduced_costs(
                        aux->get_coeff() * dual[aux->get_row()],
                        aux->get_high());
                    return false;
                }
                return true;
            });
        }

        /**
         * High edge calculation
         */
        auto g = n.forward_label[0].get_f() + n.get_reduced_cost()[1];
        if (g < p1.forward_label[0].get_f()) {
            p1.forward_label[0].forward_update(g, &(n.forward_label[0]), true);
        }

        /**
         * Low edge calculation
         */
        g = n.forward_label[0].get_f() + n.get_reduced_cost()[0];
        if (g < p0.forward_label[0].get_f()) {
            p0.forward_label[0].forward_update(g, n.forward_label[0]);
        }
    }
};

#endif  // FORWARD_BDD_HPP

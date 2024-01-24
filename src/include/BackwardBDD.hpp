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

#ifndef BACKWARD_BDD_HPP
#define BACKWARD_BDD_HPP
#include <array>
#include <range/v3/action/remove_if.hpp>  // for remove_if
#include <span>                           //for span
#include "BddCoeff.hpp"
#include "ModernDD/NodeBddEval.hpp"  // for Eval
#include "NodeBdd.hpp"               // for NodeBdd
#include "PricingSolution.hpp"       // for PricingSolution

class BackwardBddBase : public Eval<NodeBdd, PricingSolution> {
    std::span<const double> aux_pi{};

   public:
    BackwardBddBase() = default;

    void set_aux_pi(std::span<const double>& _pi) { aux_pi = _pi; }

    [[nodiscard]] auto get_aux_pi() const -> std::span<const double> {
        return aux_pi;
    }

    auto get_objective(NodeBdd& n) const -> PricingSolution override {
        PricingSolution sol(0.0);
        auto*           aux_label = &(n.backward_label[0]);
        auto*           table_tmp = Eval<NodeBdd, PricingSolution>::get_table();

        do {
            auto& tmp_node_id = aux_label->get_node_id();
            if (aux_label->get_high()) {
                auto& node = table_tmp->node(tmp_node_id);
                sol.push_job_back(node.get_job(), node.get_reduced_cost()[1]);
            }

            aux_label = aux_label->get_previous();
        } while (aux_label != nullptr);

        return sol;
    }

    void initialize_node(NodeBdd& n) const override = 0;
    void initialize_root_node(NodeBdd& n) const override = 0;
    void evalNode(NodeBdd& n) const override = 0;

    BackwardBddBase(const BackwardBddBase&) = default;
    BackwardBddBase(BackwardBddBase&&) noexcept = default;
    ~BackwardBddBase() override = default;

    auto operator=(const BackwardBddBase&) -> BackwardBddBase& = default;
    auto operator=(BackwardBddBase&&) noexcept -> BackwardBddBase& = default;
};

class BackwardBddSimple : public BackwardBddBase {
   public:
    BackwardBddSimple() = default;

    void evalNode(NodeBdd& n) const override {
        auto* table_tmp = Eval<NodeBdd, PricingSolution>::get_table();
        auto& p0_tmp = table_tmp->node(n[0]);
        auto& p1_tmp = table_tmp->node(n[1]);

        n.reset_reduced_costs();

        auto dual = BackwardBddBase::get_aux_pi();

        for (auto& list : n.get_coeff_list()) {
            list |= ranges::actions::remove_if([&](auto it) {
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

        auto obj0 = p0_tmp.backward_label[0].get_f() + n.get_reduced_cost()[0];
        auto obj1 = p1_tmp.backward_label[0].get_f() + n.get_reduced_cost()[1];

        if (obj0 > obj1) {
            n.backward_label[0].backward_update(&(p1_tmp.backward_label[0]),
                                                obj1, true);
        } else {
            n.backward_label[0].backward_update(&(p0_tmp.backward_label[0]),
                                                obj0, false);
        }
    }

    void initialize_node(NodeBdd& n) const override {
        n.backward_label[0].reset();
    }

    void initialize_root_node(NodeBdd& n) const override {
        n.backward_label[0].get_f() = 0;
    }

    BackwardBddSimple(const BackwardBddSimple&) = default;
    BackwardBddSimple(BackwardBddSimple&&) noexcept = default;
    ~BackwardBddSimple() override = default;
    auto operator=(const BackwardBddSimple&) -> BackwardBddSimple& = default;
    auto operator=(BackwardBddSimple&&) noexcept
        -> BackwardBddSimple& = default;
};

class BackwardBddCycle : public BackwardBddBase {
   public:
    BackwardBddCycle() = default;

    void evalNode(NodeBdd& n) const override {
        auto* tmp_j = n.get_job();
        auto* table_tmp = Eval<NodeBdd, PricingSolution>::get_table();
        auto& p0_tmp = table_tmp->node(n[0]);
        auto& p1_tmp = table_tmp->node(n[1]);
        auto dual = BackwardBddBase::get_aux_pi();

        n.reset_reduced_costs();

        n.adjust_reduced_costs(dual[tmp_j->job], true);
        n.adjust_reduced_costs(0.0, false);

        // for (auto& list : n.coeff_list) {
        //     list |= ranges::actions::remove_if([&](auto it) {
        //         auto aux = it.lock();
        //         if (aux) {
        //             n.adjust_reduced_costs(
        //                 aux->get_coeff() * dual[aux->get_row()],
        //                 aux->get_high());
        //             return false;
        //         } else {
        //             return true;
        //         }
        //     });
        // }

        auto* prev_job{p1_tmp.backward_label[0].prev_job_backward()};

        n.backward_label[0].backward_update(
            &(p0_tmp.backward_label[0]),
            p0_tmp.backward_label[0].get_f() + n.get_reduced_cost()[0], false);
        n.backward_label[1].backward_update(
            &(p0_tmp.backward_label[1]),
            p0_tmp.backward_label[1].get_f() + n.get_reduced_cost()[0], false);

        if (prev_job != tmp_j) {
            auto obj1{p1_tmp.backward_label[0].get_f() +
                      n.get_reduced_cost()[1]};
            auto obj2{p1_tmp.backward_label[1].get_f() +
                      n.get_reduced_cost()[1]};

            if (obj1 < n.backward_label[0].get_f()) {
                if (tmp_j != n.backward_label[0].prev_job_backward()) {
                    n.backward_label[1].backward_update(
                        &(p0_tmp.backward_label[0]),
                        p0_tmp.backward_label[0].get_f() +
                            n.get_reduced_cost()[0],
                        false);
                }

                n.backward_label[0].backward_update(&(p1_tmp.backward_label[0]),
                                                    obj1, true);
            } else if (obj1 < n.backward_label[1].get_f() &&
                       tmp_j != n.backward_label[0].prev_job_backward()) {
                n.backward_label[1].backward_update(&(p1_tmp.backward_label[0]),
                                                    obj1, true);
            } else if (obj2 < n.backward_label[1].get_f() &&
                       tmp_j != n.backward_label[0].prev_job_backward()) {
                n.backward_label[1].backward_update(&(p1_tmp.backward_label[1]),
                                                    obj2, true);
            }
        } else {
            auto obj1 =
                p1_tmp.backward_label[1].get_f() + n.get_reduced_cost()[1];

            if (obj1 < n.backward_label[0].get_f()) {
                if (tmp_j != n.backward_label[0].prev_job_backward()) {
                    n.backward_label[1].backward_update(
                        &(p0_tmp.backward_label[0]),
                        p0_tmp.backward_label[0].get_f() +
                            n.get_reduced_cost()[0],
                        false);
                }

                n.backward_label[0].backward_update(&(p1_tmp.backward_label[1]),
                                                    obj1, true);
            } else if (obj1 < n.backward_label[1].get_f() &&
                       tmp_j != n.backward_label[0].prev_job_backward()) {
                n.backward_label[1].backward_update(&(p1_tmp.backward_label[1]),
                                                    obj1, true);
            }
        }
    }

    void initialize_node(NodeBdd& n) const override {
        n.backward_label[0].reset();
    }

    void initialize_root_node(NodeBdd& n) const override {
        n.backward_label[0].get_f() = 0.0;
    }

    BackwardBddCycle(const BackwardBddCycle&) = default;
    BackwardBddCycle(BackwardBddCycle&&) noexcept = default;
    auto operator=(const BackwardBddCycle&) -> BackwardBddCycle& = default;
    auto operator=(BackwardBddCycle&&) noexcept -> BackwardBddCycle& = default;
    ~BackwardBddCycle() override = default;
};

#endif  // BACKWARD_BDD_HPP

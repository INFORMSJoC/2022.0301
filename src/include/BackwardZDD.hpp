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

#ifndef BACKWARD_ZDD_HPP
#define BACKWARD_ZDD_HPP
#include <array>                     // for array
#include <cassert>                   // for assert
#include <cstddef>                   // for size_t
#include <memory>                    // for shared_ptr, __shared_ptr_access
#include <span>                      // for span
#include <vector>                    // for vector
#include "Job.h"                     // for bool_diff_Fij, Job
#include "Label.hpp"                 // for Label
#include "ModernDD/NodeBddEval.hpp"  // for Eval
#include "PricingSolution.hpp"       // for PricingSolution
#include "ZddNode.hpp"  // for SubNodeZdd, NodeZdd, compare_sub_nodes

template <typename T = double>
class BackwardZDDBase : public Eval<NodeZdd<T>, PricingSolution> {
   private:
    std::span<const T> pi;

   public:
    BackwardZDDBase() = default;
    BackwardZDDBase(const BackwardZDDBase<T>&) = default;
    BackwardZDDBase(BackwardZDDBase<T>&&) noexcept = default;
    ~BackwardZDDBase() override = default;

    auto operator=(const BackwardZDDBase<T>&) -> BackwardZDDBase<T>& = default;
    auto operator=(BackwardZDDBase<T>&&) noexcept
        -> BackwardZDDBase<T>& = default;

    void initialize_pi(std::span<const T> _pi) { pi = _pi; }
    [[nodiscard]] auto get_pi() const -> std::span<const T> { return pi; }

    void initialize_node(NodeZdd<T>& n) const override = 0;
    void initialize_root_node(NodeZdd<T>& n) const override = 0;
    void evalNode(NodeZdd<T>& n) const override = 0;
};

template <typename T = double>
class BackwardZddSimple : public BackwardZDDBase<T> {
   public:
    BackwardZddSimple() = default;

    void evalNode(NodeZdd<T>& n) const override {
        auto* tmp_j = n.get_job();
        assert(tmp_j != nullptr);

        for (auto& it : n.list) {
            auto weight = it->weight;
            auto p0 = it->n;
            auto p1 = it->y;
            auto result = tmp_j->weighted_tardiness_start(weight);

            auto obj0 = p0->backward_label[0].get_f();
            auto obj1 = p1->backward_label[0].get_f() + result;

            if (obj0 > obj1) {
                it->backward_label[0].backward_update(obj1, true);
            } else {
                it->backward_label[0].backward_update(obj0, false);
            }
        }
    }

    void initialize_node(NodeZdd<T>& n) const override {
        for (auto& it : n.list) {
            it->backward_label[0].reset();
        }
    }

    void initialize_root_node(NodeZdd<T>& n) const override {
        for (auto& it : n.list) {
            it->backward_label[0].get_f() = 0.0;
        }
    }

    auto get_objective(NodeZdd<T>& n) const -> PricingSolution override {
        PricingSolution sol(0.0);

        auto m = std::min_element(n.list.begin(), n.list.end(),
                                  compare_sub_nodes<T>);
        auto ptr_node = (*m);
        auto aux_job = n.get_job();

        // NodeZdd<T> *aux_node = &n;
        // Job *aux_job =  n.get_job();
        auto dual = this->get_pi();
        while (aux_job) {
            if (ptr_node->backward_label[0].get_high()) {
                sol.push_job_back(aux_job, dual[aux_job->job]);
                ptr_node = ptr_node->y;
                aux_job = ptr_node->get_job();
            } else {
                ptr_node = ptr_node->n;
                aux_job = ptr_node->get_job();
            }
        }

        return sol;
    }
};

template <typename T = double>
class BackwardZddCycle : public BackwardZDDBase<T> {
   public:
    BackwardZddCycle() = default;

    void evalNode(NodeZdd<T>& n) const override {
        auto* tmp_j = n.get_job();
        std::span<const T> dual = this->get_pi();

        for (auto& it : n.list) {
            auto weight{it->get_weight()};
            auto p0{it->n};
            auto p1{it->y};
            T    result{tmp_j->weighted_tardiness_start(weight) -
                     dual[tmp_j->job]};

            auto* prev_job{p1->backward_label[0].prev_job_backward()};

            it->backward_label[0].backward_update(&(p0->backward_label[0]));
            it->backward_label[1].backward_update(&(p0->backward_label[1]));
            auto diff = bool_diff_Fij(weight, prev_job, tmp_j);
            auto diff1 = bool_diff_Fij(
                weight, p1->backward_label[0].prev_job_backward(), tmp_j);

            if (prev_job != tmp_j && diff) {
                T obj1{p1->backward_label[0].get_f() + result};
                T obj2{p1->backward_label[1].get_f() + result};

                if (obj1 < it->backward_label[0].get_f()) {
                    if (tmp_j != it->backward_label[0].prev_job_backward()) {
                        it->backward_label[1].backward_update(
                            &(p0->backward_label[0]));
                    }

                    it->backward_label[0].backward_update(
                        &(p1->backward_label[0]), obj1, true);
                } else if (obj1 < it->backward_label[1].get_f() &&
                           tmp_j != it->backward_label[0].prev_job_backward() &&
                           diff1) {
                    it->backward_label[1].backward_update(
                        &(p1->backward_label[0]), obj1, true);
                } else if (obj2 < it->backward_label[1].get_f() &&
                           tmp_j != it->backward_label[0].prev_job_backward()) {
                    it->backward_label[1].backward_update(
                        &(p1->backward_label[1]), obj2, true);
                }
            } else {
                T obj1 = p1->backward_label[1].get_f() + result;

                if (obj1 < it->backward_label[0].get_f()) {
                    if (tmp_j != it->backward_label[0].prev_job_backward()) {
                        it->backward_label[1].backward_update(
                            &(p0->backward_label[0]));
                    }

                    it->backward_label[0].backward_update(
                        &(p1->backward_label[1]), obj1, true);
                } else if (obj1 < it->backward_label[1].get_f() &&
                           tmp_j != it->backward_label[0].prev_job_backward()) {
                    it->backward_label[1].backward_update(
                        &(p1->backward_label[1]), obj1, true);
                }
            }
        }
    }

    void initialize_node(NodeZdd<T>& n) const override {
        for (auto& it : n.list) {
            it->backward_label[0].reset();
        }
    }

    void initialize_root_node(NodeZdd<T>& n) const override {
        for (auto& it : n.list) {
            it->backward_label[0].get_f() = 0.0;
        }
    }

    auto get_objective(NodeZdd<>& n) const -> PricingSolution override {
        auto  sol = PricingSolution(0.0);
        auto  m = std::min_element(n.list.begin(), n.list.end(),
                                   compare_sub_nodes<T>);
        auto* aux_label = &((*m)->backward_label[0]);
        std::span<const T>  dual = this->get_pi();

        while (aux_label) {
            if (aux_label->get_high()) {
                Job* aux_job = aux_label->get_job();
                sol.push_job_back(aux_job, dual[aux_job->job]);
            }

            aux_label = aux_label->get_previous();
        }

        return sol;
    }
};

#endif  // BACKWARD_ZDD_HPP

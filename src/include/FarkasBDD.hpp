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

#ifndef FARKASZDD_H
#define FARKASZDD_H

#include "BackwardBDD.hpp"           // for BackwardBddBase
#include "ModernDD/NodeBddEval.hpp"  // for Eval
#include "NodeBdd.hpp"               // for NodeBdd
#include "PricingSolution.hpp"       // for PricingSolution

class BackwardBddFarkas : public BackwardBddBase {
   public:
    BackwardBddFarkas() = default;

    void evalNode(NodeBdd& n) const override {
        n.reset_reduced_costs_farkas();

        auto dual = BackwardBddBase::get_aux_pi();
        for (auto& list : n.get_coeff_list()) {
            for (auto& it : list) {
                auto aux = it.lock();
                if (aux) {
                    n.adjust_reduced_costs(
                        aux->get_coeff() * dual[aux->get_row()],
                        aux->get_high());
                }
            }
        }

        auto* table_tmp = Eval<NodeBdd, PricingSolution>::get_table();
        auto& p0 = table_tmp->node(n[0]);
        auto& p1 = table_tmp->node(n[1]);

        auto obj0 = p0.backward_label[0].get_f() + n.get_reduced_cost()[0];
        auto obj1 = p1.backward_label[0].get_f() + n.get_reduced_cost()[1];

        if (obj0 > obj1) {
            n.backward_label[0].backward_update(&(p1.backward_label[0]), obj1,
                                                true);
        } else {
            n.backward_label[0].backward_update(&(p0.backward_label[0]), obj0,
                                                false);
        }
    }

    void initialize_node(NodeBdd& n) const override {
        n.backward_label[0].reset();
    }

    void initialize_root_node(NodeBdd& n) const override {
        n.backward_label[0].get_f() = 0.0;
    }
};
#endif  // FARKASZDD_H
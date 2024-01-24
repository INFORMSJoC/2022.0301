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

#ifndef CARDINALITYPATHS_H
#define CARDINALITYPATHS_H
#include <array>                             // for array
#include <boost/multiprecision/cpp_int.hpp>  // for cpp_int
#include <cstddef>                           // for size_t
#include "ModernDD/NodeBddEval.hpp"          // for Eval
#include "ModernDD/NodeBddTable.hpp"         // for NodeTableEntity
#include "ModernDD/NodeId.hpp"               // for NodeId
#include "NodeBdd.hpp"                       // for NodeBdd

class CardinalityPaths : public Eval<NodeBdd, boost::multiprecision::cpp_int> {
    using cpp_int = boost::multiprecision::cpp_int;

   public:
    CardinalityPaths() = default;

    auto get_objective(NodeBdd& n) const -> cpp_int override {
        return n.get_nb_paths();
    }

    void initialize_root_node(NodeBdd& n) const override { n.reset_nb_paths(); }

    void initialize_node(NodeBdd& n) const override { n.reset_nb_paths(); }

    void evalNode(NodeBdd& n) const override {
        auto* table_tmp = Eval<NodeBdd, cpp_int>::get_table();
        if (n[0] == 1) {
            n.update_nb_paths();
        } else if (n[0] > 1) {
            n.update_nb_paths(table_tmp->node(n[0]).get_nb_paths());
        }

        if (n[1] == 1) {
            n.update_nb_paths();
        } else {
            n.update_nb_paths(table_tmp->node(n[1]).get_nb_paths());
        }
    }
};

class BackwardDistance : public Eval<NodeBdd, std::array<int, 2>> {
   public:
    BackwardDistance() = default;

    auto get_objective([[maybe_unused]] NodeBdd& n) const
        -> std::array<int, 2> override {
        return {0, 0};
    }

    void initialize_root_node([[maybe_unused]] NodeBdd& n) const override {}

    void initialize_node([[maybe_unused]] NodeBdd& n) const override {}

    void evalNode(NodeBdd& n) const override {
        auto* table_tmp = Eval<NodeBdd, std::array<int, 2>>::get_table();
        for (size_t i = 0UL; i < 2; ++i) {
            auto& cur_node = table_tmp->node(n[i]);
            n.update_backward_distance(cur_node.get_backward_distance(),
                                       i == 1);
        }
    }
};
#endif  // CARDINALITYPATHS_H
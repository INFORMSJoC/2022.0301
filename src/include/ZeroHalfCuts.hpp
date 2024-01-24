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

#ifndef ZEROHALFCUTS_H
#define ZEROHALFCUTS_H

#include <gurobi_c++.h>               // for GRBVar, GRBEnv, GRBModel
#include <cstddef>                    // for size_t
#include <memory>                     // for shared_ptr, unique_ptr
#include <vector>                     // for vector
#include "ModernDD/NodeBddTable.hpp"  // for NodeTableEntity
#include "ModernDD/NodeId.hpp"        // for NodeId
#include "NodeBdd.hpp"                // for NodeBdd
class ConstraintGeneric;
class ReformulationModel;
class ZeroHalfCuts {
   public:
    ZeroHalfCuts(size_t                    _nb_jobs,
                 size_t                    _nb_machines,
                 ReformulationModel*       _rmp_model,
                 NodeId const&             _root,
                 NodeTableEntity<NodeBdd>* _table);
    ZeroHalfCuts(ZeroHalfCuts&&) = default;
    ZeroHalfCuts(const ZeroHalfCuts&) = delete;
    ~ZeroHalfCuts() = default;

    auto operator=(ZeroHalfCuts&&) -> ZeroHalfCuts& = default;
    auto operator=(const ZeroHalfCuts&) -> ZeroHalfCuts& = delete;

    auto add_cuts() -> bool;
    auto get_cut_list() -> std::vector<std::shared_ptr<ConstraintGeneric>>;

    void generate_cuts();

    static constexpr int ALIGN = 40;

   private:
    std::unique_ptr<GRBEnv>                         env;
    std::unique_ptr<GRBModel>                       model;
    size_t                                          nb_jobs;
    size_t                                          nb_machines;
    ReformulationModel*                             rmp_model;
    NodeId                                          root;
    NodeTableEntity<NodeBdd>*                       table;
    std::vector<NodeId>                             node_ids{};
    std::vector<NodeId>                             node_ids_lift{};
    std::vector<GRBVar>                             jobs_var;
    GRBVar                                          q;
    std::vector<std::shared_ptr<ConstraintGeneric>> cut_list;

    static constexpr double EPS_CUT = 1e-6;
    static constexpr double HALF = 2.0;
    static constexpr double TIMELIMIT = 50.0;

    void generate_model();
    void init_table();
    void init_coeff_cut();
    void init_coeff_node(NodeBdd* node);
    void construct_cut();
    void lift_operator();
    void dfs(const NodeId& v);
    void dfs_lift(const NodeId& v);
    void dfs_evaluate(const NodeId& v);
};

#endif  // ZEROHALFCUTS_H

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
#include "BranchBoundTree.hpp"
#include <memory>                       // for make_unique, unique_ptr
#include <utility>                      // for move
#include "BranchNode.hpp"               // for BranchNodeBase
#include "NodeData.h"                   // for NodeData
#include "Parms.h"                      // for Parms, bb_bfs_strategy, bb_br...
#include "Solution.hpp"                 // for Sol
#include "branch-and-bound/bfstree.h"   // for BFSTree
#include "branch-and-bound/brfstree.h"  // for BrFSTree
#include "branch-and-bound/cbfstree.h"  // for CBFSTree
#include "branch-and-bound/dfstree.h"   // for DFSTree

BranchBoundTree::BranchBoundTree(std::unique_ptr<NodeData> _data,
                                 int                       _probType,
                                 bool                      _isIntProb) {
    const auto& parms = _data->parms;
    switch (parms.bb_explore_strategy.value()) {
        case min_bb_explore_strategy:
            tree = std::make_unique<DFSTree>(_probType, _isIntProb);
            break;
        case bb_bfs_strategy:
            tree = std::make_unique<BFSTree>(_probType, _isIntProb);
            break;
        case bb_brfs_strategy:
            tree = std::make_unique<BrFSTree>(_probType, _isIntProb);
            break;
        case bb_cbfs_strategy:
            tree = std::make_unique<CBFSTree>(_probType, _isIntProb);
            break;
    }

    tree->set_global_ub(static_cast<double>(_data->upper_bound));
    tree->set_retain_states(false);
    tree->set_final_test_usage(parms.pruning_test.value());
    tree->set_time_limit(static_cast<int>(parms.branching_cpu_limit.value()));
    tree->set_node_limit(parms.bb_node_limit.value());
    tree->set_only_root_node(parms.bb_node_limit.value() == 1);
    auto aux = _data->depth;
    auto node = std::make_unique<BranchNodeBase>(std::move(_data), true);
    node->set_depth(aux);
    tree->process_state(std::move(node), true);
}

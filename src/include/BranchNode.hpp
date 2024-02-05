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

#ifndef BRANCHNODE_H
#define BRANCHNODE_H

#include <array>                     // for array
#include <limits>                    // for numeric_limits
#include <memory>                    // for unique_ptr
#include "branch-and-bound/state.h"  // for State

class BTree;
struct Instance;
struct NodeData;
struct PricerSolverBase;

class BranchNodeBase : public State {
   public:
    explicit BranchNodeBase(std::unique_ptr<NodeData> pd, bool isRoot = false);
    BranchNodeBase(BranchNodeBase&&) = default;
    BranchNodeBase(const BranchNodeBase&) = delete;
    auto operator=(BranchNodeBase&&) -> BranchNodeBase& = default;
    auto operator=(const BranchNodeBase&) -> BranchNodeBase& = delete;
    ~BranchNodeBase() override = default;
    [[nodiscard]] auto clone() const -> std::unique_ptr<State> override {
        return nullptr;
    }

    auto is_terminal_state() -> bool final;

    void branch(BTree* bt) override;
    void compute_bounds(BTree* bt) final;
    void assess_dominance(State* otherState) final;
    void apply_final_pruning_tests(BTree* bt) final;
    void update_data(double upper_bound) final;
    void print(const BTree* bt) const override;

    [[nodiscard]] auto get_data_ptr() const -> NodeData* { return pd.get(); }
    [[nodiscard]] auto get_instance_info() const -> const Instance&;
    [[nodiscard]] auto get_pricersolver() const -> PricerSolverBase*;

    static constexpr auto NbCandidates = 16UL;
    static constexpr auto EPS = 1e-4;

   private:
    std::unique_ptr<NodeData> pd;

    static constexpr auto ERROR = 1e-12;
    static constexpr auto IntegerTolerance = 1e-3;
    static constexpr auto MAX_NB_PATHS = 240000000000;
    static constexpr auto MAX_DEPTH = 7;
    static constexpr auto BRANCH_POINT_VALUE = 0.5;
};

class BranchNodeRelBranching : public BranchNodeBase {
   public:
    explicit BranchNodeRelBranching(std::unique_ptr<NodeData> _node,
                                    bool                      isRoot = false);
    BranchNodeRelBranching(BranchNodeRelBranching&& _branch_node) = default;
    BranchNodeRelBranching(const BranchNodeRelBranching& _branch_node) = delete;
    auto operator=(BranchNodeRelBranching&& _branch_node)
        -> BranchNodeRelBranching& = default;
    auto operator=(const BranchNodeRelBranching& _branch_node)
        -> BranchNodeRelBranching& = delete;
    ~BranchNodeRelBranching() override = default;

    void branch(BTree* bt) override;
};

struct BranchCandidate {
    double score{std::numeric_limits<double>::lowest()};
    bool   empty{true};
    std::array<std::unique_ptr<NodeData>, 2> data_child_nodes;

    BranchCandidate() = default;
    BranchCandidate(const BranchCandidate& _branch_cand) = delete;
    BranchCandidate(BranchCandidate&& _branch_cand) = default;
    auto operator=(const BranchCandidate& _branch_cand)
        -> BranchCandidate& = delete;
    auto operator=(BranchCandidate&& _branch_cand)
        -> BranchCandidate& = default;
    ~BranchCandidate() = default;

    BranchCandidate(double                                     _score,
                    std::array<std::unique_ptr<NodeData>, 2>&& child_nodes);
    BranchCandidate(double _score);
};
#endif  // BRANCHNODE_H
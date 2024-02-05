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

#ifndef NODEDATA_H
#define NODEDATA_H

#include <OsiSolverInterface.hpp>  //for OsiSolverInterface
#include <array>                   // for array
#include <cstddef>                 // for size_t
#include <functional>              // for function
#include <memory>                  // for unique_ptr, shared_ptr
#include <span>                    // for span
#include <string>                  // for string
#include <vector>                  // for vector
#include "or-utils/lp.h"           // for lp_interface

class PricingStabilizationBase;  // lines 14-14
class Problem;                   // lines 15-15
struct Column;                   // lines 19-19
struct Instance;                 // lines 17-17
struct Parms;                    // lines 18-18
struct PricerSolverBase;         // lines 16-16
struct Sol;
struct Statistics;

struct NodeData {
    enum NodeDataStatus {
        initialized = 0,
        LP_bound_estimated = 1,
        LP_bound_computed = 2,
        submitted_for_branching = 3,
        infeasible = 4,
        finished = 5,
    };

    size_t depth;

    NodeDataStatus status;

    // The instance information
    const Parms&    parms;
    const Instance& instance;
    Statistics&     stat;
    Sol&            opt_sol;

    size_t nb_jobs;
    size_t nb_machines;

    // The column generation lp information
    std::unique_ptr<OsiSolverInterface> osi_rmp;
    std::vector<int>                    row_status;

    std::span<const double> lambda{};
    std::span<const double> pi{};
    std::span<const double> rhs{};
    std::vector<double>     lhs_coeff{};
    std::vector<int>        id_row{};
    std::vector<double>     coeff_row{};

    // cut generation information
    size_t max_nb_cuts;
    int    id_convex_constraint;
    int    id_assignment_constraint;
    int    id_valid_cuts;

    size_t id_art_var_convex;
    int    id_art_var_assignment;
    int    id_art_var_cuts;
    int    id_next_var_cuts;
    int    id_pseudo_schedules;

    // PricerSolver
    std::unique_ptr<PricerSolverBase> solver{};

    // Columns
    int                                  zero_count;
    std::vector<int>                     column_status{};
    std::vector<std::shared_ptr<Column>> localColPool{};

    int lower_bound;
    int upper_bound;

    double LP_lower_bound;
    double LP_lower_bound_dual;
    double LP_lower_bound_BB;
    double LP_lower_min;

    int    nb_non_improvements;
    size_t iterations;

    /** Wentges smoothing technique */
    std::unique_ptr<PricingStabilizationBase> solver_stab;

    // maxiterations and retireage
    int retirementage;

    /** Branching strategies */
    size_t branch_job;
    int    completiontime;
    bool   less;

    /**
     * ptr to the parent node
     */
    explicit NodeData(Problem* problem);
    NodeData(const NodeData& src);
    NodeData(NodeData&&) = delete;
    ~NodeData();

    auto operator=(const NodeData&) -> NodeData& = delete;
    auto operator=(NodeData&&) -> NodeData& = delete;

    void prune_duplicated_sets();
    void add_solution_to_colpool(const Sol& sol);

    auto build_rmp() -> int;
    /** lowerbound.cpp */
    auto delete_unused_rows() -> int;
    auto delete_old_columns() -> int;
    auto delete_infeasible_columns() -> int;
    auto compute_objective() -> int;
    auto solve_relaxation() -> int;
    auto compute_lower_bound() -> int;
    auto estimate_lower_bound(size_t _iter) -> int;
    auto refinement() -> bool;
    void make_pi_feasible_farkas_pricing();

    /** PricerSolverWrappers.cpp */
    void build_solve_mip() const;
    void construct_lp_sol_from_rmp();
    void generate_cuts() const;
    auto delete_unused_rows_range(int first, int last) -> int;
    void call_update_rows_coeff() const;

    /** small getters */
    [[nodiscard]] auto get_score_value() const -> double;

    /** StabilizationWrappers.cpp */
    auto solve_pricing() -> int;
    void solve_farkas_dbl();

    [[nodiscard]] auto clone() const -> std::unique_ptr<NodeData>;
    [[nodiscard]] auto clone(size_t _j, int _t, bool _left) const
        -> std::unique_ptr<NodeData>;
    [[nodiscard]] auto create_child_nodes(size_t _j, int _t) const
        -> std::array<std::unique_ptr<NodeData>, 2>;

    [[nodiscard]] auto create_child_nodes(size_t _j, long _t) const
        -> std::array<std::unique_ptr<NodeData>, 2>;

   private:
    auto add_lhs_column_to_rmp(double cost) -> int;
    auto add_lhs_column_to_rmp(double cost, const std::vector<double>& _lhs)
        -> int;

    /**
     * @brief Model construction methods
     *
     */
    void create_assignment_constraints();
    void create_convex_contraint();
    void create_artificial_cols();
    void add_cols_local_pool();

    /** lowerbound.cpp */
    auto grow_ages() -> int;
    void print_ages();

    static constexpr auto CLEANUP_ITERATION = 30;
    static constexpr auto EPS = 1e-6;
    static constexpr auto EPS_BOUND = 1e-9;
    static constexpr auto min_nb_del_row_ratio = 0.9;
    static constexpr auto NB_CG_ITERATIONS = 1000000UL;
    static constexpr auto NB_CUTS = 2000;
    static constexpr auto NB_ESTIMATE_IT = 20;
    static constexpr auto NB_PATHS = 1000000000;
    static constexpr auto NB_NON_IMPROVEMENTS = 5;
};

#endif  // NODEDATA_H
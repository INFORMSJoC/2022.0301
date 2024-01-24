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

#ifndef PRICER_SOLVER_BASE_HPP
#define PRICER_SOLVER_BASE_HPP
#include <gurobi_c++.h>                      // for GRBModel
#include <boost/multiprecision/cpp_int.hpp>  // for cpp_int
#include <cstddef>                           // for size_t
#include <functional>                        // for reference_wrapper, ref
#include <memory>                            // for shared_ptr, unique_ptr
#include <string>                            // for string
#include <utility>                           // for ref
#include <vector>                            // for vector
#include "Instance.h"                        // for Instance
#include "MIP_defs.hpp"                      // for MIP_Attr
#include "ModelInterface.hpp"                // for BddCoeff, ReformulationModel
#include "PricingSolution.hpp"               // for PricingSolution
#include "Solution.hpp"                      // for Sol

/** Forward declarations **/
struct Job;
struct NodeData;  // lines 15-15
struct Column;    // lines 16-16

struct PricerSolverBase {
    using cpp_int = boost::multiprecision::cpp_int;

   public:
    const std::vector<std::shared_ptr<Job>>& jobs;

    size_t convex_constr_id;
    size_t convex_rhs;

    std::string problem_name{};

    std::shared_ptr<GRBEnv> env;
    GRBModel                model;

    ReformulationModel reformulation_model;

    bool is_integer_solution;

    double constLB;
    double UB;

    std::vector<std::vector<BddCoeff>> x_bar;

    static const std::shared_ptr<GRBEnv> genv;

    static constexpr double EPS_SOLVER = 1e-4;
    static constexpr double RC_FIXING = 1e-4;
    static constexpr int    ALIGN = 60;
    static constexpr int    ALIGN_HALF = 60;

    /**
     * @brief Construct a new Pricer Solver Base object
     *
     * @param instance
     */
    explicit PricerSolverBase(const Instance& instance);
    PricerSolverBase(const PricerSolverBase& other);
    PricerSolverBase(PricerSolverBase&& other) noexcept;
    virtual ~PricerSolverBase();

    auto operator=(const PricerSolverBase& other) -> PricerSolverBase&;
    auto operator=(PricerSolverBase&& other) noexcept -> PricerSolverBase&;

    /**
     * @brief clone a PricerSolverBase
     *
     * @return std::unique_ptr<PricerSolverBase>
     */
    [[nodiscard]] virtual auto clone() const
        -> std::unique_ptr<PricerSolverBase> = 0;

    /**
     * @brief pricing algorithms
     *
     * @param _pi
     * @return PricingSolution
     */
    virtual auto pricing_algorithm(std::span<const double>& pi)
        -> PricingSolution = 0;
    virtual auto farkas_pricing(std::span<const double>& pi)
        -> PricingSolution = 0;

    /**
     * @brief Functions for the evaluation of reduced cost fixing
     *
     * @param pi
     * @return bool true if the graph is reduced
     */
    virtual auto evaluate_nodes(std::span<const double>& pi) -> bool = 0;

    /**
     * @brief Refinement of the graph structure
     *
     * @param paths that should be removed from the graph
     * @return true paths removed
     * @return false paths not removed
     */
    virtual auto refinement_structure(
        [[maybe_unused]] const std::vector<std::shared_ptr<Column>>& paths)
        -> bool {
        return false;
    };

    /**
     * @brief Enumeration of the columns of the graph structure
     *
     */
    virtual void enumerate_columns(){};
    virtual void enumerate_columns(
        [[maybe_unused]] std::span<const double>& _pi){};

    /**
     * @brief Build mip over the original variables
     *
     */
    virtual void build_mip() = 0;
    auto         evaluate_mip_model() -> bool;

    /**
     * @brief Construction of the LP solutions over the original variables
     *
     * @param lambda values lambda corresponding columns
     * @param columns characteristics of the corresponding columns
     */

    virtual void construct_lp_sol_from_rmp(
        const std::span<const double>&              lambda,
        const std::vector<std::shared_ptr<Column>>& columns) = 0;

    /**
     * @brief Computation of sub-optimal duals
     *
     * @param lambda
     * @param columns
     * @return bool
     */
    auto compute_sub_optimal_duals(
        const std::span<const double>&              lambda,
        const std::vector<std::shared_ptr<Column>>& columns) -> bool;

    virtual void project_sol_on_original_variables(const Sol& _sol) {
        _sol.print_solution();
    };

    /**
     * @brief Calculation of constLB over pi
     *
     * @param pi
     */
    void calculate_constLB(std::span<const double>& pi);

    /**
     * @brief Branching operator on the original formulation
     *
     * @param _job
     * @param _time
     * @param _left
     */
    virtual void split_job_time([[maybe_unused]] size_t _job,
                                [[maybe_unused]] int    _time,
                                [[maybe_unused]] bool   _left) {}

    /**
     * @brief Compute the reduced cost function of a column
     *
     * @param sol characteristics of the column
     * @param pi the dual prices
     * @param lhs left hand side coefficients of the master program
     * @return double the reduced cost of sol
     */
    virtual auto compute_reduced_cost(const PricingSolution&   sol,
                                      std::span<const double>& pi,
                                      double*                  lhs) -> double;

    /**
     * @brief Compute the left hand side of solution sol in the master program
     *
     * @param sol chareteristics of the column
     * @param lhs coefficients of the left hand side
     */
    virtual void compute_lhs(const PricingSolution& sol, double* lhs);
    virtual void compute_lhs(const Column& sol, double* lhs);

    /**
     * @brief Compute only the reduced cost of a solution
     *
     * @param sol charecteristics of the column
     * @param pi the dual prices
     * @return double reduced cost associated with column sol and dual prices pi
     */
    virtual auto compute_reduced_cost_simple(const PricingSolution&   sol,
                                             std::span<const double>& pi)
        -> double;

    /**
     * @brief Compute Lagrangian bound of pi when a optimal solution is known
     * for pi
     *
     * @param sol optimal solution for pi
     * @param pi  dual prices
     * @return double Lagrangian bound
     */
    virtual auto compute_lagrange(const PricingSolution&     sol,
                                  const std::vector<double>& pi) -> double;
    virtual auto compute_lagrange(const PricingSolution&         sol,
                                  const std::span<const double>& pi) -> double;

    virtual auto compute_subgradient(const PricingSolution& sol,
                                     double* subgradient) -> double;

    inline void set_is_integer_solution(bool _is_solution) {
        is_integer_solution = _is_solution;
    }

    /**
     * Constraint on the solver
     */

    virtual void insert_constraints_lp(NodeData* pd) = 0;
    virtual auto add_constraints() -> int;
    virtual void remove_constraints(int first, int nb_del);
    virtual void update_rows_coeff(size_t first);
    virtual void update_coeff_constraints() = 0;
    virtual void update_constraints() = 0;

    /**
     * Some getters
     */
    virtual auto       print_num_paths() -> cpp_int = 0;
    [[nodiscard]] auto get_UB() const -> double;
    void               update_UB(double _ub);

    virtual auto get_nb_vertices() -> size_t = 0;
    virtual auto get_nb_edges() -> size_t = 0;
    virtual auto structure_feasible() -> bool { return true; }
    virtual auto check_column([[maybe_unused]] Column const* set) -> bool {
        return true;
    };

    virtual auto get_size_data() -> size_t { return 0UL; };
    virtual auto calculate_job_time() -> std::vector<std::vector<BddCoeff>>& {
        return x_bar;
    };

    [[nodiscard]] inline auto get_is_integer_solution() const -> bool {
        return is_integer_solution;
    }

    /**
     * Some printing functions
     */
    [[maybe_unused]] virtual auto get_int_attr_model(enum MIP_Attr c) -> int;
    virtual auto                  get_dbl_attr_model(enum MIP_Attr c) -> double;
};

#endif  // PRICER_SOLVER_BASE_HPP

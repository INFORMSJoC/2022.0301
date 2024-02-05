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

#ifndef PRICER_SOLVER_SIMPLE_DP_HPP
#define PRICER_SOLVER_SIMPLE_DP_HPP
#include <cstddef>               // for size_t
#include <memory>                // for make_unique, unique_ptr, shared_ptr
#include <vector>                // for vector
#include "Instance.h"            // for Instance
#include "PricerSolverBase.hpp"  // for PricerSolverBase
#include "PricingSolution.hpp"   // for PricingSolution
#include "gurobi_c++.h"          // for GRBVar
struct Job;
struct NodeData;
struct Column;
class PricerSolverSimpleDp : public PricerSolverBase {
   private:
    size_t                         Hmax;
    size_t                         size_graph;
    std::vector<Job*>              A;
    std::vector<double>            F;
    std::vector<double>            backward_F;
    std::vector<std::vector<Job*>> backward_graph;
    std::vector<std::vector<Job*>> forward_graph;
    std::vector<GRBVar>            TI_x;
    std::vector<bool>              take;
    std::vector<double>            lp_x;
    std::vector<double>            solution_x;

   public:
    PricerSolverSimpleDp(const Instance& instance);
    PricerSolverSimpleDp(PricerSolverSimpleDp&&) = default;
    ~PricerSolverSimpleDp() override = default;
    auto operator=(PricerSolverSimpleDp&&) -> PricerSolverSimpleDp& = default;
    auto operator=(const PricerSolverSimpleDp&)
        -> PricerSolverSimpleDp& = default;

    PricerSolverSimpleDp(const PricerSolverSimpleDp& src)
        : PricerSolverBase(src),
          Hmax(src.Hmax),
          size_graph(src.size_graph),
          A(Hmax + 1),
          F(Hmax + 1),
          backward_F(Hmax + 1),
          TI_x(convex_constr_id * (Hmax + 1), GRBVar()),
          take(convex_constr_id * (Hmax + 1)),
          lp_x(convex_constr_id * (Hmax + 1), 0.0),
          solution_x(convex_constr_id * (Hmax + 1)) {
        init_table();
    };

    [[nodiscard]] auto clone() const
        -> std::unique_ptr<PricerSolverBase> override {
        return std::make_unique<PricerSolverSimpleDp>(*this);
    };

    void init_table();

    auto evaluate_nodes([[maybe_unused]] std::span<const double>& pi)
        -> bool override;
    void build_mip() override;

    void construct_lp_sol_from_rmp(
        const std::span<const double>&              lambda,
        const std::vector<std::shared_ptr<Column>>& columns) override;
    auto get_nb_edges() -> size_t override;
    auto get_nb_vertices() -> size_t override;
    auto print_num_paths() -> cpp_int override;

    auto check_column(Column const* set) -> bool override;

    auto pricing_algorithm(std::span<const double>& _pi)
        -> PricingSolution override;
    auto farkas_pricing(std::span<const double>& _pi)
        -> PricingSolution override;
    void forward_evaluator(double* _pi);
    void forward_evaluator(std::span<const double>& _pi);
    void backward_evaluator(double* _pi);
    void backward_evaluator(std::span<const double>& _pi);

    void update_constraints() override {}

    void insert_constraints_lp([[maybe_unused]] NodeData* pd) override {}

    void update_coeff_constraints() override {}
};

#endif  // PRICER_SOLVER_SIMPLE_DP_HPP

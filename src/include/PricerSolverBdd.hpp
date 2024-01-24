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

#ifndef PRICER_SOLVER_BDD_HPP
#define PRICER_SOLVER_BDD_HPP
#include <OsiGrbSolverInterface.hpp>      // for OsiGrbSolverInterface
#include <cstddef>                        // for size_t
#include <memory>                         // for unique_ptr, shared_ptr
#include <span>                           // for span
#include <utility>                        // for pair
#include <vector>                         // for vector
#include "Instance.h"                     // for Instance
#include "MipGraph.hpp"                   // for MipGraph
#include "ModelInterface.hpp"             // for OriginalModel
#include "ModernDD/NodeBddStructure.hpp"  // for DdStructure
#include "ModernDD/NodeBddTable.hpp"      // for NodeTableEntity, TableHandler
#include "NodeBdd.hpp"                    // for NodeBdd
#include "PricerSolverBase.hpp"           // for PricerSolverBase
#include "PricingSolution.hpp"            // for PricingSolution
struct Interval;
struct Job;
struct NodeData;
struct Column;
struct Sol;

class PricerSolverBdd : public PricerSolverBase {
    DdStructure<NodeBdd> decision_diagram;

    size_t size_graph{};
    size_t nb_edges{};
    size_t nb_vertices{};
    size_t nb_removed_edges{};
    size_t nb_removed_nodes{};

    std::vector<std::pair<Job*, Interval*>> ordered_jobs_new;

    MipGraph mip_graph;

    OriginalModel<> original_model;

    // std::unordered_map<int, std::vector<std::weak_ptr<NodeId>>> t_in;
    // std::unordered_map<int, std::vector<std::weak_ptr<NodeId>>> t_out;

    size_t H_min;
    int    H_max;

   public:
    explicit PricerSolverBdd(const Instance& instance);

    PricerSolverBdd(const PricerSolverBdd& src);
    PricerSolverBdd(PricerSolverBdd&&) = default;
    ~PricerSolverBdd() override;

    auto operator=(PricerSolverBdd&&) -> PricerSolverBdd& = default;
    auto operator=(const PricerSolverBdd&) -> PricerSolverBdd& = delete;

    [[maybe_unused]] void check_infeasible_arcs();
    void                  topdown_filtering();
    void                  bottom_up_filtering();
    void                  equivalent_paths_filtering();
    [[maybe_unused]] void print_representation_file();
    void                  cleanup_arcs();

    auto refinement_structure(const std::vector<std::shared_ptr<Column>>& paths)
        -> bool override;
    void enumerate_columns() override;
    void enumerate_columns(std::span<const double>& _pi) override;
    auto evaluate_nodes(std::span<const double>& _pi) -> bool override;

    auto               check_column(Column const* col) -> bool override;
    [[nodiscard]] auto clone() const
        -> std::unique_ptr<PricerSolverBase> override = 0;
    virtual auto evaluate_rc_arc(NodeBdd& n) -> double = 0;
    auto         evaluate_rc_low_arc(NodeBdd& n) -> double;
    virtual void compute_labels(std::span<const double>& _pi) = 0;

    void                  remove_layers();
    void                  remove_edges();
    [[maybe_unused]] void remove_layers_init();
    void                  construct_mipgraph();
    void                  init_coeff_constraints();
    void                  init_table();

    void build_mip() override;
    void construct_lp_sol_from_rmp(
        const std::span<const double>&              lambda,
        const std::vector<std::shared_ptr<Column>>& columns) override;

    void project_sol_on_original_variables(const Sol& _sol) override;
    auto calculate_job_time() -> std::vector<std::vector<BddCoeff>>& override;
    void split_job_time(size_t _job, int _time, bool _left) override;
    auto print_num_paths() -> cpp_int override;
    void remove_constraints(int first, int nb_del) override;
    void update_rows_coeff(size_t first) override;
    void insert_constraints_lp(NodeData* pd) override;

    auto add_constraints() -> int override;

    auto get_nb_edges() -> size_t override;
    auto get_nb_vertices() -> size_t override;
    auto structure_feasible() -> bool override;
    auto get_size_data() -> size_t override {
        return decision_diagram.getDiagram()->totalSize();
    };

    inline auto get_decision_diagram() -> DdStructure<NodeBdd>& {
        return decision_diagram;
    }

    [[nodiscard]] inline auto get_nb_removed_edges() const {
        return nb_removed_edges;
    }

    inline void add_nb_removed_edges() { nb_removed_edges++; }

   private:

    auto compute_reduced_cost(const PricingSolution&   sol,
                              std::span<const double>& pi,
                              double*                  lhs) -> double override;

    void compute_lhs(const PricingSolution& sol, double* lhs) override;
    void compute_lhs(const Column& sol, double* lhs) override;

    auto compute_lagrange(const PricingSolution&     sol,
                          const std::vector<double>& pi) -> double override;
    auto compute_lagrange(const PricingSolution&         sol,
                          const std::span<const double>& pi) -> double override;

    auto compute_subgradient(const PricingSolution& sol, double* sub_gradient)
        -> double override;

    void update_constraints() override {}

    void update_coeff_constraints() override;
    auto build_model() -> std::unique_ptr<OsiSolverInterface>;
};

#endif  // PRICER_SOLVER_BDD_HPP

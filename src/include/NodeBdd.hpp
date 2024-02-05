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

#ifndef NODE_DURATION_HPP
#define NODE_DURATION_HPP

#include <gurobi_c++.h>                             // for GRBVar
#include <array>                                    // for array
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <boost/multiprecision/cpp_int.hpp>         // for cpp_int
#include <cstddef>                                  // for size_t
#include <memory>                                   // for allocator, shared...
#include <ostream>                                  // for operator<<, ostream
#include <range/v3/algorithm/any_of.hpp>            // for any_of
#include <range/v3/algorithm/max.hpp>               // for max, max_fn
#include <range/v3/functional/arithmetic.hpp>       // for plus
#include <range/v3/functional/comparisons.hpp>      // for less
#include <range/v3/functional/identity.hpp>         // for identity
#include <range/v3/numeric/accumulate.hpp>          // for accumulate, accum...
#include <vector>                                   // for vector
#include "BddCoeff.hpp"                             // for BddCoeff
#include "Job.h"                                    // for Job
#include "Label.hpp"                                // for Label
#include "ModernDD/NodeBase.hpp"                    // for NodeBase
#include "ModernDD/NodeId.hpp"                      // for NodeId

class NodeBdd : public NodeBase {
   public:
    using dbl_array = std::array<double, 2>;
    using int_array = std::array<int, 2>;
    using size_t_array [[maybe_unused]] = std::array<size_t, 2>;
    using bool_array = std::array<bool, 2>;
    using grb_array = std::array<GRBVar, 2>;
    using big_int = boost::multiprecision::cpp_int;
    using weak_ptr_nodeid = std::weak_ptr<NodeId>;
    using weak_ptr_bddcoeff = std::weak_ptr<BddCoeff>;

   private:
    /** Job chatacteristics */
    Job*   job{nullptr};
    size_t weight{};
    size_t key{};
    size_t key_model{};

    bool   visited{false};
    NodeId ptr_node_id{};

    /** Calculation jobs */
    boost::dynamic_bitset<> all{};

    /** original model data */
    dbl_array lp_x{0.0, 0.0};
    dbl_array cost{0.0, 0.0};
    dbl_array best_sol_x{0.0, 0.0};
    dbl_array reduced_cost{0.0, 0.0};

    std::array<std::vector<weak_ptr_bddcoeff>, 2> coeff_list{};

    /** Evaluation of backward distance */
    int_array backward_distance{};

    /** Evaluation nb_paths */
    big_int nb_paths{};

    /** Reducing size of decision diagram */
    bool_array calc{true, true};

    /** Searching equivalent paths */
    int_array in_degree{};

    /** Zero-Half cuts formulation variables */
    bool                               lp_visited{false};
    std::array<std::vector<size_t>, 2> in_edges{};
    std::array<size_t, 2>              key_edges{};
    dbl_array                          coeff_cut{0.0, 0.0};
    grb_array                          y{};
    grb_array                          r{};
    GRBVar                             sigma{};

   public:
    std::array<Label<NodeBdd>, 2> forward_label{};
    std::array<Label<NodeBdd>, 2> backward_label{};

    /**
     * Constructor
     */
    NodeBdd() = default;
    NodeBdd(size_t i, size_t j) : NodeBase(i, j) {}

    NodeBdd(const NodeBdd& src) = default;
    NodeBdd(NodeBdd&& src) noexcept = default;
    ~NodeBdd() = default;

    auto operator=(const NodeBdd& src) -> NodeBdd& = default;
    auto operator=(NodeBdd&& src) noexcept -> NodeBdd& = default;

    void               set_weight(int _weight);
    [[nodiscard]] auto get_weight() const -> size_t;

    [[nodiscard]] auto get_job() const -> Job*;
    [[nodiscard]] auto get_nb_job() const -> size_t;
    void               set_job(Job* _job);

    void add_coeff_list(const std::shared_ptr<BddCoeff>& ptr, bool high);
    void add_coeff_list_clear();
    auto get_coeff_list() -> std::array<std::vector<weak_ptr_bddcoeff>, 2>&;

    auto operator!=(NodeBdd const& o) const -> bool;

    friend auto operator<<(std::ostream& os, NodeBdd const& o) -> std::ostream&;

    void init_node(size_t                _weight,
                   [[maybe_unused]] bool _root_node = false,
                   bool                  _terminal_node = false);

    void set_node_id_label(const NodeId& _id);

    void set_job_label(Job* _job);

    auto operator<=>(const NodeBdd& rhs);

    /** Functions concerning reduced_cost */
    auto get_reduced_cost() -> dbl_array&;
    void reset_reduced_costs();
    void reset_reduced_costs_farkas();
    void adjust_reduced_costs(double _x, bool high);

    /** Functions concerning ptr_node_id */
    void set_ptr_node_id(size_t i, size_t j);
    void set_ptr_node_id(const NodeId& _node_id);
    auto get_ptr_node_id() -> NodeId&;
    void reset_ptr_node_id();

    /** Functions for manipulation of nb_paths */
    auto get_nb_paths() -> big_int&;
    void reset_nb_paths();
    void update_nb_paths(const big_int& x = 1);

    /** Functions for manipulation of visited */
    void               reset_visited();
    void               update_visited();
    [[nodiscard]] auto get_visited() const -> bool;

    /** Functions for manipulation of key */
    void set_key(const size_t& _key);
    auto get_key() -> size_t&;

    void set_key_model(const size_t& _key);
    auto get_key_model() -> size_t&;
    /** Functions for manipulation of lp_visited */
    void               update_lp_visited(bool _update);
    [[nodiscard]] auto get_lp_visited() const -> bool;

    /** Functions for manipulation of lp_x */
    [[nodiscard]] auto get_lp_x() -> dbl_array&;
    auto               get_lp_x(bool _high) -> double&;

    void update_lp_x(double _x, bool _high);

    void reset_lp_x();

    /** Functions for manipulation of cost */
    void set_cost(double _cost);

    /** Functions for manipulation of best_sol_x */
    void update_best_sol_x(bool _high);
    auto get_best_sol_x(bool _high) -> double;

    /** Functions for manipulation of backward distance */
    void reset_backward_distance(int _x = 0);
    void update_backward_distance(const std::array<int, 2>& _backward_distance,
                                  bool                      _high);

    auto get_backward_distance() -> int_array&;

    /** Functions for manipulation of all */
    void reset_all(size_t _nb_elements);

    auto get_all() -> boost::dynamic_bitset<>&;

    void intersect_all(const boost::dynamic_bitset<>& _set);
    void union_all(const boost::dynamic_bitset<>& _set);
    void add_element(size_t _element);
    auto is_element(size_t _element) -> bool;
    auto all_is_empty() -> bool;
    auto get_first_all() -> size_t;
    auto get_next_all(size_t _index) -> size_t;

    /** Functions for manipulation of calc */
    void reset_calc();
    auto any_of_calc() -> bool;

    auto get_calc(bool _high) -> bool;

    void update_calc(bool _high, bool _update = false);

    /** Functions for manipulation of in_degree */
    auto get_in_degree(bool _high) -> int;

    void update_in_degree(bool _high);

    auto alternative_paths() -> bool;

    void reset_in_degree();

    /** Functions for manipulation of Zero-Half cuts formulation cut problem */
    void reset_in_edges();

    auto get_in_edges(bool _high) -> std::vector<size_t>&;
    void add_in_edge(bool _high, size_t _key_parent);

    void set_key_edge(bool _high, size_t _key);

    auto get_key_edge(bool _high) -> size_t&;

    void reset_coeff_cut();
    void reset_coeff_cut(bool _high);

    void update_coeff_cut(bool _high);

    void reduce_coeff_cut(bool _high);

    auto get_coeff_cut(bool _high) -> double;

    auto get_y() -> grb_array&;
    auto get_r() -> grb_array&;
    void set_sigma(GRBVar&& _sigma);
    auto get_sigma() -> GRBVar&;
};

#endif  // NODE_DURATION_HPP

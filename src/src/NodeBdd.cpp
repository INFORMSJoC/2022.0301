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

#include "NodeBdd.hpp"
#include <range/v3/view/zip.hpp>

void NodeBdd::set_weight(int _weight) {
    weight = _weight;
}
[[nodiscard]] auto NodeBdd::get_weight() const -> size_t {
    return weight;
}

[[nodiscard]] auto NodeBdd::get_job() const -> Job* {
    return job;
}
[[nodiscard]] auto NodeBdd::get_nb_job() const -> size_t {
    return job->job;
}
void NodeBdd::set_job(Job* _job) {
    job = _job;
}

void NodeBdd::add_coeff_list(const std::shared_ptr<BddCoeff>& ptr, bool high) {
    if (high) {
        coeff_list[1].push_back(ptr);
    } else {
        coeff_list[0].push_back(ptr);
    }
}

void NodeBdd::add_coeff_list_clear() {
    for (auto& it : coeff_list) {
        it.clear();
    }
}

auto NodeBdd::get_coeff_list()
    -> std::array<std::vector<NodeBdd::weak_ptr_bddcoeff>, 2>& {
    return coeff_list;
}

auto NodeBdd::operator!=(NodeBdd const& o) const -> bool {
    return !(*this == (o));
}

auto operator<<(std::ostream& os, NodeBdd const& o) -> std::ostream& {
    os << "(" << o[0];

    for (int i = 1; i < 2; ++i) {
        os << "," << o[i];
    }

    return os << ")";
}

void NodeBdd::init_node(size_t                _weight,
                        [[maybe_unused]] bool _root_node,
                        bool                  _terminal_node) {
    if (!_terminal_node) {
        weight = _weight;
    } else {
        set_job(nullptr);
        weight = -1;
    }
}

void NodeBdd::set_node_id_label(const NodeId& _id) {
    for (auto&& [b, f] : ranges::views::zip(backward_label, forward_label)) {
        f.set_node_id(_id);
        b.set_node_id(_id);
    }
}

void NodeBdd::set_job_label(Job* _job) {
    for (auto&& [f, b] : ranges::views::zip(forward_label, backward_label)) {
        f.set_job(_job);
        b.set_job(_job);
    }
}

auto NodeBdd::operator<=>(const NodeBdd& rhs) {
    return this->forward_label[0].get_f() <=> rhs.forward_label[0].get_f();
}

/** Functions concerning reduced_cost */
auto NodeBdd::get_reduced_cost() -> NodeBdd::dbl_array& {
    return reduced_cost;
}

void NodeBdd::reset_reduced_costs() {
    reduced_cost = cost;
}

void NodeBdd::reset_reduced_costs_farkas() {
    reduced_cost = {0.0, 0.0};
}

void NodeBdd::adjust_reduced_costs(double _x, bool high) {
    if (high) {
        reduced_cost[1] -= _x;
    } else {
        reduced_cost[0] -= _x;
    }
}

/** Functions concerning ptr_node_id */
void NodeBdd::set_ptr_node_id(size_t i, size_t j) {
    ptr_node_id = NodeId(i, j);
}

void NodeBdd::set_ptr_node_id(const NodeId& _node_id) {
    ptr_node_id = _node_id;
}

auto NodeBdd::get_ptr_node_id() -> NodeId& {
    return ptr_node_id;
}
void NodeBdd::reset_ptr_node_id() {
    ptr_node_id = 0;
}

/** Functions for manipulation of nb_paths */
auto NodeBdd::get_nb_paths() -> NodeBdd::big_int& {
    return nb_paths;
}
void NodeBdd::reset_nb_paths() {
    nb_paths = 0UL;
}
void NodeBdd::update_nb_paths(const big_int& x) {
    nb_paths += x;
}

/** Functions for manipulation of visited */
void NodeBdd::reset_visited() {
    visited = false;
}
void NodeBdd::update_visited() {
    visited = true;
}
auto NodeBdd::get_visited() const -> bool {
    return visited;
}

/** Functions for manipulation of key */
void NodeBdd::set_key(const size_t& _key) {
    key = _key;
}
auto NodeBdd::get_key() -> size_t& {
    return key;
}

void NodeBdd::set_key_model(const size_t& _key) {
    key_model = _key;
}
auto NodeBdd::get_key_model() -> size_t& {
    return key_model;
}

/** Functions for manipulation of lp_visited */
void NodeBdd::update_lp_visited(bool _update) {
    lp_visited = _update;
}
auto NodeBdd::get_lp_visited() const -> bool {
    return lp_visited;
}

/** Functions for manipulation of lp_x */
[[nodiscard]] auto NodeBdd::get_lp_x() -> NodeBdd::dbl_array& {
    return lp_x;
}
auto NodeBdd::get_lp_x(bool _high) -> double& {
    return _high ? lp_x[1] : lp_x[0];
}

void NodeBdd::update_lp_x(double _x, bool _high) {
    if (_high) {
        lp_x[1] += _x;
    } else {
        lp_x[0] += _x;
    }
}

void NodeBdd::reset_lp_x() {
    lp_x = {0.0, 0.0};
    lp_visited = false;
}

/** Functions for manipulation of cost */
void NodeBdd::set_cost(double _cost) {
    cost = {0.0, _cost};
}

/** Functions for manipulation of best_sol_x */
void NodeBdd::update_best_sol_x(bool _high) {
    if (_high) {
        best_sol_x[1] += 1.0;
    } else {
        best_sol_x[0] += 1.0;
    }
}
auto NodeBdd::get_best_sol_x(bool _high) -> double {
    if (_high) {
        return best_sol_x[1];
    }
    return best_sol_x[0];
}

/** Functions for manipulation of backward distance */
void NodeBdd::reset_backward_distance(int _x) {
    backward_distance = {_x, _x};
}
void NodeBdd::update_backward_distance(
    const std::array<int, 2>& _backward_distance,
    bool                      _high) {
    if (_high) {
        backward_distance[1] =
            ranges::max(_backward_distance) + get_job()->processing_time;
    } else {
        backward_distance[0] = ranges::max(_backward_distance);
    }
}

auto NodeBdd::get_backward_distance() -> NodeBdd::int_array& {
    return backward_distance;
}

/** Functions for manipulation of all */
void NodeBdd::reset_all(size_t _nb_elements) {
    all = boost::dynamic_bitset<>{_nb_elements, 0};
}

auto NodeBdd::get_all() -> boost::dynamic_bitset<>& {
    return all;
}

void NodeBdd::intersect_all(const boost::dynamic_bitset<>& _set) {
    all &= _set;
}
void NodeBdd::union_all(const boost::dynamic_bitset<>& _set) {
    all |= _set;
}
void NodeBdd::add_element(size_t _element) {
    all[_element] = true;
}
auto NodeBdd::is_element(size_t _element) -> bool {
    return all[_element];
}
auto NodeBdd::all_is_empty() -> bool {
    return all.empty();
}
auto NodeBdd::get_first_all() -> size_t {
    return all.find_first();
}
auto NodeBdd::get_next_all(size_t _index) -> size_t {
    return all.find_next(_index);
}

/** Functions for manipulation of calc */
void NodeBdd::reset_calc() {
    calc = {true, true};
}
auto NodeBdd::any_of_calc() -> bool {
    return ranges::any_of(calc, [](auto& it) { return it; });
}

auto NodeBdd::get_calc(bool _high) -> bool {
    return _high ? calc[1] : calc[0];
}

void NodeBdd::update_calc(bool _high, bool _update) {
    if (_high) {
        calc[1] = _update;
    } else {
        calc[0] = _update;
    }
}

/** Functions for manipulation of in_degree */
auto NodeBdd::get_in_degree(bool _high) -> int {
    return _high ? in_degree[1] : in_degree[0];
}

void NodeBdd::update_in_degree(bool _high) {
    if (_high) {
        ++in_degree[1];
    } else {
        ++in_degree[0];
    }
}

auto NodeBdd::alternative_paths() -> bool {
    return (ranges::accumulate(in_degree, 0, ranges::plus{}) >= 2);
}

void NodeBdd::reset_in_degree() {
    in_degree = {0, 0};
}

/** Functions for manipulation of Zero-Half cuts formulation cut problem */
void NodeBdd::reset_in_edges() {
    for (auto& it : in_edges) {
        it.clear();
    }
}

auto NodeBdd::get_in_edges(bool _high) -> std::vector<size_t>& {
    return _high ? in_edges[1] : in_edges[0];
}
void NodeBdd::add_in_edge(bool _high, size_t _key_parent) {
    if (_high) {
        in_edges[1].emplace_back(_key_parent);
    } else {
        in_edges[0].emplace_back(_key_parent);
    }
}

void NodeBdd::set_key_edge(bool _high, size_t _key) {
    if (_high) {
        key_edges[1] = _key;
    } else {
        key_edges[0] = _key;
    }
}

auto NodeBdd::get_key_edge(bool _high) -> size_t& {
    return _high ? key_edges[1] : key_edges[0];
}

void NodeBdd::reset_coeff_cut() {
    coeff_cut = {0.0, 0.0};
}
void NodeBdd::reset_coeff_cut(bool _high) {
    if (_high) {
        coeff_cut[1] = 0.0;
    } else {
        coeff_cut[0] = 0.0;
    }
}

void NodeBdd::update_coeff_cut(bool _high) {
    if (_high) {
        coeff_cut[1] += 1.0;
    } else {
        coeff_cut[0] += 1.0;
    }
}

void NodeBdd::reduce_coeff_cut(bool _high) {
    if (_high) {
        coeff_cut[1] -= 1.0;
    } else {
        coeff_cut[0] -= 1.0;
    }
}

auto NodeBdd::get_coeff_cut(bool _high) -> double {
    return _high ? coeff_cut[1] : coeff_cut[0];
}

auto NodeBdd::get_y() -> NodeBdd::grb_array& {
    return y;
}
auto NodeBdd::get_r() -> NodeBdd::grb_array& {
    return r;
}
void NodeBdd::set_sigma(GRBVar&& _sigma) {
    sigma = _sigma;
}
auto NodeBdd::get_sigma() -> GRBVar& {
    return sigma;
}

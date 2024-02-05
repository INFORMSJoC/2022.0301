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

#ifndef LOCALSEARCH_HPP
#define LOCALSEARCH_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <list>
#include <random>
#include <vector>
#include "Solution.hpp"

struct SlopeType {
    int b1{};    /* begin of segment*/
    int b2{};    /* end of segment*/
    int c{};     /* value of the function at b1*/
    int alpha{}; /* slope of the function */

    SlopeType() = default;
    SlopeType(const SlopeType&) = default;
    SlopeType(SlopeType&&) = default;
    auto operator=(SlopeType&&) -> SlopeType& = default;
    auto operator=(const SlopeType&) -> SlopeType& = default;
    ~SlopeType() = default;

    SlopeType(int _b1, int _b2, int _c, int _alpha)
        : b1(_b1),
          b2(_b2),
          c(_c),
          alpha(_alpha) {}

    auto operator()(int t) const -> int { return c + alpha * (t - b1); };
    [[nodiscard]] auto in_interval(int _c) const -> bool {
        return (b1 <= _c) && (_c < b2);
    };
};

struct ProcessingListData {
    std::size_t pos{};
    int         p{};

    ProcessingListData() = default;
    ProcessingListData(std::size_t _pos, int _p) : pos(_pos), p(_p){};
    ProcessingListData(const ProcessingListData&) = default;
    ProcessingListData(ProcessingListData&&) = default;
    auto operator=(const ProcessingListData&) -> ProcessingListData& = default;
    auto operator=(ProcessingListData&&) -> ProcessingListData& = default;
    ~ProcessingListData() = default;
};

using MatrixProcessingList = std::vector<std::vector<ProcessingListData>>;
using VectorProcessingList = std::vector<ProcessingListData>;
using MatrixInt = std::vector<std::vector<int>>;
using VectorInt = std::vector<int>;
using SlopeList = std::vector<SlopeType>;
using SlopeListIt = std::vector<SlopeType>::iterator;

struct LocalSearchData {
    LocalSearchData(size_t _nb_jobs, size_t _nb_machines);
    void RVND(Sol& sol);

   private:
    size_t nb_jobs;
    size_t nb_machines;
    bool   updated{};
    int    iterations{};

    constexpr static int nb_moves = 15;

    std::array<std::function<void(Sol&)>, nb_moves> moves{
        [&](Sol& sol) { insertion_operator(sol, 1); },
        [&](Sol& sol) { swap_operator(sol, 0, 1); },
        [&](Sol& sol) { insertion_operator(sol, 2); },
        [&](Sol& sol) { swap_operator(sol, 0, 2); },
        [&](Sol& sol) { swap_operator(sol, 1, 1); },
        [&](Sol& sol) { insertion_operator_inter(sol, 1); },
        [&](Sol& sol) { insertion_operator_inter(sol, 2); },
        [&](Sol& sol) { swap_operator_inter(sol, 1, 1); },
        [&](Sol& sol) { swap_operator_inter(sol, 1, 2); },
        [&](Sol& sol) { swap_operator_inter(sol, 1, 3); },
        [&](Sol& sol) { swap_operator_inter(sol, 2, 2); },
        [&](Sol& sol) { swap_operator_inter(sol, 2, 3); },
        [&](Sol& sol) { swap_operator_inter(sol, 3, 3); },
        [&](Sol& sol) { swap_operator_inter(sol, 3, 4); },
        [&](Sol& sol) { swap_operator_inter(sol, 4, 4); }};

    MatrixInt                           W;
    std::vector<std::vector<SlopeList>> g;

    std::array<MatrixProcessingList, 2> processing_list;

    MatrixInt                G;
    MatrixInt                H;
    VectorInt                GG;
    MatrixInt                HH;
    std::vector<SlopeListIt> begin;
    std::vector<SlopeListIt> end;
    MatrixInt                B2_1;
    MatrixInt                B2_2;
    MatrixInt                B3_1;
    MatrixInt                B3_2;
    MatrixInt                B4_1;
    MatrixInt                B4_2;
    VectorInt                B3_1_;
    MatrixInt                B5_1;
    MatrixInt                B5_2;
    MatrixInt                B6_1;

    void calculate_W(const Sol& sol);
    void calculate_g(Sol& sol);
    void create_processing_insertion_operator(const Sol& sol, size_t l);
    void create_processing_insertion_operator_inter(const Sol& sol, size_t l);
    void create_processing_list_swap(const Sol& sol, size_t l1, size_t l2);
    void create_processing_list_swap_inter(const Sol& sol,
                                           size_t     l1,
                                           size_t     l2);
    void swap_operator(Sol& sol, size_t l1, size_t l2);
    void swap_operator_inter(Sol& sol, size_t l1, size_t l2);
    void insertion_operator(Sol& sol, size_t l);
    void insertion_operator_inter(Sol& sol, size_t l);
};

struct PerturbOperator {
    constexpr static int max_nb_perturbations = 8;
    std::mt19937         gen = std::mt19937(0);

    std::array<std::function<void(Sol&)>, max_nb_perturbations>
                                    perturbation_list;
    std::uniform_int_distribution<> dist;

    PerturbOperator()
        : perturbation_list{
              [&](Sol& sol) { sol.perturb_swap_inter(1UL, 2UL, gen); },
              [&](Sol& sol) { sol.perturb_swap_inter(1UL, 2UL, gen); },
              [&](Sol& sol) { sol.perturb_swap_inter(1UL, 3UL, gen); },
              [&](Sol& sol) { sol.perturb_swap_inter(1UL, 3UL, gen); },
              [&](Sol& sol) { sol.perturb_swap_inter(2UL, 2UL, gen); },
              [&](Sol& sol) { sol.perturb_swap_inter(2UL, 2UL, gen); },
              [&](Sol& sol) { sol.perturb_swap_inter(2UL, 3UL, gen); },
              [&](Sol& sol) { sol.perturb_swap_inter(2UL, 3UL, gen); }},
          dist(0, max_nb_perturbations - 1) {}
    void operator()(Sol& sol) {
        std::ranges::shuffle(perturbation_list.begin(), perturbation_list.end(),
                             gen);
        int nb = dist(gen);
        std::ranges::for_each(perturbation_list.begin(),
                              perturbation_list.begin() + nb,
                              [&sol](auto& operation) { operation(sol); });
    }
};

#endif  // LOCALSEARCH_HPP
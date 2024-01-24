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

#ifndef SOLUTION_NEW_H
#define SOLUTION_NEW_H

#include <cstddef>  // for size_t
#include <limits>   // for numeric_limits
#include <memory>   // for shared_ptr
#include <random>   // for mt19937
#include <vector>   // for vector

/** Forward declaration */
struct Instance;  // lines 60-60
struct Interval;  // lines 19-19
struct Job;       // lines 17-17

using VecJobPtr = std::vector<std::shared_ptr<Job>>;
using VecJobRawPtr = std::vector<Job*>;
using VecIntervalPtr = std::vector<std::shared_ptr<Interval>>;

template <typename IT>
void swap_ranges(IT start_a, IT end_a, IT start_b, IT end_b) {
    auto it = std::rotate(start_a, start_b, end_b);
    auto new_start_a = (end_a - start_a) + it;
    std::rotate(it, new_start_a, end_b);
}

struct Machine {
    VecJobRawPtr job_list{};

    int total_processing_time{};
    int total_weighted_tardiness{};

    bool updated{true};

    Machine() = default;
    Machine(Machine&&) = default;
    Machine(const Machine&) = default;
    ~Machine() = default;

    auto operator=(Machine&&) -> Machine& = default;
    auto operator=(const Machine&) -> Machine& = default;

    void add_job(Job* job);
    void reset_machine(std::vector<int>& c);
};

struct Sol {
    std::vector<Machine> machines{};

    size_t nb_jobs{};
    size_t nb_machines{};

    std::vector<int>    c{};
    std::vector<size_t> u{};

    int tw{std::numeric_limits<int>::max()};
    int off{};

    Sol() = default;
    explicit Sol(const Instance& instance);
    Sol(const Sol&) = default;
    Sol(Sol&&) = default;
    ~Sol() = default;

    auto operator=(const Sol&) -> Sol& = default;
    auto operator=(Sol&&) -> Sol& = default;

    void construct_edd(VecJobPtr& v);
    void construct_spt(const VecJobPtr& v);
    void construct_random_fisher_yates(const VecJobPtr& v);
    void construct_random_shuffle(const VecJobPtr& v);
    void canonical_order(const VecIntervalPtr& intervals);
    void perturb_swap_inter(size_t l1, size_t l2, std::mt19937& mt);

    void update_insertion_move(size_t i, size_t j, size_t k, size_t l);
    void update_swap_move(size_t i, size_t j, size_t k, size_t l1, size_t l2);
    void update_insertion_move_inter(size_t i,
                                     size_t j,
                                     size_t k,
                                     size_t l,
                                     size_t m);
    void update_swap_move_inter(size_t i,
                                size_t j,
                                size_t k1,
                                size_t k2,
                                size_t l1,
                                size_t l2);
    void print_solution() const;

   private:
    void add_job_front_machine(Job* job);
    void calculate_partition(const VecIntervalPtr& v);
};

#endif  // SOLUTION_NEW_H
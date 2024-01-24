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

#ifndef SCHEDULESET_H
#define SCHEDULESET_H

#include <ostream>              // for operator<<, ostream, basic_ostream::o...
#include <vector>               // for vector
#include "Job.h"                // for Job
#include "PricingSolution.hpp"  // for PricingSolution
struct Machine;

struct Column {
    int    age{};
    bool   del{};
    int    total_processing_time{};
    double total_weighted_completion_time{};

    std::vector<Job*> job_list{};

    Column() = default;
    explicit Column(const Machine& m);
    explicit Column(PricingSolution&& pricing_solution);
    ~Column() = default;
    Column(Column&& _col) = default;
    Column(const Column& _col) = default;

    auto operator<(Column const& other) -> bool;
    auto operator==(Column const& other) const -> bool;
    auto operator<=>(Column const& other) -> bool;
    auto operator=(const Column& other) -> Column& = default;
    auto operator=(Column&& other) -> Column& = default;

    friend auto operator<<(std::ostream& os, const Column& set)
        -> std::ostream& {
        for (const auto& it : set.job_list) {
            os << it->job << " ";
        }
        os << "C = " << set.total_processing_time
           << " cost = " << set.total_weighted_completion_time << '\n';
        return os;
    }
};

// namespace std {
// template <>
// struct less<std::shared_ptr<ScheduleSet>> {
//     constexpr bool operator()(auto const& lhs, auto const& rhs) {
//         return (*lhs) < (*rhs);  // or use boost::hash_combine
//     }
// };
// template <>
// struct equal_to<std::shared_ptr<ScheduleSet>> {
//     constexpr bool operator()(auto const& lhs, auto const& rhs) {
//         return (*lhs) == (*rhs);  // or use boost::hash_combine
//     }
// };

// }  // namespace std

#endif  // SCHEDULESET_H

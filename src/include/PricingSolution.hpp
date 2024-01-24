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

#ifndef OPTIMAL_SOLUTION_HPP
#define OPTIMAL_SOLUTION_HPP
#include <iostream>
#include <range/v3/algorithm/reverse.hpp>
#include <vector>
#include "Job.h"

class PricingSolution {
   public:
    double            obj{};
    int               cost{};
    int               C_max{};
    std::vector<Job*> jobs{};

    /** Default constructor */
    PricingSolution() = default;

    explicit PricingSolution(double _obj) : obj(_obj) {}

    /** Copy constructor */
    PricingSolution(const PricingSolution& other) = delete;

    /** Move constructor */
    PricingSolution(PricingSolution&& other) noexcept = default;

    /** Copy assignment operator */
    auto operator=(const PricingSolution& other) -> PricingSolution& = delete;

    /** Move assignment operator */
    auto operator=(PricingSolution&& other) noexcept
        -> PricingSolution& = default;

    /** Destructor */
    ~PricingSolution() = default;

    friend auto operator<<(std::ostream& os, PricingSolution const& o)
        -> std::ostream& {
        os << "obj = " << o.obj << "," << std::endl
           << "cost = " << o.cost << " C_max = " << o.C_max << std::endl;

        for (const auto& it : o.jobs) {
            os << it->job << ' ';
        }
        return os;
    };

    inline void push_job_back(Job* _job, double _pi) {
        jobs.emplace_back(_job);
        C_max += _job->processing_time;
        cost += _job->weighted_tardiness(C_max);
        obj += _pi;
    }

    inline void push_job_back(Job* _job, int C, double _pi) {
        jobs.emplace_back(_job);
        cost += _job->weighted_tardiness_start(C);
        obj += _pi;
    }

    inline void push_job_back(Job* _job, size_t C, double _pi) {
        jobs.emplace_back(_job);
        cost += _job->weighted_tardiness_start(C);
        obj += _pi;
    }

    void reverse_jobs() { ranges::reverse(jobs); }
};

#endif  // OPTIMAL_SOLUTION_HPP

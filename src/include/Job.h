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

#ifndef JOB_H
#define JOB_H

#include <cstddef>

struct Job {
    size_t job{};
    int    processing_time{};
    int    weight{};
    int    due_time{};

    Job() = default;
    Job(int p, int w, int d);
    Job(const Job&) = default;
    Job(Job&&) = default;
    auto operator=(const Job&) -> Job& = default;
    auto operator=(Job&&) -> Job& = default;
    ~Job() = default;

    [[nodiscard]] auto weighted_tardiness(int C) const -> int;
    [[nodiscard]] auto weighted_tardiness(size_t C) const -> int;
    [[nodiscard]] auto weighted_tardiness_start(int S) const -> int;
    [[nodiscard]] auto weighted_tardiness_start(size_t S) const -> int;
    [[nodiscard]] auto weighted_tardiness_start(long S) const -> int;
};

auto value_diff_Fij(int C, Job* i, Job* j) -> int;
auto value_diff_Fij(size_t C, Job* i, Job* j) -> int;
auto bool_diff_Fij(int weight, Job* _prev, Job* tmp_j) -> bool;

#endif  // JOB_H

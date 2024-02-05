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
#include "Job.h"
#include <algorithm>
#include <cstddef>

Job::Job(int p, int w, int d) : processing_time(p), weight(w), due_time(d) {}

auto Job::weighted_tardiness(int C) const -> int {
    return weight * std::max(0, C - due_time);
}

auto Job::weighted_tardiness(size_t C) const -> int {
    return weight * std::max(0, static_cast<int>(C) - due_time);
}

auto Job::weighted_tardiness_start(int S) const -> int {
    return weight * std::max(0, S + processing_time - due_time);
}

auto Job::weighted_tardiness_start(size_t S) const -> int {
    return weight *
           std::max(0, static_cast<int>(S) + processing_time - due_time);
}

auto Job::weighted_tardiness_start(long S) const -> int {
    return weight *
           std::max(0, static_cast<int>(S) + processing_time - due_time);
}

auto value_diff_Fij(int C, Job* i, Job* j) -> int {
    auto val = i->weighted_tardiness_start(C - j->processing_time);
    val += j->weighted_tardiness(C + i->processing_time);
    val -= j->weighted_tardiness(C);
    val -= i->weighted_tardiness_start(C);
    return val;
}

auto value_diff_Fij(size_t C, Job* i, Job* j) -> int {
    auto j_p = static_cast<size_t>(j->processing_time);
    auto i_p = static_cast<size_t>(i->processing_time);
    auto val = i->weighted_tardiness_start(C - j_p);
    val += j->weighted_tardiness(C + i_p);
    val -= j->weighted_tardiness(C);
    val -= i->weighted_tardiness_start(C);
    return val;
}

auto bool_diff_Fij(int weight, Job* _prev, Job* tmp_j) -> bool {
    return (_prev == nullptr) ? true
                              : (value_diff_Fij(weight + tmp_j->processing_time,
                                                _prev, tmp_j) >= 0);
}

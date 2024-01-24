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
#include "Interval.h"
#include <algorithm>  // for __sort_fn, sort
#include <cmath>
#include <functional>  // for identity, __invoke
#include <memory>      // for shared_ptr, __shared_ptr_access, swap
#include <utility>     // for move
#include "Job.h"       // for Job, value_diff_Fij

Interval::Interval(int _a, int _b, vector_ptr_jobs _jobs)
    : a(_a),
      b(_b),
      sigma(std::move(_jobs)) {
    auto cmp = compare_edd(a, b);
    std::ranges::sort(sigma, cmp);
}

IntervalPair::IntervalPair(int                         _a,
                           int                         _b,
                           const std::shared_ptr<Job>& j1,
                           const std::shared_ptr<Job>& j2,
                           Interval*                   I_)
    : left(_a),
      right(_b),
      jobs{j1, j2},
      I(I_) {}

auto IntervalPair::operator()() -> int {
    left = I->a;
    right = I->a + jobs[1]->processing_time;

    if (left > I->b - jobs[0]->processing_time) {
        return left;
    }

    if (value_diff_Fij(left, jobs[0].get(), jobs[1].get()) <= 0) {
        return left;
    }

    // for (int t = k - 1; t >= 0; t--) {
    //     tmp = (interval*)g_ptr_array_index(interval_array, t);
    //     pair->left = tmp->a + j->processing_time - i->processing_time;

    //     if (value_diff_Fij(pair->left, i, j) <= 0 && pair->left >= tmp->a
    //     &&
    //         pair->left <= tmp->b - i->processing_time) {
    //         break;
    //     }
    // }

    left = jobs[0]->due_time +
           static_cast<int>(std::ceil(
               static_cast<double>(jobs[1]->weight * jobs[0]->processing_time) /
               jobs[0]->weight)) -
           jobs[0]->processing_time;
    return left;
}

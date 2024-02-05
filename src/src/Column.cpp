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
#include "Column.h"             // for ScheduleSet
#include <utility>              // for move
#include "PricingSolution.hpp"  // for PricingSolution
#include "Solution.hpp"         // for Machine

Column::Column(const Machine& m)
    : total_processing_time(m.total_processing_time),
      total_weighted_completion_time(m.total_weighted_tardiness),
      job_list(m.job_list) {}

Column::Column(PricingSolution&& pricing_solution)
    : total_processing_time(pricing_solution.C_max),
      total_weighted_completion_time(pricing_solution.cost),
      job_list(std::move(pricing_solution.jobs)) {}

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

#include "Instance.h"
#include <fmt/core.h>                       // for print
#include <algorithm>                        // for min, max, min_element
#include <cmath>                            // for floor, ceil
#include <cstddef>                          // for size_t
#include <fstream>                          // for basic_istream::opera...
#include <memory>                           // for shared_ptr, __shared...
#include <range/v3/algorithm/for_each.hpp>  // for for_each
#include <range/v3/algorithm/sort.hpp>      // for sort
#include <range/v3/numeric/accumulate.hpp>  // for accumulate, accumula...
#include <range/v3/view/for_each.hpp>       // for for_each
#include <range/v3/view/reverse.hpp>        // for reverse_view, revers...
#include <range/v3/view/sliding.hpp>        // for sliding
#include <range/v3/view/take.hpp>           // for take_view, take, tak...
#include <string>                           // for getline, string
#include <tuple>                            // for tie, operator<=>
#include <utility>                          // for move, pair, make_pair
#include <vector>                           // for vector, erase_if
#include "Interval.h"                       // for IntervalPair, Interval
#include "Job.h"                            // for Job
#include "Parms.h"                          // for Parms
#include "nlohmann/json.hpp"                // for basic_json<>::object_t
#include <fmt/format.h>

Instance::Instance(Parms const& _parms)
    : path_to_instance(_parms.jobfile),
      pname(_parms.pname),
      nb_machines(_parms.nb_machines) {
    std::ifstream in_file{path_to_instance};

    if (in_file.is_open()) {
        std::string str;
        if (getline(in_file, str)) {
            std::istringstream ss(str);
            ss >> nb_jobs;
        }

        int  p{};
        int  d{};
        int  w{};
        auto counter{ 0UL };
        while (getline(in_file, str)) {
            std::istringstream ss(str);
            ss >> p >> d >> w;
            d = d / static_cast<int>(nb_machines);
            jobs.emplace_back(std::make_shared<Job>(p, w, d));
            auto* tmp = jobs.back().get();
            if (tmp->processing_time > tmp->due_time) {
                off += tmp->weight * (tmp->processing_time - tmp->due_time);
                tmp->due_time = tmp->processing_time;
            }

            p_sum += tmp->processing_time;
            pmin = std::min(pmin, tmp->processing_time);
            pmax = std::max(pmax, tmp->processing_time);
            dmin = std::min(dmin, tmp->due_time);
            dmax = std::max(dmax, tmp->due_time);

            tmp->job = counter++;
        }
        calculate_H_max_H_min();
        find_division();
    }

    if(_parms.use_bks.value())
    {
        auto json_filename = fmt::format("sol/{}.json", pname);
        std::ifstream json_file{ json_filename };
        if (json_file.is_open()) {
            nlohmann::json json;
            json_file >> json;
            std::string nb_machines_str = std::to_string(nb_machines);
            auto sol = json[nb_machines_str]["bks"].get<int>();
            best_known_sol = sol;
            json_file.close();
        }
        else {
            throw InstanceException("Could not open file\n");
        }
    }

    
}

void Instance::calculate_H_max_H_min() {
    auto temp = p_sum - pmax;
    auto temp_dbl = static_cast<double>(temp);
    auto tmp_m = static_cast<double>(nb_machines);
    temp_dbl = std::floor(temp_dbl / tmp_m);

    H_max = static_cast<int>(temp_dbl) + pmax;
    H_min = static_cast<int>(std::ceil(temp_dbl / tmp_m)) - pmax;

    ranges::sort(jobs, std::less{},
                 [](const auto& lhs) { return lhs->processing_time; });

    H_min = p_sum;
    auto tmp = ranges::accumulate(
        jobs | ranges::views::reverse | ranges::views::take(nb_machines - 1),
        p_sum, std::minus<>{}, [](auto& it) { return it->processing_time; });

    H_min = static_cast<int>(ceil(tmp / tmp_m));
    fmt::print(
        R"(H_max = {}, H_min = {},  pmax = {}, pmin = {}, p_sum = {}, off = {}
)",
        H_max, H_min, pmax, pmin, p_sum, off);

    auto projection = [](const std::shared_ptr<Job>& x) {
        return std::tie(x->due_time, x->processing_time, x->weight, x->job);
    };
    ranges::sort(jobs, std::less<>{}, projection);

    auto index = 0UL;
    ranges::for_each(jobs, [&index](auto& it) { it->job = index++; });
}

void Instance::find_division() {
    int prev = 0;

    std::vector<std::shared_ptr<Interval>> tmp_ptr_vec{};
    for (auto i = 0UL; i < nb_jobs && prev < H_max; ++i) {
        auto tmp = std::min(H_max, jobs[i]->due_time);
        if (prev < tmp) {
            tmp_ptr_vec.emplace_back(
                std::make_shared<Interval>(prev, tmp, jobs));
            prev = jobs[i]->due_time;
        }
    }

    if (prev < H_max) {
        tmp_ptr_vec.emplace_back(std::make_shared<Interval>(prev, H_max, jobs));
    }

    for (auto& it : tmp_ptr_vec) {
        std::vector<IntervalPair> special_pairs{};
        for (auto i = 0U; i < it->sigma.size() - 1; i++) {
            auto job1 = it->sigma[i];
            for (size_t j = i + 1; j < it->sigma.size(); j++) {
                auto job2 = it->sigma[j];
                auto intervalpair =
                    IntervalPair(it->a, it->b, job1, job2, it.get());

                if (job1->due_time < it->b && job2->due_time < it->b &&
                    !((it->a + job2->processing_time >= it->b) ||
                      (intervalpair() <= it->a))) {
                    special_pairs.push_back(intervalpair);
                }
            }
        }

        if (!special_pairs.empty()) {
            std::vector<int> t{};
            t.push_back(it->a);

            while (!special_pairs.empty()) {
                auto min = std::min_element(
                    special_pairs.begin(), special_pairs.end(),
                    [](const auto& lhs, const auto& rhs) -> bool {
                        return (lhs.right < rhs.right);
                    });

                t.push_back(min->right);
                special_pairs.erase(min);

                std::erase_if(special_pairs, [&t](auto& it_tmp) -> bool {
                    auto tmp = t.back();
                    if ((tmp >= it_tmp.left && tmp <= it_tmp.right) ||
                        (tmp + it_tmp.jobs[1]->processing_time >=
                         it_tmp.I->b)) {
                        return true;
                    }
                    it_tmp.right = tmp + it_tmp.jobs[1]->processing_time;
                    return false;
                });
            }

            t.push_back(it->b);

            for (auto&& tmp : t | ranges::views::sliding(2)) {
                intervals.push_back(
                    std::make_shared<Interval>(tmp[0], tmp[1], jobs));
            }
        } else {
            intervals.push_back(std::make_shared<Interval>(it->a, it->b, jobs));
        }
    }

    for (auto& it : intervals) {
        for (auto& job : it->sigma) {
            if (job->processing_time <= it->b) {
                vector_pair.emplace_back(std::make_pair(job.get(), it.get()));
            }
        }
    }

    fmt::print(R"(The number of layers = {}
)",
vector_pair.size());
}

void Instance::print_intervals() {
    fmt::print(R"(The number of jobs = {}
The number of machines = {}
The number of intervals = {}
)",
nb_jobs, nb_machines, intervals.size());

    for (auto& it : intervals) {
        fmt::print(R"(Interval: [{}, {}]
)",
it->a, it->b);
        for (auto& job : it->sigma) {
            fmt::print(R"(Job: {}, processing time = {}, due time = {}, weight = {}
)",
job->job, job->processing_time, job->due_time, job->weight);
        }
    }
}


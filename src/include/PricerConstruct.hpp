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

#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <cassert>                                  // for assert
#include <cstddef>                                  // for size_t
#include <memory>                                   // for allocator_traits<...
#include <utility>                                  // for pair
#include <vector>                                   // for vector, allocator
#include "Instance.h"                               // for Instance
#include "Interval.h"                               // for Interval
#include "Job.h"                                    // for Job
#include "ModernDD/NodeBddSpec.hpp"                 // for DdSpec

class conflict_state {
   public:
    boost::dynamic_bitset<> add;
    boost::dynamic_bitset<> remove;

    conflict_state(const conflict_state&) = default;
    conflict_state(conflict_state&&) = default;
    conflict_state() = default;
    ~conflict_state() = default;

    auto operator=(const conflict_state&) -> conflict_state& = default;
    auto operator=(conflict_state&&) -> conflict_state& = default;
};

class PricerConstruct : public DdSpec<PricerConstruct, int, 2> {
    const std::vector<std::pair<Job*, Interval*>>* ptr_vector{};
    size_t                                         nb_layers;

   public:
    explicit PricerConstruct(const Instance& instance)
        : ptr_vector(&instance.vector_pair),
          nb_layers(instance.vector_pair.size()) {}

    PricerConstruct(const PricerConstruct&) = default;
    PricerConstruct(PricerConstruct&&) = default;
    ~PricerConstruct() = default;
    auto operator=(const PricerConstruct&) -> PricerConstruct& = default;
    auto operator=(PricerConstruct&&) -> PricerConstruct& = default;

    auto getRoot(int& state) const -> int {
        state = 0;
        return static_cast<int>(nb_layers);
    };

    auto getChild(int& state, int level, size_t value) const -> int {
        auto layer = nb_layers - level;
        // assert(0 <= layer && layer <= nb_layers - 1);
        // auto* tmp_pair = static_cast<job_interval_pair*>(aux_list[layer]);
        const auto& tmp_pair = (*ptr_vector)[layer];
        auto*       tmp_j = tmp_pair.first;
        auto*       tmp_interval = tmp_pair.second;

        if (value == 1) {
            state = state + tmp_j->processing_time;
        }

        auto _j = min_job(layer, state, value);

        if (_j >= nb_layers) {
            if (state <= tmp_interval->b) {
                return -1;
            }
            return 0;
        }
        assert(_j < nb_layers);
        return static_cast<int>(nb_layers - _j);
    }

   private:
    [[nodiscard]] auto min_job(size_t j, int state, size_t value) const
        -> size_t {
        auto val = nb_layers;
        // auto* tmp = static_cast<job_interval_pair*>(aux_list[j])->j;
        auto* tmp = (*ptr_vector)[j].first;

        if (value != 0) {
            for (size_t i = j + 1; i < nb_layers; ++i) {
                const auto& tmp_pair = (*ptr_vector)[i];
                auto*       tmp_j = tmp_pair.first;
                auto*       tmp_interval = tmp_pair.second;
                // auto* tmp_pair =
                // static_cast<job_interval_pair*>(aux_list[i]); auto*
                // tmp_interval = static_cast<interval*>(tmp_pair->I); auto*
                // tmp_j = static_cast<Job*>(tmp_pair->j);

                if (state + tmp_j->processing_time > tmp_interval->a &&
                    state + tmp_j->processing_time <= tmp_interval->b) {
                    if (tmp == tmp_j) {
                        continue;
                    }
                    val = i;
                    break;
                }
            }
        } else {
            for (auto i = j + 1; i < nb_layers; ++i) {
                const auto& tmp_pair = (*ptr_vector)[i];
                auto*       tmp_j = tmp_pair.first;
                auto*       tmp_interval = tmp_pair.second;
                // auto* tmp_pair =
                // static_cast<job_interval_pair*>(aux_list[i]); auto*
                // tmp_interval = static_cast<interval*>(tmp_pair->I); auto*
                // tmp_j = static_cast<Job*>(tmp_pair->j);

                if (state + tmp_j->processing_time > tmp_interval->a &&
                    state + tmp_j->processing_time <= tmp_interval->b) {
                    val = i;
                    break;
                }
            }
        }

        return val;
    }

    static auto diff_obj(Job* i, Job* j, int C) -> int {
        return i->weighted_tardiness(C) +
               j->weighted_tardiness(C + j->processing_time) -
               (j->weighted_tardiness(C - i->processing_time +
                                      j->processing_time) +
                i->weighted_tardiness(C + j->processing_time));
    }
};

// class PricerConstructTI : public DdSpec<PricerConstructTI, int, 2> {
//     std::span<void*> aux_list;
//     // GPtrArray*       pair_list;
//     int* take_job;
//     int  Hmax;
//     int  nb_layers;

//    public:
//     explicit PricerConstructTI(GPtrArray* _pair_list, int* _take_job, int
//     _Hmax)
//         : aux_list{_pair_list->pdata, _pair_list->len},
//           take_job(_take_job),
//           Hmax(_Hmax),
//           nb_layers(_pair_list->len){};

//     PricerConstructTI(const PricerConstructTI&) = default;
//     PricerConstructTI(PricerConstructTI&&) = default;
//     PricerConstructTI& operator=(const PricerConstructTI&) = default;
//     PricerConstructTI& operator=(PricerConstructTI&&) = default;

//     int getRoot(int& state) {
//         state = 0;
//         return nb_layers;
//     };

//     int getChild(int& state, int level, int value) const {
//         auto layer = nb_layers - level;
//         // assert(0 <= layer && layer <= nb_layers - 1);
//         auto* tmp_pair = static_cast<job_interval_pair*>(aux_list[layer]);
//         auto* tmp_interval = static_cast<interval*>(tmp_pair->I);
//         auto* tmp_j = static_cast<Job*>(tmp_pair->j);

//         if (level - 1 == 0 && value) {
//             return (state + tmp_j->processing_time <= tmp_interval->b) ? -1 :
//             0;
//         } else if (level - 1 == 0) {
//             return (state <= tmp_interval->b) ? -1 : 0;
//         }

//         if (value) {
//             state = state + tmp_j->processing_time;
//         }

//         auto _j = min_job(layer, state, value);

//         if (!(_j < nb_layers)) {
//             if (state <= tmp_interval->b) {
//                 return -1;
//             }
//             return 0;
//         }
//         assert(_j < nb_layers);
//         return nb_layers - _j;
//     }

//     ~PricerConstructTI() = default;

//    private:
//     [[nodiscard]] int min_job(int j, int state, int value) const {
//         int                val = nb_layers;
//         job_interval_pair* tmp_pair = nullptr;
//         interval*          tmp_interval = nullptr;
//         Job*               tmp_j = nullptr;
//         auto* tmp = (static_cast<job_interval_pair*>(aux_list[j]))->j;

//         if (value) {
//             for (int i = j + 1; i < nb_layers; ++i) {
//                 tmp_pair = static_cast<job_interval_pair*>(aux_list[i]);
//                 tmp_interval = tmp_pair->I;
//                 tmp_j = tmp_pair->j;

//                 if (state + tmp_j->processing_time > tmp_interval->a &&
//                     state + tmp_j->processing_time <= tmp_interval->b &&
//                     take_job[tmp_j->job * (Hmax + 1) + state]) {
//                     if (tmp == tmp_j ||
//                         (tmp->job > tmp_j->job &&
//                          value_diff_Fij(state, tmp_j, tmp) < 0)) {
//                         continue;
//                     }
//                     val = i;
//                     break;
//                 }
//             }
//         } else {
//             for (int i = j + 1; i < nb_layers; ++i) {
//                 tmp_pair = static_cast<job_interval_pair*>(aux_list[i]);
//                 tmp_interval = tmp_pair->I;
//                 tmp_j = tmp_pair->j;

//                 if (state + tmp_j->processing_time > tmp_interval->a &&
//                     state + tmp_j->processing_time <= tmp_interval->b &&
//                     take_job[tmp_j->job * (Hmax + 1) + state]) {
//                     val = i;
//                     break;
//                 }
//             }
//         }

//         return val;
//     }

//     int diff_obj(Job* i, Job* j, int C) const {
//         return value_Fj(C, i) + value_Fj(C + j->processing_time, j) -
//                (value_Fj(C - i->processing_time + j->processing_time, j) +
//                 value_Fj(C + j->processing_time, i));
//     }
// };

// class ConflictConstraints
//     : public DdSpec<ConflictConstraints, conflict_state, 2> {
//     size_t                               nb_jobs;
//     std::vector<boost::dynamic_bitset<>> differsets;
//     std::vector<boost::dynamic_bitset<>> samesets;

//     bool takeable(size_t job, conflict_state& state) {
//         if (state.remove[job]) {
//             return false;
//         }

//         return true;
//     }

//     bool leaveable(size_t job, conflict_state& state) {
//         if (state.add[job]) {
//             return false;
//         }

//         return true;
//     }

//    public:
//     ConflictConstraints(size_t _nb_jobs,
//                         int*   elist_same,
//                         size_t edge_count_same,
//                         int*   elist_differ,
//                         size_t edge_count_differ)
//         : nb_jobs(_nb_jobs) {
//         std::span aux_same{elist_same, edge_count_same};
//         std::span aux_diff{elist_differ, edge_count_differ};

//         differsets.resize(_nb_jobs);
//         samesets.resize(_nb_jobs);

//         for (auto i = 0UL; i < _nb_jobs; i++) {
//             differsets[i].resize(_nb_jobs);
//             samesets[i].resize(_nb_jobs);
//         }

//         for (auto i = 0UL; i < edge_count_same; ++i) {
//             samesets[aux_same[2 * i]][aux_same[2 * i + 1]] = true;
//         }

//         for (auto i = 0UL; i < edge_count_differ; ++i) {
//             differsets[aux_diff[2 * i]][aux_diff[2 * i + 1]] = true;
//         }
//     };

//     ~ConflictConstraints() = default;
//     ConflictConstraints(const ConflictConstraints&) = default;
//     ConflictConstraints(ConflictConstraints&&) = default;
//     ConflictConstraints& operator=(ConflictConstraints&&) = default;
//     ConflictConstraints& operator=(const ConflictConstraints&) = default;

//     int getRoot(conflict_state& state) const {
//         state.add.resize(nb_jobs);
//         state.remove.resize(nb_jobs);
//         return nb_jobs;
//     }

//     int getChild(conflict_state& state, int level, int take) {
//         size_t job = nb_jobs - level;
//         size_t _j = 0UL;
//         assert(0UL <= job && job <= nb_jobs - 1);

//         if (samesets[job].intersects(differsets[job])) {
//             return 0;
//         }

//         if (level - 1 == 0 && take) {
//             return (!state.remove[job]) ? -1 : 0;
//         } else if (level - 1 == 0) {
//             return (!state.add[job]) ? -1 : 0;
//         }

//         if (take) {
//             if (!takeable(job, state)) {
//                 return 0;
//             }

//             state.add |= samesets[job];
//             state.remove |= differsets[job];
//         } else {
//             if (!leaveable(job, state)) {
//                 return 0;
//             }

//             state.remove |= samesets[job];
//         }

//         _j = min_job(job, &state);

//         if (_j == nb_jobs && take) {
//             return (!state.remove[job]) ? -1 : 0;
//         } else if (_j == nb_jobs) {
//             return (!state.add[job]) ? -1 : 0;
//         }

//         assert(_j < nb_jobs);
//         return nb_jobs - _j;
//     }

//     [[nodiscard]] bool equalTo(conflict_state const& state1,
//                                conflict_state const& state2) const {
//         if (state2.add != state1.add) {
//             return false;
//         }

//         if (state2.remove != state1.remove) {
//             return false;
//         }

//         return true;
//     }

//     [[nodiscard]] size_t hashCode(conflict_state const& state) const {
//         size_t val = 0;
//         size_t it = state.add.find_first();

//         while (it != boost::dynamic_bitset<>::npos) {
//             // val += 1213657 * it;
//             boost::hash_combine(val, it);
//             it = state.add.find_next(it);
//         }

//         it = state.remove.find_first();

//         while (it != boost::dynamic_bitset<>::npos) {
//             // val += 487239 * it;
//             boost::hash_combine(val, it);
//             it = state.remove.find_next(it);
//         }

//         return val;
//     }

//     size_t min_job(size_t j, conflict_state* state) const {
//         auto val = nb_jobs;

//         for (auto i = j + 1; i < nb_jobs; ++i) {
//             if (!state->remove[i]) {
//                 val = i;
//                 break;
//             }
//         }

//         return val;
//     }
// };

// class scheduling : public DdSpec<scheduling, int, 2> {
//     Job*                                           job;
//     const std::vector<std::pair<Job*, Interval*>>* ptr_vector;
//     // GPtrArray*                                     list_layers;
//     // std::span<void*>                               aux_list;
//     int nb_layers;
//     int order;

//    public:
//     explicit scheduling(Job*                                     _job,
//                         std::vector<std::pair<Job*, Interval*>>&
//                         _vector_pair, int _order)
//         : job(_job),
//           ptr_vector(&_vector_pair),
//           nb_layers(_vector_pair.size()),
//           order(_order) {}

//     scheduling(const scheduling&) = default;
//     scheduling(scheduling&&) = default;
//     scheduling& operator=(const scheduling&) = default;
//     scheduling& operator=(scheduling&&) = default;
//     ~scheduling() = default;

//     int getRoot(int& state) const {
//         state = 0;
//         return nb_layers;
//     }

//     int getChild(int& state, int level, int take) {
//         auto j = nb_layers - level;
//         // assert(0 <= j && j <= nb_layers - 1);
//         auto& tmp = (*ptr_vector)[j];
//         auto* tmp_j = tmp.first;

//         if (level - 1 == 0 && take) {
//             if (tmp_j != job) {
//                 return -1;
//             } else {
//                 if (state == 0) {
//                     return -1;
//                 }
//                 return 0;
//             }
//         } else if (level - 1 == 0) {
//             return -1;
//         }

//         if (take) {
//             if (tmp_j != job) {
//                 // state++;
//                 j++;
//                 return nb_layers - j;
//             } else {
//                 if (state == 0) {
//                     state++;
//                     j++;
//                     return nb_layers - j;
//                 } else if (state >= order) {
//                     state = 1;
//                     j++;
//                     return nb_layers - j;
//                 } else {
//                     state++;
//                     return 0;
//                 }
//             }
//         } else {
//             j++;
//             return nb_layers - j;
//         }
//     }
// };

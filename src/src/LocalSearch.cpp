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

#include "LocalSearch.hpp"
#include <fmt/core.h>
#include <limits>
#include <range/v3/action/shuffle.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/drop.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/take.hpp>
#include <vector>
#include "DebugLvl.hpp"
#include "Job.h"
#include "Parms.h"
#include "Solution.hpp"

LocalSearchData::LocalSearchData(size_t _nb_jobs, size_t _nb_machines)
    : nb_jobs(_nb_jobs),
      nb_machines(_nb_machines),
      W(_nb_machines, VectorInt(_nb_jobs)),
      g(_nb_machines, std::vector<SlopeList>(_nb_jobs)),
      processing_list{
          MatrixProcessingList(
              _nb_machines,
              VectorProcessingList(_nb_jobs, ProcessingListData())),
          MatrixProcessingList(
              _nb_machines,
              VectorProcessingList(_nb_jobs, ProcessingListData()))},
      G(_nb_jobs, VectorInt(_nb_jobs, 0)),
      H(_nb_jobs, VectorInt(_nb_jobs, 0)),
      GG(_nb_jobs, 0),
      HH(_nb_jobs, VectorInt(_nb_jobs, 0)),
      begin(_nb_jobs, SlopeListIt()),
      end(_nb_jobs, SlopeListIt()),
      B2_1(_nb_jobs + 1, VectorInt(_nb_jobs + 1)),
      B2_2(_nb_jobs + 1, VectorInt(_nb_jobs + 1)),
      B3_1(_nb_jobs + 1, VectorInt(_nb_jobs + 1)),
      B3_2(_nb_jobs + 1, VectorInt(_nb_jobs + 1)),
      B4_1(_nb_jobs + 1, VectorInt(_nb_jobs + 1)),
      B4_2(_nb_jobs + 1, VectorInt(_nb_jobs + 1)),
      B3_1_(_nb_jobs + 1, 0),
      B5_1(_nb_jobs + 1, VectorInt(_nb_jobs + 1)),
      B5_2(_nb_jobs + 1, VectorInt(_nb_jobs + 1)),
      B6_1(_nb_jobs + 1, VectorInt(_nb_jobs + 1)) {}

void LocalSearchData::RVND(Sol& sol) {
    calculate_W(sol);
    calculate_g(sol);
    std::random_device         rd;
    std::default_random_engine rng(rd());

    do {
        moves |= ranges::actions::shuffle(rng);
        for (auto& operator_move : moves) {
            operator_move(sol);
            if (updated) {
                break;
            }
        }
    } while (updated);
}

void LocalSearchData::insertion_operator(Sol& sol, size_t l) {
    int         max = 0;
    size_t      i_best{};
    size_t      j_best{};
    size_t      k_best{};
    int         c{};
    SlopeListIt it;
    SlopeListIt tmp_end;

    ++iterations;
    create_processing_insertion_operator(sol, l);
    updated = false;

    for (auto&& [k, m] : sol.machines | ranges::views::enumerate) {
        const auto& vec_m = m.job_list;
        const auto& n = vec_m.size();

        for (auto i = 0UL; i < n - l; i++) {
            auto& p_list = processing_list[0][k][i];
            it = g[k][p_list.pos].begin();
            tmp_end = g[k][p_list.pos].end();

            for (auto j = p_list.pos + l; j < n; j++) {
                auto start = sol.c[vec_m[j]->job] - p_list.p;
                it = std::find_if(it, tmp_end,
                                  [&start](const auto& tmp) -> bool {
                                      return tmp.in_interval(start);
                                  });
                G[p_list.pos][j] = it != tmp_end ? (*it)(start) : 0;
            }
        }

        for (auto i = 0UL; i < n - l; i++) {
            it = g[k][i + l].begin();
            tmp_end = g[k][i + l].end();

            for (auto j = i + l; j < n; j++) {
                c = sol.c[vec_m[j]->job];
                it = std::find_if(it, tmp_end, [&c](const auto& tmp) -> bool {
                    return tmp.in_interval(c);
                });
                H[i][j] = it != tmp_end ? (*it)(c) : 0;
            }
        }

        for (auto i = 0UL; i < n - l; i++) {
            it = g[k][i + l].begin();
            tmp_end = g[k][i + l].end();
            c = 0;

            if (i != 0) {
                c = sol.c[vec_m[i - 1]->job];
            }
            it = std::find_if(it, tmp_end, [&c](const auto& tmp) -> bool {
                return tmp.in_interval(c);
            });
            GG[i] = it != tmp_end ? (*it)(c) : 0;
        }

        for (auto j = l; j < n - 1; ++j) {
            begin[j] = g[k][j + 1].begin();
            end[j] = g[k][j + 1].end();
        }

        for (auto i = 0UL; i < n - l; ++i) {
            auto pos = processing_list[0][k][i].pos;
            auto p = processing_list[0][k][i].p;

            for (auto j = pos + l; j < n; ++j) {
                if (j != n - 1) {
                    c = sol.c[vec_m[j]->job] - p;
                    begin[j] = std::find_if(begin[j], end[j],
                                            [&c](const auto& tmp) -> bool {
                                                return tmp.in_interval(c);
                                            });
                    HH[pos][j] = begin[j] != end[j] ? (*begin[j])(c) : 0;
                } else {
                    HH[pos][j] = 0;
                }
            }
        }

        for (auto i = 0UL; i < n - l; i++) {
            for (auto j = i + l; j < n; j++) {
                auto tmp = 0;

                if (i != 0) {
                    tmp += W[k][i - 1];
                }

                tmp += G[i][j] - H[i][j];
                tmp += GG[i] - HH[i][j];
                tmp += W[k][n - 1] - W[k][j];

                if (m.total_weighted_tardiness - tmp > max) {
                    max = m.total_weighted_tardiness - tmp;
                    i_best = i;
                    j_best = j;
                    k_best = k;
                    updated = true;
                }
            }
        }

        // k++;
    }

    if (updated) {
        if (debug_lvl(0)) {
            sol.print_solution();
        }

        sol.update_insertion_move(i_best, j_best, k_best, l);
        calculate_W(sol);
        calculate_g(sol);

        if (debug_lvl(0)) {
            sol.print_solution();
        }
    }
}

void LocalSearchData::swap_operator(Sol& sol, size_t l1, size_t l2) {
    boost::timer::cpu_timer time_swap_operator;

    auto        max = 0;
    size_t      i_best{};
    size_t      j_best{};
    size_t      k_best{};
    SlopeListIt it;
    SlopeListIt end_tmp;
    int         c{};

    ++iterations;
    create_processing_list_swap(sol, l1, l2);
    updated = false;

    for (auto&& [k, m] : sol.machines | ranges::views::enumerate) {
        auto&       machine = m.job_list;
        const auto& n = machine.size();

        if (n + 1 <= l1 + l2) {
            continue;
        }

        /** compute g */
        for (auto& p_list :
             processing_list[0][k] | ranges::views::take(n - l1 - l2 + 1)) {
            // for (auto i = 0UL; i < n - l1 - l2 + 1; ++i) {
            // auto& p_list = it;
            it = g[k][p_list.pos].begin();
            end_tmp = g[k][p_list.pos].end();

            for (auto j = p_list.pos + l1; j < n - l2 + 1; ++j) {
                c = sol.c[machine[j + l2 - 1]->job] - p_list.p;
                it = std::find_if(it, end_tmp, [&c](const auto& tmp) -> bool {
                    return tmp.in_interval(c);
                });
                B2_1[p_list.pos][j] = it != end_tmp ? (*it)(c) : 0;
            }
        }

        /** compute h */
        for (auto i = 0UL; i + l1 + l2 < n + 1; ++i) {
            if (i + l1 >= n) {
                for (auto j = i + l1; j + l2 < n + 1; ++j) {
                    B2_2[i][j] = 0;
                }
            } else {
                it = g[k][i + l1].begin();
                end_tmp = g[k][i + l1].end();
                for (auto j = i + l1; j + l2 < n + 1; ++j) {
                    c = sol.c[machine[j + l2 - 1]->job];
                    it = std::find_if(it, end_tmp,
                                      [&c](const auto& tmp) -> bool {
                                          return tmp.in_interval(c);
                                      });
                    B2_2[i][j] = it != end_tmp ? (*it)(c) : 0;
                }
            }
        }

        /** compute B3_1 */
        for (auto i = 0UL; i + l1 + l2 < n + 1; ++i) {
            begin[i] = g[k][i + l1].begin();
            end[i] = g[k][i + l1].end();
        }

        for (auto j = l1; j + l2 < n + 1; ++j) {
            auto& p_list = processing_list[1][k][j - l1];

            for (auto i = 0UL; i + l1 < p_list.pos + 1; ++i) {
                if (i + l1 >= n) {
                    B3_1[i][p_list.pos] = 0;
                } else {
                    c = p_list.p;

                    if (i != 0) {
                        c += sol.c[machine[i - 1]->job];
                    }

                    begin[i] = std::find_if(begin[i], end[i],
                                            [&c](const auto& tmp) -> bool {
                                                return tmp.in_interval(c);
                                            });
                    B3_1[i][p_list.pos] =
                        begin[i] != end[i] ? (*begin[i])(c) : 0;
                }
            }
        }

        /** compute B3_2 */
        for (auto j = l1; j + l2 < n + 1; ++j) {
            begin[j] = g[k][j].begin();
            end[j] = g[k][j].end();
        }

        for (auto i = 0UL; i + l1 + l2 < n + 1; ++i) {
            auto& p_list = processing_list[0][k][i];

            for (auto j = p_list.pos + l1; j < n - l2 + 1; ++j) {
                if (i + l1 >= n) {
                    B3_2[p_list.pos][j] = 0;
                } else {
                    c = sol.c[machine[j + l2 - 1]->job] - p_list.p;
                    begin[j] = std::find_if(begin[j], end[j],
                                            [&c](const auto& tmp) -> bool {
                                                return tmp.in_interval(c);
                                            });
                    B3_2[p_list.pos][j] =
                        begin[j] != end[j] ? (*begin[j])(c) : 0;
                }
            }
        }

        /** compute B4_1 */
        for (auto j = l1; j + l2 < n + 1; ++j) {
            it = g[k][j].begin();
            end_tmp = g[k][j].end();

            for (auto i = 0UL; i + l1 < j + 1; ++i) {
                c = 0;

                if (i != 0) {
                    c = sol.c[machine[i - 1]->job];
                }

                it = std::find_if(it, end_tmp, [&c](const auto& tmp) -> bool {
                    return tmp.in_interval(c);
                });
                B4_1[i][j] = it != end_tmp ? (*it)(c) : 0;
            }
        }

        /** compute B4_2 */
        for (auto j = l1; j + l2 < n + 1; ++j) {
            auto& p_list = processing_list[1][k][j - l1];

            if (p_list.pos + l2 == n) {
                for (auto i = 0UL; i + l1 < p_list.pos + 1; ++i) {
                    B4_2[i][p_list.pos] = 0;
                }
            } else {
                it = g[k][p_list.pos + l2].begin();
                end_tmp = g[k][p_list.pos + l2].end();

                for (auto i = 0UL; i + l1 < p_list.pos + 1; ++i) {
                    c = p_list.p;

                    if (i != 0) {
                        c += sol.c[machine[i - 1]->job];
                    }

                    it = std::find_if(it, end_tmp,
                                      [&c](const auto& tmp) -> bool {
                                          return tmp.in_interval(c);
                                      });

                    B4_2[i][p_list.pos] = it != end_tmp ? (*it)(c) : 0;
                }
            }
        }

        for (auto i = 0UL; i + l1 + l2 < n + 1; ++i) {
            for (auto j = i + l1; j + l2 < n + 1; ++j) {
                auto t = 0;

                if (i != 0) {
                    t += W[k][i - 1];
                }

                t += B2_1[i][j] - B2_2[i][j];
                t += B3_1[i][j] - B3_2[i][j];
                t += B4_1[i][j] - B4_2[i][j];
                t += W[k][n - 1] - W[k][j + l2 - 1];

                if (m.total_weighted_tardiness - t > max) {
                    max = m.total_weighted_tardiness - t;
                    i_best = i;
                    j_best = j;
                    k_best = k;
                    updated = true;
                }
            }
        }
    }

    /** update to best improvement */
    if (updated) {
        if (debug_lvl(0)) {
            sol.print_solution();
        }

        sol.update_swap_move(i_best, j_best, k_best, l1, l2);
        calculate_W(sol);
        calculate_g(sol);

        if (debug_lvl(0)) {
            sol.print_solution();
        }
    }

    if (debug_lvl(0)) {
        constexpr int PRECISION = 5;
        fmt::print(
            R"(intra swap with l1 = {} and l2 = {}, running time = {} and improvement {} on machine {} on places {} {}
)",
            l1, l2, time_swap_operator.format(PRECISION, "%ws"), max, k_best, i_best,
            j_best);
    }
}

void LocalSearchData::swap_operator_inter(Sol& sol, size_t l1, size_t l2) {
    auto        max = 0;
    size_t      i_best{};
    size_t      j_best{};
    size_t      k_best{};
    size_t      kk_best{};
    SlopeListIt it;
    SlopeListIt end_tmp;
    int         c{};

    ++iterations;
    create_processing_list_swap_inter(sol, l1, l2);
    updated = false;

    for (auto k1 : ranges::views::iota(size_t{}, nb_machines)) {
        auto&      machine1 = sol.machines[k1].job_list;
        const auto nb_jobs1 = machine1.size();

        for (auto k2 : ranges::views::iota(size_t{}, nb_machines)) {
            auto&      machine2 = sol.machines[k2].job_list;
            const auto nb_jobs2 = machine2.size();

            if (k1 == k2 || nb_jobs1 < l1 || nb_jobs2 < l2) {
                continue;
            }

            /** compute B2_1 */
            for (auto i = 0UL; i + l1 < nb_jobs1 + 1; ++i) {
                it = g[k1][i].begin();
                end_tmp = g[k1][i].end();

                for (auto j = 0UL; j + l2 < nb_jobs2 + 1; ++j) {
                    c = 0;

                    if (j != 0) {
                        c = sol.c[machine2[j - 1]->job];
                    }

                    it = std::find_if(it, end_tmp,
                                      [&c](const auto& tmp) -> bool {
                                          return tmp.in_interval(c);
                                      });
                    B2_1[i][j] = it != end_tmp ? (*it)(c) : 0;
                }
            }

            for (auto i = 0UL; i + l1 < nb_jobs1 + 1; ++i) {
                auto& p_list = processing_list[0][k1][i];

                if (p_list.pos + l1 >= nb_jobs1) {
                    for (auto j = 0UL; j + l2 < nb_jobs2 + 1; ++j) {
                        B2_2[p_list.pos][j] = 0;
                    }
                } else {
                    it = g[k1][p_list.pos + l1].begin();
                    end_tmp = g[k1][p_list.pos + l1].end();
                    for (auto j = 0UL; j + l2 < nb_jobs2 + 1; ++j) {
                        c = p_list.p;

                        if (j != 0) {
                            c += sol.c[machine2[j - 1]->job];
                        }

                        it = std::find_if(it, end_tmp,
                                          [&c](const auto& tmp) -> bool {
                                              return tmp.in_interval(c);
                                          });
                        B2_2[p_list.pos][j] = it != end_tmp ? (*it)(c) : 0;
                    }
                }
            }

            for (auto i = 0UL; i + l1 < nb_jobs1; ++i) {
                begin[i] = g[k1][i + l1].begin();
                end[i] = g[k1][i + l1].end();
            }

            for (auto j = 0UL; j + l2 < nb_jobs2 + 1; ++j) {
                auto& p_list = processing_list[1][k2][j];

                for (auto i = 0UL; i + l1 < nb_jobs1 + 1; ++i) {
                    if (i + l1 >= nb_jobs1) {
                        B3_1[i][p_list.pos] = 0;
                    } else {
                        c = p_list.p;

                        if (i != 0) {
                            c += sol.c[machine1[i - 1]->job];
                        }

                        begin[i] = std::find_if(begin[i], end[i],
                                                [&c](const auto& tmp) -> bool {
                                                    return tmp.in_interval(c);
                                                });
                        B3_1[i][p_list.pos] =
                            begin[i] != end[i] ? (*begin[i])(c) : 0;
                    }
                }
            }

            for (auto j = 0UL; j + l2 < nb_jobs2 + 1; ++j) {
                it = g[k2][j].begin();
                end_tmp = g[k2][j].end();

                for (auto i = 0UL; i + l1 < nb_jobs1 + 1; ++i) {
                    c = 0;

                    if (i != 0) {
                        c += sol.c[machine1[i - 1]->job];
                    }

                    it = std::find_if(it, end_tmp,
                                      [&c](const auto& tmp) -> bool {
                                          return tmp.in_interval(c);
                                      });
                    B5_1[i][j] = it != end_tmp ? (*it)(c) : 0;
                }
            }

            for (auto j = 0UL; j + l2 < nb_jobs2 + 1; ++j) {
                auto& p_list = processing_list[1][k2][j];

                if (p_list.pos + l2 >= nb_jobs2) {
                    for (auto i = 0UL; i + l1 < nb_jobs1 + 1; ++i) {
                        B5_2[i][p_list.pos] = 0;
                    }
                } else {
                    it = g[k2][p_list.pos + l2].begin();
                    end_tmp = g[k2][p_list.pos + l2].end();
                    for (auto i = 0UL; i + l1 < nb_jobs1 + 1; ++i) {
                        c = p_list.p;

                        if (i != 0) {
                            c += sol.c[machine1[i - 1]->job];
                        }

                        it = std::find_if(it, end_tmp,
                                          [&c](const auto& tmp) -> bool {
                                              return tmp.in_interval(c);
                                          });
                        B5_2[i][p_list.pos] = it != end_tmp ? (*it)(c) : 0;
                    }
                }
            }

            for (auto j = 0UL; j + l2 < nb_jobs2; ++j) {
                begin[j] = g[k2][j + l2].begin();
                end[j] = g[k2][j + l2].end();
            }

            for (auto i = 0UL; i + l1 < nb_jobs1 + 1; ++i) {
                auto& p_list = processing_list[0][k1][i];

                for (auto j = 0UL; j + l2 < nb_jobs2 + 1; ++j) {
                    if (j + l2 >= nb_jobs2) {
                        B6_1[p_list.pos][j] = 0;
                    } else {
                        c = p_list.p;

                        if (j != 0) {
                            c += sol.c[machine2[j - 1]->job];
                        }

                        begin[j] = std::find_if(begin[j], end[j],
                                                [&c](const auto& tmp) -> bool {
                                                    return tmp.in_interval(c);
                                                });
                        B6_1[p_list.pos][j] =
                            begin[j] != end[j] ? (*begin[j])(c) : 0;
                    }
                }
            }

            for (auto i = 0UL; i + l1 < nb_jobs1 + 1; ++i) {
                for (auto j = 0UL; j + l2 < nb_jobs2 + 1; ++j) {
                    auto t = 0;

                    if (i != 0) {
                        t += W[k1][i - 1];
                    }

                    t += B2_1[i][j] - B2_2[i][j];
                    t += B3_1[i][j];
                    t += B5_1[i][j] - B5_2[i][j];
                    t += B6_1[i][j];

                    if (j != 0) {
                        t += W[k2][j - 1];
                    }

                    if (sol.machines[k1].total_weighted_tardiness +
                            sol.machines[k2].total_weighted_tardiness - t >
                        max) {
                        max = sol.machines[k1].total_weighted_tardiness +
                              sol.machines[k2].total_weighted_tardiness - t;
                        i_best = i;
                        j_best = j;
                        k_best = k1;
                        kk_best = k2;
                        updated = true;
                    }
                }
            }
        }
    }

    /** update to best improvement */
    if (updated) {
        if (debug_lvl(0)) {
            sol.print_solution();
        }
        sol.update_swap_move_inter(i_best, j_best, k_best, kk_best, l1, l2);
        calculate_W(sol);
        calculate_g(sol);

        if (debug_lvl(0)) {
            sol.print_solution();
        }
    }

    if (debug_lvl(0)) {
        fmt::print(
            R"(inter insertion with l1 = {} and l2 = {} and improvement {} on machines {} and {} on places {} {}
)",
            l1, l2, max, k_best, kk_best, i_best, j_best);
    }
}

void LocalSearchData::insertion_operator_inter(Sol& sol, size_t l) {
    int         max{};
    size_t      i_best{};
    size_t      j_best{};
    size_t      k_best{};
    size_t      kk_best{};
    SlopeListIt it;
    SlopeListIt tmp_end;
    int         c{};

    ++iterations;
    create_processing_insertion_operator_inter(sol, l);
    updated = false;

    for (auto&& [k1, m1] : sol.machines | ranges::views::enumerate) {
        auto&       vec_m1 = m1.job_list;
        const auto& nb_jobs1 = vec_m1.size();

        for (auto&& [k2, m2] : sol.machines | ranges::views::enumerate) {
            if (k1 == k2) {
                continue;
            }

            auto&       vec_m2 = m2.job_list;
            const auto& nb_jobs2 = vec_m2.size();
            if (nb_jobs + 1 < l) {
                continue;
            }
            for (auto i = 0UL; i < nb_jobs1 - l + 1; ++i) {
                it = g[k1][i].begin();
                tmp_end = g[k1][i].end();

                for (auto j = 0UL; j < nb_jobs2; ++j) {
                    c = 0;

                    if (j != 0) {
                        c = sol.c[vec_m2[j - 1]->job];
                    }

                    it = std::find_if(it, tmp_end,
                                      [&c](const auto& tmp) -> bool {
                                          return tmp.in_interval(c);
                                      });
                    B2_1[i][j] = (it != tmp_end) ? (*it)(c) : 0;
                }
            }

            for (auto i = 0UL; i < nb_jobs1 - l + 1; ++i) {
                auto& p_list = processing_list[0][k1][i];

                if (p_list.pos + l >= nb_jobs1) {
                    for (auto j = 0UL; j < nb_jobs2; ++j) {
                        B2_2[p_list.pos][j] = 0;
                    }
                } else {
                    it = g[k1][p_list.pos + l].begin();
                    tmp_end = g[k1][p_list.pos + l].end();
                    for (auto j = 0UL; j < nb_jobs2; ++j) {
                        c = p_list.p;

                        if (j != 0UL) {
                            c += sol.c[vec_m2[j - 1]->job];
                        }

                        it = std::find_if(it, tmp_end,
                                          [&c](const auto& tmp) -> bool {
                                              return tmp.in_interval(c);
                                          });
                        B2_2[p_list.pos][j] = it != tmp_end ? (*it)(c) : 0;
                    }
                }
            }

            for (auto i = 0UL; i < nb_jobs1 - l + 1; ++i) {
                if (i + l >= nb_jobs1) {
                    B3_1_[i] = 0;
                } else {
                    c = 0;

                    it = g[k1][i].begin();
                    tmp_end = g[k1][i].end();

                    if (i != 0) {
                        c = sol.c[vec_m1[i - 1]->job];
                    }

                    it = std::find_if(it, tmp_end,
                                      [&c](const auto& tmp) -> bool {
                                          return tmp.in_interval(c);
                                      });
                    B3_1_[i] = it != tmp_end ? (*it)(c) : 0;
                }
            }

            for (auto j = 0UL; j + 1 < nb_jobs2; ++j) {
                it = g[k2][j].begin();
                tmp_end = g[k2][j].end();

                for (auto i = 0UL; i < nb_jobs1 - l + 1; ++i) {
                    auto& p_list = processing_list[0][k1][i];
                    c = p_list.p;

                    if (j != 0) {
                        c += sol.c[vec_m2[j - 1]->job];
                    }

                    it = std::find_if(it, tmp_end,
                                      [&c](const auto& tmp) -> bool {
                                          return tmp.in_interval(c);
                                      });
                    B5_1[p_list.pos][j] = (it != tmp_end) ? (*it)(c) : 0;
                }
            }

            for (auto i = 0UL; i + l < nb_jobs1 + 1; ++i) {
                for (auto j = 0UL; j + 1 < nb_jobs2; ++j) {
                    auto t = 0;

                    if (i != 0UL) {
                        t += W[k1][i - 1];
                    }

                    t += B2_1[i][j] - B2_2[i][j];
                    t += B3_1_[i];
                    t += B5_1[i][j];

                    if (j != 0UL) {
                        t += W[k2][j - 1];
                    }

                    if (m1.total_weighted_tardiness +
                            m2.total_weighted_tardiness - t >
                        max) {
                        max = m1.total_weighted_tardiness +
                              m2.total_weighted_tardiness - t;
                        i_best = i;
                        j_best = j;
                        k_best = k1;
                        kk_best = k2;
                        updated = true;
                    }
                }
            }
            k2++;
        }
        k1++;
    }

    if (updated) {
        if (debug_lvl(0)) {
            sol.print_solution();
        }
        sol.update_insertion_move_inter(i_best, j_best, k_best, kk_best, l);
        calculate_W(sol);
        calculate_g(sol);

        if (debug_lvl(0)) {
            sol.print_solution();
        }
    }

    if (debug_lvl(0)) {
        fmt::print(
            R"(inter insertion with l = {} and improvement {} on machines {} and {} on places {} {}
)",
            l, max, k_best, kk_best, i_best, j_best);
        fmt::print("{:-^30}", "");
    }
}

void LocalSearchData::calculate_W(const Sol& sol) {
    for (auto&& [i, m] : sol.machines | ranges::views::enumerate) {
        if (!m.updated) {
            continue;
        }
        auto tmp = m.job_list.begin();
        W[i][0] = (*tmp)->weighted_tardiness(sol.c[(*tmp)->job]);
        tmp++;
        auto j = 0UL;
        for (; tmp != m.job_list.end(); tmp++) {
            W[i][j + 1] =
                W[i][j] + (*tmp)->weighted_tardiness(sol.c[(*tmp)->job]);
            j++;
        }
    }
}

void LocalSearchData::calculate_g(Sol& sol) {
    for (auto&& [i, m] : sol.machines | ranges::views::enumerate) {
        if (m.updated) {
            for (auto j = 0UL; j < nb_jobs; j++) {
                g[i][j].clear();
            }
            m.updated = false;
        } else {
            continue;
        }

        int    P{0};
        size_t j{0UL};
        for (auto it = m.job_list.begin(); it != m.job_list.end(); it++) {
            VecJobRawPtr lateness_sort(it, m.job_list.end());

            ranges::sort(lateness_sort,
                         [&sol](const auto& lhs, const auto& rhs) -> bool {
                             return (sol.c[lhs->job] - lhs->due_time >
                                     sol.c[rhs->job] - rhs->due_time);
                         });

            int  tw{0};
            int  w{0};
            int  t1{0};
            bool move{true};
            auto k = 0U;

            for (; k < lateness_sort.size() && move;) {
                move = lateness_sort[k]->weight *
                           (sol.c[lateness_sort[k]->job] - P -
                            lateness_sort[k]->due_time) >
                       0;

                if (move) {
                    tw += lateness_sort[k]->weight *
                          (sol.c[lateness_sort[k]->job] - P -
                           lateness_sort[k]->due_time);
                    w += lateness_sort[k]->weight;
                    k++;
                }
            }

            int t2 =
                lateness_sort[k]->due_time - sol.c[lateness_sort[k]->job] + P;

            g[i][j].emplace_back(t1, t2, tw, w);

            for (auto l = k; l < lateness_sort.size();) {
                tw = tw + w * (t2 - t1);
                t1 = t2;
                move = true;

                while (move) {
                    w += lateness_sort[l]->weight;
                    l++;

                    if (l == lateness_sort.size()) {
                        move = false;
                        t2 = std::numeric_limits<int>::max();
                    } else {
                        t2 = lateness_sort[l]->due_time -
                             sol.c[lateness_sort[l]->job] + P;
                        move = (t1 == t2);
                    }
                }

                g[i][j].emplace_back(t1, t2, tw, w);
            }

            P += m.job_list[j]->processing_time;
            j++;
        }
    }
}

void LocalSearchData::create_processing_insertion_operator(const Sol& sol,
                                                           size_t     l) {
    for (auto&& [i, m] : sol.machines | ranges::views::enumerate) {
        auto        C = 0;
        const auto& machine = m.job_list;

        if (machine.size() < l) {
            continue;
        }

        for (auto j = 0UL; j < l; j++) {
            C += machine[j]->processing_time;
        }

        for (size_t j = 0; j < machine.size() - l; j++) {
            processing_list[0][i][j] = {j, C};
            C = C - machine[j]->processing_time +
                machine[j + l]->processing_time;
        }

        ranges::sort(
            processing_list[0][i] | ranges::views::take(machine.size() - l),
            std::greater{}, [](const auto& tmp) { return tmp.p; });
    }
}

void LocalSearchData::create_processing_list_swap(const Sol& sol,
                                                  size_t     l1,
                                                  size_t     l2) {
    for (auto&& [i, m] : sol.machines | ranges::views::enumerate) {
        const auto& vec_m = m.job_list;
        const auto  n = vec_m.size();

        if (n < l1 + l2) {
            continue;
        }

        auto C = ranges::accumulate(
            vec_m | ranges::views::take(l1 + l2) | ranges::views::drop(l1), 0,
            std::plus<>{},
            [](const auto& tmp) { return tmp->processing_time; });

        for (std::size_t j = l1; j < n - l2; ++j) {
            processing_list[1][i][j - l1] = {j, C};
            C = C - vec_m[j]->processing_time + vec_m[j + l2]->processing_time;
        }

        processing_list[1][i][n - l2 - l1] = {n - l2, C};
        ranges::sort(
            processing_list[1][i] | ranges::views::take(n - l1 - l2 + 1),
            std::less{}, [](const auto& tmp) { return tmp.p; });

        C = ranges::accumulate(
            vec_m | ranges::views::take(l1), 0, std::plus<>{},
            [](const auto& tmp) { return tmp->processing_time; });

        for (std::size_t j = 0; j < n - l1 - l2; ++j) {
            processing_list[0][i][j] = {j, C};
            C = C - vec_m[j]->processing_time + vec_m[j + l1]->processing_time;
        }

        processing_list[0][i][n - l1 - l2] = {n - l1 - l2, C};
        ranges::sort(
            processing_list[0][i] | ranges::views::take(n - l1 - l2 + 1),
            std::less{}, [](const auto& tmp) { return tmp.p; });
    }
}

void LocalSearchData::create_processing_insertion_operator_inter(const Sol& sol,
                                                                 size_t     l) {
    for (auto&& [i, m] : sol.machines | ranges::views::enumerate) {
        const auto& nb = m.job_list.size();
        const auto& vec_m = m.job_list;

        if (nb < l) {
            continue;
        }

        auto C = ranges::accumulate(
            vec_m | ranges::views::take(l), 0, std::plus<>{},
            [](const auto& tmp) { return tmp->processing_time; });

        for (auto j = 0UL; j < nb - l; ++j) {
            auto* j1 = vec_m[j];
            auto* j2 = vec_m[j + l];
            processing_list[0][i][j] = {j, C};
            C = C - j1->processing_time + j2->processing_time;
        }

        processing_list[0][i][nb - l] = {nb - l, C};
        ranges::sort(processing_list[0][i] | ranges::views::take(nb - l + 1),
                     std::less{}, [](const auto& tmp) { return tmp.p; });
    }
}

void LocalSearchData::create_processing_list_swap_inter(const Sol& sol,
                                                        size_t     l1,
                                                        size_t     l2) {
    for (auto i = 0UL; i < nb_machines; ++i) {
        const auto& machine = sol.machines[i].job_list;
        const auto  nb = machine.size();

        if (nb < l1 || nb < l2) {
            continue;
        }

        auto C = ranges::accumulate(
            machine | ranges::views::take(l1), 0, std::plus<>{},
            [](const auto& tmp) { return tmp->processing_time; });

        for (auto j = 0UL; j < nb - l1; ++j) {
            auto* j1 = machine[j];
            auto* j2 = machine[j + l1];
            processing_list[0][i][j] = {j, C};
            C = C - j1->processing_time + j2->processing_time;
        }

        processing_list[0][i][nb - l1] = {nb - l1, C};

        ranges::sort(processing_list[0][i] | ranges::views::take(nb - l1 + 1),
                     std::less{}, [](const auto& tmp) -> int { return tmp.p; });

        C = ranges::accumulate(
            machine | ranges::views::take(l2), 0, std::plus<>(),
            [](const auto& tmp) -> int { return tmp->processing_time; });

        for (auto j = 0UL; j < nb - l2; ++j) {
            auto* j1 = machine[j];
            auto* j2 = machine[j + l2];
            processing_list[1][i][j] = {j, C};
            C = C - j1->processing_time + j2->processing_time;
        }

        processing_list[1][i][nb - l2] = {nb - l2, C};
        ranges::sort(processing_list[1][i] | ranges::views::take(nb - l2 + 1),
                     std::less{}, [](const auto& tmp) -> int { return tmp.p; });
    }
}

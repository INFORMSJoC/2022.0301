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

#ifndef STATISTICS_H
#define STATISTICS_H

#include <fmt/format.h>               // for format
#include <boost/chrono/duration.hpp>  // for den
#include <boost/timer/timer.hpp>      // for cpu_timer
#include <cstddef>                    // for size_t
#include <span>                       // for span
#include <string>                     // for string
#include "or-utils/util.h"            // for CCutil_timer
struct Parms;

class Timer : public boost::timer::cpu_timer {
    enum ClockType {
        wall_time,
        cpu_time,
    };

    std::string _name;
    ClockType   _type;

   public:
    Timer(std::string name_ = "timer", ClockType type = ClockType::wall_time)
        : _name(std::move(name_)),
          _type(type) {
        this->stop();
        this->elapsed().clear();
    };

    Timer(std::string name_ = "timer", bool type = false)
        : _name(std::move(name_)),
          _type(type ? cpu_time : wall_time) {
        this->stop();
        this->elapsed().clear();
    };

    [[nodiscard]] auto name() const -> const std::string& { return _name; };
    [[nodiscard]] auto type() const -> ClockType { return _type; }

    [[nodiscard]] auto dbl_sec() const -> double {
        auto aux = boost::chrono::nanoseconds();
        switch (_type) {
            case wall_time:

                aux = boost::chrono::nanoseconds(this->elapsed().wall);

                return static_cast<double>(aux.count()) *
                       boost::chrono::nanoseconds::period::num /
                       boost::chrono::nanoseconds::period::den;

            case cpu_time:
                aux = boost::chrono::nanoseconds(this->elapsed().user +
                                                 this->elapsed().system);

                return static_cast<double>(aux.count()) *
                       boost::chrono::nanoseconds::period::num /
                       boost::chrono::nanoseconds::period::den;
        }
        return 0.0;
    }

    [[nodiscard]] auto nano_sec() const -> boost::timer::nanosecond_type {
        switch (_type) {
            case wall_time:
                return this->elapsed().wall;

            case cpu_time:
                auto cpu_aux = boost::chrono::nanoseconds(
                    this->elapsed().user + this->elapsed().system);
                return cpu_aux.count();
        }
        return boost::timer::nanosecond_type{};
    }

    [[nodiscard]] auto str_sec(
        short precision = boost::timer::default_places) const -> std::string {
        switch (_type) {
            case cpu_time:
                return this->format(precision, "%t");
            case wall_time:
            default:
                return this->format(precision, "%w");
        }
    }
};

struct Statistics {
    int    global_upper_bound;
    int    global_lower_bound;
    double rel_error;

    int    root_upper_bound;
    int    root_lower_bound;
    double root_rel_error;

    size_t nb_generated_col;
    size_t nb_generated_col_root;
    size_t first_size_graph;
    size_t size_graph_after_reduced_cost_fixing;

    enum TimerType {
        build_dd_timer,
        cputime_timer,
        bb_timer,
        lb_root_timer,
        lb_timer,
        solve_lp_timer,
        pricing_timer,
        heuristic_timer,
        reduced_cost_fixing_timer
    };

    Timer time_build_dd;
    Timer time_total;
    Timer time_branch_and_bound;
    Timer time_strong_branching;
    Timer time_lb_root;
    Timer time_lb;
    Timer time_solve_lp;
    Timer time_pricing;
    Timer time_heuristic;
    Timer time_rc_fixing;

    double real_time_total;

    int    mip_nb_vars;
    int    mip_nb_constr;
    double mip_obj_bound;
    double mip_obj_bound_lp;
    double mip_rel_gap;
    double mip_run_time;
    int    mip_status;
    double mip_nb_iter_simplex;
    double mip_nb_nodes;
    int    mip_reduced_cost_fixing;

    std::string pname;

    void start_resume_timer(TimerType _type);
    void suspend_timer(TimerType _type);

    [[nodiscard]] auto total_time_str(TimerType _type, short precision) const
        -> std::string;

    [[nodiscard]] auto total_time_dbl(TimerType _type) const -> double;
    [[nodiscard]] auto total_time_nano_sec(TimerType _type) const
        -> boost::timer::nanosecond_type;

    Statistics(const Parms& parms);

    static constexpr auto DEFAULT_MIP_RUN = 110.0;
};

template <>
struct fmt::formatter<Statistics> {
    char presentation = 'v';
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx) -> decltype(ctx.begin()) {
        auto aux = std::span<const char>{ctx.begin(), ctx.end()};

        if (aux.empty())
        {
            return nullptr;
        }

        auto it = aux.begin();
        auto end = aux.end();

        if (it != end && (*it == 'v' || *it == 'n')) {
            presentation = *it++;
        }

        if (it != end && *it != '}') {
            throw format_error("invalid format");
        }

        return std::addressof(*it);
    }

    template <typename FormatContext>
    auto format(const Statistics& stat, FormatContext& ctx)
        -> decltype(ctx.out()) {
        if (presentation == 'n') {
            return format_to(
                ctx.out(),
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}"
                ","
                "{},{},{},{},{},{},{}",
                "tot_real_time", stat.time_total.name(),
                stat.time_branch_and_bound.name(), stat.time_lb.name(),
                stat.time_lb_root.name(), stat.time_heuristic.name(),
                stat.time_build_dd.name(), stat.time_pricing.name(),
                stat.time_rc_fixing.name(), "rel_error", "global_lower_bound",
                "global_upper_bound", "first_rel_error",
                "global_upper_bound_root", "global_lowerbound_root",
                "nb_generated_col", "nb_generated_col_root", "first_size_graph",
                "size_after_reduced_cost", "mip_nb_vars", "mip_nb_constr",
                "mip_obj_bound", "mip_obj_bound_lp", "mip_rel_gap",
                "mip_run_time", "mip_status", "mip_nb_iter_simplex",
                "mip_nb_nodes");
        }

        return format_to(
            ctx.out(),
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}"
            ",{},{},{},{},{},{}",
            stat.real_time_total,
            stat.total_time_str(Statistics::cputime_timer, PRECISION),
            stat.total_time_str(Statistics::bb_timer, PRECISION),
            stat.total_time_str(Statistics::lb_timer, PRECISION),
            stat.total_time_str(Statistics::lb_root_timer, PRECISION),
            stat.total_time_str(Statistics::heuristic_timer, PRECISION),
            stat.total_time_str(Statistics::build_dd_timer, PRECISION),
            stat.total_time_str(Statistics::pricing_timer, PRECISION),
            stat.total_time_str(Statistics::reduced_cost_fixing_timer,
                                PRECISION),
            stat.rel_error, stat.global_lower_bound, stat.global_upper_bound,
            stat.root_rel_error, stat.root_upper_bound, stat.root_lower_bound,
            stat.nb_generated_col, stat.nb_generated_col_root,
            stat.first_size_graph, stat.size_graph_after_reduced_cost_fixing,
            stat.mip_nb_vars, stat.mip_nb_constr, stat.mip_obj_bound,
            stat.mip_obj_bound_lp, stat.mip_rel_gap, stat.mip_run_time,
            stat.mip_status, stat.mip_nb_iter_simplex, stat.mip_nb_nodes);
    }

    static constexpr auto PRECISION = 5;
};

#endif  // STATISTICS_H

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

#include "Statistics.h"
#include <limits>
#include "Parms.h"  // for Parms

Statistics::Statistics(const Parms& _parms)
    : global_upper_bound(std::numeric_limits<int>::max()),
      global_lower_bound(0),
      rel_error(std::numeric_limits<double>::max()),
      root_upper_bound(std::numeric_limits<int>::max()),
      root_lower_bound(0),
      root_rel_error(std::numeric_limits<double>::max()),
      nb_generated_col(0),
      nb_generated_col_root(0),
      first_size_graph(0),
      size_graph_after_reduced_cost_fixing(0),
      time_build_dd{"tot_build_dd", _parms.use_cpu_time.value()},
      time_total{"tot_cputime", _parms.use_cpu_time.value()},
      time_branch_and_bound{"tot_bb", _parms.use_cpu_time.value()},
      time_strong_branching{"tot_strong_branching",
                            _parms.use_cpu_time.value()},
      time_lb_root{"tot_lb_root", _parms.use_cpu_time.value()},
      time_lb{"tot_lb", _parms.use_cpu_time.value()},
      time_solve_lp{"tot_solve_lp", _parms.use_cpu_time.value()},
      time_pricing{"tot_pricing", _parms.use_cpu_time.value()},
      time_heuristic{"tot_heuristic", _parms.use_cpu_time.value()},
      time_rc_fixing{"total_reduce_cost_fixing", _parms.use_cpu_time.value()},
      real_time_total(getRealTime()),
      mip_nb_vars(0),
      mip_nb_constr(0),
      mip_obj_bound(0.0),
      mip_obj_bound_lp(0.0),
      mip_rel_gap(0.0),
      mip_run_time(DEFAULT_MIP_RUN),
      mip_status(0),
      mip_nb_iter_simplex(),
      mip_nb_nodes(0),
      mip_reduced_cost_fixing(),
      pname(_parms.pname) {
    start_resume_timer(cputime_timer);
}

void Statistics::start_resume_timer(TimerType _type) {
    switch (_type) {
        case build_dd_timer:
            time_build_dd.resume();
            break;
        case cputime_timer:
            time_total.resume();
            break;
        case bb_timer:
            time_branch_and_bound.resume();
            break;
        case lb_root_timer:
            time_lb_root.resume();
            break;
        case lb_timer:
            time_lb.resume();
            break;
        case solve_lp_timer:
            time_solve_lp.resume();
            break;
        case pricing_timer:
            time_pricing.resume();
            break;
        case heuristic_timer:
            time_heuristic.resume();
            break;
        case reduced_cost_fixing_timer:
            time_rc_fixing.resume();
            break;
    }
}

void Statistics::suspend_timer(TimerType _type) {
    switch (_type) {
        case build_dd_timer:
            time_build_dd.stop();
            break;
        case cputime_timer:
            time_total.stop();
            break;
        case bb_timer:
            time_branch_and_bound.stop();
            break;
        case lb_root_timer:
            time_lb_root.stop();
            break;
        case lb_timer:
            time_lb.stop();
            break;
        case solve_lp_timer:
            time_solve_lp.stop();
            break;
        case pricing_timer:
            time_pricing.stop();
            break;
        case heuristic_timer:
            time_heuristic.stop();
            break;
        case reduced_cost_fixing_timer:
            time_rc_fixing.stop();
            break;
    }
}

auto Statistics::total_time_dbl(TimerType _type) const -> double {
    switch (_type) {
        case build_dd_timer:
            return time_build_dd.dbl_sec();
        case cputime_timer:
            return time_total.dbl_sec();
        case bb_timer:
            return time_branch_and_bound.dbl_sec();
        case lb_root_timer:
            return time_lb_root.dbl_sec();
        case lb_timer:
            return time_lb.dbl_sec();
        case solve_lp_timer:
            return time_solve_lp.dbl_sec();
        case pricing_timer:
            return time_pricing.dbl_sec();
        case heuristic_timer:
            return time_heuristic.dbl_sec();
        case reduced_cost_fixing_timer:
            return time_rc_fixing.dbl_sec();
    }
    return 0.0;
}

auto Statistics::total_time_nano_sec(TimerType _type) const
    -> boost::timer::nanosecond_type {
    switch (_type) {
        case build_dd_timer:
            return time_build_dd.nano_sec();
        case cputime_timer:
            return time_total.nano_sec();
        case bb_timer:
            return time_branch_and_bound.nano_sec();
        case lb_root_timer:
            return time_lb_root.nano_sec();
        case lb_timer:
            return time_lb.nano_sec();
        case solve_lp_timer:
            return time_solve_lp.nano_sec();
        case pricing_timer:
            return time_pricing.nano_sec();
        case heuristic_timer:
            return time_heuristic.nano_sec();
        case reduced_cost_fixing_timer:
            return time_rc_fixing.nano_sec();
    }
    return boost::timer::nanosecond_type{};
}

auto Statistics::total_time_str(TimerType _type, short precision) const
    -> std::string {
    switch (_type) {
        case build_dd_timer:
            return time_build_dd.str_sec(precision);
        case cputime_timer:
            return time_total.str_sec(precision);
        case bb_timer:
            return time_branch_and_bound.str_sec(precision);
        case lb_root_timer:
            return time_lb_root.str_sec(precision);
        case lb_timer:
            return time_lb.str_sec(precision);
        case solve_lp_timer:
            return time_solve_lp.str_sec(precision);
        case pricing_timer:
            return time_pricing.str_sec(precision);
        case heuristic_timer:
            return time_heuristic.str_sec(precision);
        case reduced_cost_fixing_timer:
            return time_rc_fixing.str_sec(precision);
    }
    return "";
}

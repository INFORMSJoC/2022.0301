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

#ifndef WCT_PRIVATE_H
#define WCT_PRIVATE_H
#include <array>         // for array
#include <cstddef>       // for size_t
#include <exception>     // for exception
#include <functional>    // for function
#include <memory>        // for unique_ptr, shared_ptr
#include <string>        // for string
#include <vector>        // for vector
#include "Instance.h"    // for Instance
#include "Parms.h"       // for Parms
#include "Solution.hpp"  // for Sol
#include "Statistics.h"  // for Statistics
class BranchBoundTree;   // lines 21-21
class PricingStabilizationBase;
struct NodeData;  // lines 19-19
struct PricerSolverBase;
struct Column;  // lines 18-18

/**
 * wct data types nodes of branch and bound tree
 */
/**
 *  CONSTANTS NODEDATA STRUCTURE
 *
 */

enum problem_status {
    no_sol = 0,
    lp_feasible = 1,
    feasible = 2,
    meta_heuristic = 3,
    optimal = 4
};

/**
 * problem data
 */
class Problem {
   private:
    /** Different Parameters */
    Parms parms;
    /*Cpu time measurement + Statistics*/
    Statistics stat;
    /** Instance data*/
    Instance instance;

    std::unique_ptr<BranchBoundTree> tree{};
    std::unique_ptr<NodeData>        root_pd;

    problem_status status;

    /* Best Solution*/
    Sol opt_sol{};

    static constexpr auto EPS = 1e-6;

   public:
    /** All methods of problem class */
    auto to_screen() -> int;
    void to_csv();
    void solve();
    /** Heuristic related */
    void heuristic();
    /** Constructors */
    Problem(int argc, const char** argv);
    Problem(const Problem&) = delete;
    Problem(Problem&&) = delete;
    ~Problem();

    auto operator=(const Problem&) -> Problem& = delete;
    auto operator=(Problem&&) -> Problem& = delete;

    friend NodeData;

    class ProblemException : public std::exception {
       public:
        ProblemException(const char* const msg = nullptr) : errmsg(msg) {}

        [[nodiscard]] auto what() const noexcept -> const char* override {
            return (errmsg);
        }

       private:
        const char* errmsg;
    };
};

#endif

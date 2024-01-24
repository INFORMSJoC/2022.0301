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

#ifndef PRICER_SOLVER_ZDD_FORWARD_HPP
#define PRICER_SOLVER_ZDD_FORWARD_HPP
#include "Instance.h"             // for Instance
#include "PricerEvaluateZdd.hpp"  // for BackwardZddCycleDouble, BackwardZdd...
#include "PricerSolverZdd.hpp"    // for PricerSolverZdd
#include "PricingSolution.hpp"    // for PricingSolution
class PricerSolverSimple : public PricerSolverZdd {
   private:
    ForwardZddSimpleDouble  evaluator;
    BackwardZddSimpleDouble reversed_evaluator;

   public:
    explicit PricerSolverSimple(const Instance& instance);
    PricerSolverSimple(PricerSolverSimple&&) = default;
    PricerSolverSimple(const PricerSolverSimple&) = default;
    ~PricerSolverSimple() override = default;

    auto operator=(const PricerSolverSimple&) -> PricerSolverSimple& = delete;
    auto operator=(PricerSolverSimple&&) -> PricerSolverSimple& = delete;

    auto pricing_algorithm(std::span<const double>& _pi)
        -> PricingSolution override;
    auto evaluate_nodes(std::span<const double>& pi) -> bool final;
    void compute_labels(std::span<const double>& _pi);
    void compute_labels(double* _pi);
};

class PricerSolverZddCycle : public PricerSolverZdd {
   private:
    ForwardZddCycleDouble  evaluator;
    BackwardZddCycleDouble reversed_evaluator;

   public:
    explicit PricerSolverZddCycle(const Instance& instance);
    auto pricing_algorithm(std::span<const double>& _pi)
        -> PricingSolution override;
    PricerSolverZddCycle(const PricerSolverZddCycle&) = default;
    PricerSolverZddCycle(PricerSolverZddCycle&&) = default;
    ~PricerSolverZddCycle() override = default;

    auto operator=(PricerSolverZddCycle&&) -> PricerSolverZddCycle& = default;
    auto operator=(const PricerSolverZddCycle&)
        -> PricerSolverZddCycle& = delete;
    auto evaluate_nodes(std::span<const double>& pi) -> bool final;
    void compute_labels(double* _pi);
    void compute_labels(std::span<const double>& _pi);
};

#endif  // PRICER_SOLVER_ZDD_FORWARD_HPP

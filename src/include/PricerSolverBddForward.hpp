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

#ifndef PRICER_SOLVER_BDD_FORWARD_HPP
#define PRICER_SOLVER_BDD_FORWARD_HPP
#include <memory>                 // for make_unique, unique_ptr
#include "Instance.h"             // for Instance
#include "NodeBdd.hpp"            // for NodeBdd
#include "PricerEvaluateBdd.hpp"  // for BackwardBddFarkasDouble, BackwardBd...
#include "PricerSolverBdd.hpp"    // for PricerSolverBdd
#include "PricingSolution.hpp"    // for PricingSolution
struct PricerSolverBase;

class PricerSolverBddSimple : public PricerSolverBdd {
   private:
    ForwardBddSimpleDouble  evaluator;
    BackwardBddSimpleDouble reversed_evaluator;
    BackwardBddFarkasDouble farkas_evaluator;

   public:
    explicit PricerSolverBddSimple(const Instance& instance);

    PricerSolverBddSimple(const PricerSolverBddSimple&) = default;
    PricerSolverBddSimple(PricerSolverBddSimple&&) = default;
    ~PricerSolverBddSimple() override = default;

    auto operator=(PricerSolverBddSimple&&) -> PricerSolverBddSimple& = default;
    auto operator=(const PricerSolverBddSimple&)
        -> PricerSolverBddSimple& = delete;

    [[nodiscard]] auto clone() const
        -> std::unique_ptr<PricerSolverBase> override {
        return std::make_unique<PricerSolverBddSimple>(*this);
    };

    auto pricing_algorithm(std::span<const double>& _pi)
        -> PricingSolution override;
    auto farkas_pricing(std::span<const double>& _pi)
        -> PricingSolution override;
    auto evaluate_rc_arc(NodeBdd& n) -> double override;

    void compute_labels(std::span<const double>& _pi) override;
};

class PricerSolverBddCycle : public PricerSolverBdd {
   private:
    ForwardBddCycleDouble   evaluator;
    BackwardBddCycleDouble  reversed_evaluator;
    BackwardBddFarkasDouble farkas_evaluator;

   public:
    explicit PricerSolverBddCycle(const Instance& instance);

    PricerSolverBddCycle(const PricerSolverBddCycle&) = default;
    PricerSolverBddCycle(PricerSolverBddCycle&&) = default;
    ~PricerSolverBddCycle() override = default;

    auto operator=(PricerSolverBddCycle&&) -> PricerSolverBddCycle& = default;
    auto operator=(const PricerSolverBddCycle&)
        -> PricerSolverBddCycle& = delete;

    [[nodiscard]] auto clone() const
        -> std::unique_ptr<PricerSolverBase> override {
        return std::make_unique<PricerSolverBddCycle>(*this);
    }

    auto pricing_algorithm(std::span<const double>& _pi)
        -> PricingSolution override;
    auto farkas_pricing(std::span<const double>& _pi)
        -> PricingSolution override;
    auto evaluate_rc_arc(NodeBdd& n) -> double override;

    void compute_labels(std::span<const double>& _pi) override;
};

#endif  // PRICER_SOLVER_BDD_FORWARD_HPP

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

#include <memory>                 // for make_unique, unique_ptr
#include "Instance.h"             // for Instance
#include "NodeBdd.hpp"            // for NodeBdd
#include "PricerEvaluateBdd.hpp"  // for BackwardBddFarkasDouble, BackwardBd...
#include "PricerSolverBdd.hpp"    // for PricerSolverBdd
#include "PricingSolution.hpp"    // for PricingSolution
struct PricerSolverBase;

class PricerSolverBddBackwardSimple : public PricerSolverBdd {
   private:
    BackwardBddSimpleDouble evaluator;
    ForwardBddSimpleDouble  reversed_evaluator;
    BackwardBddFarkasDouble farkas_evaluator;

   public:
    explicit PricerSolverBddBackwardSimple(const Instance& instance);

    PricerSolverBddBackwardSimple(const PricerSolverBddBackwardSimple&) =
        default;
    PricerSolverBddBackwardSimple(PricerSolverBddBackwardSimple&&) = default;
    ~PricerSolverBddBackwardSimple() override = default;

    auto operator=(PricerSolverBddBackwardSimple&&)
        -> PricerSolverBddBackwardSimple& = default;
    auto operator=(const PricerSolverBddBackwardSimple&)
        -> PricerSolverBddBackwardSimple& = delete;

    [[nodiscard]] auto clone() const
        -> std::unique_ptr<PricerSolverBase> override {
        return std::make_unique<PricerSolverBddBackwardSimple>(*this);
    };

    auto pricing_algorithm(std::span<const double>& _pi)
        -> PricingSolution override;
    auto farkas_pricing(std::span<const double>& _pi)
        -> PricingSolution override;
    auto evaluate_rc_arc(NodeBdd& n) -> double override;

    void compute_labels(std::span<const double>& _pi) override;
};

class PricerSolverBddBackwardCycle : public PricerSolverBdd {
   private:
    BackwardBddCycleDouble  evaluator;
    ForwardBddCycleDouble   reversed_evaluator;
    BackwardBddFarkasDouble farkas_evaluator;

   public:
    explicit PricerSolverBddBackwardCycle(const Instance& instance);
    PricerSolverBddBackwardCycle(PricerSolverBddBackwardCycle&&) = default;
    PricerSolverBddBackwardCycle(const PricerSolverBddBackwardCycle&) = default;
    ~PricerSolverBddBackwardCycle() override = default;

    auto operator=(const PricerSolverBddBackwardCycle&)
        -> PricerSolverBddBackwardCycle& = delete;
    auto operator=(PricerSolverBddBackwardCycle&&)
        -> PricerSolverBddBackwardCycle& = default;

    [[nodiscard]] auto clone() const
        -> std::unique_ptr<PricerSolverBase> override {
        return std::make_unique<PricerSolverBddBackwardCycle>(*this);
    };

    auto pricing_algorithm(std::span<const double>& _pi)
        -> PricingSolution override;
    auto farkas_pricing(std::span<const double>& _pi)
        -> PricingSolution override;
    auto evaluate_rc_arc(NodeBdd& n) -> double override;

    void compute_labels(std::span<const double>& _pi) override;
};
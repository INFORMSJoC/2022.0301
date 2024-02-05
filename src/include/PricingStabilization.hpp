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

#ifndef PRICINGSTABILIZATION_H
#define PRICINGSTABILIZATION_H

#include <cstddef>              // for size_t
#include <memory>               // for unique_ptr
#include <vector>               // for vector
#include "PricingSolution.hpp"  // for PricingSolution
struct PricerSolverBase;        // lines 9-9
class PricingStabilizationBase {
   public:
    explicit PricingStabilizationBase(PricerSolverBase*        _solver,
                                      std::span<const double>& _pi_out);
    PricingStabilizationBase(PricingStabilizationBase&&) = default;
    PricingStabilizationBase(const PricingStabilizationBase&) = delete;
    virtual ~PricingStabilizationBase() = default;

    auto operator=(PricingStabilizationBase&&)
        -> PricingStabilizationBase& = delete;
    auto operator=(const PricingStabilizationBase&)
        -> PricingStabilizationBase& = delete;

    virtual auto clone(PricerSolverBase* _solver, std::span<const double>& _pi)
        -> std::unique_ptr<PricingStabilizationBase>;

    static constexpr double EPS_RC = -1e-10;
    static constexpr double ETA_DIFF = 1e-4;
    static constexpr double EPS_STAB = 1e-9;
    static constexpr int    RC_FIXING_RATE = 50;

    PricerSolverBase* solver;
    PricingSolution   sol;

    double reduced_cost{};
    double eta_in{};
    double eta_sep{};
    double eta_out{};

    std::span<const double>& pi_out;
    std::vector<double>      pi_in;
    std::vector<double>      pi_sep;

    bool   update_stab_center{};
    size_t iterations{};
    size_t previous_fix{};

    bool continueLP{};

    auto               get_sol() -> PricingSolution&;
    [[nodiscard]] auto get_update_stab_center() const -> bool;
    [[nodiscard]] auto get_reduced_cost() const -> double;

    static auto get_eps_stab_solver() -> double;

    virtual auto get_eta_in() -> double;
    virtual auto get_eta_sep() -> double;
    virtual auto stopping_criteria() -> bool;
    virtual auto reduced_cost_fixing() -> bool;
    virtual auto do_reduced_cost_fixing() -> bool;

    virtual void solve(double eta_out, double* _lhs);
    virtual void update_duals();
    virtual void remove_constraints(int first, int nb_del);
    virtual void update_continueLP(double obj);
    virtual void set_alpha([[maybe_unused]] double _alpha){};
};

class PricingStabilizationStat : public PricingStabilizationBase {
   public:
    explicit PricingStabilizationStat(PricerSolverBase*        _solver,
                                      std::span<const double>& _pi_out);
    PricingStabilizationStat(PricingStabilizationStat&&) = default;
    PricingStabilizationStat(const PricingStabilizationStat&) = delete;
    ~PricingStabilizationStat() override = default;

    auto operator=(PricingStabilizationStat&&)
        -> PricingStabilizationStat& = delete;
    auto operator=(const PricingStabilizationStat&)
        -> PricingStabilizationStat& = delete;

    auto clone(PricerSolverBase* _solver, std::span<const double>& _pi)
        -> std::unique_ptr<PricingStabilizationBase> override;

    static constexpr double ALPHA = 0.8;
    static constexpr double ALPHA_CHG = 0.1;
    static constexpr double ALPHA_MAX = 0.99;

    double alpha{ALPHA};
    double alphabar{};
    int    update{};
    int    k{};
    bool   hasstabcenter{};
    void   compute_pi_eta_sep(double _alpha_bar);

    auto stopping_criteria() -> bool final;
    auto get_eta_in() -> double final;

    void solve(double _eta_out, double* _lhs_coeff) override;
    void update_duals() override;
    void remove_constraints(int first, int nb_del) override;
    void set_alpha(double _alpha) override;
};

class PricingStabilizationDynamic : public PricingStabilizationStat {
   public:
    explicit PricingStabilizationDynamic(PricerSolverBase*        _solver,
                                         std::span<const double>& _pi_out);
    PricingStabilizationDynamic(PricingStabilizationDynamic&&) = default;
    PricingStabilizationDynamic(const PricingStabilizationDynamic&) = delete;
    ~PricingStabilizationDynamic() override = default;

    auto operator=(PricingStabilizationDynamic&&)
        -> PricingStabilizationDynamic& = delete;
    auto operator=(const PricingStabilizationDynamic&)
        -> PricingStabilizationDynamic& = delete;

    auto clone(PricerSolverBase* _solver, std::span<const double>& _pi)
        -> std::unique_ptr<PricingStabilizationBase> override;

    std::vector<double> subgradient;

    double subgradientnorm{};
    double subgradientproduct{};
    double dualdiffnorm{};

    void solve(double _eta_out, double* _lhs_coeff) override;
    void update_duals() override;
    void remove_constraints(int first, int nb_del) override;

    void compute_subgradient(const PricingSolution& _sol);
    void adjust_alpha();
};

class PricingStabilizationHybrid : public PricingStabilizationDynamic {
   public:
    explicit PricingStabilizationHybrid(PricerSolverBase*        pricer_solver,
                                        std::span<const double>& _pi_out);

    PricingStabilizationHybrid(PricingStabilizationHybrid&&) = default;
    PricingStabilizationHybrid(const PricingStabilizationHybrid&) = delete;
    ~PricingStabilizationHybrid() override = default;

    auto operator=(PricingStabilizationHybrid&&)
        -> PricingStabilizationHybrid& = delete;
    auto operator=(const PricingStabilizationHybrid&)
        -> PricingStabilizationHybrid& = delete;

    auto clone(PricerSolverBase* _solver, std::span<const double>& _pi)
        -> std::unique_ptr<PricingStabilizationBase> override;

    std::vector<double> subgradient_in;

    double beta{};
    double hybridfactor{};
    bool   in_mispricing_schedule{};

    void update_duals() override;
    void remove_constraints(int first, int nb_del) override;

   private:
    void solve(double _eta_out, double* _lhs) override;

    void update_hybrid() {
        if (hasstabcenter && !in_mispricing_schedule && alpha > 0.0) {
            calculate_dualdiffnorm();
            calculate_beta();
            calculate_hybridfactor();
        }
    }

    auto is_stabilized() -> bool;
    auto compute_dual(auto i) -> double;

    void calculate_dualdiffnorm();
    void calculate_beta();
    void calculate_hybridfactor();
    void update_stabcenter(const PricingSolution& _sol);
    void compute_subgradient_norm(const PricingSolution& _sol);
    void update_subgradientproduct();
    void update_alpha_misprice();
    void update_alpha();
};

#endif  // PRICINGSTABILIZATION_H
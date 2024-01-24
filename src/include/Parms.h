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

#ifndef PARMS_H
#define PARMS_H

#include <fmt/format.h>
#include <array>                  // for array
#include <boost/timer/timer.hpp>  // for boost::timer::nanosecond_type
#include <cstddef>                // for size_t
#include <functional>             // for function
#include <nlohmann/json.hpp>      // for json
#include <span>
#include <string>  // for allocator, string

enum BBNodeSelection {
    min_search_strategy = 0,
    no_branching = min_search_strategy,
    min_lb_strategy = 1,
    dfs_strategy = 2,
    max_strategy = 3,
};

enum PricingSolver {
    bdd_solver_simple = 0,
    bdd_solver_cycle = 1,
    bdd_solver_backward_simple = 2,
    bdd_solver_backward_cycle = 3,
    dp_solver = 4,
    ati_solver = 5,
    dp_bdd_solver = 6
};

enum StabTechniques {
    no_stab = 0,
    stab_wentgnes = 1,
    stab_dynamic = 2,
    stab_hybrid = 3,
    min_stab = stab_wentgnes,
};

enum BBExploreStrategy {
    min_bb_explore_strategy = 0,
    bb_dfs_strategy = min_bb_explore_strategy,
    bb_bfs_strategy = 1,
    bb_brfs_strategy = 2,
    bb_cbfs_strategy = 3,
};

enum Scoring_Parameter {
    min_scoring_parameter = 0,
    product_scoring_parameter = min_scoring_parameter,
    min_function_scoring_parameter = 1,
    weighted_sum_scoring_parameter = 2,
    weighted_product_scoring_parameter = 3,
    max_function_scoring_parameter = 4,
};

enum Scoring_Value {
    min_scoring_value = 0,
    lb_scoring_value = min_scoring_value,
    size_scoring_value = 1,
    nb_paths_scoring_value = 2,
};

#define NLOHMANN_JSON_SERIALIZE_ENUM_ARRAY(ENUM_TYPE, NB_ENUM_TYPE, ...)       \
    template <typename BasicJsonType>                                          \
    inline void to_json(BasicJsonType& j, const ENUM_TYPE& e) {                \
        static_assert(std::is_enum<ENUM_TYPE>::value,                          \
                      #ENUM_TYPE " must be an enum!");                         \
        static const std::array<std::pair<ENUM_TYPE, BasicJsonType>,           \
                                NB_ENUM_TYPE>                                  \
             m = {__VA_ARGS__};                                                \
        auto it = std::find_if(                                                \
            std::begin(m), std::end(m),                                        \
            [e](const std::pair<ENUM_TYPE, BasicJsonType>& ej_pair) -> bool {  \
                return ej_pair.first == e;                                     \
            });                                                                \
        j = ((it != std::end(m)) ? it : std::begin(m))->second;                \
    }                                                                          \
    template <typename BasicJsonType>                                          \
    inline void from_json(const BasicJsonType& j, ENUM_TYPE& e) {              \
        static_assert(std::is_enum<ENUM_TYPE>::value,                          \
                      #ENUM_TYPE " must be an enum!");                         \
        static const std::array<std::pair<ENUM_TYPE, BasicJsonType>,           \
                                NB_ENUM_TYPE>                                  \
             m = {__VA_ARGS__};                                                \
        auto it = std::find_if(                                                \
            std::begin(m), std::end(m),                                        \
            [&j](const std::pair<ENUM_TYPE, BasicJsonType>& ej_pair) -> bool { \
                return ej_pair.second == j;                                    \
            });                                                                \
        e = ((it != std::end(m)) ? it : std::begin(m))->first;                 \
    }

template <>
struct fmt::formatter<PricingSolver> : formatter<string_view> {
    // parse is inherited from formatter<string_view>.
    template <typename FormatContext>
    auto format(PricingSolver _solver, FormatContext& ctx) {
        string_view name = "unknown";
        switch (_solver) {
            case PricingSolver::bdd_solver_simple:
                name = "BddForward";
                break;
            case PricingSolver::bdd_solver_cycle:
                name = "BddForwardCycle";
                break;
            case PricingSolver::bdd_solver_backward_simple:
                name = "BddBackward";
                break;
            case PricingSolver::bdd_solver_backward_cycle:
                name = "BddBackwardCycle";
                break;
            case PricingSolver::dp_solver:
                name = "TimeIndexed";
                break;
            case PricingSolver::ati_solver:
                name = "ArcTimeIndexed";
                break;
            case PricingSolver::dp_bdd_solver:
                name = "Hybrid";
                break;
        }
        return formatter<string_view>::format(name, ctx);
    }
};

NLOHMANN_JSON_SERIALIZE_ENUM_ARRAY(PricingSolver,
                                   11,
                                   {{bdd_solver_simple, "BddForward"},
                                    {bdd_solver_cycle, "BddForwardCycle"},
                                    {bdd_solver_backward_simple, "BddBackward"},
                                    {bdd_solver_backward_cycle,
                                     "BddBackwardCycle"},
                                    {dp_solver, "TimeIndexed"},
                                    {ati_solver, "ArcTimeIndexed"},
                                    {dp_bdd_solver, "Hybrid"}})

NLOHMANN_JSON_SERIALIZE_ENUM_ARRAY(StabTechniques,
                                   4,
                                   {{no_stab, "NoStabilization"},
                                    {stab_wentgnes, "WentgnesStab"},
                                    {stab_dynamic, "DynamicStab"},
                                    {stab_hybrid, "HybridStab"}})

template <>
struct fmt::formatter<StabTechniques> : formatter<string_view> {
    // parse is inherited from formatter<string_view>.
    template <typename FormatContext>
    auto format(StabTechniques _solver, FormatContext& ctx) {
        string_view name = "unknown";
        switch (_solver) {
            case StabTechniques::no_stab:
                name = "NoStabilization";
                break;
            case StabTechniques::stab_wentgnes:
                name = "WentgnesStab";
                break;
            case StabTechniques::stab_dynamic:
                name = "DynamicStab";
                break;
            case StabTechniques::stab_hybrid:
                name = "HybridStab";
                break;
        }
        return formatter<string_view>::format(name, ctx);
    }
};

NLOHMANN_JSON_SERIALIZE_ENUM_ARRAY(BBExploreStrategy,
                                   5,
                                   {{min_bb_explore_strategy, "dfs"},
                                    {bb_dfs_strategy, "dfs"},
                                    {bb_bfs_strategy, "bfs"},
                                    {bb_brfs_strategy, "brfs"},
                                    {bb_cbfs_strategy, "cbfs"}})

template <>
struct fmt::formatter<BBExploreStrategy> : formatter<string_view> {
    // parse is inherited from formatter<string_view>.
    template <typename FormatContext>
    auto format(BBExploreStrategy _solver, FormatContext& ctx) {
        string_view name = "unknown";
        switch (_solver) {
            case BBExploreStrategy::min_bb_explore_strategy:
                name = "dfs";
                break;
            case BBExploreStrategy::bb_bfs_strategy:
                name = "bfs";
                break;
            case BBExploreStrategy::bb_brfs_strategy:
                name = "brfs";
                break;
            case BBExploreStrategy::bb_cbfs_strategy:
                name = "cbfs";
                break;
        }
        return formatter<string_view>::format(name, ctx);
    }
};

NLOHMANN_JSON_SERIALIZE_ENUM_ARRAY(
    Scoring_Parameter,
    6,
    {{min_scoring_parameter, "ProductScoring"},
     {product_scoring_parameter, "ProductScoring"},
     {min_function_scoring_parameter, "MinFunction"},
     {max_function_scoring_parameter, "MaxFunction"},
     {weighted_sum_scoring_parameter, " WeightedSum"},
     {weighted_product_scoring_parameter, "WeightedProduct"}})

template <>
struct fmt::formatter<Scoring_Parameter> : formatter<string_view> {
    // parse is inherited from formatter<string_view>.
    template <typename FormatContext>
    auto format(Scoring_Parameter _solver, FormatContext& ctx) {
        string_view name = "unknown";
        switch (_solver) {
            case Scoring_Parameter::min_scoring_parameter:
                name = "ProductScoring";
                break;
            case Scoring_Parameter::min_function_scoring_parameter:
                name = "MinFunction";
                break;
            case Scoring_Parameter::max_function_scoring_parameter:
                name = "MaxFunction";
                break;
            case Scoring_Parameter::weighted_sum_scoring_parameter:
                name = " WeightedSum";
                break;
            case Scoring_Parameter::weighted_product_scoring_parameter:
                name = "WeightedProduct";
                break;
        }
        return formatter<string_view>::format(name, ctx);
    }
};

template <>
struct fmt::formatter<Scoring_Value> : formatter<string_view> {
    // parse is inherited from formatter<string_view>.
    template <typename FormatContext>
    auto format(Scoring_Value _solver, FormatContext& ctx) {
        string_view name = "unknown";
        switch (_solver) {
          case Scoring_Value::lb_scoring_value:
                name = "LbScoring";
                break;
          case Scoring_Value::size_scoring_value:
                name = "SizeScoringValue";
                break;
          case Scoring_Value::nb_paths_scoring_value:
                name = "NbPathsScoringValue";
                break;
        }
        return formatter<string_view>::format(name, ctx);
    }
};
struct Parms {
    template <typename T>
    class Parameter {
        std::string _name;
        T           _value;

       public:
        Parameter(std::string name, const T& value)
            : _name(std::move(name)),
              _value(value) {}

        Parameter(const Parameter<T>&) = default;
        Parameter(Parameter<T>&&) noexcept = default;
        auto operator=(Parameter<T>&&) noexcept -> Parameter<T>& = default;
        auto operator=(const Parameter<T>&) -> Parameter<T>& = default;
        ~Parameter() = default;

        [[nodiscard]] inline auto value() const -> const T& { return _value; }
        [[nodiscard]] auto name() const -> const std::string& { return _name; }

        inline auto value() -> T& { return _value; }
        auto        name() -> std::string& { return _name; }

        void set_value(const T& value) { _value = value; }
    };

    Parameter<BBExploreStrategy>                        bb_explore_strategy;
    Parameter<StabTechniques>                           stab_technique;
    Parameter<PricingSolver>                            pricing_solver;
    Parameter<Scoring_Parameter>                        scoring_parameter;
    Parameter<Scoring_Value>                            scoring_value;
    Parameter<int>                                      strong_branching;
    Parameter<int>                                      bb_node_limit;
    Parameter<int>                                      nb_iterations_rvnd;
    Parameter<boost::timer::nanosecond_type>            branching_cpu_limit;
    Parameter<bool>                                     use_cpu_time;
    Parameter<double>                                   alpha;
    Parameter<double>                                   branching_point;
    Parameter<bool>                                     use_heuristic;
    std::function<double(const std::array<double, 2>&)> scoring_function;
    Parameter<bool>                                     use_mip_solver;
    Parameter<bool>                                     refine_bdd;
    Parameter<bool>                                     enumerate;
    Parameter<bool>                                     pruning_test;
    Parameter<bool>                                     suboptimal_duals;
    Parameter<bool>                                     reduce_cost_fixing;
    Parameter<bool>                                     print_csv;
    Parameter<bool>                                     use_bks;

    std::string jobfile{};
    std::string pname{};

    int    nb_jobs{};
    size_t nb_machines{};

    Parms();
    Parms(int argc, const char** argv);

    Parms(const Parms&) = default;
    Parms(Parms&&) = default;
    auto operator=(const Parms&) -> Parms& = default;
    auto operator=(Parms&&) noexcept -> Parms& = default;
    ~Parms() = default;

    void        parms_set_scoring_function(int scoring);
    void        parse_cmd(int argc, const char** argv);
    friend void from_json(const nlohmann::json& j, Parms& p);
    friend void to_json(nlohmann::json& j, const Parms& p);

   private:
    static constexpr double                mu = 5.0 / 6.0;
    static constexpr std::array<double, 2> beta = {1.5, 0.5};
    static constexpr double                TargetBrTimeValue = 0.45;
    static constexpr auto                  EPS = 1e-6;
};

template <typename T>
struct fmt::formatter<Parms::Parameter<T>> {
    char presentation = 'v';

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
        auto aux = std::span<const char>{ctx.begin(), ctx.end()};

        if (aux.empty()) {
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
    constexpr auto format(const Parms::Parameter<T>& p, FormatContext& ctx)
        -> decltype(ctx.out()) {
        if (presentation == 'n') {
            return format_to(ctx.out(), "{}", p.name());
        }
        return format_to(ctx.out(), "{}", p.value());
    }
};

template <>
struct fmt::formatter<Parms> {
    char presentation = 'v';

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
        auto aux = std::span<const char>{ctx.begin(), ctx.end()};

        if (aux.empty()) {
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
    constexpr auto format(const Parms& parms, FormatContext& ctx) -> decltype(ctx.out()) {
        if (presentation == 'n') {
            return format_to(
                ctx.out(),
                "{:n},{:n},{:n},{:n},{:n},{:n},{:n},{:n},{:n},{:n},{:n},{:n},{:n}",
                parms.nb_iterations_rvnd, parms.stab_technique, parms.alpha,
                parms.pricing_solver, parms.strong_branching,
                parms.branching_point, parms.refine_bdd, parms.pruning_test,
                parms.suboptimal_duals, parms.scoring_parameter,
                parms.scoring_value, parms.use_cpu_time,parms.use_bks);
        }
        return format_to(
            ctx.out(), "{},{},{},{},{},{},{},{},{},{},{},{},{}",
            parms.nb_iterations_rvnd, parms.stab_technique, parms.alpha,
            parms.pricing_solver, parms.strong_branching, parms.branching_point,
            parms.refine_bdd, parms.pruning_test, parms.suboptimal_duals,
            parms.scoring_parameter, parms.scoring_value, parms.use_cpu_time,parms.use_bks);
    }
};

#endif  // PARMS_H

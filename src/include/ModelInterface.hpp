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

#ifndef MODEL_INTERFACE
#define MODEL_INTERFACE

#include <fmt/core.h>                     // for print
#include <cmath>                          // for fabs
#include <cstddef>                        // for size_t
#include <iostream>                       // for operator<<, size_t
#include <list>                           // for list
#include <memory>                         // for shared_ptr, __shared_...
#include <range/v3/range/conversion.hpp>  // for to
#include <range/v3/range_fwd.hpp>         // for views
#include <range/v3/view/iota.hpp>         // for iota, iota_fn
#include <range/v3/view/transform.hpp>    // for transform
#include <unordered_map>                  // for operator==, _Node_ite...
#include <utility>                        // for move, pair
#include <variant>                        // for hash
#include <vector>                         // for vector

#include "BddCoeff.hpp"
#include "VariableKeyBase.hpp"

namespace vs = ranges::views;

class ConstraintBase {
    char   sense;
    double rhs;
    bool   can_be_deleted;

   public:
    [[nodiscard]] inline auto get_rhs() const -> double { return rhs; }

    [[nodiscard]] inline auto get_sense() const -> char { return sense; }

    [[nodiscard]] inline auto get_can_be_deleted() const -> bool {
        return can_be_deleted;
    }

    ConstraintBase(const ConstraintBase&) = default;
    ConstraintBase(ConstraintBase&&) = default;
    virtual ~ConstraintBase() = default;

    auto operator=(const ConstraintBase&) -> ConstraintBase& = default;
    auto operator=(ConstraintBase&&) -> ConstraintBase& = default;

    ConstraintBase(char _sense, double _rhs, bool _can_be_delete = false)
        : sense(_sense),
          rhs(_rhs),
          can_be_deleted(_can_be_delete) {}

    virtual auto operator()(const VariableKeyBase&) -> double = 0;
};

class ConstraintAssignment : public ConstraintBase {
   private:
    size_t row;

   public:
    explicit ConstraintAssignment(size_t _row)
        : ConstraintBase('>', 1.0),
          row(_row) {}

    auto operator()(const VariableKeyBase& key) -> double override {
        if (key.get_j() == row && key.get_high()) {
            return 1.0;
        }

        return 0.0;
    }

    ConstraintAssignment(const ConstraintAssignment&) = default;
    ConstraintAssignment(ConstraintAssignment&&) = default;
    ~ConstraintAssignment() override = default;
    auto operator=(const ConstraintAssignment&)
        -> ConstraintAssignment& = default;
    auto operator=(ConstraintAssignment&&) -> ConstraintAssignment& = default;
};

class ConstraintConvex : public ConstraintBase {
   public:
    explicit ConstraintConvex(double _rhs) : ConstraintBase('>', _rhs) {}

    auto operator()(const VariableKeyBase& key) -> double override {
        if (key.get_t() == 0) {
            return -1.0;
        }
        return 0.0;
    }
};

class ReformulationModel : public std::vector<std::shared_ptr<ConstraintBase>> {
   public:
    ReformulationModel(size_t nb_assignments, size_t nb_machines);
    ReformulationModel(ReformulationModel&&) = default;
    ReformulationModel(const ReformulationModel&) = default;
    ~ReformulationModel() = default;
    auto operator=(ReformulationModel&&) -> ReformulationModel& = default;
    auto operator=(const ReformulationModel&) -> ReformulationModel& = default;

    inline void delete_constraint(auto c) {
        if ((*this)[c]->get_can_be_deleted()) {
            (*this)[c].reset();
        }
    }

    inline void delete_constraints(int first, int nb_del) {
        auto it = this->begin() + first;
        this->erase(it, it + nb_del);
    }
};

class GenericData : public std::unordered_map<VariableKeyBase, double> {
   private:
    static constexpr double EPS_GENERIC_DATA = 1e-6;

   public:
    GenericData() = default;
    ~GenericData() = default;
    GenericData(GenericData&&) = default;  // movable and noncopyable
    GenericData(const GenericData&) = default;
    auto operator=(GenericData&&) -> GenericData& = default;
    auto operator=(const GenericData&) -> GenericData& = default;

    void add_coeff_hash_table(size_t _j, size_t _t, bool _high, double _coeff) {
        VariableKeyBase key(_j, _t, _high);

        auto it = this->find(key);
        if (it == this->end()) {
            (*this)[key] = _coeff;
        } else {
            (*this)[key] += _coeff;
        }
    }

    void list_coeff() {
        for (auto& it : (*this)) {
            fmt::print("{} ({},{},{})", it.second, it.first.get_j(),
                       it.first.get_t(), it.first.get_high());
        }
        fmt::print("\n");
    }

    friend auto operator==(const GenericData& lhs, const GenericData& rhs)
        -> bool {
        if (lhs.size() != rhs.size()) {
            return false;
        }

        for (const auto& it1 : lhs) {
            const auto it2 = rhs.find(it1.first);
            if (it2 == rhs.end()) {
                return false;
            }

            if (fabs(it1.second - (*it2).second) > EPS_GENERIC_DATA) {
                return false;
            }
        }


        return true;
    };

    friend auto operator!=(const GenericData& lhs, const GenericData& rhs)
        -> bool {
        return !(lhs == rhs);
    }
};

class ConstraintGeneric : public ConstraintBase {
   private:
    std::shared_ptr<GenericData> data;

   public:
    ConstraintGeneric(GenericData* _data,
                      double       _rhs,
                      char         _sense = '>',
                      bool         _can_be_deleted = true)
        : ConstraintBase(_sense, _rhs, _can_be_deleted),
          data(_data) {}

    explicit ConstraintGeneric(double _rhs,
                               char   _sense = '>',
                               bool   _can_be_deleted = true)
        : ConstraintBase(_sense, _rhs, _can_be_deleted),
          data(nullptr) {}

    ~ConstraintGeneric() override = default;
    ConstraintGeneric(ConstraintGeneric&& op) = default;
    ConstraintGeneric(const ConstraintGeneric&) = default;

    auto operator=(ConstraintGeneric&& op) -> ConstraintGeneric& = default;
    auto operator=(const ConstraintGeneric&) -> ConstraintGeneric& = default;

    auto operator()(const VariableKeyBase& key) -> double override {
        auto it = data->find(key);
        if (it == data->end()) {
            return 0.0;
        }
        return (*it).second;
    }

    friend auto operator!=(const ConstraintGeneric& lhs,
                           const ConstraintGeneric& rhs) -> bool {
        return !(*lhs.data == *rhs.data);
    }

    friend auto operator==(const ConstraintGeneric& lhs,
                           const ConstraintGeneric& rhs) -> bool {
        return (*lhs.data == *rhs.data);
    }

    void list_coeff() { data->list_coeff(); }
};

template <typename T = BddCoeff>
class OriginalConstraint {
   private:
    size_t                        id_constr{};
    std::weak_ptr<ConstraintBase> constr;
    std::list<std::shared_ptr<T>> coeff_list;

   public:
    explicit OriginalConstraint(const std::shared_ptr<ConstraintBase>& _constr)
        : constr(_constr){};
    OriginalConstraint() : constr(){};
    OriginalConstraint(OriginalConstraint&&) noexcept = default;
    OriginalConstraint(const OriginalConstraint<T>&) = default;
    ~OriginalConstraint() = default;
    auto operator=(OriginalConstraint&&) noexcept
        -> OriginalConstraint& = default;
    auto operator=(const OriginalConstraint<T>&)
        -> OriginalConstraint& = default;

    inline auto get_coeff_list() -> std::list<std::shared_ptr<T>>& {
        return coeff_list;
    }

    inline auto get_constr() -> ConstraintBase* {
        auto aux = constr.lock();
        if (aux) {
            return aux.get();
        }
        return nullptr;
    }

    inline void add_coeff_to_list(std::shared_ptr<T> _coeff) {
        coeff_list.push_back(_coeff);
    }

    void clear_coeff() { coeff_list.clear(); }
};

template <typename T = BddCoeff>
class OriginalModel : public std::vector<OriginalConstraint<T>> {
   public:
    explicit OriginalModel(const ReformulationModel& model)
        : std::vector<OriginalConstraint<T>>(
              vs::iota(0UL, model.size()) | vs::transform([&](auto i) {
                  return OriginalConstraint<T>(model[i]);
              }) |
              ranges::to<std::vector<OriginalConstraint<T>>>()) {}

    OriginalModel(OriginalModel&& op) noexcept = default;
    OriginalModel(const OriginalModel&) = default;
    ~OriginalModel() = default;
    auto operator=(OriginalModel&& op) noexcept -> OriginalModel& = default;
    auto operator=(const OriginalModel&) -> OriginalModel& = default;

    void add_coeff_list(int c, std::shared_ptr<T> coeff) {
        (*this)[c].add_coeff_to_list(coeff);
    }

    auto get_constraint(int c) -> ConstraintBase* {
        return (*this)[c].get_constr();
    }

    inline auto get_coeff_list(size_t c) -> std::list<std::shared_ptr<T>>& {
        return (*this)[c].get_coeff_list();
    }

    inline void add_constraint(const std::shared_ptr<ConstraintBase>& _constr) {
        this->push_back(OriginalConstraint<>(_constr));
    }

    inline auto get_nb_constraints() -> size_t { return this->size(); }

    inline void delete_constraints(int first, int nb_del) {
        auto it = this->begin() + first;
        this->erase(it, it + nb_del);
    }

    void clear_all_coeff() {
        for (auto& it : *this) {
            it.clear_coeff();
        }
    }
};

#endif
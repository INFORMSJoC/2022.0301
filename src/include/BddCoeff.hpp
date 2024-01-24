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

//
// Created by daniel on 10/6/21.
//
#include <boost/container_hash/extensions.hpp>  // for hash_combine
#include <cstddef>
#include <ostream>
#include "VariableKeyBase.hpp"
#ifndef PM_BDDCOEFF_H
#define PM_BDDCOEFF_H

class BddCoeff : public VariableKeyBase {
   private:
    size_t row;
    double coeff;
    double value;
    double cum_value;

   public:
    BddCoeff(size_t _j,
             size_t _t,
             double _coeff,
             double _value = 0.0,
             size_t _row = 0UL,
             bool   _high = true,
             bool   _root = false,
             double _cum_value = 0.0);

    BddCoeff(const BddCoeff&) = default;
    BddCoeff(BddCoeff&& op) noexcept;
    ~BddCoeff() override = default;

    auto operator=(const BddCoeff& _src) -> BddCoeff&;
    auto operator=(BddCoeff&& op) noexcept -> BddCoeff& = default;

    [[nodiscard]] auto get_coeff() const -> double;
    [[nodiscard]] auto get_value() const -> double;
    [[nodiscard]] auto get_cum_value() const -> double;
    [[nodiscard]] auto get_row() const -> size_t;

    auto operator==(const BddCoeff& other) -> bool;

    void set_row(size_t _row);
    void update_cum_value(double _value);

    friend auto operator<<(std::ostream& os, const BddCoeff& object)
        -> std::ostream&;
    friend auto operator==(const BddCoeff& lhs, const BddCoeff& rhs) -> bool;
};

namespace std {
template <>
struct hash<BddCoeff> {
    auto operator()(auto const& s) const noexcept -> std::size_t {
        std::size_t seed = 0;
        boost::hash_combine(seed, s.get_j());
        boost::hash_combine(seed, s.get_t());
        boost::hash_combine(seed, s.get_high());

        return seed;  // or use boost::hash_combine
    }
};
}  // namespace std
#endif  // PM_BDDCOEFF_H

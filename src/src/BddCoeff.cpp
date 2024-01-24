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

#include "BddCoeff.hpp"

BddCoeff::BddCoeff(size_t _j,
                   size_t _t,
                   double _coeff,
                   double _value,
                   size_t _row,
                   bool   _high,
                   bool   _root,
                   double _cum_value)
    : VariableKeyBase(_j, _t, _high, _root),
      row(_row),
      coeff(_coeff),
      value(_value),
      cum_value(_cum_value) {}

auto BddCoeff::operator=(const BddCoeff&) -> BddCoeff& = default;
BddCoeff::BddCoeff(BddCoeff&& op) noexcept = default;

auto operator<<(std::ostream& os, const BddCoeff& object) -> std::ostream& {
    return os << "(j = " << object.get_j() << ", t = " << object.get_t()
              << ", x = " << object.get_value()
              << ", high = " << object.get_high() << " )\n";
}

auto operator==(const BddCoeff& lhs, const BddCoeff& rhs) -> bool {
    return lhs.get_j() == rhs.get_j() && lhs.get_t() == rhs.get_t() &&
           lhs.get_high() == rhs.get_high();
}

auto BddCoeff::operator==(const BddCoeff& other) -> bool {
    return get_j() == other.get_j() && get_t() == other.get_t() &&
           get_high() == other.get_high();
}

[[nodiscard]] auto BddCoeff::get_coeff() const -> double {
    return coeff;
}

[[nodiscard]] auto BddCoeff::get_value() const -> double {
    return value;
}

[[nodiscard]] auto BddCoeff::get_cum_value() const -> double {
    return cum_value;
}

[[nodiscard]] auto BddCoeff::get_row() const -> size_t {
    return row;
}

void BddCoeff::set_row(size_t _row) {
    row = _row;
}

void BddCoeff::update_cum_value(double _value) {
    cum_value += _value;
}

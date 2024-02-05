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

#include "VariableKeyBase.hpp"
VariableKeyBase::VariableKeyBase(size_t _j, size_t _t, bool _high, bool _root)
    : j(_j),
      t(_t),
      high(_high),
      root(_root) {}

VariableKeyBase::VariableKeyBase() = default;

[[nodiscard]] auto VariableKeyBase::get_j() const -> size_t {
    return j;
}
[[nodiscard]] auto VariableKeyBase::get_t() const -> size_t {
    return t;
}
[[nodiscard]] auto VariableKeyBase::get_root() const -> bool {
    return root;
}
[[nodiscard]] auto VariableKeyBase::get_high() const -> bool {
    return high;
}

void VariableKeyBase::set_j(size_t _j) {
    j = _j;
}
void VariableKeyBase::set_t(size_t _t) {
    t = _t;
}
void VariableKeyBase::set_root(bool _root) {
    root = _root;
}
void VariableKeyBase::set_high(bool _high) {
    high = _high;
}

VariableKeyBase::VariableKeyBase(const VariableKeyBase& src) = default;
VariableKeyBase::VariableKeyBase(VariableKeyBase&& src) noexcept = default;

auto VariableKeyBase::operator=(const VariableKeyBase& src)
    -> VariableKeyBase& = default;
auto VariableKeyBase::operator=(VariableKeyBase&& src) noexcept
    -> VariableKeyBase& = default;

auto VariableKeyBase::operator==(const VariableKeyBase& other) const -> bool {
    return (j == other.j && t == other.t && high == other.high);
}

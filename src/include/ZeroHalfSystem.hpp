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

#ifndef ZEROHALFSYSTEM_H
#define ZEROHALFSYSTEM_H

#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <cstddef>
#include <vector>
#include "CoinPackedVector.hpp"
#include "ModernDD/NodeBddStructure.hpp"
#include "NodeBdd.hpp"

class ZeroHalfSystem {
    using MatrixDouble = std::vector<std::vector<double>>;
    using VectorDouble = std::vector<double>;
    struct ConstraintData {
        CoinPackedVector        constraint{};
        double                  slack{};
        boost::dynamic_bitset<> row_indices{};
        int                     rhs{};
    };

   public:
    ZeroHalfSystem() = default;
    ZeroHalfSystem(const MatrixDouble& A,
                   const VectorDouble& b,
                   const VectorDouble& _x);
    ZeroHalfSystem(const DdStructure<NodeBdd>& bdd);
    ZeroHalfSystem(ZeroHalfSystem&&) = default;
    ZeroHalfSystem(const ZeroHalfSystem&) = default;
    ~ZeroHalfSystem() = default;

    auto operator=(ZeroHalfSystem&&) -> ZeroHalfSystem& = default;
    auto operator=(const ZeroHalfSystem&) -> ZeroHalfSystem& = default;

   private:
    std::vector<CoinPackedVector>        assignment_constraints;
    std::vector<double>                  assignment_slack;
    std::vector<boost::dynamic_bitset<>> assignment_index;
    std::vector<int>                     assignment_rhs;

    std::vector<CoinPackedVector>        flow_constraints;
    std::vector<double>                  flow_slack;
    std::vector<boost::dynamic_bitset<>> flow_index;
    std::vector<int>                     flow_rhs;

    int nb_assignments{};
    int nb_nz_vertices{};
    int nb_nz_edges{};

    std::vector<std::vector<int>> A_bar;
    std::vector<int>              b_bar;
    std::vector<double>           x_star;
    std::vector<double>           slack;

    std::vector<boost::dynamic_bitset<>> row_index;

    size_t nb_rows{};
    size_t nb_columns{};

    static constexpr double HALF = 2.0;
    static constexpr double EPS = 1e-6;

    void remove_row(size_t _row);
    void remove_col(size_t _col);

    void reduce_system();
    void reduce_gauss();
    void evaluate_rows(const std::vector<int>& _rows);

    void add_to_row(size_t i, size_t j);

   public:
    void generate_cuts();
};

#endif  // ZEROHALFSYSTEM_H
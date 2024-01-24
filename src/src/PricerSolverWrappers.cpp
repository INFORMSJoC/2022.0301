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

#include <cassert>                   // for assert
#include <memory>                    // for unique_ptr
#include <vector>                    // for vector
#include "NodeData.h"                // for NodeData
#include "PricerSolverBase.hpp"      // for PricerSolverBase
#include "PricingStabilization.hpp"  // for PricingStabilizationBase
#include "or-utils/lp.h"             // for lp_interface_deleterows, lp_inte...

void NodeData::build_solve_mip() const {
    solver->build_mip();
}

void NodeData::construct_lp_sol_from_rmp() {
    assert(osi_rmp->getNumCols() - id_pseudo_schedules ==
           static_cast<int>(localColPool.size()));

    lambda = std::span<const double>{osi_rmp->getColSolution(),
                                     static_cast<size_t>(osi_rmp->getNumCols())};
    solver->construct_lp_sol_from_rmp(lambda.subspan(id_pseudo_schedules),
                                      localColPool);
}

void NodeData::generate_cuts() const {
    // 1. add cuts to reformulation model

    solver->add_constraints();
    // solver->insert_constraints_lp(this);
    // solver->update_coeff_constraints();
    // 2. add cuts to lp relaxation wctlp
    // 3. adjust the pricing solver (add constraints to original model)
}

auto NodeData::delete_unused_rows_range(int first, int last) -> int {
    int val = 0;

    // lp_interface_deleterows(RMP.get(), first, last);
    solver->remove_constraints(first, last - first + 1);
    solver_stab->remove_constraints(first, last - first + 1);
    // pi.erase(pi.begin(), pi.begin() + last - first + 1);
    // rhs.erase(rhs.begin(), rhs.begin() + last - first + 1);
    lhs_coeff.erase(lhs_coeff.begin(), lhs_coeff.begin() + last - first + 1);

    return val;
}

void NodeData::call_update_rows_coeff() const  {

    solver->update_rows_coeff(nb_jobs + 1);

}

//
// Created by sam on 18/10/24.
//

#ifndef REWEIGHT_H
#define REWEIGHT_H

#include "aligned_vec.h"

namespace recombine {


using index_iterator = typename VECTORI::iterator;
using data_iterator = typename VECTORD::iterator;

void reweight(
        index_iterator rIDXb, //->an index of length rP for the rows in P and K
        index_iterator cIDXb, //->an index of length cK for the columns in K
        data_iterator P, //->a positive measure on rP points with strictly positive mass at at least one of the points
        data_iterator K, //->cK independent real column vectors on the same rP points - given column by column - (and
                         //tested only in the case with sum over each row is 1)
        integer rP,
        integer ldk, // the stride in K
        integer cK,
        doublereal probability_zero_tolerance // UNUSED an empirical measure of the errors rounding gives to probability
                                              // measure calculations
);


} // namespace recombine

#endif // REWEIGHT_H

/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
  ============================================================================
  Copyright (C) 2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  ============================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ============================================================================
  Implementation of simple histograms.
  ==========================================================================*/

// TODO(holtgrew): Add iterators interface?

#ifndef READ_ANALYZER_HISTOGRAM_H_
#define READ_ANALYZER_HISTOGRAM_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

// ============================================================================
// Enums, Classes
// ============================================================================

struct Dense_;
typedef Tag<Dense_> Dense;

/*
.Class.Histogram
..summary:Stores number of occurences for a contiguous interval of integers beginning in 0.
..remark:The size is automatically updated.
*/
template <typename TSpec>
struct Histogram;

template <>
struct Histogram<Dense>
{
public:
    String<double> counters;

    Histogram()
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

void tally(Histogram<Dense> & histogram, size_t value, double weight)
{
    if (value >= length(histogram.counters))
        resize(histogram.counters, value + 1, 0);
    histogram.counters[value] += weight;
}

void tally(Histogram<Dense> & histogram, size_t value)
{
    tally(histogram, value, 1.0);
}

#endif  // #ifndef READ_ANALYZER_HISTOGRAM_H_

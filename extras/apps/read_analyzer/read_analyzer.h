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
  Shared code for the Read Analyzer program.
  ==========================================================================*/

#ifndef READ_ANALYZER_READ_ANALYZER_H_
#define READ_ANALYZER_READ_ANALYZER_H_

#include <algorithm>
#include <sstream>

using namespace seqan;

// ============================================================================
// Enums, Classes
// ============================================================================

/*
.Class.ReadEvaluationResult
..summary:Stores base and quality counts for reads.
*/
template <typename TSpec>
struct ReadEvaluationResult;

/*
.Class.AlignmentEvaluationResult
..summary:Stores base and quality counts for reads.
*/
template <typename TSpec>
struct AlignmentEvaluationResult;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================


#endif  // READ_ANALYZER_READ_ANALYZER_H_

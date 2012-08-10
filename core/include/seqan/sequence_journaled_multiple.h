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
  Author: Rene Maerker <rene.maerker@fu-berlin.de>
  20.12.2010
  ==========================================================================*/

#ifndef SEQUENCE_JOURNALED_MULTIPLE_H_
#define SEQUENCE_JOURNALED_MULTIPLE_H_

// ============================================================================
// Prerequisites.
// ============================================================================

#include <algorithm>
#include <map>
#include <ostream>
#include <string>
#include <sstream>
#include <queue>

#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include <seqan/sequence.h>

// ============================================================================
// Forwards.
// ============================================================================

//#include <seqan/sequence_journaled/sequence_journaled_forwards.h>
//
//#ifdef SEQAN_SWITCH_USE_FORWARDS
//#include <seqan/sequence_journaled/sequence_journaled_generated_forwards.h>
//#endif

// ============================================================================
// Journaled Sequences.
// ============================================================================

#include <seqan/sequence_journaled_multiple/string_set_journaled_group.h>
#include <seqan/sequence_journaled_multiple/string_set_journaled_basic.h>
#include <seqan/sequence_journaled_multiple/score_journal.h>

#include <seqan/sequence_journaled_multiple/string_set_journaled_join_trace.h>
#include <seqan/sequence_journaled_multiple/string_set_journaled_join_base.h>

#include <seqan/sequence_journaled_multiple/string_set_journaled_join_simple.h>
#include <seqan/sequence_journaled_multiple/string_set_journaled_join_max_compression.h>
#include <seqan/sequence_journaled_multiple/string_set_journaled_join_synchronization.h>
#include <seqan/sequence_journaled_multiple/string_set_journaled_join_any_alignment.h>
#include <seqan/sequence_journaled_multiple/string_set_journaled_join_chaining.h>



#endif /* SEQUENCE_JOURNALED_MULTIPLE_H_ */

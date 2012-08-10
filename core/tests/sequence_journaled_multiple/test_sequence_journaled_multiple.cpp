/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de
  ===========================================================================
  Copyright (C) 2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  ===========================================================================
  Author: Rene Maerker <rene.maerker@fu-berlin.de>
  ===========================================================================
  Tests for the sequence_journaled_multiple module.
  ===========================================================================
*/

#include <seqan/basic.h>
#include <seqan/file.h>

#include "test_sequence_journaled_multiple.h"


SEQAN_BEGIN_TESTSUITE(test_sequence_journaled) {
    // Call tests of the sequence journal with unbalanced tree journal.
	SEQAN_CALL_TEST(test_string_set_journaled_group_basic);
	SEQAN_CALL_TEST(test_string_set_journaled_basic);
	SEQAN_CALL_TEST(test_sequence_journaled_is_empty);
	SEQAN_CALL_TEST(test_string_set_journaled_join);
    // Verify checkpoints.
	SEQAN_VERIFY_CHECKPOINTS("projects\\library\\seqan\\sequence_journaled_multiple\\string_set_journaled_basic.h");
	SEQAN_VERIFY_CHECKPOINTS("projects\\library\\seqan\\sequence_journaled_multiple\\string_set_journaled_group.h");
	SEQAN_VERIFY_CHECKPOINTS("projects\\library\\seqan\\sequence_journaled_multiple\\score_journal.h");
	SEQAN_VERIFY_CHECKPOINTS("projects\\library\\seqan\\sequence_journaled_multiple\\string_set_journaled_join_base.h");
	SEQAN_VERIFY_CHECKPOINTS("projects\\library\\seqan\\sequence_journaled_multiple\\string_set_journaled_join_max_compression.h");
//	SEQAN_VERIFY_CHECKPOINTS("projects\\library\\seqan\\sequence_journaled_multiple\\string_set_journaled_util.h");

}
SEQAN_END_TESTSUITE


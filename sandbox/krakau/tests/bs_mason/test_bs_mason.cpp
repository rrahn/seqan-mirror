// ==========================================================================
//                                  bs_mason
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Your Name <your.email@example.net>
// ==========================================================================

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include <seqan/basic.h>
#include <seqan/find.h>
#include <seqan/misc/misc_cmdparser.h>
#include <seqan/misc/misc_random.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <seqan/store.h>
#include <seqan/stream.h>
#include <seqan/file.h>
#include <seqan/refinement.h>

#include "../../apps/bs_mason/simulate_illumina.h"
#include "../../apps/bs_mason/simulate_454.h"
#include "../../apps/bs_mason/simulate_sanger.h"
#include "../../apps/bs_mason/simulate_illumina_bs.h"
#include "../../apps/bs_mason/bs_mason.h"

#include "test_bs_mason.h"


SEQAN_BEGIN_TESTSUITE(test_bs_mason)
{
    // No given methylation rates, simulate methylation positions

	SEQAN_CALL_TEST(test_bs_mason_simulateMethylPositions);

	SEQAN_CALL_TEST(test_bs_mason_buildHaplotype);

    SEQAN_CALL_TEST(test_bs_mason_buildMethylHaplotypeUseMethylPositions);

	SEQAN_CALL_TEST(test_bs_mason_buildSimulationInstructions);

	SEQAN_CALL_TEST(test_bs_mason_applySimulationInstructions);

	SEQAN_CALL_TEST(test_bs_mason_checkOutputSimMR);

	SEQAN_CALL_TEST(test_bs_mason_checkOutputSimMR_PE);

    // Use given methylation rates

    SEQAN_CALL_TEST(test_bs_mason_buildMethylHaplotypeUseMethylRates);

	SEQAN_CALL_TEST(test_bs_mason_checkOutputUseMR);

	SEQAN_CALL_TEST(test_bs_mason_checkOutputUseMR_PE);
}
SEQAN_END_TESTSUITE

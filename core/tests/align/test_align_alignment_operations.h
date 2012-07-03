// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Author: Birte Kehr <birte.kehr@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_CORE_TESTS_ALIGN_TEST_ALIGN_ALIGMENT_OPERATIONS_H_
#define SEQAN_CORE_TESTS_ALIGN_TEST_ALIGN_ALIGMENT_OPERATIONS_H_

SEQAN_DEFINE_TEST(test_align_integrate_align)
{
    using namespace seqan;

    // We do not test the first interface since the second one calls the first
    // one implicitely.

    // Case: Integration left of existing with leading gaps. Test second
    // specialization: One Align of Sequence, one Align of Infix of Sequence.
    {
        Dna5String seqH = "AAAACGATAAAACGAT";
        Dna5String seqV =     "CGAT"  "CGAT";

        // Alignment align:
        //
        // AAAACGATAAAACGAT
        //     CGAT----CGAT

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), seqH);
        assignSource(row(align, 1), seqV);
        insertGaps(row(align, 1), 4, 4);

        // Alignment alignInf
        //
        // AAAACGAT
        // ----CGAT

        Align<typename Infix<Dna5String>::Type> alignInf;
        resize(rows(alignInf), 2);
        assignSource(row(alignInf, 0), infix(seqH, 0, 8));
        assignSource(row(alignInf, 1), infix(seqV, 0, 4));
        insertGaps(row(alignInf, 1), 0, 4);

        integrateAlign(align, alignInf);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAACGATAAAACGAT");
        SEQAN_ASSERT_EQ(ssV.str(), "----CGAT----CGAT");
    }

    // Case: Integration left of existing with trailing gaps. Test second
    // specialization: One Align of Sequence, one Align of Infix of Sequence.
    {
        Dna5String seqH = "CGATATAAAAAACGAT";
        Dna5String seqV = "CGATT"     "CGAT";

        // Alignment align:
        //
        // CGATATAAAAAACGAT
        //    CGATT----CGAT

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), seqH);
        assignSource(row(align, 1), seqV);
        insertGaps(row(align, 1), 5, 4);

        // Alignment alignInf
        //
        // CGATATAA
        // CGAT-T--

        Align<typename Infix<Dna5String>::Type> alignInf;
        resize(rows(alignInf), 2);
        assignSource(row(alignInf, 0), infix(seqH, 0, 8));
        assignSource(row(alignInf, 1), infix(seqV, 0, 5));
        insertGaps(row(alignInf, 1), 5, 2);
        insertGaps(row(alignInf, 1), 4, 1);

        integrateAlign(align, alignInf);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "CGATATAAAAAACGAT");
        SEQAN_ASSERT_EQ(ssV.str(), "CGAT-T------CGAT");
    }

    // Case: Integration right of existing with leading gaps. Test second
    // specialization: One Align of Sequence, one Align of Infix of Sequence.
    {
        Dna5String seqH = "AAAACGATAAAACGAT";
        Dna5String seqV =     "CGAT"  "CGAT";

        // Alignment align:
        //
        // AAAACGATAAAACGAT
        // ----CGATCGAT

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), seqH);
        assignSource(row(align, 1), seqV);
        insertGaps(row(align, 1), 0, 4);

        // Alignment alignInf
        //
        // AAAACGAT
        // ----CGAT

        Align<typename Infix<Dna5String>::Type> alignInf;
        resize(rows(alignInf), 2);
        assignSource(row(alignInf, 0), infix(seqH, 8, 16));
        assignSource(row(alignInf, 1), infix(seqV, 4, 8));
        insertGaps(row(alignInf, 1), 0, 4);

        integrateAlign(align, alignInf);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAACGATAAAACGAT");
        SEQAN_ASSERT_EQ(ssV.str(), "----CGAT----CGAT");
    }

    // Case: Integration right of existing with trailing gaps. Test second
    // specialization: One Align of Sequence, one Align of Infix of Sequence.
    {
        Dna5String seqH = "AAAACGATCGATATAA";
        Dna5String seqV =     "CGATCGATT";

        // Alignment align:
        //
        // AAAACGATCGATATAA
        // ----CGATCGATT

        Align<Dna5String> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), seqH);
        assignSource(row(align, 1), seqV);
        insertGaps(row(align, 1), 0, 4);

        // Alignment alignInf
        //
        // CGATATAA
        // CGAT-T--

        Align<typename Infix<Dna5String>::Type> alignInf;
        resize(rows(alignInf), 2);
        assignSource(row(alignInf, 0), infix(seqH, 8, 16));
        assignSource(row(alignInf, 1), infix(seqV, 4, 9));
        insertGaps(row(alignInf, 1), 5, 2);
        insertGaps(row(alignInf, 1), 4, 1);

        integrateAlign(align, alignInf);

        std::stringstream ssH, ssV;
        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAAACGATCGATATAA");
        SEQAN_ASSERT_EQ(ssV.str(), "----CGATCGAT-T--");
    }
}

#endif  // #ifndef SEQAN_CORE_TESTS_ALIGN_TEST_ALIGN_ALIGMENT_OPERATIONS_H_

// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_TRACEBACK_ADAPTOR_H_
#define CORE_TESTS_ALIGN_TEST_ALIGNMENT_TRACEBACK_ADAPTOR_H_

#include <sstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/align/gaps_base.h>
#include <seqan/align/gaps_iterator_base.h>
#include <seqan/align/gaps_array.h>
#include <seqan/align/gaps_iterator_array.h>
#include <seqan/align/gaps_anchor.h>
#include <seqan/align/gaps_iterator_anchor.h>

#include <seqan/align/align_base.h>
#include <seqan/align/alignment_base.h>
#include <seqan/align/alignment_traceback_tracesegment.h>
#include <seqan/align/alignment_traceback_adaptor.h>

using namespace seqan;

void testAlignmentTracebackAdaptorAdaptFile()
{
    std::stringstream stream;
    typedef TraceSegment<size_t, size_t> TTraceSegment;
    String<TTraceSegment> traceSegments;

    appendValue(traceSegments, TTraceSegment(12,8, 4, +TraceBitMask::VERTICAL));
    appendValue(traceSegments, TTraceSegment(8,4, 4, +TraceBitMask::DIAGONAL));
    appendValue(traceSegments, TTraceSegment(8,3, 1, +TraceBitMask::VERTICAL));
    appendValue(traceSegments, TTraceSegment(4,3, 4, +TraceBitMask::HORIZONTAL));
    appendValue(traceSegments, TTraceSegment(1,0, 3, +TraceBitMask::DIAGONAL));
    appendValue(traceSegments, TTraceSegment(0,0, 1, +TraceBitMask::HORIZONTAL));

    String<char> seq0 ="AAAACCCCGGGG";
    String<char> seq1 ="AAAACCCCGGGG";

    adapt(stream, seq0, seq1, traceSegments);
    String<char> result = stream.str();

    SEQAN_ASSERT_EQ(result[1], 'A');
    SEQAN_ASSERT_EQ(result[3], gapValue<char>());

    SEQAN_ASSERT_EQ(result[7], 'A');
    SEQAN_ASSERT_EQ(result[9], 'A');

    SEQAN_ASSERT_EQ(result[13], 'A');
    SEQAN_ASSERT_EQ(result[15], 'A');

    SEQAN_ASSERT_EQ(result[19], 'A');
    SEQAN_ASSERT_EQ(result[21], 'A');

    SEQAN_ASSERT_EQ(result[25], 'C');
    SEQAN_ASSERT_EQ(result[27], gapValue<char>());

    SEQAN_ASSERT_EQ(result[31], 'C');
    SEQAN_ASSERT_EQ(result[33], gapValue<char>());

    SEQAN_ASSERT_EQ(result[37], 'C');
    SEQAN_ASSERT_EQ(result[39], gapValue<char>());

    SEQAN_ASSERT_EQ(result[43], 'C');
    SEQAN_ASSERT_EQ(result[45], gapValue<char>());

    SEQAN_ASSERT_EQ(result[49], gapValue<char>());
    SEQAN_ASSERT_EQ(result[51], 'A');

    SEQAN_ASSERT_EQ(result[55], 'G');
    SEQAN_ASSERT_EQ(result[57], 'C');

    SEQAN_ASSERT_EQ(result[61], 'G');
    SEQAN_ASSERT_EQ(result[63], 'C');

    SEQAN_ASSERT_EQ(result[67], 'G');
    SEQAN_ASSERT_EQ(result[69], 'C');

    SEQAN_ASSERT_EQ(result[73], 'G');
    SEQAN_ASSERT_EQ(result[75], 'C');

    SEQAN_ASSERT_EQ(result[79], gapValue<char>());
    SEQAN_ASSERT_EQ(result[81], 'G');

    SEQAN_ASSERT_EQ(result[85], gapValue<char>());
    SEQAN_ASSERT_EQ(result[87], 'G');

    SEQAN_ASSERT_EQ(result[91], gapValue<char>());
    SEQAN_ASSERT_EQ(result[93], 'G');

    SEQAN_ASSERT_EQ(result[97], gapValue<char>());
    SEQAN_ASSERT_EQ(result[99], 'G');
}

template <typename TGapsSpec>
void testAlignmentTracebackAdaptorAdaptGaps(TGapsSpec const &)
{

    typedef TraceSegment<size_t, size_t> TTraceSegment;
    typedef Align<CharString, TGapsSpec > TAlign;


    {
        String<TTraceSegment> traceSegments;
        appendValue(traceSegments, TTraceSegment(12,8, 4, +TraceBitMask::VERTICAL));
        appendValue(traceSegments, TTraceSegment(8,4, 4, +TraceBitMask::DIAGONAL));
        appendValue(traceSegments, TTraceSegment(8,3, 1, +TraceBitMask::VERTICAL));
        appendValue(traceSegments, TTraceSegment(4,3, 4, +TraceBitMask::HORIZONTAL));
        appendValue(traceSegments, TTraceSegment(1,0, 3, +TraceBitMask::DIAGONAL));
        appendValue(traceSegments, TTraceSegment(0,0, 1, +TraceBitMask::HORIZONTAL));

        CharString seq0 ="AAAACCCCGGGG";
        CharString seq1 ="AAAACCCCGGGG";

        TAlign align;

        resize(rows(align),2);
        assignSource(row(align, 0),seq0);
        assignSource(row(align, 1),seq1);

        // prepare align object to test the correctness of the adapt function.
        TAlign compareAlign;
        resize(rows(compareAlign),2);
        assignSource(row(compareAlign, 0),seq0);
        assignSource(row(compareAlign, 1),seq1);

        insertGap(row(compareAlign,0),8);
        insertGaps(row(compareAlign,0),13, 4);

        insertGap(row(compareAlign,1),0);
        insertGaps(row(compareAlign,1),4,4);

        adapt(row(align,0), row(align,1), traceSegments);

        SEQAN_ASSERT_EQ(row(align,0), row(compareAlign, 0));
        SEQAN_ASSERT_EQ(row(align,1), row(compareAlign, 1));
    }

    {
        String<TTraceSegment> traceSegments;
        appendValue(traceSegments, TTraceSegment(8,12, 4, +TraceBitMask::HORIZONTAL));
        appendValue(traceSegments, TTraceSegment(4,8, 4, +TraceBitMask::DIAGONAL));
        appendValue(traceSegments, TTraceSegment(3,8, 1, +TraceBitMask::HORIZONTAL));
        appendValue(traceSegments, TTraceSegment(3,4, 4, +TraceBitMask::VERTICAL));
        appendValue(traceSegments, TTraceSegment(0,1, 3, +TraceBitMask::DIAGONAL));
        appendValue(traceSegments, TTraceSegment(0,0, 1, +TraceBitMask::VERTICAL));

        CharString seq0 ="AAAACCCCGGGG";  //-AAA----ACCCCGGGG
        CharString seq1 ="AAAACCCCGGGG";  //AAAACCCC-GGGG----

        TAlign align;

        resize(rows(align),2);
        assignSource(row(align, 0),seq0);
        assignSource(row(align, 1),seq1);

        // prepare align object to test the correctness of the adapt function.
        TAlign compareAlign;
        resize(rows(compareAlign),2);
        assignSource(row(compareAlign, 0),seq0);
        assignSource(row(compareAlign, 1),seq1);

        insertGap(row(compareAlign,1),8);
        insertGaps(row(compareAlign,1),13, 4);

        insertGap(row(compareAlign,0),0);
        insertGaps(row(compareAlign,0),4,4);

        adapt(row(align,0), row(align,1), traceSegments);

        SEQAN_ASSERT_EQ(row(align,0), row(compareAlign, 0));
        SEQAN_ASSERT_EQ(row(align,1), row(compareAlign, 1));
    }

    {
        String<TTraceSegment> traceSegments;
        appendValue(traceSegments, TTraceSegment(8,4, 4, +TraceBitMask::DIAGONAL));
        appendValue(traceSegments, TTraceSegment(8,3, 1, +TraceBitMask::VERTICAL));
        appendValue(traceSegments, TTraceSegment(4,3, 4, +TraceBitMask::HORIZONTAL));

        CharString seq0 ="AAAACCCCGGGG";
        CharString seq1 ="AAAACCCCGGGG";

        //01234567890123456
        //AAAACCCC-GGGG----
        //-AAA----ACCCCGGGG
        //10123456780123456
        //-
        //    XXXXXXXXX

        TAlign align;

        resize(rows(align),2);
        assignSource(row(align, 0),seq0);
        assignSource(row(align, 1),seq1);

        // prepare align object to test the correctness of the adapt function.
        TAlign compareAlign;
        resize(rows(compareAlign),2);
        assignSource(row(compareAlign, 0),seq0);
        assignSource(row(compareAlign, 1),seq1);

        insertGap(row(compareAlign,0),8);
        insertGaps(row(compareAlign,0),13, 4);

        insertGap(row(compareAlign,1),0);
        insertGaps(row(compareAlign,1),4,4);

        // TODO(rmaerker): leading gaps in anchor gaps are enumerated with negative integers but should have the same interface as array gaps. Dave needs to fix this.
        setClippedEndPosition(row(compareAlign, 0), 13 + clippedBeginPosition(row(compareAlign, 0)));
        setClippedBeginPosition(row(compareAlign, 0), 4 + clippedBeginPosition(row(compareAlign, 0)));

        // TODO(rmaerker): leading gaps in anchor gaps are enumerated with negative integers but should have the same interface as array gaps. Dave needs to fix this.
        setClippedEndPosition(row(compareAlign, 1), 13 + clippedBeginPosition(row(compareAlign, 1)));
        setClippedBeginPosition(row(compareAlign, 1), 4 + clippedBeginPosition(row(compareAlign, 1)));


        adapt(row(align,0), row(align,1), traceSegments);

        SEQAN_ASSERT_EQ(row(align,0), row(compareAlign, 0));
        SEQAN_ASSERT_EQ(row(align,1), row(compareAlign, 1));
    }

}


void testAlignmentTracebackAdaptorAdaptFragments()
{
    // TODO (rmaerker): write me!
}

void testAlignmentTracebackAdaptorAdaptAlignmentGraph()
{
    // TODO (rmaerker): write me!
}

void testAlignmentTracebackAdaptorAdaptVertexDescriptor()
{
    // TODO (rmaerker): write me!
}

SEQAN_DEFINE_TEST(test_alignment_traceback_adaptor_adapt_file)
{
    testAlignmentTracebackAdaptorAdaptFile();
}
SEQAN_DEFINE_TEST(test_alignment_traceback_adaptor_adapt_gaps)
{
    testAlignmentTracebackAdaptorAdaptGaps(ArrayGaps());
    testAlignmentTracebackAdaptorAdaptGaps(AnchorGaps<>());
}
SEQAN_DEFINE_TEST(test_alignment_traceback_adaptor_adapt_fragments)
{
    testAlignmentTracebackAdaptorAdaptFragments();
}
SEQAN_DEFINE_TEST(test_alignment_traceback_adaptor_adapt_alignment_graph)
{
    testAlignmentTracebackAdaptorAdaptAlignmentGraph();
}
SEQAN_DEFINE_TEST(test_alignment_traceback_adaptor_adapt_vertex_descriptor)
{
    testAlignmentTracebackAdaptorAdaptVertexDescriptor();
}

#endif  // #ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_TRACEBACK_ADAPTOR_H_

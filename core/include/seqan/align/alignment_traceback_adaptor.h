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

#ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_TRACEBACK_ADAPTOR_H_
#define CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_TRACEBACK_ADAPTOR_H_

namespace seqan {

// ----------------------------------------------------------------------------
// fucntion _writeTraceSegmentToFile()
// ----------------------------------------------------------------------------

template <typename TFile, typename TSeq0Value, typename TSeq1Value>
inline void _writeTraceSegmentToFile(TFile & file, TSeq0Value const & seq0Val, TSeq1Value const & seq1Val)
{
    _streamPut(file, '(');
    _streamPut(file, seq0Val);
    _streamPut(file, ',');
    _streamPut(file, seq1Val);
    _streamPut(file, ')');
    _streamPut(file, '\n');
}

// ----------------------------------------------------------------------------
// function adapt()
// ----------------------------------------------------------------------------

template <typename TSourceHorizontal, typename TGapsSpecHorizontal, typename TSourceVertical,
        typename TGapsSpecVertical, typename TPosition, typename TSize, typename TStringSpec>
void
adapt(Gaps<TSourceHorizontal, TGapsSpecHorizontal> & gapsHorizontal,
      Gaps<TSourceVertical, TGapsSpecVertical> & gapsVertical,
      String<TraceSegment<TPosition, TSize>, TStringSpec > const & traceSegments)
{
    typedef Gaps<TSourceHorizontal, TGapsSpecHorizontal> TGapsHorizontal;
    typedef Gaps<TSourceVertical, TGapsSpecVertical> TGapsVertical;
    typedef typename Iterator<TGapsHorizontal>::Type TIteratorHorizontal;
    typedef typename Iterator<TGapsVertical>::Type TIteratorVertical;
    typedef TraceSegment<TPosition, TSize> TTraceSegment;
    typedef typename Iterator<String<TTraceSegment, TStringSpec> const>::Type TTraceIterator;

    // TODO (rmaerker): specialize for different Align structures
    clearGaps(gapsHorizontal);
    clearClipping(gapsHorizontal);

    clearGaps(gapsVertical);
    clearClipping(gapsVertical);

    TTraceIterator srcIter = end(traceSegments) - 1;
    TTraceIterator srcEnd = begin(traceSegments) - 1;

    // we build the gap structure here.
    // set the clipped begin position of the alignment.
    setBeginPosition(gapsHorizontal, getBeginHorizontal(value(srcIter))); // begin of source
    setBeginPosition(gapsVertical, getBeginVertical(value(srcIter)));

    TIteratorHorizontal it0 = begin(gapsHorizontal);
    TIteratorVertical it1 = begin(gapsVertical);

    while (srcIter != srcEnd)
    {
        TSize segmentSize = value(srcIter)._length;
        switch (value(srcIter)._traceValue)
        {
            case TraceBitMask::HORIZONTAL:
                insertGaps(it1, segmentSize);
                break;

            case TraceBitMask::VERTICAL:
                insertGaps(it0, segmentSize);
                break;
        }
        goFurther(it0, segmentSize);
        goFurther(it1, segmentSize);
        --srcIter;
    }
    setClippedEndPosition(gapsHorizontal, position(it0) + clippedBeginPosition(gapsHorizontal));
    setClippedEndPosition(gapsVertical, position(it1) + clippedBeginPosition(gapsVertical));
}

template <typename TStringSet, typename TCargo, typename TSpec, typename TSequenceIdH, typename TSequenceIdV,
        typename TPosition, typename TSize, typename TStringSpec>
inline void
adapt(Graph<Alignment<TStringSet, TCargo, TSpec> > & g,
        TSequenceIdH const & seqHId,
        TSequenceIdV const & seqVId,
        String<TraceSegment<TPosition, TSize>, TStringSpec > const & traceSegments)
{
    typedef TraceSegment<TPosition, TSize> TTraceSegment;
    // check begin and end positions of the graph.

    // insert leading gaps
    TTraceSegment traceBegin = traceSegments[length(traceSegments) - 1];
    if (getBeginVertical(traceBegin) != 0)
    {
        addVertex(g, seqVId, 0, getBeginVertical(traceBegin));
    }
    if (getBeginHorizontal(traceBegin) != 0)
    {
        addVertex(g, seqHId, 0, getBeginHorizontal(traceBegin));
    }


    for (TSize i = 0; i < length(traceSegments); ++i)
    {

        switch (traceSegments[i]._traceValue)
        {
        case TraceBitMask::DIAGONAL:
            addEdge(g, addVertex(g, seqHId, traceSegments [i]._horizontalBeginPos, traceSegments [i]._length),
                    addVertex(g, seqVId, traceSegments [i]._verticalBeginPos, traceSegments [i]._length));
            break;
        case TraceBitMask::VERTICAL:
            addVertex(g, seqVId, traceSegments [i]._verticalBeginPos, traceSegments [i]._length);
            break;
        case TraceBitMask::HORIZONTAL:
            addVertex(g, seqHId, traceSegments [i]._horizontalBeginPos, traceSegments [i]._length);
            break;
        }
    }

    // insert trailing gaps
    TTraceSegment traceEnd = traceSegments [0];
    if (getEndVertical(traceEnd) != length(value(stringSet(g), idToPosition(stringSet(g), seqVId))))
    {
        addVertex(g, seqVId, getEndVertical(traceEnd),
                  length(value(stringSet(g), idToPosition(stringSet(g), seqVId))) - getEndVertical(traceEnd));
    }
    if (getEndHorizontal(traceEnd) != length(value(stringSet(g), idToPosition(stringSet(g), seqHId))))
    {
        addVertex(g, seqHId, getEndHorizontal(traceEnd),
                  length(value(stringSet(g), idToPosition(stringSet(g), seqHId))) - getEndHorizontal(traceEnd));
    }
}

template <typename TFile, typename TSequenceH, typename TSequenceV, typename TPosition, typename TSize,
        typename TStringSpec>
inline void
adapt(TFile & file,
      TSequenceH const & seqH,
      TSequenceV const & seqV,
      String<TraceSegment<TPosition, TSize>, TStringSpec> const & traceSegments)
{
    for (TSize k = length(traceSegments); k > (TSize) 0; --k)
    {
        switch (traceSegments[k-1]._traceValue)
        {
            case TraceBitMask::DIAGONAL:
            {
                int j = traceSegments[k-1]._verticalBeginPos;
                for (int i = traceSegments[k-1]._horizontalBeginPos; i < (int) (traceSegments[k-1]._horizontalBeginPos + traceSegments[k-1]._length); ++i)
                {
                    _writeTraceSegmentToFile(file, seqH[i], seqV[j]);
                    ++j;
                }
                break;
            }
            case TraceBitMask::VERTICAL:
            {
                for (int i = traceSegments[k-1]._verticalBeginPos; i < (int) (traceSegments[k-1]._verticalBeginPos + traceSegments[k-1]._length); ++i)
                {
                    _writeTraceSegmentToFile(file, gapValue<char>(), seqV[i]);
                }
                break;
            }
            case TraceBitMask::HORIZONTAL:
            {
                for (int i = traceSegments[k-1]._horizontalBeginPos; i < (int) (traceSegments[k-1]._horizontalBeginPos + traceSegments[k-1]._length) ; ++i)
                {
                    _writeTraceSegmentToFile(file, seqH[i], gapValue<char>());
                }
                break;
            }
        }
    }
}

template <typename TSize, typename TFragmentSpec, typename TStringSpec, typename TSequenceH, typename TSequenceV,
         typename TPosition, typename TSize2, typename TStringSpec2>
inline void
adapt(String<Fragment<TSize, TFragmentSpec>, TStringSpec>& matches,
        TSequenceH const & seqHId,
        TSequenceV const & seqVId,
        String<TraceSegment<TPosition, TSize2>, TStringSpec2 > const & traceSegments)
{
    typedef Fragment<TSize, TFragmentSpec> TFragment;

    for (TSize2 i = 0; i < length(traceSegments); ++i)
    {
        if (traceSegments[i]._traceValue == TraceBitMask::DIAGONAL)
        {
            appendValue(
                    matches,
                    TFragment(seqHId, traceSegments [i]._horizontalBeginPos, seqVId,
                              traceSegments [i]._verticalBeginPos, traceSegments [i]._length),
                    Generous());
        }
    }
}


// TODO (rmaerker): Write me!
template <typename TVertexDescriptor, typename TSpec, typename TSequence0, typename TSequence1, typename TPosition,
        typename TSize, typename TStringSpec>
inline void
adapt(String<String<TVertexDescriptor, TSpec> >& /*nodeString*/,
        TSequence0 const & /*seq0*/,
        TSequence1 const & /*seq1*/,
        String<TraceSegment<TPosition, TSize>, TStringSpec > const & /*traceSegments*/)
{
//    typedef String<TVertexDescriptor, TSpec> TVertexDescriptorString;
//    typedef typename Size<TSource>::Type TSize;
//    typedef typename Iterator<TVertexDescriptorString>::Type TStringIter;
//    TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();


     // TODO (rmaerker): see how to adapt this code here for the new structure.
//        // TraceBack values
//    TTraceValue Diagonal = 0; TTraceValue Horizontal = 1; TTraceValue Vertical = 2;
//
//    if (segLen == 0) return;
//        // Number of vertex descriptors in the first string at any position (e.g., group of 5 sequences = group of 5 vertex descriptors)
//    TSize len1 = length(getValue(getValue(str,0), 0));
//        // Number of vertex descriptors in the second string at any position (e.g., group of 5 sequences = group of 5 vertex descriptors)
//    TSize len2 = length(getValue(getValue(str,1), 0));
//
//        // Resize the node string
//    TSize index = length(nodeString);
//    resize(nodeString, index + segLen);
//
//    if (tv == Horizontal) {
//        for (int i = pos1 + segLen - 1; i>= (int) pos1;--i) {
//            resize(value(nodeString, index), len1 + len2, nilVertex);
//            TStringIter it = beSEQAN_CHECKPOINTgin(value(nodeString, index));
//            for(TPos all = 0;all<len1;++all) {
//                *it = getValue(getValue(getValue(str,0),i), all);
//                goNext(it);
//            }
//            ++index;
//        }
//    }
//    else if (tv == Vertical) {
//        for (int i = pos2 + segLen - 1; i>= (int) pos2;--i) {
//            resize(value(nodeString, index), len1 + len2, nilVertex);
//            TStringIter it = begin(value(nodeString, index));
//            it+=len1;
//            for(TPos all = 0;all<len2;++all) {
//                *it = getValue(getValue(getValue(str,1),i), all);
//                goNext(it);
//            }
//            ++index;
//        }
//    }
//    else if (tv == Diagonal) {
//        int j = pos2 + segLen - 1;
//        for (int i = pos1 + segLen - 1; i>= (int) pos1;--i) {
//            resize(value(nodeString, index), len1 + len2);
//            TStringIter it = begin(value(nodeString, index));
//            for(TPos all = 0;all<len1;++all) {
//                *it = getValue(getValue(getValue(str,0),i), all);
//                goNext(it);
//            }
//            for(TPos all = 0;all<len2;++all) {
//                *it = getValue(getValue(getValue(str,1),j), all);
//                goNext(it);
//            }
//            ++index;
//            --j;
//        }
//    }
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_ALIGN_ALIGNMENT_TRACEBACK_ADAPTOR_H_

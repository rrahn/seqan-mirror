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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Code for loading and writing FAI FASTA index files.
// ==========================================================================

// TODO(holtgrew): Include in SeqAn library.
// TODO(holtgrew): Write tests.
// TODO(holtgrew): Document.

#include <iostream>

#include <seqan/sequence.h>
#include <seqan/store.h>
#include <seqan/file.h>
#include <seqan/stream.h>

#ifndef SEQAN_CORE_APPS_RABEMA_FAI_INDEX_H_
#define SEQAN_CORE_APPS_RABEMA_FAI_INDEX_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Fai_;
typedef Tag<Fai_> Fai;

class FaiIndexEntry_
{
public:
    // Name of reference sequence.
    CharString name;
    // Number of nucleotides in sequence.
    __uint64 sequenceLength;
    // Offset in the file.
    __uint64 offset;
    // Number of sequence characters per line.
    unsigned lineLength;
    // Number of overall characters per line, including newline character(s).
    unsigned overallLineLength;

    FaiIndexEntry_() :
        sequenceLength(0), offset(0), lineLength(0), overallLineLength(0)
    {}
};

class FaiIndex
{
public:
    CharString fastaFilename;
    CharString faiFilename;

    String<FaiIndexEntry_> indexEntryStore;
    StringSet<CharString> refNameStore;
    NameStoreCache<StringSet<CharString> > refNameStoreCache;

    FaiIndex() :
        refNameStoreCache(refNameStore)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

inline void clear(FaiIndex & index)
{
    clear(index.fastaFilename);
    clear(index.faiFilename);
    clear(index.indexEntryStore);
    clear(index.refNameStore);
    refresh(index.refNameStoreCache);
}

// ----------------------------------------------------------------------------
// Function getIdByName()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Fix parameter order when getIdByName() has good parameter order.

inline bool getIdByName(FaiIndex const & index, CharString const & name, unsigned & id)
{
    return getIdByName(index.refNameStore, name, id, index.refNameStoreCache);
}

// ----------------------------------------------------------------------------
// Function sequenceLength()
// ----------------------------------------------------------------------------

inline __uint64 sequenceLength(FaiIndex const & index, unsigned refId)
{
    return index.indexEntryStore[refId].sequenceLength;
}

// ----------------------------------------------------------------------------
// Function getSequenceInfix()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline int getSequenceInfix(String<TValue, TSpec> & str, FaiIndex const & index, unsigned refId, unsigned beginPos, unsigned endPos)
{
    // TODO(holtgrew): Keep file open?
    String<char, MMap<> > mmapString;
    if (!open(mmapString, toCString(index.fastaFilename), OPEN_RDONLY))
        return 1;  // Could not open file.

    typedef typename Iterator<String<char, MMap<> >, Standard>::Type TSourceIter;
    typedef typename Iterator<String<TValue, TSpec>, Standard>::Type TTargetIter;
    TSourceIter itSource = begin(mmapString, Standard());
    __uint64 offset = index.indexEntryStore[refId].offset;
    // First, compute offset of the completely filled lines.
    unsigned numLines = beginPos / index.indexEntryStore[refId].lineLength;
    unsigned numBytes = numLines * index.indexEntryStore[refId].overallLineLength;
    // Then, compute overall offset by adding remaining bytes, too.
    numBytes += beginPos % index.indexEntryStore[refId].lineLength;
    offset += numBytes;
    // Advance iterator in MMap file.
    itSource += offset;

    // Give target string appropriate size.
    unsigned toRead = endPos - beginPos;

    // Copy out the characters from FASTA file and convert via iterator assignment to target string's type.
    resize(str, toRead, TValue());
    TTargetIter itTarget = begin(str, Standard());
    for (unsigned i = 0; i < toRead; )
    {
        if (isspace(*itSource))
        {
            ++itSource;
            continue;  // Skip spaces.
        }
        *itTarget = *itSource;
        ++itTarget;
        ++itSource;
        ++i;
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function getSequence()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline int getSequence(String<TValue, TSpec> & str, FaiIndex const & index, unsigned refId)
{
    return getSequenceInfix(str, index, refId, 0, sequenceLength(index, refId));
}

// ----------------------------------------------------------------------------
// Function load()
// ----------------------------------------------------------------------------

// Load index from file, return 0 on success and 1 on errors.  Remove any existing index data.

inline int load(FaiIndex & index, char const * fastaFilename_, char const * faiFilename_)
{
    CharString fastaFilename = fastaFilename_;
    CharString faiFilename = faiFilename_;
    clear(index);  // Also clears filename, thus backup above and restore below.
    index.fastaFilename = fastaFilename;
    index.faiFilename = faiFilename;

    // Open file.
    std::ifstream faiStream(toCString(faiFilename), std::ios::binary | std::ios::in);
    if (!faiStream.good())
        return 1;

    // Read FAI file.
    RecordReader<std::ifstream, SinglePass<> > reader(faiStream);
    CharString buffer;
    while (!atEnd(reader))
    {
        FaiIndexEntry_ entry;

        // Read REF_NAME.
        if (readUntilTabOrLineBreak(entry.name, reader) != 0)
            return 1;

        appendValue(index.refNameStore, entry.name);
        if (atEnd(reader) || value(reader) != '\t')
            return 1;  // Must be on tab.

        skipChar(reader, '\t');  // Must have been on tab, no checking.

        // Read SEQ_LENGTH.
        clear(buffer);
        if (readUntilTabOrLineBreak(buffer, reader) != 0)
            return 1;

        if (!lexicalCast2(entry.sequenceLength, buffer))
            return 1;  // Could not cast to integer.

        if (atEnd(reader) || value(reader) != '\t')
            return 1;  // Must be on tab.

        skipChar(reader, '\t');  // Must have been on tab, no checking.

        // Read OFFSET.
        clear(buffer);
        if (readUntilTabOrLineBreak(buffer, reader) != 0)
            return 1;

        if (!lexicalCast2(entry.offset, buffer))
            return 1;  // Could not cast to integer.

        if (atEnd(reader) || value(reader) != '\t')
            return 1;  // Must be on tab.

        skipChar(reader, '\t');  // Must have been on tab, no checking.

        // Read LINE_LENGTH.
        clear(buffer);
        if (readUntilTabOrLineBreak(buffer, reader) != 0)
            return 1;

        if (!lexicalCast2(entry.lineLength, buffer))
            return 1;  // Could not cast to integer.

        if (atEnd(reader) || value(reader) != '\t')
            return 1;  // Must be on tab.

        skipChar(reader, '\t');  // Must have been on tab, no checking.

        // Read OVERALL_LINE_LENGTH.
        clear(buffer);
        if (readUntilTabOrLineBreak(buffer, reader) != 0)
            return 1;

        if (!lexicalCast2(entry.overallLineLength, buffer))
            return 1;  // Could not cast to integer.

        if (!atEnd(reader) && value(reader) != '\r' && value(reader) != '\n')
            return 1;  // Must be on end of line or file.

        if (!atEnd(reader))
            skipLine(reader);  // Skip over line ending.

        appendValue(index.indexEntryStore, entry);
    }

    // Refresh name store cache.
    refresh(index.refNameStoreCache);

    return 0;
}

inline int load(FaiIndex & index, char const * fastaFilename)
{
    char buffer[1000];
    snprintf(buffer, 999, "%s.fai", toCString(fastaFilename));
    return load(index, fastaFilename, buffer);
}

inline int load(FaiIndex & index)
{
    // Cannot load if FAI filename is empty.
    if (empty(index.faiFilename))
        return 1;

    return load(index, toCString(index.fastaFilename), toCString(index.faiFilename));
}

// ---------------------------------------------------------------------------
// Function buildIndex()
// ---------------------------------------------------------------------------

/**
.Function.FaiIndex#buildIndex
..summary:Build an index file for a sequence file.
..signature:buildIndex(seqFilename[, faiFilename], Fai())
..param.seqFilename:Name of sequence file to build an index for.
...type:Shortcut.CharString
..param.faiFilename:Target name of the FAI file.
...type:Shortcut.CharString
..tag:Select the index type to build.  Must be $Fai()$.
..returns:$int$, equal to 0 on success, != 0 otherwise.
..remarks:The name of the output file will be derived from $seqFilename$.
..include:seqan/stream.h
 */

inline int buildIndex(CharString const & seqFilename, CharString const & faiFilename, Fai const & /*tag*/)
{
    // Open sequence file and create RecordReader.
    typedef String<char, MMap<> > TMMapString;
    TMMapString mmapString;
    if (!open(mmapString, toCString(seqFilename), OPEN_RDONLY))
        return 1;  // Could not open file.

    RecordReader<TMMapString, SinglePass<Mapped> > reader(mmapString);
    // Get file format, must be FASTA for FAI.
    AutoSeqStreamFormat tagSelector;
    if (!checkStreamFormat(reader, tagSelector))
        return 1;  // Invalid format.

    if (tagSelector.tagId != 1)
        return 1;  // Invalid format, not FASTA.

    // Open index files.
    std::ofstream indexOut(toCString(faiFilename), std::ios::binary | std::ios::out);

    // Re-using the FASTA/FASTQ parsing code from read_fasta_fastq is not really feasible here.  We roll our own
    // mini-parser from scratch.
    CharString line;
    CharString readName;
    __uint32 seqLength = 0;
    __uint64 seqOffset = 0;
    __uint32 lineLength = 0;
    __uint32 lineSize = 0;
    while (!atEnd(reader))
    {
        clear(line);
        clear(readName);

        if (value(reader) != '>')
            return 1;  // Must be >.

        goNext(reader);

        int res = readUntilWhitespace(readName, reader);
        if (res != 0)
            return res;  // Error reading.

        res = skipLine(reader);
        if (res != 0)
            return res;  // Error reading.

        seqOffset = reader._current - begin(reader._string, Standard());

        res = readLine(line, reader);
        if (res != 0 && res != EOF_BEFORE_SUCCESS)
            return res;  // Error reading.

        lineSize = reader._current - begin(reader._string, Standard()) - seqOffset;
        lineLength = length(line);
        seqLength = lineLength;

        while (!atEnd(reader))
        {
            char c = value(reader);
            if (c == '>')
                break;
            if (!isspace(c))
                seqLength += 1;
            goNext(reader);
        }

        indexOut << readName << '\t' << seqLength << '\t' << seqOffset << '\t'
                 << lineLength << '\t' << lineSize << '\n';
    }

    return 0;
}

inline int buildIndex(CharString const & seqFilename, Fai const & tag)
{
    CharString faiFilename(seqFilename);
    append(faiFilename, ".fai");
    return buildIndex(seqFilename, faiFilename, tag);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_APPS_RABEMA_FAI_INDEX_H_

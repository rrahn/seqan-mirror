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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de.>
// ==========================================================================
// Sample NUM reads from input file FILE.
// ==========================================================================

#include <iostream>
#include <set>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/random.h>
#include <seqan/file.h>
#include <seqan/stream.h>

const unsigned SEED = 42;

int main(int argc, char const ** argv)
{
    using namespace seqan;

    // Check command line count.
    if (argc != 3)
    {
        std::cerr << "Invalid number of arguments.\n"
                  << "USAGE: sample_seqs IN.{fasta,fastq} NUM\n";
        return 1;
    }

    // Get number of sequences to sample.
    unsigned num = atoi(argv[2]);

    // Open file.
    std::cerr << "Opening file..." << std::endl;
    MultiSeqFile multiSeqFile;
    if (!open(multiSeqFile.concat, argv[1], OPEN_RDONLY))
    {
        std::cerr << "Could not open file " << argv[1] << "\n";
        return 1;
    }

    // Guess format, split files, now, the length is available.
    std::cerr << "Splitting records..." << std::endl;
    AutoSeqFormat format;
    guessFormat(multiSeqFile.concat, format);
    split(multiSeqFile, format);

    // Sanity check on number to sample.
    if (num > length(multiSeqFile))
    {
        std::cerr << "Request to sample more reads than there actually are!" << std::endl;
        return 0;
    }
    if (2 * num > length(multiSeqFile))
    {
        std::cerr << "WARNING There are " << length(multiSeqFile)
                  << " reads, we want to sample " << num
                  << ", it make quite some time!" << std::endl;
    }

    // Now, sample reads to pick.
    Rng<MersenneTwister> rng(SEED);
    Pdf<Uniform<unsigned> > pdf(0, length(multiSeqFile) - 1);
    std::cerr << "Sampling ids..." << std::endl;
    std::set<unsigned> sampledIds;
    while (sampledIds.size() < num)
    {
        unsigned x = pickRandomNumber(rng, pdf);
        sampledIds.insert(x);
    }

    // Finally, sample reads.
    CharString id, qual;
    String<Dna5Q> seq;
    std::set<unsigned>::iterator it = sampledIds.begin();
    for(; it != sampledIds.end(); ++it)
    {
        assignSeq(seq, multiSeqFile[*it], format);
        assignQual(qual, multiSeqFile[*it], format);
        assignSeqId(id, multiSeqFile[*it], format);
        assignQualities(seq, qual);
        writeRecord(std::cout, id, seq, Fastq());
    }

    return 0;
}
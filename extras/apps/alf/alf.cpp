// ==========================================================================
//                  ALF - Alignment free sequence comparison
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: Jonathan Goeke <goeke@molgen.mpg.de>
// ==========================================================================
// Alignment free sequence comparison.
//
// This application can be used to calculate pairwise scores of DNA Sequences
// without alignments.
//
// The following scores are implemented: N2, D2, D2Star, D2z.
// ==========================================================================

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/misc/edit_environment.h>
#include <seqan/misc/misc_cmdparser.h>
#include <seqan/alignment_free.h>

using namespace seqan;
using namespace std;

int main(int argc, const char * argv[])
{
    CommandLineParser parser("alf");
    addTitleLine(parser, "******************************************");
    addTitleLine(parser, "***                 ALF                ***");
    addTitleLine(parser, "*** Alignment free sequence comparison ***");
    addTitleLine(parser, "******************************************");
    addUsageLine(parser, "-i <inputFile.fasta> -o <outputFile.txt> [Options]");

    addOption(parser, CommandLineOption('i', "inputFile", "Name of the multifasta input file",  (OptionType::String | OptionType::Mandatory),""));
    addOption(parser, CommandLineOption('o', "outputFile", "Name of the file to which the tab delimited matrix with pairwise scores will be written" , OptionType::String,""));
    addOption(parser, CommandLineOption('m', "method", "[N2, D2, D2Star, D2z]", OptionType::String));
    addOption(parser, CommandLineOption('k', "kmerSize", "[integer] Size of the k-mers", OptionType::Int));
    addOption(parser, CommandLineOption("mo", "bgModelOrder", "[integer] Order of background markov model", OptionType::Int));
    addOption(parser, CommandLineOption("rc", "reverseComplement", "['','bothStrands','mean','min','max'] N2 only. Default: only input strand. Select 'bothStrands' to score both strands simultaneously)", OptionType::String));
    addOption(parser, CommandLineOption("mm", "mismatches", "[0,1] N2 only. Default '0' (off), '1': Calculate N2 using the kmer-neighbourhood with one mismatch", OptionType::Int));
    addOption(parser, CommandLineOption("mmw", "mismatchWeight", "[Double] N2 only.  Weight of counts for words with mismatches", OptionType::Double));
    addOption(parser, CommandLineOption("kwf", "kmerWeightsFile", "N2 only. Print kmerWeights for every sequence to this file.", OptionType::String,""));
    addOption(parser, CommandLineOption('v', "verbose", "[true, false] Print details on progress to the screen", OptionType::Boolean));

    if (!parse(parser, argc, argv)) {
        return 1;
    } else if(isSetShort(parser,'h')) {
        help(parser,std::cerr);
        return 0;
    }

    // Declare all parameters
    String<char> kmerWeightsFileTmp;
    String<char> inFileTmp;
    String<char> outFileTmp;

    if (isSetShort (parser, 'i'))
        getOptionValueShort (parser,'i',inFileTmp);
    String<char, CStyle> inFile = inFileTmp;

    if (isSetShort (parser, 'o'))
        getOptionValueShort (parser,'o',outFileTmp);
    String<char, CStyle> outFile = outFileTmp;

    String<char> method = "N2";
    if (isSetShort (parser, 'm'))
        getOptionValueShort (parser,'m',method);

    int kmerSize = 4;
    if (isSetShort (parser, 'k'))
        getOptionValueShort (parser,'k',kmerSize);

    int bgModelOrder = 1;
    if (isSetShort (parser, "mo"))
        getOptionValueShort (parser,"mo",bgModelOrder);

    String<char>  revCom;
    if (isSetShort (parser, "rc"))
        getOptionValueShort (parser,"rc",revCom);

    unsigned mismatches = 0;
    if (isSetShort (parser, "mm"))
        getOptionValueShort (parser,"mm",mismatches);

    double  mismatchWeight = 0.1;
    if (isSetShort (parser, "mmw"))
        getOptionValueShort (parser,"mmw",mismatchWeight);

    String<char, CStyle> kmerWeightsFile;
    if (isSetShort (parser, "kwf"))
        getOptionValueShort (parser,"kwf",kmerWeightsFileTmp);
    kmerWeightsFile=kmerWeightsFileTmp;

    bool verbose = false;
    if (isSetShort (parser, 'v'))
        getOptionValueShort (parser,'v',verbose);

    // Definition of type DNA string sets
    typedef String<Dna5>        TText;
    typedef StringSet<TText>    TStringSet;

    // Definition of mxn two-dimensional matrix
    typedef Matrix<double, 2> TMatrix;

    TMatrix myMatrix;  // myMatrix stores pairwise kmerScores

    TStringSet mySequenceSet;  // mySequenceSet stores all sequences from the multi-fasta file
    if (inFile != "")  // read in file
    {
        MultiSeqFile multiSeqFile;
        open(multiSeqFile.concat, inFile, OPEN_RDONLY);
        AutoSeqFormat format;
        guessFormat(multiSeqFile.concat, format);
        split(multiSeqFile, format);
        unsigned seqCount = length(multiSeqFile);
        StringSet<String<Dna5Q> > seqs;
        StringSet<CharString> seqIDs;
        reserve(mySequenceSet, seqCount, Exact());
        reserve(seqIDs, seqCount, Exact());
        String<Dna5Q> seq;
        TText sequenceTmp;
        CharString qual;
        CharString id;
        for (unsigned i = 0; i < seqCount; ++i)
        {
            assignSeq(sequenceTmp, multiSeqFile[i], format);        // read sequence
            assignSeqId(id, multiSeqFile[i], format);       // read sequence id
            appendValue(mySequenceSet, sequenceTmp, Generous());
            appendValue(seqIDs, id, Generous());
        }
    }

    // Dispatch to alignment free comparisons with different scores.
    if (method == "D2")
    {
        AFScore<D2> myScoreD2(kmerSize, verbose);
        alignmentFreeComparison(myMatrix, mySequenceSet, myScoreD2);
    }
    else if (method == "D2z")
    {
        AFScore<D2z> myScoreD2z(kmerSize, bgModelOrder, verbose);
        alignmentFreeComparison(myMatrix, mySequenceSet, myScoreD2z);
    }
    else if (method == "D2Star")
    {
        AFScore<D2Star> myScoreD2Star(kmerSize, bgModelOrder, verbose);
        alignmentFreeComparison(myMatrix, mySequenceSet, myScoreD2Star);
    }
    else if (method == "N2")
    {
        AFScore<N2> myScoreN2(kmerSize, bgModelOrder, revCom, mismatches, mismatchWeight, kmerWeightsFile, verbose);
        alignmentFreeComparison(myMatrix, mySequenceSet, myScoreN2);
    }

    // Write out resulting matrix; to file if file name was given, to stdout otherwise.
    if (outFile != "")
    {
        ofstream myfile(outFile, std::ios::binary | std::ios::out);
        myfile << myMatrix;
        myfile.close();
    }
    else
    {
        std::cout << "\n" << myMatrix;
    }
    return 0;
}

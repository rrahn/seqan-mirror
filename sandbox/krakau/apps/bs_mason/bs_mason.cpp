
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

#include "simulate_illumina.h"
#include "simulate_454.h"
#include "simulate_sanger.h"

#include "simulate_illumina_bs.h"

#include "bs_mason.h"

using namespace seqan;


/* Print global help.  This function is called when the user specifies
 * a wrong simulation command.
 */
void printHelpGlobal() {
    std::cerr << "Mason - A Read Simulator" << std::endl
              << "(c) 2010 by Manuel Holtgrewe" << std::endl
              << std::endl
              << "Usage: mason illumina [OPTIONS] [SEQUENCE.fasta]" << std::endl
              << "       mason illumina_bs [OPTIONS] [SEQUENCE.fasta]" << std::endl
              << "       mason 454 [OPTIONS] [SEQUENCE.fasta]" << std::endl
              << "       mason sanger [OPTIONS] [SEQUENCE.fasta]" << std::endl
              << "       mason sanger [OPTIONS] [SEQUENCE.fasta]" << std::endl
              << std::endl
              << "Call with 'mason READS-TYPE --help' to get detailed help." << std::endl;
}

int parseOptions(Options<Global> & options, const int argc, const char * argv[]) {
    if (argc == 1 ||
        (argc >= 2 && CharString(argv[1]) == "--help") ||
        (argc >= 2 && CharString(argv[1]) == "-h") ||
        (CharString(argv[1]) != "illumina" && CharString(argv[1]) != "illumina_bs" &&   CharString(argv[1]) != "454" &&
         CharString(argv[1]) != "sanger")) {
        printHelpGlobal();
        return 1;
    }

    if (CharString(argv[1]) == "illumina") {
        options.readsType = READS_TYPE_ILLUMINA;
    } else if (CharString(argv[1]) == "illumina_bs") {
        options.readsType = READS_TYPE_ILLUMINA_BS;
    } else if (CharString(argv[1]) == "454") {
        options.readsType = READS_TYPE_454;
    } else if (CharString(argv[1]) == "sanger") {
        options.readsType = READS_TYPE_SANGER;
    } else {
        printHelpGlobal();
        return 1;
    }

    return 0;
}


int main(const int argc, const char * argv[]) {
    // Switch command type (which reads to simulate) or show global help.
    Options<Global> globalOptions;
    int ret = parseOptions(globalOptions, argc, argv);
    if (ret != 0)
        return ret;

    // Kick off read simulation, depending on the chosen command.  We
    // use a type-based compile-time dispatch for selecting the
    // appropriate function for simulation.  This if/then/else switch
    // selects the right code path.
    if (globalOptions.readsType == READS_TYPE_ILLUMINA) {
        CommandLineParser parser;
        setUpCommandLineParser(parser, IlluminaReads());
        Options<IlluminaReads> options;
        CharString referenceFilename;
        int ret = parseCommandLineAndCheck(options, referenceFilename, parser, argc, argv);
        if (options.showHelp)
            return 0;
        if (ret != 0)
            return ret;
        return simulateReads(options, referenceFilename, IlluminaReads());
    } else if (globalOptions.readsType == READS_TYPE_ILLUMINA_BS) {     // bs_change:
        CommandLineParser parser;
        setUpCommandLineParser(parser, IlluminaReadsBS());
        Options<IlluminaReadsBS> options;
        CharString referenceFilename;
        int ret = parseCommandLineAndCheck(options, referenceFilename, parser, argc, argv);
        if (options.showHelp)
            return 0;
        if (ret != 0)
            return ret;
        // Bs_change
        return simulateReads(options, referenceFilename, IlluminaReadsBS());
    } else if (globalOptions.readsType == READS_TYPE_454) {
        CommandLineParser parser;
        setUpCommandLineParser(parser, LS454Reads());
        Options<LS454Reads> options;
        CharString referenceFilename;
        int ret = parseCommandLineAndCheck(options, referenceFilename, parser, argc, argv);
        if (options.showHelp)
            return 0;
        if (ret != 0)
            return ret;
        return simulateReads(options, referenceFilename, LS454Reads());
    } else if (globalOptions.readsType == READS_TYPE_SANGER) {
        CommandLineParser parser;
        setUpCommandLineParser(parser, SangerReads());
        Options<SangerReads> options;
        CharString referenceFilename;
        int ret = parseCommandLineAndCheck(options, referenceFilename, parser, argc, argv);
        if (options.showHelp)
            return 0;
        if (ret != 0)
            return ret;
        return simulateReads(options, referenceFilename, SangerReads());
    } else {
        SEQAN_ASSERT_FAIL("Invalid reads type!");
    }

    return 0;
}

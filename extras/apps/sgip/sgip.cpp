// ===========================================================================
//                 SGIP - Solution of Graph Isomorphism Problem
// ===========================================================================
// Copyright (C) 2012 by Jialu Hu
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your options) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ===========================================================================
// Author: Jialu Hu <Jialu.Hu@fu-berlin.de>
// ===========================================================================
// This application is used to determine whether two given graphs are
// isomorphic or not by using heuristic approach.
// ===========================================================================

#include <time.h>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/graph_types.h>
#include <seqan/misc/misc_cmdparser.h>

#include "sgip.h"
#include "sgip_base.h"
#include "sgip_output.h"

using namespace seqan;

// ==========================================================================
// Forwards
// ==========================================================================

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// --------------------------------------------------------------------------
// Enum FileOption
// --------------------------------------------------------------------------

enum FileOption
{
    ZERO,                  // graph given through parameter
    FIRST,                 // graph created from  orignalFile
    SECOND                 // graph created from  comparisive file
};

// --------------------------------------------------------------------------
// Enum SearchingType
// --------------------------------------------------------------------------

enum SearchingType
{
    HEURISTIC,            // heuristic searching approach
    BRUTEFORTH
};

// --------------------------------------------------------------------------
// Enum SgipOption
// --------------------------------------------------------------------------

// Option for sgip.
struct SgipOption
{
    // I/O options.
    CharString orginalFile;        // name of orginal file(first graph)
    CharString comparFile;         // name of comparisive file(second graph)
    CharString outputFile;         // name of result file
    CharString compareFolder;

    // More options.
    bool autoMetric;               // search Metric dimension if true
    FileOption activeFile;
    SearchingType searchingType;
    CharString algorithm;          // Search strategy for metric dimension,e.g. Greedy, genetic etc.
    unsigned odimension;           // metric dimension of orginal graph specified by user
    unsigned cdimension;           // metric dimension of comparitive graph specified by user
    bool showHelp;
    bool showVersion;
    bool isoCheck;                //to check whether two input graphs are isomorphic
    bool isPrintFile;             //print result if outputFile is specified
    bool isAllatOnce;
    int verbose;

    SgipOption()
    {
        algorithm = "greedy";
        activeFile = FIRST;
        searchingType = HEURISTIC;
        showHelp = 0;
        showVersion = 0;
        isoCheck = 0;
        isPrintFile = false;
        odimension = 3;
        cdimension = 3;
        verbose = 0;
        isAllatOnce = false;
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function _sgip()
// --------------------------------------------------------------------------

// Test sgip with Options.
template <typename TOption>
int _sgip(TOption & options)
{
    typedef Graph<Directed<> > TGraph;
    typedef String<bool>       TMat;
    TMat leastmat1, leastmat2;
    TGraph g1, g2;
    char const * file1, * file2;
    if (!options.isoCheck)
    {
        file1 = toCString(options.orginalFile);
        if (!_createGraph(g1, SivaLab(), file1))
            return 1;

        getCanonicalLabel(leastmat1, g1);
        if (options.verbose > 1)
            outputLabel(leastmat1, FFFF());
    }
    else
    {
        file1 = toCString(options.orginalFile);
        file2 = toCString(options.comparFile);
        if (!_createGraph(g1, SivaLab(), file1))
            return 1;

        if (!_createGraph(g2, SivaLab(), file2))
            return 1;

        if (checkIsomorphic(g1, g2))
        {
            if (options.verbose > 0)
                std::cout << "They are isomorphic!" << std::endl;
        }
        else
        {
            if (options.verbose > 0)
                std::cout << "They are not isomorphic!" << std::endl;
        }
    }
    return 0;
}

// --------------------------------------------------------------------------
// Function setupParser()
// --------------------------------------------------------------------------

// Set up parser.
template <typename TParser>
void _setupParser(TParser & parser)
{
    addVersionLine(parser, "Version 1.1 (Nov. 16th 2011) SeqAn Revision:");

    addTitleLine(parser, "***********************************************");
    addTitleLine(parser, "* SGIP - Solution of Graph Isomorphism Problem*");
    addTitleLine(parser, "*      (c) Copyright 2011 by Jialu HU         *");
    addTitleLine(parser, "*       email<Jialu.Hu@fu-berlin.de>          *");
    addTitleLine(parser, "***********************************************");
    addTitleLine(parser, "");

    addUsageLine(parser, "-o <orginal graph> [Option]");
    addSection(parser, "Mandatory Options");
    addOption(parser, CommandLineOption("o", "orginal", "File containing orginal graph", OptionType::String | OptionType::Mandatory));
    addSection(parser, "Main Options");
    addOption(parser, CommandLineOption("a", "algorithm", "Algorithm used for searching metric dimension", OptionType::String, "greedy"));
    addOption(parser, CommandLineOption("s", "searching", "Searching algorithm used for detecting resolving set,heuristic 0 bruteforce 1.", OptionType::Int, 0));
    addOption(parser, CommandLineOption("i", "isomorphism", "To check whether two given graphs are isomorphic", OptionType::Bool));
    addOption(parser, CommandLineOption("c", "comparitive", "File containing comparitive graph", OptionType::String));
    addOption(parser, CommandLineOption("od", "odimension", "Specified initial dimension of orginal graph by user", OptionType::Int, 3));
    addOption(parser, CommandLineOption("cd", "cdimension", "Specified initial dimension of comparitive graph by user", OptionType::Int, 3));
    addOption(parser, CommandLineOption("ad", "directory", "test for all graphs in a specified directory", OptionType::String));
    addOption(parser, CommandLineOption('v', "verbose", "control the level of output files ", OptionType::Int, 0));

    addSection(parser, "Output Option");
    addOption(parser, CommandLineOption("r", "result", "File containing canonical label of input graphs", OptionType::String));
}

// --------------------------------------------------------------------------
// Function _parseOptions()
// --------------------------------------------------------------------------

template <typename TOption, typename TParser>
int _parseOptions(TOption & options, TParser & parser)
{
    getOptionValueShort(parser, 'o', options.orginalFile);
    if (isSetShort(parser, 'i'))
    {
        if (!isSetShort(parser, 'c'))
        {
            if (!isSetShort(parser, "ad"))
            {
                std::cerr << "sgip" << ":comparitive file has not been specified!" << std::endl;
                shortHelp(parser, std::cerr);
                return 0;
            }
        }
        options.isoCheck = true;
        getOptionValueShort(parser, 'c', options.comparFile);
    }
    if (isSetShort(parser, 'r'))
    {
        options.isPrintFile = true;
        getOptionValueShort(parser, 'r', options.outputFile);
    }
    if (isSetShort(parser, "od"))
        getOptionValueShort(parser, "od", options.odimension);
    if (isSetShort(parser, "cd"))
        getOptionValueShort(parser, "cd", options.cdimension);
    if (isSetShort(parser, 'a'))
        getOptionValueShort(parser, 'a', options.algorithm);
    if (isSetShort(parser, 'v'))
        getOptionValueShort(parser, 'v', options.verbose);
    if (isSetShort(parser, 's'))
    {
        int s;
        getOptionValueShort(parser, 's', s);
        options.searchingType = static_cast<SearchingType>(s);
    }
    _printTitle(parser, std::cerr);
    return 1;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.
int main(int argc, char const ** argv)
{
    // Setup command line parser.
    CharString        appName("sgip");
    CommandLineParser parser(appName);
    SgipOption options;
    time_t start, end;
    double dif;

    time(&start);
    _setupParser(parser);
    // Then, parse the command line and handle the cases where help display
    // is requested or erroneous parameters were given.
    if (!parse(parser, argc, argv))
    {
        if (isSetShort(parser, 'h') || isSetShort(parser, 'v'))
            return 0;

        shortHelp(parser, std::cerr);
        return 1;
    }
    if (!_parseOptions(options, parser))
        return 1;
    // Finally, launch the program.
    int ret = 0;
    if (!options.isAllatOnce)
        ret = _sgip(options);
    time(&end);
    dif = difftime(end, start);
    if (options.verbose > 2)
        std::cout << std::setprecision(10) << "escaped time:" << dif << "seconds" << std::endl;
    return ret;
}

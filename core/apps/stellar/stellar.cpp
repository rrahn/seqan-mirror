// ==========================================================================
//                    STELLAR - SwifT Exact LocaL AligneR
//                   http://www.seqan.de/projects/stellar/
// ==========================================================================
// Copyright (C) 2010-2012 by Birte Kehr
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
//
// ==========================================================================
// Author: Birte Kehr <birte.kehr@fu-berlin.de>
// ==========================================================================

#define SEQAN_PROFILE

#include <seqan/index.h>
#include <seqan/misc/misc_cmdparser.h>
#include "stellar.h"
#include "stellar_output.h"

using namespace seqan;

///////////////////////////////////////////////////////////////////////////////
// Initializes a Finder object for a database sequence,
//  calls stellar, and writes matches to file
template<typename TSequence, typename TId, typename TPattern, typename TMatches>
inline bool
_stellarOnOne(TSequence & database,
			  TId & databaseID,
			  TPattern & swiftPattern,
			  bool databaseStrand,
			  TMatches & matches,
			  StellarOptions & options) {
	std::cout << "  " << databaseID;
	if (!databaseStrand) std::cout << ", complement";
  std::cout << std::flush;

	// finder
    typedef Finder<TSequence, Swift<SwiftLocal> > TFinder;
	TFinder swiftFinder(database, options.minRepeatLength, options.maxRepeatPeriod);

	// stellar
	if (options.fastOption == CharString("exact"))
		stellar(swiftFinder, swiftPattern, options.epsilon, options.minLength, options.xDrop, 
				options.disableThresh, options.compactThresh, options.numMatches, options.verbose,
				databaseID, databaseStrand, matches, AllLocal());
	else if (options.fastOption == "bestLocal")
		stellar(swiftFinder, swiftPattern, options.epsilon, options.minLength, options.xDrop, 
				options.disableThresh, options.compactThresh, options.numMatches, options.verbose,
				databaseID, databaseStrand, matches, BestLocal());
	else if (options.fastOption == "bandedGlobal")
		stellar(swiftFinder, swiftPattern, options.epsilon, options.minLength, options.xDrop, 
				options.disableThresh, options.compactThresh, options.numMatches, options.verbose,
				databaseID, databaseStrand, matches, BandedGlobal());
	else if (options.fastOption == "bandedGlobalExtend")
		stellar(swiftFinder, swiftPattern, options.epsilon, options.minLength, options.xDrop, 
				options.disableThresh, options.compactThresh, options.numMatches, options.verbose,
				databaseID, databaseStrand, matches, BandedGlobalExtend());
	else {
		std::cerr << "\nUnknown verification strategy: " << options.fastOption << std::endl;
		return false;
	}

	std::cout << std::endl;
	return true;
}

//////////////////////////////////////////////////////////////////////////////
namespace SEQAN_NAMESPACE_MAIN
{

template <typename TStringSet, typename TShape, typename TSpec>
struct Cargo<Index<TStringSet, IndexQGram<TShape, TSpec> > > {
	typedef struct {
		double		abundanceCut;
	} Type;
};

//////////////////////////////////////////////////////////////////////////////
// Repeat masker
template <typename TStringSet, typename TShape, typename TSpec>
inline bool _qgramDisableBuckets(Index<TStringSet, IndexQGram<TShape, TSpec> > &index) 
{
	typedef Index<TStringSet, IndexQGram<TShape, TSpec> >	TReadIndex;
	typedef typename Fibre<TStringSet, QGramDir>::Type		TDir;
	typedef typename Iterator<TDir, Standard>::Type			TDirIterator;
	typedef typename Value<TDir>::Type						TSize;

	TDir &dir    = indexDir(index);
	bool result  = false;
	unsigned counter = 0;
	TSize thresh = (TSize)(length(index) * cargo(index).abundanceCut);
	if (thresh < 100) thresh = 100;

	TDirIterator it = begin(dir, Standard());
	TDirIterator itEnd = end(dir, Standard());
	for (; it != itEnd; ++it)
		if (*it > thresh) 
		{
			*it = (TSize)-1;
			result = true;
			++counter;
		}

	if (counter > 0)
		std::cerr << "Removed " << counter << " k-mers" << ::std::endl;

	return result;
}

}

///////////////////////////////////////////////////////////////////////////////
// Initializes a Pattern object with the query sequences, 
//  and calls _stellarOnOne for each database sequence
template<typename TSequence, typename TId>
inline bool
_stellarOnAll(StringSet<TSequence> & databases,
			  StringSet<TId> & databaseIDs,
			  StringSet<TSequence> & queries,
			  StringSet<TId> & queryIDs,
			  StellarOptions & options) {
    // pattern
    typedef Index<StringSet<TSequence, Dependent<> >, IndexQGram<SimpleShape, OpenAddressing> > TQGramIndex;
    TQGramIndex qgramIndex(queries);
    resize(indexShape(qgramIndex), options.qGram);
	cargo(qgramIndex).abundanceCut = options.qgramAbundanceCut;
	Pattern<TQGramIndex, Swift<SwiftLocal> > swiftPattern(qgramIndex);
	
	if (options.verbose) swiftPattern.params.printDots = true;

	// Construct index
	std::cout << "Constructing index..." << std::endl;
	indexRequire(qgramIndex, QGramSADir());
	std::cout << std::endl;

    // container for eps-matches
	StringSet<QueryMatches<StellarMatch<TSequence, TId> > > matches;
    resize(matches, length(queries));

	std::cout << "Aligning all query sequences to database sequence..." << std::endl;
	for(unsigned i = 0; i < length(databases); ++i) {
		// positive database strand
		if (options.forward) {
			if (!_stellarOnOne(databases[i], databaseIDs[i], swiftPattern, true, matches, options))
				return 1;
		}
		// negative (reverse complemented) database strand
		if (options.reverse) {
			reverseComplement(databases[i]);
			if (!_stellarOnOne(databases[i], databaseIDs[i], swiftPattern, false, matches, options))
				return 1;
			reverseComplement(databases[i]);
		}
	}
	std::cout << std::endl;
	
	// file output
	if (options.disableThresh != (unsigned)-1) {
		if (!_outputMatches(matches, queries, queryIDs, databases, options.verbose,
		                      options.outputFile, options.outputFormat, options.disabledQueriesFile)) return 1;
	}
	else {
		if (!_outputMatches(matches, queryIDs, databases, options.verbose,
		                      options.outputFile, options.outputFormat)) return 1;
	}
	
	return 0;
}

template<typename TId>
inline bool
_checkUniqueId(std::set<TId> & uniqueIds, TId & id)
{
	TId shortId;
	typedef typename Iterator<TId>::Type TIterator;

	TIterator it = begin(id);
	TIterator itEnd = end(id);
	while (it != itEnd && *it > 32)
	{
		appendValue(shortId, *it);
		++it;
	}
	
	if (uniqueIds.count(shortId) == 0)
	{
		uniqueIds.insert(shortId);
		return 1;
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Imports sequences from a file, 
//  stores them in the StringSet seqs and their identifiers in the StringSet ids
template<typename TSequence, typename TId>
inline bool
_importSequences(CharString const & fileName,
				 CharString const & name,
                 StringSet<TSequence> & seqs,
                 StringSet<TId> & ids)
{
    MultiSeqFile multiSeqFile;
	if (!open(multiSeqFile.concat, toCString(fileName), OPEN_RDONLY))
	{
		std::cerr << "Failed to open " << name << " file." << std::endl;
        return false;
	}

    AutoSeqFormat format;
    guessFormat(multiSeqFile.concat, format);
    split(multiSeqFile, format);
    
    unsigned seqCount = length(multiSeqFile);
    reserve(seqs, seqCount, Exact());
    reserve(ids, seqCount, Exact());

	std::set<TId> uniqueIds; // set of short IDs (cut at first whitespace)
	bool idsUnique = true;

    TSequence seq;
    TId id;
    for(unsigned i = 0; i < seqCount; ++i)
	{
        assignSeq(seq, multiSeqFile[i], format);
        assignSeqId(id, multiSeqFile[i], format);

		idsUnique &= _checkUniqueId(uniqueIds, id);

        appendValue(seqs, seq, Generous());
        appendValue(ids, id, Generous());
    }

	std::cout << "Loaded " << seqCount << " " << name << " sequence" << ((seqCount>1)?"s.":".") << std::endl;
	if (!idsUnique) std::cerr << "WARNING: Non-unique " << name << " ids. Output can be ambigous.\n";
    return true;
}

///////////////////////////////////////////////////////////////////////////////
// Calculates parameters from parameters in options object and from sequences and writes them to std::cout
template<typename TStringSet>
void _writeMoreCalculatedParams(StellarOptions & options, TStringSet & databases, TStringSet & queries) {
//IOREV _notio_
	typedef typename Size<TStringSet>::Type TSize;

	if (options.qgramAbundanceCut != 1) {
		std::cout << "Calculated parameters:" << std::endl;
	}

	TSize queryLength = length(concat(queries));
	if (options.qgramAbundanceCut != 1) {
		std::cout << "  q-gram expected abundance : ";
		std::cout << queryLength/(double)((long)1<<(options.qGram<<1)) << std::endl;
		std::cout << "  q-gram abundance threshold: ";
		std::cout << _max(100,(int)(queryLength*options.qgramAbundanceCut)) << std::endl;
		std::cout << std::endl;
	}

	// Computation of maximal E-value for this search

	TSize maxLengthQueries = 0;
	TSize maxLengthDatabases = 0;

	typename Iterator<TStringSet>::Type dbIt = begin(databases);
	typename Iterator<TStringSet>::Type dbEnd = end(databases);
	while (dbIt != dbEnd) {
		if (length(*dbIt) > maxLengthDatabases) {
			maxLengthDatabases = length(*dbIt);
		}
		++dbIt;
	}

	typename Iterator<TStringSet>::Type queriesIt = begin(queries);
	typename Iterator<TStringSet>::Type queriesEnd = end(queries);
	while (queriesIt != queriesEnd) {
		if (length(*queriesIt) > maxLengthQueries) {
			maxLengthQueries = length(*queriesIt);
		}
		++queriesIt;
	}

	TSize errors = static_cast<TSize>(options.minLength * options.epsilon);
	TSize minScore = options.minLength - 3*errors; // #matches - 2*#errors // #matches = minLenght - errors, 

	std::cout << "All matches matches resulting from your search have an E-value of: " << std::endl;
	std::cout << "        " << _computeEValue(minScore, maxLengthQueries, maxLengthDatabases) << " or smaller";
	std::cout << "  (match score = 1, error penalty = -2)" << std::endl;

	std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// Calculates parameters from parameters in options object and writes them to std::cout
void _writeCalculatedParams(StellarOptions & options) {
//IOREV _notio_
	int errMinLen = (int) floor(options.epsilon * options.minLength);
	int n = (int) ceil((errMinLen + 1) / options.epsilon);
	int errN = (int) floor(options.epsilon * n);
	unsigned smin = (unsigned) _min(ceil((double)(options.minLength-errMinLen)/(errMinLen+1)),
		                            ceil((double)(n-errN)/(errN+1)));

	std::cout << "Calculated parameters:" << std::endl;
	if (options.qGram == (unsigned)-1) {
		options.qGram = (unsigned)_min(smin, 32u);
		std::cout << "  k-mer length: " << options.qGram << std::endl;
	}

	int threshold = (int) _max(1, (int) _min((n + 1) - options.qGram * (errN + 1),
											 (options.minLength + 1) - options.qGram * (errMinLen + 1)));
	int overlap = (int) floor((2 * threshold + options.qGram - 3) / (1 / options.epsilon - options.qGram));
	int distanceCut = (threshold - 1) + options.qGram * overlap + options.qGram;
	int logDelta = _max(4, (int) ceil(log((double)overlap + 1) / log(2.0)));
	int delta = 1 << logDelta;

	std::cout << "  s^min       : " << smin << std::endl;
	std::cout << "  threshold   : " << threshold << std::endl;
	std::cout << "  distance cut: " << distanceCut << std::endl;
	std::cout << "  delta       : " << delta << std::endl;
	std::cout << "  overlap     : " << overlap << std::endl;
	std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// Writes user specified parameters from options object to std::cout
template<typename TOptions>
void
_writeSpecifiedParams(TOptions & options) {
//IOREV _notio_
	// Output user specified parameters
	std::cout << "User specified parameters:" << std::endl;
	std::cout << "  minimal match length             : " << options.minLength << std::endl;
	std::cout << "  maximal error rate (epsilon)     : " << options.epsilon << std::endl;
	std::cout << "  maximal x-drop                   : " << options.xDrop << std::endl;
	if (options.qGram != (unsigned)-1)
		std::cout << "  k-mer (q-gram) length            : " << options.qGram << std::endl;
	std::cout << "  search forward strand            : " << ((options.forward)?"yes":"no") << std::endl;
	std::cout << "  search reverse complement        : " << ((options.reverse)?"yes":"no") << std::endl;
	std::cout << std::endl;

	std::cout << "  verification strategy            : " << options.fastOption << std::endl;
	if (options.disableThresh != (unsigned)-1) {
		std::cout << "  disable queries with more than   : " << options.disableThresh << " matches" << std::endl;
	}
	std::cout << "  maximal number of matches        : " << options.numMatches << std::endl;
	std::cout << "  duplicate removal every          : " << options.compactThresh << std::endl;
	if (options.maxRepeatPeriod != 1 || options.minRepeatLength != 1000) {
		std::cout << "  max low complexity repeat period : " << options.maxRepeatPeriod << std::endl;
		std::cout << "  min low complexity repeat length : " << options.minRepeatLength << std::endl;
	}
	if (options.qgramAbundanceCut != 1) {
		std::cout << "  q-gram abundance cut ratio       : " << options.qgramAbundanceCut << std::endl;
	}
	std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// Writes file name from options object to std::cout
template<typename TOptions>
void
_writeFileNames(TOptions & options) {
//IOREV _notio_
	std::cout << "Database file   : " << options.databaseFile << std::endl;
	std::cout << "Query file      : " << options.queryFile << std::endl;
	std::cout << "Output file     : " << options.outputFile << std::endl;
	std::cout << "Output format   : " << options.outputFormat << std::endl;
	if (options.disableThresh != (unsigned)-1) {
		std::cout << "Disabled queries: " << options.disabledQueriesFile << std::endl;
	}
	std::cout << std::endl;
}

inline void
_addVersion(CommandLineParser& parser) {
	::std::string rev = "$Revision: 11692 $";
	addVersionLine(parser, "Version 1.2 (April 20th 2012) SeqAn Revision: " + rev.substr(11, 5) + "");
}

///////////////////////////////////////////////////////////////////////////////
// Parses options from command line parser and writes them into options object
template<typename TParser, typename TOptions>
bool
_parseOptions(TParser & parser, TOptions & options) {
//IOREV _notio_
    // i/o options
	getOptionValueShort(parser, 'd', options.databaseFile);
    getOptionValueShort(parser, 'q', options.queryFile);
    if (isSetShort(parser, 'o')) getOptionValueShort(parser, 'o', options.outputFile);
    if (isSetShort(parser, "od")) getOptionValueShort(parser, "od", options.disabledQueriesFile);
	if (isSetShort(parser, "of")) getOptionValueShort(parser, "of", options.outputFormat);
    if (isSetLong(parser, "no-rt")) options.noRT = true;

	// main options
	if (isSetLong(parser, "kmer")) getOptionValueLong(parser, "kmer", options.qGram);
    if (isSetLong(parser, "minLength")) getOptionValueLong(parser, "minLength", options.minLength);
	if (isSetShort(parser, 'e')) getOptionValueShort(parser, 'e', options.epsilon);
    if (isSetShort(parser, 'x')) getOptionValueShort(parser, 'x', options.xDrop);

	if (isSetShort(parser, 'f')) if (!isSetShort(parser, 'r')) options.reverse = false;
	if (isSetShort(parser, 'r')) if (!isSetShort(parser, 'f')) options.forward = false;

	if (isSetShort(parser, "vs")) getOptionValueShort(parser, "vs", options.fastOption);
	if (isSetShort(parser, "dt")) getOptionValueShort(parser, "dt", options.disableThresh);
	if (isSetShort(parser, 'n')) getOptionValueShort(parser, 'n', options.numMatches);
	if (isSetShort(parser, 's')) getOptionValueShort(parser, 's', options.compactThresh);
	if (isSetShort(parser, "rp")) getOptionValueShort(parser, "rp", options.maxRepeatPeriod);
	if (isSetShort(parser, "rl")) getOptionValueShort(parser, "rl", options.minRepeatLength);
	if (isSetShort(parser, 'a')) getOptionValueShort(parser, 'a', options.qgramAbundanceCut);

	if (isSetShort(parser, 'v')) options.verbose = 1;

	if (options.outputFormat != "gff" && options.outputFormat != "text") {
		std::cerr << "Invalid parameter value: Unknown output format." << std::endl;
		return 0;
	}

	if (options.fastOption != "exact" && options.fastOption != "bestLocal"
		 && options.fastOption != "bandedGlobal") {
		std::cerr << "Invalid parameter value: Unknown verification strategy." << std::endl;
		return 0;
	}

	if (isSetShort(parser, 'k') && options.qGram < 1) {
		std::cerr << "Invalid parameter value: Please choose a greater k-mer length." << std::endl;
		return 0;
	}

	if (isSetShort(parser, 'k') && options.qGram > 32) {
		std::cerr << "Invalid parameter value: Please choose a smaller k-mer length." << std::endl;
		return 0;
	}

	if (options.epsilon < 0.0) {
		std::cerr << "Invalid parameter value: Please choose a greater error rate." << std::endl;
		return 0;
	}

	if (options.epsilon > 0.25) {
		std::cerr << "Invalid parameter value: Please choose a smaller error rate." << std::endl;
		return 0;
	}
	if (isSetShort(parser, 'k') && options.qGram >= 1/options.epsilon) {
		std::cerr << "Invalid parameter value: Please choose q-gram length lower than 1/epsilon." << std::endl; 
		return 0;
	}

	if (options.qgramAbundanceCut > 1 || options.qgramAbundanceCut < 0) {
		std::cerr << "Invalid parameter value: Please choose a k-mer overabundance cut ration between 0 and 1.\n";
		return 0;
	}

	if (options.numMatches > options.compactThresh) {
		std::cerr << "Invalid parameter values: Please choose numMatches <= sortThresh." << std::endl;
		return 0;
	}
	return 1;
}

///////////////////////////////////////////////////////////////////////////////
// Set-Up of Command Line Parser
template<typename TParser>
void
_setParser(TParser & parser) {
	_addVersion(parser);

    addTitleLine(parser, "*******************************************");
	addTitleLine(parser, "* STELLAR - the SwifT Exact LocaL AligneR *");
	addTitleLine(parser, "* (c) Copyright 2010-2012 by Birte Kehr   *");
	addTitleLine(parser, "*******************************************");

	addUsageLine(parser, "-d <FASTA sequence file> -q <FASTA sequence file> [Options]");

	addLine(parser, "");
    addLine(parser, "An implementation of the SWIFT filter algorithm (Rasmussen et al., 2006)");
    addLine(parser, "and subsequent verification of the SWIFT hits using local alignment,");
    addLine(parser, "gapped X-drop extension, and extraction of the longest epsilon-match.");

	addSection(parser, "Non-optional Arguments:");
    addOption(parser, CommandLineOption('d', "database", "Fasta file containing the database sequences",
              (OptionType::String | OptionType::Mandatory)));
    addOption(parser, CommandLineOption('q', "query", "Fasta file containing the query sequences", 
              (OptionType::String | OptionType::Mandatory)));

	addSection(parser, "Main Options:");
    addOption(parser, CommandLineOption('e', "epsilon", "Maximal error rate (max 0.25)", OptionType::Double, "0.05"));
    addOption(parser, CommandLineOption('l', "minLength", "Minimal length of epsilon-matches", OptionType::Int, 100));
	addOption(parser, CommandLineOption('f', "forward", "Search only in forward strand of database",
		OptionType::Boolean, "both"));
	addOption(parser, CommandLineOption('r', "reverse", "Search only in reverse complement of database",
		OptionType::Boolean, "both"));
    addOption(parser, CommandLineOption('v', "verbose", "Verbosity mode.", OptionType::Bool, "false"));
    
	addSection(parser, "Filtering Options:");
    addOption(parser, CommandLineOption('k', "kmer", "Length of the q-grams (max 32)", OptionType::Int, "smin"));
    addOption(parser, CommandLineOption("rp", "repeatPeriod",
		"Maximal period of low complexity repeats to be filtered", OptionType::Int, 1));
    addOption(parser, CommandLineOption("rl", "repeatLength",
		"Minimal length of low complexity repeats to be filtered", OptionType::Int, 1000));
    addOption(parser, CommandLineOption('a', "abundanceCut",
		"k-mer overabundance cut ratio", OptionType::Double, "1"));

	addSection(parser, "Verification Options:");
    addOption(parser, CommandLineOption('x', "xDrop", "Maximal x-drop for extension", OptionType::Double, 5));
	addOption(parser, CommandLineOption("vs", "verification", "Verification strategy", OptionType::String, "exact"));
	addHelpLine(parser, "exact        = compute and extend all local alignments in SWIFT hits");
	addHelpLine(parser, "bestLocal    = compute and extend only best local alignment in SWIFT hits");
	addHelpLine(parser, "bandedGlobal = banded global alignment on SWIFT hits");
	addOption(parser, CommandLineOption("dt", "disableThresh",
		"Maximal number of verified matches before disabling verification", OptionType::Int));
	addHelpLine(parser, "for one query sequence (default infinity)");
	addOption(parser, CommandLineOption('n', "numMatches",
		"Maximal number of kept matches per query and database", OptionType::Int, 50));
	addHelpLine(parser, "If there are more matches, only the longest ones are kept.");
	addOption(parser, CommandLineOption('s', "sortThresh",
		"Number of matches triggering removal of duplicates", OptionType::Int, 500));
	addHelpLine(parser, "Choose a smaller value for saving space.");

	addSection(parser, "Output Options:");
    addOption(parser, CommandLineOption('o', "out", "Name of output file", OptionType::String, "stellar.gff"));
	addOption(parser, CommandLineOption("of", "outFormat", "Output format", OptionType::String, "gff"));
	addHelpLine(parser, "Possible formats: gff, text");
	addOption(parser, CommandLineOption("od", "outDisabled",
		"Name of output file containing disabled query sequences", OptionType::String));
	addHelpLine(parser, "(default stellar.disabled.fasta)");
    addOption(parser, CommandLineOption("t", "no-rt", "Suppress printing running time.", OptionType::Bool | OptionType::Hidden, false));
}

// TODO(holtgrew): Move this into a SeqAn misc module.

/*
.Class.ScientificNotationExponentOutputNormalizer
..summary:RAII class to normalize the output of exponents in scientific notation.
..signature:ScientificNotationExponentOutputNormalizer
..remarks:
The C/C++ standard only specifies that exponents in scientific notation have to be printed with at least two characters.
MSVC prints three and thus the output of otherwise identical programs differs.
This RAII class fixes the problem.
..example.text:Use this class in your $main()$ function, for example.
..example.code:
int main()
{
	// Make sure doubles are printed consistently on all platforms.
	ScientificNotationExponentOutputNormalizer normalizer;

	std::cout << 1e-8 << std::endl;
    return 0;
}
*/

class ScientificNotationExponentOutputNormalizer
{
public:
	unsigned _oldExponentFormat;

	ScientificNotationExponentOutputNormalizer() : _oldExponentFormat(0)
	{
#ifdef PLATFORM_WINDOWS_VS
		// Set scientific format to print two places.
		unsigned _oldExponentFormat = _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
	}

	~ScientificNotationExponentOutputNormalizer()
	{
#ifdef PLATFORM_WINDOWS_VS
		// Enable old exponent format.
		_set_output_format(_oldExponentFormat);
#endif
	}
};

///////////////////////////////////////////////////////////////////////////////
// Parses and outputs parameters, calls _stellarOnAll
int main(int argc, const char *argv[]) {
	// Makes sure that printing doubles in scientific notation is normalized on all platforms.
	ScientificNotationExponentOutputNormalizer scientificNotationNormalizer; 

    // command line parsing
    CommandLineParser parser("stellar");

    _setParser(parser);
    if (!parse(parser, argc, argv)) {
		if (isSetShort(parser, 'h') || isSetShort(parser, 'v')) return 0; 
        shortHelp(parser, std::cerr);
        return 1;
    }

	StellarOptions options = StellarOptions();
	if (!_parseOptions(parser, options)) {
		return 1;
	}

	typedef String<Dna5> TSequence;

	// output header
	_printTitle(parser, std::cout);
	std::cout << std::endl;

	// output file names
	_writeFileNames(options);

	// output parameters
	_writeSpecifiedParams(options);
	_writeCalculatedParams(options);

    // import query sequences
    StringSet<TSequence> queries;
    StringSet<CharString> queryIDs;
	if (!_importSequences(options.queryFile, "query", queries, queryIDs)) return 1;

    // import database sequence
    StringSet<TSequence > databases;
    StringSet<CharString> databaseIDs;
	if (!_importSequences(options.databaseFile, "database", databases, databaseIDs)) return 1;

    std::cout << std::endl;
	_writeMoreCalculatedParams(options, databases, queries);

    // open output files
    std::ofstream file;
    file.open(toCString(options.outputFile));
	if (!file.is_open()) {
		std::cerr << "Could not open output file." << std::endl;
		return 1;
	}
    file.close();

	if(options.disableThresh != (unsigned)-1) {
		std::ofstream daFile;
		daFile.open(toCString(options.disabledQueriesFile));
		if (!daFile.is_open()) {
			std::cerr << "Could not open file for disabled queries." << std::endl;
			return 1;
		}
		daFile.close();
	}

	// stellar on all databases and queries writing results to file
    SEQAN_PROTIMESTART(timeStellar);
	if (!_stellarOnAll(databases, databaseIDs, queries, queryIDs, options)) return 1;

    if (options.verbose > 0 && options.noRT == false) 
		std::cout << "Running time: " << SEQAN_PROTIMEDIFF(timeStellar) << "s" << std::endl;

	return 0;
}

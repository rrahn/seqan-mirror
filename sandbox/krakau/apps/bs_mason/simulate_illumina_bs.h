// ==========================================================================
//                           simulate_illumina_bs.h
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
//
//
// SIMULATE METHYLATION PATTERNS DEPENDING ON CONTEXT 
// OR
// BASED ON GIVEN METHYLATION RATES
//
// TODO
// Lister protocol
// 

#ifndef SANDBOX_KRAKAU_APPS_BS_MASON_SIMULATE_ILLUMINA_BS_H_
#define SANDBOX_KRAKAU_APPS_BS_MASON_SIMULATE_ILLUMINA_BS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/misc/misc_cmdparser.h>
#include <seqan/store.h>

#include "bs_mason.h"

using namespace seqan;

// ============================================================================
// Enums, Tags, Classes.
// ============================================================================

struct IlluminaReadsBS_;
typedef Tag<IlluminaReadsBS_> IlluminaReadsBS;

template<>
struct Options<IlluminaReadsBS> : public Options<IlluminaReads>
{
    double probabilityMethylationCG;
    double probabilityMethylationCHG;
    double probabilityMethylationCHH;
    // Protocol: Lister (1) 
    unsigned sequencingProtocol;
    // Rate of unmethylated Cs converted to Ts
    double conversionRate;
    // For haplotype creation: 
    // rate of unmethylated Cs becoming methylated
    double haplotypeMethylInsRate;
    // rate of methylated Cs becoming unmethylated
    double haplotypeMethylDelRate;

    // For import of given methylation rates
    bool useMethylRates;
    CharString methylRatesFile;    

    Options():
                // Context dependent methylation probabilities:
                probabilityMethylationCG(0.25),
                probabilityMethylationCHG(0.06),
                probabilityMethylationCHH(0.01),
                // 
                sequencingProtocol(1),
                //
                conversionRate(0.90),
                //
                haplotypeMethylInsRate(0.001),
                haplotypeMethylDelRate(0.001),
                //
                useMethylRates(false)
    {}                
};

template<>
struct ModelParameters<IlluminaReadsBS> : public ModelParameters<IlluminaReads>
{
    // For import of given methylation rates:
    StringSet<String<double> > topMethylRates; 
    StringSet<String<double> > bottomMethylRates;
    
    // Indicates methylation positions on contig 
    // StringSet: corresponding to contigs in fragmentStore
    // True if methylation at corresponding position
    StringSet<String<bool> > topMethylPositions;
    StringSet<String<bool, Journaled<Alloc<> > > > topMethylPositionsHaplotype;  
    StringSet<String<bool> > bottomMethylPositions;            
    StringSet<String<bool, Journaled<Alloc<> > > > bottomMethylPositionsHaplotype;

};

template <>
struct ReadSimulationInstruction<IlluminaReadsBS> : public ReadSimulationInstruction<IlluminaReads> 
{
    // Bs conversion String: regarding original read before base call errors (editString) are simulated
    // At the beginning: contains FALSE where methylations in reference prevent bs conversions
    // After applySimulationInstructions: contains TRUE where bs conversions where applied (only C or G positions)
    String<bool> bsConversionString;
    // Wheather or not the read is coming from the original top strand 
    bool isFromTopStrand;

};

// ============================================================================
// Metafunctions.
// ============================================================================

// ============================================================================
// Functions.
// ============================================================================

// Extract Options<IlluminaReads> out of Options<IlluminaReadsBS>
void extractOptionsBase(Options<IlluminaReads> & options, Options<IlluminaReadsBS> const & optionsBS)
{
    assign(options.readLength, optionsBS.readLength);
    assign(options.probabilityInsert, optionsBS.probabilityInsert);
    assign(options.probabilityDelete, optionsBS.probabilityDelete);
    assign(options.probabilityMismatchFromFile, optionsBS.probabilityMismatchFromFile);
    assign(options.probabilityMismatchScale, optionsBS.probabilityMismatchScale);
    assign(options.probabilityMismatch, optionsBS.probabilityMismatch);
    assign(options.probabilityMismatchBegin, optionsBS.probabilityMismatchBegin);
    assign(options.probabilityMismatchEnd, optionsBS.probabilityMismatchEnd);
    assign(options.positionRaise, optionsBS.positionRaise);
    assign(options.illuminaNoN, optionsBS.illuminaNoN);
    assign(options.meanQualityBegin, optionsBS.meanQualityBegin);
    assign(options.meanQualityEnd, optionsBS.meanQualityEnd);
    assign(options.stdDevQualityBegin, optionsBS.stdDevQualityBegin);
    assign(options.stdDevQualityEnd, optionsBS.stdDevQualityEnd);
    assign(options.meanMismatchQualityBegin, optionsBS.meanMismatchQualityBegin);
    assign(options.meanMismatchQualityEnd, optionsBS.meanMismatchQualityEnd);
    assign(options.stdDevMismatchQualityBegin, optionsBS.stdDevMismatchQualityBegin);
    assign(options.stdDevMismatchQualityEnd, optionsBS.stdDevMismatchQualityEnd);
}
// Extract ModelParameters <IlluminaReads> out of ModelParameters <IlluminaReadsBS>
void extractModelParametersBase(ModelParameters<IlluminaReads> & parameters, ModelParameters<IlluminaReadsBS> const & parametersBS)
{
    assign(parameters.mismatchProbabilities, parametersBS.mismatchProbabilities);
    assign(parameters.mismatchQualityMeans, parametersBS.mismatchQualityMeans);
    assign(parameters.mismatchQualityStdDevs, parametersBS.mismatchQualityStdDevs);
    assign(parameters.qualityMeans, parametersBS.qualityMeans);
    assign(parameters.qualityStdDevs, parametersBS.qualityStdDevs);
}
// Merge Options<IlluminaReads> into Options<IlluminaReadsBS>
void mergeBaseIntoBSOptions(Options<IlluminaReadsBS> & optionsBS, Options<IlluminaReads> const & options)
{
    assign(optionsBS.readLength, options.readLength);
    assign(optionsBS.probabilityInsert, options.probabilityInsert);
    assign(optionsBS.probabilityDelete, options.probabilityDelete);
    assign(optionsBS.probabilityMismatchFromFile, options.probabilityMismatchFromFile);
    assign(optionsBS.probabilityMismatchScale, options.probabilityMismatchScale);
    assign(optionsBS.probabilityMismatch, options.probabilityMismatch);
    assign(optionsBS.probabilityMismatchBegin, options.probabilityMismatchBegin);
    assign(optionsBS.probabilityMismatchEnd, options.probabilityMismatchEnd);
    assign(optionsBS.positionRaise, options.positionRaise);
    assign(optionsBS.illuminaNoN, options.illuminaNoN);
    assign(optionsBS.meanQualityBegin, options.meanQualityBegin);
    assign(optionsBS.meanQualityEnd, options.meanQualityEnd);
    assign(optionsBS.stdDevQualityBegin, options.stdDevQualityBegin);
    assign(optionsBS.stdDevQualityEnd, options.stdDevQualityEnd);
    assign(optionsBS.meanMismatchQualityBegin, options.meanMismatchQualityBegin);
    assign(optionsBS.meanMismatchQualityEnd, options.meanMismatchQualityEnd);
    assign(optionsBS.stdDevMismatchQualityBegin, options.stdDevMismatchQualityBegin);
    assign(optionsBS.stdDevMismatchQualityEnd, options.stdDevMismatchQualityEnd);
}
// Merge ModelParameters<IlluminaReads> into ModelParameters<IlluminaReadsBS>
void mergeBaseIntoBSModelParameters(ModelParameters<IlluminaReadsBS> & parametersBS, ModelParameters<IlluminaReads> & parameters)
{
    assign(parametersBS.mismatchProbabilities, parameters.mismatchProbabilities);
    assign(parametersBS.mismatchQualityMeans, parameters.mismatchQualityMeans);
    assign(parametersBS.mismatchQualityStdDevs, parameters.mismatchQualityStdDevs);
    assign(parametersBS.qualityMeans, parameters.qualityMeans);
    assign(parametersBS.qualityStdDevs, parameters.qualityStdDevs);
}

template <typename TStream>
TStream & operator<<(TStream & stream, Options<IlluminaReadsBS> const & options) {
    stream << static_cast<Options<IlluminaReads> >(options);
    stream << "bisulfite-options {" << std::endl
           << "  probabilityMethylationCG;               " << options.probabilityMethylationCG << std::endl
           << "  probabilityMethylationCHG;              " << options.probabilityMethylationCHG << std::endl           
           << "  probabilityMethylationCHH;              " << options.probabilityMethylationCHH<< std::endl
           << "  sequencingProtocol;                     " << options.sequencingProtocol << std::endl
           << "  conversionRate                          " << options.conversionRate << std::endl
           << "  useMethylRates:                    " << (options.useMethylRates ? "true" : "false") << std::endl
           << "  methylRatesFile:                 \"" << options.methylRatesFile << "\"" << std::endl
           << "}" << std::endl;
    return stream;
}

void setUpCommandLineParser(CommandLineParser & parser,
                            IlluminaReadsBS const &)
{
    setUpCommandLineParser(parser, IlluminaReads());

    addSection(parser, "Illumina bisulfite parameters");

    addOption(parser, CommandLineOption("pcg",  "prob-cg", "The probability of an 'C' beeing methylated in CG-context.  Default: ?.", OptionType::Integer | OptionType::Label));
    addOption(parser, CommandLineOption("pchg",  "prob-chg", "The probability of an 'C' beeing methylated in CHG-context.  Default: ?.", OptionType::Integer | OptionType::Label));
    addOption(parser, CommandLineOption("pchh",  "prob-chh", "The probability of an 'C' beeing methylation in CHH-context.  Default: ?.", OptionType::Integer | OptionType::Label));

    addOption(parser, CommandLineOption("p",  "sequence-protocol", "The sequence protocol used (Lister or Cokus).  Default: Cokus.", OptionType::Integer | OptionType::Label));
    addOption(parser, CommandLineOption("c", "conversion-rate", "The C-T conversion rate.  Default: ?.", OptionType::Double));

    addOption(parser, CommandLineOption("umr", "use-methyl-rates", "Use given methylation rates.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("mr",  "methyl-rates-file", "Input file with given methylation rates.  Default: \"\".", OptionType::String));
}

int parseCommandLineAndCheckModelSpecific(Options<IlluminaReadsBS> & options,
                                          CommandLineParser & parser)
{
    Options<IlluminaReads> optionsBase;
    extractOptionsBase(optionsBase, options);  
    int ret = parseCommandLineAndCheckModelSpecific(optionsBase, parser);
    mergeBaseIntoBSOptions(options, optionsBase);

    if (ret != 0)
        return ret;

    if (isSetLong(parser, "prob-cg"))
        getOptionValueLong(parser, "prob-cg", options.probabilityMethylationCG);
    if (isSetLong(parser, "prob-chg"))
        getOptionValueLong(parser, "prob-chg", options.probabilityMethylationCG);
    if (isSetLong(parser, "prob-chh"))
        getOptionValueLong(parser, "prob-chh", options.probabilityMethylationCG);

    
    if (isSetLong(parser, "sequence-protocol"))
        getOptionValueLong(parser, "sequence-protocol", options.sequencingProtocol);
    if (isSetLong(parser, "conversion-rate"))
        getOptionValueLong(parser, "conversion-rate", options.conversionRate);

    if (isSetLong(parser, "use-methyl-rates"))
        options.useMethylRates = true;
    if (isSetLong(parser, "methyl-rates-file"))
        getOptionValueLong(parser, "methyl-rates-file", options.methylRatesFile);

    return 0;
}

// Load given methylation rates from a file.
int loadMethylRates(ModelParameters<IlluminaReadsBS> & parameters, 
                         Options<IlluminaReadsBS> const & options, 
                         FragmentStore<MyFragmentStoreConfig> & fragmentStore)
{
    typedef Stream<std::fstream> TStream;
    typedef RecordReader<std::fstream, SinglePass<> > TRecordReader;

    std::fstream inputFile(toCString(options.methylRatesFile), std::ios::binary | std::ios::in);
    TRecordReader reader(inputFile);

    // Resize methylation rate strings corresponding contig lengths   
    clear(parameters.topMethylRates);
    clear(parameters.bottomMethylRates);
    resize(parameters.topMethylRates, length(fragmentStore.contigStore), Exact());
    resize(parameters.bottomMethylRates, length(fragmentStore.contigStore), Exact());
    for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i){
        resize(value(parameters.topMethylRates, i), length(fragmentStore.contigStore[i].seq), 0.0, Exact());      
        resize(value(parameters.bottomMethylRates, i), length(fragmentStore.contigStore[i].seq), 0.0, Exact());     
    }

    CharString contigName;
    CharString contigName_Prev = "";
    CharString position;
    CharString strand;
    CharString methylationRate;

    size_t contigId = 0;
    refresh(fragmentStore.contigNameStoreCache);
    while(!atEnd(reader))
    {
        clear(contigName);
        clear(position);
        clear(strand);
        clear(methylationRate);
        // Read contigName
        readUntilWhitespace(contigName, reader);
        skipWhitespaces(reader);
        // Read position
        readUntilWhitespace(position, reader);
        skipWhitespaces(reader);
        // Read strand
        readUntilWhitespace(strand, reader);
        skipWhitespaces(reader);
        // Read methylationRate
        readUntilWhitespace(methylationRate, reader);
        skipLine(reader);

        // If 1st occurence of contigName:
        if (contigName != contigName_Prev){
            assign(contigName_Prev, contigName);
            contigId = 0;
            // Get contigId for current contigName out of contigNameStore 
            if (!getIdByName(fragmentStore.contigNameStore, contigName, contigId, fragmentStore.contigNameStoreCache)) {
                        std::cerr << "ERROR: Could not find contig with name \"" << contigName << "\" (from methylation rates file) in contigs." << std::endl;
                        return 1;
            }         
        }
        // Assign methylation rate to top/bottomMethylRate strings
        if (strand == '+'){
            SEQAN_ASSERT_EQ(fragmentStore.contigStore[contigId].seq[lexicalCast<size_t>(position)], 'C'); 
            assignValue(value(parameters.topMethylRates, contigId), lexicalCast<size_t>(position), lexicalCast<double>(methylationRate));   
        } else {
            SEQAN_ASSERT_EQ(fragmentStore.contigStore[contigId].seq[lexicalCast<size_t>(position)], 'G'); 
            assignValue(value(parameters.bottomMethylRates, contigId), lexicalCast<size_t>(position), lexicalCast<double>(methylationRate)); 
        }
    }
    return 0;
}

// Calls original functions for non-BS tags
template<typename TModelParameters, typename TOptions>
int simulateReadsSetupModelSpecificData(TModelParameters & parameters,
                                        TOptions const & options,
                                        FragmentStore<MyFragmentStoreConfig> & /*NOP*/)
{
    int ret = simulateReadsSetupModelSpecificData(parameters, options);
    if (ret != 0)
        return ret;

    return 0;
}

// Bs_change: need fragmentStore for contig size -> size of methylationRates strings
int simulateReadsSetupModelSpecificData(ModelParameters<IlluminaReadsBS> & parameters,
                                        Options<IlluminaReadsBS> const & options,
                                        FragmentStore<MyFragmentStoreConfig> & fragmentStore)
{
    // Build standard illumina parameters:
    ModelParameters<IlluminaReads> parametersBase;
    // Extract IlluminaReads parameters out of IlluminaReadsBS parameters 
    extractModelParametersBase(parametersBase, parameters); 
    int ret = simulateReadsSetupModelSpecificData(parametersBase, static_cast<Options<IlluminaReads> >(options) );
    // Merge IlluminaReads parameters again into IlluminaReadsBS parameters
    mergeBaseIntoBSModelParameters(parameters, parametersBase);
    if (ret != 0)
        return ret;
   
    // Read methylation rates if given
    if (options.useMethylRates){
        std::cerr << "Loading methylation rates from \"" << options.methylRatesFile << "\"" << std::endl;
        ret = loadMethylRates(parameters, options, fragmentStore);
        if (ret != 0)
            return ret;
    }

    return 0;
}


// Does actualy nothing for non-BS read tags
template<typename TReadsTag, typename TRNG>
int simulateMethylPositions(ModelParameters<TReadsTag> & /*NOP*/, TRNG & /*NOP*/, Options<TReadsTag> const & /*NOP*/, FragmentStore<MyFragmentStoreConfig> & /*NOP*/)
{
    return 0;
}

// Simulates methylation state of each C in reference dependent on context and its given methylation probability
template<typename TRNG>
int simulateMethylPositions(ModelParameters<IlluminaReadsBS> & parameters, TRNG & rng, Options<IlluminaReadsBS> const & options, FragmentStore<MyFragmentStoreConfig> & fragmentStore)
{
    // Only if no given methylation rates
    // For given methylation rates, methylation position will be simulated later for each Haplotype corresponding methylation rates 
    if (!options.useMethylRates){
        clear(parameters.topMethylPositions);
        clear(parameters.bottomMethylPositions);
        resize(parameters.topMethylPositions, length(fragmentStore.contigStore), Exact());
        resize(parameters.bottomMethylPositions, length(fragmentStore.contigStore), Exact());

        // For Cs in topStrand:
        for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i) {
            String<Dna5> const & contig = fragmentStore.contigStore[i].seq;
            clear(value(parameters.topMethylPositions, i));
            resize(value(parameters.topMethylPositions, i), length(contig), false, Exact());
            for (size_t j = 0; j < length(contig)-2; ++j){  // not to the end because of context
                // CG context:
                if (contig[j] == 'C' && contig[j+1] == 'G') {
                    double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
                    if (x < options.probabilityMethylationCG){
                        assignValue(value(parameters.topMethylPositions, i), j, true);           
                    }
                }
                // CHG context:
                else if (contig[j] == 'C' && contig[j+2] == 'G') {
                    double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
                    if (x < options.probabilityMethylationCHG){
                        assignValue(value(parameters.topMethylPositions, i), j, true);               
                    }
                }
                // CHH context:
                else if (contig[j] == 'C') {
                    double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
                    if (x < options.probabilityMethylationCHH){
                        assignValue(value(parameters.topMethylPositions, i), j, true);
                    }
                }
            }
        }
        // For Cs in bottomStrand:
        for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i) {
            String<Dna5> const & contig = fragmentStore.contigStore[i].seq;
            Dna5StringReverseComplement revCompl(contig);
            clear(value(parameters.bottomMethylPositions, i));
            resize(value(parameters.bottomMethylPositions, i), length(revCompl), false, Exact());
            for (size_t j = 0; j < length(revCompl)-2; ++j){  // Not to the end because of context length
                // CG context:
                if (revCompl[j] == 'C' && revCompl[j+1] == 'G') {
                    double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
                    if (x < options.probabilityMethylationCG){
                        assignValue(value(parameters.bottomMethylPositions, i), j, true);
                    }
                }
                // CHG context:
                else if (revCompl[j] == 'C' && revCompl[j+2] == 'G') {
                    double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
                    if (x < options.probabilityMethylationCHG){
                        assignValue(value(parameters.bottomMethylPositions, i), j, true);
                    }
                }
                // CHH context:
                else if (revCompl[j] == 'C') {
                    double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
                    if (x < options.probabilityMethylationCHH){
                        assignValue(value(parameters.bottomMethylPositions, i), j, true);
                    }
                }
            }
            // Reverse: methylations regarding to positions on top forward strand  
            reverse(value(parameters.bottomMethylPositions, i));
        }
    }
    return 0;
}


template <typename TRNG>
unsigned pickReadLength(TRNG const &, Options<IlluminaReadsBS> const & options)
{
    return options.readLength;
}


// Build a haplotype, based on the contigs from the given fragment store.
// And adjust methylation patterns ( forwardMethylationPositions, reverse ...)
template <typename TRNG>
void buildHaplotype(StringSet<String<Dna5, Journaled<Alloc<> > > > & haplotype,
                    ModelParameters<IlluminaReadsBS> & parameters,
                    FragmentStore<MyFragmentStoreConfig> & fragmentStore,
                    TRNG & rng,
                    Options<IlluminaReadsBS> const & options) {
    resize(haplotype, length(fragmentStore.contigStore), Exact());
    String<Dna5> buffer;
    reserve(buffer, options.haplotypeIndelRangeMax);

    // Bs_change:
    // Without given methylation rates: use already simulated methyl positions and adjust them
    // With given methylation rates: methylation position strings are empty, resize will happen later
    /*
    assign(parameters.topMethylPositionsHaplotype, parameters.topMethylPositions);
    assign(parameters.bottomMethylPositionsHaplotype, parameters.bottomMethylPositions);
    */
    // Journaled strings: store only differences in haplotypes
    clear(parameters.topMethylPositionsHaplotype);
    clear(parameters.bottomMethylPositionsHaplotype);
    resize(parameters.topMethylPositionsHaplotype, length(parameters.topMethylPositions));
    resize(parameters.bottomMethylPositionsHaplotype, length(parameters.bottomMethylPositions));
    for (unsigned i = 0; i < length(parameters.topMethylPositions); ++i){
        setHost(parameters.topMethylPositionsHaplotype[i], parameters.topMethylPositions[i]);
        setHost(parameters.bottomMethylPositionsHaplotype[i], parameters.bottomMethylPositions[i]);
    }

    String<bool> bs_buffer;

    // If simulated methylation positions: 
    // apply haplotype methylation changes 
    if (!options.useMethylRates)
        buildMethylHaplotypeUseMethylPositions(parameters, fragmentStore, rng, options);
    // If given methylation rates: 
    // simulate haplotype methylation positions corresponding given methylation rates
    else
        buildMethylHaplotypeUseMethylRates(parameters, rng);
    
    // Apply general haplotype changes and adjust methylation positions if necessary
    for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i) {
        std::cout << "    contig # " << i+1 << "/" << length(fragmentStore.contigStore) << std::endl;
        clear(haplotype[i]);
        setHost(haplotype[i], fragmentStore.contigStore[i].seq);
        String<Dna5> const & contig = fragmentStore.contigStore[i].seq;
        String<Dna5, Journaled<Alloc<> > > & haplotypeContig = haplotype[i];

        // Only generate Ns in the haplotype if allowed.
        int maxOrdValue = options.haplotypeNoN ? 3 : 4;

        // j is position in original sequence, k is position in haplotype
        for (size_t j = 0, k = 0; j < length(contig);) {
            double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
            if (x < options.haplotypeSnpRate) {
                // SNP
                Dna5 c = Dna5(pickRandomNumber(rng, Pdf<Uniform<int> >(0, maxOrdValue - 1)));
                if (c == contig[j])
                    c = Dna5(ordValue(c) + 1);
                if (options.haplotypeNoN)
                    SEQAN_ASSERT(c != Dna5('N'));
                assignValue(haplotypeContig, k, c);
                // Bs_change: adjust methylation state if C is changed to different base
                if (contig[j] == 'C'){
                    assignValue(value(parameters.topMethylPositionsHaplotype, i), k, false);
                } else if (contig[j] == 'G'){
                    assignValue(value(parameters.bottomMethylPositionsHaplotype, i), k, false);
                }
                // TODO if convertion into 'C': maybe methylate it, depending on context and probabilities ?
                j += 1;
                k += 1;
            } else if (x < options.haplotypeSnpRate + options.haplotypeIndelRate) {
                // Indel of random length.
                unsigned rangeLen = options.haplotypeIndelRangeMax - options.haplotypeIndelRangeMin;
                unsigned indelLen = options.haplotypeIndelRangeMin + static_cast<unsigned>(pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1)) * rangeLen);
                if (pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1)) < 0.5) {
                    // Insertion.
                    clear(buffer);
                    for (unsigned ii = 0; ii < indelLen; ++ii)
                        appendValue(buffer, Dna5(pickRandomNumber(rng, Pdf<Uniform<int> >(0, maxOrdValue))));
                    insert(haplotypeContig, k, buffer);
                    // Bs_change: adjust length of methylation string
                    clear(bs_buffer);
                    resize(bs_buffer, indelLen, false, Exact());
                    insert(value(parameters.topMethylPositionsHaplotype, i), k, bs_buffer);
                    insert(value(parameters.bottomMethylPositionsHaplotype, i), k, bs_buffer);
                    // TODO if insertion of 'C': maybe methylate it, depending on context and probabilities ?
                    k += indelLen;
                } else {
                    // Deletion.
                    indelLen = _min(indelLen, length(haplotypeContig) - k);
                    erase(haplotypeContig, k, k + indelLen);
                    // Bs_change: adjust length of methylation string
                    erase(value(parameters.topMethylPositionsHaplotype, i), k, k + indelLen);
                    erase(value(parameters.bottomMethylPositionsHaplotype, i), k, k + indelLen); 
                    j += indelLen;
                }
            } else {
                // Match.
                j += 1;
                k += 1;
            }
        }
        SEQAN_ASSERT_EQ(length(haplotypeContig), length(parameters.topMethylPositionsHaplotype[i]));
    }
}

// Build haplotype methylation changes
// Based on simulated methylation positions
// No methylation rates
template <typename TRNG>
void buildMethylHaplotypeUseMethylPositions(ModelParameters<IlluminaReadsBS> & parameters,
                    FragmentStore<MyFragmentStoreConfig> & fragmentStore,
                    TRNG & rng,
                    Options<IlluminaReadsBS> const & options) {
    // Apply methylation changes for haplotype
    // just corresponding to given probability for haplotype methylation changes, independent of context etc.
    SEQAN_ASSERT_EQ(length(fragmentStore.contigStore), length(parameters.topMethylPositionsHaplotype));
    for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i) {
        String<Dna5> const & contig = fragmentStore.contigStore[i].seq;
        SEQAN_ASSERT_EQ(length(contig), length(parameters.topMethylPositionsHaplotype[i]));
        for (unsigned j = 0; j < length(parameters.topMethylPositionsHaplotype[i]); ++j){
            // topStrand
            double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
            if (x < options.haplotypeMethylInsRate && contig[j] == 'C'){
                // Methylation insertion
                if (getValue(parameters.topMethylPositionsHaplotype[i],j) == false)
                    assignValue(value(parameters.topMethylPositionsHaplotype, i), j, true);
            } else if (x < options.haplotypeMethylInsRate + options.haplotypeMethylDelRate){
                // Methylation deletion
                if (getValue(parameters.topMethylPositionsHaplotype[i], j) == true)
                    assignValue(value(parameters.topMethylPositionsHaplotype, i), j, false);
            } 
            // bottomStrand 
            x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
            if (x < options.haplotypeMethylInsRate && contig[j] == 'G'){
                // Methylation insertion
                if (getValue(parameters.bottomMethylPositionsHaplotype[i], j) == false)
                    assignValue(value(parameters.bottomMethylPositionsHaplotype, i), j, true);
            } else if (x < options.haplotypeMethylInsRate + options.haplotypeMethylDelRate){
                // Methylation deletion
                if (getValue(parameters.bottomMethylPositionsHaplotype[i], j) == true)
                    assignValue(value(parameters.bottomMethylPositionsHaplotype, i), j, false);
            }  
        }
    }
}

// Build methylation haplotype
// Based on given methylation rates
template <typename TRNG>
void buildMethylHaplotypeUseMethylRates(ModelParameters<IlluminaReadsBS> & parameters,
                    TRNG & rng) {
    // Simulate haplotype methyl positions based on given methylation rate
    // Methylation positions strings are not yet initialised (function simulateMethylPositions not called for given methylation rates)
    // Use of journaled string: 
    // At this point parameters.topMethylPositions is still is not yet filled
    // We store the methyl positions in parameters.topMethylPositions instead of in parameters.topMethylPositionsHaplotype
    // -> indirect stored in parameters.topMethylPositionsHaplotype
    // (no differences, no extra space)
    clear(parameters.topMethylPositions);
    clear(parameters.bottomMethylPositions);
    resize(parameters.topMethylPositions, length(parameters.topMethylRates), Exact());
    resize(parameters.bottomMethylPositions, length(parameters.bottomMethylRates), Exact());
    resize(parameters.topMethylPositionsHaplotype, length(parameters.topMethylRates), Exact());
    resize(parameters.bottomMethylPositionsHaplotype, length(parameters.bottomMethylRates), Exact());

    for (unsigned i = 0; i < length(parameters.topMethylPositions); ++i) {
        clear(value(parameters.topMethylPositions, i));
        clear(value(parameters.bottomMethylPositions, i));
        resize(value(parameters.topMethylPositions, i), length(parameters.topMethylRates[i]), false, Exact());
        resize(value(parameters.bottomMethylPositions, i), length(parameters.bottomMethylRates[i]), false, Exact());
        for (unsigned j = 0; j < length(parameters.topMethylPositions[i]); ++j){
            // topStrand
            double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
            if (x < getValue(parameters.topMethylRates[i], j)){
                // Methylation insertion
                assignValue(value(parameters.topMethylPositions, i), j, true);
            } 
            x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
            if (x < getValue(parameters.bottomMethylRates[i], j)){
                // Methylation insertion
                assignValue(value(parameters.bottomMethylPositions, i), j, true);
            } 
        }
        // Transfer content to haplotype contig
        clear(parameters.topMethylPositionsHaplotype[i]);
        clear(parameters.bottomMethylPositionsHaplotype[i]);
        setHost(parameters.topMethylPositionsHaplotype[i], parameters.topMethylPositions[i]);
        setHost(parameters.bottomMethylPositionsHaplotype[i], parameters.bottomMethylPositions[i]);
    }
}

template<typename TRNG, typename TReadSimulationInstruction>
void pickOriginalStrand(TRNG & /*NOP*/, TReadSimulationInstruction & /*NOP*/){}

template<typename TRNG>
void pickOriginalStrand(TRNG & rng, ReadSimulationInstruction<IlluminaReadsBS> & inst){
    inst.isFromTopStrand = pickRandomNumber(rng, Pdf<Uniform<int> >(0, 1));
}

template <typename TRNG, typename TContig>
void buildSimulationInstructions(ReadSimulationInstruction<IlluminaReadsBS> & inst, TRNG & rng, unsigned readLength, TContig const & contig, ModelParameters<IlluminaReadsBS> const & parameters, Options<IlluminaReadsBS> const & options) {
    clear(inst.editString);
    reserve(inst.editString, static_cast<size_t>(1.2 * readLength), Generous());
    inst.delCount = 0;
    inst.insCount = 0;

    //
    // Build Edit String.    //
    for (unsigned i = 0; i < readLength; /*NOP*/) {
        double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
        double pMismatch = parameters.mismatchProbabilities[i];
        double pInsert   = options.probabilityInsert;
        double pDelete   = options.probabilityDelete;
        double pMatch    = 1.0 - pMismatch - pInsert - pDelete;
        if (x < pMatch) {
            // match
            i += 1;
            appendValue(inst.editString, ERROR_TYPE_MATCH);
        } else if (x < pMatch + pMismatch) {
            // mismatch
            i += 1;
            appendValue(inst.editString, ERROR_TYPE_MISMATCH);
        } else if (x < pMatch + pMismatch + pInsert) {
            // insert
            if (length(inst.editString) > 0 && back(inst.editString == ERROR_TYPE_DELETE)) {
                inst.delCount -= 1;
                eraseBack(inst.editString);
            } else {
                i += 1;
                inst.insCount += 1;
                appendValue(inst.editString, ERROR_TYPE_INSERT);
            }
        } else {
            // Decrement string size, do not add a delete if string is
            // too short, possibly remove insert from edit string.
            if (length(inst.editString) > 0) {
                if (back(inst.editString == ERROR_TYPE_INSERT)) {
                    i -= 1;
                    inst.insCount -= 1;
                    eraseBack(inst.editString);
                } else {
                    inst.delCount += 1;
                    appendValue(inst.editString, ERROR_TYPE_DELETE);
                }
            }
        }
    }
    SEQAN_ASSERT_EQ(readLength, length(inst.editString) - inst.delCount);

    //
    // Adjust Positions.
    //

    // If the number of deletions does not equal the number of inserts
    // then we have to adjust the read positions.
    if (inst.delCount != inst.insCount) {
        int delta = static_cast<int>(inst.delCount) - static_cast<int>(inst.insCount);
        inst.endPos += delta;
        if (inst.endPos > length(contig)) {
            delta = inst.endPos - length(contig);
            inst.endPos -= delta;
            inst.beginPos -= delta;
        }
        SEQAN_ASSERT_EQ(inst.endPos - inst.beginPos + inst.insCount - inst.delCount,
                        readLength);
    }

    //
    // Build bsConversionPositions String (regarding adjusted positions in forward haplotype contig sequence)
    // Apart of begin and end positions, this is independently of the edit string, 
    // because later this bs conversions will be applied before the edit string 
    //
    clear(inst.bsConversionString);
    resize(inst.bsConversionString, inst.endPos - inst.beginPos, false, Exact());
    if (inst.isFromTopStrand) {
        for (unsigned i = 0; i < (inst.endPos - inst.beginPos); ++i) {
            if (getValue(parameters.topMethylPositionsHaplotype[inst.contigId], inst.beginPos + i) == false) {  
                // Only methylated positions are protected from conversion
                // remember later to check if position is C or G
                double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
                if (x < options.conversionRate) {
                    // Bs conversion
                    assignValue(inst.bsConversionString, i, true);
                }
            } else
                SEQAN_ASSERT_EQ(inst.bsConversionString[i], false);                
        }
    } else {
        for (unsigned i = 0; i < (inst.endPos - inst.beginPos); ++i) {
            if (getValue(parameters.bottomMethylPositionsHaplotype[inst.contigId], inst.beginPos + i) == false) {
                double x = pickRandomNumber(rng, Pdf<Uniform<double> >(0, 1));
                if (x < options.conversionRate) {
                    // Bs conversion
                    assignValue(inst.bsConversionString, i, true);
                }
            } else
                SEQAN_ASSERT_EQ(inst.bsConversionString[i], false);
        }
    }
    SEQAN_ASSERT_EQ((inst.endPos - inst.beginPos), length(inst.bsConversionString));

    //
    // Quality Simulation.
    //
    SEQAN_ASSERT_GT(length(inst.editString), 0u);
    if (options.simulateQualities) {
        clear(inst.qualities);
        resize(inst.qualities, length(inst.editString), 0, Exact());

        for (unsigned i = 0, j = 0; i < length(inst.editString); i++) {
            SEQAN_ASSERT_LEQ(j, inst.endPos - inst.beginPos + inst.delCount);
            if (inst.editString[i] == ERROR_TYPE_MISMATCH) {
                // std::cout << "i == " << i << ", j == " << j << ", parameters.mismatchQualityMeans[j] == " << parameters.mismatchQualityMeans[j] << ", parameters.mismatchQualityStdDevs[j] == " << parameters.mismatchQualityStdDevs[j] << std::endl;
                Pdf<Normal> pdf(parameters.mismatchQualityMeans[j], parameters.mismatchQualityStdDevs[j]);
                inst.qualities[i] = static_cast<int>(pickRandomNumber(rng, pdf));
            } else {
                Pdf<Normal> pdf(parameters.qualityMeans[j], parameters.qualityStdDevs[j]);
                inst.qualities[i] = static_cast<int>(pickRandomNumber(rng, pdf));
            }

            if (inst.editString[i] == ERROR_TYPE_MISMATCH || inst.editString[i] == ERROR_TYPE_MATCH)
                j += 1;
        }
    }

}


template <typename TRNG, typename TString>
void applySimulationInstructions(TString & read, TRNG & rng, ReadSimulationInstruction<IlluminaReadsBS> & inst, Options<IlluminaReadsBS> const & options)
{
    typedef typename Value<TString>::Type TAlphabet;

    if (options.simulateQualities)
        SEQAN_ASSERT_EQ(length(inst.qualities), length(inst.editString));
    
    //
    // apply bs conversions
    //    

    // Top strand, forward strand sequence of read is edited
    // For bottom strand and/or reverse strand: modifying will happen later
    if (inst.isFromTopStrand){
        for (unsigned i = 0; i < length(inst.bsConversionString); ++i){
            if (inst.bsConversionString[i] == true && read[i] == 'C'){
                assignValue(read, i, 'T');
            } else {
                // Modify bsConversionString on the fly to a indicator string, if bs conversion was applied or not
                // TRUE only if current base is C and this C is converted
                assignValue(inst.bsConversionString, i, false);
            }
        }
    }else{
        for (unsigned i = 0; i < length(inst.bsConversionString); ++i){
            if (inst.bsConversionString[i] == true && read[i] == 'G'){
                assignValue(read, i, 'A');
            } else {
                assignValue(inst.bsConversionString, i, false);
            }

        }
    }
    
    //
    // apply errors
    //
    TString tmp;
    reserve(tmp, length(read) + inst.insCount - inst.delCount);
    unsigned j = 0;
    for (unsigned i = 0; i < length(inst.editString); ++i) {
        SEQAN_ASSERT_LEQ(j, i);

        TAlphabet c;
        //int x, xold;
        switch (inst.editString[i]) {
            case ERROR_TYPE_MATCH:
                SEQAN_ASSERT_LT_MSG(j, length(read), "i = %u", i);
                appendValue(tmp, read[j]);
                if (options.simulateQualities)
                    assignQualityValue(back(tmp), inst.qualities[i]);
                // std::cout << i << " " << getQualityValue(back(tmp)) << " " << inst.qualities[i] << " " << convert<char>(back(tmp)) << " match" << std::endl;
                //std::cout << back(tmp) << " " << read[j] << " " << inst.qualities[i] << std::endl;
                j += 1;
                break;
            case ERROR_TYPE_MISMATCH:
                if (options.illuminaNoN) {
                    c = TAlphabet(pickRandomNumber(rng, Pdf<Uniform<int> >(0, ValueSize<TAlphabet>::VALUE - 3)));  // -3, N not allowed
                } else {
                    c = TAlphabet(pickRandomNumber(rng, Pdf<Uniform<int> >(0, ValueSize<TAlphabet>::VALUE - 2)));  // -2, N allowed
                }
                //xold = ordValue(c);
                SEQAN_ASSERT_LT_MSG(j, length(read), "i = %u", i);
                if (ordValue(c) >= ordValue(read[j]))
                    c = TAlphabet(ordValue(c) + 1);
                if (options.illuminaNoN)
                    SEQAN_ASSERT(c != TAlphabet('N'));
                //x = ordValue(c);
                appendValue(tmp, c);
                if (options.simulateQualities) {
                    if (options.illuminaNoN)  // Ns can be introduced through quality, too.
                        assignQualityValue(back(tmp), _max(1, inst.qualities[i]));
                    else
                        assignQualityValue(back(tmp), inst.qualities[i]);
                }
                // std::cout << i << " q(q_i)=" << getQualityValue(back(tmp)) << " q(i)=" << inst.qualities[i] << " char=" << convert<char>(back(tmp)) << " c_old=" << xold << " c=" << x << " r_j=" << ordValue(read[j]) << std::endl;
                // std::cout << i << " " << getQualityValue(back(tmp)) << " " << inst.qualities[i] << " " << convert<char>(back(tmp)) << " mismatch" << std::endl;
                //std::cout << "MM " << c << " " << back(tmp) << " " << inst.qualities[i] << std::endl;
                j += 1;
                break;
            case ERROR_TYPE_INSERT:
                if (options.illuminaNoN)
                    appendValue(tmp, TAlphabet(pickRandomNumber(rng, Pdf<Uniform<int> >(0, ValueSize<TAlphabet>::VALUE - 2))));  // -2 == no N
                else
                    appendValue(tmp, TAlphabet(pickRandomNumber(rng, Pdf<Uniform<int> >(0, ValueSize<TAlphabet>::VALUE - 1))));  // -1 == N allowed
                if (options.simulateQualities) {
                    if (options.illuminaNoN)  // Ns can be introduced through quality, too.
                        assignQualityValue(back(tmp), _max(1, inst.qualities[i]));
                    else
                        assignQualityValue(back(tmp), inst.qualities[i]);
                }
                // std::cout << i << " " << getQualityValue(back(tmp)) << " " << inst.qualities[i] << " " << convert<char>(back(tmp)) << " insertion" << std::endl;
                break;
            case ERROR_TYPE_DELETE:
                j += 1;
                break;
            default:
                SEQAN_ASSERT_FAIL("Invalid error type.");
        }
    }
    SEQAN_ASSERT_EQ(j, length(read));
    SEQAN_ASSERT_GEQ(length(tmp), options.readLength);

    //std::cout << "tmp == " << tmp << std::endl;
    resize(tmp, options.readLength, Exact());
    move(read, tmp);
}

// Performs the actual read simulation.
template <typename TRNG, typename TOptions>
int simulateReadsMain(FragmentStore<MyFragmentStoreConfig> & fragmentStore,
                      TRNG & rng,
                      TOptions const & options,
                      ModelParameters<IlluminaReadsBS> & parameters) {        // bs_change: not const, need to adjust a few parameters
    typedef FragmentStore<MyFragmentStoreConfig> TFragmentStore;
    typedef Value<TFragmentStore::TMatePairStore>::Type TMatePairStoreElement;
    typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type	TAlignedElement;
    typedef typename TAlignedElement::TGapAnchors								TReadGapAnchors;

    if (options.verbose)
        std::cerr << "Simulating reads..." << std::endl;

    typedef Position<CharString>::Type TPos;

    // Number of reads comes from command line by default.  If sample
    // counts are given, we compute it from this instead.
    size_t numReads = options.numReads;
    if (length(parameters.sampleCounts) > 0u) {
        numReads = 0;
        for (unsigned i = 0; i < length(parameters.sampleCounts); ++i)
            numReads += parameters.sampleCounts[i];
    }

    // First, we randomly pick the haplotype for each read to be
    // simulated or read it from the sample counts in parameters.
    String<unsigned> haplotypeIds;
    // Pick random haplotype origin.
    reserve(haplotypeIds, numReads);
    for (size_t i = 0; i < numReads; ++i)
        appendValue(haplotypeIds, pickRandomNumber(rng, Pdf<Uniform<unsigned> >(0, options.numHaplotypes - 1)));

    // Maybe pick contig ids to sample from.
    String<unsigned> contigIds;
    if (length(parameters.sampleCounts) > 0) {
        for (unsigned i = 0; i < length(fragmentStore.contigNameStore); ++i) {
            resize(contigIds, length(contigIds) + parameters.sampleCounts[i], i);
            if (options.veryVerbose)
                std::cerr << parameters.sampleCounts[i] << " reads from haplotype " << i << "..." << std::endl;
        }
        shuffle(contigIds, rng);
    }

    // We do not build all haplotypes at once since this could cost a
    // lot of memory.
    //
    // TODO(holtgrew): Would only have to switch pointers to journals which is possible.
    //
    // for each haplotype id
    //   simulate haplotype
    //   for each simulation instruction for this haplotype:
    //     build simulated read
    reserve(fragmentStore.readSeqStore, numReads, Exact());
    reserve(fragmentStore.readNameStore, numReads, Exact());
    CharString readNameBuf;
    char outFileName[151];
    snprintf(outFileName, 150, "%s", toCString(options.readNamePrefix));
    for (unsigned haplotypeId = 0; haplotypeId < options.numHaplotypes; ++haplotypeId) {
        std::cerr << "Simulating for haplotype #" << haplotypeId << "..." << std::endl;
        std::cout << "  Building haplotype..." << std::endl;
        StringSet<String<Dna5, Journaled<Alloc<> > > > haplotypeContigs;
        // bs_change
        // need of parameters, adjust corresponding to haplotype changes
        buildHaplotype(haplotypeContigs, parameters, fragmentStore, rng, options);
        // TODO(holtgrew): Assigning of string set with compatible string should be possible.
        StringSet<String<Dna5> > haplotypeContigsCopy;
        for (unsigned i = 0; i < length(haplotypeContigs); ++i)
            appendValue(haplotypeContigsCopy, haplotypeContigs[i]);

        // Build partial sums over relative contig lengths so we can pick the contigs later on.
        size_t totalLength = 0;
        for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i)
            totalLength += length(fragmentStore.contigStore[i].seq);
        String<double> relativeContigLengths;
        resize(relativeContigLengths, length(fragmentStore.contigStore) + 1, Exact());
        front(relativeContigLengths) = 0.0;
        for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i) {
            double l = static_cast<double>(length(fragmentStore.contigStore[i].seq));
            relativeContigLengths[i + 1] = l / totalLength;
        }
        std::partial_sum(begin(relativeContigLengths), end(relativeContigLengths), begin(relativeContigLengths));
        back(relativeContigLengths) = 1.0;

        // Simulate the reads...
        std::cerr << "  Simulating reads for haplotype #" << haplotypeId << "..." << std::endl;

//         std::cerr << "Journal: " << haplotypeContigs[0]._journalEntries << std::endl;

        for (unsigned j = 0; j < length(haplotypeIds); ++j) {
            if (haplotypeIds[j] != haplotypeId)
                continue;  // Guard against instructions on wrong haplotype.

            // Build simulation instructions.
            String<ReadSimulationInstruction<IlluminaReadsBS> > instructions;
            // TODO(holtgrew): Pick contig id outside of instructions.
            size_t contigId = 0;
            bool fixedContigId = false;
            if (length(parameters.sampleCounts) > 0) {
                fixedContigId = true;
                if (options.generateMatePairs)
                    contigId = contigIds[length(fragmentStore.readSeqStore) / 2];
                else
                    contigId = contigIds[length(fragmentStore.readSeqStore)];
            }
            int res = buildReadSimulationInstruction(instructions, rng, haplotypeId, haplotypeContigsCopy, relativeContigLengths, contigId, fixedContigId, parameters, options);
            if (res != 0)
                return res;

            for (unsigned k = 0; k < length(instructions); ++k) {
                ReadSimulationInstruction<IlluminaReadsBS> & inst = instructions[k];
                // Reverse edit string and qualities string if this instruction simulates from reverse strand.
                if (!inst.isForward)
                {
                    reverse(inst.editString);
                    reverse(inst.qualities);
                }
                // Apply simulation instructions.
                SEQAN_ASSERT_EQ(length(fragmentStore.readSeqStore), length(fragmentStore.readNameStore));
                // Cut out segment from haplotype.
                String<Dna5Q> read = infix(haplotypeContigsCopy[inst.contigId], inst.beginPos, inst.endPos);
                String<Dna5Q> haplotypeInfix = read;  // Copy for printing later on.
                applySimulationInstructions(read, rng, inst, options);
                if (!inst.isForward)
                {
                    reverseComplement(read);
                    // Reconstruct edit string and qualities in sequencing direction.
                    reverse(inst.editString);
                    reverse(inst.qualities);
                }
                // Append read sequence to read seq store and mate pair to read name store.  This also yields the read
                // id.  We will generate and append the read name below, depending on the read id.
                unsigned readId;
                if (options.generateMatePairs)
                    readId = appendRead(fragmentStore, read, length(fragmentStore.matePairStore));
                else
                    readId = appendRead(fragmentStore, read);

                // Get expected begin/end position in the original sequence.
                TPos origBeginPos = virtualToHostPosition(haplotypeContigs[inst.contigId], inst.beginPos);
                TPos origEndPos = virtualToHostPosition(haplotypeContigs[inst.contigId], inst.endPos);

                // Generate read name.
                // TODO(holtgrew): Remove mateNum, not required?
                if (options.generateMatePairs)
                {
                    // Generate the mate num \in {1, 2}, randomly but consistent so two entries belonging together have
                    // different nums.
                    if (options.includeReadInformation)
                    {
                        resize(readNameBuf, 1024 + length(haplotypeInfix) + length(inst.editString) + length(inst.bsConversionString));
                        sprintf(&readNameBuf[0], "%s.%09u contig=%s haplotype=%u length=%lu orig_begin=%lu orig_end=%lu haplotype_infix=%s edit_string=", outFileName, readId / 2, toCString(fragmentStore.contigNameStore[inst.contigId]), haplotypeId, static_cast<long unsigned>(length(read)), static_cast<long unsigned>(origBeginPos), static_cast<long unsigned>(origEndPos), toCString(CharString(haplotypeInfix)));
                    }
                    else
                    {
                        resize(readNameBuf, 1024);
                        sprintf(&readNameBuf[0], "%s.%09u", outFileName, readId / 2);
                    }
                } else {
                    if (options.includeReadInformation)
                    {
                        resize(readNameBuf, 1024 + length(haplotypeInfix) + length(inst.editString) + length(inst.bsConversionString));
                        sprintf(&readNameBuf[0], "%s.%09u contig=%s haplotype=%u length=%lu orig_begin=%lu orig_end=%lu haplotype_infix=%s edit_string=", outFileName, readId, toCString(fragmentStore.contigNameStore[inst.contigId]), haplotypeId, static_cast<long unsigned>(length(read)), static_cast<long unsigned>(origBeginPos), static_cast<long unsigned>(origEndPos), toCString(CharString(haplotypeInfix)));
                    }
                    else
                    {
                        resize(readNameBuf, 1024);
                        sprintf(&readNameBuf[0], "%s.%09u", outFileName, readId);
                    }
                }
                if (options.includeReadInformation) {
                    for (unsigned i = 0; i < length(inst.editString); ++i) {
                        char buffer[2] = "*";
                        buffer[0] = "MEID"[static_cast<int>(inst.editString[i])];
                        strcat(&readNameBuf[0], buffer);
                    }
                    // bs_change:
                    strcat(&readNameBuf[0], " bs_conversion_string=");
                    for (unsigned i = 0; i < length(inst.bsConversionString); ++i) {
                        char buffer[2] = "*";
                        buffer[0] = "01"[static_cast<int>(inst.bsConversionString[i])];
                        strcat(&readNameBuf[0], buffer);
                    }

                }
                // bs_change:
                CharString readName(&readNameBuf[0]);
                if (inst.isFromTopStrand)
                    append(readName, " original_strand=top");
                else
                    append(readName, " original_strand=bottom");

                if (inst.isForward)
                    append(readName, " strand=forward");
                else
                    append(readName, " strand=reverse");
                appendValue(fragmentStore.readNameStore, readName);

                // Print info about read and haplotype.
                if (options.veryVerbose) {
                    std::cout << ",-- Read #" << readId << std::endl
                              << "| inst.beginPos    " << inst.beginPos << std::endl
                              << "| inst.endPos      " << inst.endPos << std::endl
                              << "| origBeginPos     " << origBeginPos << std::endl
                              << "| origEndPos       " << origEndPos << std::endl
                              << "| isgapinhost      " << isGapInHost(haplotypeContigs[inst.contigId], inst.beginPos-1) << std::endl
                              << "| isgapinhost      " << isGapInHost(haplotypeContigs[inst.contigId], inst.beginPos) << std::endl
                              << "| isgapinhost      " << isGapInHost(haplotypeContigs[inst.contigId], inst.beginPos+1) << std::endl
                              << "| name:            " << readName << std::endl
                              << "| original infix:  " << infix(fragmentStore.contigStore[inst.contigId].seq, origBeginPos, origEndPos) << std::endl
                              << "| haplotype infix: " << infix(haplotypeContigsCopy[inst.contigId], inst.beginPos, inst.endPos) << std::endl
                              << "| read:            " << read << std::endl
                              << "`-- " << std::endl;
                }

                // Swap original begin and end position if from reverse strand.
                if (!inst.isForward)
                    std::swap(origBeginPos, origEndPos);

                // Add matches to aligned read store.
                if (options.generateMatePairs)
                    appendAlignedRead(fragmentStore, readId, inst.contigId, origBeginPos, origEndPos, length(fragmentStore.matePairStore));
                else
                    appendAlignedRead(fragmentStore, readId, inst.contigId, origBeginPos, origEndPos);

                // Adding mate pair information.
                if (options.generateMatePairs && readId % 2 == 1) {  // Only append mate pair info after simulating second mate.
                    // Append mate pair element to fragment store's mate pair store.
                    TMatePairStoreElement matePair;
                    matePair.readId[0] = readId - 1;
                    matePair.readId[1] = readId;
                    appendValue(fragmentStore.matePairStore, matePair);
                }
            }
			if (options.generateMatePairs)
            {
                // When generating mate pairs, an even number of reads is generated in each step.
 				SEQAN_ASSERT_EQ(length(fragmentStore.alignedReadStore) % 2, 0u);
				SEQAN_ASSERT_EQ(length(fragmentStore.readNameStore) % 2, 0u);
				SEQAN_ASSERT_EQ(length(fragmentStore.readSeqStore) % 2, 0u);
			}
		}
    }

    // Last but not least, convert the matches collected before to a global alignment.
    convertMatchesToGlobalAlignment(fragmentStore, Score<int, EditDistance>(), True());
	
	// AlignedReadLayout layout;
	// layoutAlignment(layout, fragmentStore);
	// printAlignment(std::cout, Raw(), layout, fragmentStore, 0, 0, 300, 0, 100);
    
    if (options.verbose)
        std::cerr << "Simulated " << length(fragmentStore.readSeqStore) << " reads" << std::endl;

    return 0;
}



#endif  // #ifndef SANDBOX_KRAKAU_APPS_BS_MASON_SIMULATE_ILLUMINA_BS_H_

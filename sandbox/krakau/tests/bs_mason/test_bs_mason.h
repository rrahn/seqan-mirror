// ==========================================================================
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

#ifndef SANDBOX_KRAKAU_TESTS_BS_MASON_TEST_BS_MASON_H_
#define SANDBOX_KRAKAU_TESTS_BS_MASON_TEST_BS_MASON_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/misc/misc_cmdparser.h>

//CharString path = "/home/takifugu/sabrina7/";
CharString path = "/Volumes/sabrina7/";


void testSetUpCG_simulateMethylPositions(FragmentStore<MyFragmentStoreConfig> & fragmentStore, Options<IlluminaReadsBS> & options, String<bool> & methylString1, String<bool> & methylString2){
    // Set methylation probabilities
    options.probabilityMethylationCG = 1.0;
    options.probabilityMethylationCHG = 0.0;
    options.probabilityMethylationCHH = 0.0;

    // Only for 1 contig
    resize(fragmentStore.contigStore, 1);   

    // "AAAANAAAANAAAANAAAANAAAAN"
    //
    String<Dna5Q> contig = "AGTTNACGANCTAANACTGNACTAN";
    assign(fragmentStore.contigStore[0].seq, contig); 

    clear(methylString1);
    resize(methylString1, length(contig), false);
    assignValue(methylString1, 6, true);
    clear(methylString2);
    resize(methylString2, length(contig), false, Exact());
    assignValue(methylString2, 7, true);

    SEQAN_ASSERT_EQ(length(methylString1), length(contig));
    SEQAN_ASSERT_EQ(length(methylString2), length(contig));
    SEQAN_ASSERT_EQ(length(fragmentStore.contigStore[0].seq), length(contig));
 }

void testSetUpCHG_simulateMethylPositions(FragmentStore<MyFragmentStoreConfig> & fragmentStore, Options<IlluminaReadsBS> & options, String<bool> & methylString1, String<bool> & methylString2){
    // Set methylation probabilities
    options.probabilityMethylationCG = 0.0;
    options.probabilityMethylationCHG = 1.0;
    options.probabilityMethylationCHH = 0.0;

    // Only for 1 contig
    resize(fragmentStore.contigStore, 1);   

    // "AAAANAAAANAAAANAAAANAAAAN"
    //
    String<Dna5Q> contig = "ACGANACTTNAAAANACAGNACATN";    
    assign(fragmentStore.contigStore[0].seq, contig); 

    clear(methylString1);
    resize(methylString1, length(contig), false, Exact());
    assignValue(methylString1, 16, true);
    clear(methylString2);
    resize(methylString2, length(contig), false, Exact());
    assignValue(methylString2, 18, true);
    
    SEQAN_ASSERT_EQ(length(methylString1), length(contig));
    SEQAN_ASSERT_EQ(length(methylString2), length(contig));
    SEQAN_ASSERT_EQ(length(fragmentStore.contigStore[0].seq), length(contig));
}

void testSetUpCHH_simulateMethylPositions(FragmentStore<MyFragmentStoreConfig> & fragmentStore, Options<IlluminaReadsBS> & options, String<bool> & methylString1, String<bool> & methylString2){
    // Set methylation probabilities
    options.probabilityMethylationCG = 0.0;
    options.probabilityMethylationCHG = 0.0;
    options.probabilityMethylationCHH = 1.0;

    // Only for 1 contig
    resize(fragmentStore.contigStore, 1);   

    // "AAAANAAAANAAAANAAAANAAAAN"
    //
    String<Dna5Q> contig = "ACGANACTTNAGAANACAGNACATN";
    assign(fragmentStore.contigStore[0].seq, contig); 

    clear(methylString1);
    resize(methylString1, length(contig), false, Exact());
    assignValue(methylString1, 6, true);
    assignValue(methylString1, 21, true);
    clear(methylString2);
    resize(methylString2, length(contig), false, Exact());
    assignValue(methylString2, 11, true);
    
    SEQAN_ASSERT_EQ(length(methylString1), length(contig));
    SEQAN_ASSERT_EQ(length(methylString2), length(contig));
    SEQAN_ASSERT_EQ(length(fragmentStore.contigStore[0].seq), length(contig));
}


// A test for strings.
SEQAN_DEFINE_TEST(test_bs_mason_simulateMethylPositions)
{
    using namespace seqan;

    Options<IlluminaReadsBS> options;
    ModelParameters<IlluminaReadsBS> parameters;
    FragmentStore<MyFragmentStoreConfig> fragmentStore;
    Rng<MersenneTwister> rng(options.seed);
    String<bool> methylString1;
    String<bool> methylString2;

    //
    // Only CG context: set methylation probability in CG context to 1.0, others to 0.0
    //
    // Setup contig in fragmentsStore, methylation probabilities in options and build correct methylation strings for both strands
    testSetUpCG_simulateMethylPositions(fragmentStore, options, methylString1, methylString2);
    simulateMethylPositions(parameters, rng, options, fragmentStore);
    SEQAN_ASSERT_EQ(parameters.topMethylPositions[0], methylString1);
    SEQAN_ASSERT_EQ(parameters.bottomMethylPositions[0], methylString2);

    //
    // Only CHG context: set methylation probability in CHG context to 1.0, others to 0.0
    //
    testSetUpCHG_simulateMethylPositions(fragmentStore, options, methylString1, methylString2);
    simulateMethylPositions(parameters, rng, options, fragmentStore);
    SEQAN_ASSERT_EQ(parameters.topMethylPositions[0], methylString1);
    SEQAN_ASSERT_EQ(parameters.bottomMethylPositions[0], methylString2);

    //
    // Only CHH context: set methylation probability in CHH context to 1.0, others to 0.0
    //
    testSetUpCHH_simulateMethylPositions(fragmentStore, options, methylString1, methylString2);
    simulateMethylPositions(parameters, rng, options, fragmentStore);
    /*for (unsigned i = 0; i < length(methylString1); ++i){
        std::cout << "i: " << i << "   " << getValue(parameters.topMethylPositions[0] , i)  << "  :  " << getValue(methylString1 , i) << std::endl;
    }*/
    SEQAN_ASSERT_EQ(parameters.topMethylPositions[0], methylString1);
    SEQAN_ASSERT_EQ(parameters.bottomMethylPositions[0], methylString2);
}

void testSetUp_buildHaplotype(FragmentStore<MyFragmentStoreConfig> & fragmentStore, Options<IlluminaReadsBS> & options, ModelParameters<IlluminaReadsBS> & parameters){
   
    resize(fragmentStore.contigStore, 2, Exact());   
    clear(parameters.topMethylPositions);
    clear(parameters.bottomMethylPositions);
    resize(parameters.topMethylPositions, 2, Exact());
    resize(parameters.bottomMethylPositions, 2, Exact());
    // Create haplotype parameters
    assign(options.haplotypeSnpRate, 0.3);
    assign(options.haplotypeIndelRate, 0.3);
    assign(options.haplotypeIndelRangeMin, 1); 
    assign(options.haplotypeIndelRangeMax, 6);
    //
    // Contig1: contains no Cs -> no methyl positions at all
    //
    // Create original reference contig
    String<Dna5Q> contig = "AAAANAAAANAAAANAAAANAAAAN";
    assign(fragmentStore.contigStore[0].seq, contig); 
    // Create methylation positions for original reference (all c positions in this case)
    resize(parameters.topMethylPositions[0], length(contig), false, Exact());
    resize(parameters.bottomMethylPositions[0], length(contig), false, Exact());
    //
    // Contig2: contains Cs -> check if methylations are only at methyl positions
    //
    // Create original reference contig
    contig = "CCCCNGGGGNCCCCNAACANACGANACACNCCCCNAAGCNAACANCCCCNGGGGNCCCCN";
    assign(fragmentStore.contigStore[1].seq, contig); 
    // Create methylation positions for original reference (all c positions in this case)
    resize(parameters.topMethylPositions[1], length(contig), false, Exact());
    resize(parameters.bottomMethylPositions[1], length(contig), false, Exact());
    for(unsigned i = 0; i < length(contig); ++i){
        // ... true at all positions with C in top strand
        if (contig[i] == 'C'){
            assignValue(parameters.topMethylPositions[1], i, true);
        }
        // ... true at all positions with C in bottom strand
        if (contig[i] == 'G'){
            assignValue(parameters.bottomMethylPositions[1], i, true);
        }
    }

}

// Check for case that all original Cs are methylated
void check_buildHaplotype(StringSet<String<Dna5, Journaled<Alloc<> > > > & haplotypeContigs, ModelParameters<IlluminaReadsBS> & parameters){
    //
    // Contig1: contains no Cs -> no methyl positions at all
    //
    SEQAN_ASSERT_EQ(length(haplotypeContigs[0]), length(parameters.topMethylPositionsHaplotype[0]));
    SEQAN_ASSERT_EQ(length(haplotypeContigs[0]), length(parameters.bottomMethylPositionsHaplotype[0]));
 
    for (unsigned i = 0; i < length(haplotypeContigs[0]); ++i){
        // Check if no methyl positions
        SEQAN_ASSERT_NOT(getValue(parameters.topMethylPositionsHaplotype[0], i));
        SEQAN_ASSERT_NOT(getValue(parameters.bottomMethylPositionsHaplotype[0], i));
     }  
    /*std::cout << "haplo    :" << haplotypeContigs[0] << std::endl;
    std::cout << "top      :";
    for (unsigned i = 0; i < length(parameters.topMethylPositionsHaplotype[0]); ++i){
        std::cout << getValue(parameters.topMethylPositionsHaplotype[0], i);     
    }
    std::cout << std::endl;
    std::cout << "bottom   :";
    for (unsigned i = 0; i < length(parameters.bottomMethylPositionsHaplotype[0]); ++i){
        std::cout << getValue(parameters.bottomMethylPositionsHaplotype[0], i);     
    }
    std::cout << std::endl;*/
    //
    // Contig2: contains Cs -> check if methylations are only at methyl positions
    //
    SEQAN_ASSERT_EQ(length(haplotypeContigs[1]), length(parameters.topMethylPositionsHaplotype[1]));
    SEQAN_ASSERT_EQ(length(haplotypeContigs[1]), length(parameters.bottomMethylPositionsHaplotype[1]));
 
    for (unsigned i = 0; i < length(haplotypeContigs[1]); ++i){
        // Check if topMethylPositions is only at C positions true
        if (getValue(haplotypeContigs[1], i) != 'C')
            SEQAN_ASSERT_NOT(getValue(parameters.topMethylPositionsHaplotype[1], i));
        // Check if bottomMethylPositions is only at C positions true
        if (getValue(haplotypeContigs[1], i) != 'G')
            SEQAN_ASSERT_NOT(getValue(parameters.bottomMethylPositionsHaplotype[1], i));
     }  
}


SEQAN_DEFINE_TEST(test_bs_mason_buildHaplotype)
{
    using namespace seqan;

    Options<IlluminaReadsBS> options;
    ModelParameters<IlluminaReadsBS> parameters;
    FragmentStore<MyFragmentStoreConfig> fragmentStore;
    Rng<MersenneTwister> rng(options.seed);
    StringSet<String<Dna5, Journaled<Alloc<> > > > haplotypeContigs;
    // 
    // Check case: all original Cs are methylated
    //
    // Setup contig in fragmentStore, haplotyeSnpRate etc. in options  and methylation positions in parameters  
    testSetUp_buildHaplotype(fragmentStore, options, parameters);
    buildHaplotype(haplotypeContigs, parameters, fragmentStore, rng, options);
    // Check if haplotype methyl positions suit to haplotype C positions
    check_buildHaplotype(haplotypeContigs, parameters);

}


void testSetUp_buildMethylHaplotypeUseMethylPositions(FragmentStore<MyFragmentStoreConfig> & fragmentStore, ModelParameters<IlluminaReadsBS> & parameters){
   
    resize(fragmentStore.contigStore, 2, Exact());   
    clear(parameters.topMethylPositions);
    clear(parameters.bottomMethylPositions);
    resize(parameters.topMethylPositions, 2, Exact());
    resize(parameters.bottomMethylPositions, 2, Exact());
    clear(parameters.topMethylPositionsHaplotype);
    clear(parameters.bottomMethylPositionsHaplotype);
    resize(parameters.topMethylPositionsHaplotype, 2, Exact());
    resize(parameters.bottomMethylPositionsHaplotype, 2, Exact());
    //
    // Contig1: contains no Cs -> no methyl positions at all
    //
    // Create original reference contig
    String<Dna5Q> contig = "AAAANAAAANAAAANAAAANAAAAN";
    assign(fragmentStore.contigStore[0].seq, contig); 
    // Create methylation positions for original reference (all c positions in this case)
    resize(parameters.topMethylPositions[0], length(contig), false, Exact());
    resize(parameters.bottomMethylPositions[0], length(contig), false, Exact());
    clear(parameters.topMethylPositionsHaplotype[0]);
    clear(parameters.bottomMethylPositionsHaplotype[0]);
    setHost(parameters.topMethylPositionsHaplotype[0], parameters.topMethylPositions[0]);
    setHost(parameters.bottomMethylPositionsHaplotype[0], parameters.bottomMethylPositions[0]);
    //
    // Contig2: contains Cs -> check if methylations are only at methyl positions
    //
    // Create original reference contig
    contig = "CCCCNGGGGNCCCCNAACANACGANACACNCCCCNAAGCNAACANCCCCNGGGGNCCCCN";
    assign(fragmentStore.contigStore[1].seq, contig); 
    // Create methylation positions for original reference (all c positions in this case)
    resize(parameters.topMethylPositions[1], length(contig), false, Exact());
    resize(parameters.bottomMethylPositions[1], length(contig), false, Exact());
    clear(parameters.topMethylPositionsHaplotype[1]);
    clear(parameters.bottomMethylPositionsHaplotype[1]);
    setHost(parameters.topMethylPositionsHaplotype[1], parameters.topMethylPositions[1]);
    setHost(parameters.bottomMethylPositionsHaplotype[1], parameters.bottomMethylPositions[1]);   
    for(unsigned i = 0; i < length(contig); ++i){
        // ... true at all positions with C in top strand
        if (contig[i] == 'C'){
            assignValue(parameters.topMethylPositions[1], i, true);
        }
        // ... true at all positions with C in bottom strand
        if (contig[i] == 'G'){
            assignValue(parameters.bottomMethylPositions[1], i, true);
        }
    }
}


SEQAN_DEFINE_TEST(test_bs_mason_buildMethylHaplotypeUseMethylPositions)
{
    using namespace seqan;

    Options<IlluminaReadsBS> options;
    ModelParameters<IlluminaReadsBS> parameters;
    FragmentStore<MyFragmentStoreConfig> fragmentStore;
    Rng<MersenneTwister> rng(options.seed);
    // 1)
    testSetUp_buildMethylHaplotypeUseMethylPositions(fragmentStore, parameters);
    options.haplotypeMethylInsRate = 0.0;
    options.haplotypeMethylDelRate = 1.0;
    buildMethylHaplotypeUseMethylPositions(parameters, fragmentStore, rng, options);
    // Contig1: contains no Cs -> no methyl positions at all
    SEQAN_ASSERT_EQ(length(fragmentStore.contigStore[0].seq), length(parameters.topMethylPositionsHaplotype[0]));
    SEQAN_ASSERT_EQ(length(fragmentStore.contigStore[0].seq), length(parameters.bottomMethylPositionsHaplotype[0]));
    for (unsigned i = 0; i < length(fragmentStore.contigStore[0].seq); ++i){
        // Check if no methyl positions
        SEQAN_ASSERT_NOT(getValue(parameters.topMethylPositionsHaplotype[0], i));
        SEQAN_ASSERT_NOT(getValue(parameters.bottomMethylPositionsHaplotype[0], i));
     }  
    // Contig2: contains Cs -> check if no methyl positions
    SEQAN_ASSERT_EQ(length(fragmentStore.contigStore[1].seq), length(parameters.topMethylPositionsHaplotype[1]));
    SEQAN_ASSERT_EQ(length(fragmentStore.contigStore[1].seq), length(parameters.bottomMethylPositionsHaplotype[1]));
    for (unsigned i = 0; i < length(fragmentStore.contigStore[1].seq); ++i){
        // Check if no methyl positions
        SEQAN_ASSERT_NOT(getValue(parameters.topMethylPositionsHaplotype[1], i));
        SEQAN_ASSERT_NOT(getValue(parameters.bottomMethylPositionsHaplotype[1], i));
    }
    // 2)
    options.haplotypeMethylInsRate = 1.0;
    options.haplotypeMethylDelRate = 0.0;
    buildMethylHaplotypeUseMethylPositions(parameters, fragmentStore, rng, options);
    // Contig1: contains no Cs -> no methyl positions at all
    SEQAN_ASSERT_EQ(length(fragmentStore.contigStore[0].seq), length(parameters.topMethylPositionsHaplotype[0]));
    SEQAN_ASSERT_EQ(length(fragmentStore.contigStore[0].seq), length(parameters.bottomMethylPositionsHaplotype[0]));
    for (unsigned i = 0; i < length(fragmentStore.contigStore[0].seq); ++i){
        // Check if no methyl positions
        SEQAN_ASSERT_NOT(getValue(parameters.topMethylPositionsHaplotype[0], i));
        SEQAN_ASSERT_NOT(getValue(parameters.bottomMethylPositionsHaplotype[0], i));
     }  
    // Contig2: contains Cs -> check if all Cs are methylated
    SEQAN_ASSERT_EQ(length(fragmentStore.contigStore[1].seq), length(parameters.topMethylPositionsHaplotype[1]));
    SEQAN_ASSERT_EQ(length(fragmentStore.contigStore[1].seq), length(parameters.bottomMethylPositionsHaplotype[1]));
    for (unsigned i = 0; i < length(fragmentStore.contigStore[1].seq); ++i){
        // Check if topMethylPositions is (only) true at C positions
        if (getValue(fragmentStore.contigStore[1].seq, i) != 'C')
            SEQAN_ASSERT_NOT(getValue(parameters.topMethylPositionsHaplotype[1], i));
        else
            SEQAN_ASSERT(getValue(parameters.topMethylPositionsHaplotype[1], i));
        // Check if bottomMethylPositions is (only) at C positions true
        if (getValue(fragmentStore.contigStore[1].seq, i) != 'G')
            SEQAN_ASSERT_NOT(getValue(parameters.bottomMethylPositionsHaplotype[1], i));
        else
            SEQAN_ASSERT(getValue(parameters.bottomMethylPositionsHaplotype[1], i));
     }  

}

void testSetUp_buildSimulationInstructions(String<bool> & check_bsConversionString, ReadSimulationInstruction<IlluminaReadsBS> & inst, unsigned & readLength, String<Dna5> & contig, ModelParameters<IlluminaReadsBS> & parameters, Options<IlluminaReadsBS> & options){
    
    assign(contig, "AAAANAAAANAAAANAAAANAAAANAAAANAAAAN");
    assign(readLength, length(contig));
    assign(inst.contigId, 0);

    options.probabilityInsert = 0.0;
    options.probabilityDelete = 0.0;
    options.conversionRate = 1.0;

    clear(parameters.mismatchProbabilities);
    resize(parameters.mismatchProbabilities, readLength, 0.0, Exact());

    if (inst.isFromTopStrand){
        // Create c and methylation patterns for top strand
        clear(parameters.topMethylPositions);
        resize(parameters.topMethylPositions, 1, Exact());
        clear(parameters.topMethylPositionsHaplotype);
        resize(parameters.topMethylPositionsHaplotype, 1, Exact());
        resize(parameters.topMethylPositions[0], length(contig), false, Exact());    
        setHost(parameters.topMethylPositionsHaplotype[0], parameters.topMethylPositions[0]);

        // Default: no methylation -> all Cs converted to Ts
        assignValue(contig, 1, 'C');     
        assignValue(contig, 3, 'C');
        assignValue(contig, 7, 'C');
        assignValue(contig, 8, 'C');
        assignValue(contig, 21, 'C');
        assignValue(contig, 23, 'C');
        assignValue(contig, 33, 'C');    

        // Methylation: no C-T conversion
        assignValue(parameters.topMethylPositionsHaplotype[0], 1, true);
        assignValue(parameters.topMethylPositionsHaplotype[0], 33, true);

        inst.beginPos = 0;
        inst.endPos = inst.beginPos + readLength;

        clear(check_bsConversionString);
        resize(check_bsConversionString, readLength, true, Exact());
        assignValue(check_bsConversionString, 1, false);
        assignValue(check_bsConversionString, 33, false); 
    } else {
        // Create c and methylation patterns for bottom strand
        clear(parameters.bottomMethylPositions);
        resize(parameters.bottomMethylPositions, 1, Exact());
        clear(parameters.bottomMethylPositionsHaplotype);
        resize(parameters.bottomMethylPositionsHaplotype, 1, Exact());
        resize(parameters.bottomMethylPositions[0], length(contig), false, Exact());    
        setHost(parameters.bottomMethylPositionsHaplotype[0], parameters.bottomMethylPositions[0]);

        assignValue(contig, 12, 'G');
        assignValue(contig, 15, 'G');   
        assignValue(contig, 18, 'G');   
        assignValue(contig, 26, 'G');
        assignValue(contig, 28, 'G');
        assignValue(contig, 31, 'G');
        assignValue(contig, 32, 'G');

        // Methylation: no G-A conversion
        assignValue(parameters.bottomMethylPositionsHaplotype[0], 15, true);
        assignValue(parameters.bottomMethylPositionsHaplotype[0], 28, true);

        inst.beginPos = 0;
        inst.endPos = inst.beginPos + readLength;

        clear(check_bsConversionString);
        resize(check_bsConversionString, readLength, true, Exact());
        assignValue(check_bsConversionString, 15, false);
        assignValue(check_bsConversionString, 28, false);    }

}

SEQAN_DEFINE_TEST(test_bs_mason_buildSimulationInstructions)
{
    using namespace seqan;
    
    ModelParameters<IlluminaReadsBS> parameters;
    Options<IlluminaReadsBS> options;
    ReadSimulationInstruction<IlluminaReadsBS> inst;
    Rng<MersenneTwister> rng(options.seed);
    unsigned readLength;
    String<Dna5> contig;
    String<bool> check_bsConversionString;
    
    // From top strand
    inst.isFromTopStrand = true;
    testSetUp_buildSimulationInstructions(check_bsConversionString, inst, readLength, contig, parameters, options);
    buildSimulationInstructions(inst, rng, readLength, contig, parameters, options);
    SEQAN_ASSERT_EQ(inst.bsConversionString, check_bsConversionString);

    // From bottom strand:
    inst.isFromTopStrand = false;
    testSetUp_buildSimulationInstructions(check_bsConversionString, inst, readLength, contig, parameters, options);
    buildSimulationInstructions(inst, rng, readLength, contig, parameters, options);
    /*std::cout << "inst.bsConversionString   :";
    for (unsigned i = 0; i < length(inst.bsConversionString ); ++i){
        std::cout << inst.bsConversionString[i];     
    }
    std::cout << std::endl;
    std::cout << "check_bsConversionString  :";
    for (unsigned i = 0; i < length(check_bsConversionString); ++i){
        std::cout << check_bsConversionString[i];     
    }
    std::cout << std::endl;*/
    SEQAN_ASSERT_EQ(inst.bsConversionString, check_bsConversionString);
}


void testSetUp_applySimulationInstructions(String<Dna5Q> & check_read, String<Dna5Q> & read, ReadSimulationInstruction<IlluminaReadsBS> & inst, Options<IlluminaReadsBS> & options){
    
    assign(inst.insCount, 0);
    assign(inst.delCount, 0);

    read = "AACANACGANACAGNAACANAGAANAGGANAAAAN";
    assign(options.readLength, length(read));

    clear(inst.bsConversionString);
    resize(inst.bsConversionString, length(read), true, Exact());

    if (inst.isFromTopStrand)  
        check_read = "AATANATGANATAGNAATANAGAANAGGANAAAAN";
    else 
        check_read = "AACANACAANACAANAACANAAAANAAAANAAAAN";
    
    clear(inst.editString);
    resize(inst.editString, length(read),  ERROR_TYPE_MATCH, Exact());
}

SEQAN_DEFINE_TEST(test_bs_mason_applySimulationInstructions)
{
    using namespace seqan;

    Options<IlluminaReadsBS> options;
    ReadSimulationInstruction<IlluminaReadsBS> inst;
    Rng<MersenneTwister> rng(options.seed);
    String<Dna5Q> read;
    String<Dna5Q> check_read;

    // Top strand:
    assign(inst.isForward, true);
    assign(inst.isFromTopStrand, true);
    testSetUp_applySimulationInstructions(check_read, read, inst, options);

    applySimulationInstructions(read, rng, inst, options); 
    SEQAN_ASSERT_EQ(read, check_read);

    // Bottom strand:
    assign(inst.isForward, true);
    assign(inst.isFromTopStrand, false);
    testSetUp_applySimulationInstructions(check_read, read, inst, options);

    applySimulationInstructions(read, rng, inst, options); 
    SEQAN_ASSERT_EQ(read, check_read);
}


void testCheck_checkOutput(CharString & outputFilename){
    
    typedef Stream<std::fstream> TStream;
    typedef RecordReader<std::fstream, SinglePass<> > TRecordReader;
    std::fstream file(toCString(outputFilename), std::ios::binary | std::ios::in);
    TRecordReader reader(file);

    CharString str;
    CharString beginPos_Str;
    CharString endPos_Str;
    size_t beginPos;
    size_t endPos;
    Dna5String haplotypeInfix;
    CharString editString;
    CharString bsConversionString;
    CharString hString;
    bool isFromTopStrand;
    bool isForward;

    Dna5String simulatedRead;

    while(!atEnd(reader))
    {
        clear(str);
        clear(beginPos_Str);
        clear(endPos_Str);
        clear(haplotypeInfix);
        clear(editString);
        clear(bsConversionString);
        clear(hString);
        clear(simulatedRead);
        //
        // Read input line
        //
        readNChars(str, reader, 1);
        if (str[0] != '>'){
            skipLine(reader);
            std::cerr << "ERROR: Expected '>', got instead " << str[0] << std::endl;
            continue;
        }        
        skipUntilWhitespace(reader);    // readname
        skipUntilWhitespace(reader);    // contig
        skipUntilWhitespace(reader);    // haplotype
        skipUntilWhitespace(reader);    // length
        // beginPos       
        skipWhitespaces(reader);
        skipUntilChar(reader, '=');
        readDigits(beginPos_Str, reader);
        lexicalCast2(beginPos, beginPos_Str); 
        // endPos
        skipWhitespaces(reader);
        skipUntilChar(reader, '=');
        readDigits(endPos_Str, reader);
        lexicalCast2(endPos, endPos_Str); 
        // haplotypeInfix
        skipWhitespaces(reader);
        skipUntilChar(reader, '=');
        readUntilWhitespace(haplotypeInfix, reader);
        // editString
        skipWhitespaces(reader);
        skipUntilChar(reader, '=');
        readUntilWhitespace(editString, reader);
        // bsConversionString
        skipWhitespaces(reader);
        skipUntilChar(reader, '=');
        readUntilWhitespace(bsConversionString, reader);
        // Get "top" or "bottom"
        skipWhitespaces(reader);
        skipUntilChar(reader, '=');
        readUntilWhitespace(hString, reader);
        if (hString == "top")
            isFromTopStrand = true;
        else if (hString == "bottom")
            isFromTopStrand = false;
        // Get "forward" or "reverse"
        skipWhitespaces(reader);
        skipUntilChar(reader, '=');
        readUntilWhitespace(hString, reader);
        if (hString == "forward")
            isForward = true;
        else if (hString == "reverse")
            isForward = false;
        skipLine(reader);
        // simulatedRead
        readUntilWhitespace(simulatedRead, reader);
        skipLine(reader);

        // 
        // Check simulatedRead
        //
        // Build revCompl of simulatedRead if reverse -> we just have to check forward directions, easy comparison to editString etc.
        if (!isForward)
            reverseComplement(simulatedRead);
        // For each position:
        // if bsConversionString indicates conversion: 
        //  check if 'C' or 'G' in original sequence (haplotypeInfix)
        // if no error (match) and if bsConversionString indicates no conversion: 
        //  check if simulatedReads is equal to haplotypeInfix at current position
        unsigned b = 0; // position in bsConversionString
        unsigned i = 0; // position in haplotype infix
        unsigned j = 0; // position in simulated read
        for (unsigned e = 0; e < length(editString); ++e){
            if (editString[e] == 'M'){
                // MATCH
                if (bsConversionString[b] == '0')
                    SEQAN_ASSERT_EQ(simulatedRead[j], haplotypeInfix[i]);
                else {
                    if (isFromTopStrand){
                        SEQAN_ASSERT_EQ(haplotypeInfix[i], 'C');
                        SEQAN_ASSERT_EQ(simulatedRead[j], 'T');    
                    } else {
                        SEQAN_ASSERT_EQ(haplotypeInfix[i], 'G');
                        SEQAN_ASSERT_EQ(simulatedRead[j], 'A'); 
                    }
                }
                ++b;
                ++i;
                ++j;
            } else if (editString[e] == 'E'){
                // MISMATCH
                if (bsConversionString[b] == '1'){
                    if (isFromTopStrand)
                        SEQAN_ASSERT_EQ(haplotypeInfix[i], 'C'); 
                    else
                        SEQAN_ASSERT_EQ(haplotypeInfix[i], 'G'); 
                }
                ++b;
                ++i;
                ++j;
            } else if (editString[e] == 'I'){
                // INSERTION
                ++j;
            } else if (editString[e] == 'D'){
                // DELETION
                if (bsConversionString[b] == '1'){
                    if (isFromTopStrand)
                        SEQAN_ASSERT_EQ(haplotypeInfix[i], 'C'); 
                    else
                        SEQAN_ASSERT_EQ(haplotypeInfix[i], 'G'); 
                }
                ++b;
                ++i;
            }
        }
    }
}
// Simulate methylation rates
SEQAN_DEFINE_TEST(test_bs_mason_checkOutputSimMR)
{
    using namespace seqan;

    CharString referenceFilename = path;
    append(referenceFilename, "Development/workspace_seqan/seqan-trunk/sandbox/krakau/apps/bs_mason/tests/adeno-genome.fa");
    CharString outputFilename = path;
    append(outputFilename, "Development/workspace_seqan/seqan-trunk/sandbox/krakau/apps/bs_mason/tests/adeno-genome.fa.fasta");

    Options<IlluminaReadsBS> options;
    // 1)
    options.includeReadInformation = true; 
    options.numReads = 1000;
    simulateReads(options, referenceFilename, IlluminaReadsBS());      
    testCheck_checkOutput(outputFilename);

    // 2)
    options.includeReadInformation = true; 
    options.numHaplotypes = 2;
    simulateReads(options, referenceFilename, IlluminaReadsBS());      
    testCheck_checkOutput(outputFilename);

    // 3)
    options.includeReadInformation = true; 
    options.probabilityMismatch = 0.1;
    options.probabilityMismatchBegin = 0.05;
    options.probabilityMismatchEnd = 0.2;
    simulateReads(options, referenceFilename, IlluminaReadsBS());      
    testCheck_checkOutput(outputFilename);

}

// check output additional for matepairs
void testCheck_checkOutput_PE(CharString & outputFilename){
    
    typedef Stream<std::fstream> TStream;
    typedef RecordReader<std::fstream, SinglePass<> > TRecordReader;
    CharString filename = outputFilename;
    append(filename, "_1");
    std::fstream file1(toCString(filename), std::ios::binary | std::ios::in);
    TRecordReader reader1(file1);
    clear(filename);
    append(filename, "_2");
    std::fstream file2(toCString(filename), std::ios::binary | std::ios::in);
    TRecordReader reader2(file2);

    CharString str;
    CharString beginPos_Str;
    CharString endPos_Str;
    size_t beginPos;
    size_t endPos;
    Dna5String haplotypeInfix;
    CharString editString;
    CharString bsConversionString;
    CharString hString;
    bool isFromTopStrand1;
    bool isForward1;
    bool isFromTopStrand2;
    bool isForward2;

    Dna5String simulatedRead;

    while(!atEnd(reader1) && !atEnd(reader2))
    {
        //////////////////////////////////////////////////////////////////
        // mate 1
        clear(str);
        clear(beginPos_Str);
        clear(endPos_Str);
        clear(haplotypeInfix);
        clear(editString);
        clear(bsConversionString);
        clear(hString);
        clear(simulatedRead);
        //
        // Read input line
        //
        readNChars(str, reader1, 1);
        if (str[0] != '>'){
            skipLine(reader1);
            std::cerr << "ERROR: Expected '>', got instead " << str[0] << std::endl;
            continue;
        }        
        skipUntilWhitespace(reader1);    // readname
        skipUntilWhitespace(reader1);    // contig
        skipUntilWhitespace(reader1);    // haplotype
        skipUntilWhitespace(reader1);    // length
        // beginPos       
        skipWhitespaces(reader1);
        skipUntilChar(reader1, '=');
        readDigits(beginPos_Str, reader1);
        lexicalCast2(beginPos, beginPos_Str); 
        // endPos
        skipWhitespaces(reader1);
        skipUntilChar(reader1, '=');
        readDigits(endPos_Str, reader1);
        lexicalCast2(endPos, endPos_Str); 
        // haplotypeInfix
        skipWhitespaces(reader1);
        skipUntilChar(reader1, '=');
        readUntilWhitespace(haplotypeInfix, reader1);
        // editString
        skipWhitespaces(reader1);
        skipUntilChar(reader1, '=');
        readUntilWhitespace(editString, reader1);
        // bsConversionString
        skipWhitespaces(reader1);
        skipUntilChar(reader1, '=');
        readUntilWhitespace(bsConversionString, reader1);
        // Get "top" or "bottom"
        skipWhitespaces(reader1);
        skipUntilChar(reader1, '=');
        readUntilWhitespace(hString, reader1);
        if (hString == "top")
            isFromTopStrand1 = true;
        else if (hString == "bottom")
            isFromTopStrand1 = false;
        // Get "forward" or "reverse"
        skipWhitespaces(reader1);
        skipUntilChar(reader1, '=');
        readUntilWhitespace(hString, reader1);
        if (hString == "forward")
            isForward1 = true;
        else if (hString == "reverse")
            isForward1 = false;
        skipLine(reader1);
        // simulatedRead
        readUntilWhitespace(simulatedRead, reader1);
        skipLine(reader1);

        //////////////////////////////////////////////////////////////////
        // mate 2
        clear(str);
        clear(beginPos_Str);
        clear(endPos_Str);
        clear(haplotypeInfix);
        clear(editString);
        clear(bsConversionString);
        clear(hString);
        clear(simulatedRead);
        //
        // Read input line
        //
        readNChars(str, reader2, 1);
        if (str[0] != '>'){
            skipLine(reader2);
            std::cerr << "ERROR: Expected '>', got instead " << str[0] << std::endl;
            continue;
        }        
        skipUntilWhitespace(reader2);    // readname
        skipUntilWhitespace(reader2);    // contig
        skipUntilWhitespace(reader2);    // haplotype
        skipUntilWhitespace(reader2);    // length
        // beginPos       
        skipWhitespaces(reader2);
        skipUntilChar(reader2, '=');
        readDigits(beginPos_Str, reader2);
        lexicalCast2(beginPos, beginPos_Str); 
        // endPos
        skipWhitespaces(reader2);
        skipUntilChar(reader2, '=');
        readDigits(endPos_Str, reader2);
        lexicalCast2(endPos, endPos_Str); 
        // haplotypeInfix
        skipWhitespaces(reader2);
        skipUntilChar(reader2, '=');
        readUntilWhitespace(haplotypeInfix, reader2);
        // editString
        skipWhitespaces(reader2);
        skipUntilChar(reader2, '=');
        readUntilWhitespace(editString, reader2);
        // bsConversionString
        skipWhitespaces(reader2);
        skipUntilChar(reader2, '=');
        readUntilWhitespace(bsConversionString, reader2);
        // Get "top" or "bottom"
        skipWhitespaces(reader2);
        skipUntilChar(reader2, '=');
        readUntilWhitespace(hString, reader2);
        if (hString == "top")
            isFromTopStrand2 = true;
        else if (hString == "bottom")
            isFromTopStrand2 = false;
        // Get "forward" or "reverse"
        skipWhitespaces(reader2);
        skipUntilChar(reader2, '=');
        readUntilWhitespace(hString, reader2);
        if (hString == "forward")
            isForward2 = true;
        else if (hString == "reverse")
            isForward2 = false;
        skipLine(reader2);
        // simulatedRead
        readUntilWhitespace(simulatedRead, reader2);
        skipLine(reader2);

        /////////////////////////////////////////////////////////
        // check orientation of matepairs
        
        SEQAN_ASSERT_EQ(isFromTopStrand1, isFromTopStrand2);
        SEQAN_ASSERT_NEQ(isForward1, isForward2);
    }
}
// Simulate methylation rates for matepairs
SEQAN_DEFINE_TEST(test_bs_mason_checkOutputSimMR_PE)
{
    using namespace seqan;

    CharString referenceFilename = path;
    append(referenceFilename, "Development/workspace_seqan/seqan-trunk/sandbox/krakau/apps/bs_mason/tests/adeno-genome.fa");
    CharString outputFilename = path;
    append(outputFilename, "Development/workspace_seqan/seqan-trunk/sandbox/krakau/apps/bs_mason/tests/adeno-genome.fa.fasta");

    Options<IlluminaReadsBS> options;
    // 1)
    options.includeReadInformation = true; 
    options.numReads = 1000;
    options.generateMatePairs = true;
    simulateReads(options, referenceFilename, IlluminaReadsBS());      
    testCheck_checkOutput_PE(outputFilename);
}

////////////////////////////////////////////////////////////////////////////////////////
// Use given methylation rates
////////////////////////////////////////////////////////////////////////////////////////


SEQAN_DEFINE_TEST(test_bs_mason_buildMethylHaplotypeUseMethylRates){

    using namespace seqan;

    ModelParameters<IlluminaReadsBS> parameters;
    Options<IlluminaReadsBS> options;
    Rng<MersenneTwister> rng(options.seed);

    clear(parameters.topMethylRates);
    clear(parameters.bottomMethylRates);
    resize(parameters.topMethylRates, 2, Exact());
    resize(parameters.bottomMethylRates, 2, Exact());

    resize(parameters.topMethylRates[0], 10, 0.0, Exact());
    resize(parameters.bottomMethylRates[0], 10, 0.0, Exact());
    resize(parameters.topMethylRates[1], 15, 0.0, Exact());
    resize(parameters.bottomMethylRates[1], 15, 0.0, Exact());

    value(parameters.topMethylRates[0], 2) = 1.0;
    value(parameters.topMethylRates[0], 5) = 1.0;
    value(parameters.topMethylRates[0], 9) = 1.0;
    value(parameters.bottomMethylRates[0], 3) = 1.0;
    value(parameters.bottomMethylRates[0], 8) = 1.0;

    value(parameters.topMethylRates[1], 1) = 1.0;
    value(parameters.topMethylRates[1], 2) = 1.0;
    value(parameters.bottomMethylRates[1], 11) = 1.0;
    value(parameters.bottomMethylRates[1], 12) = 1.0;

    buildMethylHaplotypeUseMethylRates(parameters,rng);

    SEQAN_ASSERT_EQ(length(parameters.topMethylRates), length(parameters.topMethylPositionsHaplotype));
    SEQAN_ASSERT_EQ(length(parameters.bottomMethylRates), length(parameters.bottomMethylPositionsHaplotype));
    SEQAN_ASSERT_EQ(length(parameters.topMethylRates), length(parameters.bottomMethylRates));
    for (unsigned i = 0; i < length(parameters.topMethylRates); ++i){
        SEQAN_ASSERT_EQ(length(parameters.topMethylRates[i]), length(parameters.topMethylPositionsHaplotype[i]));
        for (unsigned j = 0; j < length(parameters.topMethylRates[i]); ++j){
            if (getValue(parameters.topMethylRates[i], j) == 1.0 )
               SEQAN_ASSERT(getValue(parameters.topMethylPositionsHaplotype[i], j));
            else
               SEQAN_ASSERT_NOT(getValue(parameters.topMethylPositionsHaplotype[i], j));       
            if (getValue(parameters.bottomMethylRates[i], j) == 1.0 )
               SEQAN_ASSERT(getValue(parameters.bottomMethylPositionsHaplotype[i], j));
            else
               SEQAN_ASSERT_NOT(getValue(parameters.bottomMethylPositionsHaplotype[i], j));         
        }
    }
}

SEQAN_DEFINE_TEST(test_bs_mason_checkOutputUseMR)
{
    using namespace seqan;

    CharString referenceFilename = path;
    append(referenceFilename, "RAID/human_methylome/h1/h1_c_basecalls/simReads_21_22/hg18_21_22_small.fa");
    CharString outputFilename = path;
    append(outputFilename, "RAID/human_methylome/h1/h1_c_basecalls/simReads_21_22_small/reads.fasta");
    CharString mrFilename = path;
    append(mrFilename, "RAID/human_methylome/h1/h1_c_basecalls/simReads_21_22/h1_21_22_methLevel_s2.txt");
    
    Options<IlluminaReadsBS> options;
    options.includeReadInformation = true; 
    options.numReads = 1000;
    options.useMethylRates = true;
    options.methylRatesFile = mrFilename;
    simulateReads(options, referenceFilename, IlluminaReadsBS());      
    testCheck_checkOutput(outputFilename);
}

// check output additional for matepair orientations
SEQAN_DEFINE_TEST(test_bs_mason_checkOutputUseMR_PE)
{
    using namespace seqan;

    CharString referenceFilename = path;
    append(referenceFilename, "RAID/human_methylome/h1/h1_c_basecalls/simReads_21_22/hg18_21_22_small.fa");
    CharString outputFilename = path;
    append(outputFilename, "RAID/human_methylome/h1/h1_c_basecalls/simReads_21_22_small/reads.fasta");
    CharString mrFilename = path;
    append(mrFilename, "RAID/human_methylome/h1/h1_c_basecalls/simReads_21_22/h1_21_22_methLevel_s2.txt");
    
    Options<IlluminaReadsBS> options;
    options.includeReadInformation = true; 
    options.numReads = 1000;
    options.useMethylRates = true;
    options.methylRatesFile = mrFilename;
    options.generateMatePairs = true;
    simulateReads(options, referenceFilename, IlluminaReadsBS());      
    testCheck_checkOutput_PE(outputFilename);
}


#endif  // SANDBOX_KRAKAU_TESTS_BS_MASON_TEST_BS_MASON_H_

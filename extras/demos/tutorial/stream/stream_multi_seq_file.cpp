// FRAGMENT(includes)
#include <seqan/file.h>
#include <iostream>

using namespace seqan;

// FRAGMENT(open_file)
int main (int argc, char const ** argv)
{
// FRAGMENT(open-guess-split)
    MultiSeqFile multiSeqFile;
    if (argc < 2 || !open(multiSeqFile.concat, argv[1], OPEN_RDONLY))
        return 1;

    AutoSeqFormat format;
    guessFormat(multiSeqFile.concat, format);
    split(multiSeqFile, format);

// FRAGMENT(load)
    unsigned seqCount = length(multiSeqFile);
    StringSet<String<Dna5Q> > seqs;
    StringSet<CharString> seqIDs;

    reserve(seqs, seqCount, Exact());
    reserve(seqIDs, seqCount, Exact());

    String<Dna5Q> seq;
    CharString qual;
    CharString id;

// FRAGMENT(output)
    for (unsigned i = 0; i < seqCount; ++i)
    {
        assignSeq(seq, multiSeqFile[i], format);    // read sequence
        assignQual(qual, multiSeqFile[i], format);  // read ascii quality values
        assignSeqId(id, multiSeqFile[i], format);   // read sequence id

        // Convert ascii to values from 0..62.
        // We store DNA and quality together in Dna5Q.
        for (unsigned j = 0; j < length(qual) && j < length(seq); ++j)
            assignQualityValue(seq[j], static_cast<int>(ordValue(qual[j]) - 33));

        // We use reserve and append, as assign is not supported
        // by StringSet<..., Owner<ConcatDirect<> > >
        appendValue(seqs, seq, Generous());
        appendValue(seqIDs, id, Generous());
    }
// FRAGMENT(return)
    return 0;
}

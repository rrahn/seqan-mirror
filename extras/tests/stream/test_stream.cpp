// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Tests for the stream module.
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/stream.h>

#include "test_stream_char_array.h"
#if SEQAN_HAS_ZLIB
#include "test_stream_gz_file.h"
#include "test_stream_bgzf.h"
#endif  // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
#include "test_stream_bz2_file.h"
#endif  // #if SEQAN_HAS_BZIP2
#include "test_stream_adapt_cstdio.h"
#include "test_stream_adapt_fstream.h"
#include "test_stream_adapt_sstream.h"
#include "test_stream_adapt_mmap.h"
#include "test_stream_tokenize.h"
#include "test_stream_lexical_cast.h"
#include "test_stream_record_reader_fasta.h"
#include "test_stream_record_reader_fastq.h"
#include "test_stream_guess_stream_format.h"
#include "test_stream_write_fasta.h"
#include "test_stream_read_embl.h"
#include "test_stream_read_genbank.h"
#include "test_stream_read_auto_format.h"

SEQAN_BEGIN_TESTSUITE(test_stream)
{
    // Tests for Char Array Stream.
    SEQAN_CALL_TEST(test_stream_char_array_metafunctions);
    SEQAN_CALL_TEST(test_stream_char_array_read_simple_usage);
    SEQAN_CALL_TEST(test_stream_char_array_read_complex_usage);
    SEQAN_CALL_TEST(test_stream_char_array_write_simple_usage);
    SEQAN_CALL_TEST(test_stream_char_array_write_complex_usage);
    SEQAN_CALL_TEST(test_stream_char_array_eof);
    SEQAN_CALL_TEST(test_stream_char_array_peek);
    SEQAN_CALL_TEST(test_stream_char_array_read_char);
    SEQAN_CALL_TEST(test_stream_char_array_read_block);
    SEQAN_CALL_TEST(test_stream_char_array_write_block);
    SEQAN_CALL_TEST(test_stream_char_array_streamPut);
    SEQAN_CALL_TEST(test_stream_char_array_write_char);
    SEQAN_CALL_TEST(test_stream_char_array_flush);
    SEQAN_CALL_TEST(test_stream_char_array_seek);
    SEQAN_CALL_TEST(test_stream_char_array_tell);

#if SEQAN_HAS_ZLIB  // Enable tests for Stream<GZFile>, Stream<Bgzf> if available.
    // Tests for BZ2 File Stream.
    SEQAN_CALL_TEST(test_stream_gz_file_metafunctions);
    SEQAN_CALL_TEST(test_stream_gz_file_read_simple_usage);
    SEQAN_CALL_TEST(test_stream_gz_file_read_complex_usage);
    SEQAN_CALL_TEST(test_stream_gz_file_write_simple_usage);
    SEQAN_CALL_TEST(test_stream_gz_file_write_complex_usage);
    SEQAN_CALL_TEST(test_stream_gz_file_eof);
    SEQAN_CALL_TEST(test_stream_gz_file_peek);
    SEQAN_CALL_TEST(test_stream_gz_file_read_char);
    SEQAN_CALL_TEST(test_stream_gz_file_read_block);
    SEQAN_CALL_TEST(test_stream_gz_file_write_block);
    SEQAN_CALL_TEST(test_stream_gz_file_streamPut);
    SEQAN_CALL_TEST(test_stream_gz_file_write_char);
    SEQAN_CALL_TEST(test_stream_gz_file_flush);
    SEQAN_CALL_TEST(test_stream_gz_file_seek);
    SEQAN_CALL_TEST(test_stream_gz_file_tell);

    // Test Stream<Bgzf>.
    SEQAN_CALL_TEST(test_stream_bgzf_metafunctions);
    SEQAN_CALL_TEST(test_stream_bgzf_read_simple_usage);
    SEQAN_CALL_TEST(test_stream_bgzf_read_complex_usage);
    SEQAN_CALL_TEST(test_stream_bgzf_write_simple_usage);
    SEQAN_CALL_TEST(test_stream_bgzf_write_complex_usage);
    SEQAN_CALL_TEST(test_stream_bgzf_eof);
    SEQAN_CALL_TEST(test_stream_bgzf_peek);
    SEQAN_CALL_TEST(test_stream_bgzf_read_char);
    SEQAN_CALL_TEST(test_stream_bgzf_read_block);
    SEQAN_CALL_TEST(test_stream_bgzf_write_block);
    SEQAN_CALL_TEST(test_stream_bgzf_streamPut);
    SEQAN_CALL_TEST(test_stream_bgzf_write_char);
    SEQAN_CALL_TEST(test_stream_bgzf_flush);
    SEQAN_CALL_TEST(test_stream_bgzf_seek);
    SEQAN_CALL_TEST(test_stream_bgzf_tell);

    SEQAN_CALL_TEST(test_stream_bgzf_write_large_and_compare_with_file);
    SEQAN_CALL_TEST(test_stream_bgzf_from_file_and_compare);
#endif  // #if SEQAN_HAS_ZLIB

#if SEQAN_HAS_BZIP2  // Enable tests for Stream<BZ2File> if available.
    // Tests for BZ2 File Stream.
    SEQAN_CALL_TEST(test_stream_bz2_file_metafunctions);
    SEQAN_CALL_TEST(test_stream_bz2_file_read_simple_usage);
    SEQAN_CALL_TEST(test_stream_bz2_file_read_complex_usage);
    SEQAN_CALL_TEST(test_stream_bz2_file_write_simple_usage);
    SEQAN_CALL_TEST(test_stream_bz2_file_write_complex_usage);
    SEQAN_CALL_TEST(test_stream_bz2_file_eof);
    SEQAN_CALL_TEST(test_stream_bz2_file_read_char);
    SEQAN_CALL_TEST(test_stream_bz2_file_read_block);
    SEQAN_CALL_TEST(test_stream_bz2_file_write_block);
    SEQAN_CALL_TEST(test_stream_bz2_file_write_char);
    SEQAN_CALL_TEST(test_stream_bz2_file_streamPut);
    SEQAN_CALL_TEST(test_stream_bz2_file_flush);
#endif  // #if SEQAN_HAS_BZIP2

    // Tests for cstdio.
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_metafunctions);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_read_simple_usage);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_read_complex_usage);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_write_simple_usage);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_write_complex_usage);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_eof);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_peek);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_read_char);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_read_block);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_write_block);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_write_char);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_streamPut);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_flush);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_seek);
    SEQAN_CALL_TEST(test_stream_adapt_cstdio_tell);

    // Tests for std::stringstream adaptions.
    SEQAN_CALL_TEST(test_stream_adapt_sstream_metafunctions);
    SEQAN_CALL_TEST(test_stream_adapt_sstream_read_simple_usage);
    SEQAN_CALL_TEST(test_stream_adapt_sstream_read_complex_usage);
    SEQAN_CALL_TEST(test_stream_adapt_sstream_write_simple_usage);
    SEQAN_CALL_TEST(test_stream_adapt_sstream_write_complex_usage);
    SEQAN_CALL_TEST(test_stream_adapt_sstream_eof);
    SEQAN_CALL_TEST(test_stream_adapt_sstream_peek);
    SEQAN_CALL_TEST(test_stream_adapt_sstream_read_char);
    SEQAN_CALL_TEST(test_stream_adapt_sstream_read_block);
    SEQAN_CALL_TEST(test_stream_adapt_sstream_write_char);
    SEQAN_CALL_TEST(test_stream_adapt_sstream_write_block);
    SEQAN_CALL_TEST(test_stream_adapt_sstream_streamPut);
    SEQAN_CALL_TEST(test_stream_adapt_sstream_flush);
    SEQAN_CALL_TEST(test_stream_adapt_sstream_seek);
    SEQAN_CALL_TEST(test_stream_adapt_sstream_tell);

    // Tests for std::istringstream adaptions.
    SEQAN_CALL_TEST(test_stream_adapt_istringstream_metafunctions);
    SEQAN_CALL_TEST(test_stream_adapt_istringstream_read_simple_usage);
    SEQAN_CALL_TEST(test_stream_adapt_istringstream_read_complex_usage);
    SEQAN_CALL_TEST(test_stream_adapt_istringstream_eof);
    SEQAN_CALL_TEST(test_stream_adapt_istringstream_peek);
    SEQAN_CALL_TEST(test_stream_adapt_istringstream_read_char);
    SEQAN_CALL_TEST(test_stream_adapt_istringstream_read_block);
    SEQAN_CALL_TEST(test_stream_adapt_istringstream_seek);
    SEQAN_CALL_TEST(test_stream_adapt_istringstream_tell);

    // Tests for std::ostringstream adaptions.
    SEQAN_CALL_TEST(test_stream_adapt_ostringstream_metafunctions);
    SEQAN_CALL_TEST(test_stream_adapt_ostringstream_write_simple_usage);
    SEQAN_CALL_TEST(test_stream_adapt_ostringstream_write_complex_usage);
    SEQAN_CALL_TEST(test_stream_adapt_ostringstream_write_char);
    SEQAN_CALL_TEST(test_stream_adapt_ostringstream_write_block);
    SEQAN_CALL_TEST(test_stream_adapt_ostringstream_streamPut);
    SEQAN_CALL_TEST(test_stream_adapt_ostringstream_flush);
    SEQAN_CALL_TEST(test_stream_adapt_ostringstream_seek);
    SEQAN_CALL_TEST(test_stream_adapt_ostringstream_tell);

    // Tests for std::fstream adaptions.
    SEQAN_CALL_TEST(test_stream_adapt_fstream_metafunctions);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_read_simple_usage);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_read_complex_usage);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_write_simple_usage);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_write_complex_usage);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_eof);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_peek);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_read_char);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_read_block);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_write_block);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_write_char);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_streamPut);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_flush);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_seek);
    SEQAN_CALL_TEST(test_stream_adapt_fstream_tell);

    // Tests for std::ifstream adaptions.
    SEQAN_CALL_TEST(test_stream_adapt_ifstream_metafunctions);
    SEQAN_CALL_TEST(test_stream_adapt_ifstream_read_simple_usage);
    SEQAN_CALL_TEST(test_stream_adapt_ifstream_read_complex_usage);
    SEQAN_CALL_TEST(test_stream_adapt_ifstream_eof);
    SEQAN_CALL_TEST(test_stream_adapt_ifstream_peek);
    SEQAN_CALL_TEST(test_stream_adapt_ifstream_read_char);
    SEQAN_CALL_TEST(test_stream_adapt_ifstream_read_block);
    SEQAN_CALL_TEST(test_stream_adapt_ifstream_seek);
    SEQAN_CALL_TEST(test_stream_adapt_ifstream_tell);

    // Tests for std::ofstream adaptions.
    SEQAN_CALL_TEST(test_stream_adapt_ofstream_metafunctions);
    SEQAN_CALL_TEST(test_stream_adapt_ofstream_write_simple_usage);
    SEQAN_CALL_TEST(test_stream_adapt_ofstream_write_complex_usage);
    SEQAN_CALL_TEST(test_stream_adapt_ofstream_write_block);
    SEQAN_CALL_TEST(test_stream_adapt_ofstream_write_char);
    SEQAN_CALL_TEST(test_stream_adapt_ofstream_streamPut);
    SEQAN_CALL_TEST(test_stream_adapt_ofstream_flush);

    // Tests for mmap-stream adaptation
    SEQAN_CALL_TEST(test_stream_adapt_mmap_streamPut);

    // Tests for tokenize.h
    SEQAN_CALL_TEST(test_stream_tokenizing_readUntil);
    SEQAN_CALL_TEST(test_stream_tokenizing_readNChars);
    SEQAN_CALL_TEST(test_stream_tokenizing_readIgnoring);
    SEQAN_CALL_TEST(test_stream_tokenizing_readLine);
    SEQAN_CALL_TEST(test_stream_tokenizing_skipUntil);
    SEQAN_CALL_TEST(test_stream_tokenizing_skipWhile);
    SEQAN_CALL_TEST(test_stream_tokenizing_skipLine);
    SEQAN_CALL_TEST(test_stream_tokenizing_skipUntilString);
    SEQAN_CALL_TEST(test_stream_tokenizing_skipUntilLineBeginsWithChar);
    SEQAN_CALL_TEST(test_stream_tokenizing_skipUntilLineBeginsWithStr);
    SEQAN_CALL_TEST(test_stream_tokenizing_skipUntilLineBeginsWithOneCharOfStr);
    SEQAN_CALL_TEST(test_stream_tokenizing_read_until_tab_or_line_break);
    SEQAN_CALL_TEST(test_stream_tokenizing_read_until_one_of);

    // Tests for lexical_cast
    SEQAN_CALL_TEST(test_stream_lexical_cast_1_stdstring);
    SEQAN_CALL_TEST(test_stream_lexical_cast_1_chararray);
    SEQAN_CALL_TEST(test_stream_lexical_cast_1_seqanstring);
    SEQAN_CALL_TEST(test_stream_lexical_cast_2_stdstring);
    SEQAN_CALL_TEST(test_stream_lexical_cast_2_chararray);
    SEQAN_CALL_TEST(test_stream_lexical_cast_2_seqanstring);

    // -------------- File format specific code ------------------

    /* RecordReader is tested for each Format with fstream and mmap-backend.
     * Every backend is tested with single-pass record-reading, those that have
     * double-pass implementations are also tested for double-pass
     * record-reading and batch reading (double-pass whole document). For the
     * last case an additional test for concat-direct-strings is performed. */

    // Tests for FASTA
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_single_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_double_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_batch_single_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_batch_double_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_single_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_double_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_batch_single_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_batch_single_concat_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_batch_double_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_batch_double_concat_mmap);

    // Tests for FASTQ
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_single_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_double_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_batch_single_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_batch_double_fstream);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_single_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_double_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_batch_single_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_batch_single_concat_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_batch_double_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_batch_double_concat_mmap);
    SEQAN_CALL_TEST(test_stream_record_reader_fastq_check_stream_format);

    // Tests for EMBL
    SEQAN_CALL_TEST(test_stream_read_embl_single_char_array_stream);
    SEQAN_CALL_TEST(test_stream_read_embl_record_char_array_stream);
    SEQAN_CALL_TEST(test_stream_read_embl_batch_char_array_stream);
    SEQAN_CALL_TEST(test_stream_read_embl_single_mmap);
    SEQAN_CALL_TEST(test_stream_read_embl_single_batch_mmap);
    SEQAN_CALL_TEST(test_stream_read_embl_single_batch_concat_mmap);

    // Tests for GenBank
    SEQAN_CALL_TEST(test_stream_read_genbank_single_char_array_stream);
    SEQAN_CALL_TEST(test_stream_read_genbank_record_char_array_stream);
    SEQAN_CALL_TEST(test_stream_read_genbank_batch_char_array_stream);
    SEQAN_CALL_TEST(test_stream_read_genbank_single_mmap);
    SEQAN_CALL_TEST(test_stream_read_genbank_single_batch_mmap);
    SEQAN_CALL_TEST(test_stream_read_genbank_single_batch_concat_mmap);

    // TODO Tests for other formats once they are supported

    // Tests for file format auto-detection
    SEQAN_CALL_TEST(test_stream_guess_stream_format_auto_fasta);
    SEQAN_CALL_TEST(test_stream_guess_stream_format_auto_fastq);
    SEQAN_CALL_TEST(test_stream_guess_stream_format_auto_bogus);

    // Tests for reading with automatic file format detection.
    SEQAN_CALL_TEST(test_stream_read_record_auto_format_quals_fasta);
    SEQAN_CALL_TEST(test_stream_read_record_auto_format_quals_fastq);
    SEQAN_CALL_TEST(test_stream_read_record_auto_format_no_quals_fasta);
    SEQAN_CALL_TEST(test_stream_read_record_auto_format_no_quals_fastq);

    SEQAN_CALL_TEST(test_stream_read_auto_format_quals_fasta);
    SEQAN_CALL_TEST(test_stream_read_auto_format_quals_fastq);
    SEQAN_CALL_TEST(test_stream_read_auto_format_no_quals_fasta);
    SEQAN_CALL_TEST(test_stream_read_auto_format_no_quals_fastq);

    // Tests for FASTA with Amino Acid alphabet.
    SEQAN_CALL_TEST(test_stream_record_reader_fasta_protein_single_fstream);

    /* Tests for writing file formats.
     * Only fstream used as a representative */

    // Tests for FASTA
    SEQAN_CALL_TEST(test_stream_write_record_fasta_default);
    SEQAN_CALL_TEST(test_stream_write_record_fasta_nolinebreaks);
    SEQAN_CALL_TEST(test_stream_write_record_fastq_default_separate_qual);
    SEQAN_CALL_TEST(test_stream_write_record_fastq_default_qual_in_seq);
    SEQAN_CALL_TEST(test_stream_write_record_fastq_linebreaks_qualmeta);
    SEQAN_CALL_TEST(test_stream_write2_fasta_default);
    SEQAN_CALL_TEST(test_stream_write2_fastq_default_separate_qual);
    SEQAN_CALL_TEST(test_stream_write2_fastq_default_qual_in_seq);
}
SEQAN_END_TESTSUITE


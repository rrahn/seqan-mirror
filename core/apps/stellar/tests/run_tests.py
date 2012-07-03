#!/usr/bin/env python
"""Execute the tests for stellar.

The golden test outputs are generated by the script generate_outputs.sh.

You have to give the root paths to the source and the binaries as arguments to
the program.  These are the paths to the directory that contains the 'projects'
directory.

Usage:  run_tests.py SOURCE_ROOT_PATH BINARY_ROOT_PATH
"""
import logging
import os
import os.path
import sys

# Automagically add util/py_lib to PYTHONPATH environment variable.
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..',
                                    '..', '..', 'util', 'py_lib'))
sys.path.insert(0, path)

import seqan.app_tests as app_tests

def main(source_base, binary_base):
    """Main entry point of the script."""

    print 'Executing test for stellar'
    print '========================='
    print
    
    ph = app_tests.TestPathHelper(
        source_base, binary_base,
        'core/apps/stellar/tests')  # tests dir

    # ============================================================
    # Auto-detect the binary path.
    # ============================================================

    path_to_program = app_tests.autolocateBinary(
      binary_base, 'core/apps/stellar', 'stellar')

    # ============================================================
    # Built TestConf list.
    # ============================================================

    # Build list with TestConf objects, analoguely to how the output
    # was generated in generate_outputs.sh.
    conf_list = []

    # We prepare a list of transforms to apply to the output files.  This is
    # used to strip the input/output paths from the programs' output to
    # make it more canonical and host independent.
    ph.outFile('-')  # To ensure that the out path is set.
    transforms = [
        app_tests.ReplaceTransform(os.path.join(ph.source_base_path, 'core/apps/stellar/tests') + os.sep, '', right=True),
        app_tests.ReplaceTransform(ph.temp_dir + os.sep, '', right=True),
        app_tests.NormalizeScientificExponentsTransform(),
        ]

    # ============================================================
    # Run STELLAR.
    # ============================================================

    # Error rate 0.1:
    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('e-1.stdout'),
        args=['-d', ph.inFile('512_simSeq1_e-1.fa'),
              '-q', ph.inFile('512_simSeq2_e-1.fa'),
              '-e', '0.1', '-l', '50', '-x', '10', '-k', '7', '-n', '5000',
              '-s', '10000', '-f', '-v', '-t', '-of', 'gff',
              '-o', ph.outFile('e-1.gff')],
        to_diff=[(ph.inFile('e-1.stdout'),
                  ph.outFile('e-1.stdout'),
                  transforms),
                 (ph.inFile('e-1.gff'),
                  ph.outFile('e-1.gff'),
                  transforms)])
    conf_list.append(conf)

    # Error rate 0.05:
    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('5e-2.stdout'),
        args=['-d', ph.inFile('512_simSeq1_5e-2.fa'),
              '-q', ph.inFile('512_simSeq2_5e-2.fa'),
              '-e', '0.05', '-l', '50', '-x', '10', '-k', '7', '-n', '5000',
              '-s', '10000', '-f', '-v', '-t', '-of', 'gff',
              '-o', ph.outFile('5e-2.gff')],
        to_diff=[(ph.inFile('5e-2.stdout'),
                  ph.outFile('5e-2.stdout'),
                  transforms),
                 (ph.inFile('5e-2.gff'),
                  ph.outFile('5e-2.gff'),
                  transforms)])
    conf_list.append(conf)

    # Error rate 0.25:
    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('25e-3.stdout'),
        args=['-d', ph.inFile('512_simSeq1_25e-3.fa'),
              '-q', ph.inFile('512_simSeq2_25e-3.fa'),
              '-e', '0.025', '-l', '50', '-x', '10', '-k', '7', '-n', '5000',
              '-s', '10000', '-f', '-v', '-t', '-of', 'gff',
              '-o', ph.outFile('25e-3.gff')],
        to_diff=[(ph.inFile('25e-3.stdout'),
                  ph.outFile('25e-3.stdout'),
                  transforms),
                 (ph.inFile('25e-3.gff'),
                  ph.outFile('25e-3.gff'),
                  transforms)])
    conf_list.append(conf)

    # Error rate 0.75:
    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('75e-3.stdout'),
        args=['-d', ph.inFile('512_simSeq1_75e-3.fa'),
              '-q', ph.inFile('512_simSeq2_75e-3.fa'),
              '-e', '0.075', '-l', '50', '-x', '10', '-k', '7', '-n', '5000',
              '-s', '10000', '-f', '-v', '-t', '-of', 'gff',
              '-o', ph.outFile('75e-3.gff')],
        to_diff=[(ph.inFile('75e-3.stdout'),
                  ph.outFile('75e-3.stdout'),
                  transforms),
                 (ph.inFile('75e-3.gff'),
                  ph.outFile('75e-3.gff'),
                  transforms)])
    conf_list.append(conf)

    # Error rate 0.0001:
    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('e-4.stdout'),
        args=['-d', ph.inFile('512_simSeq1_e-4.fa'),
              '-q', ph.inFile('512_simSeq2_e-4.fa'),
              '-e', '0.0001', '-l', '50', '-x', '10', '-k', '7', '-n', '5000',
              '-s', '10000', '-f', '-v', '-t', '-of', 'gff',
              '-o', ph.outFile('e-4.gff')],
        to_diff=[(ph.inFile('e-4.stdout'),
                  ph.outFile('e-4.stdout'),
                  transforms),
                 (ph.inFile('e-4.gff'),
                  ph.outFile('e-4.gff'),
                  transforms)])
    conf_list.append(conf)

    # Minimal length: 20, Error rate 0.05:
    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('minLen20.stdout'),
        args=['-d', ph.inFile('512_simSeq1_5e-2.fa'),
              '-q', ph.inFile('512_simSeq2_5e-2.fa'), '-e', '0.05', '-l', '20',
              '-x', '10', '-k', '7', '-n', '5000', '-s', '10000', '-f', '-v',
              '-t', '-of', 'gff',
              '-o', ph.outFile('minLen20.gff')],
        to_diff=[(ph.inFile('minLen20.stdout'),
                  ph.outFile('minLen20.stdout'),
                  transforms),
                 (ph.inFile('minLen20.gff'),
                  ph.outFile('minLen20.gff'),
                  transforms)])
    conf_list.append(conf)

    # Minimal length: 150, Error rate 0.05:
    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('minLen150.stdout'),
        args=['-d', ph.inFile('512_simSeq1_5e-2.fa'),
              '-q', ph.inFile('512_simSeq2_5e-2.fa'),
              '-e', '0.05', '-l', '150', '-x', '10', '-k', '7', '-n', '5000',
              '-s', '10000', '-f', '-v', '-t', '-of', 'gff',
              '-o', ph.outFile('minLen150.gff')],
        to_diff=[(ph.inFile('minLen150.stdout'),
                  ph.outFile('minLen150.stdout'),
                  transforms),
                 (ph.inFile('minLen150.gff'),
                  ph.outFile('minLen150.gff'),
                  transforms)])
    conf_list.append(conf)

    # Output format text:
    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('5e-2txt.stdout'),
        args=['-d', ph.inFile('512_simSeq1_5e-2.fa'),
              '-q', ph.inFile('512_simSeq2_5e-2.fa'),
              '-e', '0.05', '-l', '50', '-x', '10', '-k', '7', '-n', '5000',
              '-s', '10000', '-f', '-v', '-t', '-of', 'text',
              '-o', ph.outFile('5e-2.txt')],
        to_diff=[(ph.inFile('5e-2txt.stdout'),
                  ph.outFile('5e-2txt.stdout'),
                  transforms),
                 (ph.inFile('5e-2.txt'),
                  ph.outFile('5e-2.txt'),
                  transforms)])
    conf_list.append(conf)

    # ============================================================
    # Execute the tests.
    # ============================================================
    failures = 0
    for conf in conf_list:
        print ' '.join(['stellar'] + conf.args),
        sys.stdout.flush()
        res = app_tests.runTest(conf)
        # Output to the user.
        print ' '.join(['stellar'] + conf.args),
        if res:
             print 'OK'
        else:
            failures += 1
            print 'FAILED'

    print '=============================='
    print '     total tests: %d' % len(conf_list)
    print '    failed tests: %d' % failures
    print 'successful tests: %d' % (len(conf_list) - failures)
    print '=============================='
    # Compute and return return code.
    return failures != 0


if __name__ == '__main__':
    sys.exit(app_tests.main(main))

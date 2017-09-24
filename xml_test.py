#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""A test script for XML importing development.
"""

# Ensure Python 2 and 3 compatibility
from __future__ import division
from __future__ import print_function

import os

from xena_dataset import read_clinical, read_biospecimen

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    test_dir = os.path.join(script_dir, 'gitignore', 'test')
    os.chdir(test_dir)
    
#    test_file = os.path.join(
#            test_dir, 'nationwidechildrens.org_biospecimen.TCGA-A7-A13E.xml'
#        )
#    with open(test_file, 'r') as f:
#        test = read_biospecimen(f)
#    test.to_csv('biospecimen.test.tsv', sep='\t')
#    return

    test_file = os.path.join(
            test_dir, 'nationwidechildrens.org_clinical.TCGA-A7-A13E.xml'
        )
    with open(test_file, 'r') as f:
        test = read_clinical(f)
    test.T.to_csv('clinical.test.tsv', sep='\t')
    return

if __name__ == '__main__':
    main()

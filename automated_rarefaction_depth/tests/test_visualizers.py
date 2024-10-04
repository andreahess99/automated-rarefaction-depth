import os.path

from skbio.alignment import TabularMSA
from skbio.sequence import DNA

from qiime2.plugin.testing import TestPluginBase

from q2_dwq2._visualizers import summarize_alignment

#example test
class SummarizeAlignmentTests(TestPluginBase):
    package = 'q2_dwq2.tests'

    def test_simple1(self):
        aligned_sequence1 = DNA('AAAAAAAAGGTGGCCTTTTTTTT')
        aligned_sequence2 = DNA('AAAAAAAAGG-GGCCTTTTTTTT')
        msa = TabularMSA([aligned_sequence1, aligned_sequence2])

        with self.temp_dir as output_dir:
            summarize_alignment(output_dir, msa)

            with open(os.path.join(output_dir, 'index.html'), 'r') as fh:
                observed = fh.read()

            self.assertIn('AAAAAAAAGGTGGCCTTTTTTTT', observed)
            self.assertIn('AAAAAAAAGG-GGCCTTTTTTTT', observed)

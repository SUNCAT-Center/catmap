from __future__ import absolute_import

import unittest

import numpy as np
import os
import shutil

from catmap.model import ReactionModel
from catmap.analyze.mechanism import MechanismAnalysis

ref_dir = os.path.join('tutorials', '3-refining_microkinetic_model')

class MechanismTest(unittest.TestCase):
    def setUp(self):
        shutil.copy(os.path.join(ref_dir, 'CO_oxidation.mkm'), '.')
        shutil.copy(os.path.join(ref_dir, 'energies.txt'), '.')
        self.rxn_model = ReactionModel(setup_file='CO_oxidation.mkm') 
        self.mechanism = MechanismAnalysis(self.rxn_model)

    def test_graph(self):
        graph = self.mechanism.create_graph()
        graph_steps = self.mechanism.create_graph(mechanism='steps')

    def tearDown(self):
        os.remove('CO_oxidation.mkm')
        os.remove('energies.txt')

if __name__ == '__main__':
    unittest.main()

"""Test solvers of CatMAP."""
import os
import json
import pytest
import pickle
from catmap import ReactionModel
import mpmath as mp
# Decimal precision, note that CatMAP buffers the
# provided coverage with an extra 25 points
_DECIMAL_PRECISION = 100 + 25
# Decimal precision to feed into mpmath 
mp.mp.dps = _DECIMAL_PRECISION
# Compare the coverages up to this precision
_COMPARE_PRECISION = 30

@pytest.fixture
def generate_mkm_coverages(reaction='CO_oxidation', solver='coverages'):
    """Generate a set of inputs for an mkm to run tests."""

    if reaction == 'CO_oxidation':

        mkm_file = os.path.join(os.path.dirname(__file__), 'input', reaction, 'CO_oxidation.mkm')
        model = ReactionModel(setup_file=mkm_file)
        model.output_variables += ['coverage', 'production_rate']
        model.data_file = os.path.join(os.path.dirname(__file__), 'solution', reaction, f'{reaction}_{solver}.pkl')
        model.log_file = os.path.join(os.path.dirname(__file__), 'solution', reaction, f'{reaction}_{solver}.log')

    return model

@pytest.fixture
def generate_solution_coverages(reaction='CO_oxidation', benchmark_solver='coverages'):
    """Read in the solution coverages for a given reaction."""

    if reaction == 'CO_oxidation':
        solution_file = os.path.join(os.path.dirname(__file__), 'solution', reaction, f'{reaction}.json')
        with open(solution_file, 'r') as f:
            data = json.load(f)
        return data 


def test_setup_solver_coverages(generate_mkm_coverages):
    """Test solver coverages."""
    model = generate_mkm_coverages
    model.use_numbers_solver = False 
    model.fix_x_star = True

    # Check inputs are appropriate
    assert not model.use_numbers_solver 
    assert model.fix_x_star

    model.run()

    # Make sure that a coverage map is produced
    assert model.coverage_map

def test_solver_coverages(generate_mkm_coverages, generate_solution_coverages):
    """Test the solver coverages based on the coverage solver."""
    model = generate_mkm_coverages
    model.use_numbers_solver = False 
    model.fix_x_star = True
    model.run()
    data_calculation = model.coverage_map
    coverages_calculation = data_calculation[0][1]

    # Make sure that the coverages are correct
    data_benchmark = generate_solution_coverages
    # Check the coverages element by element
    for i, element_calc in enumerate(coverages_calculation):
        element_benchmark = data_benchmark['coverages'][i]
        assert element_calc == mp.mpmathify(element_benchmark)

def test_fix_x_coverages(generate_mkm_coverages, generate_solution_coverages):
    """Test the solver fix-x based on the coverages.""" 
    model = generate_mkm_coverages
    model.use_numbers_solver = True
    model.fix_x_star = True
    model.run()
    data_calculation = model.coverage_map
    coverages_calculation = data_calculation[0][1]
    # If the numbers solver is used, it will always 
    # have one more than the number of coverages as previous
    # so take only the first two in the comparison
    coverages_calculation = coverages_calculation[:2]

    # Make sure that the coverages are correct
    data_benchmark = generate_solution_coverages
    # Check the coverages element by element
    for i, element_calc in enumerate(coverages_calculation):
        element_benchmark = data_benchmark['coverages'][i]
        assert str(element_calc)[:_COMPARE_PRECISION] == str(element_benchmark)[:_COMPARE_PRECISION]

def test_free_x_coverages(generate_mkm_coverages, generate_solution_coverages):
    """Test the solver free-x based on the coverages.""" 
    model = generate_mkm_coverages
    model.use_numbers_solver = True
    model.fix_x_star = False
    model.run()
    data_calculation = model.coverage_map
    coverages_calculation = data_calculation[0][1]
    # If the numbers solver is used, it will always 
    # have one more than the number of coverages as previous
    # so take only the first two in the comparison
    coverages_calculation = coverages_calculation[:2]

    # Make sure that the coverages are correct
    data_benchmark = generate_solution_coverages
    # Check the coverages element by element
    for i, element_calc in enumerate(coverages_calculation):
        element_benchmark = data_benchmark['coverages'][i]
        assert str(element_calc)[:_COMPARE_PRECISION] == str(element_benchmark)[:_COMPARE_PRECISION]
"""Test solvers of CatMAP."""
import json
import os, shutil
from pathlib import Path

import numpy as np
import mpmath as mp
import pytest
from catmap import ReactionModel

# Decimal precision, note that CatMAP buffers the
# provided coverage with an extra 25 points
_DECIMAL_PRECISION = 100 + 25
# Decimal precision to feed into mpmath
mp.mp.dps = _DECIMAL_PRECISION
# Compare the coverages up to this precision
_COMPARE_PRECISION = 4


@pytest.fixture
def generate_mkm_coverages(tmp_path):
    def _generate_mkm_coverages(reaction):
        mkm_file = Path(__file__).parent / "config" / reaction / "input.mkm"
        energy_file = Path(__file__).parent / "config" / reaction / "energies.txt"
        shutil.copy(energy_file, tmp_path)
        os.chdir(tmp_path)
        model = ReactionModel(setup_file=str(mkm_file))
        model.output_variables += ["coverage", "production_rate"]
        data_file = tmp_path / f"coverages.pkl"
        model.data_file = str(data_file)
        log_file = tmp_path / f"coverages.log"
        model.log_file = str(log_file)
        return model

    return _generate_mkm_coverages


@pytest.fixture
def generate_solution_coverages():
    def _generate_solution_coverages(reaction):
        solution_file = Path(__file__).parent / "config" / reaction / f"data.json"
        with open(solution_file, "r") as f:
            data = json.load(f)
        return data

    return _generate_solution_coverages


@pytest.mark.parametrize("reaction", ["CO_oxidation", "CO_oxidation_ads_ads"])
@pytest.mark.parametrize("use_numbers_solver", [True, False])
@pytest.mark.parametrize("fix_x_star", [True, False])
@pytest.mark.parametrize("DEBUG", [True, False])
@pytest.mark.parametrize("numbers_type", ["squared", "exponential"])
def test_solver_setup(
    reaction: str, use_numbers_solver: bool, fix_x_star: bool, DEBUG:bool, numbers_type: str, generate_mkm_coverages
):
    model = generate_mkm_coverages(reaction)
    model.use_numbers_solver = use_numbers_solver
    model.fix_x_star = fix_x_star
    model.DEBUG = DEBUG
    model.numbers_type = numbers_type
    model.run()
    assert model.fix_x_star == fix_x_star
    assert model.use_numbers_solver == use_numbers_solver
    assert model.coverage_map


@pytest.mark.parametrize("reaction", ["CO_oxidation", "CO_oxidation_ads_ads"])
@pytest.mark.parametrize("use_numbers_solver", [True, False])
@pytest.mark.parametrize("fix_x_star", [True, False])
@pytest.mark.parametrize("DEBUG", [True, False])
def test_solver_coverages(
    reaction: str,
    use_numbers_solver: bool,
    fix_x_star: bool,
    DEBUG:bool,
    generate_mkm_coverages,
    generate_solution_coverages,
):
    model = generate_mkm_coverages(reaction)
    model.use_numbers_solver = use_numbers_solver
    model.fix_x_star = fix_x_star
    model.DEBUG = DEBUG
    model.run()
    data_calculation = model.coverage_map
    coverages_calculation = data_calculation[0][1]
    data_benchmark = generate_solution_coverages(reaction)
    coverages_calculation = np.array(coverages_calculation)
    coverage_calculation = [float(i) for i in coverages_calculation]
    coverage_benchmark = data_benchmark["coverages"]
    coverage_benchmark = np.array(coverage_benchmark, dtype=float)
    n_species = coverage_benchmark.shape[0]
    coverage_calculation = np.array(coverage_calculation)[:n_species]
    np.testing.assert_almost_equal(
        coverage_calculation,
        coverage_benchmark,
        decimal=_COMPARE_PRECISION,
    )
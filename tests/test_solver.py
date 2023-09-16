"""Test solvers of CatMAP."""
import json
import os
from pathlib import Path

import mpmath as mp
import pytest
from catmap import ReactionModel

# Decimal precision, note that CatMAP buffers the
# provided coverage with an extra 25 points
_DECIMAL_PRECISION = 100 + 25
# Decimal precision to feed into mpmath
mp.mp.dps = _DECIMAL_PRECISION
# Compare the coverages up to this precision
_COMPARE_PRECISION = 10


@pytest.fixture
def generate_mkm_coverages():
    def _generate_mkm_coverages(reaction):
        mkm_file = Path(__file__).parent / "input" / reaction / "input.mkm"
        model = ReactionModel(setup_file=str(mkm_file))
        model.output_variables += ["coverage", "production_rate"]
        data_file = Path(__file__).parent / "solution" / reaction / f"coverages.pkl"
        model.data_file = str(data_file)
        log_file = Path(__file__).parent / "solution" / reaction / f"coverages.log"
        model.log_file = str(log_file)
        return model

    return _generate_mkm_coverages


@pytest.fixture
def generate_solution_coverages():
    def _generate_solution_coverages(reaction):
        solution_file = Path(__file__).parent / "solution" / reaction / f"data.json"
        with open(solution_file, "r") as f:
            data = json.load(f)
        return data

    return _generate_solution_coverages


@pytest.mark.parametrize("reaction", ["CO_oxidation"])
@pytest.mark.parametrize("use_numbers_solver", [True, False])
@pytest.mark.parametrize("fix_x_star", [True, False])
def test_solver_setup(
    reaction: str, use_numbers_solver: bool, fix_x_star: bool, generate_mkm_coverages
):
    model = generate_mkm_coverages(reaction)
    model.use_numbers_solver = use_numbers_solver
    model.fix_x_star = fix_x_star
    model.run()
    assert model.fix_x_star == fix_x_star
    assert model.use_numbers_solver == use_numbers_solver
    assert model.coverage_map


@pytest.mark.parametrize("reaction", ["CO_oxidation"])
@pytest.mark.parametrize("use_numbers_solver", [True, False])
@pytest.mark.parametrize("fix_x_star", [True, False])
def test_solver_coverages(
    reaction: str,
    use_numbers_solver: bool,
    fix_x_star: bool,
    generate_mkm_coverages,
    generate_solution_coverages,
):
    model = generate_mkm_coverages(reaction)
    model.use_numbers_solver = use_numbers_solver
    model.fix_x_star = fix_x_star
    model.run()
    data_calculation = model.coverage_map
    coverages_calculation = data_calculation[0][1]
    data_benchmark = generate_solution_coverages(reaction)
    for idx, element_benchmark in enumerate(data_benchmark["coverages"]):
        element_calc = coverages_calculation[idx]
        assert (
            str(element_calc)[:_COMPARE_PRECISION]
            == str(element_benchmark)[:_COMPARE_PRECISION]
        )

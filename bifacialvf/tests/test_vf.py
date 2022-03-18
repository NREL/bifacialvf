"""
Tests of the view factors module.  
py.test --cov-report term-missing --cov=bifacialvf
"""
import pytest
import numpy as np
from bifacialvf.vf import getSkyConfigurationFactors
from bifacialvf.tests import (
    SKY_BETA160_C05_D1, SKY_BETA20_C05_D1, SKY_BETA20_C0_D1, SKY_BETA160_C0_D1,
    SKY_BETA160_C1_D1, SKY_BETA20_C1_D1, SKY_BETA20_C1_D0, SKY_BETA160_C1_D0,
    SKY_BETA160_C05_D0, SKY_BETA20_C05_D0, SKY_20_1_1_first, SKY_20_1_1_last,
    SKY_20_1_1_single)


@pytest.mark.parametrize('beta, C, D, expected',
    [(160, 0.5, 1, SKY_BETA160_C05_D1), (20, 0.5, 1, SKY_BETA20_C05_D1),
     (20, 0, 1, SKY_BETA20_C0_D1), (160, 0, 1, SKY_BETA160_C0_D1),
     (160, 1, 1, SKY_BETA160_C1_D1), (20, 1, 1, SKY_BETA20_C1_D1),
     (20, 1, 0, SKY_BETA20_C1_D0), (160, 1, 0, SKY_BETA160_C1_D0),
     (160, 0.5, 0, SKY_BETA160_C05_D0), (20, 0.5, 0, SKY_BETA20_C05_D0)])
def test_getSkyConfigurationFactors(beta, C, D, expected):
    """
    Benchmark against to the master branch on 2018-08-20 at 91e785d.
    """
    assert np.allclose(
        getSkyConfigurationFactors("interior", beta, C, D), expected)

@pytest.mark.parametrize('rowtype, expected',
    [('first', SKY_20_1_1_first), ('last', SKY_20_1_1_last),
     ('single', SKY_20_1_1_single)])
def test_getSkyConfigurationsRowtype(rowtype, expected):
    """
    Benchmark against to the master branch on 2018-08-20 at 91e785d.
    """
    assert np.allclose(
        getSkyConfigurationFactors(rowtype, beta=20, C=1, D=1), expected)
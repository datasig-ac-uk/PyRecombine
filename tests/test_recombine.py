import pytest
import numpy as np

from numpy.linalg import norm
from pyrecombine import recombine


@pytest.fixture
def dimension():
    return 30

@pytest.fixture
def no_points():
    return 10000

def test_recombine_1(dimension, no_points):
    data = np.random.default_rng(12345).random(size=(no_points, dimension), dtype=np.float64)

    ## test 1
    selected_points, new_weights = recombine(data)  ## degree = 1

    ## check mean preserved
    old_average = np.sum(data, 0)
    new_average = new_weights.dot(np.take(data, selected_points, 0))
    normalised_error = norm(old_average - new_average) / (norm(old_average) + norm(new_average))

    ## report
    assert len(selected_points) <= dimension + 1
    assert normalised_error <= 1e-13


def test_recombine_2(dimension, no_points):
    rng = np.random.default_rng(12345)

    data = rng.random(size=(no_points, dimension))
    ## test2
    ### the points are not spanning the full space and so the minimal set should have cardinality less than or equal rank + 1
    matrix = rng.random(size=(dimension, dimension + 20))
    new_data = data.dot(matrix)
    selected_points, new_weights = recombine(new_data)  ## degree = 1

    ## check mean preserved
    old_average = np.sum(data, 0)
    new_average = new_weights.dot(np.take(data, selected_points, 0))
    normalised_error = norm(old_average - new_average) / (norm(old_average) + norm(new_average))

    ## report
    assert len(selected_points) <= dimension + 1
    assert normalised_error <= 1e-12

def test_recombine_3():
    rng = np.random.default_rng(12345)

    # test3
    ## test the degree > 1 case - match second moments
    dimension = 10
    no_points = 1000
    data = rng.random(size=(no_points, dimension))

    selected_points, new_weights = recombine(data, degree=2)

    old_average = np.sum(data, 0)
    new_average = new_weights.dot(np.take(data, selected_points, 0))
    normalised_error_in_mean = norm(old_average - new_average) / (norm(old_average) + norm(new_average))

    new_cov = np.cov(np.take(data, selected_points, 0), rowvar=False, bias=True, aweights=new_weights)
    old_cov = np.cov(data, rowvar=False, bias=True, aweights=np.full(1000, 1.))
    normalised_error_in_cov = norm(old_cov - new_cov) / (norm(old_cov) + norm(new_cov))

    assert normalised_error_in_mean <= 1e-13
    assert normalised_error_in_cov <= 1e-13
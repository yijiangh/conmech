import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal_nulp, assert_array_almost_equal
from pyconmech.frame_analysis import create_local_stiffness_matrix
# from pyconmech.frame_analysis import numpy_stiffness

@pytest.fixture
def tol():
    return 1e-8

def assert_aleq_gmz_array(a, b, precision=2, np_test=False):
    # 2digit precision is used in [GMZ] results
    if np_test:
        assert_array_almost_equal(a, b)
    else:
        assert a.shape == b.shape
        a_iter = np.nditer(a, flags=['multi_index'])
        b_iter = np.nditer(b, flags=['multi_index'])
        while not a_iter.finished:
            assert np.format_float_scientific(a_iter[0], precision=precision) == \
                np.format_float_scientific(b_iter[0], precision=precision), \
                    '{} : a {} | b {}'.format(a_iter.multi_index, a_iter[0], b_iter[0])
            a_iter.iternext()
            b_iter.iternext()

@pytest.mark.local_stiff
def test_2Dbeam_stiffness_matrix(tol):
    # base on example 4.8, p.79 [MGZ2000]
    # assume no bending normal to the plane of the paper
    # notice that in [MGZ2000], unit convention:
    # length : mm (angle: rad)
    # force : kN
    # pressure : kN / mm2 = 1e9 Pa = 1e3 MPa
    # moment : kN m
    # material properties
    E = 200000. # MPa
    E *= 1e-3 # => kN / mm2
    mu = 0.3 # Poisson ratio

    # area and inertia
    A = 6*1e3 #mm2
    Iz = 200*1e6 # mm4
    Jx = 300*1e3 #mm4

    L = 8 #m
    L *= 1e3 # => mm
    Iy = Iz # not used

    K_full = create_local_stiffness_matrix(L, A, Jx, Iy, Iz, E, mu)
    # omit out-of-plane shear and bending dof 
    # delete (My1, My2, w1, \theta_y1, w2, \theta_y1)
    ids = list(range(6))
    ids.remove(2) # Fz1
    ids.remove(4) # My1
    ids.extend([i+6 for i in ids])
    K = K_full[np.ix_(ids,ids)]

    K_ans = np.zeros([8,8])
    K_ans[0,:] = [0.75,0,0,0,-0.75,0,0,0]
    K_ans[1,1:] = [0.00469, 0, 18.750, 0, -0.00469, 0, 18.750]
    K_ans[2,2:] = [14.423, 0, 0, 0, -14.423, 0]
    K_ans[3,3:] = [1.0*1e5, 0, -18.750, 0, 0.5*1e5]
    K_ans[4,4:] = [0.75, 0, 0, 0]
    K_ans[5,5:] = [0.00469, 0, -18.750]
    K_ans[6,6:] = [14.423, 0]
    K_ans[7,7]  = 1.0*1e5
    K_ans = K_ans + K_ans.T
    K_ans[range(8), range(8)] = np.diag(K_ans) / 2
    K_ans *= 200

    assert_aleq_gmz_array(K_ans, K)


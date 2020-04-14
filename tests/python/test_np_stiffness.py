import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal_nulp, assert_array_almost_equal
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.sparse import find
import scipy.sparse.linalg as SPLIN
from pyconmech.frame_analysis import create_local_stiffness_matrix, global2local_transf_matrix, \
    assemble_global_stiffness_matrix, get_element_shape_fn

# from pyconmech.frame_analysis import numpy_stiffness

def isPSD(A, tol = 1e-8):
    vals, _ = SPLIN.eigsh(A, k = 2, which = 'BE') 
    # return the ends of spectrum of A
    return np.all(vals > -tol)

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

def assert_aleq_sparse(A, B, atol = 1e-14):
    assert np.array_equal(A.shape, B.shape)
    r1,c1,v1 = find(A)
    r2,c2,v2 = find(B)
    index_match = np.array_equal(r1,r2) & np.array_equal(c1,c2)
    assert index_match!=0, 'number of nonzero id not matched: {}'.format(index_match)
    assert np.allclose(v1,v2, atol=atol)

@pytest.mark.np_wip
def test_2Dbeam_stiffness_matrix(viewer):
    # base on example 4.8-4.12, p.79 [MGZ2000]
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
    # Iy = Iz # not used
    Iy = 0 # not used

    K_full = create_local_stiffness_matrix(L, A, Jx, Iy, Iz, E, mu)
    assert (12,12) == K_full.shape
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

    # assembl & solve stiffness matrix
    A_bc = 4000 # mm2
    Iz_bc = 50*1e6 # mm4
    J_bc = 100*1e3 # mm4
    L_bc = 5 # m
    L_bc *= 1e3 # => mm

    K_bc = create_local_stiffness_matrix(L_bc, A_bc, J_bc, 0, Iz_bc, E, mu)
    assert (12,12) == K_full.shape
    assert (12,12) == K_bc.shape
    k_list = [K_full, K_bc]

    nV = 3
    id_map = np.vstack([np.array(range(0, 12)), np.array(range(6, 18))])
    Ksp = assemble_global_stiffness_matrix(k_list, nV, id_map)
    assert_aleq_sparse(Ksp, Ksp.T)
    assert isPSD(Ksp)

    # fixities
    # boundary condition: 
    # end a fully fixed, v_b = 0, theta_b = theta_c = 0
    v_map = np.vstack([np.array(range(v*6, (v+1)*6)) for v in range(3)])
    enabled_dofs = []
    for v in range(3):
        enabled_dofs.extend([0+v*6,1+v*6,3+v*6,5+v*6])

    dof_ids = [v_map[1,0], v_map[2,0], v_map[2,1], v_map[1,5], v_map[2,5]]
    
    Pc = 5/np.sqrt(2) # kN
    load_ff = np.array([0, Pc, -Pc, 0, 0])

    Ksp_ff = Ksp[np.ix_(dof_ids, dof_ids)]
    # print(Ksp_ff / 200)
    u = SPLIN.spsolve(Ksp_ff, load_ff)

    # u_b, theta_zb, u_c, v_c, theta_zc
    # u_ans = np.array([0.024, -0.00088, 0.046, -19.15,  -0.00530])
    # u_b, u_c, v_c, theta_zb, theta_zc
    u_ans = np.array([0.024, 0.046, -19.15, -0.00088, -0.00530])
    assert_aleq_gmz_array(u_ans, u, precision=1)

    # compute reactions
    # R_xa, R_ya, R_mza, R_yb
    fixed_dofs = [0,1,5,7]
    Ksp_fs = Ksp[np.ix_(fixed_dofs, dof_ids)]
    R = Ksp_fs.dot(u)
    assert R.shape[0] == 4

    # R_xa, R_ya, Rm_za, R_yb
    # textbook moment reaction unit kN m
    R_ans = np.array([-3.6, -3.3, -8.8*1e3, 6.85])
    # (God knows what type of hand-calc precision the authors were using...)
    assert_aleq_gmz_array(R_ans, R, precision=0)

    # * plot moment diagram

    # * plot deflected structure
    end_pts = np.array([[0.,0.,0], [0., 8, 0], [0, 13., 0]])
    end_pts *= 1e3 # m => mm
    d_a = np.array([0,0,0,0,0,0])
    d_b = np.array([0,u[1],0,0,0,u[3]])
    shape_fn = get_element_shape_fn(end_pts[0, :], end_pts[1, :], d_a, d_b)

    if viewer:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        x_plt = np.linspace(0., 0.008, 50) 
        u_plt = np.array([shape_fn(x) for x in x_plt])
        # print('s: ', u_plt.shape)
        # https://matplotlib.org/3.1.1/gallery/subplots_axes_and_figures/axis_equal_demo.html

        ax.plot(u_plt[:,0], u_plt[:,1], u_plt[:,2])
        ax.legend()
        plt.show()

    assert_array_almost_equal(np.zeros((3)), shape_fn(0.0))
    assert_array_almost_equal(end_pts[1,:] + np.array([0, 0.024, 0]), shape_fn(8*1e3))


def test_transf_matrix():
    pt1 = np.array([0,0,0])
    pt2 = np.array([1,0,0])
    R = global2local_transf_matrix(pt1, pt2)
    assert_array_almost_equal(np.eye(3), R)
import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal_nulp, assert_array_almost_equal
from matplotlib import pyplot as plt
from scipy.sparse import find
import scipy.sparse.linalg as SPLIN
from pyconmech.frame_analysis import create_local_stiffness_matrix, global2local_transf_matrix, \
    assemble_global_stiffness_matrix, get_element_shape_fn, get_internal_reaction_fn
from pyconmech.frame_analysis import bending_stiffness_matrix

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
        # assert type(a) == type(b)
        if isinstance(a, float) and isinstance(b, float):
            a = np.array([a])
            b = np.array([b])
        assert a.shape == b.shape
        a_iter = np.nditer(a, flags=['multi_index'])
        b_iter = np.nditer(b, flags=['multi_index'])
        while not a_iter.finished:
            if abs(a_iter[0]-b_iter[0]) > 1e-10:
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

# @pytest.mark.skip('not fully developed')
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

    K_ab_full = create_local_stiffness_matrix(L, A, Jx, Iy, Iz, E, mu)
    assert (12,12) == K_ab_full.shape
    # omit out-of-plane shear and bending dof 
    # delete (My1, My2, w1, \theta_y1, w2, \theta_y1)
    ids = list(range(6))
    ids.remove(2) # Fz1
    ids.remove(4) # My1
    ids.extend([i+6 for i in ids])
    K = K_ab_full[np.ix_(ids,ids)]

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
    assert (12,12) == K_ab_full.shape
    assert (12,12) == K_bc.shape
    k_list = [K_ab_full, K_bc]

    nV = 3
    # node id : dof id map
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
    print(Ksp_ff.todense())
    # print(Ksp_ff / 200)
    # TODO: should always solve for -P here?
    u = SPLIN.spsolve(Ksp_ff, load_ff)

    # u_b, theta_zb, u_c, v_c, theta_zc
    # u_ans = np.array([0.024, -0.00088, 0.046, -19.15,  -0.00530])
    # u_b, u_c, v_c, theta_zb, theta_zc
    u_ans = np.array([0.024, 0.046, -19.15, -0.00088, -0.00530])
    assert_aleq_gmz_array(u_ans, u, precision=1)

    d_a = np.zeros(6)
    d_b = np.array([u[0],0,0,    0,0,u[3]])
    d_c = np.array([u[1],u[2],0, 0,0,u[4]])

    # * compute reactions
    # - fixities reaction
    # R_xa, R_ya, R_mza, R_yb
    fixed_dofs = [0,1,5,7]
    Ksp_fs = Ksp[np.ix_(fixed_dofs, dof_ids)]
    R = Ksp_fs.dot(u)

    # R_xa, R_ya, Rm_za, R_yb
    # textbook moment reaction unit kN m
    # convert kN m to kN mm here
    R_ans = np.array([-3.6, -3.3, -8.8*1e3, 6.85])
    # (God knows what type of hand-calc precision the authors were using...)
    assert_aleq_gmz_array(R_ans, R, precision=0)

    # fixities at node b
    # node_a_fr = np.array([R[0],R[1],0,0,0,R[2]])
    # node_b_fr = np.array([0,R[3],0,0,0,0])

    # - internal reaction
    # no need to do transformation here
    # beam_ab_f = K_ab_full.dot(np.hstack([d_a, d_b]))
    # beam_bc_f = K_bc.dot(np.hstack([d_b, d_c]))
    # node_a_f = beam_ab_f[:6]
    # assert_array_almost_equal(node_a_fr, node_a_f)

    # node_b_f = beam_ab_f[6:12]
    # node_c_f = beam_bc_f[6:12]
    # # assert_array_almost_equal(node_b_f + node_b_fr + beam_bc_f[:6], np.zeros(6))

    # # * plot deflected structure
    # end_pts = np.array([[0.,0.,0], [8., 0., 0], [13., 0., 0]])
    # end_pts *= 1e3 # m => mm

    # beam_a_shape_fn = get_element_shape_fn(end_pts[0, :], end_pts[1, :], d_a, d_b)
    # assert_array_almost_equal(np.zeros((3)), beam_a_shape_fn(0.0))
    # assert_array_almost_equal(end_pts[1,:] + d_b[:3], beam_a_shape_fn(1.))

    # beam_b_shape_fn = get_element_shape_fn(end_pts[1, :], end_pts[2, :], d_b, d_c)
    # assert_array_almost_equal(end_pts[1,:] + d_b[:3], beam_b_shape_fn(0.0))
    # assert_array_almost_equal(end_pts[2,:] + d_c[:3], beam_b_shape_fn(1.))

    # beam_a_reaction_fn = get_internal_reaction_fn(node_a_f, node_b_f)
    # beam_b_reaction_fn = get_internal_reaction_fn(node_b_f, node_c_f)
    # assert_aleq_gmz_array(8.8*1e3, beam_a_reaction_fn(0.)[-1], precision=1)
    # assert_aleq_gmz_array(-17.68*1e3, beam_a_reaction_fn(1.)[-1], precision=1)
    # assert_aleq_gmz_array(-17.68*1e3, beam_b_reaction_fn(0.)[-1], precision=1)
    # assert_aleq_gmz_array(0.0, beam_b_reaction_fn(1.)[-1])

    # if viewer:
    #     # fig = plt.figure()
    #     # ax = fig.gca(projection='2d')
    #     # https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.subplots.html#examples-using-matplotlib-pyplot-subplots
    #     fig, axes = plt.subplots(1, 1+6)

    #     ax = axes[0]
    #     # original nodal position
    #     ax.plot(end_pts[:,0], end_pts[:,1], marker='o')
    #     # ax.scatter(end_pts[:,0], end_pts[:,1], marker='o')

    #     t_plt = np.linspace(0., 1., 50) 
    #     beam_a_u_plt = np.array([beam_a_shape_fn(t) for t in t_plt])
    #     ax.plot(beam_a_u_plt[:,0], beam_a_u_plt[:,1], beam_a_u_plt[:,2])

    #     beam_b_u_plt = np.array([beam_b_shape_fn(t) for t in t_plt])
    #     ax.plot(beam_b_u_plt[:,0], beam_b_u_plt[:,1], beam_b_u_plt[:,2])

    #     ax.set_xlabel('x / mm')
    #     ax.set_ylabel('y / mm')
    #     # ax.set_zlabel('z / mm')
    #     ax.legend()
    #     ax.set_title('Deflected shape')
    #     # https://matplotlib.org/3.1.1/gallery/subplots_axes_and_figures/axis_equal_demo.html
    #     # ax.set_aspect('equal', 'box')

    #     titles = ['Nx / kN', 'Fy / kN', 'Fz / kN', 'Mx / kN-m', 'My / kN-m', 'Mz / kN-m']
    #     for i in range(len(axes)-1):
    #         ax_R = axes[i+1]
    #         ax_R.plot(t_plt, [beam_a_reaction_fn(t)[i]*1e-3 for t in t_plt])
    #         ax_R.plot(t_plt+np.ones(t_plt.shape), [beam_b_reaction_fn(t)[i]*1e-3 for t in t_plt])

    #         dt_plt = np.linspace(0., 2., 100) 
    #         ax_R.plot(dt_plt, np.zeros(dt_plt.shape), linewidth=0.5)
    #         # ax_R.scatter([0,1,2], [-8.8, 17.68, 0])
    #         ax_R.set_xlabel('t')
    #         # ax_R.set_ylabel(labels[i])
    #         ax_R.set_title(titles[i])
    #         ax_R.legend()

    #     plt.show()


def test_transf_matrix():
    pt1 = np.array([0,0,0])
    pt2 = np.array([1,0,0])
    R = global2local_transf_matrix(pt1, pt2)
    assert_array_almost_equal(np.eye(3), R)

@pytest.mark.rot_stiff
def test_rotational_stiffness():
    Kz = bending_stiffness_matrix(1.0, 1.0, 1.0, cr1=0, cr2=0)
    print(Kz)

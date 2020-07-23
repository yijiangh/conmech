import pytest
from numpy.testing import assert_array_almost_equal_nulp, assert_array_almost_equal

import numpy as np
from numpy.linalg import norm
import scipy
import scipy.sparse.linalg as SPLA
from scipy.sparse import csc_matrix
from matplotlib import pyplot as plt
from scipy.sparse import find
from pyconmech.frame_analysis import create_local_stiffness_matrix, global2local_transf_matrix, \
    assemble_global_stiffness_matrix
from pyconmech.frame_analysis.visualization import get_element_shape_fn, get_internal_reaction_fn
from pyconmech.frame_analysis.numpy_stiffness import bending_stiffness_matrix, axial_stiffness_matrix, mu2G, G2mu, \
    turn_diagblock
from termcolor import cprint

# from pyconmech.frame_analysis import numpy_stiffness

def isPSD(A, tol = 1e-8):
    vals, _ = SPLA.eigsh(A, k = 2, which = 'BE') 
    # return the ends of spectrum of A
    return np.all(vals > -tol)

def isPD(A, tol = 1e-8):
    vals, _ = SPLA.eigsh(A, k = 2, which = 'BE') 
    return np.all(vals > tol)

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

@pytest.mark.flat_beam
def test_2Dbeam_stiffness_matrix(viewer):
    # base on example 4.8-4.12, p.79 [MGZ2000]
    # assume no bending normal to the plane of the paper
    # notice that in [MGZ2000], unit convention:
    # length : mm (angle: rad)
    # force : kN
    # pressure : kN / mm2 = 1e9 Pa = 1e3 MPa
    # moment : kN m
    # material properties

    # TODO: use global stiffness matrix and dof map

    E = 200000. # MPa
    E *= 1e-3 # => kN / mm2
    mu = 0.3 # Poisson ratio

    # area and inertia
    cross_sec = {}
    cross_sec[0] = {
        'A' : 6*1e3, #mm2
        'Iz' : 200*1e6, # mm4
        'Jx' : 300*1e3, #mm4
    }
    cross_sec[1] = {
        'A' : 4000, # mm2
        'Iz' : 50*1e6, # mm4
        'Jx' : 100*1e3, # mm4
    }

    nV = 3
    nE = 2

    V = np.array([[0.,0,0], [8,0,0], [13.,0,0]])*1e3 # m => mm
    T = np.array([[0,1], [1,2]], dtype=int)

    # element id : dof id map
    id_map = np.vstack([np.array(range(i*6, (i+2)*6)) for i in range(nE)])

    def compute_R(e):
        R_3 = global2local_transf_matrix(V[e[0],:], V[e[1],:]) #, rot_y2x=np.pi)
        return turn_diagblock(R_3)

    def compute_Ke(e, A, Jx, Iz, return_global=True):
        Iy = 0 # not used
        Le = norm(V[e[0],:] - V[e[1],:])
        Ke_l = create_local_stiffness_matrix(Le, A, Jx, Iy, Iz, E, mu)
        if not return_global:
            return Ke_l
        else:
            R_f = compute_R(e)
            Ke = ((R_f.T).dot(Ke_l)).dot(R_f)
            return Ke

    k_list = [compute_Ke(T[i,:], cross_sec[i]['A'], cross_sec[i]['Jx'], cross_sec[i]['Iz']) for i in range(nE)]
    K_full = assemble_global_stiffness_matrix(k_list, nV, id_map)

    # compute permutation map
    total_dof = nV*6
    dof_stat = np.ones(total_dof, dtype=int)*-1
    dof_stat[id_map[0,:6]] = [1,1,1,1,1,1]
    dof_stat[id_map[0,6:]] = [0,1,1,1,1,0]
    dof_stat[id_map[1,6:]] = [0,0,1,1,1,0]
    cprint('dof specs:', 'green')
    print(dof_stat[id_map[0,:6]])
    print(dof_stat[id_map[0,6:]])
    print(dof_stat[id_map[1,6:]])

    n_fixed_dof = np.sum(dof_stat)
    n_free_dof = total_dof - n_fixed_dof

    # K_ab = create_local_stiffness_matrix(L, A, Jx, Iy, Iz, E, mu)
    K_ab = k_list[0]
    assert (12,12) == K_ab.shape
    # omit out-of-plane shear and bending dof 
    # delete (My1, My2, w1, \theta_y1, w2, \theta_y1)
    ids = list(range(6))
    ids.remove(2) # Fz1
    ids.remove(4) # My1
    ids.extend([i+6 for i in ids])
    K = K_ab[np.ix_(ids,ids)]

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

    # node id : dof id map
    assert_aleq_sparse(K_full, K_full.T)
    assert isPSD(K_full)

    # load vector
    P = np.zeros(total_dof)
    Pc = 5/np.sqrt(2) # kN
    P[id_map[1,6:]] = np.array([Pc, -Pc, 0, 0, 0, 0])

    free_tail = 0
    fix_tail = n_free_dof
    id_map_RO = np.array(range(total_dof))
    for i in range(total_dof):
        if 0 == dof_stat[i]:
            id_map_RO[free_tail] = i
            free_tail += 1
        elif 1 == dof_stat[i]:
            id_map_RO[fix_tail] = i
            fix_tail += 1
        else:
            raise ValueError

    # a row permuatation matrix (multiply left)
    perm_data = []
    perm_row = []
    perm_col= []
    for i in range(total_dof):
        perm_row.append(i)
        perm_col.append(id_map_RO[i])
        perm_data.append(1)
    Perm = csc_matrix((perm_data, (perm_row, perm_col)), shape=(total_dof, total_dof))
    PermT = Perm.T

    # permute the full stiffness matrix & carve the needed portion out
    K_perm = (Perm.dot(K_full)).dot(PermT)
    K_mm = K_perm[0:n_free_dof, 0:n_free_dof]
    K_fm = K_perm[n_free_dof:n_free_dof+n_fixed_dof, 0:n_free_dof]

    # eigenvalues, eigenvectors = SPLA.eigsh(K_mm, k=5)
    eigenvalues, eigenvectors = scipy.linalg.eigh(K_mm.toarray())
    print('Eigen values: ', eigenvalues)
    for val, vec in zip(eigenvalues, eigenvectors):
        if abs(val) < 1e-8:
            vv = np.zeros(total_dof)
            vv[0:n_free_dof] = vec
            vv_p = PermT.dot(vv)
            print('Eigen val: {} | Eigen vec: {}'.format(val, vv_p))
    assert isPD(K_mm)

    P_perm = Perm.dot(P) 
    P_m = P_perm[0:n_free_dof]
    P_f = P_perm[n_free_dof:n_free_dof+n_fixed_dof]

    # D_diag = np.array([1/np.sqrt(a) for a in K_mm.diagonal()])
    # D = scipy.sparse.diags(D_diag, format='csc')
    # # D_inv = scipy.sparse.diags([np.sqrt(a) for a in K_mm.diagonal()], format='csc')
    # K_mm_precond = D.dot(K_mm)

    # K_LU = SPLA.splu(K_mm_precond, diag_pivot_thresh=0, options=dict(SymmetricMode=True)) 
    # U_m_pred = K_LU.solve(D.dot(P_m))
    # U_m = U_m_pred
    # K_LU = SPLA.splu(K_mm, diag_pivot_thresh=0, options=dict(SymmetricMode=True)) 
    # U_m = K_LU.solve(P_m)

    U_m = SPLA.spsolve(K_mm, P_m)

    print('P_m: ', P_m)
    print('U_m: ', U_m)
    assert_array_almost_equal(K_mm.dot(U_m), P_m)

    # displacement
    Utmp = np.zeros(total_dof)
    Utmp[0:n_free_dof] = U_m
    U = PermT.dot(Utmp)
    cprint('Displacement:', 'green')
    for i in range(nV):
        print('N#{}: {} | {}'.format(i, U[6*i:6*i+3], U[6*i+3:6*i+6]))

    # u_b, u_c, v_c, theta_zb, theta_zc
    u_ans = np.array([0.024, 0.046, -19.15, -0.00088, -0.00530])
    assert_aleq_gmz_array(u_ans, U[[6, 6*2, 6*2+1, 6+5, 6*2+5]], precision=1)

    # reaction
    # print('Pf: ', P_f)
    Rtmp = np.zeros(total_dof)
    Rtmp[n_free_dof:n_free_dof+n_fixed_dof] = K_fm.dot(U_m) - P_f
    R = PermT.dot(Rtmp)
    cprint('Reaction:', 'green')
    for i in range(nV):
        print('R#{}: {} | {}'.format(i, R[6*i:6*i+3], R[6*i+3:6*i+6]))

    # R_xa, R_ya, Rm_za, R_yb
    # textbook moment reaction unit kN m
    # convert kN m to kN mm here
    R_ans = np.array([-3.6, -3.3, -8.8*1e3, 6.85])
    # (God knows what type of hand-calc precision the authors were using...)
    assert_aleq_gmz_array(R_ans, R[[0, 1, 5, 6+1]], precision=0)

    Ues = []
    Res = []
    for i in range(nE):
        e = T[i,:]
        # Ke = compute_Ke(e, *e_cr_y[i])
        # Ke = compute_Ke(e, return_global=False)
        Ke = k_list[i]
        R_f = compute_R(e)
        print('local axis E#{}:\n{}'.format(i, R_f[:3, :3]))
        print('dof id: ', id_map[i, :])
        Ue = U[id_map[i, :]]
        # Ues.append(R_f.dot(Ue))
        Ues.append(Ue)
        Res.append(R_f.dot(Ke.dot(Ue)))
        # Res.append(Ke.dot(Ue))

    cprint('element displacement in global axis', 'green')
    for i in range(nE):
        Ue = Ues[i]
        print('-'*5)
        print('E#{}:'.format(i))
        print('N1: {} | {}\n'.format(Ue[:3], Ue[3:6]))
        print('N2: {} | {}\n'.format(Ue[6:9], Ue[9:12]))

    cprint('local reaction', 'green')
    for i in range(nE):
        Re = Res[i]
        print('-'*5)
        print('E#{}:'.format(i))
        print('Reaction E#{}:'.format(i))
        print('N1: {} | {}\n'.format(-Re[:3], -Re[3:6]))
        print('N2: {} | {}\n'.format(Re[6:9], Re[9:12]))

    # fixities at node b
    node_a_fr = R[:6]
    assert_array_almost_equal(node_a_fr, Res[0][:6])

    node_b_fr = R[6:12]
    assert_array_almost_equal(node_b_fr, Res[0][6:12] + Res[1][:6])

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
    print('*'*10)
    L = 10.0
    E = 300.0
    Iz = 2.0

    Kz = bending_stiffness_matrix(L, E, Iz, cr1=0, cr2=0)
    assert_array_almost_equal(Kz, np.zeros((4,4)))
    print('Both released passed.')

    K = np.zeros((4,4))
    for axis in [1,2]:
        sign = 1. if axis==2 else -1.
        K[0,:] = np.array([12.,      sign*6*L, -12.,      sign*6*L])
        K[1,:] = np.array([sign*6*L, 4*(L**2), sign*-6*L, 2*(L**2)])
        K[2,:] = -K[0,:]
        K[3,:] = np.array([sign*6*L, 2*(L**2), -sign*6*L, 4*(L**2)])
        K *= (E*Iz/L**3)
        Kz = bending_stiffness_matrix(L, E, Iz, axis=axis, cr1=np.inf, cr2=np.inf)
        assert_array_almost_equal(Kz, K)
    print('Both fixed passed.')

    K = np.array([[0,3/L,-3/L**2,0], [0,0,-3/L,0], [0,0,0,0], [0,0,0,0]])
    K += K.T
    K += np.diag([3/L**2,3,3/L**2, 0])
    K *= E*Iz/L
    Kz = bending_stiffness_matrix(L, E, Iz, cr1=np.inf, cr2=0)
    assert_array_almost_equal(Kz, K)

    K = np.array([[0,0,-3/L**2,3/L], [0,0,0,0], [0,0,0,-3/L], [0,0,0,0]])
    K += K.T
    K += np.diag([3/L**2,0,3/L**2, 3])
    K *= E*Iz/L
    Kz = bending_stiffness_matrix(L, E, Iz, cr1=0, cr2=np.inf)
    assert_array_almost_equal(Kz, K)
    print('One fixed one released passed.')

    Kz = bending_stiffness_matrix(L, E, Iz, cr1=0, cr2=0.4)
    assert_array_almost_equal(Kz[1,:], np.zeros(4))
    assert_array_almost_equal(Kz[:,1], np.zeros(4))
    print(Kz)
    Kz = bending_stiffness_matrix(L, E, Iz, cr1=0.4, cr2=0.)
    assert_array_almost_equal(Kz[3,:], np.zeros(4))
    assert_array_almost_equal(Kz[0,3], np.zeros(4))
    print(Kz)
    print('One rot stiff one released passed.')

    Kz = bending_stiffness_matrix(L, E, Iz, cr1=np.inf, cr2=0.4)
    print(Kz)
    Kz = bending_stiffness_matrix(L, E, Iz, cr1=0.4, cr2=np.inf)
    print(Kz)
    print('One rot stiff one fixed passed.')

    Kz = bending_stiffness_matrix(L, E, Iz, cr1=0.4, cr2=0.4)
    print(Kz)
    print('Both rot stiff passed.')

@pytest.mark.axial_stiff
def test_axial_stiffness():
    print('*'*10)
    L = 1.0
    A = 0.02
    E = 2e4

    Ka = axial_stiffness_matrix(L, A, E, cr1=0, cr2=0)
    assert_array_almost_equal(Ka, np.zeros((2,2)))
    print('Both released passed.')

    K = np.ones([2,2])
    K[0,1] = -1.
    K[1,0] = -1.
    K *= (E*A/L)
    Ka = axial_stiffness_matrix(L, A, E, cr1=np.inf, cr2=np.inf)
    assert_array_almost_equal(Ka, K)
    print('Both fixed passed.')

    Ka = axial_stiffness_matrix(L, A, E, cr1=0, cr2=np.inf)
    assert_array_almost_equal(Ka, np.zeros((2,2)))
    Ka = axial_stiffness_matrix(L, A, E, cr1=np.inf, cr2=0)
    assert_array_almost_equal(Ka, np.zeros((2,2)))
    Ka = axial_stiffness_matrix(L, A, E, cr1=0, cr2=0.4)
    assert_array_almost_equal(Ka, np.zeros((2,2)))
    Ka = axial_stiffness_matrix(L, A, E, cr1=0.4, cr2=0)
    assert_array_almost_equal(Ka, np.zeros((2,2)))
    print('One joint released passed.')

    Ka = axial_stiffness_matrix(L, A, E, cr1=0.4, cr2=0.5)
    K = np.ones([2,2])
    K[0,1] = -1.
    K[1,0] = -1.
    K *= 1/(L/(E*A)+ 1/0.4 + 1/0.5)
    assert_array_almost_equal(Ka, K)
    print('Both joints rotational stiffness passed.')

@pytest.mark.pointup_truss
def test_joint_release_solve():
    A =  0.002 # m2
    Jx = 3.6630108993439872e-08 # m4
    Iy = 3.3194666666666672e-06
    Iz = 1.3341666666666669e-06

    E = 210000000.0 # kN/m2
    G12 = 80760000.0 # kN/m2 
    mu = G2mu(G12, E)
    assert mu < 0.5 and mu > 0

    V = np.array([[-10.,0,0], [0,0,5.], [10., 0., 0.]])
    T = np.array([[0,1], [1,2]], dtype=int)
    nV = 3
    nE = 2

    # V = np.array([[-10.,0,0], [0,0,5.], [10.,0,0]])
    # T = np.array([[0,1], [1,2]], dtype=int)
    # nV = 2
    # nE = 1

    # ebending = [True, True]
    ebending = [False, False]
    e_cr_y = [[np.inf, np.inf], [np.inf, np.inf]]
    # e_cr_y = [[np.inf, 0.004], [0.008, np.inf]]
    # e_cr_y = [[np.inf, 0.004], [np.inf, 0.008]]
    # e_cr_y = [[np.inf, 0.0], [0.0, np.inf]]
    # e_cr_y = [[np.inf, 0.75], [0.5, np.inf]]
    for i in range(2):
        if not ebending[i]:
            e_cr_y[i] = [0, 0]

    # element id : dof id map
    id_map = np.vstack([np.array(range(i*6, (i+2)*6)) for i in range(nE)])
    # id_map = np.vstack([np.array(range(0, 12))])

    def compute_R(e):
        R_3 = global2local_transf_matrix(V[e[0],:], V[e[1],:]) #, rot_y2x=np.pi)
        return turn_diagblock(R_3)

    def compute_Ke(e, cr_y1=np.inf, cr_y2=np.inf, return_global=True):
        Le = norm(V[e[0],:] - V[e[1],:])
        Ke_l = create_local_stiffness_matrix(Le, A, Jx, Iy, Iz, E, mu, cr1=[np.inf, cr_y1, np.inf], cr2=[np.inf, cr_y2, np.inf])
        if not return_global:
            return Ke_l
        else:
            R_f = compute_R(e)
            Ke = ((R_f.T).dot(Ke_l)).dot(R_f)
            return Ke

    k_list = [compute_Ke(T[i,:], *e_cr_y[i]) for i in range(nE)]
    K_full = assemble_global_stiffness_matrix(k_list, nV, id_map)

    # compute permutation map
    total_dof = nV*6
    dof_stat = np.ones(total_dof, dtype=int)*-1
    dof_stat[id_map[0,:6]] = [1,1,1,1,int(not ebending[0]),1]
    dof_stat[id_map[0,6:]] = [0,1,0,1,int(not (ebending[0] or ebending[1])),1]
    dof_stat[id_map[1,6:]] = [1,1,1,1,int(not ebending[1]),1]
    # dof_stat[id_map[0,:6]] = [1,1,1,1,0,1]
    # dof_stat[id_map[0,6:]] = [0,1,0,1,0,1]
    # dof_stat[id_map[1,6:]] = [1,1,1,1,0,1]
    cprint('dof specs:', 'green')
    print(dof_stat[id_map[0,:6]])
    print(dof_stat[id_map[0,6:]])
    print(dof_stat[id_map[1,6:]])

    n_fixed_dof = np.sum(dof_stat)
    n_free_dof = total_dof - n_fixed_dof

    # load vector
    P = np.zeros(total_dof)
    P[2] = -10.0 #kN
    P[6+2] = -10.0 #kN
    P[12+2] = -10.0 #kN

    free_tail = 0
    fix_tail = n_free_dof
    id_map_RO = np.array(range(total_dof))
    for i in range(total_dof):
        if 0 == dof_stat[i]:
            id_map_RO[free_tail] = i
            free_tail += 1
        elif 1 == dof_stat[i]:
            id_map_RO[fix_tail] = i
            fix_tail += 1
        else:
            raise ValueError

    # a row permuatation matrix (multiply left)
    perm_data = []
    perm_row = []
    perm_col= []
    for i in range(total_dof):
        perm_row.append(i)
        perm_col.append(id_map_RO[i])
        perm_data.append(1)
    Perm = csc_matrix((perm_data, (perm_row, perm_col)), shape=(total_dof, total_dof))
    PermT = Perm.T

    # permute the full stiffness matrix & carve the needed portion out
    K_perm = (Perm.dot(K_full)).dot(PermT)
    K_mm = K_perm[0:n_free_dof, 0:n_free_dof]
    K_fm = K_perm[n_free_dof:n_free_dof+n_fixed_dof, 0:n_free_dof]

    # eigenvalues, eigenvectors = SPLA.eigsh(K_mm, k=5)
    eigenvalues, eigenvectors = scipy.linalg.eigh(K_mm.toarray())
    print('Eigen values: ', eigenvalues)
    for val, vec in zip(eigenvalues, eigenvectors):
        if abs(val) < 1e-8:
            vv = np.zeros(total_dof)
            vv[0:n_free_dof] = vec
            vv_p = PermT.dot(vv)
            print('Eigen val: {} | Eigen vec: {}'.format(val, vv_p))
    # assert isPD(K_mm)

    P_perm = Perm.dot(P) 
    P_m = P_perm[0:n_free_dof]
    P_f = P_perm[n_free_dof:n_free_dof+n_fixed_dof]

    # D_diag = np.array([1/np.sqrt(a) for a in K_mm.diagonal()])
    # D = scipy.sparse.diags(D_diag, format='csc')
    # # D_inv = scipy.sparse.diags([np.sqrt(a) for a in K_mm.diagonal()], format='csc')
    # K_mm_precond = D.dot(K_mm)

    # K_LU = SPLA.splu(K_mm_precond, diag_pivot_thresh=0, options=dict(SymmetricMode=True)) 
    # U_m_pred = K_LU.solve(D.dot(P_m))
    # U_m = U_m_pred
    # K_LU = SPLA.splu(K_mm, diag_pivot_thresh=0, options=dict(SymmetricMode=True)) 
    # U_m = K_LU.solve(P_m)

    U_m = SPLA.spsolve(K_mm, P_m)

    print('P_m: ', P_m)
    print('U_m: ', U_m)
    assert_array_almost_equal(K_mm.dot(U_m), P_m)
    
    # displacement
    Utmp = np.zeros(total_dof)
    Utmp[0:n_free_dof] = U_m
    U = PermT.dot(Utmp)
    cprint('Displacement:', 'green')
    for i in range(nV):
        print('N#{}: {} | {}'.format(i, U[6*i:6*i+3], U[6*i+3:6*i+6]))

    # reaction
    # print('Pf: ', P_f)
    Rtmp = np.zeros(total_dof)
    Rtmp[n_free_dof:n_free_dof+n_fixed_dof] = K_fm.dot(U_m) - P_f
    R = PermT.dot(Rtmp)
    cprint('Reaction:', 'green')
    for i in range(nV):
        print('R#{}: {} | {}'.format(i, R[6*i:6*i+3], R[6*i+3:6*i+6]))

    Ues = []
    Res = []
    for i in range(nE):
        e = T[i,:]
        # Ke = compute_Ke(e, *e_cr_y[i])
        # Ke = compute_Ke(e, return_global=False)
        Ke = k_list[i]
        R_f = compute_R(e)
        print('local axis E#{}:\n{}'.format(i, R_f[:3, :3]))
        print('dof id: ', id_map[i, :])
        Ue = U[id_map[i, :]]
        # Ues.append(R_f.dot(Ue))
        Ues.append(Ue)
        Res.append(R_f.dot(Ke.dot(Ue)))
        # Res.append(Ke.dot(Ue))

    cprint('element displacement in global axis', 'green')
    for i in range(nE):
        Ue = Ues[i]
        print('-'*5)
        print('E#{}:'.format(i))
        print('N1: {} | {}\n'.format(Ue[:3], Ue[3:6]))
        print('N2: {} | {}\n'.format(Ue[6:9], Ue[9:12]))

    cprint('local reaction', 'green')
    for i in range(nE):
        Re = Res[i]
        print('-'*5)
        print('E#{}:'.format(i))
        print('Reaction E#{}:'.format(i))
        print('N1: {} | {}\n'.format(-Re[:3], -Re[3:6]))
        print('N2: {} | {}\n'.format(Re[6:9], Re[9:12]))

# def pointup_truss_truth(setup):
#     if setup == 'beam-beam-f-f-f':
#         beam_disp = np.array([[0,0,-0.000539, 0,0,0, 
#                                0,0,0, 0,-0.000065,0],
#                               [0,0,0, 0,0.000065,0, 
#                                0,0,-0.000539, 0,0,0]])
# 
#         beam_reaction = np.array([[-11.178299,0,0.001021, 0,0,0, 
#                                    -11.178299,0,-0.001021, 0,0.01141,0, ],
#                                   [ 11.178299,0,0.001021, 0,0,0, 
#                                    -11.178299,0,0.001021, 0,0.01141,0,  ]])

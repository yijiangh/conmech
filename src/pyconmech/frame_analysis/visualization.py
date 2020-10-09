import numpy as np
import scipy
import scipy.sparse.linalg as SPLA
import numpy.linalg as LA
from numpy.linalg import norm, solve
from numpy.polynomial.polynomial import polyvander, polyval, Polynomial, polyfit

def interp_poly(d_u, d_v, L):
    """compute shape polynomial coeff for local x,y,z based on end nodal displacement value

    Parameters
    ----------
    d_u : 6x1 np array
        nodal u's displacement in local coordinate
    d_v : 6x1 np array
    L : float
        element length

    Return
    ------
    u_poly : (3x4) np array
        each row is polynomail coeff vector
    """
    D_local = np.hstack([d_u, d_v])
    dx_u, dx_v = (D_local[0], D_local[6])
    # * linear monomials 1+x
    A_l = polyvander([dx_u, L+dx_v], 1)

    # * cubic monomials 1, x, x**2, x**3
    A_c = np.zeros((4, 4))
    A_c[[0,2], :] = polyvander([dx_u, L+dx_v], 3)
    # c = np.hstack([Polynomial([1]*4).deriv().coef, [0]])
    # A[[1,3], :] = polyvander([0, L], 3) * c # apply entrywise mul on each row
    # * derivative of the monomials
    c = np.array([1,2,3])
    A_c[[1,3], 1:] = polyvander([dx_u, L+dx_v], 2) * c

    # polynomial evaluated at node pts
    # notice the negative sign here to make sign convention agree
    # coefficents of the polynomial functions under monomial basis
    d_x = [dx_u, dx_v]
    d_y = [D_local[1], D_local[5], D_local[7], D_local[11]]
    d_z = [D_local[2], -D_local[4], D_local[8], -D_local[10]]
    u_poly = [solve(A_l, d_x),
              solve(A_c, d_y),
              solve(A_c, d_z)]
    assert np.allclose(u_poly[0], polyfit([dx_u, L+dx_v], d_x, 1))
    return u_poly
    

def get_element_shape_fn(end_u, end_v, d_u, d_v, exagg=1.0):
    """cubic polynomial interpolation given its boundary condition:
        d = [u1, du1/dx, u2, du2/dx]
    
    Parameters
    ----------
    d : [type]
        [description]
    L : float
        element length

    Return
    ------
    poly_eval_fn : function handle
        shape_fn(t) -> np array of size (3,), node pts in global coordinate
        time parameter reparameterized to [0,1.0]
    """
    L = norm(end_u - end_v)
    R3 = global2local_transf_matrix(end_u, end_v)
    assert np.allclose((R3.T).dot(R3), np.eye(3))
    assert np.allclose((R3).dot(R3.T), np.eye(3))

    R = turn_diagblock(R3)
    # compute end deflection in local ordinates
    D_global = np.hstack([d_u, d_v])
    D_local = exagg * R.dot(D_global)

    dx_u, dx_v = (D_local[0], D_local[6])
    u_poly = interp_poly(D_local[:6], D_local[6:12], L)

    def shape_fn(t):
        assert t>=0 and t<=1.0
        d = np.array([polyval(dx_u + t*(L+dx_v-dx_u), u_poly[i]) for i in range(3)])
        pt_t = end_u + np.array([t*L, 0, 0])
        return (pt_t + R3.T.dot(d))
    return shape_fn


def get_internal_reaction_fn(fu, fv):
    # TODO cubic interpolation when in-span load is applied
    # * linear monomials 1+x
    A_l = polyvander([0., 1.], 1)
    # * cubic monomials 1, x, x**2, x**3
    A_c = np.zeros((4, 4))
    A_c[[0,2], :] = polyvander([0., 1.], 3)
    # * derivative of the monomials
    c = np.array([1,2,3])
    A_c[[1,3], 1:] = polyvander([0., 1.], 2) * c

    # polynomial evaluated at node pts
    # notice the negative sign here to make sign convention agree
    # coefficents of the polynomial functions under monomial basis
    # - linear interpolation for Nx, Fy, Fz, Tx
    u_poly = [solve(A_l, [fu[i], fv[i]]) for i in range(6)]
    # - cubic interpolation for moments, shear as tangents
    # Mz = [fu[5], -fu[1], fv[5], -fv[1]]
    # My = [fu[4], fu[2], fv[4], fv[2]]
    # u_poly.extend([solve(A_c, My), solve(A_c, Mz)])

    # Vz = dMy/dx 
    # assert np.allclose(-u_poly[1], Polynomial(u_poly[5]).deriv().coef)
    print((-fu[1], Polynomial(u_poly[5]).deriv().coef))
    # -Vy = dMz/dx 
    # assert np.allclose(u_poly[2], Polynomial(u_poly[4]).deriv().coef)
    print((-u_poly[2], Polynomial(u_poly[4]).deriv().coef))

    def local_force_fn(t):
        assert t>=0 and t<=1.0
        return np.array([polyval(t, u_poly[i]) for i in range(6)])
    return local_force_fn

###############################################
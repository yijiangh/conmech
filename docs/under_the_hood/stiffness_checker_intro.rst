What exactly is StiffnessChecker tring to solve?
------------------------------------------------

Elastic structure deforms under load, and by deforming themselves they develop resistance 
(or reaction) force to balance the external load. In a nutshell, the stiffness checker 
calculates elastic deformation and corresponding reaction force of a frame structure 
under given load. For example, in the context of construction sequencing, we mainly 
consider the gravity load induced by elements' self-weight in construction sequencing.

Conceptually, the solver tries to piece many elements' unit behavior together to 
reach equilibrium with the external force. Each element obeys both Hooke's law and 
the `Beam equations <https://en.wikipedia.org/wiki/Euler%E2%80%93Bernoulli_beam_theory>`__, 
which tells us how does an element **develops internal force** to balance external load 
via **deformation**.

So locally in **each beam's own local coordinate system**, we have the local elastic equation:

.. math::

  \begin{pmatrix} 
  F^{e-n1}_{L} \\ 
  --- \\ 
  F^{e-n2}_{L} 
  \end{pmatrix} := 
  \begin{pmatrix} 
  F^{n1}_{Lx} \\ 
  F^{n1}_{Ly} \\ 
  F^{n1}_{Lz} \\ 
  M^{n1}_{Lx} \\ 
  M^{n1}_{Ly} \\ 
  M^{n1}_{Lz} \\ 
  --- \\ 
  F^{n2}_{Lx} \\ 
  F^{n2}_{Ly} \\ 
  F^{n2}_{Lz} \\ 
  M^{n2}_{Lx} \\ 
  M^{n2}_{Ly} \\ 
  M^{n2}_{Lz}
  \end{pmatrix} = 
  \mathbf{K_e} \cdot 
  \begin{pmatrix} 
  d^{n1}_{Lx} \\ 
  d^{n1}_{Ly} \\ 
  d^{n1}_{Lz} \\ 
  \theta^{n1}_{Lx} \\ 
  \theta^{n1}_{Ly} \\ 
  \theta^{n1}_{Lz} \\ 
  --- \\ 
  d^{n2}_{Lx} \\ 
  d^{n2}_{Ly} \\ 
  d^{n2}_{Lz} \\ 
  \theta^{n2}_{Lx} \\ 
  \theta^{n2}_{Ly} \\ 
  \theta^{n2}_{Lz} 
  \end{pmatrix} = 
  \mathbf{K_e} 
  \begin{pmatrix} 
  u^{e-n1}_{L} \\ 
  --- \\ 
  u^{e-n2}_{L} 
  \end{pmatrix}

Here you can conceptually think about this :math:`12 \times 12` element stiffness matrix 
:math:`\mathbf{K_e}` as the stiffness factor :math:`k` in Hooke's law 
:math:`\Delta{F} = k \Delta{x}` in a string system. 
The only difference is that it's capturing the shear, bending, and torsion effect as well, 
not only the axial elongation (see picture below).

.. image:: ../images/frame_equation.png

image source: `MIT 1.571 lecture note <https://github.com/yijiangh/conmech/blob/master/docs/literature/MIT_1.571_L11_Displacement_Method.pdf>`__, page 11 (Pierre Ghisbain)

But since some elements are sharing a node, these elements' reaction must relate to 
each other so that 

1. the deformation at the shared node (:math:`u^{e-v}_{L}`) is the same
2. the reaction forces of these elements (:math:`F^{e-v}_{L}`) reach equilibrium 
   at the shared node.

Thus, at each node :math:`v`, we have first the equilibrium equation:

.. math::

  \sum_{e \in \{e | e \sim v\}} 
  \begin{pmatrix}
  R_{e, GL} & 0 \\ 
  0 & R_{e, GL}\\
  \end{pmatrix}
  \begin{pmatrix} 
  F^{v}_{e, Lx} \\ 
  F^{v}_{e, Ly} \\ 
  F^{v}_{e, Lz} \\ 
  M^{v}_{e, Lx} \\ 
  M^{v}_{e, Lx} \\ 
  M^{v}_{e, Lz} 
  \end{pmatrix} = 
  \sum_{e \in \{e | e \sim v\}} (F_{v, \textrm{e self-w load}}) + 
  F_{v, \textrm{pt load}}

Where :math:`R_{e, GL}` represents the element's local-to-global :math:`3 \times 3`
transformation matrix.

The RHS of the equation above represents all the loads (gravity, external point load) 
at node :math:`v`. I will come back to loads in the section "The relationship between 
`fixities_reaction` and `element_reaction`" in later sections. Notice that we have to 
apply this local to global rotation matrix to transform all the element internal force 
to the global coordinate system.

Then, we also have to make sure the shared nodal deformation is the same for all the 
connected elements, so by plugging into the local elastic equation above, the 
equilibrium equation becomes:

.. math::

  (\sum_{e \in \{e | e \sim v\}} 
  \begin{pmatrix}
  R_{e, GL} & 0 & 0 & 0\\ 
  0 & R_{e, GL} & 0 & 0\\
  0 & 0 & R_{e, GL} & 0\\ 
  0 & 0 & 0 & R_{e, GL}
  \end{pmatrix} \mathbf{K_e} 
  \begin{pmatrix}
  R_{e, GL} & 0 & 0 & 0\\ 
  0 & R_{e, GL} & 0 & 0\\
  0 & 0 & R_{e, GL} & 0\\ 
  0 & 0 & 0 & R_{e, GL}
  \end{pmatrix}^T) 
  \begin{pmatrix} 
  u^{n1}_{G}\\
  u^{n2}_{G}\\
  \end{pmatrix} = 
  \sum_{e \in \{e | e \sim v\}} 
  (F_{v, \textrm{e self-w load}}) + F_{v, \textrm{pt load}}

Notice that we are enforcing the nodal deformation compatibility by having all the 
connected elements share the same nodal deformation 
:math:`\begin{pmatrix} u^{n1}_{G}, u^{n2}_{G} \end{pmatrix}` **in global coordinate**.

The equation above must be satisfied by all the nodes in the structure, and we have to solve all of them together by *assembling* the global stiffness matrix.

PS: we can really see the essence of FEM here: first we have the physics model for 
one single element, then we try to enforce: 

1. internal reaction equibilirium at nodes (shared element boundary) 
2. compability on the deformation at the shared element boundary. 
3. we assembly these nodal equations together into a giant linear system and we solve.
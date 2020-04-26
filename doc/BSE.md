**Abstract** Apart from the formal paper, this document contains very detailed derivation such that it is easier to check the validity. It is organized as seven parts:

1. TDHF and its eigenvalue formulation;
2. TDHF with static $W$ and its eigenvalue formulation;
3. TDHF with dynamic $W$ and its eigenvalue formulation;
4. Hubbard model introduction;
5. Hubbard model: Hartree-Fock calculations;
6. Hubbard model: eigenvalue formulation for both static and dynamic $W$;
7. Hubbard model: time-dependent formulation for both static and first-order $W$.

# TDHF and its Eigenvalue Formulation

When a 1-particle wave function propagates under a Hamiltonian
$$
h_{\rm HF}=t+v+v_{\rm H}+k_{\rm X}
$$
where (closed-shell molecules)
$$
v_{\rm H}(r,t)\psi(r,t)=\int\mathrm dr'n(r',t)v(r,r')\psi(r',t)
$$

$$
k_{\rm X}(r,t)\psi(r,t)=-\frac12\int\mathrm dr'\rho(r,r',t)v(r,r')\psi(r',t)
$$

Then the electron-hole excitation will satisfy a eigenvalue equation (often called Casida equation), assuming all orbitals are real:
$$
\begin{pmatrix}A&B\\-B&-A\end{pmatrix}
\begin{pmatrix}X\\Y\end{pmatrix}=\omega
\begin{pmatrix}X\\Y\end{pmatrix}
$$
Note: The size of this matrix is $2N\times 2N$, where $N=N_{\rm occ}\times N_{\rm vir}$. We use $i,j,k$ for occupied and $a,b,c$ for virtual orbitals.

Specifically, $A=D+2K^X+K^D$, and $B=A-D$. Then matrix elements are given by
$$
\begin{aligned}
D_{ia,jb}&=(\varepsilon_a-\varepsilon_i)\delta_{ia,jb}\\
K^X_{ia,jb}&=\langle ia|v|jb\rangle\\
K^D_{ia,jb}&=\langle ab|v|ij\rangle
\end{aligned}
$$
(Notation: for a 2-particle operator $O_2(1,2)$, $\langle ab|O_2|cd\rangle$ means:
$$
\int \varphi_a(1)\varphi_b(1)O_2(1,2)\varphi_c(2)\varphi_d(2)\mathrm d1\mathrm d2
$$
)

# TDHF with static $W$ and its eigenvalue formulation

It is well known that the original BSE take the form
$$
L(1,2,1',2')=L_0(1,2,1',2')+L_0(1,3,1',3') \Xi(3',4',3,4)L(4,2,4',2')
$$
In the context of $GW$ approximation, we approximate $\Xi$ by
$$
\Xi(3',4',3,4)=-i\delta(3',3)\delta(4',4)v(3',4')+i\delta(3',4)\delta(4',3)W(3',4')
$$

After we transform the BSE to frequency domain, it becomes
$$
L(\omega_1,\omega_2)=L_0(\omega_1,\omega_2)+\frac{L_0(\omega_1,\omega_2)}{2\pi}\int\mathrm d\omega_3[v-W(\omega_2-\omega_3)]L(\omega_1,\omega_3)
$$
Due to the coupling in frequency, it cannot be further simplified. However, if we take $W(\omega)\approx W(\omega=0)$, we can integrate out the second frequency argument of $L$, and obtain
$$
L(\omega)=L_0(\omega)-iL_0(\omega)vL(\omega)+iL_0(\omega)WL(\omega)
$$
Then, by inserting complete sets of quasi-particle states, we are able to formulate a eigenvalue problem:
$$
\begin{pmatrix}A&B\\-B&-A\end{pmatrix}
\begin{pmatrix}X\\Y\end{pmatrix}=\omega
\begin{pmatrix}X\\Y\end{pmatrix}
$$
... while the only difference is, the direct interaction matrix elements are now given by
$$
K^D_{ia,jb}=\langle ab|W(\omega=0)|ij\rangle
$$
i.e., simply replace $v\to W(0)$.

# TDHF with static $W$ and its eigenvalue formulation

## Direct interaction matrix elements

We use $i,j,k$ for occupied and $a,b,c$ for virtual orbitals. Going beyond the static screening $W(\omega=0)$, we approximate the interaction matrix elements by
$$
K^d_{ia,jb}(\omega)=\int\mathrm dr\mathrm dr'a(r)b(r)i(r')j(r')\times\frac{i}{2\pi}\\
\int\mathrm d\omega'e^{-i\omega' 0^+}W(r,r',\omega')\left[\frac1{A-\omega'+i0^+}+\frac1{B+\omega'+i0^+}\right]
$$
Where $A=\omega-(\varepsilon_b-\varepsilon_i)$ and $B=\omega-(\varepsilon_a-\varepsilon_j)$.

## Expansion of $W(t)$

We assume that the retarded screened interaction (which has only $t>0$ parts) can be expanded into
$$
W^r(r,r',t)=\sum_l W_l(r,r')\theta(t)e^{-\omega_lt}
$$
Where $\omega_l$ is complex resonance frequency with $\Re{\omega_l}>0$. Since
$$
\int_{0}^{+\infty}e^{-\omega_lt}e^{i\omega t}=\frac{i}{\omega+i\omega_l}
$$
We have
$$
W^r(r,r',\omega)=\sum_l W_l(r,r')\frac{i}{\omega+i\omega_l}
$$
Thus the time-ordered one is
$$
W(r,r',\omega)=\begin{cases}
\sum_l W_l(r,r')\frac{i}{\omega+i\omega_l}&\omega>0\\
\sum_l W_l(r,r')\frac{-i}{\omega-i\omega_l}&\omega<0\\
\end{cases}
$$

## Expansion of $W(\omega)$

We can expand the screened interaction in plasmon frequency:
$$
W(r,r', \omega)= \sum_{l} \frac{W_{l}(r, r')}2 \left(\frac{\omega_{l}}{\omega-\omega_{l}}-\frac{\omega_{l}}{\omega+\omega_{l}}\right)
$$
We insert this equation to evaluate the frequency integral analytically, to obtain
$$
K^d_{ia,bj}(\omega)=-\int\mathrm dr\mathrm dr'a(r)b(r)i(r')j(r')\left[\sum_l\frac{W_{l}(r, r')}2 \left(\frac{\omega_{l}}{\omega_{l}-A}+\frac{\omega_{l}}{\omega_{l}-B}\right)\right]
$$
A reasonable approximation will be $A\ll\omega_l$, thus
$$
\begin{aligned}
&\sum_l-\frac{W_{l}(r, r')}2 \left(\frac{\omega_{l}}{\omega_{l}-A}+\frac{\omega_{l}}{\omega_{l}-B}\right)\\
&\approx\sum_l-\frac{W_{l}(r, r')}2\left(1+\frac{A}{\omega_l}+\frac{A^2}{\omega_l^2}+1+\frac{B}{\omega_l}+\frac{B^2}{\omega_l^2}\right)\\
&=W(r,r',0)+\frac{A+B}2X(r,r')+\frac{A^2+B^2}2Y(r,r')
\end{aligned}
$$
Where the summation of zero-order term recover the static W:
$$
W(r,r',0)=-\sum_{l}W_{l}(r, r')
$$
And the first- and second-order term
$$
\begin{aligned}
X(r,r')&=-\sum_{l}\frac{W_{l}(r, r')}{\omega_l}\\
Y(r,r')&=-\sum_{l}\frac{W_{l}(r, r')}{\omega_l^2}
\end{aligned}
$$
What can we say about $X$ and $Y$? Well, for $Y$ it is simply
$$
Y(r,r')=\frac12\frac{\partial^2W(r,r',\omega)}{\partial\omega^2}\bigg |_{\omega=0}
$$
Which can be calculated from $W(t)$; however, (we made a mistake in the afternoon such that)
$$
\frac{\partial W(r,r',\omega)}{\partial\omega}\bigg |_{\omega=0}=0
$$
Thus we cannot directly obtain $X(r,r')$, and generally $X(r,r')\ne0$. We may need specific information about plasmon frequency.

## Time domain approach

We thus consider only the first-order term and investigate how to produce $(A+B)/2$ term.
$$
\frac{A+B}2=\omega+\frac{\varepsilon_i+\varepsilon_j-\varepsilon_a-\varepsilon_b}2
$$
By mapping
$$
\omega\to i\frac{\partial}{\partial t}
$$
and
$$
\varepsilon\to \hat h_{qp}
$$
We are able to create a modified exchange operator
$$
\hat k=\hat k_0+\hat k_1+\hat k_2+\hat k_3
$$
Such that it reproduce the desired matrix elements through Casida transformation.

## $\hat k_0$

$$
\hat k_0\varphi_i(r,t)=-\frac12\int\mathrm dr'W(r,r',0)\rho(r,r',t)\varphi_i(r',t)
$$

Which have commonly been implemented in TDBSE.

## $\hat k_1$

$$
\hat k_1\varphi_i(r,t)=-\frac12\int\mathrm dr'X(r,r')i\frac{\partial\rho}{\partial t}(r,r',t)\varphi_i(r',t)
$$

Which clearly contribute to $\omega$ term. By definition,
$$
\frac{\partial\rho(r,r',t)}{\partial t}=\sum_{i\le N_{occ}}\frac{\partial}{\partial t}\varphi_i(r,t)\varphi_i^*(r',t)+\varphi_i(r,t)\frac{\partial}{\partial t}\varphi_i^*(r',t)
$$
And the time derivatives can be calculated either by backward difference or by $i\partial_t=h_{BSE}$, depending on which is easier.

## $\hat k_2$

$$
\hat k_2\varphi_i(r,t)=-\frac12\int\mathrm dr'X(r,r')\tilde\rho(r,r',t)\varphi_i(r',t)
$$

where
$$
\tilde\rho(r,r',t)=\langle r|[\hat\rho(t),\hat h_{qp}]|r'\rangle
$$
Then, we consider the excitation part $\rho(r,r',t)=\rho_0(r,r')+\rho_1(r,r',t)$, and assume the time dependence as $\rho_1(r,r',t)=z(r,r')e^{-i\omega t}+z^*(r,r')e^{i\omega t}$.
$$
\begin{aligned}
&2\langle a|\hat k_2[\rho_1]|i\rangle\\
&=\int\mathrm dr\mathrm dr'a(r)X(r,r')\tilde\rho_1(r,r',t)i(r')\\
&=\int\mathrm dr\mathrm dr'a(r)X(r,r')\langle r|[\hat z, \hat h_{qp}]|r'\rangle e^{-i\omega t} i(r')+c.c.\\
&=\sum_{bj}\int\mathrm dr\mathrm dr'(\varepsilon_j-\varepsilon_b)X(r,r')a(r)b(r)i(r')j(r') e^{-i\omega t} \langle b|\hat z|j\rangle\\
&+\cdots+c.c.
\end{aligned}
$$

Where the $\cdots$ stands for the terms contain $\langle j|\hat z|b\rangle$ and is omitted for clarity.

## $\hat k_3$

$$
\hat k_3\varphi_i(r'',t)=-\frac12\sum_nn(r'')\int\mathrm dr\mathrm dr'X(r,r')\rho(r,r',t)n(r)[(h_{qp}-\varepsilon_n)\varphi_i](r',t)
$$

So clearly,
$$
\begin{aligned}
&2\langle a|\hat k_3[\rho_1]|i\rangle\\
&=\sum_{bj}\int\mathrm dr\mathrm dr'(\varepsilon_i-\varepsilon_a)X(r,r')a(r)b(r)i(r')j(r') e^{-i\omega t} \langle b|\hat z|j\rangle\\
&+\cdots+c.c.
\end{aligned}
$$

# Introduction of the Hubbard Model

We consider an one-dimensional lattice with $n$ atoms, with atom indices $l=1\sim n$. In such a periodic system, the Wannier functions $w_l(r)$ serve as a orthonormal complete basis, thus we are able to write a second quantization form:
$$
H=\sum_{lm\sigma}T_{lm}c_{l\sigma}^{\dagger}c_{m\sigma}+\frac12\sum_{lmpq\sigma\sigma'}U_{lmpq}c_{l\sigma}^{\dagger}c_{m\sigma'}^{\dagger}c_{q\sigma'}c_{p\sigma}
$$

In Hubbard model, we make several approximations. Because ($h$ is the one-body Hamiltonian)
$$
T_{lm}=\int\mathrm dr w_l^*(r)h(r)w_m(r)
$$
Assuming that the overlap only happen when $l,m$ are adjacent atoms, and shift the energy zero s.t. $T_{ll}=0$. Similarly, because
$$
U_{lmpq}=\int \frac{w_l^*(r)w_m^*(r')w_p(r)w_q(r')}{\left|r-r^{\prime}\right|}\mathrm dr\mathrm dr'
$$
is non-vanishing only when the four grids are actually the same one, we approximate $U_{lmpq}$ by $U\delta_{lm}\delta_{lp}\delta_{lq}$, and WLOG we set $U=1$. In this way,
$$
H=\sum_{\langle lm\rangle,\sigma}T_{lm}c_{l\sigma}^{\dagger}c_{m\sigma}+\frac12 \sum_{l\sigma\sigma'}c_{l\sigma}^{\dagger}c_{l\sigma'}^{\dagger}c_{l\sigma'}c_{l\sigma}
$$
We are not allowed to create two electron of the same $i$ and $\sigma$, so $\sigma'=\sigma^*$.
$$
H=\sum_{\langle lm\rangle,\sigma}T_{lm}c_{l\sigma}^{\dagger}c_{m\sigma}+ \sum_{l}n_{l\uparrow}n_{l\downarrow}
$$

Furthermore, we use an non-periodic and alternating hopping matrix:

$$
T_{l,l+1}=
\begin{cases}
\alpha&l=1,3,5,7\\
\beta&l=2,4,6\\
\end{cases}
$$

# Hubbard Model: Hartree-Fock Calculations

## Hartree-Fock Equation for Hubbard Model

Suppose we have a closed-shell $N$-electron system (thus with occupied orbitals $N_{\rm occ}=N/2$ and virtual orbitals $N_{\rm vir}=n^2-N_{\rm occ}$), the general closed-shell Hartree-Fock equation is
$$
f(r)\varphi_i(r)=\varepsilon_i\varphi_i(r)
$$
where
$$
f(r)=h(r)+\sum_{j=1}^{N_{\rm occ}}2J_j(r)-K_j(r)
$$
and
$$
\begin{aligned}
J_j(r)\varphi_i(r)&=\int\mathrm dr'\varphi_j^*(r')\varphi_j(r')v(r,r')\varphi_i(r)\\
K_j(r)\varphi_i(r)&=\int\mathrm dr'\varphi_j^*(r')\varphi_i(r')v(r,r')\varphi_j(r)
\end{aligned}
$$

Now, for Hubbard model we introduce a Wannier basis $w_l(r)$, where $l$ is a vector from $(1, 1)$ to $(n, n)$, and write
$$
\varphi_i=\sum_{l=1}^{n^2}C_{l i}w_{l}
$$
With Fock matrix
$$
F_{lm}=\int\mathrm dr w_l^*(r)f(r)w_m(r)
$$
we are going to get
$$
\mathbf F\mathbf C=\mathbf C\boldsymbol\varepsilon
$$
this is a $n^2\times n^2$ eigenvalue problem. For convenience, we also introduce a density matrix $\mathbf P$:
$$
P_{lm}=2 \sum_{i}^{N_{\rm occ}} C_{li} C_{mi}^{*}
$$

## Evaluation of Fock matrix

The Fock matrix have one-electron and two-electron part:
$$
F_{lm}=T_{lm}+G_{lm}
$$

The $T_{lm}$ is shown above, and $G_{lm}$ is also easy to calculate:
$$
\begin{aligned}
G_{lm}&=\sum_{pq}P_{pq}\left[U_{lqmp}-\frac12U_{lqpm}\right]\\
&=\sum_{pq}P_{pq}\left[U\delta_{lm}\delta_{lq}\delta_{lp}-\frac12U\delta_{lp}\delta_{lq}\delta_{lm}\right]\\
&=\frac12U\delta_{lm}P_{ll}
\end{aligned}
$$

## SCF Procedure

1. Prepare $\mathbf T$;
2. Guess $\mathbf P=\mathbf G=0$, so $\mathbf F=\mathbf T$;
3. Solve $\mathbf F\mathbf C=\mathbf C\boldsymbol\varepsilon$;
4. Calculate $\mathbf P$ from $\mathbf C$;
5. Repeat step 3 and 4 until $|\mathbf P^{(k+1)}-\mathbf P^{(k)}|_F$ (Frobenius norm) less than some tolerance.

# Hubbard Model: Eigenvalue Formulation for Both Static and Dynamic $W$

## Static case

We adopt the RPA expression for $W$:
$$
W(1,2,\omega)=v(1,2)+\int\mathrm d(3,4)v(1,3)\chi_0(3,4,\omega)W(4,2)
$$
where
$$
\chi_0(3,4,\omega)=2\sum_{ia}\varphi_i^*(3)\varphi_a(3)\varphi_i(4)\varphi_a^*(4)f_{ia}(\omega)
$$
in which $i$ is an occupied orbital and $a$ is a virtual one; and the $f_{ia}(\omega)$ is some function that we will explicitly indicate later. But let's do the matrix elements first:
$$
\begin{aligned}
W_{lmpq}&=\int w_l^*(1)w_m^*(2)w_p(1)w_q(2)W(1,2,\omega)\mathrm dr\mathrm dr'\\
&=U\delta_{lmpq}+2\sum_{ia}f_{ia}(\omega)\int\mathrm d(1,2)w_l^*(1)\varphi_i^*(3)v(1,3)w_p(1)\varphi_a(3)\\
&\times \int\mathrm d(2,4)\varphi_a^*(4)w_m^*(2)W(4,2)\varphi_i(4)w_q(2)\\
&=U\delta_{lmpq}+2U\delta_{lp}\sum_{ia}f_{ia}(\omega)C_{li}^*C_{pa}\sum_{st}C_{sa}^*C_{ti}W_{smtq}
\end{aligned}
$$
It seems complicated, but we can reformulate it by introducing a rank-4 tensor $Q$:
$$
W_{lmpq}=U_{lmpq}+\sum_{st}Q_{lspt}W_{smtq}
$$
and formally we regard $\alpha=(l,p)$, $\beta=(m,q)$, $\gamma=(s,t)$, we can write
$$
W_{\alpha\beta}=U_{\alpha\beta}+\sum_\gamma Q_{\alpha\gamma}W_{\gamma\beta}
$$
and even more clearly,
$$
\mathbf W=\mathbf U+\mathbf Q\mathbf W
$$
i.e. $\mathbf W=\mathbf U(1-\mathbf Q)^{-1}$.

If we define $C_{li}C_{pa}=Z_{lipa}=Z_{\alpha\lambda}$, then
$$
Q_{\alpha\gamma}=2U\delta_{lp}\sum_{\lambda}Z_{\alpha\lambda}F_{\lambda}Z_{\gamma\lambda}
$$

This can be calculated easily. Now we expand the $f_{ia}$:
$$
f_{ia}=\frac1{\omega-\omega_{ia}+i0^+}-\frac1{\omega+\omega_{ia}-i0^+}
$$
and $\omega_{ia}=\varepsilon_a-\varepsilon_i$ is a (always positive) excitation energy. In the static case, it reduces to
$$
f_{ia}(0)=-\frac 2{\omega_{ia}}
$$

To summarize, the calculation is performed as following:

1. Calculate $\mathbf Q(0)$;
2. Inverse and obtain $\mathbf W(0)$;
3. Do a 4-tensor multiplication (practically two matrix-matrix product) to obtain Casida matrix;
4. Solve the eigenvalue problem.

In the Casida matrix, the exchange block is:
$$
\begin{aligned}
K^X_{ia,jb}&=\int\mathrm d(1,2)\varphi_i(1)\varphi_j(2)v(1,2)\varphi_a(1)\varphi_b(2)\\
&=\sum_{lmpq}C_{li}C_{mj}C_{pa}C_{qb}\int\mathrm d(1,2)w_l(1)w_m(2)v(1,2)w_p(1)w_q(2)\\
&=\sum_lC_{li}C_{lj}C_{la}C_{lb}U
\end{aligned}
$$

and the direct block:
$$
\begin{aligned}
K^D_{ia,jb}&=\int\mathrm d(1,2)\varphi_a(1)\varphi_i(2)W(1,2)\varphi_b(1)\varphi_j(2)\\
&=\sum_{lmpq}C_{la}C_{mi}C_{pb}C_{qj}W_{lmpq}
\end{aligned}
$$

They can also be converted to efficient matrix-form evaluation, but it is too complicated to write them out.

## Dynamic case (single plasmon)

elements by
$$
\begin{aligned}
K^D_{ia,bj}(\omega)&=\int\mathrm d(1,2)a(1)b(1)i(2)j(2)\left[\frac{W(1,2,0)}2 \left(\frac{\omega_{0}}{\omega_{0}-A}+\frac{\omega_{0}}{\omega_{0}-B}\right)\right]\\
&=K^D_{ia,bj}(0)\frac{\omega_0}{2}\left(\frac1{\omega_0-A}+\frac1{\omega_0-B}\right)
\end{aligned}
$$
Where $A=\omega-(\varepsilon_b-\varepsilon_i)$ and $B=\omega-(\varepsilon_a-\varepsilon_j)$.

## Dynamic case

The strategy above is not straightforward to be generalized to dynamic case, because a matrix inverse spoiled the plasmon structure. However, we can first get
$$
\chi(3,4,\omega)=L(3,4,3,4,\omega)=\\\sum_{\lambda}\sum_{iajb}\frac{A_{\lambda}^{(ia)}A_{\lambda}^{(jb)}\varphi_i(3)\varphi_a(3)\varphi_j(4)\varphi_b(4)}{\omega-\omega_{\lambda}+i\eta}-\frac{A_{\lambda}^{(ia)}A_{\lambda}^{(jb)}\varphi_i(3)\varphi_a(3)\varphi_j(4)\varphi_b(4)}{\omega-\omega_{\lambda}+i\eta}
$$


However, we are more interested in the case that we approximate the dynamic effect by a convolution:
$$
\begin{aligned}
\tilde f_{ia}(\omega)&=\frac{i}{2\pi}\int\mathrm d\omega' e^{-i\omega 0^+}f_{ia}(\omega')\left[\frac1{\omega-\omega'-(\varepsilon_l-\varepsilon_p)+i0^+}+\frac1{\omega+\omega'-(\varepsilon_m-\varepsilon_q)+i0^+}\right]\\
&=\frac1{\omega-\omega_{ia}+\varepsilon_l-\varepsilon_p}+\frac1{\omega-\omega_{ia}+\varepsilon_m-\varepsilon_q}
\end{aligned}
$$

Eventually, we get
$$
\begin{aligned}
K^D_{ia,jb}&=\sum_lC_{li}C_{lj}C_{la}C_{lb}+\sum_{\lambda}\sum_{kchd}\sum_{lm}C_{la}C_{lb}C_{lk}C_{lc}C_{mh}C_{md}C_{mi}C_{mj}\\&\left(\frac{A_{\lambda}^{kc}A_{\lambda}^{ld}}{\omega-\omega_{\lambda}-(\varepsilon_b-\varepsilon_i)}+\frac{B_{\lambda}^{ck}B_{\lambda}^{dl}}{\omega-\omega_{\lambda}-(\varepsilon_a-\varepsilon_j)}\right)\\&=\sum_lZ_{l\alpha}Z_{l\beta}+\sum_{\lambda}\sum_{\gamma\delta}Z_{l\alpha}Z_{l\gamma}Z_{m\delta}Z_{m\beta}\left(\frac{A_{\lambda}^{kc}A_{\lambda}^{ld}}{\omega-\omega_{\lambda}-(\varepsilon_b-\varepsilon_i)}+\frac{B_{\lambda}^{ck}B_{\lambda}^{dl}}{\omega-\omega_{\lambda}-(\varepsilon_a-\varepsilon_j)}\right)\\
&=\mathbf Z^T\mathbf Z+\sum_{\lambda}\mathbf Z^T(\mathbf Z\mathbf a_\lambda)(\mathbf Z\mathbf a_{\lambda})^T\mathbf Z\otimes\mathbf F(\lambda,b,i)+\mathbf Z^T(\mathbf Z\mathbf b_\lambda)(\mathbf Z\mathbf b_{\lambda})^T\mathbf Z\otimes\mathbf F(\lambda,a,j)
\end{aligned}
$$
The calculation is performed as following:

1. From static case obtain eigenvectors;
2. Choose a initial guess $\omega=\omega_0$, solve the equation;
3. Obtain the lowest eigenvalue $\omega_0^{(1)}$.

# Hubbard Model: Time-Dependent Formulation for Both Static and First-Order Dynamic $W$

## Preparation

The algorithm can be described as:

1. Perturb all occupied states;
2. Propagate with various Hamiltonian;
3. Calculate the density deviation and do Fourier transform.

### Perturbation

Suppose we have occupied states $\varphi_i$, where $i=1\sim N_{\rm occ}$. Within the Wannier representation, one can write
$$
\varphi_i(0)=\sum_{l}C_{li}(0)w_l
$$
The perturbation is defined as:
$$
\varphi_i^{\gamma}(0)=\sum_le^{-i\gamma l}C_{li}(0)w_l
$$
We can define $\boldsymbol C^{\gamma}$ and related $\boldsymbol P^{\gamma}$ as the perturbed orbitals and densities.
$$
C^{\gamma}_{li}(0)=e^{-i\gamma l}C_{li}(0)
$$


### Propagation

Splitting pattern:
$$
e^{ih\Delta t}=e^{iv_H\Delta t/2}e^{it\Delta t/2}e^{ik\Delta t}e^{it\Delta t/2}e^{iv_H\Delta t/2}
$$
#### One Body

$$
\begin{aligned}
\langle w_m|t|\varphi_i(t)\rangle
&=\sum_lC_{li}(t)\langle w_m|t|w_l\rangle\\
&=\sum_lC_{li}(t)T_{ml}\\
\end{aligned}
$$
#### Hartree

$$
\begin{aligned}
\langle w_m|v_H|\varphi_i(t)\rangle&=\int\mathrm drw_m^*(r)\langle r|v_H|\varphi_i(t)\rangle\\
&=\int\mathrm dr\mathrm dr'w_m^*(r)\rho(r')v(r,r')\varphi_i(r,t)\\
&=\sum_j\int\mathrm dr\mathrm dr'w_m^*(r)2\varphi_j^*(r')\varphi_j(r')v(r,r')\varphi_i(r,t)\\
&=2\sum_j\int\mathrm dr\mathrm dr'w_m^*(r)C_{mj}^*w_m^*(r')C_{mj}w_m(r')v(r,r')C_{mi}(t)w_m(r)\\
&=2\sum_jC_{mj}^*C_{mj}C_{mi}(t)\\
&=P_{mm}(t)C_{mi}(t)\\
\end{aligned}
$$

#### Exchange-Correlation (Fock)

$$
\begin{aligned}
\langle w_m|k|\varphi_i(t)\rangle&=\int\mathrm drw_m^*(r)\langle r|k|\varphi_i(t)\rangle\\
&=-\frac12\int\mathrm dr\mathrm dr'w_m^*(r)\rho(r,r',t)v(r,r')\varphi_i(r',t)\\
&=-\frac12\sum_j\int\mathrm dr\mathrm dr'w_m^*(r)2\varphi_j^*(r')\varphi_j(r)v(r,r')\varphi_i(r',t)\\
&=-\sum_j\int\mathrm dr\mathrm dr'w_m^*(r)C_{mj}^*w_m^*(r')C_{mj}w_m(r)v(r,r')C_{mi}(t)w_m(r')\\
&=-\sum_jC_{mj}^*C_{mj}C_{mi}(t)\\
&=-\frac12P_{mm}(t)C_{mi}(t)\\
\end{aligned}
$$

#### Exchange-Correlation (Static BSE)

$$
\begin{aligned}
\langle w_l|k|\varphi_i(t)\rangle
&=-\frac12\sum_j\int\mathrm dr\mathrm dr'w_l^*(r)2\varphi_j^*(r')\varphi_j(r)W(r,r',0)\varphi_i(r',t)\\
&=-\sum_j\sum_{mpq}C_{mj}^*(t)C_{pj}(t)C_{qi}(t)\int\mathrm dr\mathrm dr'w_l^*(r)w_m^*(r')W(r,r',0)w_p(r)w_q(r')\\
&=-\sum_j\sum_{mpq}C_{mj}^*C_{pj}C_{qi}W_{lmpq}\\
&=-\frac12\sum_{mpq}P_{pm}C_{qi}W_{lmpq}\\
&=-\frac12\sum_{q}F_{lq}C_{qi}\\
\end{aligned}
$$



#### Exchange-Correlation (Dynamic BSE)



### Correlation Function

$$
\begin{aligned}
d(t)&=\frac1\gamma\int\mathrm drz(n^\gamma(r,t)-n^0(r,t))\\
&=\frac1\gamma\int\mathrm drz2\sum_{i=1}^{N_{\rm occ}}\varphi_i^{\gamma,*}(r,t)\varphi_i^{\gamma}(r,t)-\varphi_i^{0,*}(r,t)\varphi_i^{0}(r,t)\\
&=\frac2\gamma\sum_{i=1}^{N_{\rm occ}}\sum_{lm}\int\mathrm drz\left[C^{\gamma,*}_{li}(t)C^{\gamma}_{mi}(t)-C^{0,*}_{li}(t)C^{0}_{mi}(t)\right]w_l(r)w_m(r)\\
&=\frac1\gamma\sum_{lm}\left[P_{lm}^{\gamma}(t)-P_{lm}^0(t)\right]l\delta_{lm}\\
\end{aligned}
$$

# Hubbard model: CI

## One-body and Two-body matrix

Suppose we have some Hartree-Fock orbital labeled $i,j,k,l$. Then
$$
\langle i|h|j\rangle=\sum_{lm}C_{li}^*C_{mj}\langle l|h|m\rangle=\sum_{lm}C_{li}^*C_{mj}T_{lm}
$$
So,
$$
\boldsymbol H_1=\boldsymbol C^{\dagger}\boldsymbol T\boldsymbol C
$$
As well as
$$
\langle ij|kh\rangle=\sum_{lmpq}C_{li}C_{mj}C_{pk}C_{qh}\langle lm|pq\rangle=\sum_lC_{li}C_{lj}C_{lk}C_{lh}
$$

## Check CI result

A valid optical excitation needs
$$
\langle m|\hat Z|n\rangle\ne0
$$

$$
|n\rangle=\sum_cD_{cn}|\Phi_c\rangle
$$

so we basically need
$$
\sum_{bc}D_{bm}^*D_{cn}\langle\Phi_b|\hat Z|\Phi_c\rangle\ne0
$$
Now,
$$
\langle i|\hat z|j\rangle=\sum_{lm}C_{li}^*C_{mj}\langle w_l|\hat z|w_m\rangle=\sum_lC_{li}^*C_{lj}l
$$
This matrix is $\boldsymbol C^\dagger\boldsymbol L\boldsymbol C$.




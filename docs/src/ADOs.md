# [Auxiliary Density Operators](@id doc-ADOs)

## Introduction
The auxiliary density operators ([`ADOs`](@ref)) ``\rho_{\textbf{j}\vert\textbf{q}}^{(m,n,p)}(t)`` encode environmental effects related to different exponential terms ([`Exponent`](@ref)) present in the [`Bosonic Bath`](@ref doc-Bosonic-Bath) and [`Fermionic Bath`](@ref doc-Fermionic-Bath) correlation functions and provide an iterative description of high-order system-baths memory effects.

In ``\rho_{\textbf{j}\vert\textbf{q}}^{(m,n,p)}(t)``, the tuple ``(m, n, p)`` represents the ``m``th-level-bosonic-and-``n``th-level-fermionic ADO with parity ``p``, and ``\textbf{j}`` (``\textbf{q}``) denotes a vector ``[j_1,\cdot\cdot\cdot,j_m]`` (``[q_1,\cdot\cdot\cdot,q_n]``) where each ``j`` (``q``) represents a specific multi-index ensemble ``\{\beta, u, l, \sigma_\textrm{b}\}`` (``\{\alpha, \nu, h, \sigma_\textrm{f}\}``) with
 - ``\beta`` : denotes the index of bosonic bath
 - ``\alpha`` : denotes the index of fermionic bath
 - ``u`` : denotes the real (``u=\textrm{R}``) and imaginary (``u=\textrm{I}``) part of bosonic bath correlation function
 - ``\nu`` : denotes the fermionic bath correaltion function corresponding to absorption (``\nu=+``) and emission (``\nu=-``) process
 - ``l`` : denotes the index of exponent in the bosonic bath correlation function ``C_\beta^u(\tau)=\sum_l \eta_{\beta, l}^u \exp(-\gamma_{\beta,l}^u)(\tau)``
 - ``h`` : denotes the index of exponent in the fermionic bath correlation function ``C_\alpha^\nu(\tau)=\sum_h \eta_{\alpha, h}^\nu \exp(-\gamma_{\alpha,h}^\nu)(\tau)``
 - ``\sigma_\textrm{b}`` : specify other quantum numbers (such as energy or spin) of the system interacting with a bosonic bath
 - ``\sigma_\textrm{f}`` : specify other quantum numbers (such as energy or spin) of the system interacting with a fermionic bath

## Methods
using HierarchicalEOM
using SparseArrays

e = -5
U = 10
d_up = kron(      [0 1; 0 0], [1 0; 0 1])
d_dn = kron(-1 * [1 0; 0 -1], [0 1; 0 0])
iden = kron(      [1 0; 0 1], [1 0; 0 1])

H0 = e * (d_up' * d_up + d_dn' * d_dn)
H1 = U * (d_up' * d_up * d_dn' * d_dn)
Hsys = H0 + H1

λ   =  1
μ_l =  1
μ_r = -1
W   = 10
kT  =  0.5
N = 5
fuL = Fermion_Lorentz_Pade(d_up, λ, μ_l, W, kT, N)
fdL = Fermion_Lorentz_Pade(d_dn, λ, μ_l, W, kT, N)
fuR = Fermion_Lorentz_Pade(d_up, λ, μ_r, W, kT, N)
fdR = Fermion_Lorentz_Pade(d_dn, λ, μ_r, W, kT, N)

tier = 2
Le = M_Fermion(Hsys, tier, [fuL, fdL, fuR, fdR]; verbose=false)
Lo = M_Fermion(Hsys, tier, [fuL, fdL, fuR, fdR], ODD; verbose=false)

ados_s = SteadyState(Le; verbose=false)
ωlist = -20:2:20

if isfile("DOS.txt")
    rm("DOS.txt")
end
d_up_normal = HEOMSuperOp(d_up, ODD, Le)
dos1 = DensityOfStates(Lo, ados_s, d_up, ωlist; verbose=false, filename="DOS")
dos2 = PowerSpectrum(Lo, ados_s, d_up_normal, d_up', ωlist, true; verbose=false) .+ PowerSpectrum(Lo, ados_s, d_up', d_up_normal, ωlist, false; verbose=false)
dos3 = [
    0.0007920428534358747,
    0.0012795202828027256,
    0.0022148985361417936,
    0.004203651086852703,
    0.009091831276694192,
    0.024145570990254425,
    0.09111929563034224,
    0.2984323182857298,
    0.1495897397559431,
    0.12433521300531171,
    0.17217519700362036,
    0.1243352130053117,
    0.14958973975594306,
    0.29843231828572964,
    0.09111929563034214,
    0.024145570990254432,
    0.009091831276694183,
    0.004203651086852702,
    0.0022148985361417923,
    0.0012795202828027236,
    0.0007920428534358735
]
for i in 1:length(ωlist)
    @test dos1[i] ≈ dos3[i] atol=1.0e-10
    @test dos2[i] ≈ dos3[i] atol=1.0e-10
end

mat  = spzeros(ComplexF64, 2, 2)
mat2 = spzeros(ComplexF64, 3, 3)
bathb = Boson_DrudeLorentz_Pade(mat, 1, 1, 1, 2)
@test_throws ErrorException spectrum(Lo, ados_s, d_up, ωlist; verbose=false)
@test_throws ErrorException DensityOfStates(Lo, ados_s, mat2, ωlist; verbose=false)
@test_throws ErrorException DensityOfStates(Lo, ADOs(zeros(8), 2), d_up, ωlist; verbose=false)
@test_throws ErrorException DensityOfStates(Lo, ADOs(zeros(32), 2), d_up, ωlist; verbose=false)
@test_throws ErrorException DensityOfStates(Lo, ADOs(ados_s.data, ados_s.N, ODD), d_up, ωlist; verbose=false)
@test_throws ErrorException DensityOfStates(M_Boson(mat, 2, bathb; verbose=false), mat, mat, [0])
@test_throws ErrorException DensityOfStates(M_Fermion(mat, 2, fuL; verbose=false), mat, mat, [0])

# remove all the temporary files
rm("DOS.txt")
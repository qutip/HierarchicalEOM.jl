export BosonDynamicalField
export AbstractBosonDynamicalField, AbstractBosonFunctionField,
    bosonInputFunction, bosonOutputLeft, bosonOutputRight, bosonOutputFunctionLeft, bosonOutputFunctionRight

abstract type AbstractBosonDynamicalField <: AbstractBosonBath end
abstract type AbstractBosonFunctionField <: AbstractBosonDynamicalField end

function _check_dynamical_field_function(ηlist::Vector{Function})
    Nterm = length(ηlist)
    ηnew = Vector{ScalarOperator}(undef, Nterm)
    for i in 1:Nterm
        f_method = methods(ηlist[i], [Any, Real])
        length(f_method.ms) == 0 &&
            error("The following function must accept two arguments: `$(f_method.mt.name)(p, t)` with t<:Real")

        update_func = (a, u, p, t) -> ηlist[i](p, t)
        ηnew[i] = ScalarOperator(zero(ComplexF64), update_func)
    end
    return ηnew
end

function BosonDynamicalField(
    op::QuantumObject;
    η_in = nothing,
    η_out_L = [],
    γ_out_L = [],
    η_out_R = [],
    γ_out_R = [],
    η_out_fn_L = nothing,
    η_out_fn_R = nothing,
    δ::Number = 0.0,
)
    bath_list = AbstractBosonBath[]

    (η_in isa Nothing) || push!(bath_list, bosonInputFunction(op, η_in))

    N_out_L = length(η_out_L)
    if N_out_L == length(γ_out_L)
        (N_out_L != 0) && push!(bath_list, bosonOutputLeft(op, η_out_L, γ_out_L))
    else
        error("The length of \'η_out_L\' and \'γ_out_L\' should be the same.")
    end

    N_out_R = length(η_out_R)
    if N_out_R == length(γ_out_R)
        (N_out_R != 0) && push!(bath_list, bosonOutputRight(op, η_out_R, γ_out_R))
    else
        error("The length of \'η_out_R\' and \'γ_out_R\' should be the same.")
    end

    (η_out_fn_L isa Nothing) || push!(bath_list, bosonOutputFunctionLeft(op, η_out_fn_L))
    (η_out_fn_R isa Nothing) || push!(bath_list, bosonOutputFunctionRight(op, η_out_fn_R))

    Nterm = (length(bath_list) == 0) ? 0 : sum(getfield.(bath_list, :Nterm))
    return BosonBath(bath_list, op, Nterm, δ)
end

struct bosonInputFunction <: AbstractBosonFunctionField
    Comm::SparseMatrixCSC{ComplexF64,Int64}
    dimensions::Dimensions
    η::Vector{ScalarOperator}
    γ::AbstractVector
    Nterm::Int

    function bosonInputFunction(op::QuantumObject, η::Vector{Function})
        _η = _check_dynamical_field_function(η)
        _op = _check_bosonic_coupling_operator(op)

        Nterm = length(_η)
        γ = zeros(Nterm)
        Id_cache = I(size(_op, 1))
        return new(_spre(_op.data, Id_cache) - _spost(_op.data, Id_cache), _op.dimensions, _η, γ, Nterm)
    end
end
bosonInputFunction(op::QuantumObject, η::Function) = bosonInputFunction(op, Function[η])

struct bosonOutputLeft <: AbstractBosonDynamicalField
    spre::SparseMatrixCSC{ComplexF64,Int64}
    dimensions::Dimensions
    η::AbstractVector
    γ::AbstractVector
    Nterm::Int

    function bosonOutputLeft(op::QuantumObject, η::Vector{Ti}, γ::Vector{Tj}) where {Ti<:Number,Tj<:Number}
        _op = _check_bosonic_coupling_operator(op)

        N_exp_term = length(η)
        (N_exp_term == length(γ)) || error("The length of \'η\' and \'γ\' should be the same.")

        Id_cache = I(size(_op, 1))
        return new(_spre(_op.data, Id_cache), _op.dimensions, η, γ, N_exp_term)
    end
end

struct bosonOutputRight <: AbstractBosonDynamicalField
    spost::SparseMatrixCSC{ComplexF64,Int64}
    dimensions::Dimensions
    η::AbstractVector
    γ::AbstractVector
    Nterm::Int

    function bosonOutputRight(op::QuantumObject, η::Vector{Ti}, γ::Vector{Tj}) where {Ti<:Number,Tj<:Number}
        _op = _check_bosonic_coupling_operator(op)

        N_exp_term = length(η)
        (N_exp_term == length(γ)) || error("The length of \'η\' and \'γ\' should be the same.")

        Id_cache = I(size(_op, 1))
        return new(_spost(_op.data, Id_cache), _op.dimensions, η, γ, N_exp_term)
    end
end

struct bosonOutputFunctionLeft <: AbstractBosonFunctionField
    spre::SparseMatrixCSC{ComplexF64,Int64}
    dimensions::Dimensions
    η::Vector{ScalarOperator}
    γ::AbstractVector
    Nterm::Int

    function bosonOutputFunctionLeft(op::QuantumObject, η::Vector{Function})
        _η = _check_dynamical_field_function(η)
        _op = _check_bosonic_coupling_operator(op)

        Nterm = length(_η)
        γ = zeros(Nterm)
        Id_cache = I(size(_op, 1))
        return new(_spre(_op.data, Id_cache), _op.dimensions, _η, γ, Nterm)
    end
end
bosonOutputFunctionLeft(op::QuantumObject, η::Function) = bosonOutputFunctionLeft(op, Function[η])

struct bosonOutputFunctionRight <: AbstractBosonFunctionField
    spost::SparseMatrixCSC{ComplexF64,Int64}
    dimensions::Dimensions
    η::Vector{ScalarOperator}
    γ::AbstractVector
    Nterm::Int

    function bosonOutputFunctionRight(op::QuantumObject, η::Vector{Function})
        _η = _check_dynamical_field_function(η)
        _op = _check_bosonic_coupling_operator(op)

        Nterm = length(_η)
        γ = zeros(Nterm)
        Id_cache = I(size(_op, 1))
        return new(_spost(_op.data, Id_cache), _op.dimensions, _η, γ, Nterm)
    end
end
bosonOutputFunctionRight(op::QuantumObject, η::Function) = bosonOutputFunctionRight(op, Function[η])

function Base.getproperty(b::BType, key::Symbol) where {BType<:AbstractBosonDynamicalField}
    # a comment here to avoid bad render by JuliaFormatter
    if key === :dims
        return dimensions_to_dims(getfield(b, :dimensions))
    else
        return getfield(b, key)
    end
end

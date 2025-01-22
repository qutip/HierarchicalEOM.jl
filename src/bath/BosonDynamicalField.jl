export AbstractBosonDynamicalField,
    bosonInputFunction, bosonOutputLeft, bosonOutputRight, bosonOutputFunctionLeft, bosonOutputFunctionRight

abstract type AbstractBosonDynamicalField <: AbstractBosonBath end

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

struct bosonInputFunction <: AbstractBosonDynamicalField
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

struct bosonOutputFunctionLeft <: AbstractBosonDynamicalField
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

struct bosonOutputFunctionRight <: AbstractBosonDynamicalField
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

function Base.getproperty(b::BType, key::Symbol) where {BType<:AbstractBosonDynamicalField}
    # a comment here to avoid bad render by JuliaFormatter
    if key === :dims
        return dimensions_to_dims(getfield(b, :dimensions))
    else
        return getfield(b, key)
    end
end

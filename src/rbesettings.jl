struct RBESettings{T}
    g_tol::Real
    x_tol::Real
    f_tol::Real
    iterations::Int
    store_trace::Bool
    extended_trace::Bool
    show_trace::Bool
    memopt::Bool
    init::Vector{T}
    postopt::Bool
    vlm::Real
    maxopttry::Int
    rhoadjstep::Real
    rholink::Symbol
    singlim::Real
end

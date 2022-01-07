struct RBESettings{T}
    g_tol::Float64
    x_tol::Float64
    f_tol::Float64
    iterations::Int
    store_trace::Bool
    extended_trace::Bool
    show_trace::Bool
    memopt::Bool
    init::Vector{T}
    postopt::Bool
    vlm::Float64
    maxopttry::Int
    rholink::Symbol
    singlim::Float64
    function RBESettings(;g_tol = 1e-12, x_tol = 0.0, f_toll = 0.0, iterations = 100,
        store_trace = false, extended_trace = false, show_trace = false,
        memopt = true,
        init = zeros(0),
        postopt = false, vlm = 1.0, maxopttry = 100,
        rholink = :psigmoid,
        singlim = 1e-10)
        new{typeof(init)}(g_tol, x_tol, f_tol, iterations,
            store_trace, extended_trace, show_trace,
            memopt,
            init,
            postopt, vlm, maxopttry,
            rholink,
            singlim)
    end
end

ENV["GKSwstype"] = "100"

using Dash
using Plots
using Base64
using LaTeXStrings

gr()

# Load the module
include("qwz_pbc_module.jl")
using .QWZ_Model

# Helper function to convert Julia plot to base64 image string for the web
function plot_to_b64(p)
    io = IOBuffer()
    show(io, MIME("image/png"), p)
    return "data:image/png;base64," * base64encode(take!(io))
end

app = dash()

app.layout = html_div(style=Dict("display" => "flex", "font-family" => "Arial", "minHeight" => "100vh", "margin" => "0")) do
    
    # Left Sidebar (1/5 Width ~ 20%)
    html_div(style=Dict(
        "width" => "20%",
        "backgroundColor" => "#f8f9fa",
        "padding" => "20px",
        "boxSizing" => "border-box",
        "borderRight" => "1px solid #e0e0e0",
        "display" => "flex",
        "flexDirection" => "column",
        "gap" => "15px"
    )) do
        html_h3("Parameters", style=Dict("marginTop" => "0", "marginBottom" => "10px")),
        
        # Grid Resolution Inputs
        html_div() do
            html_label("Nkx (kx Resolution)", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_input(id="input-Nkx", type="number", value=101, min=10, max=300, step=1, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,

        html_div() do
            html_label("Nky (ky Resolution)", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_input(id="input-Nky", type="number", value=101, min=10, max=300, step=1, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,

        html_hr(style=Dict("margin" => "5px 0", "border" => "none", "borderTop" => "1px solid #ddd")),

        # Model Parameters
        html_div() do
            html_label("A", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_input(id="input-A", type="number", value=1.0, step=0.1, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,
        
        html_div() do
            html_label("B", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_input(id="input-B", type="number", value=1.0, step=0.1, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,
        
        html_div() do
            html_label("m", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_input(id="input-m", type="number", value=-1.0, step=0.1, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,

        html_div() do
            html_label("Gamma", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_input(id="input-gamma", type="number", value=0.0, step=0.1, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,

        html_div() do
            html_label("Chern crit", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_input(id="input-contour_value", type="number", value=0.5, step=0.1, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,

        html_div() do
            html_label("Distortion Type", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_dropdown(
                id="dropdown-pert",
                options=[
                    Dict("label" => "None", "value" => "none"),
                    Dict("label" => "Symmetric", "value" => "symmetric"),
                    Dict("label" => "Tilt", "value" => "tilt")
                ],
                value="symmetric",
                style=Dict("marginTop" => "4px")
            )
        end
    end,

    # Right Main Panel (4/5 Width ~ 80%)
    html_div(style=Dict(
        "width" => "80%",
        "padding" => "20px",
        "boxSizing" => "border-box",
        "display" => "flex",
        "flexDirection" => "column",
        "alignItems" => "center",
        "gap" => "20px"
    )) do
        html_h2("QWZ Model - Topological Band structure & Berry Curvature", style=Dict("marginTop" => "0")),

        # Plots area wrapped in loaders
        dcc_loading(id="loading-1", type="circle", children=html_img(id="plot-bandstruct", style=Dict("maxWidth" => "100%", "height" => "auto"))),
        dcc_loading(id="loading-2", type="circle", children=html_img(id="plot-fxy", style=Dict("maxWidth" => "100%", "height" => "auto"))),
        dcc_loading(id="loading-3", type="circle", children=html_img(id="plot-cum", style=Dict("maxWidth" => "100%", "height" => "auto")))
    end
end

# Callback linking typed inputs to the plots
callback!(app,
    Output("plot-bandstruct", "src"),
    Output("plot-fxy", "src"),
    Output("plot-cum", "src"),
    Input("input-Nkx", "value"),
    Input("input-Nky", "value"),
    Input("input-A", "value"),
    Input("input-B", "value"),
    Input("input-m", "value"),
    Input("input-gamma", "value"),
    Input("input-contour_value", "value"),
    Input("dropdown-pert", "value")
) do Nkx, Nky, A, B, m, gamma, contour_value, pert
    
    # Guard against invalid or missing input values during typing
    nkx_val = isnothing(Nkx) || Nkx < 2 ? 20 : Int(Nkx)
    nky_val = isnothing(Nky) || Nky < 2 ? 20 : Int(Nky)
    a_val = isnothing(A) ? 1.0 : Float64(A)
    b_val = isnothing(B) ? 1.0 : Float64(B)
    m_val = isnothing(m) ? 0.0 : Float64(m)
    g_val = isnothing(gamma) ? 0.0 : Float64(gamma)
    c_val = isnothing(contour_value) ? 0.5 : Float64(contour_value)
    pert_type = isnothing(pert) ? :none : Symbol(pert)

    # 1. Run the bulk computation
    data = compute_bulk_band_berry_data(
        Nkx=nkx_val,
        Nky=nky_val, 
        A=a_val, 
        B=b_val, 
        m=m_val, 
        gamma=g_val, 
        perturbation_type=pert_type
    )

    title_str = LaTeXString("\\gamma = $g_val")
    
    # 2. Render plots
    p1 = plt_bandstructure_heatmap(data.kx_vals, data.ky_vals, data.energies; title=title_str)
    p2 = plt_k_resolved_F_xy_heatmaps(data.kx_vals, data.ky_vals, data.berry_curvature; title=title_str)
    p3 = plt_accumulated_chern_heatmaps(data.kx_vals, data.ky_vals, data.cum_chern_per_band; title=title_str, contour_value=c_val)

    # 3. Convert to base64 images for HTML rendering
    return plot_to_b64(p1), plot_to_b64(p2), plot_to_b64(p3)
end

println("Starting Interactive Dashboard... Open http://127.0.0.1:8050 in your browser.")
run_server(app, "0.0.0.0", 8050)
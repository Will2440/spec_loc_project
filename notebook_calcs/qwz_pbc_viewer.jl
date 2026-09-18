ENV["GKSwstype"] = "100"

using Dash
using Plots
using Base64
using LaTeXStrings
using Printf

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
            dcc_input(id="input-Nkx", type="number", value=101, min=10, max=501, step=1, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,

        html_div() do
            html_label("Nky (ky Resolution)", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_input(id="input-Nky", type="number", value=101, min=10, max=501, step=1, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,

        html_hr(style=Dict("margin" => "5px 0", "border" => "none", "borderTop" => "1px solid #ddd")),

        # Model Parameters
        html_div() do
            html_label("A", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_input(id="input-A", type="number", value=1.0, step=0.01, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,
        
        html_div() do
            html_label("B", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_input(id="input-B", type="number", value=1.0, step=0.01, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,
        
        html_div() do
            html_label("m", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_input(id="input-m", type="number", value=-1.0, step=0.01, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,

        html_div() do
            html_label("Gamma", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_input(id="input-gamma", type="number", value=0.0, step=0.01, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,

        html_div() do
            html_label("Indirect Distortion Type (H=H+g(k)I)", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),            dcc_dropdown(
                id="dropdown-indirect-distortion",
                options=[
                    Dict("label" => "Symmetric (cos)", "value" => "", "disabled" => true),
                    Dict("label" => "  cos(kx) + cos(ky)", "value" => "sym_cos_sum"),
                    Dict("label" => "  cos(kx) - cos(ky)", "value" => "sym_cos_diff"),
                    Dict("label" => "  cos(kx + ky)", "value" => "sym_cos_add"),
                    Dict("label" => "  cos(kx - ky)", "value" => "sym_cos_sub"),
                    Dict("label" => "Antisymmetric (sin)", "value" => "", "disabled" => true),
                    Dict("label" => "  sin(kx) + sin(ky)", "value" => "asym_sin_sum"),
                    Dict("label" => "  sin(kx) - sin(ky)", "value" => "asym_sin_diff"),
                    Dict("label" => "  sin(kx + ky)", "value" => "asym_sin_add"),
                    Dict("label" => "  sin(kx - ky)", "value" => "asym_sin_sub")
                ],
                value="sym_cos_sum",
                style=Dict("marginTop" => "4px")
            )
        end,

        html_hr(style=Dict("margin" => "5px 0", "border" => "none", "borderTop" => "1px solid #ddd")),

        html_div() do
            html_label("Orbital Embedding Angle (φ=aπ)", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_input(id="input-phi", type="number", value=0.0, step=0.05, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,

        html_div() do
            html_label("Orbital Displacement (d=a*sqrt(2))", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_input(id="input-d", type="number", value=0.0, step=1.0, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,

        html_hr(style=Dict("margin" => "5px 0", "border" => "none", "borderTop" => "1px solid #ddd")),

        html_div() do
            html_label("Fermi Energy (E_f)", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_input(id="input-fermi-energy", type="number", value=0.0, step=0.05, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,

        html_div() do
            html_label("Fermi Surface Type", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_dropdown(
                id="dropdown-fermi-mode",
                options=[
                    Dict("label" => "None", "value" => "none"),
                    Dict("label" => "Contour", "value" => "contour"),
                    Dict("label" => "Mask", "value" => "mask")
                ],
                value="none",
                style=Dict("marginTop" => "4px")
            )
        end,

        html_div() do
            dcc_checklist(
                id="checklist-fermi-surface",
                options=[
                    Dict("label" => " Show on Bandstructure", "value" => "bandstructure"),
                    Dict("label" => " Show on Berry Flux", "value" => "berry_flux")
                ],
                value=["bandstructure"],
                style=Dict("marginTop" => "8px", "fontSize" => "13px")
            )
        end,

        html_hr(style=Dict("margin" => "5px 0", "border" => "none", "borderTop" => "1px solid #ddd")),

        html_div() do
            html_label("Chern crit 1", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_input(id="input-contour_value_1", type="number", value=0.5, step=0.01, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,

        html_div() do
            html_label("Chern crit 2", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_input(id="input-contour_value_2", type="number", value=0.0, step=0.01, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,

        html_div() do
            html_label("Chern crit tol", style=Dict("fontWeight" => "bold", "fontSize" => "14px")),
            dcc_input(id="input-contour_tol", type="number", value=0.01, step=0.001, style=Dict("width" => "100%", "padding" => "6px", "marginTop" => "4px"))
        end,


        html_hr(style=Dict("margin" => "15px 0", "border" => "none", "borderTop" => "2px solid #999")),

        html_h3("Properties", style=Dict("marginTop" => "10px", "marginBottom" => "10px")),
        
        html_div() do
            dcc_checklist(
                id="checklist-properties",
                options=[
                    Dict("label" => " Euler Characteristic", "value" => "euler_char"),
                    Dict("label" => " Accumulated Chern C(E_F)", "value" => "chern_accum"),
                    Dict("label" => " Energy of Crit Chern", "value" => "energy_crit_chern")
                ],
                value=[],
                style=Dict("fontSize" => "13px")
            )
        end,

        html_div(id="properties-info", style=Dict("marginTop" => "15px", "fontSize" => "12px", "color" => "#666"))
    end,

    # Center Main Panel (60% Width ~ 60%)
    html_div(style=Dict(
        "width" => "60%",
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
    end,

    # Right Properties Panel (20% Width ~ 20%)
    html_div(style=Dict(
        "width" => "20%",
        "backgroundColor" => "#f0f4f8",
        "padding" => "20px",
        "boxSizing" => "border-box",
        "borderLeft" => "1px solid #e0e0e0",
        "display" => "flex",
        "flexDirection" => "column",
        "gap" => "15px",
        "overflowY" => "auto"
    )) do
        html_h3("Calculated Properties", style=Dict("marginTop" => "0", "marginBottom" => "15px")),
        html_div(id="right-properties-panel", style=Dict("fontSize" => "13px", "lineHeight" => "1.6"))
    end
end

# Callback linking typed inputs to the plots
callback!(app,
    Output("plot-bandstruct", "src"),
    Output("plot-fxy", "src"),
    Output("plot-cum", "src"),
    Output("right-properties-panel", "children"),
    Input("input-Nkx", "value"),
    Input("input-Nky", "value"),
    Input("input-A", "value"),
    Input("input-B", "value"),
    Input("input-m", "value"),
    Input("input-gamma", "value"),
    Input("input-phi", "value"),
    Input("input-d", "value"),
    Input("input-fermi-energy", "value"),
    Input("dropdown-fermi-mode", "value"),
    Input("checklist-fermi-surface", "value"),
    Input("input-contour_value_1", "value"),
    Input("input-contour_value_2", "value"),
    Input("input-contour_tol", "value"),
    Input("dropdown-indirect-distortion", "value"),
    Input("checklist-properties", "value")
) do Nkx, Nky, A, B, m, gamma, phi, d, fermi_energy, fermi_mode, fermi_surface_checks, contour_value_1, contour_value_2, contour_tol, distortion_type, properties_checked
    
    # Guard against invalid or missing input values during typing
    nkx_val = isnothing(Nkx) || Nkx < 2 ? 20 : Int(Nkx)
    nky_val = isnothing(Nky) || Nky < 2 ? 20 : Int(Nky)
    a_val = isnothing(A) ? 1.0 : Float64(A)
    b_val = isnothing(B) ? 1.0 : Float64(B)
    m_val = isnothing(m) ? 0.0 : Float64(m)
    g_val = isnothing(gamma) ? 0.0 : Float64(gamma)
    phi_input = isnothing(phi) ? 0.0 : Float64(phi)
    phi_val = phi_input * π  # Convert to radians
    d_input = isnothing(d) ? 0.0 : Float64(d)
    d_val = d_input * 1
    ef_val = isnothing(fermi_energy) ? 0.0 : Float64(fermi_energy)
    c_val_1 = isnothing(contour_value_1) ? 0.5 : Float64(contour_value_1)
    c_val_2 = isnothing(contour_value_2) ? 0.5 : Float64(contour_value_2)
    c_tol = isnothing(contour_tol) ? 0.01 : Float64(contour_tol)
    
    # Convert distortion type string to symbol
    pert_type = isnothing(distortion_type) || distortion_type == "" ? :sym_cos_sum : Symbol(distortion_type)
    
    properties_list = isnothing(properties_checked) ? [] : properties_checked

    show_contour = (fermi_mode == "contour")
    mask_surface = (fermi_mode == "mask")
    
    # Determine which plots show fermi surface
    checks = isnothing(fermi_surface_checks) ? [] : fermi_surface_checks
    show_fermi_bandstruct = "bandstructure" in checks
    show_fermi_berry = "berry_flux" in checks

    # 1. Run the bulk computation
    data = compute_bulk_band_berry_data(
        Nkx=nkx_val,
        Nky=nky_val, 
        A=a_val, 
        B=b_val, 
        m=m_val, 
        gamma=g_val,
        phi=phi_val,
        orbital_displacement=d_val,
        perturbation_type=pert_type
    )

    title_str = if pert_type == :sym_cos_sum || pert_type == :sym_cos_diff || pert_type == :sym_cos_add || pert_type == :sym_cos_sub
        func_type = split(string(pert_type), "_")[3]  # Extract "sum", "diff", "add", or "sub"
        LaTeXString("\\gamma = $g_val, cos($func_type), \\phi = $(phi_input)\\pi, d = $d_val")
    elseif pert_type == :asym_sin_sum || pert_type == :asym_sin_diff || pert_type == :asym_sin_add || pert_type == :asym_sin_sub
        func_type = split(string(pert_type), "_")[3]  # Extract "sum", "diff", "add", or "sub"
        LaTeXString("\\gamma = $g_val, sin($func_type), \\phi = $(phi_input)\\pi, d = $d_val")
    else
        LaTeXString("\\phi = $(phi_input)\\pi, d = $d_val")
    end
    
    # 2. Render plots
    p1 = plt_bandstructure_fermi_surface_heatmap(
        data.kx_vals, data.ky_vals, data.energies; 
        title=title_str,
        fermi_energy=ef_val,
        show_fermi_contour=(show_fermi_bandstruct && show_contour),
        mask_fermi_surface=(show_fermi_bandstruct && mask_surface)
    )
    p2 = plt_k_resolved_F_xy_heatmaps(
        data.kx_vals, data.ky_vals, data.berry_curvature; 
        title=title_str,
        energies=data.energies,
        fermi_energy=ef_val,
        show_fermi_contour=(show_fermi_berry && show_contour),
        mask_fermi_surface=(show_fermi_berry && mask_surface)
    )
    p3 = plt_accumulated_chern_heatmaps(data.kx_vals, data.ky_vals, data.cum_chern_per_band; title=title_str, contour_value_1=c_val_1, contour_value_2=c_val_2, tol=c_tol)

    # 3. Build properties panel content
    properties_content = Vector{Any}()
    
    if "euler_char" in properties_list
        # Compute Euler characteristic at the Fermi energy
        E_fermi = ef_val
        
        # Create mask for all occupied states (E <= E_F) across all bands
        chi_total = 0
        pockets_total = 0
        holes_total = 0
        
        for b in 1:size(data.plaquette_energies, 1)
            M = @view(data.plaquette_energies[b, :, :]) .<= E_fermi
            chi_b = compute_euler_characteristic(M; periodic=true)
            pockets_b, holes_b = compute_pockets_and_holes(M; periodic=true)
            chi_total += chi_b
            pockets_total += pockets_b
            holes_total += holes_b
        end
        
        push!(properties_content, html_div(style=Dict("backgroundColor" => "#ffffff", "padding" => "12px", "borderRadius" => "4px", "border" => "1px solid #ddd")) do
            html_h4("Euler Characteristic (χ)", style=Dict("marginTop" => "0", "marginBottom" => "10px", "color" => "#333")),
            html_p("χ = " * string(chi_total), style=Dict("fontSize" => "16px", "fontWeight" => "bold", "margin" => "8px 0")),
            html_p("N_pockets = " * string(pockets_total), style=Dict("fontSize" => "14px", "margin" => "6px 0")),
            html_p("N_holes = " * string(holes_total), style=Dict("fontSize" => "14px", "margin" => "6px 0")),
            html_hr(style=Dict("margin" => "8px 0", "borderColor" => "#ddd")),
            html_p("Relation: χ = N_pockets - N_holes", style=Dict("fontSize" => "12px", "fontStyle" => "italic", "color" => "#666", "margin" => "8px 0"))
        end)
    end
    
    if "chern_accum" in properties_list
        # Compute accumulated Chern number at the Fermi energy
        E_fermi = ef_val
        
        # Compute per-band accumulated Chern
        nbands = size(data.berry_curvature, 1)
        chern_accum_bands = zeros(Float64, nbands)
        
        for b in 1:nbands
            # Sum berry curvature for all k-points where E(k) <= E_F
            for ix in 1:size(data.plaquette_energies, 2), iy in 1:size(data.plaquette_energies, 3)
                if data.plaquette_energies[b, ix, iy] <= E_fermi
                    chern_accum_bands[b] += data.berry_curvature[b, ix, iy] / (2π)
                end
            end
        end
        
        # Total accumulated Chern across all bands
        chern_accum_total = sum(chern_accum_bands)
        
        push!(properties_content, html_div(style=Dict("backgroundColor" => "#ffffff", "padding" => "12px", "borderRadius" => "4px", "border" => "1px solid #ddd")) do
            html_h4("Accumulated Chern C(E_F)", style=Dict("marginTop" => "0", "marginBottom" => "10px", "color" => "#333")),
            html_p("C(E_F) = " * @sprintf("%.4f", chern_accum_total), style=Dict("fontSize" => "16px", "fontWeight" => "bold", "margin" => "8px 0", "color" => "#1f77b4")),
            html_hr(style=Dict("margin" => "8px 0", "borderColor" => "#ddd")),
            html_p("Per-band contributions:", style=Dict("fontSize" => "12px", "fontWeight" => "bold", "margin" => "8px 0", "color" => "#666")),
            [html_p("Band " * string(b) * ": " * @sprintf("%.4f", chern_accum_bands[b]), style=Dict("fontSize" => "12px", "margin" => "4px 0 4px 12px")) for b in 1:nbands]...
        end)
    end
    
    if "energy_crit_chern" in properties_list
        # Compute band-resolved energy corresponding to ± critical Chern values
        # Using 2D contour extraction for robustness to non-monotonic C(E) curves
        nbands = size(data.berry_curvature, 1)
        
        # Tolerance for contour matching (should match contour plotting)
        contour_tol = c_tol
        
        # For each band, compute the energies at ± critical Chern values from 2D data
        band_crit_energies = []
        
        for b in 1:nbands
            # Get full 2D arrays for this band (NOT flattened or sorted)
            E_band_2d = @view(data.plaquette_energies[b, :, :])
            C_band_2d = @view(data.cum_chern_per_band[b, :, :])
            
            # Find energies for both ± critical Chern 1 from 2D contour
            E_crit_1_pos = get_energy_at_chern_2d(E_band_2d, C_band_2d, c_val_1; tol=contour_tol)
            E_crit_1_neg = get_energy_at_chern_2d(E_band_2d, C_band_2d, -c_val_1; tol=contour_tol)
            
            # Find energies for both ± critical Chern 2 from 2D contour
            E_crit_2_pos = get_energy_at_chern_2d(E_band_2d, C_band_2d, c_val_2; tol=contour_tol)
            E_crit_2_neg = get_energy_at_chern_2d(E_band_2d, C_band_2d, -c_val_2; tol=contour_tol)
            
            push!(band_crit_energies, (E_crit_1_pos, E_crit_1_neg, E_crit_2_pos, E_crit_2_neg))
        end
        
        # Build display content
        band_contents = []
        for b in 1:nbands
            E_1_pos, E_1_neg, E_2_pos, E_2_neg = band_crit_energies[b]
            
            # Build the band's content dynamically based on which values exist
            band_items = []
            
            # Header
            push!(band_items, html_p("Band " * string(b) * ":", style=Dict("fontSize" => "12px", "fontWeight" => "bold", "margin" => "0 0 6px 0", "color" => "#333")))
            
            # Critical Chern 1
            push!(band_items, html_p("Crit 1:", style=Dict("fontSize" => "11px", "margin" => "4px 0 2px 12px", "color" => "#666", "fontWeight" => "bold")))
            
            if !isnan(E_1_pos)
                push!(band_items, html_p("C = +" * @sprintf("%.4f", c_val_1) * ": E = " * @sprintf("%.4f", E_1_pos), style=Dict("fontSize" => "10px", "margin" => "2px 0 0 24px", "fontWeight" => "bold", "color" => "#1f77b4")))
            else
                push!(band_items, html_p("C = +" * @sprintf("%.4f", c_val_1) * ": N/A", style=Dict("fontSize" => "10px", "margin" => "2px 0 0 24px", "color" => "#ccc")))
            end
            
            if !isnan(E_1_neg)
                push!(band_items, html_p("C = -" * @sprintf("%.4f", c_val_1) * ": E = " * @sprintf("%.4f", E_1_neg), style=Dict("fontSize" => "10px", "margin" => "2px 0 0 24px", "fontWeight" => "bold", "color" => "#ff7f0e")))
            else
                push!(band_items, html_p("C = -" * @sprintf("%.4f", c_val_1) * ": N/A", style=Dict("fontSize" => "10px", "margin" => "2px 0 0 24px", "color" => "#ccc")))
            end
            
            # Critical Chern 2
            push!(band_items, html_p("Crit 2:", style=Dict("fontSize" => "11px", "margin" => "4px 0 2px 12px", "color" => "#666", "fontWeight" => "bold")))
            
            if !isnan(E_2_pos)
                push!(band_items, html_p("C = +" * @sprintf("%.4f", c_val_2) * ": E = " * @sprintf("%.4f", E_2_pos), style=Dict("fontSize" => "10px", "margin" => "2px 0 0 24px", "fontWeight" => "bold", "color" => "#d62728")))
            else
                push!(band_items, html_p("C = +" * @sprintf("%.4f", c_val_2) * ": N/A", style=Dict("fontSize" => "10px", "margin" => "2px 0 0 24px", "color" => "#ccc")))
            end
            
            if !isnan(E_2_neg)
                push!(band_items, html_p("C = -" * @sprintf("%.4f", c_val_2) * ": E = " * @sprintf("%.4f", E_2_neg), style=Dict("fontSize" => "10px", "margin" => "2px 0 0 24px", "fontWeight" => "bold", "color" => "#2ca02c")))
            else
                push!(band_items, html_p("C = -" * @sprintf("%.4f", c_val_2) * ": N/A", style=Dict("fontSize" => "10px", "margin" => "2px 0 0 24px", "color" => "#ccc")))
            end
            
            band_div = html_div(band_items, style=Dict("marginBottom" => "12px", "paddingBottom" => "12px", "borderBottom" => "1px solid #eee"))
            
            push!(band_contents, band_div)
        end
        
        push!(properties_content, html_div([
            html_h4("Energy of Critical Chern (per-band ±)", style=Dict("marginTop" => "0", "marginBottom" => "12px", "color" => "#333")),
            band_contents...
        ], style=Dict("backgroundColor" => "#ffffff", "padding" => "12px", "borderRadius" => "4px", "border" => "1px solid #ddd")))
    end
    
    if isempty(properties_content)
        push!(properties_content, html_p("Select properties above to display", style=Dict("color" => "#999", "fontStyle" => "italic")))
    end

    # 4. Convert to base64 images for HTML rendering
    return plot_to_b64(p1), plot_to_b64(p2), plot_to_b64(p3), html_div(properties_content)
end

println("Starting Interactive Dashboard... Open http://127.0.0.1:8050 in your browser.")
println("To broadcast run this command in new termal: ")
println("npx localtunnel --port 8050 --subdomain my-qwz-dashboard")
run_server(app, "0.0.0.0", 8050)
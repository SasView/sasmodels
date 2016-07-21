CONVERSION_TABLE = {
    "adsorbed_layer": [
        "Core2ndMomentModel",
        {
            "scale": "scale",
            "second_moment": "second_moment",
            "density_shell": "density_poly",
            "sld_solvent": "sld_solv",
            "radius": "radius_core",
            "volfraction": "volf_cores",
            "background": "background",
            "adsorbed_amount": "ads_amount",
            "sld_shell": "sld_poly"
        }
    ],
    "barbell": [
        "BarBellModel",
        {
            "sld": "sld_barbell",
            "length": "len_bar",
            "bell_radius": "rad_bell",
            "radius": "rad_bar",
            "sld_solvent": "sld_solv"
        }
    ],
    "bcc_paracrystal": [
        "BCCrystalModel",
        {
            "sld": "sldSph",
            "sld_solvent": "sldSolv"
        }
    ],
    "be_polyelectrolyte": [
        "BEPolyelectrolyte",
        {
            "ionization_degree": "alpha",
            "polymer_concentration": "c",
            "salt_concentration": "cs",
            "virial_param": "h",
            "background": "background",
            "contrast_factor": "k",
            "bjerrum_length": "lb",
            "monomer_length": "b"
        }
    ],
    "binary_hard_sphere": [
        "BinaryHSModel",
        {
            "sld_sm": "ss_sld",
            "sld_lg": "ls_sld",
            "volfraction_sm": "vol_frac_ss",
            "radius_lg": "l_radius",
            "radius_sm": "s_radius",
            "volfraction_lg": "vol_frac_ls",
            "sld_solvent": "solvent_sld"
        }
    ],
    "broad_peak": [
        "BroadPeakModel",
        {
            "peak_pos": "q_peak",
            "scale": None,
            "lorentz_length": "length_l",
            "porod_scale": "scale_p",
            "lorentz_exp": "exponent_l",
            "lorentz_scale": "scale_l",
            "porod_exp": "exponent_p"
        }
    ],
    "capped_cylinder": [
        "CappedCylinderModel",
        {
            "sld": "sld_capcyl",
            "length": "len_cyl",
            "cap_radius": "rad_cap",
            "radius": "rad_cyl",
            "sld_solvent": "sld_solv"
        }
    ],
    "core_multi_shell": [
        "CoreMultiShellModel",
        {
            "thickness": "thick_shell",
            "sld": "sld_shell",
            "radius": "rad_core0",
            "sld_core": "sld_core0",
            "sld_solvent": "sld_solv",
            "n": "n_shells"
        }
    ],
    "core_shell_bicelle": [
        "CoreShellBicelleModel",
        {
            "phi": "axis_phi",
            "sld_core": "core_sld",
            "sld_rim": "rim_sld",
            "face_thickness": "face_thick",
            "sld_solvent": "solvent_sld",
            "rim_thickness": "rim_thick",
            "sld_face": "face_sld",
            "theta": "axis_theta"
        }
    ],
    "core_shell_cylinder": [
        "CoreShellCylinderModel",
        {
            "theta": "axis_theta",
            "phi": "axis_phi",
            "sld_shell": "shell_sld",
            "sld_solvent": "solvent_sld",
            "sld_core": "core_sld"
        }
    ],
    "core_shell_ellipsoid": [
        "CoreShellEllipsoidModel",
        {
            "phi": "axis_phi",
            "sld_core": "sld_core",
            "polar_shell": "polar_shell",
            "sld_solvent": "sld_solvent",
            "equat_shell": "equat_shell",
            "equat_core": "equat_core",
            "theta": "axis_theta",
            "polar_core": "polar_core",
            "sld_shell": "sld_shell"
        }
    ],
    "core_shell_ellipsoid_xt": [
        "CoreShellEllipsoidXTModel",
        {
            "phi": "axis_phi",
            "sld_core": "sld_core",
            "x_core": "X_core",
            "sld_solvent": "sld_solvent",
            "t_shell": "T_shell",
            "x_polar_shell": "XpolarShell",
            "theta": "axis_theta",
            "sld_shell": "sld_shell"
        }
    ],
    "core_shell_parallelepiped": [
        "CSParallelepipedModel",
        {
            "phi": "parallel_phi",
            "psi": "parallel_psi",
            "sld_core": "sld_pcore",
            "sld_c": "sld_rimC",
            "sld_b": "sld_rimB",
            "sld_solvent": "sld_solv",
            "a_side": "shortA",
            "sld_a": "sld_rimA",
            "b_side": "midB",
            "crim_thickness": "rimC",
            "theta": "parallel_theta",
            "arim_thickness": "rimA",
            "c_side": "longC",
            "brim_thickness": "rimB"
        }
    ],
    "core_shell_sphere": [
        "CoreShellModel",
        {
            "sld_shell": "shell_sld",
            "sld_solvent": "solvent_sld",
            "sld_core": "core_sld"
        }
    ],
    "correlation_length": [
        "CorrLengthModel",
        {
            "porod_scale": "scale_p",
            "lorentz_scale": "scale_l",
            "exponent_p": "exponent_p",
            "exponent_l": "exponent_l",
            "cor_length": "length_l"
        }
    ],
    "cylinder": [
        "CylinderModel",
        {
            "sld": "sldCyl",
            "theta": "cyl_theta",
            "phi": "cyl_phi",
            "sld_solvent": "sldSolv",
            "M0:sld": "M0_sld_cyl",
            "mtheta:sld": "M_theta_cyl",
            "mphi:sld": "M_phi_cyl",
            "M0:sld_solvent": "M0_sld_solv",
            "mtheta:sld_solvent": "M_theta_solv",
            "mphi:sld_solvent": "M_phi_solv",
            "up:frac_i": "Up_frac_i",
            "up:frac_f": "Up_frac_f",
            "up:angle": "Up_theta"
        }
    ],
    "dab": [
        "DABModel",
        {
            "length": "length"
        }
    ],
    "ellipsoid": [
        "EllipsoidModel",
        {
            "phi": "axis_phi",
            "r_equatorial": "radius_b",
            "sld": "sldEll",
            "theta": "axis_theta",
            "r_polar": "radius_a",
            "sld_solvent": "sldSolv"
        }
    ],
    "elliptical_cylinder": [
        "EllipticalCylinderModel",
        {
            "phi": "cyl_phi",
            "psi": "cyl_psi",
            "theta": "cyl_theta",
            "sld": "sldCyl",
            "axis_ratio": "r_ratio",
            "sld_solvent": "sldSolv"
        }
    ],
    "fcc_paracrystal": [
        "FCCrystalModel",
        {
            "sld": "sldSph",
            "sld_solvent": "sldSolv"
        }
    ],
    "flexible_cylinder": [
        "FlexibleCylinderModel",
        {
            "sld": "sldCyl",
            "sld_solvent": "sldSolv"
        }
    ],
    "flexible_cylinder_elliptical": [
        "FlexCylEllipXModel",
        {
            "sld": "sldCyl",
            "sld_solvent": "sldSolv"
        }
    ],
    "fractal": [
        "FractalModel",
        {
            "sld_block": "sldBlock",
            "radius": "radius",
            "cor_length": "cor_length",
            "sld_solvent": "sldSolv",
            "fractal_dim": "fractal_dim"
        }
    ],
    "fractal_core_shell": [
        "FractalCoreShellModel",
        {
            "sld_shell": "shell_sld",
            "sld_solvent": "solvent_sld",
            "sld_core": "core_sld"
        }
    ],
    "fuzzy_sphere": [
        "FuzzySphereModel",
        {
            "sld": "sldSph",
            "fuzziness": "fuzziness",
            "radius": "radius",
            "sld_solvent": "sldSolv"
        }
    ],
    "gauss_lorentz_gel": [
        "GaussLorentzGelModel",
        {
            "gauss_scale_factor": "scale_g",
            "dynamic_cor_length": "dyn_colength",
            "static_cor_length": "stat_colength",
            "background": "background",
            "lorentz_scale_factor": "scale_l"
        }
    ],
    "gaussian_peak": [
        "PeakGaussModel",
        {
            "sigma": "B"
        }
    ],
    "gel_fit": [
        "GelFitModel",
        {
            "gyration_radius": "radius",
            "lorentzian_scale": "lScale",
            "fractal_exp": "FractalExp",
            "cor_length": "zeta",
            "guinier_scale": "gScale"
        }
    ],
    "guinier": [
        "GuinierModel",
        {
            "rg": "rg"
        }
    ],
    "guinier_porod": [
        "GuinierPorodModel",
        {
            "s": "dim",
            "rg": "rg",
            "m": "m",
            "scale": "scale",
            "background": "background"
        }
    ],
    "hardsphere": [
        "HardsphereStructure",
        {
            "radius_effective_pd": "effect_radius_pd",
            "radius_effective": "effect_radius",
            "radius_effective_pd_n": "effect_radius_pd_n"
        }
    ],
    "hardsphere2": [
        "HardsphereStructure",
        {}
    ],
    "hardsphere3": [
        "HardsphereStructure",
        {}
    ],
    "hardsphere4": [
        "HardsphereStructure",
        {}
    ],
    "hayter_msa": [
        "HayterMSAStructure",
        {
            "salt_concentration": "saltconc",
            "radius_effective_pd": "effect_radius_pd",
            "radius_effective": "effect_radius",
            "radius_effective_pd_n": "effect_radius_pd_n"
        }
    ],
    "hollow_cylinder": [
        "HollowCylinderModel",
        {
            "phi": "axis_phi",
            "scale": "scale",
            "core_radius": "core_radius",
            "sld_solvent": "sldSolv",
            "length": "length",
            "radius": "radius",
            "background": "background",
            "sld": "sldCyl",
            "theta": "axis_theta"
        }
    ],
    "hollow_rectangular_prism": [
        "RectangularHollowPrismModel",
        {
            "b2a_ratio": "b2a_ratio",
            "a_side": "short_side",
            "sld": "sldPipe",
            "c_side": "c2a_ratio",
            "sld_solvent": "sldSolv",
            "thickness": "thickness"
        }
    ],
    "hollow_rectangular_prism_thin_walls": [
        "RectangularHollowPrismInfThinWallsModel",
        {
            "sld": "sldPipe",
            "b2a_ratio": "b2a_ratio",
            "a_side": "short_side",
            "c_side": "c2a_ratio",
            "sld_solvent": "sldSolv"
        }
    ],
    "lamellar": [
        "LamellarModel",
        {
            "sld": "sld_bi",
            "sld_solvent": "sld_sol",
            "thickness": "bi_thick"
        }
    ],
    "lamellar_hg": [
        "LamellarFFHGModel",
        {
            "sld": "sld_tail",
            "sld_solvent": "sld_solvent",
            "sld_head": "sld_head",
            "tail_length": "t_length",
            "head_length": "h_thickness"
        }
    ],
    "lamellar_hg_stack_caille": [
        "LamellarPSHGModel",
        {
            "Caille_parameter": "caille",
            "Nlayers": "n_plates",
            "sld_head": "sld_head",
            "tail_length": "deltaT",
            "head_length": "deltaH",
            "sld": "sld_tail",
            "sld_solvent": "sld_solvent"
        }
    ],
    "lamellar_stack_caille": [
        "LamellarPSModel",
        {
            "sld": "sld_bi",
            "Caille_parameter": "caille",
            "Nlayers": "N_plates",
            "sld_solvent": "sld_sol",
            "thickness": "delta"
        }
    ],
    "lamellar_stack_paracrystal": [
        "LamellarPCrystalModel",
        {
            "sld": "sld_layer",
            "spacing_polydisp": "pd_spacing",
            "sld_solvent": "sld_solvent"
        }
    ],
    "line": [
        "LineModel",
        {
            "slope": "B",
            "scale": None,
            "background": None,
            "intercept": "A"
        }
    ],
    "linear_pearls": [
        "LinearPearlsModel",
        {
            "sld": "sld_pearl",
            "sld_solvent": "sld_solv",
            "edge_sep": "edge_separation"
        }
    ],
    "lorentz": [
        "LorentzModel",
        {
            "cor_length": "length"
        }
    ],
    "mass_fractal": [
        "MassFractalModel",
        {
            "cutoff_length": "co_length",
            "radius": "radius",
            "mass_dim": "mass_dim"
        }
    ],
    "mass_surface_fractal": [
        "MassSurfaceFractal",
        {
            "cluster_rg": "cluster_rg",
            "mass_dim": "mass_dim",
            "radius": "radius",
            "surface_dim": "surface_dim",
            "primary_rg": "primary_rg"
        }
    ],
    "mono_gauss_coil": [
        "DebyeModel",
        {
            "radius_gyration": "rg",
            "scale": "scale",
            "background": "background"
        }
    ],
    "multilayer_vesicle": [
        "MultiShellModel",
        {
            "sld": "shell_sld",
            "thick_solvent": "w_thickness",
            "thick_shell": "s_thickness",
            "radius": "core_radius",
            "sld_solvent": "core_sld"
        }
    ],
    "onion": [
        "OnionExpShellModel",
        {
            "n_shells": "n_shells",
            "A": "A_shell",
            "sld_core": "sld_core0",
            "core_radius": "rad_core0",
            "sld_solvent": "sld_solv",
            "thickness": "thick_shell",
            "sld_in": "sld_in_shell",
            "sld_out": "sld_out_shell"
        }
    ],
    "parallelepiped": [
        "ParallelepipedModel",
        {
            "phi": "parallel_phi",
            "psi": "parallel_psi",
            "sld_solvent": "sldSolv",
            "a_side": "short_a",
            "b_side": "short_b",
            "sld": "sldPipe",
            "theta": "parallel_theta",
            "c_side": "long_c"
        }
    ],
    "peak_lorentz": [
        "PeakLorentzModel",
        {
            "peak_pos": "q0",
            "peak_hwhm": "B"
        }
    ],
    "pearl_necklace": [
        "PearlNecklaceModel",
        {
            "scale": "scale",
            "string_thickness": "thick_string",
            "sld_string": "sld_string",
            "sld_solvent": "sld_solv",
            "edge_separation": "edge_separation",
            "number_of_pearls": "num_pearls",
            "radius": "radius",
            "background": "background",
            "sld": "sld_pearl"
        }
    ],
    "poly_gauss_coil": [
        "Poly_GaussCoil",
        {
            "radius_gyration": "rg",
            "polydispersity": "poly_m",
            "scale": "scale",
            "background": "background"
        }
    ],
    "polymer_excl_volume": [
        "PolymerExclVolume",
        {
            "rg": "rg",
            "scale": "scale",
            "background": "background",
            "porod_exp": "m"
        }
    ],
    "polymer_micelle": [
        "MicelleSphCoreModel",
        {
            "sld_corona": "rho_corona",
            "sld_solvent": "rho_solv",
            "sld_core": "rho_core"
        }
    ],
    "porod": [
        "PorodModel",
        {
            "scale": "scale",
            "background": "background"
        }
    ],
    "power_law": [
        "PowerLawAbsModel",
        {
            "scale": "scale",
            "background": "background",
            "power": "m"
        }
    ],
    "pringle": [
        "PringlesModel",
        {
            "scale": "scale",
            "solvent_sld": "sld_solvent",
            "thickness": "thickness",
            "beta": "beta",
            "radius": "radius",
            "background": "background",
            "alpha": "alpha",
            "pringle_sld": "sld_pringle"
        }
    ],
    "raspberry": [
        "RaspBerryModel",
        {
            "volfraction_lg": "volf_Lsph",
            "volfraction_sm": "volf_Ssph",
            "radius_sm": "radius_Ssph",
            "radius_lg": "radius_Lsph",
            "sld_lg": "sld_Lsph",
            "sld_sm": "sld_Ssph",
            "sld_solvent": "sld_solv",
            "surface_fraction": "surfrac_Ssph",
            "penetration": "delta_Ssph"
        }
    ],
    "rectangular_prism": [
        "RectangularPrismModel",
        {
            "sld": "sldPipe",
            "b2a_ratio": "b2a_ratio",
            "a_side": "short_side",
            "c_side": "c2a_ratio",
            "sld_solvent": "sldSolv"
        }
    ],
    "rpa": [
        "RPA10Model",
        {
            "K12": "Kab", "K13": "Kac", "K14": "Kad",
            "K23": "Kbc", "K24": "Kbd", "K34": "Kcd",
            "N1": "Na", "N2": "Nb", "N3": "Nc", "N4": "Nd",
            "L1": "La", "L2": "Lb", "L3": "Lc", "L4": "Ld",
            "v1": "va", "v2": "vb", "v3": "vc", "v4": "vd",
            "b1": "ba", "b2": "bb", "b3": "bc", "b4": "bd",
            "Phi1": "Phia", "Phi2": "Phib", "Phi3": "Phic", "Phi4": "Phid",
            "case_num": "lcase_n"
        }
    ],
    "sc_paracrystal": [
        "SCCrystalModel",
        {
            "sld": "sldSph",
            "sld_solvent": "sldSolv"
        }
    ],
    "sphere": [
        "SphereModel",
        {
            "sld": "sldSph",
            "radius": "radius",
            "sld_solvent": "sldSolv",
            "M0:sld": "M0_sld_sph",
            "mtheta:sld": "M_theta_sph",
            "mphi:sld": "M_phi_sph",
            "M0:sld_solvent": "M0_sld_solv",
            "mtheta:sld_solvent": "M_theta_solv",
            "mphi:sld_solvent": "M_phi_solv",
            "up:frac_i": "Up_frac_i",
            "up:frac_f": "Up_frac_f",
            "up:angle": "Up_theta"
        }
    ],
    "magsphere": [
        "SphereModel",
        {
            "sld": "sldSph",
            "radius": "radius",
            "sld_solvent": "sldSolv",
            "sld_M0": "M0_sld_sph",
            "sphere_theta": "M_theta_sph",
            "sphere_phi": "M_phi_sph",
            "sld_solvent_M0": "M0_sld_solv",
            "solvent_theta": "M_theta_solv",
            "solvent_phi": "M_phi_solv",
            "in_spin": "Up_frac_i",
            "out_spin": "Up_frac_f",
            "spin_theta": "Up_theta"
        }
    ],
    "_spherepy": [
        "SphereModel",
        {
            "sld": "sldSph",
            "radius": "radius",
            "sld_solvent": "sldSolv"
        }
    ],
    "spherical_sld": [
        "SphericalSLDModel",
        {
            "radius_core": "rad_core0",
            "sld_core": "sld_core0",
            "sld_solvent": "sld_solv"
        }
    ],
    "squarewell": [
        "SquareWellStructure",
        {
            "radius_effective_pd": "effect_radius_pd",
            "radius_effective": "effect_radius",
            "radius_effective_pd_n": "effect_radius_pd_n"
        }
    ],
    "stacked_disks": [
        "StackedDisksModel",
        {
            "phi": "axis_phi",
            "sld_layer": "layer_sld",
            "sld_core": "core_sld",
            "theta": "axis_theta",
            "sld_solvent": "solvent_sld",
            "n_stacking": "n_stacking"
        }
    ],
    "star_polymer": [
        "StarPolymer",
        {
            "arms": "arms",
            "rg_squared": "R2"
        }
    ],
    "stickyhardsphere": [
        "StickyHSStructure",
        {
            "radius_effective_pd": "effect_radius_pd",
            "radius_effective": "effect_radius",
            "radius_effective_pd_n": "effect_radius_pd_n"
        }
    ],
    "surface_fractal": [
        "SurfaceFractalModel",
        {
            "cutoff_length": "co_length",
            "radius": "radius",
            "surface_dim": "surface_dim"
        }
    ],
    "teubner_strey": [
        "TeubnerStreyModel",
        {
            "a2": "scale"
        }
    ],
    "triaxial_ellipsoid": [
        "TriaxialEllipsoidModel",
        {
            "phi": "axis_phi",
            "req_minor": "semi_axisA",
            "rpolar": "semi_axisC",
            "req_major": "semi_axisB",
            "solvent_sld": "sldSolv",
            "psi": "axis_psi",
            "sld": "sldEll",
            "theta": "axis_theta"
        }
    ],
    "two_lorentzian": [
        "TwoLorentzianModel",
        {
            "lorentz_scale_1": "scale_1",
            "lorentz_scale_2": "scale_2",
            "lorentz_exp_1": "exponent_1",
            "lorentz_exp_2": "exponent_2",
            "lorentz_length_2": "length_2",
            "lorentz_length_1": "length_1",
            "background": "background"
        }
    ],
    "two_power_law": [
        "TwoPowerLawModel",
        {
            "coefficent_1": "coef_A",
            "power_2": "power2",
            "power_1": "power1",
            "background": "background",
            "crossover": "qc"
        }
    ],
    "unified_power_Rg": [
        "UnifiedPowerRgModel",
        {
        }
    ],
    "vesicle": [
        "VesicleModel",
        {
            "sld": "shell_sld",
            "sld_solvent": "solv_sld"
        }
    ]
}

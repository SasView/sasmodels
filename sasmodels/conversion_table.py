"""
Parameter conversion table

*CONVERSION_TABLE* gives the old model name and a dictionary of old parameter
names for each parameter in sasmodels.  This is used by :mod:`convert` to
determine the equivalent parameter set when comparing a sasmodels model to
the models defined in previous versions of SasView and sasmodels. This is now
versioned based on the version number of SasView.

When any sasmodels parameter or model name is changed, this must be modified to
account for that.

Usage::

    <old_Sasview_version> : {
        <new_model_name> : [
            <old_model_name> ,
            {
                <new_param_name_1> : <old_param_name_1>,
                ...
                <new_param_name_n> : <old_param_name_n>
            }
        ]
    }

Any future parameter and model name changes can and should be given in this
table for future compatibility.
"""

CONVERSION_TABLE = {
    (3, 1, 2) : {
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
                "radius_bell": "rad_bell",
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
                "radius_cap": "rad_cap",
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
                "n": "n_shells",
                "M0:sld_core": "M0_sld_core0",
                "mtheta:sld_core": "M_theta_core0",
                "mphi:sld_core": "M_phi_core0",
                "M0:sld1": "M0_sld_shell1",
                "mtheta:sld1": "M_theta_shell1",
                "mphi:sld1": "M_phi_shell1",
                "M0:sld2": "M0_sld_shell2",
                "mtheta:sld2": "M_theta_shell2",
                "mphi:sld2": "M_phi_shell2",
                "M0:sld3": "M0_sld_shell3",
                "mtheta:sld3": "M_theta_shell3",
                "mphi:sld3": "M_phi_shell3",
                "M0:sld4": "M0_sld_shell4",
                "mtheta:sld4": "M_theta_shell4",
                "mphi:sld4": "M_phi_shell4",
                "M0:sld_solvent": "M0_sld_solv",
                "mtheta:sld_solvent": "M_theta_solv",
                "mphi:sld_solvent": "M_phi_solv",
                "up:frac_i": "Up_frac_i",
                "up:frac_f": "Up_frac_f",
                "up:angle": "Up_theta",
            }
        ],
        "core_shell_bicelle": [
            "CoreShellBicelleModel",
            {
                "phi": "axis_phi",
                "sld_core": "core_sld",
                "sld_rim": "rim_sld",
                "thick_face": "face_thick",
                "sld_solvent": "solvent_sld",
                "thick_rim": "rim_thick",
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
        "core_shell_ellipsoid:1": [
            "CoreShellEllipsoidModel",
            {
                "sld_core": "sld_core",
                "sld_shell": "sld_shell",
                "sld_solvent": "sld_solvent",
                "radius_equat_core": "equat_core",
                "x_core": "polar_core",
                "thick_shell": "equat_shell",
                "x_polar_shell": "polar_shell",
                "theta": "axis_theta",
                "phi": "axis_phi",
            }
        ],
        "core_shell_ellipsoid": [
            "CoreShellEllipsoidXTModel",
            {
                "sld_core": "sld_core",
                "sld_shell": "sld_shell",
                "sld_solvent": "sld_solvent",
                "radius_equat_core": "equat_core",
                "thick_shell": "T_shell",
                "x_core": "X_core",
                "x_polar_shell": "XpolarShell",
                "theta": "axis_theta",
                "phi": "axis_phi",
            }
        ],
        "core_shell_parallelepiped": [
            "CSParallelepipedModel",
            {
                "sld_core": "sld_pcore",
                "sld_a": "sld_rimA",
                "sld_b": "sld_rimB",
                "sld_c": "sld_rimC",
                "sld_solvent": "sld_solv",
                "length_a": "shortA",
                "length_b": "midB",
                "length_c": "longC",
                "thick_rim_a": "rimA",
                "thick_rim_c": "rimC",
                "thick_rim_b": "rimB",
                "theta": "parallel_theta",
                "phi": "parallel_phi",
                "psi": "parallel_psi",
            }
        ],
        "core_shell_sphere": [
            "CoreShellModel",
            {
                "sld_core": "core_sld",
                "sld_shell": "shell_sld",
                "sld_solvent": "solvent_sld",
                "M0:sld_core": "M0_sld_core",
                "mtheta:sld_core": "M_theta_core",
                "mphi:sld_core": "M_phi_core",
                "M0:sld_shell": "M0_sld_shell",
                "mtheta:sld_shell": "M_theta_shell",
                "mphi:sld_shell": "M_phi_shell",
                "M0:sld_solvent": "M0_sld_solv",
                "mtheta:sld_solvent": "M_theta_solv",
                "mphi:sld_solvent": "M_phi_solv",
                "up:frac_i": "Up_frac_i",
                "up:frac_f": "Up_frac_f",
                "up:angle": "Up_theta"
            }
        ],
        "correlation_length": [
            "CorrLength",
            {
                "porod_scale": "scale_p",
                "lorentz_scale": "scale_l",
                "porod_exp": "exponent_p",
                "lorentz_exp": "exponent_l",
                "cor_length": "length_l"
            },
            "CorrLengthModel"
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
                "cor_length": "length"
            }
        ],
        "ellipsoid": [
            "EllipsoidModel",
            {
                "phi": "axis_phi",
                "radius_equatorial": "radius_b",
                "sld": "sldEll",
                "theta": "axis_theta",
                "radius_polar": "radius_a",
                "sld_solvent": "sldSolv"
            }
        ],
        "elliptical_cylinder": [
            "EllipticalCylinderModel",
            {
                "axis_ratio": "r_ratio",
                "radius_minor": "r_minor",
                "sld": "sldCyl",
                "sld_solvent": "sldSolv",
                "theta": "cyl_theta",
                "phi": "cyl_phi",
                "psi": "cyl_psi",
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
            "FractalCoreShell",
            {
                "sld_core": "core_sld",
                "sld_shell": "shell_sld",
                "sld_solvent": "solvent_sld",
                "radius": "radius",
                "thickness": "thickness",
                "fractal_dim": "frac_dim",
                "cor_length": "cor_length",
                "volfraction": "volfraction",
            },
            "FractalCoreShellModel"
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
            "GaussLorentzGel",
            {
                "gauss_scale": "scale_g",
                "cor_length_dynamic": "dyn_colength",
                "cor_length_static": "stat_colength",
                "background": "background",
                "lorentz_scale": "scale_l"
            },
            "GaussLorentzGelModel"
        ],
        "gaussian_peak": [
            "Peak Gauss Model",
            {
                "peak_pos": "q0",
                "sigma": "B",
            },
            "PeakGaussModel",
        ],
        "gel_fit": [
            "GelFitModel",
            {
                "rg": "radius",
                "lorentz_scale": "lScale",
                "guinier_scale": "gScale",
                "fractal_dim": "FractalExp",
                "cor_length": "zeta",
            }
        ],
        "guinier": [
            "Guinier",
            {
                "rg": "rg"
            },
            "GuinierModel",
        ],
        "guinier_porod": [
            "GuinierPorod",
            {
                "s": "dim",
                "rg": "rg",
                "porod_exp": "m",
                "scale": "scale",
                "background": "background"
            },
            "GuinierPorodModel",
        ],
        "hardsphere": [
            "HardsphereStructure",
            {
                "scale": "scale_factor",
                "radius_effective": "effect_radius",
            }
        ],
        "hayter_msa": [
            "HayterMSAStructure",
            {
                "scale": "scale_factor",
                "radius_effective": "effect_radius",
                "volfraction": "volfraction",
                "charge": "charge",
                "temperature": "temperature",
                "concentration_salt": "saltconc",
                "dielectconst": "dielectconst",
            }
        ],
        "hollow_cylinder": [
            "HollowCylinderModel",
            {
                "sld": "sldCyl",
                "sld_solvent": "sldSolv",
                "radius": "core_radius",
                "thickness": "radius",
                "length": "length",
                "theta": "axis_theta",
                "phi": "axis_phi",
            }
        ],
        "hollow_rectangular_prism": [
            "RectangularHollowPrismModel",
            {
                "sld": "sldPipe",
                "sld_solvent": "sldSolv",
                "length_a": "short_side",
                "b2a_ratio": "b2a_ratio",
                "c2a_ratio": "c2a_ratio",
                "thickness": "thickness",
            }
        ],
        "hollow_rectangular_prism_thin_walls": [
            "RectangularHollowPrismInfThinWallsModel",
            {
                "sld": "sldPipe",
                "sld_solvent": "sldSolv",
                "length_a": "short_side",
                "b2a_ratio": "b2a_ratio",
                "c2a_ratio": "c2a_ratio",
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
                "length_tail": "t_length",
                "length_head": "h_thickness"
            }
        ],
        "lamellar_hg_stack_caille": [
            "LamellarPSHGModel",
            {
                "sld": "sld_tail",
                "sld_head": "sld_head",
                "sld_solvent": "sld_solvent",
                "length_tail": "deltaT",
                "length_head": "deltaH",
                "d_spacing": "spacing",
                "Caille_parameter": "caille",
                "Nlayers": "n_plates",
            }
        ],
        "lamellar_stack_caille": [
            "LamellarPSModel",
            {
                "sld": "sld_bi",
                "sld_solvent": "sld_sol",
                "thickness": "delta",
                "d_spacing": "spacing",
                "Caille_parameter": "caille",
                "Nlayers": "n_plates",
            }
        ],
        "lamellar_stack_paracrystal": [
            "LamellarPCrystalModel",
            {
                "sld": "sld_layer",
                "sld_solvent": "sld_solvent",
                "thickness": "thickness",
                "d_spacing": "spacing",
                "sigma_d": "pd_spacing",
                "Nlayers": "Nlayers",
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
            "Lorentz",
            {
                "cor_length": "length"
            },
            "LorentzModel",
        ],
        "mass_fractal": [
            "MassFractalModel",
            {
                "cutoff_length": "co_length",
                "radius": "radius",
                "fractal_dim_mass": "mass_dim"
            }
        ],
        "mass_surface_fractal": [
            "MassSurfaceFractal",
            {
                "rg_cluster": "cluster_rg",
                "fractal_dim_mass": "mass_dim",
                "radius": "radius",
                "fractal_dim_surf": "surface_dim",
                "rg_primary": "primary_rg"
            }
        ],
        "mono_gauss_coil": [
            "Debye",
            {
                "rg": "rg",
                "i_zero": "scale",
                "background": "background",
            },
            "DebyeModel",
        ],
        "multilayer_vesicle": [
            "MultiShellModel",
            {
                "radius": "core_radius",
                "sld_solvent": "core_sld",
                "n_shells": "n_pairs",
                "thick_shell": "s_thickness",
                "sld": "shell_sld",
                "thick_solvent": "w_thickness",
            }
        ],
        "onion": [
            "OnionExpShellModel",
            {
                "n_shells": "n_shells",
                "A": "A_shell",
                "sld_core": "sld_core0",
                "radius_core": "rad_core0",
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
                "length_a": "short_a",
                "length_b": "short_b",
                "sld": "sldPipe",
                "theta": "parallel_theta",
                "length_c": "long_c",
                "M0:sld": "M0_sld_pipe",
                "mtheta:sld": "M_theta_pipe",
                "mphi:sld": "M_phi_pipe",
                "M0:sld_solvent": "M0_sld_solv",
                "mtheta:sld_solvent": "M_theta_solv",
                "mphi:sld_solvent": "M_phi_solv",
                "up:frac_i": "Up_frac_i",
                "up:frac_f": "Up_frac_f",
                "up:angle": "Up_theta",
            }
        ],
        "peak_lorentz": [
            "Peak Lorentz Model",
            {
                "peak_pos": "q0",
                "peak_hwhm": "B"
            },
            "PeakLorentzModel",
        ],
        "pearl_necklace": [
            "PearlNecklaceModel",
            {
                "scale": "scale",
                "thick_string": "thick_string",
                "sld_string": "sld_string",
                "sld_solvent": "sld_solv",
                "edge_sep": "edge_separation",
                "num_pearls": "num_pearls",
                "radius": "radius",
                "background": "background",
                "sld": "sld_pearl"
            }
        ],
        "poly_gauss_coil": [
            "Poly_GaussCoil",
            {
                "rg": "rg",
                "polydispersity": "poly_m",
                "i_zero": "scale",
                "background": "background",
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
                "sld_core": "rho_core",
                "ndensity": "ndensity",
                "v_core": "v_core",
                "v_corona": "v_corona",
                "radius_core": "radius_core",
                "rg": "radius_gyr",
                "d_penetration": "d_penetration",
                "n_aggreg": "n_aggreg",
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
                "sld_solvent": "sld_solvent",
                "thickness": "thickness",
                "beta": "beta",
                "radius": "radius",
                "background": "background",
                "alpha": "alpha",
                "sld": "sld_pringle"
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
                "length_a": "short_side",
                "b2a_ratio": "b2a_ratio",
                "c2a_ratio": "c2a_ratio",
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
        "spherical_sld": [
            "SphericalSLDModel",
            # Be lazy and use a generator expression to define
            #    sld1: sld_flat0, ...
            #    thickness1: thick_flat0, ...
            #    interface1: thick_inter0, ...
            #    shape1: func_inter0, ...
            #    nu1: nu_inter0, ...
            # but override thickness1 => rad_cor0 and sld1 => sld_core0.
            # Note: explicit key,value pairs given by **{...} override the
            # keys from the gnerator expression ((k,v) for k,v in seq) when
            # used as dict((generator), **{...})
            dict(((field_new+str(index+1), field_old+str(index))
                  for field_new, field_old in [("sld", "sld_flat"),
                                               ("thickness", "thick_flat"),
                                               ("interface", "thick_inter"),
                                               ("shape", "func_inter"),
                                               ("nu", "nu_inter"),]
                  for index in range(11)),
                 **{
                     "n_shells": "n_shells",
                     "n_steps": "npts_inter",
                     "sld_solvent": "sld_solv",
                     "thickness1": "rad_core0",
                     "sld1": "sld_core0",
                 })
        ],
        "squarewell": [
            "SquareWellStructure",
            {
                "scale": "scale_factor",
                "radius_effective": "effect_radius",
                "wellwidth": "wellwidth",
                "welldepth": "welldepth",
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
                "n_stacking": "n_stacking",
                "thick_layer": "layer_thick",
                "thick_core": "core_thick",
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
                "scale": "scale_factor",
                "radius_effective": "effect_radius",
            }
        ],
        "surface_fractal": [
            "SurfaceFractalModel",
            {
                "cutoff_length": "co_length",
                "radius": "radius",
                "fractal_dim_surf": "surface_dim"
            }
        ],
        "teubner_strey": [
            "TeubnerStreyModel",
            {
                # Note: parameters are completely rewritten in convert.py
                "volfraction_a": "volfraction_a",
                "sld_a": "sld_a",
                "sld_b": "sld_b",
                "d": "d",
                "xi": "xi",
            }
        ],
        "triaxial_ellipsoid": [
            "TriaxialEllipsoidModel",
            {
                "phi": "axis_phi",
                "radius_equat_minor": "semi_axisA",
                "radius_polar": "semi_axisC",
                "radius_equat_major": "semi_axisB",
                "sld_solvent": "sldSolv",
                "psi": "axis_psi",
                "sld": "sldEll",
                "theta": "axis_theta"
            }
        ],
        "two_lorentzian": [
            "TwoLorentzian",
            {
                "lorentz_scale_1": "scale_1",
                "lorentz_scale_2": "scale_2",
                "lorentz_exp_1": "exponent_1",
                "lorentz_exp_2": "exponent_2",
                "lorentz_length_2": "length_2",
                "lorentz_length_1": "length_1",
                "background": "background"
            },
            "TwoLorentzianModel",
        ],
        "two_power_law": [
            "TwoPowerLaw",
            {
                "coefficent_1": "coef_A",
                "power_2": "power2",
                "power_1": "power1",
                "background": "background",
                "crossover": "qc"
            },
            "TwoPowerLawModel",
        ],
        "unified_power_Rg": [
            "UnifiedPowerRg",
            dict(((field_new+str(index), field_old+str(index))
                  for field_new, field_old in [("rg", "Rg"),
                                               ("power", "power"),
                                               ("G", "G"),
                                               ("B", "B"),]
                  for index in range(11)),
                 **{
                     "background": "background",
                     "scale": "scale",
                 }),
            "UnifiedPowerRgModel",
        ],
        "vesicle": [
            "VesicleModel",
            {
                "sld": "shell_sld",
                "sld_solvent": "solv_sld"
            }
        ]
    }
}

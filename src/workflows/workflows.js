module.exports = {
    workflowData: {
        subworkflows: {
            espresso: {
                band_gap: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Band Gap",
                    units: [
                        {
                            config: { execName: "pw.x", flavorName: "pw_scf", name: "pw_scf" },
                            functions: { head: true },
                            type: "executionBuilder",
                        },
                        {
                            config: { execName: "pw.x", flavorName: "pw_nscf", name: "pw_nscf" },
                            type: "executionBuilder",
                        },
                    ],
                },
                band_structure: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Band Structure",
                    units: [
                        {
                            config: { execName: "pw.x", flavorName: "pw_scf", name: "pw_scf" },
                            functions: { head: true },
                            type: "executionBuilder",
                        },
                        {
                            config: { execName: "pw.x", flavorName: "pw_bands", name: "pw_bands" },
                            type: "executionBuilder",
                        },
                        {
                            config: { execName: "bands.x", flavorName: "bands", name: "bands" },
                            type: "executionBuilder",
                        },
                    ],
                },
                band_structure_average_esp: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Band Structure + averaged ESP",
                    units: [
                        {
                            config: {
                                name: "Set Material Index",
                                operand: "MATERIAL_INDEX",
                                value: 0,
                            },
                            type: "assignment",
                        },
                        {
                            config: { execName: "pw.x", flavorName: "pw_scf", name: "pw_scf" },
                            type: "executionBuilder",
                        },
                        {
                            attributes: { results: [{ name: "band_gaps" }] },
                            config: {
                                execName: "pw.x",
                                flavorName: "pw_bands",
                                flowchartId: "pw-bands-calculate-band-gap",
                                name: "pw_bands",
                            },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                input: [
                                    { name: "band_gaps", scope: "pw-bands-calculate-band-gap" },
                                ],
                                name: "Select indirect band gap",
                                operand: "BAND_GAP_INDIRECT",
                                value: '[bandgap for bandgap in band_gaps["values"] if bandgap["type"] == "indirect"][0]',
                            },
                            type: "assignment",
                        },
                        {
                            config: {
                                name: "Set Valence Band Maximum",
                                operand: "VBM",
                                value: 'BAND_GAP_INDIRECT["eigenvalueValence"]',
                            },
                            type: "assignment",
                        },
                        {
                            config: { execName: "bands.x", flavorName: "bands", name: "bands" },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "pp.x",
                                flavorName: "pp_electrostatic_potential",
                                name: "Electrostatic Potential (ESP)",
                            },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "average.x",
                                flavorName: "average",
                                flowchartId: "average-electrostatic-potential",
                                name: "average ESP",
                            },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                input: [
                                    {
                                        name: "averaged_potential_profile",
                                        scope: "average-electrostatic-potential",
                                    },
                                ],
                                name: "Set Macroscopically Averaged ESP Data",
                                operand: "AVG_DATA",
                                value: 'averaged_potential_profile["yDataSeries"][1]',
                            },
                            type: "assignment",
                        },
                    ],
                },
                band_structure_dos: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Band Structure + Density of States",
                    units: [
                        {
                            config: { execName: "pw.x", flavorName: "pw_scf", name: "pw_scf" },
                            type: "executionBuilder",
                        },
                        {
                            config: { execName: "pw.x", flavorName: "pw_bands", name: "pw_bands" },
                            type: "executionBuilder",
                        },
                        {
                            config: { execName: "bands.x", flavorName: "bands", name: "bands" },
                            type: "executionBuilder",
                        },
                        {
                            config: { execName: "pw.x", flavorName: "pw_nscf", name: "pw_nscf" },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "projwfc.x",
                                flavorName: "projwfc",
                                name: "projwfc",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                calc_valence_band_offset: {
                    application: { name: "python", version: "3.8.6" },
                    method: { name: "UnknownMethod" },
                    model: { name: "UnknownModel" },
                    name: "Calculate VBO",
                    units: [
                        {
                            config: {
                                name: "Difference of valence band maxima",
                                operand: "VBM_DIFF",
                                value: "VBM_LEFT - VBM_RIGHT",
                            },
                            type: "assignment",
                        },
                        {
                            config: {
                                name: "Difference of macroscopically averaged ESP in bulk",
                                operand: "AVG_ESP_DIFF",
                                value: "AVG_ESP_LEFT - AVG_ESP_RIGHT",
                            },
                            type: "assignment",
                        },
                        {
                            config: {
                                name: "Lineup of macroscopically averaged ESP in interface",
                                operand: "ESP_LINEUP",
                                value: "np.abs(AVG_ESP[0] - AVG_ESP[1])",
                            },
                            type: "assignment",
                        },
                        {
                            config: {
                                name: "Valence Band Offset",
                                operand: "VALENCE_BAND_OFFSET",
                                results: [{ name: "valence_band_offset" }],
                                value: "VBM_DIFF - AVG_ESP_DIFF + (np.sign(AVG_ESP_DIFF) * ESP_LINEUP)",
                            },
                            type: "assignment",
                        },
                    ],
                },
                dos: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Density of States",
                    units: [
                        {
                            config: { execName: "pw.x", flavorName: "pw_scf", name: "pw_scf" },
                            type: "executionBuilder",
                        },
                        {
                            config: { execName: "pw.x", flavorName: "pw_nscf", name: "pw_nscf" },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "projwfc.x",
                                flavorName: "projwfc",
                                name: "projwfc",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                electronic_density_mesh: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Electronic Density Mesh",
                    units: [
                        {
                            config: { execName: "pw.x", flavorName: "pw_scf", name: "pw_scf" },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "pp.x",
                                flavorName: "pp_density",
                                name: "pp_density",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                esm: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Effective Screening Medium (ESM)",
                    units: [
                        {
                            config: { execName: "pw.x", flavorName: "pw_esm", name: "pw_esm" },
                            type: "executionBuilder",
                        },
                    ],
                },
                esm_relax: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Effective Screening Medium (ESM) Relax",
                    units: [
                        {
                            config: {
                                execName: "pw.x",
                                flavorName: "pw_esm_relax",
                                name: "pw_esm_relax",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                espresso_xml_get_qpt_irr: {
                    application: { name: "python", version: "2.7.5" },
                    dynamicSubworkflow: { name: "getQpointIrrep", subfolder: "espresso" },
                    method: { name: "UnknownMethod" },
                    model: { name: "UnknownModel" },
                    name: "espresso-xml-get-qpt-irr",
                    units: [],
                },
                fixed_cell_relaxation: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Fixed-cell Relaxation",
                    units: [
                        {
                            config: { execName: "pw.x", flavorName: "pw_relax", name: "pw_relax" },
                            functions: { head: true },
                            type: "executionBuilder",
                        },
                    ],
                },
                gw_band_structure_band_gap_full_frequency: {
                    application: { name: "espresso", version: "6.3" },
                    method: { name: "PseudopotentialMethod", setSearchText: ".*dojo-oncv.*" },
                    model: { name: "DFTModel" },
                    name: "Full Frequency GW Band Structure + Band Gap",
                    units: [
                        {
                            config: { execName: "pw.x", flavorName: "pw_scf", name: "pw_scf" },
                            functions: { head: true },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "gw.x",
                                flavorName: "gw_bands_full_frequency",
                                name: "gw_bands_full_frequency",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                gw_band_structure_band_gap_plasmon_pole: {
                    application: { name: "espresso", version: "6.3" },
                    method: { name: "PseudopotentialMethod", setSearchText: ".*dojo-oncv.*" },
                    model: { name: "DFTModel" },
                    name: "Plasmon-Pole GW Band Structure + Band Gap",
                    units: [
                        {
                            config: { execName: "pw.x", flavorName: "pw_scf", name: "pw_scf" },
                            functions: { head: true },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "gw.x",
                                flavorName: "gw_bands_plasmon_pole",
                                name: "gw_bands_plasmon_pole",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                neb: {
                    application: { name: "espresso", version: "5.4.0" },
                    config: { isMultiMaterial: true },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Nudged Elastic Band (NEB)",
                    units: [
                        {
                            config: { execName: "neb.x", flavorName: "neb", name: "neb" },
                            type: "executionBuilder",
                        },
                    ],
                },
                ph_init_qpoints: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "ph-init-qpoints",
                    units: [
                        {
                            config: {
                                execName: "ph.x",
                                flavorName: "ph_init_qpoints",
                                name: "ph_init_qpoints",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                ph_single_irr_qpt: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "ph-single-irr-qpt",
                    units: [
                        {
                            config: {
                                execName: "ph.x",
                                flavorName: "ph_single_irr_qpt",
                                name: "ph_single_irr_qpt",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                phonon_dispersions: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Phonon Dispersions",
                    units: [
                        {
                            config: { execName: "pw.x", flavorName: "pw_scf", name: "pw_scf" },
                            type: "executionBuilder",
                        },
                        {
                            config: { execName: "ph.x", flavorName: "ph_grid", name: "ph_grid" },
                            type: "executionBuilder",
                        },
                        {
                            config: { execName: "q2r.x", flavorName: "q2r", name: "q2r" },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "matdyn.x",
                                flavorName: "matdyn_path",
                                name: "matdyn_path",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                phonon_dos: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Phonon Density of States",
                    units: [
                        {
                            config: { execName: "pw.x", flavorName: "pw_scf", name: "pw_scf" },
                            type: "executionBuilder",
                        },
                        {
                            config: { execName: "ph.x", flavorName: "ph_grid", name: "ph_grid" },
                            type: "executionBuilder",
                        },
                        {
                            config: { execName: "q2r.x", flavorName: "q2r", name: "q2r" },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "matdyn.x",
                                flavorName: "matdyn_grid",
                                name: "matdyn_grid",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                phonon_dos_dispersion: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Phonon Density of States + Dispersions",
                    units: [
                        {
                            config: { execName: "pw.x", flavorName: "pw_scf", name: "pw_scf" },
                            type: "executionBuilder",
                        },
                        {
                            config: { execName: "ph.x", flavorName: "ph_grid", name: "ph_grid" },
                            type: "executionBuilder",
                        },
                        {
                            config: { execName: "q2r.x", flavorName: "q2r", name: "q2r" },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "matdyn.x",
                                flavorName: "matdyn_grid",
                                name: "matdyn_grid",
                            },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "matdyn.x",
                                flavorName: "matdyn_path",
                                name: "matdyn_path",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                phonon_reduce: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "reduce",
                    units: [
                        {
                            config: {
                                execName: "ph.x",
                                flavorName: "ph_grid_restart",
                                name: "ph_grid_restart",
                            },
                            type: "executionBuilder",
                        },
                        {
                            config: { execName: "q2r.x", flavorName: "q2r", name: "q2r" },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "matdyn.x",
                                flavorName: "matdyn_grid",
                                name: "matdyn_grid",
                            },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "matdyn.x",
                                flavorName: "matdyn_path",
                                name: "matdyn_path",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                post_processor: {
                    application: { name: "shell", version: "4.2.46" },
                    method: { name: "UnknownMethod" },
                    model: { name: "UnknownModel" },
                    name: "post-processor",
                    units: [
                        {
                            config: {
                                execName: "sh",
                                flavorName: "espresso_collect_dynmat",
                                name: "shell",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                pre_processor: {
                    application: { name: "shell", version: "4.2.46" },
                    method: { name: "UnknownMethod" },
                    model: { name: "UnknownModel" },
                    name: "pre-processor",
                    units: [
                        {
                            config: {
                                execName: "sh",
                                flavorName: "espresso_link_outdir_save",
                                name: "shell",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                processing_find_minima: {
                    application: { name: "python", version: "3.8.6" },
                    method: { name: "UnknownMethod" },
                    model: { name: "UnknownModel" },
                    name: "Find ESP Value",
                    units: [
                        {
                            config: {
                                execName: "python",
                                flavorName: "processing:find_extrema",
                                flowchartId: "python-find-extrema",
                                name: "Find Extrema",
                            },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                input: [{ name: "STDOUT", scope: "python-find-extrema" }],
                                name: "Set Averaged ESP Value",
                                operand: "AVG_ESP",
                                value: 'json.loads(STDOUT)["minima"]',
                            },
                            type: "assignment",
                        },
                    ],
                },
                pw_scf: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "pw-scf",
                    units: [
                        {
                            config: { execName: "pw.x", flavorName: "pw_scf", name: "pw_scf" },
                            functions: { head: true },
                            type: "executionBuilder",
                        },
                    ],
                },
                recalculate_bands: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Recalculate Bands",
                    units: [
                        {
                            config: { execName: "pw.x", flavorName: "pw_bands", name: "pw_bands" },
                            type: "executionBuilder",
                        },
                        {
                            config: { execName: "bands.x", flavorName: "bands", name: "bands" },
                            type: "executionBuilder",
                        },
                    ],
                },
                surface_energy: {
                    application: { name: "espresso", version: "5.4.0" },
                    dynamicSubworkflow: { name: "surfaceEnergy" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Surface Energy",
                    units: [
                        {
                            config: { execName: "pw.x", flavorName: "pw_scf", name: "pw_scf" },
                            type: "executionBuilder",
                        },
                    ],
                },
                total_energy: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Total Energy",
                    units: [
                        {
                            config: { execName: "pw.x", flavorName: "pw_scf", name: "pw_scf" },
                            functions: { head: true },
                            type: "executionBuilder",
                        },
                    ],
                },
                variable_cell_relaxation: {
                    application: { name: "espresso", version: "5.4.0" },
                    config: { systemName: "espresso-variable-cell-relaxation" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Variable-cell Relaxation",
                    units: [
                        {
                            config: {
                                execName: "pw.x",
                                flavorName: "pw_vc-relax",
                                name: "pw_vc-relax",
                            },
                            functions: { head: true },
                            type: "executionBuilder",
                        },
                    ],
                },
                zero_point_energy: {
                    application: { name: "espresso", version: "5.4.0" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Zero Point Energy",
                    units: [
                        {
                            config: { execName: "pw.x", flavorName: "pw_scf", name: "pw_scf" },
                            type: "executionBuilder",
                        },
                        {
                            config: { execName: "ph.x", flavorName: "ph_gamma", name: "ph_zpe" },
                            type: "executionBuilder",
                        },
                    ],
                },
            },
            exabyteml: {
                classification_tail: {
                    application: { name: "python", version: "3.8.6" },
                    method: { name: "UnknownMethod" },
                    model: { name: "UnknownModel" },
                    name: "Machine Learning",
                    units: [
                        {
                            attributes: { enableRender: true },
                            config: {
                                execName: "python",
                                flavorName: "pyml:setup_variables_packages",
                                name: "Setup Variables and Packages",
                            },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "python",
                                flavorName: "pyml:data_input:read_csv:pandas",
                                name: "Data Input",
                            },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "python",
                                flavorName: "pyml:data_input:train_test_split:sklearn",
                                name: "Train Test Split",
                            },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "python",
                                flavorName: "pyml:pre_processing:standardization:sklearn",
                                name: "Data Standardize",
                            },
                            type: "executionBuilder",
                        },
                        {
                            attributes: {
                                results: [{ name: "workflow:pyml_predict" }],
                                tags: [
                                    "remove-all-results",
                                    "creates-predictions-csv-during-predict-phase",
                                ],
                            },
                            config: {
                                execName: "python",
                                flavorName: "pyml:model:random_forest_classification:sklearn",
                                name: "Model Train and Predict",
                            },
                            type: "executionBuilder",
                        },
                        {
                            attributes: {
                                postProcessors: [{ name: "remove_virtual_environment" }],
                                results: [
                                    {
                                        basename: "my_roc_plot.png",
                                        filetype: "image",
                                        name: "file_content",
                                    },
                                ],
                                tags: ["remove-all-results"],
                            },
                            config: {
                                execName: "python",
                                flavorName: "pyml:post_processing:roc_curve:sklearn",
                                name: "ROC Curve Plot",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                clustering_tail: {
                    application: { name: "python", version: "3.8.6" },
                    method: { name: "UnknownMethod" },
                    model: { name: "UnknownModel" },
                    name: "Machine Learning",
                    units: [
                        {
                            attributes: { enableRender: true },
                            config: {
                                execName: "python",
                                flavorName: "pyml:setup_variables_packages",
                                name: "Setup Variables and Packages",
                            },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "python",
                                flavorName: "pyml:data_input:read_csv:pandas",
                                name: "Data Input",
                            },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "python",
                                flavorName: "pyml:data_input:train_test_split:sklearn",
                                name: "Train Test Split",
                            },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "python",
                                flavorName: "pyml:pre_processing:standardization:sklearn",
                                name: "Data Standardize",
                            },
                            type: "executionBuilder",
                        },
                        {
                            attributes: {
                                results: [{ name: "workflow:pyml_predict" }],
                                tags: [
                                    "remove-all-results",
                                    "creates-predictions-csv-during-predict-phase",
                                ],
                            },
                            config: {
                                execName: "python",
                                flavorName: "pyml:model:k_means_clustering:sklearn",
                                name: "Model Train and Predict",
                            },
                            type: "executionBuilder",
                        },
                        {
                            attributes: {
                                postProcessors: [{ name: "remove_virtual_environment" }],
                                results: [
                                    {
                                        basename: "train_test_split.png",
                                        filetype: "image",
                                        name: "file_content",
                                    },
                                    {
                                        basename: "train_clusters.png",
                                        filetype: "image",
                                        name: "file_content",
                                    },
                                    {
                                        basename: "test_clusters.png",
                                        filetype: "image",
                                        name: "file_content",
                                    },
                                ],
                                tags: ["remove-all-results"],
                            },
                            config: {
                                execName: "python",
                                flavorName: "pyml:post_processing:pca_2d_clusters:matplotlib",
                                name: "2D PCA Clusters Plot",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                regression_tail: {
                    application: { name: "python", version: "3.8.6" },
                    method: { name: "UnknownMethod" },
                    model: { name: "UnknownModel" },
                    name: "Machine Learning",
                    units: [
                        {
                            attributes: { enableRender: true },
                            config: {
                                execName: "python",
                                flavorName: "pyml:setup_variables_packages",
                                name: "Setup Variables and Packages",
                            },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "python",
                                flavorName: "pyml:data_input:read_csv:pandas",
                                name: "Data Input",
                            },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "python",
                                flavorName: "pyml:data_input:train_test_split:sklearn",
                                name: "Train Test Split",
                            },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "python",
                                flavorName: "pyml:pre_processing:standardization:sklearn",
                                name: "Data Standardize",
                            },
                            type: "executionBuilder",
                        },
                        {
                            attributes: {
                                results: [{ name: "workflow:pyml_predict" }],
                                tags: [
                                    "remove-all-results",
                                    "creates-predictions-csv-during-predict-phase",
                                ],
                            },
                            config: {
                                execName: "python",
                                flavorName: "pyml:model:multilayer_perceptron:sklearn",
                                name: "Model Train and Predict",
                            },
                            type: "executionBuilder",
                        },
                        {
                            attributes: {
                                postProcessors: [{ name: "remove_virtual_environment" }],
                                results: [
                                    {
                                        basename: "my_parity_plot.png",
                                        filetype: "image",
                                        name: "file_content",
                                    },
                                ],
                                tags: ["remove-all-results"],
                            },
                            config: {
                                execName: "python",
                                flavorName: "pyml:post_processing:parity_plot:matplotlib",
                                name: "Parity Plot",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                train_head: {
                    application: { name: "python", version: "3.8.6" },
                    method: { name: "UnknownMethod" },
                    model: { name: "UnknownModel" },
                    name: "Set Up the Job",
                    units: [
                        {
                            config: {
                                flowchartId: "head-set-predict-status",
                                name: "Set Workflow Mode",
                                operand: "IS_WORKFLOW_RUNNING_TO_PREDICT",
                                tags: ["pyml:workflow-type-setter"],
                                value: "False",
                            },
                            type: "assignment",
                        },
                        {
                            config: {
                                enableRender: true,
                                flowchartId: "head-fetch-training-data",
                                input: [
                                    {
                                        basename: "{{DATASET_BASENAME}}",
                                        objectData: {
                                            CONTAINER: "",
                                            NAME: "{{DATASET_FILEPATH}}",
                                            PROVIDER: "",
                                            REGION: "",
                                        },
                                    },
                                ],
                                name: "Fetch Dataset",
                                source: "object_storage",
                            },
                            type: "io",
                        },
                        {
                            config: {
                                else: "end-of-ml-train-head",
                                flowchartId: "head-branch-on-predict-status",
                                input: [
                                    { name: "IS_WORKFLOW_RUNNING_TO_PREDICT", scope: "global" },
                                ],
                                name: "Train or Predict?",
                                statement: "IS_WORKFLOW_RUNNING_TO_PREDICT",
                                then: "head-fetch-trained-model",
                            },
                            type: "condition",
                        },
                        {
                            config: {
                                enableRender: true,
                                flowchartId: "head-fetch-trained-model",
                                input: [
                                    {
                                        basename: "",
                                        objectData: {
                                            CONTAINER: "",
                                            NAME: "",
                                            PROVIDER: "",
                                            REGION: "",
                                        },
                                    },
                                ],
                                name: "Fetch Trained Model as file",
                                source: "object_storage",
                                tags: ["set-io-unit-filenames"],
                            },
                            type: "io",
                        },
                        {
                            config: {
                                flowchartId: "end-of-ml-train-head",
                                name: "End Setup",
                                operand: "IS_SETUP_COMPLETE",
                                value: "True",
                            },
                            type: "assignment",
                        },
                    ],
                },
                train_kernel_ridge: {
                    application: { name: "exabyteml", version: "0.2.0" },
                    config: { isDraft: true, isMultiMaterial: true },
                    method: {
                        config: { data: {}, subtype: "least_squares", type: "kernel_ridge" },
                        name: "Method",
                    },
                    model: { config: { subtype: "re", type: "ml" }, name: "Model" },
                    name: "ML: Kernel Ridge Regression Train Model",
                    units: [
                        {
                            config: {
                                flowchartId: "io",
                                head: true,
                                input: [
                                    {
                                        endpoint: "data-frame",
                                        endpoint_options: {
                                            features: [
                                                "elemental_ratio:Si",
                                                "elemental_ratio:Ge",
                                                "ionization_potential:Ge",
                                                "ionization_potential:Si",
                                            ],
                                            ids: ["KuAsBRwofzGfHPWiT"],
                                            targets: ["band_gaps:indirect", "band_gaps:direct"],
                                        },
                                    },
                                ],
                                name: "input",
                                next: "data_transformation_manipulation",
                                source: "api",
                                status: "idle",
                                subtype: "dataFrame",
                                type: "io",
                            },
                            type: "io",
                        },
                        {
                            config: {
                                flowchartId: "data_transformation_manipulation",
                                inputData: {
                                    cleanMissingData: true,
                                    removeDuplicateRows: true,
                                    replaceNoneValuesWith: 0,
                                },
                                name: "clean data",
                                next: "data_transformation_scale_and_reduce",
                                operation: "data_transformation",
                                operationType: "manipulation",
                                status: "idle",
                                type: "processing",
                            },
                            type: "processing",
                        },
                        {
                            config: {
                                flowchartId: "data_transformation_scale_and_reduce",
                                inputData: { scaler: "standard_scaler" },
                                name: "scale and reduce",
                                next: "feature_selection_filter_based",
                                operation: "data_transformation",
                                operationType: "scale_and_reduce",
                                status: "idle",
                                type: "processing",
                            },
                            type: "processing",
                        },
                        {
                            config: {
                                flowchartId: "feature_selection_filter_based",
                                inputData: { algorithm: "f_regression", nFeatures: 0 },
                                name: "select features",
                                next: "train",
                                operation: "feature_selection",
                                operationType: "filter_based",
                                status: "idle",
                                type: "processing",
                            },
                            type: "processing",
                        },
                        {
                            attributes: { results: [{ name: "workflow:ml_predict" }] },
                            config: { execName: "train", flavorName: "train", name: "train" },
                            type: "executionBuilder",
                        },
                    ],
                },
                train_linear_least_squares: {
                    application: { name: "exabyteml", version: "0.2.0" },
                    config: { isDraft: true, isMultiMaterial: true },
                    method: {
                        config: { data: {}, subtype: "least_squares", type: "linear" },
                        name: "Method",
                    },
                    model: { config: { subtype: "re", type: "ml" }, name: "Model" },
                    name: "ML: Linear Least Squares Train Model",
                    units: [
                        {
                            config: {
                                flowchartId: "io",
                                head: true,
                                input: [
                                    {
                                        endpoint: "data-frame",
                                        endpoint_options: {
                                            features: [
                                                "elemental_ratio:Si",
                                                "elemental_ratio:Ge",
                                                "ionization_potential:Ge",
                                                "ionization_potential:Si",
                                            ],
                                            ids: ["KuAsBRwofzGfHPWiT"],
                                            targets: ["band_gaps:indirect", "band_gaps:direct"],
                                        },
                                    },
                                ],
                                name: "input",
                                next: "data_transformation_manipulation",
                                source: "api",
                                status: "idle",
                                subtype: "dataFrame",
                                type: "io",
                            },
                            type: "io",
                        },
                        {
                            config: {
                                flowchartId: "data_transformation_manipulation",
                                inputData: {
                                    cleanMissingData: true,
                                    removeDuplicateRows: true,
                                    replaceNoneValuesWith: 0,
                                },
                                name: "clean data",
                                next: "data_transformation_scale_and_reduce",
                                operation: "data_transformation",
                                operationType: "manipulation",
                                status: "idle",
                                type: "processing",
                            },
                            type: "processing",
                        },
                        {
                            config: {
                                flowchartId: "data_transformation_scale_and_reduce",
                                inputData: { scaler: "standard_scaler" },
                                name: "scale and reduce",
                                next: "feature_selection_filter_based",
                                operation: "data_transformation",
                                operationType: "scale_and_reduce",
                                status: "idle",
                                type: "processing",
                            },
                            type: "processing",
                        },
                        {
                            config: {
                                flowchartId: "feature_selection_filter_based",
                                inputData: { algorithm: "f_regression", nFeatures: 0 },
                                name: "select features",
                                next: "train",
                                operation: "feature_selection",
                                operationType: "filter_based",
                                status: "idle",
                                type: "processing",
                            },
                            type: "processing",
                        },
                        {
                            attributes: { results: [{ name: "workflow:ml_predict" }] },
                            config: { execName: "train", flavorName: "train", name: "train" },
                            type: "executionBuilder",
                        },
                    ],
                },
            },
            jupyterLab: {
                jupyter_notebook: {
                    application: { name: "jupyterLab", version: "3.0.3" },
                    method: { name: "UnknownMethod" },
                    model: { name: "UnknownModel" },
                    name: "Jupyter Notebook",
                    units: [
                        {
                            attributes: { preProcessors: [{ name: "record_python_environment" }] },
                            config: {
                                execName: "jupyter",
                                flavorName: "notebook",
                                name: "notebook",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
            },
            nwchem: {
                total_energy: {
                    application: { name: "nwchem", version: "7.0.2" },
                    method: { name: "LocalOrbitalMethod" },
                    model: { name: "DFTModel" },
                    name: "Total Energy",
                    units: [
                        {
                            config: {
                                execName: "nwchem",
                                flavorName: "nwchem_total_energy",
                                name: "nwchem_total_energy",
                            },
                            functions: { head: true },
                            type: "executionBuilder",
                        },
                    ],
                },
            },
            python: {
                python_script: {
                    application: { name: "python", version: "2.7.5" },
                    method: { name: "UnknownMethod" },
                    model: { name: "UnknownModel" },
                    name: "Python Script",
                    units: [
                        {
                            config: {
                                execName: "python",
                                flavorName: "hello_world",
                                name: "python",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
            },
            shell: {
                batch_espresso_pwscf: {
                    application: { name: "shell", version: "4.2.46" },
                    method: { name: "UnknownMethod" },
                    model: { name: "UnknownModel" },
                    name: "Shell Batch Job (Espresso PWSCF)",
                    units: [
                        {
                            config: {
                                execName: "sh",
                                flavorName: "job_espresso_pw_scf",
                                name: "shell",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                hello_world: {
                    application: { name: "shell", version: "4.2.46" },
                    method: { name: "UnknownMethod" },
                    model: { name: "UnknownModel" },
                    name: "Shell Hello World",
                    units: [
                        {
                            config: { execName: "sh", flavorName: "hello_world", name: "shell" },
                            type: "executionBuilder",
                        },
                    ],
                },
            },
            vasp: {
                band_gap: {
                    application: { name: "vasp", version: "5.3.5" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Band Gap",
                    units: [
                        {
                            config: { execName: "vasp", flavorName: "vasp", name: "vasp" },
                            functions: { head: true },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "vasp",
                                flavorName: "vasp_nscf",
                                name: "vasp_nscf",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                band_structure: {
                    application: { name: "vasp", version: "5.3.5" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Band Structure",
                    units: [
                        {
                            config: { execName: "vasp", flavorName: "vasp", name: "vasp" },
                            functions: { head: true },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "vasp",
                                flavorName: "vasp_bands",
                                name: "vasp_bands",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                band_structure_dos: {
                    application: { name: "vasp", version: "5.3.5" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Band Structure + Density of States",
                    units: [
                        {
                            config: { execName: "vasp", flavorName: "vasp", name: "vasp" },
                            functions: { addResults: ["density_of_states"], head: true },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "vasp",
                                flavorName: "vasp_bands",
                                name: "vasp_bands",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                dos: {
                    application: { name: "vasp", version: "5.3.5" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Density of States",
                    units: [
                        {
                            config: { execName: "vasp", flavorName: "vasp", name: "vasp" },
                            functions: { addResults: ["density_of_states"], head: true },
                            type: "executionBuilder",
                        },
                    ],
                },
                fixed_cell_relaxation: {
                    application: { name: "vasp", version: "5.3.5" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Fixed-cell Relaxation",
                    units: [
                        {
                            config: {
                                execName: "vasp",
                                flavorName: "vasp_relax",
                                name: "vasp_relax",
                            },
                            functions: { head: true },
                            type: "executionBuilder",
                        },
                    ],
                },
                initial_final_total_energies: {
                    application: { name: "vasp", version: "5.3.5" },
                    config: { functions: { setDefaultCompute: null }, isMultiMaterial: true },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Initial/Final Total Energies",
                    units: [
                        {
                            config: {
                                execName: "vasp",
                                flavorName: "vasp_neb_initial",
                                name: "vasp_neb_initial",
                            },
                            functions: { head: true },
                            type: "executionBuilder",
                        },
                        {
                            config: {
                                execName: "vasp",
                                flavorName: "vasp_neb_final",
                                name: "vasp_neb_final",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                neb_subworkflow: {
                    application: { name: "vasp", version: "5.3.5" },
                    config: { isMultiMaterial: true },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Nudged Elastic Band (NEB)",
                    units: [
                        {
                            config: { execName: "vasp", flavorName: "vasp_neb", name: "vasp_neb" },
                            type: "executionBuilder",
                        },
                    ],
                },
                prepare_images: {
                    application: { name: "shell", version: "4.2.46" },
                    config: { isMultiMaterial: true },
                    method: { name: "Method" },
                    model: { name: "Model" },
                    name: "Prepare Directories",
                    units: [
                        {
                            config: {
                                execName: "sh",
                                flavorName: "bash_vasp_prepare_neb_images",
                                name: "prepare-neb-images",
                            },
                            type: "executionBuilder",
                        },
                    ],
                },
                recalculate_bands: {
                    application: { name: "vasp", version: "5.3.5" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Recalculate Bands",
                    units: [
                        {
                            config: {
                                execName: "vasp",
                                flavorName: "vasp_bands",
                                name: "vasp_bands",
                            },
                            functions: { head: true },
                            type: "executionBuilder",
                        },
                    ],
                },
                surface_energy: {
                    application: { name: "vasp", version: "5.3.5" },
                    dynamicSubworkflow: { name: "surfaceEnergy" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Surface Energy",
                    units: [
                        {
                            config: { execName: "vasp", flavorName: "vasp", name: "vasp" },
                            functions: { head: true },
                            type: "executionBuilder",
                        },
                    ],
                },
                total_energy: {
                    application: { name: "vasp", version: "5.3.5" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Total Energy",
                    units: [
                        {
                            config: { execName: "vasp", flavorName: "vasp", name: "vasp" },
                            functions: { head: true },
                            type: "executionBuilder",
                        },
                    ],
                },
                variable_cell_relaxation: {
                    application: { name: "vasp", version: "5.3.5" },
                    config: { systemName: "vasp-variable-cell-relaxation" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Variable-cell Relaxation",
                    units: [
                        {
                            config: {
                                execName: "vasp",
                                flavorName: "vasp_vc_relax",
                                name: "vasp_vc_relax",
                            },
                            functions: { head: true },
                            type: "executionBuilder",
                        },
                    ],
                },
                zero_point_energy: {
                    application: { name: "vasp", version: "5.3.5" },
                    method: { name: "PseudopotentialMethod" },
                    model: { name: "DFTModel" },
                    name: "Zero Point Energy",
                    units: [
                        {
                            config: { execName: "vasp", flavorName: "vasp_zpe", name: "vasp_zpe" },
                            functions: { head: true },
                            type: "executionBuilder",
                        },
                    ],
                },
            },
        },
        workflows: {
            espresso: {
                band_gap: { name: "Band Gap", units: [{ name: "band_gap", type: "subworkflow" }] },
                band_structure: {
                    name: "Band Structure",
                    units: [{ name: "band_structure", type: "subworkflow" }],
                },
                band_structure_dos: {
                    name: "Band Structure + Density of States",
                    units: [{ name: "band_structure_dos", type: "subworkflow" }],
                },
                dos: { name: "Density of States", units: [{ name: "dos", type: "subworkflow" }] },
                electronic_density_mesh: {
                    name: "Electronic Density Mesh",
                    units: [{ name: "electronic_density_mesh", type: "subworkflow" }],
                },
                esm: {
                    name: "Effective Screening Medium (ESM)",
                    units: [{ name: "esm", type: "subworkflow" }],
                },
                esm_relax: {
                    name: "Effective Screening Medium (ESM) Relax",
                    units: [{ name: "esm_relax", type: "subworkflow" }],
                },
                fixed_cell_relaxation: {
                    name: "Fixed-cell Relaxation",
                    units: [{ name: "fixed_cell_relaxation", type: "subworkflow" }],
                },
                gw_band_structure_band_gap_full_frequency: {
                    name: "Full Frequency GW Band Structure + Band Gap",
                    units: [
                        { name: "gw_band_structure_band_gap_full_frequency", type: "subworkflow" },
                    ],
                },
                gw_band_structure_band_gap_plasmon_pole: {
                    name: "Plasmon-Pole GW Band Structure + Band Gap",
                    units: [
                        { name: "gw_band_structure_band_gap_plasmon_pole", type: "subworkflow" },
                    ],
                },
                neb: {
                    name: "Nudged Elastic Band (NEB)",
                    units: [{ name: "neb", type: "subworkflow" }],
                },
                phonon_dispersions: {
                    name: "Phonon Dispersions",
                    units: [{ name: "phonon_dispersions", type: "subworkflow" }],
                },
                phonon_dos: {
                    name: "Phonon Density of States",
                    units: [{ name: "phonon_dos", type: "subworkflow" }],
                },
                phonon_dos_dispersion: {
                    name: "Phonon Density of States + Dispersions",
                    units: [{ name: "phonon_dos_dispersion", type: "subworkflow" }],
                },
                phonon_map: {
                    name: "Phonon Map",
                    units: [
                        {
                            name: "phononMap",
                            type: "workflow",
                            units: [
                                { name: "pw_scf", type: "subworkflow" },
                                { name: "ph_init_qpoints", type: "subworkflow" },
                                { name: "espresso_xml_get_qpt_irr", type: "subworkflow" },
                            ],
                        },
                        {
                            config: {
                                functions: { setDefaultCompute: null },
                                input: { name: "Q_POINTS" },
                                mapUnit: true,
                            },
                            name: "phonon_map_workflow",
                            type: "workflow",
                            units: [
                                { name: "pre_processor", type: "subworkflow" },
                                { name: "ph_single_irr_qpt", type: "subworkflow" },
                                { name: "post_processor", type: "subworkflow" },
                            ],
                        },
                        { name: "phonon_reduce", type: "subworkflow" },
                    ],
                },
                recalculate_bands: {
                    name: "Recalculate Bands",
                    units: [{ name: "recalculate_bands", type: "subworkflow" }],
                },
                surface_energy: {
                    name: "Surface Energy",
                    units: [{ name: "surface_energy", type: "subworkflow" }],
                },
                total_energy: {
                    name: "Total Energy",
                    units: [{ name: "total_energy", type: "subworkflow" }],
                },
                valence_band_offset: {
                    name: "Valence Band Offset (2D)",
                    units: [
                        {
                            config: { name: "BS + Avg ESP (Interface)" },
                            name: "band_structure_average_esp",
                            type: "subworkflow",
                            unitConfigs: [
                                {
                                    config: { attributes: { operand: "INTERFACE", value: 0 } },
                                    index: 0,
                                    type: "assignment",
                                },
                            ],
                        },
                        {
                            config: { name: "Find ESP Values (Interface)" },
                            name: "processing_find_minima",
                            type: "subworkflow",
                            unitConfigs: [
                                {
                                    config: { attributes: { operand: "AVG_ESP_INTERFACE" } },
                                    index: 1,
                                    type: "assignment",
                                },
                            ],
                        },
                        {
                            config: { name: "BS + Avg ESP (interface left)" },
                            name: "band_structure_average_esp",
                            type: "subworkflow",
                            unitConfigs: [
                                {
                                    config: { attributes: { operand: "INTERFACE_LEFT", value: 1 } },
                                    index: 0,
                                    type: "assignment",
                                },
                                {
                                    config: { attributes: { operand: "VBM_LEFT" } },
                                    index: 3,
                                    type: "assignment",
                                },
                            ],
                        },
                        {
                            config: { name: "Find ESP Value (Interface left)" },
                            name: "processing_find_minima",
                            type: "subworkflow",
                            unitConfigs: [
                                {
                                    config: { attributes: { operand: "AVG_ESP_LEFT" } },
                                    index: 1,
                                    type: "assignment",
                                },
                            ],
                        },
                        {
                            config: { name: "BS + Avg ESP (interface right)" },
                            name: "band_structure_average_esp",
                            type: "subworkflow",
                            unitConfigs: [
                                {
                                    config: {
                                        attributes: { operand: "INTERFACE_RIGHT", value: 2 },
                                    },
                                    index: 0,
                                    type: "assignment",
                                },
                                {
                                    config: { attributes: { operand: "VBM_RIGHT" } },
                                    index: 3,
                                    type: "assignment",
                                },
                            ],
                        },
                        {
                            config: { name: "Find ESP Value (Interface right)" },
                            name: "processing_find_minima",
                            type: "subworkflow",
                            unitConfigs: [
                                {
                                    config: { attributes: { operand: "AVG_ESP_RIGHT" } },
                                    index: 1,
                                    type: "assignment",
                                },
                            ],
                        },
                        { name: "calc_valence_band_offset", type: "subworkflow" },
                    ],
                },
                variable_cell_relaxation: {
                    name: "Variable-cell Relaxation",
                    units: [{ name: "variable_cell_relaxation", type: "subworkflow" }],
                },
                zero_point_energy: {
                    name: "Zero Point Energy",
                    units: [{ name: "zero_point_energy", type: "subworkflow" }],
                },
            },
            exabyteml: {
                classification_workflow: {
                    config: { attributes: { isUsingDataset: true } },
                    name: "Python ML Train Classification",
                    units: [
                        { name: "train_head", type: "subworkflow" },
                        { name: "classification_tail", type: "subworkflow" },
                    ],
                },
                clustering_workflow: {
                    config: { attributes: { isUsingDataset: true } },
                    name: "Python ML Train Clustering",
                    units: [
                        { name: "train_head", type: "subworkflow" },
                        { name: "classification_tail", type: "subworkflow" },
                    ],
                },
                regression_workflow: {
                    config: { attributes: { isUsingDataset: true } },
                    name: "Python ML Train Regression",
                    units: [
                        { name: "train_head", type: "subworkflow" },
                        { name: "regression_tail", type: "subworkflow" },
                    ],
                },
                train_kernel_ridge: {
                    name: "ML: Kernel Ridge Regression Train Model",
                    units: [{ name: "train_kernel_ridge", type: "subworkflow" }],
                },
                train_linear_least_squares: {
                    name: "ML: Linear Least Squares Train Model",
                    units: [{ name: "train_linear_least_squares", type: "subworkflow" }],
                },
            },
            jupyterLab: {
                jupyter_notebook: {
                    name: "Jupyter Notebook",
                    units: [{ name: "jupyter_notebook", type: "subworkflow" }],
                },
            },
            nwchem: {
                total_energy: {
                    name: "Total Energy",
                    units: [{ name: "total_energy", type: "subworkflow" }],
                },
            },
            python: {
                python_script: {
                    name: "Python Script",
                    units: [{ name: "python_script", type: "subworkflow" }],
                },
            },
            shell: {
                batch_espresso_pwscf: {
                    name: "Shell Batch Job (Espresso PWSCF)",
                    units: [{ name: "batch_espresso_pwscf", type: "subworkflow" }],
                },
                hello_world: {
                    name: "Shell Hello World",
                    units: [{ name: "hello_world", type: "subworkflow" }],
                },
            },
            vasp: {
                band_gap: { name: "Band Gap", units: [{ name: "band_gap", type: "subworkflow" }] },
                band_structure: {
                    name: "Band Structure",
                    units: [{ name: "band_structure", type: "subworkflow" }],
                },
                band_structure_dos: {
                    name: "Band Structure + Density of States",
                    units: [{ name: "band_structure_dos", type: "subworkflow" }],
                },
                dos: { name: "Density of States", units: [{ name: "dos", type: "subworkflow" }] },
                fixed_cell_relaxation: {
                    name: "Fixed-cell Relaxation",
                    units: [{ name: "fixed_cell_relaxation", type: "subworkflow" }],
                },
                neb: {
                    name: "Nudged Elastic Band (NEB)",
                    units: [
                        { name: "initial_final_total_energies", type: "subworkflow" },
                        { name: "prepare_images", type: "subworkflow" },
                        { name: "neb_subworkflow", type: "subworkflow" },
                    ],
                },
                recalculate_bands: {
                    name: "Recalculate Bands",
                    units: [{ name: "recalculate_bands", type: "subworkflow" }],
                },
                surface_energy: {
                    name: "Surface Energy",
                    units: [{ name: "surface_energy", type: "subworkflow" }],
                },
                total_energy: {
                    name: "Total Energy",
                    units: [{ name: "total_energy", type: "subworkflow" }],
                },
                variable_cell_relaxation: {
                    name: "Variable-cell Relaxation",
                    units: [{ name: "variable_cell_relaxation", type: "subworkflow" }],
                },
                zero_point_energy: {
                    name: "Zero Point Energy",
                    units: [{ name: "zero_point_energy", type: "subworkflow" }],
                },
            },
        },
    },
};

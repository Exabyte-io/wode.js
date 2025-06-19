import context from "./context";

const {
    BoundaryConditionsFormDataProvider,
    MLSettingsContextProvider,
    MLTrainTestSplitContextProvider,
    NEBFormDataProvider,
    PlanewaveCutoffsContextProvider,
    PointsGridFormDataProvider,
    PointsPathFormDataProvider,
    ExplicitPointsPathFormDataProvider,
    ExplicitPointsPath2PIBAFormDataProvider,
    HubbardJContextProvider,
    HubbardUContextProvider,
    HubbardVContextProvider,
    HubbardContextProviderLegacy,
    IonDynamicsContextProvider,
    CollinearMagnetizationContextProvider,
    NonCollinearMagnetizationContextProvider,
    VASPContextProvider,
    VASPNEBContextProvider,
    QEPWXContextProvider,
    QENEBContextProvider,
    NWChemTotalEnergyContextProvider,
} = context;

const CONTEXT_DOMAINS = {
    important: "important", // used to generate `ImportantSettings` form
};

const _makeImportant = (config) => Object.assign(config, { domain: CONTEXT_DOMAINS.important });

/** ********************************
 * Method-based context providers *
 ********************************* */

export const wodeProviders = {
    // NOTE: subworkflow-level data manager. Will override the unit-level data with the same name via subworkflow context.
    PlanewaveCutoffDataManager: {
        providerCls: PlanewaveCutoffsContextProvider,
        config: _makeImportant({ name: "cutoffs", entityName: "subworkflow" }),
    },
    KGridFormDataManager: {
        providerCls: PointsGridFormDataProvider,
        config: _makeImportant({ name: "kgrid" }),
    },
    QGridFormDataManager: {
        providerCls: PointsGridFormDataProvider,
        config: _makeImportant({ name: "qgrid", divisor: 5 }), // Using less points for Qgrid by default
    },
    IGridFormDataManager: {
        providerCls: PointsGridFormDataProvider,
        config: _makeImportant({ name: "igrid", divisor: 0.2 }), // Using more points for interpolated grid by default
    },
    QPathFormDataManager: {
        providerCls: PointsPathFormDataProvider,
        config: _makeImportant({ name: "qpath" }),
    },
    IPathFormDataManager: {
        providerCls: PointsPathFormDataProvider,
        config: _makeImportant({ name: "ipath" }),
    },
    KPathFormDataManager: {
        providerCls: PointsPathFormDataProvider,
        config: _makeImportant({ name: "kpath" }),
    },
    ExplicitKPathFormDataManager: {
        providerCls: ExplicitPointsPathFormDataProvider,
        config: _makeImportant({ name: "explicitKPath" }),
    },
    ExplicitKPath2PIBAFormDataManager: {
        providerCls: ExplicitPointsPath2PIBAFormDataProvider,
        config: _makeImportant({ name: "explicitKPath2PIBA" }),
    },
    HubbardJContextManager: {
        providerCls: HubbardJContextProvider,
        config: _makeImportant({ name: "hubbard_j" }),
    },
    HubbardUContextManager: {
        providerCls: HubbardUContextProvider,
        config: _makeImportant({ name: "hubbard_u" }),
    },
    HubbardVContextManager: {
        providerCls: HubbardVContextProvider,
        config: _makeImportant({ name: "hubbard_v" }),
    },
    HubbardContextManagerLegacy: {
        providerCls: HubbardContextProviderLegacy,
        config: _makeImportant({ name: "hubbard_legacy" }),
    },
    // NEBFormDataManager context is stored under the same key (`input`) as InputDataManager contexts.
    NEBFormDataManager: {
        providerCls: NEBFormDataProvider,
        config: _makeImportant({ name: "neb" }),
    },
    BoundaryConditionsFormDataManager: {
        providerCls: BoundaryConditionsFormDataProvider,
        config: _makeImportant({ name: "boundaryConditions" }),
    },
    MLSettingsDataManager: {
        providerCls: MLSettingsContextProvider,
        config: _makeImportant({ name: "mlSettings" }),
    },
    MLTrainTestSplitDataManager: {
        providerCls: MLTrainTestSplitContextProvider,
        config: _makeImportant({ name: "mlTrainTestSplit" }),
    },
    IonDynamicsContextProvider: {
        providerCls: IonDynamicsContextProvider,
        config: _makeImportant({ name: "dynamics" }),
    },
    CollinearMagnetizationDataManager: {
        providerCls: CollinearMagnetizationContextProvider,
        config: _makeImportant({ name: "collinearMagnetization" }),
    },
    NonCollinearMagnetizationDataManager: {
        providerCls: NonCollinearMagnetizationContextProvider,
        config: _makeImportant({ name: "nonCollinearMagnetization" }),
    },
    QEPWXInputDataManager: {
        providerCls: QEPWXContextProvider,
        config: { name: "input" },
    },
    QENEBInputDataManager: {
        providerCls: QENEBContextProvider,
        config: { name: "input" },
    },
    VASPInputDataManager: {
        providerCls: VASPContextProvider,
        config: { name: "input" },
    },
    VASPNEBInputDataManager: {
        providerCls: VASPNEBContextProvider,
        config: { name: "input" },
    },
    NWChemInputDataManager: {
        providerCls: NWChemTotalEnergyContextProvider,
        config: { name: "input" },
    },
};

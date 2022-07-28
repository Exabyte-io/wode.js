import { context } from "@exabyte-io/mode.js";
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
} = context;

const CONTEXT_DOMAINS = {
    important: "important",  // used to generate `ImportantSettings` form
}

const _makeImportant = (config) => Object.assign(config, {domain: CONTEXT_DOMAINS.important});

/**********************************
 * Method-based context providers *
 **********************************/

// NOTE: subworkflow-level data manager. Will override the unit-level data with the same name via subworkflow context.
export const PlanewaveCutoffDataManager = PlanewaveCutoffsContextProvider.getConstructorConfig(_makeImportant({
    name: "cutoffs",
    entityName: "subworkflow"
}));
export const KGridFormDataManager = PointsGridFormDataProvider.getConstructorConfig(_makeImportant({name: "kgrid"}));
export const QGridFormDataManager = PointsGridFormDataProvider.getConstructorConfig(_makeImportant({
    name: "qgrid",
    divisor: 5,  // Using less points for Qgrid by default
}));
export const IGridFormDataManager = PointsGridFormDataProvider.getConstructorConfig(_makeImportant({
    name: "igrid",
    divisor: 0.2,  // Using more points for interpolated grid by default
}));

export const QPathFormDataManager = PointsPathFormDataProvider.getConstructorConfig(_makeImportant({name: "qpath"}));
export const IPathFormDataManager = PointsPathFormDataProvider.getConstructorConfig(_makeImportant({name: "ipath"}));
export const KPathFormDataManager = PointsPathFormDataProvider.getConstructorConfig(_makeImportant({name: "kpath"}));
export const ExplicitKPathFormDataManager = ExplicitPointsPathFormDataProvider.getConstructorConfig(_makeImportant({name: "explicitKPath"}));
export const ExplicitKPath2PIBAFormDataManager = ExplicitPointsPath2PIBAFormDataProvider.getConstructorConfig(_makeImportant({name: "explicitKPath2PIBA"}));

// NEBFormDataManager context is stored under the same key (`input`) as InputDataManager contexts.
export const NEBFormDataManager = NEBFormDataProvider.getConstructorConfig(_makeImportant({name: "neb"}));

export const BoundaryConditionsFormDataManager = BoundaryConditionsFormDataProvider.getConstructorConfig(_makeImportant({name: "boundaryConditions"}));

export const MLSettingsDataManager = MLSettingsContextProvider.getConstructorConfig(_makeImportant({name: "mlSettings"}));
export const MLTrainTestSplitDataManager = MLTrainTestSplitContextProvider.getConstructorConfig(_makeImportant({name: "mlTrainTestSplit"}));


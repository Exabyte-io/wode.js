import type { AssignmentUnitSchema } from "@mat3ra/esse/dist/js/types";

import type { ApplicationExternalContext } from "../mixins/ApplicationContextMixin";
import type { MaterialExternalContext } from "../mixins/MaterialContextMixin";
import type { MaterialsExternalContext } from "../mixins/MaterialsContextMixin";
import type { MaterialsSetExternalContext } from "../mixins/MaterialsSetContextMixin";
import type {
    JobExternalContext,
    UnitContext,
    WorkflowExternalContext,
} from "./base/ContextProvider";
import type { JinjaExternalContext } from "./base/JSONSchemaDataProvider";
import BoundaryConditionsFormDataManager from "./BoundaryConditionsFormDataManager";
import QENEBInputDataManager from "./by_application/espresso/QENEBInputDataManager";
import QEPWXInputDataManager, {
    type MethodDataExternalContext,
} from "./by_application/espresso/QEPWXInputDataManager";
import NWChemInputDataManager from "./by_application/nwchem/NWChemInputDataManager";
import VASPInputDataManager from "./by_application/vasp/VASPInputDataManager";
import VASPNEBInputDataManager from "./by_application/vasp/VASPNEBInputDataManager";
import CollinearMagnetizationDataManager from "./CollinearMagnetizationDataManager";
import HubbardContextManagerLegacy from "./Hubbard/HubbardContextManagerLegacy";
import HubbardJContextManager from "./Hubbard/HubbardJContextManager";
import HubbardUContextManager from "./Hubbard/HubbardUContextManager";
import HubbardVContextManager from "./Hubbard/HubbardVContextManager";
import IonDynamicsDataManager from "./IonDynamicsDataManager";
import MLSettingsDataManager from "./MLSettingsDataManager";
import MLTrainTestSplitDataManager from "./MLTrainTestSplitDataManager";
import NEBFormDataManager from "./NEBFormDataManager";
import NonCollinearMagnetizationDataManager from "./NonCollinearMagnetizationDataManager";
import PlanewaveCutoffDataManager from "./PlanewaveCutoffDataManager";
import IGridFormDataManager from "./PointsGrid/IGridFormDataManager";
import KGridFormDataManager from "./PointsGrid/KGridFormDataManager";
import QGridFormDataManager from "./PointsGrid/QGridFormDataManager";
import ExplicitKPath2PIBAFormDataManager from "./PointsPath/ExplicitKPath2PIBAFormDataManager";
import ExplicitKPathFormDataManager from "./PointsPath/ExplicitKPathFormDataManager";
import IPathFormDataManager from "./PointsPath/IPathFormDataManager";
import KPathFormDataManager from "./PointsPath/KPathFormDataManager";
import QPathFormDataManager from "./PointsPath/QPathFormDataManager";

/**
 * Registry mapping provider names (as they appear in templates) to their classes.
 * This is the single source of truth for provider mappings.
 */
export const PROVIDER_REGISTRY = {
    PlanewaveCutoffDataManager,
    KGridFormDataManager,
    QGridFormDataManager,
    IGridFormDataManager,
    QPathFormDataManager,
    IPathFormDataManager,
    KPathFormDataManager,
    ExplicitKPathFormDataManager,
    ExplicitKPath2PIBAFormDataManager,
    HubbardJContextManager,
    HubbardUContextManager,
    HubbardVContextManager,
    HubbardContextManagerLegacy,
    NEBFormDataManager,
    BoundaryConditionsFormDataManager,
    MLSettingsDataManager,
    MLTrainTestSplitDataManager,
    IonDynamicsContextProvider: IonDynamicsDataManager, // Note: name mismatch preserved from original
    CollinearMagnetizationDataManager,
    NonCollinearMagnetizationDataManager,
    QEPWXInputDataManager,
    QENEBInputDataManager,
    VASPInputDataManager,
    VASPNEBInputDataManager,
    NWChemInputDataManager,
} as const;

export type AssignmentContext = Record<
    AssignmentUnitSchema["operand"],
    AssignmentUnitSchema["value"]
>;

export type SubworkflowContext = {
    subworkflowContext: AssignmentContext;
};

export type ScopeGlobalExternalContext = {
    scopeGlobal?: Record<string, unknown>;
};

/**
 * External context type used by ExecutionUnitInput when creating providers.
 * This type is always expected to be present when providers are instantiated.
 */
export type ExternalContext = ApplicationExternalContext &
    WorkflowExternalContext &
    JobExternalContext &
    MaterialsExternalContext &
    MethodDataExternalContext &
    MaterialsSetExternalContext &
    MaterialExternalContext &
    JinjaExternalContext &
    SubworkflowContext &
    ScopeGlobalExternalContext;

/**
 * Type for provider names as they appear in templates.
 */
export type ProviderName = keyof typeof PROVIDER_REGISTRY;

/**
 * Union type of all context provider instances.
 * Derived from the registry for type safety.
 */
export type AnyContextProvider = InstanceType<(typeof PROVIDER_REGISTRY)[ProviderName]>;

/**
 * Factory function to create a context provider instance from its name.
 *
 * @param name - The provider name as it appears in templates
 * @param context - The unit context
 * @param externalContext - The external context (must match the ExternalContext type defined in this file)
 * @returns An instance of the requested context provider
 * @throws Error if the provider name is unknown
 */
export function createProvider(
    name: ProviderName,
    context: UnitContext,
    externalContext: ExternalContext,
): AnyContextProvider {
    const ProviderClass = PROVIDER_REGISTRY[name];
    if (!ProviderClass) {
        throw new Error(`Unknown provider: ${name}`);
    }
    // The full ExternalContext is a superset of each provider's expected context type,
    // so passing it to every provider's createFromUnitContext is type-safe (no assertion).
    return ProviderClass.createFromUnitContext(context, externalContext);
}

// Re-export all provider classes for convenience
export {
    BoundaryConditionsFormDataManager,
    QENEBInputDataManager,
    QEPWXInputDataManager,
    NWChemInputDataManager,
    VASPInputDataManager,
    VASPNEBInputDataManager,
    CollinearMagnetizationDataManager,
    HubbardContextManagerLegacy,
    HubbardJContextManager,
    HubbardUContextManager,
    HubbardVContextManager,
    IonDynamicsDataManager,
    MLSettingsDataManager,
    MLTrainTestSplitDataManager,
    NEBFormDataManager,
    NonCollinearMagnetizationDataManager,
    PlanewaveCutoffDataManager,
    IGridFormDataManager,
    KGridFormDataManager,
    QGridFormDataManager,
    ExplicitKPath2PIBAFormDataManager,
    ExplicitKPathFormDataManager,
    IPathFormDataManager,
    KPathFormDataManager,
    QPathFormDataManager,
};

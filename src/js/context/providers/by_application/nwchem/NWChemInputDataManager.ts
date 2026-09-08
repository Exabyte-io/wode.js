import { PERIODIC_TABLE } from "@exabyte-io/periodic-table.js";
import JSONSchemasInterface from "@mat3ra/esse/dist/js/esse/JSONSchemasInterface";
import type {
    InputContextItemSchema,
    NWChemTotalEnergyContextProviderSchema,
} from "@mat3ra/esse/dist/js/types";
import type { JSONSchema7 } from "json-schema";

import materialContextMixin, {
    type MaterialContextMixin,
    type MaterialExternalContext,
} from "../../../mixins/MaterialContextMixin";
import type {
    JobExternalContext,
    UnitContext,
    WorkflowExternalContext,
} from "../../base/ContextProvider";
import JSONSchemaDataProvider, {
    type JinjaExternalContext,
} from "../../base/JSONSchemaDataProvider";

type Data = NWChemTotalEnergyContextProviderSchema;
type Schema = InputContextItemSchema & { data: Data };
type ExternalContext = JinjaExternalContext &
    MaterialExternalContext &
    WorkflowExternalContext &
    JobExternalContext;

// eslint-disable-next-line @typescript-eslint/no-empty-interface
interface NWChemInputDataManager extends MaterialContextMixin {}

const jsonSchemaId =
    "context-providers-directory/by-application/nwchem-total-energy-context-provider";

class NWChemInputDataManager extends JSONSchemaDataProvider<Schema, ExternalContext> {
    readonly name = "input" as const;

    readonly domain = "executable" as const;

    readonly entityName = "unit" as const;

    isEdited = false;

    workflowHasRelaxation: boolean;

    static createFromUnitContext(unitContext: UnitContext, externalContext: ExternalContext) {
        const contextItem = this.findContextItem<Schema>(unitContext, "input");

        return new NWChemInputDataManager(contextItem, externalContext);
    }

    readonly contextProviderName = "nwchem-total-energy" as const;

    readonly jsonSchema: JSONSchema7;

    constructor(config: Partial<Schema>, externalContext: ExternalContext) {
        super(config, externalContext);
        this.initMaterialContextMixin(externalContext);

        this.workflowHasRelaxation = externalContext.workflowHasRelaxation;

        const jsonSchema = JSONSchemasInterface.getSchemaById(jsonSchemaId);

        if (!jsonSchema) {
            throw new Error("Failed to get JSON schema");
        }

        this.jsonSchema = jsonSchema;
    }

    /*
     * TODO: Create ability for user to define CHARGE, MULT, BASIS and FUNCTIONAL parameters.
     */
    getDefaultData() {
        const basis = this.material.getBasis();
        const { workflowHasRelaxation } = this;

        const NTYP = basis.uniqueElements.length;
        const ATOMIC_POSITIONS_WITHOUT_CONSTRAINTS = basis.atomicPositions.join("\n") || "";
        const ATOMIC_SPECIES = basis.uniqueElements
            .map((symbol) => `${symbol} ${PERIODIC_TABLE[symbol].atomic_mass} `)
            .join("\n");

        basis.toCartesian();
        const atomicPositions = basis.getAtomicPositionsWithConstraintsAsStrings();

        return {
            CHARGE: 0,
            MULT: 1,
            BASIS: "6-31G",
            NAT: atomicPositions.length,
            NTYP,
            ATOMIC_POSITIONS: atomicPositions.join("\n"),
            ATOMIC_POSITIONS_WITHOUT_CONSTRAINTS,
            ATOMIC_SPECIES,
            FUNCTIONAL: "B3LYP",
            CARTESIAN: basis.toCartesian !== undefined,
            // Continue from the previous unit's RTDB, which carries its optimized geometry. Units of
            // one job share a work_dir, so nwchem.db is already there. Deliberately NOT
            // `jobHasParent || ...` as espresso does: rupy's nwchem prepare_restart symlinks only the
            // parent's perm/, while nwchem.db sits at the work-dir root, so a child job would render
            // `restart` with no database to restart from.
            RESTART: Boolean(workflowHasRelaxation),
            contextProviderName: this.contextProviderName,
        };
    }
}

materialContextMixin(NWChemInputDataManager.prototype);

export default NWChemInputDataManager;

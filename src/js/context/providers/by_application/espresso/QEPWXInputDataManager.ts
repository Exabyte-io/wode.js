import { PERIODIC_TABLE } from "@exabyte-io/periodic-table.js";
import JSONSchemasInterface from "@mat3ra/esse/dist/js/esse/JSONSchemasInterface";
import type {
    BaseMethod,
    InputContextItemSchema,
    QEPwxContextProviderSchema,
} from "@mat3ra/esse/dist/js/types";
import type { AtomicElementValue } from "@mat3ra/made/dist/js/basis/elements";
import type { JSONSchema7 } from "json-schema";
import path from "path";

import materialContextMixin, {
    type MaterialContextMixin,
    type MaterialExternalContext,
    type OrderedMaterial,
} from "../../../mixins/MaterialContextMixin";
import materialsContextMixin, {
    type MaterialsContextMixin,
    type MaterialsExternalContext,
} from "../../../mixins/MaterialsContextMixin";
import type {
    JobExternalContext,
    UnitContext,
    WorkflowExternalContext,
} from "../../base/ContextProvider";
import JSONSchemaDataProvider, {
    type JinjaExternalContext,
} from "../../base/JSONSchemaDataProvider";

// TODO: create a task to define correct type for MethodData
type MethodData = BaseMethod["data"] & {
    pseudo?: { element: AtomicElementValue; filename?: string; path?: string }[];
};

export type MethodDataExternalContext = {
    methodData?: MethodData;
};

type Data = QEPwxContextProviderSchema;
type Schema = InputContextItemSchema & { data: Data };
type ExternalContext = JinjaExternalContext &
    WorkflowExternalContext &
    MaterialExternalContext &
    JobExternalContext &
    MethodDataExternalContext &
    MaterialsExternalContext;

interface QEPWXInputDataManager extends MaterialContextMixin, MaterialsContextMixin {}

const jsonSchemaId = "context-providers-directory/by-application/qe-pwx-context-provider";

class QEPWXInputDataManager extends JSONSchemaDataProvider<Schema, ExternalContext> {
    readonly name = "input" as const;

    readonly domain = "executable" as const;

    readonly entityName = "unit" as const;

    isEdited = false;

    methodData?: MethodData;

    jobHasParent: boolean;

    workflowHasRelaxation: boolean;

    static createFromUnitContext(unitContext: UnitContext, externalContext: ExternalContext) {
        const contextItem = this.findContextItem<Schema>(unitContext, "input");

        return new QEPWXInputDataManager(contextItem, externalContext);
    }

    readonly jsonSchema: JSONSchema7;

    constructor(config: Partial<Schema>, externalContext: ExternalContext) {
        super(config, externalContext);
        this.initMaterialsContextMixin(externalContext);
        this.initMaterialContextMixin(externalContext);

        this.methodData = externalContext.methodData || {};
        this.jobHasParent = externalContext.jobHasParent;
        this.workflowHasRelaxation = externalContext.workflowHasRelaxation;

        const jsonSchema = JSONSchemasInterface.getSchemaById(jsonSchemaId);

        if (!jsonSchema) {
            throw new Error("Failed to get JSON schema");
        }

        this.jsonSchema = jsonSchema;
    }

    private buildQEPWXContext(material: OrderedMaterial): Data {
        const basis = material.getBasis();
        const lattice = material.getLattice();
        const { jobHasParent, workflowHasRelaxation } = this;

        const ATOMIC_SPECIES = basis.uniqueElements.map((symbol) => {
            const pseudo = (this.methodData?.pseudo || []).find((p) => p.element === symbol);
            return {
                X: symbol,
                Mass_X: PERIODIC_TABLE[symbol].atomic_mass,
                PseudoPot_X: pseudo?.filename || path.basename(pseudo?.path || ""),
            };
        });

        const uniqueElementsWithLabels = [...new Set(basis.elementsWithLabelsArray)];

        const ATOMIC_SPECIES_WITH_LABELS = uniqueElementsWithLabels.map((symbol) => {
            const symbolWithoutLabel = symbol.replace(/\d$/, "") as AtomicElementValue;
            const label = symbol.match(/\d$/g) ? symbol.match(/\d$/g)?.[0] : "";
            const pseudo = (this.methodData?.pseudo || []).find(
                (p) => p.element === symbolWithoutLabel,
            );
            return {
                X: `${symbolWithoutLabel}${label}`,
                Mass_X: PERIODIC_TABLE[symbolWithoutLabel].atomic_mass,
                PseudoPot_X: pseudo?.filename || path.basename(pseudo?.path || ""),
            };
        });

        const CELL_PARAMETERS = {
            v1: lattice.vectorArrays[0],
            v2: lattice.vectorArrays[1],
            v3: lattice.vectorArrays[2],
        };

        const ATOMIC_POSITIONS = basis.elementsCoordinatesConstraintsArray.map(
            ([element, label, coordinate, constraint]) => {
                return {
                    X: `${element}${label}`,
                    x: coordinate[0],
                    y: coordinate[1],
                    z: coordinate[2],
                    "if_pos(1)": constraint[0] ? 1 : 0,
                    "if_pos(2)": constraint[1] ? 1 : 0,
                    "if_pos(3)": constraint[2] ? 1 : 0,
                };
            },
        );

        return {
            IBRAV: 0,
            RESTART_MODE: jobHasParent || workflowHasRelaxation ? "restart" : "from_scratch",
            ATOMIC_SPECIES,
            ATOMIC_SPECIES_WITH_LABELS,
            NAT: basis.atomicPositions.length,
            NTYP: basis.uniqueElements.length,
            NTYP_WITH_LABELS: uniqueElementsWithLabels.length,
            ATOMIC_POSITIONS,
            ATOMIC_POSITIONS_WITHOUT_CONSTRAINTS: basis.atomicPositions.join("\n"),
            CELL_PARAMETERS,
            contextProviderName: "qe-pwx" as const,
        };
    }

    getDefaultData() {
        // the below values are read from PlanewaveDataManager instead
        // ECUTWFC = 40;
        // ECUTRHO = 200;

        return {
            ...this.buildQEPWXContext(this.material),
            perMaterial: this.materials.map((material) => this.buildQEPWXContext(material)),
        };
    }
}

materialContextMixin(QEPWXInputDataManager.prototype);
materialsContextMixin(QEPWXInputDataManager.prototype);
// methodDataContextMixin(QEPWXInputDataManager.prototype);
// workflowContextMixin(QEPWXInputDataManager.prototype);
// jobContextMixin(QEPWXInputDataManager.prototype);

export default QEPWXInputDataManager;

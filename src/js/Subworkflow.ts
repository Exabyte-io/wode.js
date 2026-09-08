import { Application } from "@mat3ra/ade";
import { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import {
    type Defaultable,
    defaultableEntityMixin,
} from "@mat3ra/code/dist/js/entity/mixins/DefaultableMixin";
import {
    type HashedEntity,
    hashedEntityMixin,
} from "@mat3ra/code/dist/js/entity/mixins/HashedEntityMixin";
import {
    type NamedEntity,
    namedEntityMixin,
} from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";
import JSONSchemasInterface from "@mat3ra/esse/dist/js/esse/JSONSchemasInterface";
import type {
    AssignmentUnitSchema,
    BaseInMemoryEntitySchema,
    JobSchema,
    SubworkflowSchema,
} from "@mat3ra/esse/dist/js/types";
import { type ComputedEntityMixin, computedEntityMixin } from "@mat3ra/ide/dist/js/compute";
import type { Material } from "@mat3ra/made";
import { type PseudopotentialMethod, Model, ModelFactory } from "@mat3ra/mode";
import type { MetaPropertyHolder } from "@mat3ra/prode";
import { ApplicationRegistry, setUnitLinks } from "@mat3ra/standata";
import { Utils } from "@mat3ra/utils";

import type { AssignmentContext, ExternalContext } from "./context/providers";
import type { WorkflowExternalContext } from "./context/providers/base/ContextProvider";
import { createConvergenceParameter } from "./convergence/factory";
import { UnitTag, UnitType } from "./enums";
import {
    type SubworkflowSchemaMixin,
    subworkflowSchemaMixin,
} from "./generated/SubworkflowSchemaMixin";
import { AssignmentUnit, ConditionUnit, SubworkflowUnit, UnitFactory } from "./units";
import type { AnySubworkflowUnit } from "./units/factory";
import type { WorkflowRenderContext } from "./Workflow";

type ConvergenceConfig = {
    parameter: "N_k" | "N_k_nonuniform";
    parameterInitial: number | [number, number, number];
    parameterIncrement: number;
    result: string;
    resultInitial: number;
    condition: string;
    operator: string;
    tolerance: number;
    maxOccurrences: number;
    externalContext: SubworkflowExternalContext;
};

interface Subworkflow
    extends Defaultable,
        NamedEntity,
        SubworkflowSchemaMixin,
        HashedEntity,
        ComputedEntityMixin<SubworkflowSchemaMixin["compute"]> {}

type SubworkflowExternalContext = WorkflowRenderContext & WorkflowExternalContext;

function isAssignmentUnitSchema(
    unit: SubworkflowSchema["units"][number],
): unit is AssignmentUnitSchema {
    return unit.type === UnitType.assignment;
}

export type SubworkflowEntity = SubworkflowSchema & BaseInMemoryEntitySchema;

type Schema = SubworkflowEntity;

class Subworkflow<S extends Schema = Schema> extends InMemoryEntity<S> {
    private ModelFactory: typeof ModelFactory;

    private applicationInstance: Application;

    unitsInstances!: AnySubworkflowUnit[];

    modelInstance: Model;

    properties: string[] = [];

    repetition = 0;

    declare static createDefault: () => Subworkflow;

    static get jsonSchema() {
        return JSONSchemasInterface.getSchemaById("workflow/subworkflow");
    }

    // NoInfer: keep default S (or an explicit type arg) instead of inferring S from the config literal.
    constructor(config: NoInfer<S>, _ModelFactory = ModelFactory) {
        super(config);
        this.ModelFactory = _ModelFactory;

        this.applicationInstance = new Application(this.application);
        this.modelInstance = this.ModelFactory.create({
            ...this.model,
            application: this.application,
        });
        this.setUnits(this.units.map((cfg) => UnitFactory.createInSubworkflow(cfg)));
    }

    static get defaultConfig() {
        const defaultName = "New Subworkflow";
        const application = new ApplicationRegistry().getDefaultApplication();

        if (!application) {
            throw new Error("No default application found");
        }

        return {
            _id: Utils.uuid.getUUID(),
            name: defaultName,
            application,
            // `Model.defaultConfig` from @mat3ra/mode may omit `functional`; ESSE subworkflow schema requires it once schemas are registered.
            model: { ...Model.defaultConfig, functional: "pbe" as const },
            properties: [],
            units: [],
        };
    }

    setRepetition(repetition: number) {
        this.repetition = repetition;
        this.unitsInstances.forEach((u) => u.setRepetition(repetition));
    }

    getAsUnit() {
        return new SubworkflowUnit({
            type: UnitType.subworkflow,
            _id: this.id,
            name: this.name,
            preProcessors: [],
            postProcessors: [],
            monitors: [],
            results: [],
            // Concrete id for graph linking ("" is a valid string but not a usable node id).
            flowchartId: Utils.uuid.getUUID(),
        });
    }

    setApplication(application: Application) {
        // TODO: adjust the logic above to take into account whether units need re-rendering after version change etc.
        // reset units if application name changes
        const previousApplicationName = this.application.name;
        this.applicationInstance = application;

        if (previousApplicationName !== application.name) {
            // TODO: figure out how to set a default unit per new application instead of removing all
            this.setUnits([]);
        } else {
            // propagate new application version to all units
            this.unitsInstances
                .filter((unit) => unit.type === UnitType.execution)
                .forEach((unit) => {
                    unit.setApplication({ application });
                });
        }

        this.application = application.toJSON();

        // set model to the default one for the application selected
        this.setModel(
            this.ModelFactory.createFromApplication({
                application: this.application,
            }),
        );
    }

    setModel(model: Model) {
        this.modelInstance = model;
        this.model = model.toJSON() as SubworkflowSchema["model"];
    }

    private buildExternalContext(context: SubworkflowExternalContext): ExternalContext {
        const subworkflowContext = this.units.filter(isAssignmentUnitSchema).reduce((acc, u) => {
            return {
                ...acc,
                [u.operand]: u.value,
            };
        }, {} as AssignmentContext);

        return {
            ...context,
            application: this.applicationInstance.toJSON(),
            methodData: this.model.method.data,
            subworkflowContext,
        };
    }

    render(context: SubworkflowExternalContext) {
        const ctx = this.buildExternalContext(context);

        this.unitsInstances.forEach((u) => {
            if (u.type === UnitType.execution) {
                u.render(ctx);
            }
        });

        // Keep serialized `units` in sync with instances after execution units update `input`.
        this.units = this.unitsInstances.map((u) => u.toJSON());
    }

    /**
     * TODO: reuse workflow function instead
     */
    addUnit(unit: AnySubworkflowUnit, index = -1) {
        const { unitsInstances } = this;

        if (unitsInstances.length === 0) {
            this.setUnits([unit]);
        } else {
            if (index >= 0) {
                unitsInstances.splice(index, 0, unit);
            } else {
                unitsInstances.push(unit);
            }
            this.setUnits(unitsInstances);
        }
    }

    private setUnits(units: AnySubworkflowUnit[]) {
        this.unitsInstances = setUnitLinks(units);
        this.units = units.map((x) => x.toJSON());
        this.properties = units.map((x) => x.resultNames).flat();
    }

    removeUnit(flowchartId: string) {
        const previousUnit = this.unitsInstances.find((x) => x.next === flowchartId);

        if (previousUnit) {
            previousUnit.unsetProp("next");
        }

        this.setUnits(this.unitsInstances.filter((x) => x.flowchartId !== flowchartId));
    }

    getUnit(flowchartId: string) {
        return this.unitsInstances.find((x) => x.flowchartId === flowchartId);
    }

    unitIndex(flowchartId: string) {
        return this.units.findIndex((unit) => {
            return unit.flowchartId === flowchartId;
        });
    }

    replaceUnit(index: number, unit: AnySubworkflowUnit) {
        this.unitsInstances[index] = unit;
        this.setUnits([...this.unitsInstances]);
    }

    setIsDraft(bool: boolean) {
        this.isDraft = bool;
    }

    get methodData() {
        return this.modelInstance.Method.data;
    }

    /**
     * @summary
     * Returns object for hashing of the workflow. Meaningful fields are units, app and model.
     * units must be sorted topologically before hashing (already sorted).
     */
    getHashObject() {
        return {
            application: this.applicationInstance.calculateHash(),
            model: this.modelInstance.calculateHash(),
            units: this.unitsInstances.map((u) => u.calculateHash()).join(),
        };
    }

    findUnitById(id: string) {
        // TODO: come back and refactor after converting flowchartId to id
        return this.units.find((u) => u.flowchartId === id);
    }

    findUnitKeyById(id: string) {
        const index = this.units.findIndex((u) => u.flowchartId === id);
        return `units.${index}`;
    }

    private findAssignmentUnitWithTag(tag: UnitTag): AssignmentUnitSchema | undefined {
        return this.units.filter(isAssignmentUnitSchema).find((unit) => unit.tags?.includes(tag));
    }

    get hasConvergence() {
        return !!this.convergenceParam && !!this.convergenceResult;
    }

    get convergenceParam() {
        return this.findAssignmentUnitWithTag(UnitTag.hasConvergenceParam)?.operand;
    }

    get convergenceResult() {
        return this.findAssignmentUnitWithTag(UnitTag.hasConvergenceResult)?.operand;
    }

    convergenceSeries(scopeTrack: JobSchema["scopeTrack"]) {
        const { convergenceParam, convergenceResult } = this;

        if (!convergenceParam || !convergenceResult || !scopeTrack?.length) {
            return [];
        }

        let prevResult: unknown;
        // `scopeTrack` stores per-repetition diffs: each item only carries the global/local keys
        // that were added or changed in that repetition (see UnitEndpoint.saveUnitStatus). Accumulate
        // the global scope across items so each iteration reads the full scope, not just its delta.
        // This also stays correct for legacy full-snapshot scopeTrack data, since re-applying a full
        // snapshot is idempotent.
        const accumulatedGlobal: Record<string, unknown> = {};

        return scopeTrack
            .map((scopeItem, i) => {
                Object.assign(accumulatedGlobal, scopeItem.scope?.global);
                return {
                    x: i,
                    param: accumulatedGlobal[convergenceParam],
                    y: accumulatedGlobal[convergenceResult],
                };
            })
            .filter(({ y }) => {
                const changed = prevResult !== y;
                prevResult = y;
                return changed;
            })
            .map((item, i) => {
                return {
                    x: i + 1,
                    param: item.param,
                    y: item.y,
                };
            });
    }

    updateMethodData(materials: Material[], metaProperties: MetaPropertyHolder[]) {
        const method = this.modelInstance.Method as PseudopotentialMethod;
        const { model } = this;
        const uniqueElements = [...new Set(materials.map((m) => m.uniqueElements).flat())];
        const appName = this.application.name;

        const methodDataItems = metaProperties
            .filter((metaProperty) => {
                return (
                    // @ts-ignore TODO: fix types
                    uniqueElements.includes(metaProperty.data.element) &&
                    metaProperty.data.apps.includes(appName)
                );
            })
            .map((metaProperty) => metaProperty.property);

        if (method.type !== "pseudopotential" || !methodDataItems.length) {
            return;
        }

        const filters = {
            appName,
            exchangeCorrelation: {
                approximation: model.subtype,
                functional: "functional" in model ? model.functional : undefined,
            },
        };

        // We cycle materials in reverse order below b/c of render(),
        // since the default state index is zero, the last material thus corresponds to index 0.
        // Without reversing, context providers in workflow consider material as changed when any update to the workflow
        // is triggered.
        // TODO: figure out how to simplify or remove the need for the above
        (materials || [])
            .concat()
            .reverse()
            .forEach((material) => {
                // updates methodData & overwrites method in subworkflow.model
                method.updateMethodDataByApplicationAndMaterials(methodDataItems, {
                    elements: material.uniqueElements,
                    ...filters,
                });

                this.modelInstance.setMethod(method);
                this.model = this.modelInstance.toJSON() as SubworkflowSchema["model"];
            });

        // TODO: Try if/else instead of running both
        if (materials.length > 1) {
            method.updateMethodDataByApplicationAndMaterials(methodDataItems, {
                elements: uniqueElements,
                ...filters,
            });

            this.modelInstance.setMethod(method);
            this.model = this.modelInstance.toJSON() as SubworkflowSchema["model"];
        }
    }

    addConvergence({
        parameter,
        parameterInitial,
        parameterIncrement,
        result,
        resultInitial,
        condition,
        operator,
        tolerance,
        maxOccurrences,
        externalContext,
    }: ConvergenceConfig) {
        // Find unit to converge: should contain passed result in its results list
        // TODO: make user to select unit for convergence explicitly
        const unitForConvergence = this.unitsInstances
            .filter((x) => x.type === UnitType.execution)
            .find((x) => {
                return x.resultNames.find((name) => name === result);
            });

        if (!unitForConvergence) {
            throw new Error(
                `Subworkflow does not contain unit with '${result}' as extracted property.`,
            );
        }

        // initialize parameter
        const convergenceParameter = createConvergenceParameter({
            name: parameter,
            initialValue: parameterInitial,
            increment: parameterIncrement,
        });

        const context = this.buildExternalContext(externalContext);

        unitForConvergence.render(context, convergenceParameter);

        const prevResult = "prev_result";
        const iteration = "iteration";

        // Assignment with result's initial value
        const prevResultInit = new AssignmentUnit({
            name: "init result",
            head: true,
            operand: prevResult,
            value: resultInitial,
        });

        // Assignment with initial value of convergence parameter
        const paramInit = new AssignmentUnit({
            name: "init parameter",
            operand: convergenceParameter.name,
            value: convergenceParameter.initialValue,
            tags: [UnitTag.hasConvergenceParam],
        });

        // Assignment with initial value of iteration counter
        const iterInit = new AssignmentUnit({
            name: "init counter",
            operand: iteration,
            value: 1,
        });

        // Assignment for storing iteration result: extracts 'result' from convergence unit scope
        const storePrevResult = new AssignmentUnit({
            name: "store result",
            input: [
                {
                    scope: unitForConvergence.flowchartId,
                    name: result,
                },
            ],
            operand: prevResult,
            value: result,
        });

        // Assignment for convergence param increase
        const nextStep = new AssignmentUnit({
            name: "update parameter",
            input: convergenceParameter.useVariablesFromUnitContext(unitForConvergence.flowchartId),
            operand: convergenceParameter.name,
            value: convergenceParameter.increment,
            next: unitForConvergence.flowchartId,
        });

        // Final step of convergence
        const exit = new AssignmentUnit({
            name: "exit",
            input: [],
            operand: convergenceParameter.name,
            value: convergenceParameter.finalValue,
        });

        // Final step of convergence
        const storeResult = new AssignmentUnit({
            name: "update result",
            input: [
                {
                    scope: unitForConvergence.flowchartId,
                    name: result,
                },
            ],
            operand: result,
            value: result,
            tags: [UnitTag.hasConvergenceResult],
        });

        // Assign next iteration value
        const nextIter = new AssignmentUnit({
            name: "update counter",
            input: [],
            operand: iteration,
            value: `${iteration} + 1`,
        });

        // Convergence condition unit
        const conditionUnit = new ConditionUnit({
            name: "check convergence",
            statement: `${condition} ${operator} ${tolerance}`,
            then: exit.flowchartId,
            else: storePrevResult.flowchartId,
            maxOccurrences,
            next: storePrevResult.flowchartId,
        });

        this.addUnit(paramInit, 0);
        this.addUnit(prevResultInit, 1);
        this.addUnit(iterInit, 2);
        this.addUnit(storeResult);
        this.addUnit(conditionUnit);
        this.addUnit(storePrevResult);
        this.addUnit(nextIter);
        this.addUnit(nextStep);
        this.addUnit(exit);

        // `addUnit` adjusts the `next` field, hence the below.
        nextStep.next = unitForConvergence.flowchartId;
    }
}

namedEntityMixin(Subworkflow.prototype);
defaultableEntityMixin(Subworkflow);
computedEntityMixin(Subworkflow.prototype);
subworkflowSchemaMixin(Subworkflow.prototype);
hashedEntityMixin(Subworkflow.prototype);

export default Subworkflow;

import { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import {
    type Defaultable,
    defaultableEntityMixin,
} from "@mat3ra/code/dist/js/entity/mixins/DefaultableMixin";
import {
    type HasDescription,
    hasDescriptionMixin,
} from "@mat3ra/code/dist/js/entity/mixins/HasDescriptionMixin";
import {
    type HashedEntity,
    hashedEntityMixin,
} from "@mat3ra/code/dist/js/entity/mixins/HashedEntityMixin";
import {
    type NamedEntity,
    namedEntityMixin,
} from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";
import { Taggable, taggableMixin } from "@mat3ra/code/dist/js/entity/mixins/TaggableMixin";
import JSONSchemasInterface from "@mat3ra/esse/dist/js/esse/JSONSchemasInterface";
import type {
    ApplicationSchema,
    BaseInMemoryEntitySchema,
    SubworkflowSchema,
} from "@mat3ra/esse/dist/js/types";
import { ComputedEntityMixin, computedEntityMixin } from "@mat3ra/ide/dist/js/compute";
import type { Material } from "@mat3ra/made";
import type { MetaPropertyHolder } from "@mat3ra/prode";
import { setUnitLinks, SubworkflowStandata } from "@mat3ra/standata";
import { Utils } from "@mat3ra/utils";

import type { MaterialExternalContext } from "./context/mixins/MaterialContextMixin";
import type { MaterialsExternalContext } from "./context/mixins/MaterialsContextMixin";
import type { MaterialsSetExternalContext } from "./context/mixins/MaterialsSetContextMixin";
import type { JobExternalContext } from "./context/providers/base/ContextProvider";
import { UnitType } from "./enums";
import { type WorkflowSchemaMixin, workflowSchemaMixin } from "./generated/WorkflowSchemaMixin";
import Subworkflow from "./Subworkflow";
import { MapUnit } from "./units";
import { type AnyWorkflowUnit, UnitFactory } from "./units/factory";
import {
    getDefaultDescription,
    getHumanReadableProperties,
    getHumanReadableUsedModels,
    getProperties,
    getSystemName,
    getUsedApplications,
    getUsedModels,
} from "./utils/workflow";
import defaultWorkflowConfig from "./workflows/default";
import type { WorkflowSchema } from "./workflows/types";

export type WorkflowEntity = WorkflowSchema & BaseInMemoryEntitySchema;

type Schema = WorkflowEntity;

interface Workflow
    extends Defaultable,
        NamedEntity,
        WorkflowSchemaMixin,
        Taggable,
        HashedEntity,
        ComputedEntityMixin<WorkflowSchema["compute"]>,
        HasDescription {}

/** Context passed to Workflow.render(); subworkflows also receive `workflowHasRelaxation`. */
export type WorkflowRenderContext = MaterialExternalContext &
    MaterialsExternalContext &
    MaterialsSetExternalContext &
    JobExternalContext & {
        scopeGlobal?: Record<string, unknown>;
    };

class Workflow<S extends Schema = Schema> extends InMemoryEntity<S> implements WorkflowSchema {
    declare createDefault: () => Workflow;

    static readonly defaultConfig = defaultWorkflowConfig;

    static get jsonSchema() {
        return JSONSchemasInterface.getSchemaById("workflow");
    }

    subworkflowInstances: Subworkflow[];

    unitInstances!: AnyWorkflowUnit[];

    workflowInstances: Workflow[];

    repetition = 0;

    totalRepetitions = 1;

    setTotalRepetitions(totalRepetition: number) {
        this.totalRepetitions = totalRepetition;
    }

    setRepetition(repetition: number) {
        this.repetition = repetition;
        this.unitInstances.forEach((u) => u.setRepetition(repetition));
        this.subworkflowInstances.forEach((sw) => sw.setRepetition(repetition));
        this.workflowInstances.forEach((wf) => wf.setRepetition(repetition));
    }

    static fromSubworkflow(subworkflow: Subworkflow) {
        const config = {
            name: subworkflow.name,
            subworkflows: [subworkflow.toJSON()],
            units: [subworkflow.getAsUnit().toJSON()],
            properties: subworkflow.properties,
            applicationName: subworkflow.application.name,
            workflows: [] as WorkflowSchema[],
        };
        return new this(config);
    }

    // NoInfer: keep default S (or an explicit type arg) instead of inferring S from the config literal.
    constructor(config: NoInfer<S> & { applicationName?: string }) {
        // Strip applicationName (fromSubworkflow helper); not part of WorkflowEntity.
        // eslint-disable-next-line @typescript-eslint/no-unused-vars
        const { applicationName, ...rest } = config;
        super({
            ...rest,
            _id: (rest as Partial<Schema>)._id || Utils.uuid.getUUID(),
        } as S);

        this.subworkflowInstances = this.subworkflows.map((x) => new Subworkflow(x));
        this.workflowInstances = this.workflows?.map((x) => new Workflow(x)) || [];
        this.setUnits(this.units.map((unit) => UnitFactory.createInWorkflow(unit)));
    }

    get workflows() {
        return this.requiredProp("workflows");
    }

    set workflows(value: WorkflowSchema[]) {
        (this._json as Schema).workflows = value;
    }

    addSubworkflow(subworkflow: Subworkflow, head = false, index = -1) {
        const subworkflowUnit = subworkflow.getAsUnit();

        if (head) {
            this.subworkflowInstances.unshift(subworkflow);
            this.addUnit(subworkflowUnit, head, index);
        } else {
            this.subworkflowInstances.push(subworkflow);
            this.addUnit(subworkflowUnit, head, index);
        }
    }

    updateMethodData(materials: Material[], metaProperties: MetaPropertyHolder[]) {
        this.subworkflowInstances.forEach((sw) => {
            sw.updateMethodData(materials, metaProperties);
        });
    }

    removeSubworkflow(id: string) {
        const subworkflowUnit = this.unitInstances.find((u) => u.id === id);

        if (subworkflowUnit) {
            this.removeUnit(subworkflowUnit.flowchartId);
        }
    }

    setUnits(arr: AnyWorkflowUnit[]) {
        this.unitInstances = setUnitLinks(arr);
        this.units = this.unitInstances.map((u) => u.toJSON());
    }

    /**
     * A top-level subworkflow branch is represented twice: a {@link UnitType.subworkflow} unit in
     * `unitInstances` (flowchart card) and a full {@link Subworkflow} in `subworkflowInstances`
     * (linked by the same `id` as {@link Subworkflow.getAsUnit}). Their display names must match
     * so `toJSON()` and editors stay consistent. Call this after mutating a subworkflow unit's
     * `name` (and updating `unitInstances` via {@link setUnits}).
     */
    syncLinkedSubworkflowNameFromUnit(unit: AnyWorkflowUnit): void {
        if (unit.type !== UnitType.subworkflow) return;
        const linked = this.subworkflowInstances.find((s) => s.id === unit.id);
        linked?.setName(unit.name);
    }

    render(context: WorkflowRenderContext) {
        this.subworkflowInstances.forEach((sw) => {
            sw.render({
                ...context,
                workflowHasRelaxation: this.hasRelaxation,
            });
        });
        this.workflowInstances.forEach((wf) => {
            wf.render(context);
        });
    }

    get usedApplications(): ApplicationSchema[] {
        return getUsedApplications(this);
    }

    // return application names
    get usedApplicationNames() {
        return this.usedApplications.map((a) => a.name);
    }

    get usedApplicationVersions() {
        return this.usedApplications.map((a) => a.version);
    }

    get usedApplicationNamesWithVersions() {
        return this.usedApplications.map((a) => `${a.name} ${a.version}`);
    }

    getUsedModels() {
        return getUsedModels(this);
    }

    getHumanReadableUsedModels() {
        return getHumanReadableUsedModels(this);
    }

    toJSON(): S {
        return {
            ...super.toJSON(),
            name: this.name,
            properties: getProperties(this),
            units: this.unitInstances.map((x) => x.toJSON()),
            subworkflows: this.subworkflowInstances.map((x) => x.toJSON()),
            workflows: this.workflowInstances.map((x) => x.toJSON()),
        } as S;
    }

    getHumanReadableProperties() {
        return getHumanReadableProperties(this);
    }

    getProperties() {
        return getProperties(this);
    }

    getSystemName() {
        return getSystemName(this);
    }

    getDefaultDescription() {
        return getDefaultDescription(this);
    }

    private addUnit(unit: AnyWorkflowUnit, head = false, index = -1) {
        const [...unitInstances] = this.unitInstances;

        if (unitInstances.length === 0) {
            this.setUnits([unit]);
        } else {
            if (head) {
                unitInstances.unshift(unit);
            } else if (index >= 0) {
                unitInstances.splice(index, 0, unit);
            } else {
                unitInstances.push(unit);
            }
            this.setUnits(unitInstances);
        }
    }

    removeUnit(flowchartId: string) {
        if (this.units.length < 2) {
            return;
        }

        const unit = this.unitInstances.find((x) => x.flowchartId === flowchartId);
        if (!unit) {
            return;
        }

        const previousUnit = this.unitInstances.find((x) => x.next === unit.flowchartId);
        if (previousUnit) {
            delete previousUnit.next;
        }

        this.subworkflowInstances = this.subworkflowInstances.filter((x) => x.id !== unit.id);

        this.setUnits(this.unitInstances.filter((x) => x.flowchartId !== flowchartId));
    }

    /*
     * @param type {String|Object} Unit type, map or subworkflow
     * @param head {Boolean}
     * @param index {Number} Index at which the unit will be added. -1 by default (ignored).
     */
    addUnitType(type: UnitType, head = false, index = -1) {
        switch (type) {
            case UnitType.map: {
                const mapWorkflowConfig = {
                    ...defaultWorkflowConfig,
                    _id: Utils.uuid.getUUID(),
                };
                const mapUnit = new MapUnit({
                    workflowId: mapWorkflowConfig._id,
                });

                this.workflows = [...(this.workflows || []), mapWorkflowConfig];
                this.workflowInstances = this.workflows.map((x) => new Workflow(x));

                this.addUnit(mapUnit, head, index);

                break;
            }
            case UnitType.subworkflow:
                this.addSubworkflow(Subworkflow.createDefault(), head, index);
                break;
            default:
                console.log(`unit_type=${type} unrecognized, skipping.`);
        }
    }

    addMapUnit(mapUnit: MapUnit, mapWorkflow: Workflow) {
        const mapWorkflowConfig = {
            _id: Utils.uuid.getUUID(),
            ...mapWorkflow.toJSON(),
        };

        mapUnit.setWorkflowId(mapWorkflowConfig._id);
        this.addUnit(mapUnit);

        this.workflows = [...(this.workflows || []), mapWorkflowConfig];
        this.workflowInstances = this.workflows.map((x) => new Workflow(x));
    }

    get allSubworkflows() {
        const subworkflowsList: Subworkflow[] = [];
        this.subworkflowInstances.forEach((sw) => subworkflowsList.push(sw));
        this.workflowInstances.forEach((workflow) => {
            return Array.prototype.push.apply(subworkflowsList, workflow.allSubworkflows);
        });
        return subworkflowsList;
    }

    get hasRelaxation() {
        return Boolean(this.getRelaxationSubworkflow());
    }

    toggleRelaxation() {
        const relaxSubworkflow = this.getRelaxationSubworkflow();
        if (relaxSubworkflow?._id) {
            this.removeSubworkflow(relaxSubworkflow._id);
        } else {
            const vcRelax = this.getStandataRelaxationSubworkflow();
            if (vcRelax) {
                const application = structuredClone(this.subworkflowInstances[0].application);
                this.addSubworkflow(
                    new Subworkflow({
                        ...vcRelax,
                        application,
                        units: vcRelax.units.map((unit) =>
                            unit.type === UnitType.execution ? { ...unit, application } : unit,
                        ),
                    }),
                    true,
                );
            }
        }
    }

    getHashObject() {
        return {
            units: this.unitInstances.map((u) => u.calculateHash()).join(),
            subworkflows: this.subworkflowInstances.map((sw) => sw.calculateHash()).join(),
            workflows: this.workflowInstances.map((w) => w.calculateHash()).join(),
        };
    }

    private getStandataRelaxationSubworkflow() {
        // TODO: fix standata type
        const subworkflow = new SubworkflowStandata().getRelaxationSubworkflowByApplication(
            this.subworkflowInstances[0].application.name,
        ) as unknown as SubworkflowSchema | undefined;

        if (!subworkflow) {
            return undefined;
        }

        const hasExecutionUnit = subworkflow.units.some((unit) => unit.type === UnitType.execution);
        if (!hasExecutionUnit) {
            throw new Error("Relaxation subworkflow is missing an execution unit");
        }

        return subworkflow;
    }

    private getRelaxationSubworkflow() {
        const standataSubworkflow = this.getStandataRelaxationSubworkflow();

        return this.subworkflowInstances.find((sw) => {
            return standataSubworkflow && standataSubworkflow.systemName === sw.systemName;
        });
    }
}

namedEntityMixin(Workflow.prototype);
workflowSchemaMixin(Workflow.prototype);
taggableMixin(Workflow.prototype);
computedEntityMixin(Workflow.prototype);
defaultableEntityMixin(Workflow);
hashedEntityMixin(Workflow.prototype);
hasDescriptionMixin(Workflow.prototype);

export default Workflow;

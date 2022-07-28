import lodash from "lodash";
import s from "underscore.string";
import { mix } from "mixwith";

import { NamedDefaultableRepetitionContextAndRenderInMemoryEntity } from "@exabyte-io/code.js/dist/entity";
import { getUUID } from "@exabyte-io/code.js/dist/utils";
import { UNIT_TYPES } from "../enums";
import { tree } from "@exabyte-io/mode.js";
import { ComputedEntityMixin, getDefaultComputeConfig } from "@exabyte-io/ide.js";

import { Subworkflow } from "../subworkflows/subworkflow";
import { UnitFactory } from "../units/factory";
import { MapUnit } from "../units";
import { setNextLinks, setUnitsHead } from "../utils";
import defaultWorkflowConfig from "./default";
import _ from "underscore";
import { calculateHashFromObject } from "@exabyte-io/code.js/dist/utils";
import { createSubworkflowByName } from "../subworkflows";
import { RelaxationLogicMixin } from "./relaxation";

const { MODEL_NAMES } = tree;

class BaseWorkflow extends mix(NamedDefaultableRepetitionContextAndRenderInMemoryEntity).with(
    ComputedEntityMixin,
    RelaxationLogicMixin,
) {
}


export class Workflow extends BaseWorkflow {
    static getDefaultComputeConfig = getDefaultComputeConfig;

    constructor(config) {
        super(config);
        this._Subworkflow = Subworkflow;
        this._UnitFactory = UnitFactory;
        this._Workflow = Workflow;
        this._MapUnit = MapUnit;
        this.initialize();
    }

    // get _allRelaxationSubworkflows() {
    //     return {
    //         espresso: createSubworkflowByName({ appName: "espresso", swfName: "variable_cell_relaxation" }),
    //         vasp: createSubworkflowByName({ appName: "vasp", swfName: "variable_cell_relaxation" }),
    //     }
    // }

    initialize() {
        const me = this;
        this._subworkflows = this.prop("subworkflows").map(x => new me._Subworkflow(x));
        this._units = this.prop("units").map(unit => me._UnitFactory.create(unit));
        this._json.workflows = (this._json.workflows || []);
        this._workflows = this.prop("workflows").map(x => new me._Workflow(x));
    }

    static get defaultConfig() {
        return defaultWorkflowConfig;
    }

    static generateWorkflowId() {
        return getUUID();
    }

    static fromSubworkflow(subworkflow, cls = Workflow) {
        const config = {
            name: subworkflow.name,
            subworkflows: [subworkflow.toJSON()],
            units: setNextLinks(setUnitsHead([subworkflow.getAsUnit().toJSON()])),
            properties: subworkflow.properties,
        };
        return new cls(config);
    }

    static fromSubworkflows(name, cls = Workflow, ...subworkflows) {
        return new cls(name, subworkflows, subworkflows.map(sw => sw.getAsUnit()));
    }

    /**
     * @summary Adds subworkflow to current workflow.
     * @param subworkflow {Subworkflow}
     * @param head {Boolean}
     */
    addSubworkflow(subworkflow, head = false, index = -1) {
        const subworkflowUnit = subworkflow.getAsUnit();
        if (head) {
            this.subworkflows.unshift(subworkflow);
            this.addUnit(subworkflowUnit, head, index);
        } else {
            this.subworkflows.push(subworkflow);
            this.addUnit(subworkflowUnit, head, index);
        }
    }

    removeSubworkflow(id) {
        const subworkflowUnit = this.units.find(u => u.id === id);
        subworkflowUnit && this.removeUnit(subworkflowUnit.flowchartId);
    }

    subworkflowId(index) {
        const sw = this.prop(`subworkflows[${index}]`);
        return sw ? sw._id : null;
    }

    replaceSubworkflowAtIndex(index, newSubworkflow) {
        this._subworkflows[index] = newSubworkflow;
        this.setUnits(setNextLinks(setUnitsHead(this._units)));
    }

    get units() {
        return this._units;
    }

    setUnits(arr) {
        this._units = arr;
    }

    // returns a list of `app` Classes
    get usedApplications() {
        const swApplications = this.subworkflows.map(sw => sw.application);
        const wfApplications = lodash.flatten(this.workflows.map(w => w.usedApplications));
        return lodash.uniqBy(swApplications.concat(wfApplications), (a) => a.name);
    }

    // return application names
    get usedApplicationNames() {
        return this.usedApplications.map(a => a.name);
    }

    get usedApplicationVersions() {
        return this.usedApplications.map(a => a.version);
    }

    get usedApplicationNamesWithVersions() {
        return this.usedApplications.map(a => `${a.name} ${a.version}`);
    }

    get usedModels() {
        return lodash.uniq(this.subworkflows.map(sw => sw.model.type));
    }

    get humanReadableUsedModels() {
        return this.usedModels.filter(m => m !== "unknown").map(m => MODEL_NAMES[m]);
    }

    toJSON(exclude = []) {
        return lodash.omit(Object.assign({}, super.toJSON(), {
            units: this._units.map(x => x.toJSON()),
            subworkflows: this._subworkflows.map(x => x.toJSON()),
            workflows: this.workflows.map(x => x.toJSON()),
            compute: this.compute
        }), exclude);
    }

    get isDefault() {
        return this.prop("isDefault", false);
    }

    // returns true if any of subworkflows are multimaterial
    get isMultiMaterial() {
        const _getMM = (x) => lodash.get(x, "isMultiMaterial", false);
        const fromSubworkflows = this.subworkflows.reduce((a, b) => _getMM(a) || _getMM(b), 0);
        return this.prop("isMultiMaterial") || fromSubworkflows;
    }

    set isUsingDataset(value){
        this.setProp("isUsingDataset", value);
    }

    get isUsingDataset(){
        return !!this.prop("isUsingDataset", false);
    }

    get properties() {
        return lodash.uniq(lodash.flatten(this._subworkflows.map(x => x.properties)));
    }

    get humanReadableProperties() {
        return this.properties.map(name => s.humanize(name));
    }

    get systemName() {
        return s.slugify(`${this.usedApplicationNames.join(":")}-${this.name.toLowerCase()}`);
    }

    get defaultDescription() {
        return `${this.usedModels.join(", ").toUpperCase()} workflow using ${this.usedApplicationNames.join(", ")}.`;
    }

    get exabyteId() {
        return this.prop("exabyteId");
    }

    get hash() {
        return this.prop("hash", "");
    }

    get isOutdated() {
        return this.prop("isOutdated", false);
    }

    get history() {
        return this.prop("history", [])
    }


    setMethodData(methodData) {
        this.subworkflows.forEach(sw => {
            const method = methodData.getMethodBySubworkflow(sw);
            method && sw.model.setMethod(method);
        });

        this.workflows.forEach(wf => {
            wf.subworkflows.forEach(sw => {
                const method = methodData.getMethodBySubworkflow(sw);
                method && sw.model.setMethod(method);
            });
        });
    }

    /**
     * @param unit {Unit}
     * @param head {Boolean}
     * @param index {Number}
     */
    addUnit(unit, head = false, index = -1) {
        const units = this.units;
        if (units.length === 0) {
            unit.head = true;
            this.setUnits([unit]);
        } else {
            if (head) {
                const first = lodash.first(units);
                units.unshift(unit);
            } else {
                const last = lodash.last(units);
                (index >= 0) ? units.splice(index, 0, unit) : units.push(unit);
            }
            this.setUnits(setNextLinks(setUnitsHead(units)));
        }
    }

    removeUnit(flowchartId) {
        if (this.units.length < 2) return;

        const unit = this.units.find(x => x.flowchartId === flowchartId);
        const previousUnit = this.units.find(x => x.next === unit.flowchartId);
        if (previousUnit) {
            delete previousUnit.next;
        }

        this._subworkflows = this._subworkflows.filter(x => x.id !== unit.id);
        this._units = setNextLinks(setUnitsHead(this._units.filter(x => x.flowchartId !== flowchartId)));

    }

    /**
     * @return Subworkflow[]
     */
    get subworkflows() {
        return this._subworkflows;
    }

    get workflows() {
        return this._workflows;
    }

    /*
     * @param type {String|Object} Unit type, map or subworkflow
     * @param head {Boolean}
     * @param index {Number} Index at which the unit will be added. -1 by default (ignored).
     */
    addUnitType(type, head = false, index = -1) {
        switch (type) {
            case UNIT_TYPES.map:
                const workflowConfig = defaultWorkflowConfig;
                workflowConfig._id = this._Workflow.generateWorkflowId();
                this.prop("workflows").push(workflowConfig);
                this._workflows = this.prop("workflows").map(x => new this._Workflow(x));
                const mapUnit = new this._MapUnit();
                mapUnit.setWorkflowId(workflowConfig._id);
                this.addUnit(mapUnit, head, index);
                break;
            case UNIT_TYPES.subworkflow:
                this.addSubworkflow(this._Subworkflow.createDefault(), head, index);
                break;
            default:
                console.log(`unit_type=${type} unrecognized, skipping.`);
        }
    }

    addMapUnit(mapUnit, mapWorkflow) {
        const mapWorkflowConfig = mapWorkflow.toJSON();
        if (!mapWorkflowConfig._id) mapWorkflowConfig._id = this._Workflow.generateWorkflowId();
        mapUnit.setWorkflowId(mapWorkflowConfig._id);
        this.addUnit(mapUnit);
        this._json.workflows.push(mapWorkflowConfig);
        const me = this;
        this._workflows = this.prop("workflows").map(x => new me._Workflow(x));
    }

    findSubworkflowById(id) {
        if (!id) return;

        const workflows = this.workflows || [];
        const subworkflows = this.subworkflows || [];

        const subworkflow = subworkflows.find(subworkflow => subworkflow.id === id);
        if (subworkflow) return subworkflow;

        const workflow = workflows.find(w => w.findSubworkflowById(id));
        if (workflow) return workflow.findSubworkflowById(id);

        console.warn("attempted to find a non-existing subworkflow");
    }

    get allSubworkflows() {
        let subworkflowsList = [];
        this.subworkflows.forEach(sw => subworkflowsList.push(sw));
        this.workflows.forEach(workflow => Array.prototype.push.apply(subworkflowsList, workflow.allSubworkflows));
        return subworkflowsList;
    }

    /**
     * @summary Calculates hash of the workflow. Meaningful fields are units and subworkflows.
     * units and subworkflows must be sorted topologically before hashing (already sorted).
     */
    calculateHash() {
        const meaningfulFields = {
            units: _.map(this.units, u => u.calculateHash()).join(),
            subworkflows: _.map(this.subworkflows, s => s.calculateHash()).join(),
            workflows: _.map(this.workflows, w => w.calculateHash()).join(),
        };
        return calculateHashFromObject(meaningfulFields);
    }

}

import lodash from "lodash";
import { mix } from "mixwith";

import { NamedDefaultableRepetitionImportantSettingsInMemoryEntity } from "@exabyte-io/code.js/dist/entity";
import { calculateHashFromObject, getUUID, removeTimestampableKeysFromConfig } from "@exabyte-io/code.js/dist/utils";

import { Application } from "@exabyte-io/ade.js";
import { UnitFactory } from "/imports/wode";
import { UNIT_TYPES } from "/imports/wode/enums";
import { Model, ModelFactory } from "@exabyte-io/mode.js";

import { setUnitsHead, setNextLinks } from "../utils";
import { ConvergenceMixin } from "./convergence";
import _ from "underscore";


class BaseSubworkflow extends mix(NamedDefaultableRepetitionImportantSettingsInMemoryEntity).with(ConvergenceMixin) {}


export class Subworkflow extends BaseSubworkflow {
    constructor(config) {
        super(config);
        this._Application = Application;
        this._ModelFactory = ModelFactory;
        this._UnitFactory = UnitFactory;
        this.initialize();
    }

    initialize() {
        this._application = new this._Application(this.prop("application"));
        this._model = this._ModelFactory.create({
            ...this.prop("model"),
            application: this.prop("application"),
        });
        this._units = setNextLinks(setUnitsHead(this.prop("units", [])), this.id).map(
            (cfg) => this._UnitFactory.create(Object.assign(cfg, { application: this.application.toJSON() }))
        );
    }

    static generateSubworkflowId() {
        return getUUID();
    }

    static get defaultConfig() {
        return {
            _id: this.generateSubworkflowId(),
            name: "New Subworkflow",
            application: Application.defaultConfig,
            model: Model.defaultConfig,
            properties: [],
            units: [],
        }
    }

    /*
     * @returns {SubworkflowUnit}
     */
    getAsUnit() {
        return this._UnitFactory.create({
            type: UNIT_TYPES.subworkflow,
            _id: this.id,
            name: this.name,
        });
    }

    /*
     * @summary Used to generate initial application tree, therefore omit setting application.
     */
    static fromArguments(application, model, method, name, units = [], config = {}, cls = Subworkflow) {
        return new cls({
            ...config,
            _id: cls.generateSubworkflowId(),
            name,
            application: application.toJSON(),
            properties: lodash.sortedUniq(lodash.flatten(units.map(x => x.resultNames))),
            model: {
                ...model.toJSON(),
                method: method.toJSON()
            },
            units,
        });
    }

    get application() {
        return this._application;
    }

    setApplication(application) {
        // TODO: adjust the logic above to take into account whether units need re-rendering after version change etc.
        // reset units if application name changes
        if (this.application.name !== application.name) this.setUnits([]);
        this._application = application;
        this.setProp('application', application.toJSON());
        // set model to the default one for the application selected
        this.setModel(this._ModelFactory.createFromApplication({application: this.prop('application')}));

    }

    get model() {
        return this._model;
    }

    setModel(model) {
        this._model = model;
    }

    get units() {
        return this._units;
    }

    setUnits(units) {
        this._units = units;
    }

    toJSON(exclude = []) {
        return Object.assign({}, super.toJSON(exclude), {
            model: this.model.toJSON(),
            units: this.units.map(x => x.toJSON()),
            compute: this.compute,
        });
    }

    get contextProviders() {
        const unitsWithContextProviders = this.units.filter(u => (u.allContextProviders && u.allContextProviders.length));
        const allContextProviders = _.flatten(unitsWithContextProviders.map(u => u.allContextProviders));
        const subworkflowContextProviders = allContextProviders.filter(p => p.isSubworkflowContextProvider);
        return _.uniq(subworkflowContextProviders, p => p.name);

    }

    /**
     * Extracts a reduced version of the entity config to be stored inside redux state.
     * This is used to track changes to context, monitors, properties, etc. when multiple materials are in state.
     */
    extractReducedExternalDependentConfig() {
        return {
            id: this.id,
            context: this.context || {},
            units: this.units.map(unit => unit.extractReducedExternalDependentConfig()),
        }
    }

    /**
     * Applies the reduced config obtained from extractReducedExternalDependentConfig on the entity.
     */
    applyReducedExternalDependentConfig(config) {
        this.context = config.context || {};
        this.units.forEach(unit => {
            const unitConfig = (config.units || []).find(c => c.id === unit.flowchartId);
            unit.applyReducedExternalDependentConfig(unitConfig || {});
        });
    }

    render(context = {}) {
        const ctx = Object.assign({},
            context,
            {
                application: this.application,
                methodData: this.model.method.data,
                model: this.model.toJSON()
            },
            // context below is assembled from context providers and passed to units to override theirs
            this.context,
        );

        this.units.forEach(u => u.render(ctx));
    }

    /**
     * TODO: reuse workflow function instead
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
                unit.next = first.flowchartId;
                units.unshift(unit);
            } else {
                const last = lodash.last(units);
                last.next = unit.flowchartId;
                (index >= 0) ? units.splice(index, 0, unit) : units.push(unit);
            }
            this.setUnits(setNextLinks(setUnitsHead(units)));
        }
    }

    removeUnit(flowchartId) {
        const previousUnit = this.units.find(x => x.next === flowchartId);
        if (previousUnit) previousUnit.unsetProp('next');
        // TODO: remove the setNextLinks and setUnitsHead and handle the logic via flowchart designer
        this.setUnits(setNextLinks(setUnitsHead(this.units.filter(x => x.flowchartId !== flowchartId))));
    }

    get properties() {
        return lodash.flatten(this.units.map(x => x.resultNames));
    }

    getUnit(flowchartId) {
        return this.units.find(x => x.flowchartId === flowchartId);
    }

    unitIndex(flowchartId) {
        return lodash.findIndex(this.units, (unit) => {
            return unit.flowchartId === flowchartId;
        })
    }

    replaceUnit(index, unit) {
        this.units[index] = unit;
        this.setUnits(setNextLinks(setUnitsHead(this.units)));
    }

    get scopeVariables() {
        return ['N_k'];//this.units.filter(x => x.type === UNIT_TYPES.assignment).map(x => x.name);
    }

    get scalarResults() {
        return ['total_energy', 'pressure'];
    }

    /**
     * @return {Model}
     */
    get isMultiMaterial() {
        return this.prop('isMultiMaterial', false);
    }

    get isDraft() {return this.prop('isDraft', false)}

    setIsDraft(bool) {
        return this.setProp('isDraft', bool);
    }

    get methodData() {
        return this.model.method.data;
    }

    /**
     * @summary Calculates hash of the subworkflow. Meaningful fields are units, app and model.
     * units must be sorted topologically before hashing (already sorted).
     */
    calculateHash() {
        const config = this.toJSON();
        const meaningfulFields = {
            application: removeTimestampableKeysFromConfig(config.application),
            model: this._calculateModelHash(),
            units: _.map(this.units, u => u.calculateHash()).join()
        };
        return calculateHashFromObject(meaningfulFields);
    }

    _calculateModelHash() {
        const model = this.toJSON().model;
        // ignore empty data object
        if (this.model.method.omitInHashCalculation) delete model.method.data;
        return calculateHashFromObject(model);
    }

    findUnitById(id) {
        // TODO: come back and refactor after converting flowchartId to id
        return this.units.find(u => u.flowchartId === id);
    }

    findUnitKeyById(id) {
        const index = this.units.findIndex(u => u.flowchartId === id);
        return `units.${index}`;
    }


}

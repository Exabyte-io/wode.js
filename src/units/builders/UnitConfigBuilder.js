import _ from "underscore";
import { getUUID } from "@exabyte-io/code.js/dist/utils";


export class UnitConfigBuilder {
    constructor({ name, type }) {
        this.type = type;
        this._name = name;
        this._head = false;
        this._results = [];
        this._monitors = [];
        this._preProcessors = [];
        this._postProcessors = [];
        this._flowchartId = this.constructor.defaultFlowchartId();
    }

    name(str) {
        this._name = str;
        return this;
    }

    head(bool) {
        this._head = bool;
        return this;
    }

    static defaultFlowchartId() {
        return getUUID();
    }

    flowchartId(flowchartId) {
        this._flowchartId = flowchartId;
        return this;
    }

    _stringArrayToNamedObject(array) {
        return array.map(name => _.isString(name) ? {name} : name);
    }

    addPreProcessors(preProcessorNames) {
        this._preProcessors = _.union(this._stringArrayToNamedObject(preProcessorNames), this._preProcessors);
        return this;
    }

    addPostProcessors(postProcessorNames) {
        this._postProcessors = _.union(this._stringArrayToNamedObject(postProcessorNames), this._postProcessors);
        return this;
    }

    addResults(resultNames) {
        this._results = _.union(this._stringArrayToNamedObject(resultNames), this._results);
        return this;
    }

    addMonitors(monitorNames) {
        this._monitors = _.union(this._stringArrayToNamedObject(monitorNames), this._monitors);
        return this;
    }

    build() {
        return {
            type: this.type,
            name: this._name,
            head: this._head,
            results: this._results,
            monitors: this._monitors,
            flowchartId: this._flowchartId,
            preProcessors: this._preProcessors,
            postProcessors: this._postProcessors,
        };
    }
}

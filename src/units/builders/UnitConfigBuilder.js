import { getUUID } from "@mat3ra/code/dist/js/utils";
import _ from "underscore";

export class UnitConfigBuilder {
    constructor({ name, type, flowchartId }) {
        this.type = type;
        this._name = name;
        this._head = false;
        this._results = [];
        this._monitors = [];
        this._preProcessors = [];
        this._postProcessors = [];
        this._flowchartId = flowchartId || this.constructor.generateFlowChartId();
    }

    name(str) {
        this._name = str;
        return this;
    }

    head(bool) {
        this._head = bool;
        return this;
    }

    static generateFlowChartId() {
        return getUUID();
    }

    flowchartId(flowchartId) {
        this._flowchartId = flowchartId;
        return this;
    }

    static _stringArrayToNamedObject(array) {
        return array.map((name) => (_.isString(name) ? { name } : name));
    }

    addPreProcessors(preProcessorNames) {
        this._preProcessors = _.union(
            this.constructor._stringArrayToNamedObject(preProcessorNames),
            this._preProcessors,
        );
        return this;
    }

    addPostProcessors(postProcessorNames) {
        this._postProcessors = _.union(
            this.constructor._stringArrayToNamedObject(postProcessorNames),
            this._postProcessors,
        );
        return this;
    }

    addResults(resultNames) {
        this._results = _.union(
            this.constructor._stringArrayToNamedObject(resultNames),
            this._results,
        );
        return this;
    }

    addMonitors(monitorNames) {
        this._monitors = _.union(
            this.constructor._stringArrayToNamedObject(monitorNames),
            this._monitors,
        );
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

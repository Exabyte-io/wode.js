import { UNIT_TYPES } from "../../enums";
import { UnitConfigBuilder } from "./UnitConfigBuilder";

export class AssignmentUnitConfigBuilder extends UnitConfigBuilder {
    constructor(name, variableName, variableValue, input = [], results = []) {
        super({ name, type: UNIT_TYPES.assignment });
        this._variableName = variableName;
        this._variableValue = variableValue;
        this._input = input;
        this._results = results;
    }

    input(arr) {
        this._input = arr;
        return this;
    }

    variableName(str) {
        this._variableName = str;
        return this;
    }

    variableValue(str) {
        this._variableValue = str;
        return this;
    }

    build() {
        return {
            ...super.build(),
            input: this._input,
            operand: this._variableName,
            value: this._variableValue,
        };
    }
}

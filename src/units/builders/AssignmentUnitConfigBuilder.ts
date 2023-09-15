import { AssignmentUnitConfig, UNIT_TYPES, UnitInput } from "../types";
import { UnitConfigBuilder } from "./UnitConfigBuilder";

export class AssignmentUnitConfigBuilder extends UnitConfigBuilder<AssignmentUnitConfig> {
    private _operand: string;
    private _value: string;
    private _input: UnitInput;

    constructor(
        name: AssignmentUnitConfig["name"],
        operand: AssignmentUnitConfig["operand"],
        value: AssignmentUnitConfig["value"],
        input: AssignmentUnitConfig["input"] = []
    ) {
        super({ name, type: UNIT_TYPES.assignment });
        this._operand = operand;
        this._value = value;
        this._input = input;
    }

    input(arr: UnitInput): this {
        this._input = arr;
        return this;
    }

    variableName(str: string): this {
        this._operand = str;
        return this;
    }

    variableValue(str: string): this {
        this._value = str;
        return this;
    }

    build() {
        return {
            ...super.build(),
            input: this._input,
            operand: this._operand,
            value: this._value,
        };
    }
}

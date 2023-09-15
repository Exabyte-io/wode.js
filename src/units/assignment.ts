import { BaseUnit } from "./base";
import { AssignmentUnitConfig, UNIT_TYPES } from "./types";

export class AssignmentUnit extends BaseUnit<AssignmentUnitConfig> {
    constructor(config: AssignmentUnitConfig) {
        super({ ...AssignmentUnit.getAssignmentConfig(), ...config });
    }

    static getAssignmentConfig() {
        return {
            name: UNIT_TYPES.assignment,
            type: UNIT_TYPES.assignment,
            operand: "X",
            value: 1,
            input: [],
        };
    }

    get operand() {
        return this.prop("operand");
    }

    get value() {
        return this.prop("value");
    }

    get input() {
        return this.prop("input");
    }

    getHashObject(): Partial<AssignmentUnitConfig> {
        return { input: this.input, operand: this.operand, value: this.value };
    }
}

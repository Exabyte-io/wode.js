import { UNIT_TYPES } from "../enums";
import { BaseUnit } from "./base";

export class AssignmentUnit extends BaseUnit {
    constructor(config) {
        super({ ...AssignmentUnit.getAssignmentConfig(), ...config });
    }

    static getAssignmentConfig() {
        return {
            name: UNIT_TYPES.assignment,
            type: UNIT_TYPES.assignment,
            operand: "X",
            value: "1",
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

    getHashObject() {
        return { input: this.input, operand: this.operand, value: this.value };
    }
}

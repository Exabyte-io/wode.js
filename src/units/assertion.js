import { UNIT_TYPES } from "../enums";
import { BaseUnit } from "./base";

export class AssertionUnit extends BaseUnit {
    constructor(config) {
        super({ ...AssertionUnit.getAssertionConfig(), ...config });
    }

    static getAssertionConfig() {
        return {
            name: UNIT_TYPES.assertion,
            type: UNIT_TYPES.assertion,
            statement: "true",
            errorMessage: "assertion failed",
        };
    }

    get statement() {
        return this.prop("statement");
    }

    get errorMessage() {
        return this.prop("errorMessage");
    }

    getHashObject() {
        return { statement: this.statement, errorMessage: this.errorMessage };
    }
}

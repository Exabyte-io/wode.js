import { BaseUnit } from "./base";
import { AssertionUnitConfig, UNIT_TYPES } from "./types";

export class AssertionUnit extends BaseUnit<AssertionUnitConfig> {
    constructor(config: AssertionUnitConfig) {
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

    getHashObject(): Partial<AssertionUnitConfig> {
        return { statement: this.statement, errorMessage: this.errorMessage };
    }
}

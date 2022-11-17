import { UNIT_TYPES } from "../enums";
import { BaseUnit } from "./base";

export class ConditionUnit extends BaseUnit {
    constructor(config) {
        super({ ...ConditionUnit.getConditionConfig(), ...config });
    }

    static getConditionConfig() {
        return {
            name: UNIT_TYPES.condition,
            type: UNIT_TYPES.condition,
            input: [],
            results: [],
            preProcessors: [],
            postProcessors: [],
            then: undefined,
            else: undefined,
            statement: "true",
            maxOccurrences: 100,
        };
    }

    get input() {
        return this.prop("input");
    }

    get then() {
        return this.prop("then");
    }

    get else() {
        return this.prop("else");
    }

    get statement() {
        return this.prop("statement");
    }

    get maxOccurrences() {
        return this.prop("maxOccurrences");
    }

    getHashObject() {
        return { statement: this.statement, maxOccurrences: this.maxOccurrences };
    }
}

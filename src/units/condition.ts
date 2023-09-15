import { BaseUnit } from "./base";
import { ConditionUnitConfig, UNIT_TYPES } from "./types";

export class ConditionUnit extends BaseUnit<ConditionUnitConfig> {
    constructor(config: ConditionUnitConfig) {
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

    get then() {
        return this.prop("then");
    }

    get else(){
        return this.prop("else");
    }

    get statement() {
        return this.prop("statement");
    }

    get maxOccurrences() {
        return this.prop("maxOccurrences");
    }

    getHashObject(): Partial<ConditionUnitConfig> {
        return { statement: this.statement, maxOccurrences: this.maxOccurrences };
    }
}

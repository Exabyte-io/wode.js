import { UNIT_TYPES } from "../enums";
import { BaseUnit } from "./base";
import { SubworkflowUnitConfig } from "./types";

export class SubworkflowUnit extends BaseUnit<SubworkflowUnitConfig> {
    constructor(config: SubworkflowUnitConfig) {
        super({ ...SubworkflowUnit.getSubworkflowConfig(), ...config });
    }

    static getSubworkflowConfig() {
        return {
            name: "New Subworkflow",
            type: UNIT_TYPES.subworkflow,
        };
    }
}

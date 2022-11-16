import { UNIT_TYPES } from "../enums";
import { BaseUnit } from "./base";

export class SubworkflowUnit extends BaseUnit {
    constructor(config) {
        super({ ...SubworkflowUnit.getSubworkflowConfig(), ...config });
    }

    static getSubworkflowConfig() {
        return {
            name: "New Subworkflow",
            type: UNIT_TYPES.subworkflow,
        };
    }
}

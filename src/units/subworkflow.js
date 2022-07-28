import { BaseUnit } from "./base";

import { UNIT_TYPES } from "/imports/wode/enums";


export class SubworkflowUnit extends BaseUnit {

    constructor(config) {
        super({ ...SubworkflowUnit.getSubworkflowConfig(), ...config });
    }

    static getSubworkflowConfig() {
        return {
            name: "New Subworkflow",
            type: UNIT_TYPES.subworkflow,
        }
    }

}

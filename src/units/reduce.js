import { UNIT_TYPES } from "../enums";
import { BaseUnit } from "./base";

export class ReduceUnit extends BaseUnit {
    constructor(unitName, mapUnit, input) {
        super({ ...ReduceUnit.getReduceConfig(unitName, mapUnit, input) });
    }

    static getReduceConfig(unitName, mapUnit, input) {
        return {
            type: UNIT_TYPES.reduce,
            name: unitName,
            mapFlowchartId: mapUnit,
            input,
        };
    }
}

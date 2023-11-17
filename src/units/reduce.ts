import { UNIT_TYPES } from "../enums";
import { BaseUnit } from "./base";
import { ReduceUnitConfig } from "./types";

export class ReduceUnit extends BaseUnit<ReduceUnitConfig> {
    constructor({name, mapFlowchartId, input, ...config}: ReduceUnitConfig) {
        super({...ReduceUnit.getReduceConfig(name, mapFlowchartId, input), ...config });
    }

    static getReduceConfig(name, mapFlowchartId, input) {
        return {
            type: UNIT_TYPES.reduce,
            name: name,
            mapFlowchartId: mapFlowchartId,
            input,
        };
    }
}

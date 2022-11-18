import { UNIT_TYPES } from "../enums";
import { AssertionUnit } from "./assertion";
import { AssignmentUnit } from "./assignment";
import { BaseUnit } from "./base";
import { ConditionUnit } from "./condition";
import { ExecutionUnit } from "./execution";
import { IOUnit } from "./io";
import { MapUnit } from "./map";
import { ProcessingUnit } from "./processing";
import { SubworkflowUnit } from "./subworkflow";

export class UnitFactory {
    static create(config) {
        switch (config.type) {
            case UNIT_TYPES.execution:
                return new ExecutionUnit(config);
            case UNIT_TYPES.assignment:
                return new AssignmentUnit(config);
            case UNIT_TYPES.condition:
                return new ConditionUnit(config);
            case UNIT_TYPES.io:
                return new IOUnit(config);
            case UNIT_TYPES.processing:
                return new ProcessingUnit(config);
            case UNIT_TYPES.map:
                return new MapUnit(config);
            case UNIT_TYPES.subworkflow:
                return new SubworkflowUnit(config);
            case UNIT_TYPES.assertion:
                return new AssertionUnit(config);
            default:
                return new BaseUnit(config);
        }
    }
}

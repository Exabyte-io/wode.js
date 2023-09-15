import { AssertionUnit } from "./assertion";
import { AssignmentUnit } from "./assignment";
import { BaseUnit } from "./base";
import { ConditionUnit } from "./condition";
import { ExecutionUnit } from "./execution";
import { IOUnit } from "./io";
import { MapUnit } from "./map";
import { ProcessingUnit } from "./processing";
import { SubworkflowUnit } from "./subworkflow";
import {
    AssertionUnitConfig,
    AssignmentUnitConfig,
    ConditionUnitConfig,
    IOUnitConfig,
    MapUnitConfig,
    ProcessingUnitConfig,
    SubworkflowUnitConfig,
    UNIT_TYPES,
    UnitConfig
} from "./types";

export class UnitFactory {
    static create(config: UnitConfig) {
        switch (config.type) {
            case UNIT_TYPES.execution:
                return new ExecutionUnit();
            case UNIT_TYPES.assignment:
                return new AssignmentUnit(config as AssignmentUnitConfig);
            case UNIT_TYPES.condition:
                return new ConditionUnit(config as ConditionUnitConfig);
            case UNIT_TYPES.io:
                return new IOUnit(config as IOUnitConfig);
            case UNIT_TYPES.processing:
                return new ProcessingUnit(config as ProcessingUnitConfig);
            case UNIT_TYPES.map:
                return new MapUnit(config as MapUnitConfig);
            case UNIT_TYPES.subworkflow:
                return new SubworkflowUnit(config as SubworkflowUnitConfig);
            case UNIT_TYPES.assertion:
                return new AssertionUnit(config as AssertionUnitConfig);
            default:
                return new BaseUnit(config as UnitConfig);
        }
    }
}

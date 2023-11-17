import { UNIT_TYPES } from "../enums";
import { BaseUnit } from "./base";
import { ProcessingUnitConfig, UnitInput } from "./types";

export class ProcessingUnit extends BaseUnit<ProcessingUnitConfig> {
    constructor(config: ProcessingUnitConfig) {
        super({ ...ProcessingUnit.getProcessingConfig(), ...config });
    }

    static getProcessingConfig() {
        return {
            name: UNIT_TYPES.processing,
            type: UNIT_TYPES.processing,
        };
    }

    setOperation(op: string) {
        this.setProp("operation", op);
    }

    setOperationType(type: string) {
        this.setProp("operationType", type);
    }

    setInput(input: UnitInput) {
        this.setProp("input", input);
    }

    get operation() {
        return this.prop("operation");
    }

    get operationType() {
        return this.prop("operationType");
    }

    get input() {
        return this.prop("input");
    }
}

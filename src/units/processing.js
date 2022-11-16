import { UNIT_TYPES } from "../enums";
import { BaseUnit } from "./base";

export class ProcessingUnit extends BaseUnit {
    constructor(config) {
        super({ ...ProcessingUnit.getProcessingConfig(), ...config });
    }

    static getProcessingConfig() {
        return {
            name: UNIT_TYPES.processing,
            type: UNIT_TYPES.processing,
        };
    }

    setOperation(op) {
        this.setProp("operation", op);
    }

    setOperationType(type) {
        this.setProp("operationType", type);
    }

    setInput(input) {
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

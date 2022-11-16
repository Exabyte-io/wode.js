import { UNIT_TYPES } from "../enums";
import { BaseUnit } from "./base";

export const defaultMapConfig = {
    name: UNIT_TYPES.map,
    type: UNIT_TYPES.map,
    workflowId: "",
    input: {
        target: "MAP_DATA",
        scope: "global",
        name: "",
        values: [],
        useValues: false,
    },
};

export class MapUnit extends BaseUnit {
    constructor(config) {
        super({ ...defaultMapConfig, ...config });
    }

    get input() {
        return this.prop("input");
    }

    get workflowId() {
        return this.prop("workflowId");
    }

    setWorkflowId(id) {
        this.setProp("workflowId", id);
    }
}

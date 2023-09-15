import { UNIT_TYPES } from "../enums";
import { BaseUnit } from "./base";
import { MapUnitConfig } from "./types";

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

export class MapUnit extends BaseUnit<MapUnitConfig> {
    constructor(config: MapUnitConfig) {
        super({ ...defaultMapConfig, ...config });
    }

    get input() {
        return this.prop("input");
    }

    get workflowId() {
        return this.prop("workflowId");
    }

    setWorkflowId(id: string) {
        this.setProp("workflowId", id);
    }
}

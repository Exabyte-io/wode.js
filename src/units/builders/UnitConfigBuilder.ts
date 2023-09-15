import { getUUID } from "@exabyte-io/code.js/dist/utils";
import _ from "underscore";
import { NamedEntity, UNIT_TYPES, UnitConfig } from "../types";

export class UnitConfigBuilder<T extends UnitConfig> {
    private _name: string;
    private _flowchartId: string;
    private _type: UNIT_TYPES;
    private _head: boolean;

    constructor({ name, type, flowchartId }: Pick<T, "type" | "name"> & Partial<Pick<T, "flowchartId">>) {
        this._type = type;
        this._name = name;
        this._head = false;
        this._flowchartId = flowchartId || UnitConfigBuilder.generateFlowChartId();
    }

    static generateFlowChartId(): string {
        return getUUID();
    }

    static _stringArrayToNamedObject(array: string[]): NamedEntity[] {
        return array.map((name) => ({ name }));
    }

    name(str: string): this {
        this._name = str;
        return this;
    }

    head(bool: boolean): this {
        this._head = bool;
        return this;
    }

    flowchartId(flowchartId: string): this {
        this._flowchartId = flowchartId;
        return this;
    }

    build(): T {
        return {
            type: this._type,
            name: this._name,
            head: this._head,
            flowchartId: this._flowchartId,
        } as T;
    }
}

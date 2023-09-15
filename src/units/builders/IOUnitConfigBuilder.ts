import { IOUnitConfig, IOUnitSubTypes, UNIT_TYPES } from "../types";
import { UnitConfigBuilder } from "./UnitConfigBuilder";

export class IOUnitConfigBuilder extends UnitConfigBuilder<IOUnitConfig> {
    private _endpointName: string;
    private _endpointOptions: string;
    private _variableName: string;
    private _subtype: IOUnitSubTypes;
    private _source: string;

    constructor(
        name: IOUnitConfig["name"],
        endpointName: string,
        endpointOptions: string,
    ) {
        super({ name, type: UNIT_TYPES.io });
        this._endpointName = endpointName;
        this._endpointOptions = endpointOptions;
        this._variableName = "DATA";
        this._subtype = "input";
        this._source = "api";
    }

    endpointName(str: string): this {
        this._endpointName = str;
        return this;
    }

    endpointOptions(options: string): this {
        this._endpointOptions = options;
        return this;
    }

    variableName(str: string): this {
        this._variableName = str;
        return this;
    }

    subtype(str: IOUnitSubTypes): this {
        this._subtype = str;
        return this;
    }

    source(str: string): this {
        this._source = str;
        return this;
    }

    build(): IOUnitConfig {
        return {
            ...super.build(),
            subtype: this._subtype,
            source: this._source,
            input: [
                {
                    endpoint: this._endpointName,
                    endpoint_options: this._endpointOptions,
                    name: this._variableName,
                },
            ],
        };
    }
}

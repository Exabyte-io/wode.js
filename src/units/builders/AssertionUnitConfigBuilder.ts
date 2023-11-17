import { AssertionUnitConfig, UNIT_TYPES } from "../types";
import { UnitConfigBuilder } from "./UnitConfigBuilder";

export class AssertionUnitConfigBuilder extends UnitConfigBuilder<AssertionUnitConfig> {
    private _statement: string;
    private _errorMessage: string;

    constructor(
        name: AssertionUnitConfig["name"],
        statement: AssertionUnitConfig["statement"],
        errorMessage: AssertionUnitConfig["errorMessage"]
    ){
        super({ name, type: UNIT_TYPES.assertion });
        this._statement = statement;
        this._errorMessage = errorMessage;
    }

    statement(str: string): this {
        this._statement = str;
        return this;
    }

    errorMessage(str: string): this {
        this._errorMessage = str;
        return this;
    }

    build(): AssertionUnitConfig {
        return {
            ...super.build(),
            statement: this._statement,
            errorMessage: this._errorMessage,
        };
    }
}

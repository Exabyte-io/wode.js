import { UNIT_TYPES } from "../../enums";
import { UnitConfigBuilder } from "./UnitConfigBuilder";

export class AssertionUnitConfigBuilder extends UnitConfigBuilder {
    constructor(name, statement, errorMessage) {
        super({ name, type: UNIT_TYPES.assertion });
        this._statement = statement;
        this._errorMessage = errorMessage;
    }

    statement(str) {
        this._statement = str;
        return this;
    }

    errorMessage(str) {
        this._errorMessage = str;
        return this;
    }

    build() {
        return {
            ...super.build(),
            statement: this._statement,
            errorMessage: this._errorMessage,
        };
    }
}

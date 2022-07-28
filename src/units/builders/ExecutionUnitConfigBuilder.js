import { Application, Executable, Flavor } from "@exabyte-io/ade.js";
import { UNIT_TYPES } from "/imports/wode/enums";

import { UnitConfigBuilder } from "./UnitConfigBuilder";


export class ExecutionUnitConfigBuilder extends UnitConfigBuilder {
    static Application = Application;
    static Executable = Executable;
    static Flavor = Flavor;

    constructor(name, application, execName, flavorName) {
        super({name, type: UNIT_TYPES.execution});

        try {
            this.initialize(application, execName, flavorName);
        } catch (e) {
            console.error(`Can't initialize executable/flavor: ${execName}/${flavorName}`);
            throw e;
        }

        // initialize runtimeItems
        this._results = this.flavor.results;
        this._monitors = this.flavor.monitors;
        this._preProcessors = this.flavor.preProcessors;
        this._postProcessors = this.flavor.postProcessors;
    }

    initialize(application, execName, flavorName) {
        this.application = application;
        this.executable = this.constructor.Executable.create({name: execName, application: this.application});
        this.flavor = this.constructor.Flavor.create({name: flavorName, executable: this.executable});
    }

    build() {
        return {
            ...super.build(),
            application: this.application.toJSON(),
            executable: this.executable.toJSON(),
            flavor: this.flavor.toJSON(),
        };
    }
}

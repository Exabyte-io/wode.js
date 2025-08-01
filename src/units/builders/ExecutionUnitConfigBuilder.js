/* eslint-disable class-methods-use-this */
import ApplicationRegistry from "@exabyte-io/ade.js/dist/js/ApplicationRegistry";

import { UNIT_TYPES } from "../../enums";
import { UnitConfigBuilder } from "./UnitConfigBuilder";

export class ExecutionUnitConfigBuilder extends UnitConfigBuilder {
    constructor(name, application, execName, flavorName, flowchartId) {
        super({ name, type: UNIT_TYPES.execution, flowchartId });

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
        this.executable = this._createExecutable(this.application, execName);
        this.flavor = this._createFlavor(this.executable, flavorName);
    }

    build() {
        return {
            ...super.build(),
            application: this.application.toJSON(),
            executable: this.executable.toJSON(),
            flavor: this.flavor.toJSON(),
        };
    }

    /**
     * Creates an executable instance. This method is intended to be overridden in subclasses.
     * @param {Application} application - The application object
     * @param {string} execName - The name of the executable
     * @returns {Executable} The created executable instance
     */
    _createExecutable(application, execName) {
        return ApplicationRegistry.getExecutableByName(application.name, execName);
    }

    /**
     * Creates a flavor instance. This method is intended to be overridden in subclasses.
     * @param {Executable} executable - The executable object
     * @param {string} flavorName - The name of the flavor
     * @returns {Flavor} The created flavor instance
     */
    _createFlavor(executable, flavorName) {
        return ApplicationRegistry.getFlavorByName(executable, flavorName);
    }
}

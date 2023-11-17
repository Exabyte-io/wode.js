import { Application, Executable, Flavor } from "@exabyte-io/ade.js";
import { UnitConfigBuilder } from "./UnitConfigBuilder";
import { ExecutionUnitConfig, UNIT_TYPES } from "../types";

export class ExecutionUnitConfigBuilder extends UnitConfigBuilder<ExecutionUnitConfig> {

    private flavor: ExecutionUnitConfig["flavor"];
    private executable: ExecutionUnitConfig["executable"];
    private application: ExecutionUnitConfig["application"];
    private _results: ExecutionUnitConfig["flavor"]["results"];
    private _monitors: ExecutionUnitConfig["flavor"]["monitors"];
    private _preProcessors: ExecutionUnitConfig["flavor"]["preProcessors"];
    private _postProcessors: ExecutionUnitConfig["flavor"]["postProcessors"];

    static Application = Application;
    static Executable = Executable;
    static Flavor = Flavor;

    constructor(
        name : ExecutionUnitConfig["name"],
        application: ExecutionUnitConfig["application"],
        execName: ExecutionUnitConfig["executable"]["name"],
        flavorName: ExecutionUnitConfig["flavor"]["name"],
        flowchartId: ExecutionUnitConfig["flowchartId"]
    ) {
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

    initialize(
        application: ExecutionUnitConfig["application"],
        execName: ExecutionUnitConfig["executable"]["name"],
        flavorName: ExecutionUnitConfig["flavor"]["name"]) {
        this.application = application;
        this.executable = ExecutionUnitConfigBuilder.Executable.create({
            name: execName,
            application: this.application,
        });
        this.flavor = ExecutionUnitConfigBuilder.Flavor.create({
            name: flavorName,
            executable: this.executable,
        });
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

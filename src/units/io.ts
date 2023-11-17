import lodash from "lodash";

import { IO_ID_COLUMN, UNIT_TYPES } from "../enums";
import { BaseUnit } from "./base";
import { IOUnitConfig } from "./types";

export class IOUnit extends BaseUnit<IOUnitConfig> {
    private _materials: any[]; // TODO: define Material type
    private _defaultTargets: string[];
    private _features: string[];
    private _targets: string[];
    private _ids: string[];
    private _jobId: string | null;
    public dataGridValues: any[];
    private clean: <T>(obj: T) => T; // TODO: should be coming from mixins
    // public toJSON: () => T;


    constructor(config: IOUnitConfig) {
        super({ ...IOUnit.getIOConfig(), ...config });
        this.initialize(config);
    }

    static getIOConfig() {
        return {
            name: UNIT_TYPES.io,
            type: UNIT_TYPES.io,
            subtype: "input",
        };
    }

    initialize(config: IOUnitConfig) {
        this._materials = [];
        this._defaultTargets = ["band_gaps:direct", "band_gaps:indirect"];
        this._features = lodash.get(config, "input.0.endpoint_options.data.features", []);
        this._targets = lodash.get(
            config,
            "input.0.endpoint_options.data.targets",
            this._defaultTargets,
        );
        this._ids = lodash.get(config, "input.0.endpoint_options.data.ids", []);
        this._jobId = null;
    }

    get materials() { // TODO: define Material type
        return this._materials || [];
    }

    get defaultTargets() {
        return this._defaultTargets;
    }

    get features() {
        return this._features;
    }

    get featuresWithoutId() {
        return this.features.filter((x) => x !== IO_ID_COLUMN);
    }

    // TODO: refactor to not use lodash
    get availableFeatures(): any[] {
        const { materials } = this;
        return lodash.uniq(
            lodash
                .flatten(materials.map((x) => lodash.keys(x.propertiesDict())))
                .concat(this.features),
        );
    }

    get availableFeaturesWithoutId() {
        return this.availableFeatures.filter((feature) => feature !== IO_ID_COLUMN);
    }

    get targets() {
        return this._targets;
    }

    /**
     * @summary Checks whether selected features contain only IO_ID_COLUMN ('exabyteId').
     * Used to identify that no features are selected yet (features set always contains ID_COLUMN)
     */
    get onlyIdFeatureSelected(): boolean {
        return lodash.isEmpty(lodash.without(this.features, IO_ID_COLUMN));
    }

    /**
     * @summary Returns object with targets as key and arrays of appropriate values.
     * E.g. {'band_gap:indirect': [0.1, 0.3], 'pressure': [100, undefined]}
     */
    get valuesByTarget() {
        const values = this.dataGridValues;
        const result = {};
        this.targets.forEach((target) => {
            result[target] = values.map((v) => v[target]);
        });
        return result;
    }

    get dataFrameConfig() {
        return {
            subtype: "dataFrame",
            source: "api",
            input: [
                {
                    endpoint: "dataframe",
                    endpoint_options: {
                        method: "POST",
                        data: {
                            targets: this._targets,
                            features: this._features,
                            ids: this._ids,
                            jobId: this._jobId,
                        },
                        headers: {},
                        params: {},
                    },
                },
            ],
        };
    }

    get isDataFrame() {
        return this.prop("subtype") === "dataFrame";
    }

    setMaterials(materials: any[]) { //TODO: define Material type
        this._materials = materials;
        this._ids = materials.map((m) => m.exabyteId);
    }

    addFeature(feature: string) {
        // only add if not already present
        if (this._features.indexOf(feature) === -1) this._features.push(feature);
    }

    removeFeature(feature: string) {
        if (this.featuresWithoutId.length === 1) {
            throw new Error("At least one feature is required");
        }
        this._features = this._features.filter((x) => feature !== x && x !== IO_ID_COLUMN);
    }

    addTarget(target: string) {
        if (this._targets.indexOf(target) === -1) this._targets.push(target);
    }

    removeTarget(target: string) {
        if (this._targets.length === 1) {
            throw new Error("At least one target is required");
        }
        this._targets = this._targets.filter((x) => target !== x);
    }

    hasFeature(feature: string) {
        return this._features.indexOf(feature) > -1;
    }

    hasTarget(target: string) {
        return this._targets.indexOf(target) > -1;
    }

    toJSON(): IOUnitConfig {
        const config = this.isDataFrame ? this.dataFrameConfig : {};
        return this.clean({ ...super.toJSON(), ...config });
    }
}

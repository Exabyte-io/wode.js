import lodash from "lodash";

import { IO_ID_COLUMN, UNIT_TYPES } from "../enums";
import { BaseUnit } from "./base";

export class IOUnit extends BaseUnit {
    /**
     * IO Unit Builder for Object Storage sources.
     *
     * @param {Object} config - config object with other parameters:
     * @param {String} config.name - the name of the unit this builder is creating
     * @param {String} config.subtype - "input", "output", or "dataframe"
     * @param {Object} config.input - input containing information on the file to download
     * @param {Boolean} config.enableRender - Whether to use Jinja templating at runtime
     */
    constructor(config) {
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

    initialize(config) {
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

    get materials() {
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

    get availableFeatures() {
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
    get onlyIdFeatureSelected() {
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

    setMaterials(materials) {
        this._materials = materials;
        this._ids = materials.map((m) => m.exabyteId);
    }

    addFeature(feature) {
        // only add if not already present
        if (this._features.indexOf(feature) === -1) this._features.push(feature);
    }

    removeFeature(feature) {
        if (this.featuresWithoutId.length === 1) {
            throw new Error("At least one feature is required");
        }
        this._features = this._features.filter((x) => feature !== x && x !== IO_ID_COLUMN);
    }

    addTarget(target) {
        if (this._targets.indexOf(target) === -1) this._targets.push(target);
    }

    removeTarget(target) {
        if (this._targets.length === 1) {
            throw new Error("At least one target is required");
        }
        this._targets = this._targets.filter((x) => target !== x);
    }

    hasFeature(feature) {
        return this._features.indexOf(feature) > -1;
    }

    hasTarget(target) {
        return this._targets.indexOf(target) > -1;
    }

    toJSON() {
        const config = this.isDataFrame ? this.dataFrameConfig : {};
        return this.clean({ ...super.toJSON(), ...config });
    }
}

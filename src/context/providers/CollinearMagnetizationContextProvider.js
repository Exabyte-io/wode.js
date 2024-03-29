import { JSONSchemaFormDataProvider, MaterialContextMixin } from "@mat3ra/code/dist/js/context";
import { Made } from "@mat3ra/made";
import lodash from "lodash";
import { mix } from "mixwith";

export class CollinearMagnetizationContextProvider extends mix(JSONSchemaFormDataProvider).with(
    MaterialContextMixin,
) {
    static Material = Made.Material;

    constructor(config) {
        super(config);
        this.firstElement =
            this.uniqueElementsWithLabels?.length > 0 ? this.uniqueElementsWithLabels[0] : "";
        this.is_constrained_magnetization = lodash.get(
            this.data,
            "is_constrained_magnetization",
            false,
        );
    }

    get uniqueElementsWithLabels() {
        const elementsWithLabelsArray = this.material?.Basis?.elementsWithLabelsArray || [];
        return [...new Set(elementsWithLabelsArray)];
    }

    indexOfElement = (element) => {
        return this.uniqueElementsWithLabels.indexOf(element) + 1;
    };

    // eslint-disable-next-line class-methods-use-this
    get defaultData() {
        return {
            starting_magnetization: [
                {
                    index: 1,
                    atomicSpecies: this.firstElement,
                    value: 0.0,
                },
            ],
            is_constrained_magnetization: false,
            total_magnetization: 0.0,
        };
    }

    transformData = (data) => {
        const startingMagnetizationWithIndex = data.starting_magnetization.map((row) => ({
            ...row,
            index: this.indexOfElement(row.atomicSpecies),
        }));

        return {
            ...data,
            starting_magnetization: startingMagnetizationWithIndex,
        };
    };

    // eslint-disable-next-line class-methods-use-this
    get uiSchemaStyled() {
        return {
            starting_magnetization: {
                items: {
                    atomicSpecies: {
                        "ui:classNames": "col-xs-3",
                    },
                    value: {
                        "ui:classNames": "col-xs-6",
                    },
                },
                "ui:readonly": this.is_constrained_magnetization,
            },
            is_constrained_magnetization: {},
            total_magnetization: {
                "ui:classNames": "col-xs-6",
                "ui:readonly": !this.is_constrained_magnetization,
            },
        };
    }

    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-07/schema#",
            title: "",
            description: "Set starting magnetization, can have values in the range [-1, +1].",
            type: "object",
            properties: {
                starting_magnetization: {
                    type: "array",
                    items: {
                        type: "object",
                        properties: {
                            atomicSpecies: {
                                type: "string",
                                title: "Atomic species",
                                enum: this.uniqueElementsWithLabels,
                                default: this.firstElement,
                            },
                            value: {
                                type: "number",
                                title: "Starting magnetization",
                                default: 0.0,
                                minimum: -1.0,
                                maximum: 1.0,
                            },
                        },
                    },
                },
                is_constrained_magnetization: {
                    type: "boolean",
                    title: "Set constrained magnetization instead",
                    default: false,
                },
                total_magnetization: {
                    type: "number",
                    title: "Total magnetization",
                    default: 0.0,
                },
            },
        };
    }
}

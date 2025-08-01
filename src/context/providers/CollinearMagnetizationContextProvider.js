import JSONSchemaFormDataProvider from "@exabyte-io/ade.js/dist/js/context/JSONSchemaFormDataProvider";
import lodash from "lodash";

import { materialContextMixin } from "../mixins/MaterialContextMixin";

export class CollinearMagnetizationContextProvider extends JSONSchemaFormDataProvider {
    constructor(config) {
        super(config);

        this.initMaterialContextMixin();

        this.firstElement =
            this.uniqueElementsWithLabels?.length > 0 ? this.uniqueElementsWithLabels[0] : "";
        this.isTotalMagnetization = lodash.get(this.data, "isTotalMagnetization", false);
    }

    get uniqueElementsWithLabels() {
        const elementsWithLabelsArray = this.material?.Basis?.elementsWithLabelsArray || [];
        return [...new Set(elementsWithLabelsArray)];
    }

    indexOfElement = (element) => {
        return this.uniqueElementsWithLabels.indexOf(element) + 1;
    };

    get defaultData() {
        return {
            startingMagnetization: [
                {
                    index: 1,
                    atomicSpecies: this.firstElement,
                    value: 0.0,
                },
            ],
            isTotalMagnetization: false,
            totalMagnetization: 0.0,
        };
    }

    transformData = (data) => {
        const startingMagnetizationWithIndex = data.startingMagnetization.map((row) => ({
            ...row,
            index: this.indexOfElement(row.atomicSpecies),
        }));

        return {
            ...data,
            startingMagnetization: startingMagnetizationWithIndex,
        };
    };

    get uiSchemaStyled() {
        return {
            startingMagnetization: {
                items: {
                    atomicSpecies: {
                        "ui:classNames": "col-xs-3",
                    },
                    value: {
                        "ui:classNames": "col-xs-6",
                    },
                },
                "ui:readonly": this.isTotalMagnetization,
            },
            isTotalMagnetization: {},
            totalMagnetization: {
                "ui:classNames": "col-xs-6",
                "ui:readonly": !this.isTotalMagnetization,
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
                startingMagnetization: {
                    type: "array",
                    maxItems: this.uniqueElementsWithLabels.length,
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
                isTotalMagnetization: {
                    type: "boolean",
                    title: "Set total magnetization instead",
                    default: false,
                },
                totalMagnetization: {
                    type: "number",
                    title: "Total magnetization",
                    default: 0.0,
                },
            },
        };
    }
}

materialContextMixin(CollinearMagnetizationContextProvider.prototype);

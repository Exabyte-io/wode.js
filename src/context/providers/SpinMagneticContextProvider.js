import { JSONSchemaFormDataProvider, MaterialContextMixin } from "@exabyte-io/code.js/dist/context";
import { Made } from "@exabyte-io/made.js";
import { mix } from "mixwith";

export class SpinMagneticContextProvider extends mix(JSONSchemaFormDataProvider).with(
    MaterialContextMixin,
) {
    static Material = Made.Material;

    constructor(config) {
        super(config);
        this.atomicLabelsArray = this.material?.Basis?.atomicLabelsArray || [];
        this.elementsArray = this.material?.Basis?.elementsArray || [];
    }

    get uniqueElementsWithLabels() {
        const elementsWithLabelsArray = [];
        this.elementsArray.forEach((item, idx) =>
            elementsWithLabelsArray.push(item + this.atomicLabelsArray[idx]),
        );
        return [...new Set(elementsWithLabelsArray)];
    }

    indexOfElement = (element) => {
        return this.uniqueElementsWithLabels.indexOf(element) + 1;
    };

    // eslint-disable-next-line class-methods-use-this
    get defaultData() {
        return [
            {
                index: 1,
                atomicSpecies:
                    this.uniqueElementsWithLabels?.length > 0
                        ? this.uniqueElementsWithLabels[0]
                        : "",
                starting_magnetization: 0.0,
            },
        ];
    }

    transformData = (data) => {
        return data.map((row) => ({
            ...row,
            index: this.indexOfElement(row.atomicSpecies),
        }));
    };

    get uiSchemaStyled() {
        return {
            "ui:options": {
                addable: true,
                orderable: false,
                removable: true,
            },
            title: {
                "ui:classNames": "col-xs-12",
            },
            items: {
                atomicSpecies: this.defaultFieldStyles,
                value: this.defaultFieldStyles,
            },
        };
    }

    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-04/schema#",
            title: "",
            description: "Set starting magnetization.",
            type: "array",
            items: {
                type: "object",
                properties: {
                    atomicSpecies: {
                        type: "string",
                        title: "Atomic species",
                        enum: this.uniqueElementsWithLabels,
                        default:
                            this.uniqueElementsWithLabels?.length > 0
                                ? this.uniqueElementsWithLabels[0]
                                : "",
                    },
                    value: {
                        type: "number",
                        title: "Starting magnetization",
                        default: 0.0,
                    },
                },
            },
        };
    }
}

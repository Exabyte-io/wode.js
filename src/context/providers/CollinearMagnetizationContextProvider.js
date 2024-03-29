import { JSONSchemaFormDataProvider, MaterialContextMixin } from "@mat3ra/code/dist/js/context";
import { Made } from "@mat3ra/made";
import { mix } from "mixwith";

export class CollinearMagnetizationContextProvider extends mix(JSONSchemaFormDataProvider).with(
    MaterialContextMixin,
) {
    static Material = Made.Material;

    get uniqueElementsWithLabels() {
        const elementsWithLabelsArray = this.material?.Basis?.elementsWithLabelsArray || [];
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
                value: 0.0,
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
            items: {
                atomicSpecies: this.defaultFieldStyles,
                value: {
                    "ui:classNames": "col-xs-6 ",
                },
            },
        };
    }

    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-04/schema#",
            title: "",
            description: "Set starting magnetization, can have values in the range [-1, +1].",
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
                        minimum: -1.0,
                        maximum: 1.0,
                    },
                },
            },
        };
    }
}
